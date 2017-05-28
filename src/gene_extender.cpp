/*****************************************************************************
 * This is the main of the CRAIG program. It predicts multi-exonic gene
 * structures. Given a gene model, which was previously trained with
 * train_craig and a fasta file it outputs a set of predictions for the
 * input fasta file
 *
 *   Copyright (C) 2002-2007  Axel E. Bernal (abernal@seas.upenn.edu)
 *
 *   This program is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU General Public License as
 *   published by the Free Software Foundation; either version 2 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful, but
 *   WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 *   02111-1307, USA.
 *
 *   The GNU General Public License is contained in the file COPYING.
 *
 ****************************************************************************/


#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "Utils.h"
#include "SequenceUtils.h"
#include "ContextIMM.h"
#include "Sequence.h"
#include "ResourceEngine.h"
#include "FeatureEngine.h"
#include "StructureCore.h"
#include "Lattice.h"
#include "GeneEvaluator.h"
#include "GeneTagPrinter.h"
#include "GeneUtils.h"
#include "InpFile.h"

#define NUM_PHASES 3
using namespace craig;

void printHelp(const char * pname) {
    cerr << "CRAIG v. " << CRAIG_VERSION << " discriminative  prediction tool for producing extended gene models(in development).\nWritten by Axel E. Bernal (abernal@seas.upenn.edu)\n\n";
    cerr << "  usage : " << pname << " [options] PARAMS_FILE FASTA_FILE > predictions\n\n";
    cerr << "  parameters is the name of the file containing the gene model parameters;\n             if file parameters is not found, craig assumes the filename\n             to be $(CRAIG_HOME)/models/PARAMS_FILE\n";
    cerr << "  fasta is the name of the file containing the input query sequence(s)\n        in fasta format\n\n";
    cerr << "  options:\n" ;
    cerr << "    --version\t\tPrint version name and license information\n";
    cerr << "    -h --help\t\tPrint this message\n";
    cerr << "    -v --verbose\tTurns on output of debugging information\n";
    cerr << "    --pfam=<pfam_input>\tFile containing pfam scores\n";
    cerr << "    --format=gtf|locs|tags\tOutput format for predictions. \"locs\" is an internal\n\t\t\tformat. \"tags\" is usually used for reranking [gtf]\n";
    cerr << "    --best=<number < 1000>\tOutputs the k-best segmentations [1]\n";
    cerr << "    -m --masked\t\tTreat the input sequences as if they had been\n\t\t\tpreviously masked. When this option is specified,\n\t\t\tlowercase letters MUST be used to represent\n\t\t\tthe masked-out regions\n";
    cerr << "    -c --complete\tAllows for prediction of complete genes only.\n\t\t\tIncomplete genes are predicted by default\n\n";
    cerr << "  predictions :";
    cerr << " The list of genes predicted in the input sequence\n\n";
    cerr << "  Report bugs to <abernal@seas.upenn.edu>\n\n";
}

bool verbose = false;

int main(int argc, char *argv[]) {
    std::string seqFileName("");
    std::string pfamFileName("");
    std::string paramFile("");
    std::string format("gtf");
    int topK = 1, i;
    std::string topology("partial");
    TStrand origStrand = BOTH_STRANDS;
    bool masked = false;

    try {

        for(i = 1; i < argc; i++) {

            if(!strncmp(argv[i], "--format=", 9)) {
                format = std::string(argv[i] + 9);
            }
            else if(!strncmp(argv[i], "--verbose", 9)
                    || !strncmp(argv[i], "-v", 2)) {
                verbose = true;
            }
            else if(!strncmp(argv[i], "--help", 6)
                    || !strncmp(argv[i], "-h", 2)) {
                printHelp("craig");
                exit(0);
            }
            else if(!strncmp(argv[i], "--version", 9)) {
                PRINT_VERSION(cerr, "craig", "discriminative gene prediction tool for eukarya");
                PRINT_DISCLAIMER(cerr, "craig");
                exit(0);
            }
            else if(!strncmp(argv[i], "--best=", 7)) {
                sscanf(argv[i] + 7, "%d", &topK);
            }
            else if(!strncmp(argv[i], "--masked", 8)
                    || !strncmp(argv[i], "-m", 2)) {
                masked = true;
            }
            else if(!strncmp(argv[i], "--pfam=", 7)) {
                pfamFileName = std::string(argv[i] + 7);
            }
            else if(!strncmp(argv[i], "--complete", 10)
                    || !strncmp(argv[i], "-c", 2)) {
                topology = std::string("complete");
            }
            else
                break;
        }

        if(argc - i < 2)
            throw EXCEPTION(BAD_USAGE, "insufficient arguments");

        paramFile = std::string(argv[i++]);
        seqFileName = std::string(argv[i++]);

        ::ifstream *paramStream,  ifstr1(paramFile.c_str()), ifstr2;
        paramStream = &ifstr1;

        if(!ifstr1.is_open()) {
            if(!getenv("CRAIG_HOME"))
                throw EXCEPTION(NOT_ANNOTATED,
                                "CRAIG_HOME must be initialized!. See README for details");

            paramFile = std::string(getenv("CRAIG_HOME")) + "/models/" + paramFile;
            ifstr2.open(paramFile.c_str());

            if(!ifstr2.is_open())
                throw EXCEPTION(FILE_UNAVAILABLE, paramFile);

            paramStream = &ifstr2;
        }

        // Retrieving resources

        ResourceEngine re(*paramStream);
        paramStream->clear();

        Sigma *sigma = (Sigma *)re.getResource(std::string("dna-alpha"));

        FSM & fsm = *(FSM *)re.getResource(topology);
        fsm.setParseStrand(origStrand);
        GeneEvaluator evaluator(fsm);

        std::string extTopology("extender-" + topology);
        FSM & bwfsm = *(FSM *)re.getResource(extTopology);
        bwfsm.setParseStrand(STRAND_COMP);
        GeneEvaluator bwevaluator(bwfsm);

        GeneTagPrinter printer;
        InpFile annotSeqFile("contigs", seqFileName, FASTA, sigma, true);

        //loading annotSeqs
        list<Sequence *> & annotSeqs = (list<Sequence *> &)annotSeqFile.sequences();

        FilterEngine fe(re, *paramStream);

        streampos fpos = paramStream->tellg();
        FeatureEngine fte(fsm, fe, *paramStream);
        paramStream->clear();
        paramStream->seekg(fpos);
        FeatureEngine bwfte(bwfsm, fe, *paramStream);

        GlobalVector params(fte.getFeatures(), paramStream);

        fte.setParamVector(&params);
        bwfte.setParamVector(&params);

        // Initializing special variables
        TypedFilter<EdgeInst> **signals = new TypedFilter<EdgeInst> * [NUM_EDGE_INSTS];
        signals[START] = (TypedFilter<EdgeInst> *)fe.getFilter("Start-Signal");
        signals[STOP] = (TypedFilter<EdgeInst> *)fe.getFilter("Stop-Signal");
        signals[DONOR] = (TypedFilter<EdgeInst> *)fe.getFilter("Donor-Signal");
        signals[ACCEPTOR] = (TypedFilter<EdgeInst> *)fe.getFilter("Acceptor-Signal");

        TypedFilter<UCHAR> *contexts = (TypedFilter<UCHAR> *)fe.findFilter("GC-Content");
        FT_HiddenSeq *vftHidden = (FT_HiddenSeq *)fte.findFeature("Hidden-Stop-Sequence");


        Lattice lattice(fe, fsm, NUM_PHASES, evaluator, MAX_NUM_EXONS,
                        INTERGENIC, signals, contexts, vftHidden);

        Lattice bwlattice(fe, bwfsm, NUM_PHASES, bwevaluator, MAX_NUM_EXONS,
                          INTERGENIC, signals, contexts, vftHidden);

        // turn OFF masking feature if input is not masked
        if(!masked) {
            Feature *maskFeature = fte.findFeature("Test-Is-Masked");
            if(maskFeature)
                maskFeature->turnOff();
            maskFeature = bwfte.findFeature("Test-Is-Masked");
            if(maskFeature)
                maskFeature->turnOff();
        }

        // Creating StructureCore object

        StructureCore core(fsm, lattice, re, fe, fte, evaluator, printer, AVG_ALL, COMB_EUCLID, (TStrand)origStrand, topK);
        StructureCore bwcore(bwfsm, bwlattice, re, fe, bwfte, bwevaluator, printer, AVG_ALL, COMB_EUCLID, STRAND_COMP, topK);
        cerr << "Done!\nProcessing annotSeqs...\n";
        printer.displayHeader(format, (std::ofstream &)cout);

        InpFile *pfamSeqFile = NULL;

        if(pfamFileName.length())
            pfamSeqFile = new InpFile("pfam", pfamFileName, EXTENDED_LDSCORE, NULL, true);

        list<Sequence *>::iterator cit = annotSeqs.begin();

        for( ; cit != annotSeqs.end(); cit++) {
            Sequence &superC = *(*cit);
            ScoreSeq<double> *pfamSeq = pfamSeqFile ?
                (ScoreSeq<double> *)pfamSeqFile->findSeq(superC.id()) :
                NULL;
            Sequence *c = NULL;

            if(pfamSeq) {
                for(int i  = 0; i < NUM_STRANDS; i++) {
                    TStrand strand = (TStrand)i;
                    double *seq = pfamSeq->getSeq(strand);
                    int pos = 1;
                    // get the start of the exon closest to the five' end of the sequence

                    for( ; pos <= pfamSeq->length(); pos++)
                        if(seq[pos])
                            break;
                    // get the three' end of the exon
                    if(pos > pfamSeq->length())
                        continue;

                    for(pos = pos + NUM_PHASES; pos <= pfamSeq->length(); pos += NUM_PHASES) {
                        if(seq[pos - NUM_PHASES] && !seq[pos]) {
                            break;
                        }
                    }
                    if(pos > pfamSeq->length())
                        continue;

                    int begin = 1;
                    int end = pos - (NUM_PHASES + 1);
                    if(strand == STRAND_COMP) {
                        begin = superC.length() - end + 1;
                        end = superC.length();
                    }

                    //          cerr << begin << " " << end << endl;
                    c = (Sequence *)superC.getSubSequence(begin, end);


                    if(strand == STRAND_FWD)
                        c->reverseComplement(fsm, NUM_PHASES);

                    //          cerr << ">" << c->length() << endl << c->getSeq(STRAND_COMP) << endl;

                    TParseNode node = bwfsm.node("LAST_EXON_B")->id();
                    bwlattice.setInitScore(100, node, 0);
                    node = bwfsm.node("INTERNAL_EXON_B")->id();
                    bwlattice.setInitScore(100, node, 0);
                    node = bwfsm.node("INIT_EXON_B")->id();
                    bwlattice.setInitScore(100, node, 0);
                    node = bwfsm.node("PEPINIT_EXON_B")->id();
                    bwlattice.setInitScore(100, node, 0);
                    node = bwfsm.node("SINGLE_EXON_B")->id();
                    bwlattice.setInitScore(100, node, 0);
                    node = bwfsm.node("PEPSINGLE_EXON_B")->id();
                    bwlattice.setInitScore(100, node, 0);

                    bwcore.predictTags(*c, PRED_SET);

                    if(strand == STRAND_FWD)
                        c->reverseComplement(fsm, NUM_PHASES, PRED_SET);

                    node = INVALID_NODE;
                    superC.appendTags(node, fsm, *c, PRED_SET,
                                      begin - 1, INT_MAX);

                    delete c;
                }
            }
            else
                core.predictTags(superC, PRED_SET);

            printer.displayTags(format, (std::ofstream &)cout, superC,
                                fsm, "PRED", PRED_SET);

        }

        paramStream->close();

        if(pfamSeqFile)
            delete pfamSeqFile;
        delete [] signals;

    } catch(exception *e) {
        perror(e->what());

        GenLibExcept *gle = (GenLibExcept *)e;
        if(gle->error() == BAD_USAGE)
            cerr << "use --help for more information\n";

        delete e;
        exit(1);
    }
}
