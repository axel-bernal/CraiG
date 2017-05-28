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
#include "SelfLattice.h"
#include "GeneEvaluator.h"
#include "GeneTagPrinter.h"
#include "GeneUtils.h"
#include "InpFile.h"
#include "ArgParseUtils.h"
#define NUM_PHASES 3

using namespace craig;

void printHelp(const char * pname, ::ofstream &fd) {
    fd << "CRAIG v. " << CRAIG_VERSION << " discriminative gene prediction tool for eukarya.\n";
    fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
    fd << "Usage : " << pname << " [options] PARAMS_FILE FASTA_FILE > PREDICTION_FILE\n\n";
    fd << "positional arguments:\n";
    ArgParseUtils::displayOption(fd, "PARAMS_FILE", "Name of the file containing the gene model parameters. If PARAMS_FILE is not found, craig assumes the filename to be $(CRAIG_HOME)/models/PARAMS_FILE");
    ArgParseUtils::displayOption(fd, "FASTA_FILE", "Name of the file containing the input query sequence(s) in fasta format");
    fd << "optional arguments:\n";
    ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
    ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
    ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
    ArgParseUtils::displayOption(fd, "--format=FORMAT", "One of gtf|locs|tags. Specifies the output format for predictions. The  format \"locs\" is an internal representation of genes.  The format \"tags\" provides greater detail and is  usually used for reranking [gtf]");
    ArgParseUtils::displayOption(fd, "--output-file=OUTPUT_FILE", "Filename for the predicted genes[STDOUT]");
    ArgParseUtils::displayOption(fd, "--strand=STRAND", "One of forward|backward|both. Specifies the strand for predicting genes [both]");
    ArgParseUtils::displayOption(fd, "--best=K", "Outputs K-best segmentations [1]");
    ArgParseUtils::displayOption(fd, "--igenics=IG_FILE", "Reads intergenic regions in \"locs\" format");
    ArgParseUtils::displayOption(fd, "-m --masked", "Treat the input sequences as if they had been  previously masked. When this option is specified, lowercase letters MUST be used to represent the masked-out regions");
    ArgParseUtils::displayOption(fd, "-c --complete", "Allows for prediction of complete genes only. Incomplete genes are predicted by default");
    ArgParseUtils::displayOption(fd, "-utr --predict-utrs", "Predictions will contain UTR regions. If not specified and paramters model UTRs, expect worse performance, so choose carefully!");
    ArgParseUtils::displayOption(fd, "-rnc --report-noncoding-genes", "Report genes with no CDS segments.");
    ArgParseUtils::displayOption(fd, "--signal-resource=RESOURCE_NAME", "Read signal locations from resource RESOURCE_NAME. Here, signal is either START, STOP, DONOR, ACCEPTOR, TSS or PAS");
    ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "prefix for fetching files used as input resources");
    ArgParseUtils::displayOption(fd, "PREDICTIONS_FILE", "The list of genes predicted in the input sequence");
    fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";
}

bool verbose = false;

int main(int argc, char *argv[]) {
    std::string seqFileName("");
    std::string paramFile("");
    std::string outputFile("");
    std::string igenicsFile(""), prefix_evidence("");
    std::string format("gtf"), strand("both");
    int topK = 1, i, num_sigfilters = 4;

    std::string **sig_info = new std::string* [3];
    for(int i = 0; i < 3; i++)
        sig_info[i] = new std::string [6];

    sig_info[0][0] = "Start-Signal"; sig_info[0][1] = "Stop-Signal";
    sig_info[0][2] = "Donor-Signal"; sig_info[0][3] = "Acceptor-Signal";
    sig_info[0][4] = "TSS-Signal"; sig_info[0][5] =  "PAS-Signal";
    sig_info[1][0] = "S"; sig_info[1][1] = "T"; sig_info[1][2] = "D";
    sig_info[1][3] = "A"; sig_info[1][4] = "X"; sig_info[1][5] = "Z";
    sig_info[2][0] = ""; sig_info[2][1] = ""; sig_info[2][2] = "";
    sig_info[2][3] = ""; sig_info[2][4] = ""; sig_info[2][5] = "";

    std::string topology("partial");
    int origStrand = BOTH_STRANDS;
    bool masked = false, predictUTRs = false, reportNonCoding = false;
    boost::RegEx rExSigRes("^\\-\\-(\\S+)\\-resource=(\\S+)$");

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
                printHelp("craig", (std::ofstream &)cout);
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
            else if(!strncmp(argv[i], "-utr", 4) ||
                    !strncmp(argv[i], "--predict-utrs", 14)) {
                predictUTRs = true;
            }
            else if(!strncmp(argv[i], "-rnc", 4) ||
                    !strncmp(argv[i], "--report-noncoding-genes", 24)) {
                reportNonCoding = true;
            }
            else if(!strncmp(argv[i], "--strand=", 9)) {
                strand = std::string(argv[i] + 9);
            }
            else if(!strncmp(argv[i], "--igenics=", 10)) {
                igenicsFile = std::string(argv[i] + 10);
            }
            else if(!strncmp(argv[i], "--output-file=", 14)) {
                outputFile = std::string(argv[i] + 14);
            }
            else if(!strncmp(argv[i], "--complete", 10)
                    || !strncmp(argv[i], "-c", 2)) {
                topology = std::string("complete");
            }
            else if(!strncmp(argv[i], "--prefix-evidence=", 18)) {
                prefix_evidence = string(argv[i] + 18);
            }
            else if(rExSigRes.Match(argv[i])) {
                sig_info[2][TypeDefs::stringToTEdgeId2(rExSigRes[1])] = rExSigRes[2];
            }
            else
                break;
        }

        if(argc - i < 2) {
            printHelp("craig", (std::ofstream &)cout);
            exit(0);
        }

        if(predictUTRs)
            num_sigfilters = 6;

        paramFile = std::string(argv[i++]);
        seqFileName = std::string(argv[i++]);

        if(strand.compare("forward") == 0)
            origStrand = STRAND_FWD;
        else if(strand.compare("backward") == 0)
            origStrand = STRAND_COMP;
        else if(strand.compare("both") != 0)
            throw EXCEPTION(BAD_USAGE,
                            string("Unrecognized argument") + strand);

        std::ifstream *paramStream,  ifstr1(paramFile.c_str()), ifstr2;
        paramStream = &ifstr1;

        std::ofstream *outputStream = (std::ofstream *)&cout, ofstr2;
        if(outputFile.length()) {
            ofstr2.open(outputFile.c_str());
            outputStream = &ofstr2;
        }

        if(!ifstr1.is_open()) {
            if(!getenv("CRAIG_HOME"))
                throw EXCEPTION(BAD_USAGE,
                                "CRAIG_HOME must be initialized!. See README for details");

            paramFile = std::string(getenv("CRAIG_HOME")) + "/models/" + paramFile;
            ifstr2.open(paramFile.c_str());

            if(!ifstr2.is_open())
                throw EXCEPTION(FILE_UNAVAILABLE, paramFile);

            paramStream = &ifstr2;
        }

        // Retrieving resources
        map<std::string, std::string> resource_subs;
        resource_subs["PREFIX_EVIDENCE"] = prefix_evidence;
        ResourceEngine re(*paramStream, &resource_subs);
        Sigma *sigma = (Sigma *)re.getResource(std::string("dna-alpha"));
        FSM & fsm = *(FSM *)re.getResource(topology);
        fsm.setParseStrand((TStrand)origStrand);

        if(!predictUTRs) {
            fsm.removeId2Node(ANY_5UTR);
            fsm.removeId2Node(ANY_3UTR);
            fsm.removeId2Node(ANY_UTR_INTRON);
        }

        GeneEvaluator evaluator(fsm);
        GeneTagPrinter printer;
        InpFile annotSeqFile("default", seqFileName, FASTA, sigma, true);
        //loading annotSeqs
        list<Sequence *> & annotSeqs = (list<Sequence *> &)annotSeqFile.sequences();
        map<string, list<BioFeature> > igenics;
        GeneUtils::loadIGenicRegions(igenics, igenicsFile);

        FilterEngine fe(re, *paramStream, sig_info, num_sigfilters);
        FeatureEngine fte(fsm, fe, *paramStream);

        GlobalVector params(fte.getFeatures(), paramStream);
        fte.setParamVector(&params);

        // Initializing special variables
        TypedFilter<EdgeInst> **signals = new TypedFilter<EdgeInst> * [NUM_EDGE_INSTS];
        for(i = 0; i < NUM_EDGE_INSTS; i++)
            signals[i] = (TypedFilter<EdgeInst> *)fe.findFilter(sig_info[0][i]);

        TypedFilter<UCHAR> *contexts = (TypedFilter<UCHAR> *)fe.findFilter("GC-Content");
        FT_HiddenSeq *vftHidden = (FT_HiddenSeq *)fte.findFeature("Hidden-Stop-Sequence");

        Lattice lattice(fe, fsm, NUM_PHASES, evaluator, MAX_NUM_EXONS,
                        INTERGENIC, signals, contexts, vftHidden);
        //    SelfLattice lattice(fe, fsm, NUM_PHASES, evaluator, MAX_NUM_EXONS,
        //                        INTERGENIC, signals, contexts);
        // turn OFF masking feature if input is not masked
        if(!masked) {
            Feature *maskFeature = fte.findFeature("Test-Is-Masked");
            if(maskFeature)
                maskFeature->turnOff();
        }

// Creating StructureCore object
        StructureCore core(fsm, lattice, re, fe, fte,
                           evaluator, printer,
                           AVG_ALL, COMB_EUCLID, (TStrand)origStrand, topK);

        cerr << "\nProcessing Fasta Sequences[PID=" << getpid() << "]...\n";
        printer.displayHeader(format, *outputStream);

        if(reportNonCoding)
            printer.reportNonCodingGenes();

        /*    int loStrand = 0, upStrand = NUM_STRANDS;
              if(origStrand != BOTH_STRANDS) {
              loStrand = origStrand;
              upStrand = loStrand + 1;
              } */

        list<Sequence *>::iterator cit = annotSeqs.begin();
        for( ; cit != annotSeqs.end(); cit++) {
            Sequence &c = *(*cit);
            //      cerr << ">" << c.id() << " " << c.length() << endl;

            map<string, list<BioFeature> >::iterator igit = igenics.find(c.id());
            if(igit != igenics.end()) {
                list<BioFeature> genics;
                GeneUtils::getGenicRegions(igit->second, genics);
                core.predictTags(c, PRED_SET, &genics);
            }
            else
                core.predictTags(c, PRED_SET);

            printer.displayTags(format, *outputStream, c,
                                fsm, "", PRED_SET);
        }

        paramStream->close();
        // releasing tags
        SeqTags::releaseMem();
        delete [] signals;

        for(int i = 0; i < 3; i++)
            delete [] sig_info[i];
        delete [] sig_info;

        cerr << "Done\n";

    } catch(exception *e) {
        perror(e->what());

        GenLibExcept *gle = (GenLibExcept *)e;
        if(gle->error() == BAD_USAGE)
            cerr << "use --help for more information\n";

        delete e;
        exit(1);
    }
}
