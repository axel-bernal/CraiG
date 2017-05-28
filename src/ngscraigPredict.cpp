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

#define NUM_PHASES 3

using namespace craig;

void printHelp(const char * pname, ::ofstream &fd) {
    fd << "CRAIG v. " << CRAIG_VERSION << " prediction tool for lternatively spliced genes in eukarya.\n";
    fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
    fd << "  usage : " << pname << " [options] PARAMS_FILE FASTA_FILE > PREDICTIONS\n\n";
    fd << "positional arguments:\n";
    ArgParseUtils::displayOption(fd, "PARAMS_FILE", "Name of the file containing the gene model parameters. If PARAMS_FILE is not found, craig assumes the filename to be $(CRAIG_HOME)/models/PARAMS_FILE");
    ArgParseUtils::displayOption(fd, "FASTA_FILE", "Name of the file containing the input query sequence(s) in fasta format");
    fd << "optional arguments:\n";
    ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
    ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
    ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
    ArgParseUtils::displayOption(fd, "--format=FORMAT", "One of gtf|locs|tags. Specifies the output format for predictions. The  format \"locs\" is an internal representation of genes.  The format \"tags\" provides greater detail and is  usually used for reranking [gtf]");
    ArgParseUtils::displayOption(fd, "--strand=STRAND", "One of forward|backward|both. Specifies the strand for predicting genes [both]");
    ArgParseUtils::displayOption(fd, "-m --masked", "Treat the input sequences as if they had been  previously masked. When this option is specified, lowercase letters MUST be used to represent the masked-out regions");
    ArgParseUtils::displayOption(fd, "-c --complete", "Allows for prediction of complete genes only. Incomplete genes are predicted by default");
    ArgParseUtils::displayOption(fd, "PREDICTIONS_FILE", "The list of genes predicted in the input sequence");
    fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";
}

bool verbose = false;

int main(int argc, char *argv[]) {
    std::string seqFileName("");
    std::string paramFile("");
    std::string format("gtf"), strand("both");
    int i;
    std::string topology("partial");
    int origStrand = BOTH_STRANDS;
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
                printHelp("craig", (std::ofstream &)cout);
                exit(0);
            }
            else if(!strncmp(argv[i], "--version", 9)) {
                PRINT_VERSION(cerr, "craig", "discriminative gene prediction tool for eukarya");
                PRINT_DISCLAIMER(cerr, "craig");
                exit(0);
            }
            else if(!strncmp(argv[i], "--masked", 8)
                    || !strncmp(argv[i], "-m", 2)) {
                masked = true;
            }
            else if(!strncmp(argv[i], "--strand=", 9)) {
                strand = std::string(argv[i] + 9);
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

        if(strand.compare("forward") == 0)
            origStrand = STRAND_FWD;
        else if(strand.compare("backward") == 0)
            origStrand = STRAND_COMP;
        else if(strand.compare("both") != 0)
            throw EXCEPTION(BAD_USAGE,
                            string("Unrecognized c.l. argument") + strand);

        std::ifstream *paramStream,  ifstr1(paramFile.c_str()), ifstr2;
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
        Sigma *sigma = (Sigma *)re.getResource(std::string("dna-alpha"));
        FSM & fsm = *(FSM *)re.getResource(topology);
        fsm.setParseStrand((TStrand)origStrand);

        GeneEvaluator evaluator(fsm);
        GeneTagPrinter printer;
        InpFile annotSeqFile("default", seqFileName, FASTA, sigma, true);
        //loading annotSeqs
        list<Sequence *> & annotSeqs = (list<Sequence *> &)annotSeqFile.sequences();

        FilterEngine fe(re, *paramStream);
        FeatureEngine fte(fsm, fe, *paramStream);

        GlobalVector params(fte.getFeatures(), paramStream);
        fte.setParamVector(&params);

        // Initializing special variables

        TypedFilter<EdgeInst> **signals = new TypedFilter<EdgeInst> * [NUM_EDGE_INSTS];
        signals[START] = (TypedFilter<EdgeInst> *)fe.getFilter("Start-Signal");
        signals[STOP] = (TypedFilter<EdgeInst> *)fe.getFilter("Stop-Signal");
        signals[DONOR] = (TypedFilter<EdgeInst> *)fe.getFilter("Donor-Signal");
        signals[ACCEPTOR] = (TypedFilter<EdgeInst> *)fe.getFilter("Acceptor-Signal");
        signals[PAS] = (TypedFilter<EdgeInst> *)fe.findFilter("PAS-Signal");
        signals[TSS] = (TypedFilter<EdgeInst> *)fe.findFilter("TSS-Signal");
        TypedFilter<UCHAR> *contexts = (TypedFilter<UCHAR> *)fe.findFilter("GC-Content");
        FT_HiddenSeq *vftHidden = (FT_HiddenSeq *)fte.findFeature("Hidden-Stop-Sequence");

        Lattice lattice(fe, fsm, NUM_PHASES, evaluator, MAX_NUM_EXONS,
                        INTERGENIC, signals, contexts, vftHidden);

        if(!masked) {
            Feature *maskFeature = fte.findFeature("Test-Is-Masked");
            if(maskFeature)
                maskFeature->turnOff();
        }

        // Creating StructureCore object

        StructureCore core(fsm, lattice, re, fe, fte,
                           evaluator, printer,
                           AVG_ALL, COMB_EUCLID, (TStrand)origStrand, 1);

        cerr << "\nProcessing Fasta Sequences[PID=" << getpid() << "]...\n";
        printer.displayHeader(format, (std::ofstream &)cout);

        list<Sequence *>::iterator cit = annotSeqs.begin();

        InpFile *cov_re = (InpFile *)re.getResource("day3-coverage");
        InpFile *junc_re = (InpFile *)re.getResource("RNAseq-day3-ev-edges");
        InpFile *sig_re = (InpFile *)re.getResource("RNAseq-day3-ev-signals");
        string bfid("BF");

        for( ; cit != annotSeqs.end(); cit++) {
            Sequence &c = *(*cit);
            cerr << "processing " <<  c.id() << endl;

            // making first pass to determine intergenic regions, i.e. more than 1 gene
            list<BioFeature> igenics;
            core.predictTags(c, EXTRA1_SET);
            SeqTags &ref_pred = c.getTags(EXTRA1_SET)[0];
            SeqTags::iterator it = ref_pred.begin();
            for(; it != ref_pred.end(); it++) {
                if((*it)->getGEClass() != NODE_INST)
                    continue;

                NodeInst *b = (NodeInst *)(*it);
                Node *node = fsm.node(b->getParseType());
                int beg = b->getPos();
                int end = b->getPos() + b->getLen() - 1;
                if(b->getType() == INTERGENIC && beg != 1 && end != c.length()) {
                    igenics.push_back(BioFeature(bfid, &c, true, STRAND_FWD));
                    igenics.back().setBegin(beg);
                    igenics.back().setEnd(end);
                }
            }

            if(igenics.size()) {
                igenics.push_front(BioFeature(bfid, &c, true, STRAND_FWD));
                igenics.front().setBegin(0);
                igenics.front().setEnd(0);
                igenics.push_back(BioFeature(bfid, &c, true, STRAND_FWD));
                igenics.back().setBegin(c.length() + 1);
                igenics.back().setEnd(c.length() + 1);
            }

            list<BioFeature> allgenics;

            if(igenics.size())
                GeneUtils::getGenicRegions(igenics, allgenics);
            else {
                allgenics.push_back(BioFeature(bfid, &c, true, STRAND_FWD));
                allgenics.back().setBegin(1); allgenics.back().setEnd(c.length());
            }

            vector<vector<SeqTags> > seq_preds(allgenics.size());
            list<BioFeature>::iterator bit = allgenics.begin();
            for( ; bit != allgenics.end(); bit++) {
                vector<SeqTags> *gene_preds = &seq_preds[std::distance(allgenics.begin(), bit)];

                for(int rounds = 0; ; rounds++) {
                    c.resetTags(PRED_SET);
                    core.predictTags(c, PRED_SET);

                    EdgeAnnotSeq *junctions = (EdgeAnnotSeq *)junc_re->findSeq(c.id());
                    Sequence *signals = (Sequence *)sig_re->findSeq(c.id());
                    ScoreSeq<double> *coverage = (ScoreSeq<double> *)cov_re->findSeq(c.id());

                    if(bit == allgenics.begin() && !rounds) {
                        cov_re->lockContents();
                        junc_re->lockContents();
                        sig_re->lockContents();
                    }

                    SeqTags &last_pred = c.getTags(PRED_SET)[0];

                    // computing update
                    double cov_upd = DOUBLE_INFINITY;
                    for(it = last_pred.begin(); it != last_pred.end() ; it++) {
                        if((*it)->getGEClass() != NODE_INST)
                            continue;

                        NodeInst b = *(NodeInst *)(*it);
                        Node *inode = fsm.node(b.getParseType());

                        if(b.getPos() < bit->begin())
                            continue;
                        if(b.getPos() + b.getLen() - 1 > bit->end())
                            break;
                        if(b.getStrand() == STRAND_COMP)
                            b.setPos(c.length() - b.getPos() - b.getLen() + 2);

                        if(b.getType() == INTRON) {
                            pair<int, int> e = TagUtils::tagCoords(inode, &b);
                            pair<int, int> v = junctions->findEdge(e, b.getStrand());

                            vector<int> *jseq = junctions->getSeq(b.getStrand());

                            if(v.second >= 0 && jseq[e.second][v.second + 3]/2 < cov_upd)
                                cov_upd = jseq[e.second][v.second + 3]/2;
                        }
                    }

                    int not_juncreads = 0;
                    if(cov_upd != DOUBLE_INFINITY) {   // most likely no introns
                        // update coverage and junctions
                        for(it = last_pred.begin(); it != last_pred.end(); it++) {
                            if((*it)->getGEClass() != NODE_INST)
                                continue;

                            NodeInst b = *(NodeInst *)(*it);
                            Node *node = fsm.node(b.getParseType());

                            if(b.getPos() < bit->begin())
                                continue;
                            if(b.getPos() + b.getLen() - 1 > bit->end())
                                break;

                            if(b.getStrand() == STRAND_COMP)
                                b.setPos(c.length() - b.getPos() - b.getLen() + 2);

                            int beg = b.getPos();
                            int end = beg + b.getLen();

                            if(b.getType() == INTRON) { // update junctions
                                pair<int, int> e = TagUtils::tagCoords(node, &b);
                                pair<int, int> v = junctions->findEdge(e, b.getStrand());
                                vector<int> *jseq = junctions->getSeq(b.getStrand());

                                if(v.second < 0) {
                                    not_juncreads++;
                                    continue;
                                }

                                jseq[e.second][v.second + 3] -= cov_upd;
                                jseq[e.first][v.first + 3] -= cov_upd;

                                if(jseq[e.second][v.second + 3] < 0 ||
                                   jseq[e.first][v.first + 3] < 0)
                                    throw EXCEPTION(BAD_USAGE, "junct spl updates are wrong");

                                // check whether edges have still read support, if not remove them
                                if(jseq[e.second][v.second + 3] != 0)
                                    continue;

                                // remove first the edges
                                vector<int>::difference_type d1 = v.second + 1;
                                vector<int>::difference_type d2 = v.second + 4;
                                vector<int>::iterator it1 = jseq[e.second].begin() + d1;
                                vector<int>::iterator it2 = jseq[e.second].begin() + d2;
                                jseq[e.second].erase(it1, it2);
                                d1 = v.first + 1;
                                d2 = v.first + 4;
                                it1 = jseq[e.first].begin() + d1;
                                it2 = jseq[e.first].begin() + d2;
                                jseq[e.first].erase(it1, it2);

                                // now remove the signals
                                if(jseq[e.second].size() == 1)
                                    (*signals)(b.getStrand(), e.second - 1) = '-';

                                if(jseq[e.first].size() == 1)
                                    (*signals)(b.getStrand(), e.first - 1) = '-';
                            }
                            else if(b.getType() == EXON || b.getType() == UTR) {
                                // type if exon or utr, update coverage
                                for(i = beg; i < end; i++) {
                                    double cov = ((ScoreSeq<double> *)coverage)->getSeq(b.getStrand())[i] - cov_upd;
                                    cov = (cov < 0 ? 0 : cov);
                                    ((ScoreSeq<double> *)coverage)->getSeq(b.getStrand())[i] = cov;
                                }
                            }
                        }
                    }

                    //	cerr << cov_upd << " " << not_juncreads << " " << endl;

                    // decide to break if more than one junction was not used
                    if(not_juncreads > 3)
                        break;

                    //check whether prediction is already in c.getTags()
                    bool already_inserted = false;
                    for(i = 0; i < gene_preds->size(); i++) {
                        if((*gene_preds)[i].equals(last_pred, INTERGENIC)) {
                            //	    cerr << "prediction already made\n";
                            already_inserted = true;
                            break;
                        }
                    }

                    if(!already_inserted)
                        gene_preds->push_back(last_pred);
                    else break;
                }
            }

            cov_re->unlockContents();
            junc_re->unlockContents();
            sig_re->unlockContents();

            for(i = 0; i < seq_preds.size(); i++) {
                list<Gene> genes;
                std::ostringstream ostr;
                ostr << "PRED" << "_" << i;
                GeneUtils::tags2Genes(c, seq_preds[i], ostr.str().c_str(), fsm, genes, 0);
                printer.displayGenes(format, (std::ofstream &)cout, genes, 0, INT_MAX);
            }
        }

        paramStream->close();
        // releasing tags
        SeqTags::releaseMem();
        delete [] signals;
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
