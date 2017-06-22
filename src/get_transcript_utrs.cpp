/*****************************************************************************
 * This program reports alternative splicing events from rna-seq coverage
 * and junction data. The type of events that can be identified are intron
 * retention and exon-skipping.
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


#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include "Utils.h"
#include "ContextIMM.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "Sequence.h"
#include "SequenceUtils.h"
#include "ResourceEngine.h"
#include "Configuration.h"
#include "FeatureEngine.h"
#include "StructureCore.h"
#include "GeneUtils.h"
#include "Lattice.h"
#include "GeneEvaluator.h"
#include "GeneTagPrinter.h"
#include "TagUtils.h"
#include "InpFile.h"
#include "RSqUtils.h"
#include "FL_EvdEdgeAligner.h"
#include "FL_ScoreFile.h"
#include <boost/regex.hpp>
#include "ArgParseUtils.h"
#include <numeric>
#include <functional>

#define NUM_PHASES 3

using namespace craig;

void printHelp(const char * pname, ::ofstream &fd) {
    fd << "CRAIG v. " << CRAIG_VERSION <<  "tool to estimate the UTRs of transcripts based on RNA-Seq evidence.\nThe tool uses some heuristics whose objective is to reduce the number of false positives.\nThe output can be used as training set for craigTrain\n\n";
    fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
    fd << "  usage :" << pname << " [options] FASTA_FILE ANNOT_FILE\n\n";
    ArgParseUtils::displayOption(fd, "FASTA_FILE", "Name of the file containing the input query sequence(s) in fasta format");
    ArgParseUtils::displayOption(fd,  "ANNOT_FILE", "Name of the file containing genes in locs format");
    fd << "optional arguments:\n";
    ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
    ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
    ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
    ArgParseUtils::displayOption(fd, "-c --conservative", "Conservatively filter out all transcripts that might have an ambiguous UTR prediction. Typically used for producing training data");
    ArgParseUtils::displayOption(fd, "-a --filter-all", "Filters out all transcripts corresponding to a contig if any transcript belonging to that contig is filtered out. Typically used for producing training data");
    ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "Prefix for associated xsignal files. According to naming convention, PREFIX_EVIDENCE.xsignal and PREFIX_EVIDENCE.xsigscores are the files containing the TSS and PSS signals and scores, respectively");
    ArgParseUtils::displayOption(fd, "-s --stranded", "RNA-Seq data is stranded");
    fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";
}

bool verbose = false;

int main(int argc, char *argv[]) {
    ::string fastaFileName(""), annotFileName("");
    ::string prefix_evidence = "";
    int i, strand = BOTH_STRANDS;
    std::string topology("partial");
    bool stranded = false, conservative = false, filter_all = false;

    try {
        if(!getenv("CRAIG_HOME"))
            throw EXCEPTION(NOT_ANNOTATED,
                            "CRAIG_HOME must be initialized!. See README for details");

        for(i = 1; i < argc; i++) {
            if(!strncmp(argv[i], "--verbose", 9)
               || !strncmp(argv[i], "-v", 2)) {
                verbose = true;
            }
            else if(!strncmp(argv[i], "--prefix-evidence=", 18)) {
                prefix_evidence = string(argv[i] + 18);
            }
            else if(!strncmp(argv[i], "--help", 6)
                    || !strncmp(argv[i], "-h", 2)) {
                printHelp("get_transcripts_utrs", (std::ofstream &)cout);
                exit(0);
            }
            else if(!strncmp(argv[i], "--version", 9)) {
                PRINT_VERSION(cerr, "get_transcripts_utrs", "tool for reporting AS events in RNA-seq data");
                PRINT_DISCLAIMER(cerr, "get_transcripts_utrs");
                exit(0);
            }
            else if(!strncmp(argv[i], "--stranded", 10)
                    || !strncmp(argv[i], "-s", 2)) {
                stranded = true;
            }
            else if(!strncmp(argv[i], "--conservative", 14)
                    || !strncmp(argv[i], "-c", 2)) {
                conservative = true;
            }
            else if(!strncmp(argv[i], "--filter-all", 12)
                    || !strncmp(argv[i], "-a", 2)) {
                filter_all = true;
            }
            else
                break;
        }

        if(argc - i < 1)
            throw EXCEPTION(BAD_USAGE, "insufficient arguments");

        fastaFileName = std::string(argv[i++]);
        annotFileName = std::string(argv[i++]);

        if(!prefix_evidence.length())
            prefix_evidence = fastaFileName;

        DNASigma sigma;
        Sigma evsigma(7, "-STDAXZ", "-STDAXZ", "- S T D A X Z");
        FilterEngine fe;

        /* Creating filters */
        InpFile junctionsFile("default", prefix_evidence+".junction.signal", EDGEANNOT_FASTA, &evsigma, true, false);
        InpFile covFile("default", prefix_evidence+".cov", EXTENDED_LDSCORE, NULL, true, false);

        FL_ScoreFile<ScoreSeq<double>, ScoreSeq<double>, double, double> coverage(0,  "RNASeq-Cov", &covFile, 1, 0, FT_DOUBLE);
        FL_EvdEdgeAligner junctions(1, "RNASeq-Edges", false);
        FL_FastaxFilter<EdgeAnnotSeq, EdgeAnnotSeq, EvdEdges, vector<int> > inp4junctions(2, "FastaxRNASeq-Edges", &junctionsFile, &junctions, 0);

        fe.setFilter(coverage.getName(), &coverage, 0, 1, 0, false, false, false);
        fe.setFilter(junctions.getName(), &junctions, 1, 1, 0, true, false, false);
        fe.setFilter(inp4junctions.getName(), &inp4junctions, 2, 1, 0, false, false, false);
        fe.allocateFilterValArrays();

        /*
         * loading annotSeqs and gene annotations.
         * It's assumed at training time that each Sequence object contains only
         * one gene.
         * One needs to preprocess the input data to fulfil this requirement.
         * The perl utility splitGenes.pl is provided along the C++ sources in
         * the craig distribution for this purpose.
         */

        InpFile xsignalsFile("default", prefix_evidence+".xsignal", XFASTA, &evsigma, true, false);
        InpFile xsigscoresFile("default", prefix_evidence+".xsigscores", EXTENDED_MULTI_INTSCORE, NULL, true, false);
        InpFile fastaFile("default", fastaFileName, FASTA, &sigma, true);
        list<Sequence *> & fastaSequences = (list<Sequence *> &)fastaFile.sequences();

        list<Gene> geneSet;
        GeneUtils::loadGenes(geneSet, annotFileName.c_str(), false, false,
                             fastaSequences, &sigma);

        list<Sequence *>::iterator cit = fastaSequences.begin();

        cerr << "Processing sequences" << endl;

        std::hash<std::string> hash_fn;

        for(; cit != fastaSequences.end(); cit++) {
            Sequence &c = *(*cit);
            srand ( hash_fn(c.id()) );//unsigned ( time (NULL) ) );
            cout << "CONTIG " <<  c.id() << "\t" << c.length() << endl;

            list<Gene *> &listGenes = (list<Gene *> &)c.getAnnotBioFeats();

            fe.setSequence(c);
            for(int s = 0; s < NUM_STRANDS; s++)
                fe.computeSeqFilters((TStrand)s);

            boost::RegEx rExscont("^(\\S+)\\.([^\\.]+)$");
            bool match = rExscont.Match(c.id());
            int offset = 1;

            if(match)
                if(!sscanf(rExscont[2].c_str(), "%d", &offset))
                    assert(0);

            int numStrands = stranded ? NUM_STRANDS : 1;
            vector<float> *raw_cov = new vector<float> [numStrands];
            vector<float> *contig_cov = new vector<float> [numStrands];
            vector<RSqInterval> introns;

            RSqUtils::computeSeqCoverage(c, &coverage,
                                         raw_cov,
                                         stranded);

            RSqUtils::computeJunctionCoverage(c, &junctions,
                                              contig_cov,
                                              stranded);

            RSqUtils::getRSqIntrons(c, &junctions,
                                    true, introns, BOTH_STRANDS);

            vector<double> intron_mean(numStrands), exon_mean(numStrands);
            vector<RSqChangePoint> *pchanges = new vector<RSqChangePoint> [numStrands];

            // find ChangePoints
            MultiScoreSeq<int> *xsigscore_c = (MultiScoreSeq<int> *)xsigscoresFile.findSeq(c.id());
            Sequence *xsignal_c = (Sequence *)xsignalsFile.findSeq(c.id());
            if(!xsigscore_c || !xsignal_c) {
                assert(0);
                throw EXCEPTION(CONTIG_UNAVAILABLE, c.id());
            }

            // for(int s = 0; s < numStrands; s++)  {
            //  cout << "xsignal strand " << s << "\n";
            //  for(i = 0; i < c.length(); i++) {
            //    char symbol = (*xsignal_c)((TStrand)s, i);
            //    if(symbol != '-')
            //      cout << s << "," << i << " " << symbol<< endl;
            //  }
            // }

            for(int s = 0; s < numStrands; s++)  {
                TStrand strand = (TStrand)s;
                pair<float, float> covmm = Utils::findMinMax(raw_cov[s], 0,
                                                             raw_cov[s].size());
                intron_mean[s] = (covmm.first > 1 ? log(covmm.first) : 1e-10);
                exon_mean[s] = 2*(covmm.second > 1 ? log(covmm.second) : 1e-10)/3.0;

                cerr << "intron mean " << intron_mean[s] << endl;
                cerr << "exon mean " << exon_mean[s] << endl;

                for(i = 1; i < contig_cov[s].size(); i++) {
                    contig_cov[s][i] += raw_cov[s][i];
                    float cov = contig_cov[s][i];
                    contig_cov[s][i] = (cov > 1 ? log(cov) : 1e-10);
                }

                for(i = 1; i < contig_cov[s].size(); i++)
                    raw_cov[s][i] += raw_cov[s][i - 1];

                for(i = 1; i < contig_cov[s].size(); i++) {
                    int pos = i;
                    int pos_c = (*xsigscore_c)(strand, 1, i);

                    if(!pos_c)
                        continue;

                    int pos_l = (*xsigscore_c)(strand, 0, i);
                    int pos_r = (*xsigscore_c)(strand, 2, i);

                    if(strand == STRAND_COMP) {
                        pos = c.length() - pos + 1;
                        int pos_tmp = pos_l;
                        pos_l = c.length() - pos_r + 2;
                        pos_c = c.length() - pos_c + 2;
                        pos_r = c.length() - pos_tmp + 2;
                    }

                    float mean_l = 0.1, mean_r = 0.1;
                    if(pos_l >= 1)
                        mean_l = RSqUtils::findMean(contig_cov[s], pos_l, pos_c);

                    if(pos_r <= contig_cov[s].size())
                        mean_r = RSqUtils::findMean(contig_cov[s], pos_c, pos_r);

                    if(mean_l > mean_r && pos <= c.length())
                        pos++;
                    //   cerr << pos << " " << mean_l << " " << mean_r << endl;
                    if(strand == STRAND_COMP)
                        pchanges[s].insert(pchanges[s].begin(), RSqChangePoint(pos, 100, mean_l, mean_r));
                    else
                        pchanges[s].push_back(RSqChangePoint(pos, 100, mean_l, mean_r));
                }
            }

            list<Transcript *> transcripts;
            GeneUtils::genes2transcripts(c, DEFAULT_SET, transcripts);
            int last_end = 1;
            int next_begin = -1;
            list<Transcript *>::iterator it = transcripts.begin(), nit = it;
            ostringstream ostr;
            bool validated = true;

            for( ; nit++, it != transcripts.end(); it++) {
                Transcript *transcript = (*it);
                TStrand gstrand = (stranded ? transcript->getStrand() : STRAND_FWD);
                int tLen = transcript->end() - transcript->begin() + 1;
                int tcod_beg = transcript->codingBegin();
                int tcod_end = transcript->codingEnd();
                bool gene_dws = (nit == transcripts.end() ? false : true);
                next_begin = (gene_dws ?
                              (*nit)->codingBegin() - 1:
                              contig_cov[gstrand].size());

                // for(i = tcod_beg - 1; i >= last_end; i--)
                //   if(contig_cov[gstrand][i] < 2) {
                //     last_end = i;
                //     break;
                //   }

                cout << "UPSTREAM of " << transcript->getId() << " " << last_end - 1 << " " << tcod_beg << " in strand " << transcript->getStrand() << endl;

                RSqChangePoint uponset;
                RSqUtils::selectOnsetPoints(pchanges[gstrand], last_end - 1, tcod_beg,
                                            contig_cov[gstrand],
                                            intron_mean[gstrand], exon_mean[gstrand],
                                            80, uponset);

                int tss, tes;

                if(!uponset.score) {
                    cerr << "couldn't find onset for " << c.id() << endl;
                    validated = false;
                    break;
                }

                RSqChangePoint *chosen_onset = &uponset;

                if(transcript->getStrand() == STRAND_FWD)  {
                    if(chosen_onset->pos == tcod_beg)
                        validated = false;
                    tss = chosen_onset->pos;
                    if((*xsignal_c)(STRAND_FWD, tss - 1) != 'X') {
                        cerr << " tss " << tss << " fwd is not X\n";
                        validated = false;
                    }
                }
                else  {
                    if(chosen_onset->pos == tcod_end)
                        validated = false;
                    tes = chosen_onset->pos;
                    if((*xsignal_c)(STRAND_COMP, c.length() - tes) != 'Z') {
                        cerr << " tes " << tes << " comp is not Z\n";
                        validated = false;
                    }
                }

                if(conservative) {
                    //if(!RSqUtils::isRealOnset(chosen_onset, pchanges[gstrand],
                    //                          last_end - 1, tcod_beg,
                    //                          intron_mean[gstrand],
                    //                          exon_mean[gstrand], 1)) {
                    if(chosen_onset->pos <= last_end) {
                        validated = false;
                        cerr << "no real onset " << chosen_onset->left_mean << " " << chosen_onset->right_mean << " " << chosen_onset->score << " " << validated << endl;
                    }
                }

                if(verbose)
                    RSqUtils::changePoints2Histogram(pchanges[gstrand], last_end - 1, tcod_beg, chosen_onset->pos);

                ostr << ">" << c.id() << "." << chosen_onset->pos
                     << "." << transcript->getStrand();
                ostr << "\t<" << c.id() << "_";

                // for(i = tcod_end; i <= next_begin; i++)
                //   if(contig_cov[gstrand][i] < 2) {
                //     next_begin = i;
                //     break;
                //   }

                cout << "DOWNSTREAM of " << transcript->getId() << " " << tcod_end << " " << next_begin << endl;

                RSqChangePoint dwoffset;
                RSqUtils::selectOffsetPoints(pchanges[gstrand], tcod_end, next_begin,
                                             contig_cov[gstrand],
                                             exon_mean[gstrand], intron_mean[gstrand],
                                             80, dwoffset);

                if(!dwoffset.score) {
                    cerr << "couldn't find offset for " << c.id() << endl;
                    validated = false;
                    break;
                }

                RSqChangePoint *chosen_offset = &dwoffset;
                last_end = transcript->end() + 1;

                if(transcript->getStrand() == STRAND_FWD)  {
                    if(chosen_offset->pos == tcod_beg)
                        validated = false;
                    tes = chosen_offset->pos - 1;
                    if((*xsignal_c)(STRAND_FWD, tes - 1) != 'Z') {
                        cerr << " tes " << tes << " fwd is not Z\n";
                        validated = false;
                    }
                }
                else  {
                    if(chosen_offset->pos == tcod_end)
                        validated = false;
                    tss = chosen_offset->pos - 1;
                    if((*xsignal_c)(STRAND_COMP, c.length() - tss) != 'X') {
                        cerr << " tss " << tss << " comp is not X\n";
                        validated = false;
                    }
                }

                if(conservative) {
                    //	  if(!RSqUtils::isRealOffset(chosen_offset, pchanges[gstrand],
                    //	  			     tcod_end, next_begin,
                    //	  			     exon_mean[gstrand],
                    //	  			     intron_mean[gstrand], 1, gene_dws)) {
                    if(chosen_offset->pos >= next_begin) {
                        validated = false;
                        cerr << "no real offset " << chosen_offset->left_mean << " " << chosen_offset->right_mean << " " << chosen_offset->score << " " << validated << endl;
                    }
                }

                //	cerr << "change point is found by substracting one from the given value\n";

                if(verbose)
                    RSqUtils::changePoints2Histogram(pchanges[gstrand], tcod_end, next_begin, chosen_offset->pos, '+');
                last_end = chosen_offset->pos + 1;
                ostr  << tss << "_";

                vector<Exon>& exons = transcript->exons();
                float t_support = 0;
                for(i = 1; i < exons.size(); i++) {
                    int i_beg = transcript->getStrand() == STRAND_FWD ?
                        exons[i - 1].end() + 1 : exons[i].begin() - 1;
                    int i_end = transcript->getStrand() == STRAND_FWD ?
                        exons[i].begin() - 2 : exons[i - 1].end() + 2;
                    cerr << "ann intron " << i_beg << " " << i_end << endl;
                }



                int lbegin = tss, lend = -1, rbegin = -1;
                vector<float> &cov = raw_cov[gstrand];

                for(i = 0; i < introns.size(); i++) {
                    RSqInterval &intron = introns[i];
                    cerr << "intron " << intron.start << " " << intron.end << " "
                         << intron.weight << " " << intron.strand << endl;

                    if(intron.strand != transcript->getStrand())
                        continue;

                    int rend = -1;
                    RSqInterval *next_i = i < introns.size() - 1 ? &introns[i + 1] : NULL;
                    if(intron.strand == STRAND_FWD)  {
                        lend = intron.start - 1;
                        rbegin = intron.end + 2;
                        if(tss < lend && rbegin < tcod_beg) {
                            rend = next_i
                                ? (next_i->start - 1 < tcod_beg
                                   ? next_i->start - 1
                                   : exons[0].end())
                                 : exons[0].end();
                        }
                        if(tcod_end < lend && rbegin < tes) {
                            rend = next_i ? next_i->start - 1 : tes;
                            lbegin = lbegin < tcod_end
                                              ? exons[0].begin()
                                              : lbegin;
                        }
                    }
                    else {
                        lend = intron.end + 1;
                        rbegin = intron.start - 2;
                        if(tss > lend && rbegin > tcod_end) {
                            rend = next_i
                                ? (next_i->end + 1 < tcod_end
                                   ? exons[0].begin()
                                   : next_i->end + 1)
                                : exons[0].begin();
                        }
                        if(tcod_beg > lend && rbegin > tes)  {
                            rend = next_i ? next_i->end + 1 : tes;
                            lbegin = lbegin > tcod_beg
                                ? exons[0].end()
                                : lbegin;
                        }
                    }

                    if(rend < 0) {
                        cerr << "couldn't find a suitable rexon for intron " << intron.start << ", " << intron.end << endl;
                        continue;
                    }

                    bool v;
                    if(strand == STRAND_FWD) {
                        v = RSqUtils::validateIntron(lbegin, lend,
                                                     intron.start, intron.end,
                                                     rbegin, rend, intron.weight,
                                                     false, cov);
                        if(v)
                            lbegin = rbegin > tcod_end ? rbegin
                                : exons[exons.size() - 1].begin();
                    }
                    else {
                        v = RSqUtils::validateIntron(lend, lbegin,
                                                     intron.start, intron.end,
                                                     rend, rbegin, intron.weight,
                                                     false, cov);
                        if(v)
                            lbegin = tcod_beg > rbegin ? exons[exons.size() - 1].end()
                                : rbegin;

                    }
                    if(v)
                        ostr << lend << ";" << c.id() << "_" << rbegin << "_";

                }

                ostr << tes << ">\n";

                if(!filter_all && validated) {
                    cout << ostr.str();
                    ostr.str("");
                }
            }

            if(filter_all && validated)
                cout << ostr.str();

            delete [] contig_cov;
            delete [] raw_cov;

            fe.deleteSeqFilters();
            fe.releaseSequence();
            junctionsFile.releaseSeqContents();
            covFile.releaseSeqContents();

        }

        cerr << "Done!\n";

    } catch(exception *e) {
        perror(e->what());

        GenLibExcept *gle = (GenLibExcept *)e;
        if(gle->error() == BAD_USAGE)
            cerr << "use --help for more information\n";
        delete e;
        exit(1);

    }
}
