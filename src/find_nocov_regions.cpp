/*****************************************************************************
 * This program reports alternative splicing events from rna-seq coverage and
 * junction data. The type of events that can be identified are intron
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


#include <string.h>
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
#include "FL_ScoreFile.h"
#include <boost/regex.hpp>
#include <numeric>
#define NUM_PHASES 3

using namespace craig;

void printHelp(const char * pname, ::ofstream &fd) {
    fd << "CRAIG v. " << CRAIG_VERSION <<  "tool for reporting AS events and the ids of the contigs that fall into each category\n\n";
    fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
    fd << "Usage :" << pname << " [options] FASTA_FILE > NOCOV_REGIONS\n\n";
    fd << "positional arguments:\n";
    ArgParseUtils::displayOption(fd, "FASTA_FILE", "Name of the file containing the input query sequence(s) in fasta format");
    fd << "optional arguments:\n";
    ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
    ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
    ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
    ArgParseUtils::displayOption(fd, "-s --stranded", "RNA-Seq data is stranded");
    ArgParseUtils::displayOption(fd, "--coverage-density=COVERAGE_DENSITY", "The coverage density factor. A number between 1 (very sparse) to 12 (very dense) [1]");
    ArgParseUtils::displayOption(fd, "--smooth-window=SMOOTH_WINDOW", "Window size for smoothing coverage data[0]");
    ArgParseUtils::displayOption(fd, "--min-igenic=MIN_IGENIC", "Minimum length of a region with no coverage (intergenic)[200]");
    ArgParseUtils::displayOption(fd, "--max-genic=MAX_GENIC", "Maximum length of a region with coverage (genic)[250000]");
    ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "Prefix for associated RNAseq files. According to naming convention, junction and small reads must be located in files PREFIX_EVIDENCE.junctions.signal and PREFIX_EVIDENCE.cov respectively. The output will be reported in PREFIX_EVIDENCE.xsignal and PREFIX_EVIDENCE.xsigscore [FASTA_FILE]");
    ArgParseUtils::displayOption(fd, "--min-coverage=MIN_COVERAGE", "minimum coverage value for a region to be considered genic.[2]");
    fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";
}


bool verbose = false;

int main(int argc, char *argv[]) {
    ::string modelPath, fastaFileName("");
    int i, strand = BOTH_STRANDS, smooth_window = 0;
    ::string prefix = "";
    bool stranded = false;
    int min_igenic = 200;
    int min_genic = 200;
    int cdfactor = 1;
    int max_genic = 250000;
    double gen_mincov = 2;
    std::string topology("partial");

    try {
        if(!getenv("CRAIG_HOME"))
            throw EXCEPTION(NOT_ANNOTATED,
                            "CRAIG_HOME must be initialized!. See README for details");

        modelPath = std::string(getenv("CRAIG_HOME")) + "/models/";

        for(i = 1; i < argc; i++) {
            if(!strncmp(argv[i], "--verbose", 9)
               || !strncmp(argv[i], "-v", 2)) {
                verbose = true;
            }
            else if(!strncmp(argv[i], "--stranded", 10)
                    || !strncmp(argv[i], "-s", 2)) {
                stranded = true;
            }
            else if(!strncmp(argv[i], "--prefix-evidence=", 18)) {
                prefix = std::string(argv[i] + 18);
            }
            else if(!strncmp(argv[i], "--smooth-window=", 16)) {
                sscanf(argv[i] + 16, "%d", &smooth_window);
            }
            else if(!strncmp(argv[i], "--min-coverage=", 15)) {
                sscanf(argv[i] + 15, "%lf", &gen_mincov);
            }
            else if(!strncmp(argv[i], "--coverage-density=", 19)) {
                sscanf(argv[i] + 19, "%d", &cdfactor);
            }
            else if(!strncmp(argv[i], "--min-igenic=", 13)) {
                sscanf(argv[i] + 13, "%d", &min_igenic);
            }
            else if(!strncmp(argv[i], "--min-genic=", 12)) {
                sscanf(argv[i] + 12, "%d", &min_genic);
            }
            else if(!strncmp(argv[i], "--max-genic=", 12)) {
                sscanf(argv[i] + 12, "%d", &max_genic);
            }
            else if(!strncmp(argv[i], "--help", 6)
                    || !strncmp(argv[i], "-h", 2)) {
                printHelp("find_nocov_regions", (std::ofstream &)cout);
                exit(0);
            }
            else if(!strncmp(argv[i], "--version", 9)) {
                PRINT_VERSION(cerr, "find_nocov_regions", "tool for finding non covered regions  in RNA-seq data");
                PRINT_DISCLAIMER(cerr, "find_nocov_regions");
                exit(0);
            }
            else if(!strncmp(argv[i], "--stranded", 10)
                    || !strncmp(argv[i], "-s", 2)) {
                stranded = true;
            }
            else
                break;

        }

        if(argc - i < 1)
            throw EXCEPTION(INSUFFICIENT_ARGS, "insufficient arguments");

        fastaFileName = std::string(argv[i++]);

        if(!prefix.length())
            prefix = fastaFileName;

        DNASigma sigma;
        Sigma evsigma(7, "-STDAXZ", "-STDAXZ", "- S T D A X Z");
        FilterEngine fe;

        /* Creating filters */
        TypedFilter<EdgeInst> **signals = new TypedFilter<EdgeInst> * [NUM_EDGE_INSTS];
        signals[START] = new FL_Signal<EdgeInst>(0, "Start-Signal", START, "ATG", &sigma);
        signals[STOP] = new FL_Signal<EdgeInst>(1, "Stop-Signal", STOP, "TAG TGA TAA", &sigma);
        signals[DONOR] = new FL_Signal<EdgeInst>(2, "Donor-Signal", DONOR, "GT", &sigma);
        signals[ACCEPTOR] = new FL_Signal<EdgeInst>(3, "Acceptor-Signal", ACCEPTOR, "AG", &sigma);
        InpFile junctionsFile("default", prefix + ".junction.signal", EDGEANNOT_FASTA, &evsigma, true, false);
        InpFile covFile("default", prefix + ".cov", EXTENDED_LDSCORE, NULL, true, false);

        FL_ScoreFile<ScoreSeq<double>, ScoreSeq<double>, double, double> coverage(4,  "RNASeq-Cov", &covFile, 1, 0, FT_DOUBLE);
        FL_EvdEdgeAligner junctions(5, "RNASeq-Edges", false);
        FL_FastaxFilter<EdgeAnnotSeq, EdgeAnnotSeq, EvdEdges, vector<int> > inp4junctions(6, "FastaxRNASeq-Edges", &junctionsFile, &junctions, 0);

        for(i = 0; i < 4; i++)
            fe.setFilter(signals[i]->getName(), signals[i], i, 1, 0, false, false, false);

        fe.setFilter(coverage.getName(), &coverage, 4, 1, 0, false, false, false);
        fe.setFilter(junctions.getName(), &junctions, 5, 1, 0, true, false, false);
        fe.setFilter(inp4junctions.getName(), &inp4junctions, 6, 1, 0, false, false, false);

        fe.allocateFilterValArrays();

        if(cdfactor < 5) {
            min_igenic = pow(5 - cdfactor, log(5)/log(4))*400;
            min_genic = pow(2, cdfactor - 1)*20;
        }
        else if(cdfactor < 13) {
            min_igenic = 200 - (cdfactor - 5)*25;
            min_genic = 200;
        }
        else
            throw EXCEPTION(BAD_USAGE, "density must be between 1 and 12");


        /*
         * loading annotSeqs and gene annotations.
         * It's assumed at training time that each Sequence object contains only
         * one gene.
         * One needs to preprocess the input data to fulfil this requirement.
         * The perl utility splitGenes.pl is provided along the C++ sources in
         * the craig distribution for this purpose.
         */

        InpFile fastaFile("default", fastaFileName, FASTA, &sigma, true);
        list<Sequence *> & fastaSequences = (list<Sequence *> &)fastaFile.sequences();

        list<Sequence *>::iterator cit = fastaSequences.begin();

        cerr << "Processing sequences"<< endl;

        for(; cit != fastaSequences.end(); cit++) {
            Sequence &c = *(*cit);
            fe.setSequence(c);
            //      cout << "CONTIG " << c.id() << endl;
            for(int s = 0; s < NUM_STRANDS; s++)
                fe.computeSeqFilters((TStrand)s);

            boost::RegEx rExscont("^(\\S+)\\.([^\\.]+)$");
            bool match = rExscont.Match(c.id());
            int offset = 1;

            if(match)
                if(!sscanf(rExscont[2].c_str(), "%d", &offset))
                    assert(0);

            int numStrands = stranded ? NUM_STRANDS : 1;
            vector<float> *contig_cov = new vector<float> [numStrands];
            vector<float> *junction_cov = new vector<float> [numStrands];

            RSqUtils::computeSeqCoverage(c, &coverage,
                                         contig_cov,
                                         stranded);

            RSqUtils::computeJunctionCoverage(c, &junctions,
                                              junction_cov,
                                              stranded);

            vector<pair<int, int> > nocov;
            vector<float> fcontig_cov(contig_cov[0].size());

            for(int s = 0; s < numStrands; s++) {
                for(i = 1; i < junction_cov[s].size(); i++) {
                    fcontig_cov[i] += contig_cov[s][i] + junction_cov[s][i];
                    //	    fcontig_cov[i] += contig_cov[s][i] + ((junction_cov[s][i] < gen_mincov && junction_cov[s][i] >= 1) ? gen_mincov : junction_cov[s][i]);
                }
            }

            //      cout << "coverage\n";
            //      for(i = 1; i < fcontig_cov.size(); i++)
            //	cout << i << " " << fcontig_cov[i] << endl;

            if(smooth_window)
                Utils::smoothScores(fcontig_cov, 1, smooth_window, c.length());

            vector<int> splpos;
            vector<float> splcov;

            RSqUtils::computeSplSCoverage(signals, fcontig_cov,
                                          -1, 1, c.length(),
                                          splpos, splcov, gen_mincov - 1,
                                          0, TR_ANYWHERE, BOTH_STRANDS);

            //      cout << "splcoverage\n";
            //      for(i = 0; i < splpos.size(); i++)
            //      	cout << i << " " << splpos[i] << " " << splcov[i] << endl;


            int ncbeg = 0, ncend = 0;
            int tolerance = 0;
            int state = 0;
            for(i = 0; i < splcov.size(); i++) {
                if(state == 0) { // initial state
                    if(splcov[i] <= tolerance) {
                        //	    cout << " tr 0->1 @ " << splcov[i] << " " << splpos[i] << endl;
                        state = 1;
                        ncbeg = i;
                    }
                    else {
                        ncend = 0;
                        //	    cout << " tr 0->2 @ " << splcov[i] << " " << splpos[i] << endl;
                        state = 2;
                    }
                }
                else if(state == 1) { // no coverage
                    if(splcov[i] > tolerance) {
                        ncend = i - 1;
                        state = 2;
                        tolerance = 0;
                        //	    cout << " tr 1->2 @ " << splcov[i] << " " << splpos[i] << endl;
                        if(splpos[ncend] > splpos[ncbeg] + min_igenic) {
                            nocov.push_back(pair<int, int>(splpos[ncbeg], splpos[ncend]));
                            //	      cout << "pushed " << tolerance << " " << nocov.back().first << " " << nocov.back().second << endl;
                        }
                    }
                }
                else if(state == 2) { // coverage
                    if(splpos[i] > splpos[ncend + 1] + (int)(0.45*max_genic)) { // too long
                        //	    cerr << "too large " << c.id() << "_" << splpos[ncend + 1] << "_" << splpos[i] << " " << tolerance << endl;
                        i = ncend + 1;
                        tolerance = (tolerance == 0 ? 2 : tolerance + 1);
                    }
                    else if(splcov[i] <= tolerance) {
                        if(splpos[i] > splpos[ncend + 1] + min_genic) {
                            ncbeg = i;
                            //	      cerr << "genic is ok " << splpos[ncend + 1] << " " << splpos[i] << " new no cov beg " << splpos[ncbeg] << " " << tolerance << "\n";
                        }
                        else {
                            if(nocov.size() && nocov.back().first == splpos[ncbeg]) {
                                nocov.pop_back();
                                //		cerr << "poped\n";
                            }
                            ncend = i;
                            //	      cerr << "genic is not ok " << splpos[ncend + 1] << " " << splpos[i] << " new no cov end " << splpos[ncend] << " " << tolerance << "\n";
                        }
                        //	    cout << " tr 2->1 @ " << splcov[i] << " " << splpos[i] << endl;
                        state = 1;
                    }
                }
            }

            assert(state != 0);
            if(state == 1 && c.length() > splpos[ncbeg] + min_igenic)
                nocov.push_back(pair<int, int>(splpos[ncbeg], c.length()));
            else
                nocov.push_back(pair<int, int>(c.length() + 1, c.length() + 1));

            //      cout << "pushed " << state << " " << tolerance << " " << nocov.back().first << " " << nocov.back().second << endl;

            //      cout << last << " " << i << " " << ncbeg << " " << ncend << " " << tolerance << " " << splcov[last] << " " << splcov[i] << endl;
            //      cout << "pushed\n";

            for(i = 0; i < nocov.size(); i++) {
                if(i == 0 && nocov[i].first > 1) {
                    cout << ">" << c.id() << ".0\t"
                         << c.id() << "_0_0" << endl;
                }
                cout << ">" << c.id() << "." << nocov[i].first
                     << "\t" << c.id() << "_" << nocov[i].first
                     << "_" << nocov[i].second << endl;
            }

            // freeing memory
            delete [] contig_cov;
            delete [] junction_cov;

            fe.deleteSeqFilters();
            fe.releaseSequence();
            junctionsFile.releaseSeqContents();
            covFile.releaseSeqContents();
        }

        for(i = 0; i < 4; i++)
            delete signals[i];
        delete [] signals;

        cerr << "Done!\n";

    } catch(exception *e) {
        GenLibExcept *gle = (GenLibExcept *)e;
        if(gle->error() == INSUFFICIENT_ARGS)
            cerr << "use --help for more information\n";
        else
            perror(e->what());
        delete e;
        exit(1);

    }
}
