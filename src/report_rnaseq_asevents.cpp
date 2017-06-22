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
#include "FL_EvdEdgeAligner.h"
#include "FL_ScoreFile.h"
#include <boost/regex.hpp>
#include "RSqUtils.h"
#include "ArgParseUtils.h"
#include <list>
#include <vector>

#define NUM_PHASES 3
#define BIN_SIZE 5
#define NUM_JUNCTION_TYPES 4
using namespace craig;

struct genloc {
    string cid;
    int start;
    int end;
    double abundance;
    double ratio;

    genloc(string & cid, int start, int end, double abundance, double ratio) {
        this->cid = cid;
        this->start = start;
        this->end = end;
        this->abundance = abundance;
        this->ratio = ratio;
    }

    genloc(Sequence &c, RSqInterval &ri, int offset = 0,
           double abundance = 0, double ratio = 0) {

        cid = c.id();
        start = ri.start;
        end = ri.end;
        this->abundance = abundance;
        this->ratio = ratio;
        shiftCoordinates(offset);
    }

    void shiftCoordinates(int offset) {
        start += offset;
        end += offset;
    }

    void swapCoordinates() {
        int tmp = start;
        start = end;
        end = tmp;
    }

};


// compute Is and V
/*
  if(annotated_junction)
  insert_interval(Is, start, end, weight, strand);
  // compute Ins
  for(int i = 1; i <= c.length(); i++) {
  if(!I[i].first)
  continue;
  if(!RSqUtils::interval_match(Is, I[i].first, I[i].second))
  RSqUtils::insert_interval(Ins, I[i].first, I[i].second, 0, strand);
  }

  float ratio1 = lit->weight*100/(rit->weight ? rit->weight : 1),
  ratio2 = rit->weight*100/(lit->weight ? lit->weight : 1);
  float ratio = ratio1 < ratio2 ? ratio1 : ratio2;
*/

void print_olap_entry(vector<genloc> &ptrS,
		      int type, int trate, int ratio,
		      bool in_phase_skips) {

    if(!ptrS.size())
        return;

    int num_ovregs = type == 3 ? 1 : (type == 2 ? 2 : 3);
    int num_olaps = ptrS.size()/num_ovregs;
    string prefix = type == 3 ? "nolap" : (type == 2 ? "alt3or5" : (type == 0 ? "exonskip" : "multiolap"));
    cout << "stats[" << type << "][" << trate << "][" << ratio << "] = " << num_olaps;

    for(int i = 0; i < ptrS.size(); i += num_ovregs) {
        if(type == 0 && in_phase_skips) {
            if((abs(ptrS[i].start - ptrS[i].end) + 2) % 3)
                continue;
        }
        cout << "\t" << ptrS[i].abundance << " " << ptrS[i].ratio;
    }
    cout << "\n";

    for(int i = 0; i < ptrS.size(); i += num_ovregs) {
        if(type == 0 && in_phase_skips) {
            if((abs(ptrS[i].start - ptrS[i].end) + 2) % 3)
                continue;
        }

        cout << ">" << ptrS[i].cid << "_" << prefix << "_" << i << "\t";

        for(int j = 0; j < num_ovregs; j++) {
            int idx = j ? num_ovregs - j : j;
            cout << ptrS[i + idx].cid << "_" << ptrS[i + idx].start   << "_" << ptrS[i + idx].end;
            if(num_ovregs > 1 && j < num_ovregs - 1)
                cout << ";";
        }
        cout << "\n";
    }
}

void interval_overlaps(Sequence &c,
		       vector<vector<vector<genloc> > > *S,
		       list<RSqInterval> &V, vector<pair<int, int> > &I,
		       int offset, TStrand strand) {

    list<RSqInterval>::iterator pivot = V.begin(), it = V.begin(),
        last_olap = V.end();
    vector<genloc> *ptrS;

    while(++it != V.end()) {
        if(pivot->end >= it->start) {
            //overlap between pivot and it
            if(last_olap != V.end()) {
                float w = last_olap->end < it->start ?
                    (last_olap->weight + it->weight)/2 :
                    (last_olap->weight > it->weight ? last_olap->weight : it->weight);
                float min_w = w < pivot->weight ? w : pivot->weight;
                float max_w = w < pivot->weight ? pivot->weight : w;
                float ratio = max_w ? min_w*100/max_w : 0;

                if(last_olap->end < it->start)  // exon skip
                    ptrS = &S[0][(int)log(max_w + 1)][ratio/BIN_SIZE];
                else // multiple overlap
                    ptrS = &S[1][(int)log(max_w + 1)][ratio/BIN_SIZE];

                ptrS->push_back(genloc(c, *pivot, offset - 1, max_w, ratio));
                ptrS->push_back(genloc(c, *it, offset - 1));
                ptrS->push_back(genloc(c, *last_olap, offset - 1));

                if(strand == STRAND_COMP) {
                    unsigned int sz = ptrS->size();
                    (*ptrS)[sz - 1].swapCoordinates();
                    (*ptrS)[sz - 2].swapCoordinates();
                    (*ptrS)[sz - 3].swapCoordinates();
                }
            }
            else {
                float min_w = it->weight < pivot->weight ? it->weight : pivot->weight;
                float max_w = it->weight < pivot->weight ? pivot->weight : it->weight;
                float ratio = max_w ? min_w*100/max_w : 0;
                ptrS = &S[2][(int)log(max_w + 1)][ratio/BIN_SIZE];

                ptrS->push_back(genloc(c, *pivot, offset - 1, max_w, ratio));
                ptrS->push_back(genloc(c, *it, offset - 1));

                if(strand == STRAND_COMP) {
                    unsigned int sz = ptrS->size();
                    (*ptrS)[sz - 1].swapCoordinates();
                    (*ptrS)[sz - 2].swapCoordinates();
                }
            }
            //      cout << "iter " << ptrS->size() - 1 << " " << ptrS->back()->start << " " << ptrS->back()->end << endl;

            if(pivot->end < it->end) {
                last_olap = pivot;
                pivot = it;
            }
            else
                last_olap = it;
        }
        else {
            ptrS = &S[3][0][0];
            ptrS->push_back(genloc(c, *pivot, offset - 1, pivot->weight, 0));
            if(strand == STRAND_COMP) {
                unsigned int sz = ptrS->size();
                (*ptrS)[sz - 1].swapCoordinates();
            }

            last_olap = V.end();
            pivot = it;
        }
    }
}


void printHelp(const char * pname, ::ofstream &fd) {
    fd << "CRAIG v. " << CRAIG_VERSION <<  "tool for reporting AS events and the ids of the contigs that fall into each category\n\n";
    fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
    fd << "Usage :" << pname << " [options] FASTA_FILE ANNOT_FILE > AS_REPORT\n\n";
    fd << "positional arguments:\n";
    ArgParseUtils::displayOption(fd, "FASTA_FILE", "Name of the file containing the input query sequence(s) in fasta format");
    ArgParseUtils::displayOption(fd,  "ANNOT_FILE", "Name of the file containing genes in locs format");
    fd << "optional arguments:\n";
    ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
    ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
    ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
    ArgParseUtils::displayOption(fd, "-a --abs-coords", "Report is made using absolute chromosome coordinates");
    ArgParseUtils::displayOption(fd, "-p --in-phase-only", "exon skip events are reported only when skipped exon is in phase");
    ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "Prefix for associated RNAseq files. According to naming convention, junction and small reads must be located in files PREFIX_EVIDENCE.junctions.signal and PREFIX_EVIDENCE.cov respectively. The output will be reported in PREFIX_EVIDENCE.xsignal and PREFIX_EVIDENCE.xsigscore [FASTA_FILE]");
    fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";
}


bool verbose = false;

int main(int argc, char *argv[]) {
    ::string fastaFileName(""), annotFileName("");
    ::string prefix_evidence("");
    int i, strand = BOTH_STRANDS;
    bool abs_coords = false, in_phase_skips = false;
    std::string topology("partial");
    std::string species = "human";

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
                printHelp("report_rnaseq_asevents", (std::ofstream &)cout);
                exit(0);
            }
            else if(!strncmp(argv[i], "--version", 9)) {
                PRINT_VERSION(cerr, "report_rnaseq_asevents", "tool for reporting AS events in RNA-seq data");
                PRINT_DISCLAIMER(cerr, "report_rnaseq_asevents");
                exit(0);
            }
            else if(!strncmp(argv[i], "--species=", 10)) {
                species = std::string(argv[i] + 10);
            }
            else if(!strncmp(argv[i], "-a", 2) || !strncmp(argv[i], "--abs-coords=", 13)) {
                abs_coords = true;
            }
            else if(!strncmp(argv[i], "-p", 2) || !strncmp(argv[i], "--in-phase-only=", 16)) {
                in_phase_skips = true;
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
        InpFile junctionsFile("default", prefix_evidence+".junction.signal", EDGEANNOT_FASTA, &evsigma, true);
        InpFile covFile("default", prefix_evidence+".cov", EXTENDED_LDSCORE, NULL, true);

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

        InpFile fastaFile("default", fastaFileName, FASTA, &sigma, true);
        list<Sequence *> & fastaSequences = (list<Sequence *> &)fastaFile.sequences();

        list<Gene> geneSet;
        GeneUtils::loadGenes(geneSet, annotFileName.c_str(), false, false,
                             fastaSequences, &sigma);

        cerr << "Processing sequences" << endl;
        int skip = 1, retained = 1, intronct = 1;

        vector<vector<vector<genloc> > > *S = new vector<vector<vector<genloc > > > [NUM_JUNCTION_TYPES];

        for(i = 0; i < NUM_JUNCTION_TYPES; i++) {
            S[i] = vector<vector<vector<genloc> > >(40);
            for(int j = 0; j < 40; j++)
                S[i][j] = vector<vector<genloc> >(100/BIN_SIZE + 1);
        }

        list<Sequence *>::iterator cit = fastaSequences.begin();

        for(; cit != fastaSequences.end(); cit++) {
            Sequence &c = *(*cit);
            list<Gene *> &listGenes = (list<Gene *> &)c.getAnnotBioFeats();

            int offset = 1;

            if(abs_coords) {
                boost::RegEx rExscont("^(\\S+)\\.([^\\.]+)$");
                bool match = rExscont.Match(c.id());

                if(match) {
                    if(!sscanf(rExscont[2].c_str(), "%d", &offset))
                        assert(0);
                }
                else cerr << "Absolute coordinates not available\nProceed with relative ones";
            }

            fe.setSequence(c);

            for(int s = 0; s < NUM_STRANDS; s++) {
                TStrand strand = (TStrand)s;
                fe.computeSeqFilters(strand);
                list<RSqInterval> V;
                vector<pair<int, int> > I(c.length() + 1);
                list<RSqInterval> Is;
                list<RSqInterval> Ins;

                list<Gene *>::iterator git = listGenes.begin();
                for( ; git != listGenes.end(); git++) {
                    Gene &gene = *(*git);
                    if(gene.getStrand() != strand)
                        continue;

                    vector<pair<int, int> > tIntrons;
                    for(int j = 0; j < gene.transcripts().size(); j++) {
                        Transcript &t = gene.transcripts()[j];
                        // extract all transcript introns to I in orderly manner.
                        t.extractIntrons(tIntrons);
                    }

                    vector<pair<int, int> >::iterator it = tIntrons.begin();
                    for( ; it != tIntrons.end(); it++) {
                        //	    cout << gene.getId() << " intron " << it->first << " " << it->second << " " << strand << endl;
                        RSqUtils::insert_hash_interval(*it, I);
                    }

                    /*	  double ovmean = 0;
                          double ovcounter = 0;

                          double params[] = {5,3,2};
                          bool intron_retain = false;
                          ovmean = ovcounter ? ovmean/ovcounter : 0;
                          cerr << "ovmean=" << ovmean << endl;
                          // looking for retained introns
                          list<RSqInterval>::iterator it = I.begin();
                          for( ; it != I.end(); it++) {
                          vector<float> cov(it->end - it->start + 1);
                          for(k = it->start; k <= it->end; k++)  {
                          cov[k - start] = coverage.value(k, strand);
                          //	      cout << cov[k - leftB] << " ";
                          }
                          //	    cout << endl;
                          it->computeCoverageStats(cov.begin(), cov.end());
                          it->.push_back(intron);
                          ovmean += intron.cov_mean;
                          ovcounter++;

                          RSqInterval &is = introns[i];
                          if(is.cov_mean > ovmean + params[ovmean < 4 ? 0 : (ovmean < 10 ? 1 : 2)]*is.cov_stdev) {
                          cout << ">" << c.id() << ".retained." << retained++ << "\t" << c.id() << "_" << is.start << "_" << is.end << " " << strand << endl;
                          intron_retain = true;
                          }
                          } */
                }

                RSqUtils::getRSqIntrons(c, &junctions, V, strand);

                // Compute alternative 3' or 5' consistent with annotation
                interval_overlaps(c, S, V, I, offset, strand);
            }

            fe.deleteSeqFilters();
            fe.releaseSequence();
            junctionsFile.releaseSeqContents();
            covFile.releaseSeqContents();
        }

        for(int i = 0; i < NUM_JUNCTION_TYPES; i++) {
            for(int j = 0; j < S[i].size(); j++) {
                for(int k = 0; k < S[i][j].size(); k++)
                    print_olap_entry(S[i][j][k], i, j, k, in_phase_skips);
            }
        }

        delete [] S;

        cerr << "Done!\n";

        //DONT FORGET
        //    pStream.close();

    } catch(exception *e) {
        perror(e->what());

        GenLibExcept *gle = (GenLibExcept *)e;
        if(gle->error() == BAD_USAGE)
            cerr << "use --help for more information\n";
        delete e;
        exit(1);

    }
}
