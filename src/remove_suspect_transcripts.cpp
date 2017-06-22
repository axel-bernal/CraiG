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
#include <algorithm>
#include <numeric>
#include <functional>

#define NUM_PHASES 3

using namespace craig;

void printHelp(const char * pname, ::ofstream &fd) {
    fd << "CRAIG v. " << CRAIG_VERSION <<  "tool to remove any transcripts from a prediction based on RNA-Seq predictions\ncomputed using a naive model that are given as parameter. The tool extracts\ndifferences between the two annotations and assesses which one is the more\nprobable one based on coverage depth of introns and exons\n\n";
    fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
    fd << "  usage :" << pname << " [options] FASTA_FILE ANNOT_FILE PRED_FILE > OUTPUT_GENES\n\n";
    ArgParseUtils::displayOption(fd, "FASTA_FILE", "Name of the file containing the input query sequence(s) in fasta format");
    ArgParseUtils::displayOption(fd,  "PRED_FILE", "Name of the file containing RNA-Seq-based naive gene predictions in locs format");
    ArgParseUtils::displayOption(fd,  "ANNOT_FILE", "Name of the file containing the annotated genes in locs format");
    fd << "optional arguments:\n";
    ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
    ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
    ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
    ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "Prefix for associated xsignal files. According to naming convention, PREFIX_EVIDENCE.xsignal and PREFIX_EVIDENCE.xsigscores are the files containing the TSS and PSS signals and scores, respectively");
    ArgParseUtils::displayOption(fd, "-s --stranded", "RNA-Seq data is stranded");
    fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";
}

bool validateTranscript(list<Transcript *>::iterator last_rit,
                        list<Transcript *>::iterator rit,
                        list<Transcript *>::iterator it,
                        list<RSqInterval> &rsq_introns, bool stranded,
                        vector<float> *contig_cov) {

    bool validated = true;
    int i;
    Transcript *a_tr = (*it), *p_tr;
    list<Transcript *>::iterator aux_it = last_rit;
    vector<pair<int, int> > vann_introns, vpred_introns;
    std::map<pair<int, int>, Transcript*> intron2tr;

    a_tr->extractIntrons(vann_introns);
    for( ; aux_it != rit; aux_it++) {
        p_tr = (*aux_it);
        int offset = vpred_introns.size();
        p_tr->extractIntrons(vpred_introns);
        for(i = offset; i < vpred_introns.size(); i++)
            intron2tr[vpred_introns[i]] = p_tr;
    }

    set<pair<int, int> > s1, s2, only_ann, only_pred;

    for(i = 0; i < vann_introns.size(); i++)
        s1.insert(vann_introns[i]);
    for(i = 0; i < vpred_introns.size(); i++)
        s2.insert(vpred_introns[i]);

    std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::inserter(only_ann, only_ann.end()));
    std::set_difference(s2.begin(), s2.end(), s1.begin(), s1.end(),
                        std::inserter(only_pred, only_pred.end()));

    // check whether there are reasons to believe this annotation
    // is missing introns
    set<pair<int, int> >::iterator diff_it = only_pred.begin();
    for( ; diff_it != only_pred.end(); diff_it++) {
        pair<int, int> p = *diff_it;
        list<RSqInterval>::iterator rsq_it = rsq_introns.begin();

        for( ; rsq_it != rsq_introns.end(); rsq_it++) {
            if(rsq_it->start == p.first && rsq_it->end == p.second)
                break;
        }

        // test whether predicted intron has strong evidence
        float junction_support = 1;
        if(rsq_it != rsq_introns.end()) // found!
            junction_support = rsq_it->weight;

        std::map<pair<int, int>, Transcript*>::iterator hit  = intron2tr.find(p);
        if(hit == intron2tr.end())  {
            cerr << "problem intron " << p.first << ", " << p.second << endl;
            throw EXCEPTION(BAD_USAGE, "couldn't find transcript corresponding to intron");
        }

        p_tr = (Transcript *)hit->second;
        cerr << "transcript " << p_tr->getId() << " for " << " intron "
             << p.first << ", " << p.second << endl;
        vector<float> &cov = contig_cov[stranded ? p_tr->getStrand() :0];
        // find the neighbooring exons
        Exon *lexon, *rexon;
        if(p_tr->getNeighboorExons4Intron(p, lexon, rexon)) {
            validated &= RSqUtils::validateIntron(lexon->begin(), lexon->end(),
                                                  p.first, p.second,
                                                  rexon->begin(), rexon->end(),
                                                  junction_support,
                                                  true, cov);
        }
        else {
            cerr << "No exons for " << p.first << " " << p.second << endl;
            validated = false;
        }
    }

    // check whether there are reasons to believe this annotations
    // has wrong intron annotations
    diff_it = only_ann.begin();
    for( ; diff_it != only_ann.end(); diff_it++) {
        pair<int, int> p = *diff_it;
        list<RSqInterval>::iterator rsq_it = rsq_introns.begin();
        for( ; rsq_it != rsq_introns.end(); rsq_it++) {
            if(rsq_it->start == p.first && rsq_it->end == p.second)
                break;
        }

        // test whether annotated intron has strong evidence
        float junction_support = 1;
        if(rsq_it != rsq_introns.end()) // found!
            junction_support = rsq_it->weight;

        vector<float> &cov = contig_cov[stranded ? a_tr->getStrand() : 0];

        // find the neighbooring exons
        Exon *lexon, *rexon;
        if(a_tr->getNeighboorExons4Intron(p, lexon, rexon)) {
            validated &= RSqUtils::validateIntron(lexon->begin(), lexon->end(),
                                                  p.first, p.second,
                                                  rexon->begin(), rexon->end(),
                                                  junction_support,
                                                  false, cov);
        }
        else {
            cerr << "No exons for " << p.first << " " << p.second << endl;
            validated = false;
        }
    }

    return validated;
}


bool verbose = false;

int main(int argc, char *argv[]) {
    ::string fastaFileName(""), predFileName(""), annotFileName("");
    ::string prefix_evidence = "";
    int i, strand = BOTH_STRANDS;
    std::string topology("partial");
    bool stranded = false;

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
            else
                break;
        }

        if(argc - i < 1)
            throw EXCEPTION(BAD_USAGE, "insufficient arguments");

        fastaFileName = std::string(argv[i++]);
        annotFileName = std::string(argv[i++]);
        predFileName = std::string(argv[i++]);

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

        InpFile fastaFile("default", fastaFileName, FASTA, &sigma, true);
        list<Sequence *> & fastaSequences = (list<Sequence *> &)fastaFile.sequences();

        list<Gene> predSet, annotSet;
        GeneUtils::loadGenes(predSet, predFileName.c_str(), false, false,
                             fastaSequences, &sigma, PRED_SET);

        GeneUtils::loadGenes(annotSet, annotFileName.c_str(), false, false,
                             fastaSequences, &sigma, TRAIN_SET);

        list<Sequence *>::iterator cit = fastaSequences.begin();

        cerr << "Processing sequences" << endl;

        std::hash<std::string> hash_fn;

        for(; cit != fastaSequences.end(); cit++) {
            Sequence &c = *(*cit);
            srand ( hash_fn(c.id()) );//unsigned ( time (NULL) ) );
            cerr << "CONTIG " <<  c.id() << "\t" << c.length() << endl;

            fe.setSequence(c);
            for(int s = 0; s < NUM_STRANDS; s++)
                fe.computeSeqFilters((TStrand)s);

            int numStrands = stranded ? NUM_STRANDS : 1;
            vector<float> *contig_cov = new vector<float> [numStrands];
            RSqUtils::computeSeqCoverage(c, &coverage,
                                         contig_cov,
                                         stranded);

            list<RSqInterval> rsq_introns;
            RSqUtils::getRSqIntrons(c, &junctions, rsq_introns, BOTH_STRANDS);

            for(int s = 0; s < numStrands; s++)  {
                TStrand strand = (TStrand)s;
                for(i = 1; i < contig_cov[s].size(); i++) {
                    contig_cov[s][i] += contig_cov[s][i - 1];

                }
            }

            list<Transcript *> pred_transcripts, ann_transcripts;
            GeneUtils::genes2transcripts(c, PRED_SET, pred_transcripts);
            GeneUtils::genes2transcripts(c, TRAIN_SET, ann_transcripts);
            list<Transcript *>::iterator it = ann_transcripts.begin(),
                rit = pred_transcripts.begin(),
                first_olap_rit = pred_transcripts.end();

            while(it != ann_transcripts.end() && rit != pred_transcripts.end()) {
                Transcript *a_tr = (*it), *p_tr = (*rit);
                //	print "testing $gene->{'Id'} (", $gene->lowestCoord(), "-", $gene->highestCoord(),  ") and $rgene->{'Id'} (", $rgene->lowestCoord(), "-", $rgene->highestCoord(), ") $b\n";
                if(p_tr->end() < a_tr->begin() || p_tr->begin() > a_tr->end()) {
                    if(p_tr->end() < a_tr->begin())
                        rit++;
                    else {
                        if(first_olap_rit != pred_transcripts.end())
                            if(validateTranscript(first_olap_rit, rit, it,
                                                  rsq_introns, stranded, contig_cov))
                                cout << ">" << a_tr->gene()->getId() << "|"
                                     << a_tr->getId() << endl;
                        it++;
                        first_olap_rit = pred_transcripts.end();
                    }
                }
                else { // overlap found!
                    cerr << "overlap " << a_tr->getId() << " " << p_tr->getId() << endl;
                    if(first_olap_rit == pred_transcripts.end())
                        first_olap_rit = rit;
                    rit++;
                }
            }

            if(first_olap_rit != pred_transcripts.end())
                if(validateTranscript(first_olap_rit, rit, it,
                                      rsq_introns, stranded, contig_cov))
                    cout << ">" << (*it)->gene()->getId() << "|"
                         << (*it)->getId() << endl;

            delete [] contig_cov;

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
