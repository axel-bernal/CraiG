/*****************************************************************************
 * This is the main of the list_unreachable_genes. It converts gene locations
 * into tags, which will be later used within the lless library. The structure
 * of ach gene is specified with the --species option
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
#include "Utils.h"
#include "ContextIMM.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "Sequence.h"
#include "SequenceUtils.h"
#include "ResourceEngine.h"
#include "Lattice.h"
#include "FeatureEngine.h"
#include "GeneUtils.h"
#include "GeneTagPrinter.h"
#include "GeneEvaluator.h"
#include "FSM.h"
#include "TagUtils.h"
#include "InpFile.h"

#define NUM_PHASES 3

using namespace craig;


void printHelp(const char * pname, ::ofstream &fd) {
    fd << "CRAIG v. " << CRAIG_VERSION << " tool for listing unreachable locations.\n";
    fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
    cerr << "Usage : " << pname << " [Options] FASTA_FILE ANNOT_FILE > UNREACHABLE_LIST\n\n";
    fd << "positional arguments:\n";
    ArgParseUtils::displayOption(fd, "FASTA_FILE", "Name of the file containing the input query sequence(s) in fasta format");
    ArgParseUtils::displayOption(fd,  "ANNOT_FILE", "Name of the file containing genes in locs format");
    fd << "optional arguments:\n";
    ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
    ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
    ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
    ArgParseUtils::displayOption(fd, "--resources=RES_FILE", "Resources file RES_FILE will be used instead of default one in model directory");
    ArgParseUtils::displayOption(fd, "-utr --have-utrs", "Genes contain UTR regions.");
    ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "prefix for fetching files used as input resources");
    ArgParseUtils::displayOption(fd, "--species=SPECIES", "Species name, given by the common suffix of file present in the models subdirectory. Possible values are: human, celegans, human-est and so on [human]");
    ArgParseUtils::displayOption(fd, "UNREACHABLE_LIST", "The list of genes that are not supported by available genomic signal the evidence");
    fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";

}

bool verbose = false;

int main(int argc, char *argv[]) {
    ::string modelPath;
std:string seqFileName = "", annotFile = "";
    ::string resourcesFile = "", prefix_evidence("");
    std::string topology = "partial";
    std::string format("locs");
    std::string species = "human";
    int i, num_sigfilters = 4;
    bool haveUTRs = false;
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
            else if(!strncmp(argv[i], "--resources=", 12)) {
                resourcesFile = ::string(argv[i] + 12);
            }
            else if(!strncmp(argv[i], "--prefix-evidence=", 18)) {
                prefix_evidence = string(argv[i] + 18);
            }
            else if(!strncmp(argv[i], "-utr", 4) ||
                    !strncmp(argv[i], "--have-utrs", 11)) {
                haveUTRs = true;
            }
            else if(!strncmp(argv[i], "--help", 6)
                    || !strncmp(argv[i], "-h", 2)) {
                printHelp("list_unreachable_genes", (std::ofstream &)cout);
                exit(0);
            }
            else if(!strncmp(argv[i], "--version", 9)) {
                PRINT_VERSION(cerr, "list_unreachable_genes", "program to list the set of genes with no supporting signal evidence");
                PRINT_DISCLAIMER(cerr, "list_unreachable_genes");
                exit(0);
            }
            else if(!strncmp(argv[i], "--species=", 10)) {
                species = std::string(argv[i] + 10);
            }
            else
                break;
        }

        if(argc - i < 2)
            throw EXCEPTION(BAD_USAGE, "insufficient arguments");

        if(haveUTRs)
            num_sigfilters = 6;

        seqFileName = std::string(argv[i++]);
        annotFile = std::string(argv[i++]);

        modelPath = species.find('/', 0) == std::string::npos ?
            modelPath + species + "." :
            species + ".";

        if(!resourcesFile.length())
            resourcesFile = modelPath + "resources";

        std::ifstream reStream(resourcesFile.c_str());
        std::ifstream flStream((modelPath + "filters").c_str());
        std::ifstream ftStream((modelPath + "features").c_str());

        map<std::string, std::string> resource_subs;
        resource_subs["PREFIX_EVIDENCE"] = prefix_evidence;
        resource_subs["~"] = modelPath;
        ResourceEngine re(reStream, &resource_subs);

        Sigma *sigma = (Sigma *)re.getResource(std::string("dna-alpha"));
        FSM & fsm = *(FSM *)re.getResource(topology);
        fsm.setParseStrand((TStrand)BOTH_STRANDS);
        GeneEvaluator evaluator(fsm);
        FilterEngine fe(re, flStream, sig_info, num_sigfilters);
        FeatureEngine fte(fsm, fe, ftStream);

        //loading annotSeqs and gene annotations
        InpFile seqFile("default", seqFileName, FASTA, sigma, true);
        list<Sequence *> & annotSeqs = (list<Sequence *> &)seqFile.sequences();

        list<Gene> geneSet;
        GeneUtils::loadGenes(geneSet, annotFile.c_str(), false, false, annotSeqs, sigma);

        // Initializing special variables
        TypedFilter<EdgeInst> **signals = new TypedFilter<EdgeInst> * [NUM_EDGE_INSTS];
        for(i = 0; i < NUM_EDGE_INSTS; i++)
            signals[i] = (TypedFilter<EdgeInst> *)fe.findFilter(sig_info[0][i]);

        TypedFilter<UCHAR> *contexts = (TypedFilter<UCHAR> *)fe.findFilter("GC-Content");

        Lattice lattice(fe, fsm, NUM_PHASES, evaluator, MAX_NUM_EXONS,
                        INTERGENIC, signals, contexts, NULL);

        for(list<Sequence *>::iterator cit = annotSeqs.begin(); cit != annotSeqs.end(); cit++) {
            Sequence &c = *(*cit);
            fe.setSequence(c);

            for(int s = 0; s < NUM_STRANDS; s++)
                fe.computeSeqFilters((TStrand)s);

            GeneUtils::annotSeq2Tags(c.getTags(), &fsm, contexts, signals, c);

            list<Gene *> &listGenes = (list<Gene *> &)c.getAnnotBioFeats();
            vector<Transcript *> transcripts;
            list<Gene *>::iterator it;

            for(it = listGenes.begin(); it != listGenes.end(); it++) {
                Gene &gene = *(*it);
                for(int j = 0; j < gene.transcripts().size(); j++) {
                    Transcript &t = gene.transcripts()[j];
                    transcripts.push_back(&t);
                }
            }

            vector<bool> reachSeqTags;
            GeneTagPrinter printer;

            cout << c.id() << " has "<< listGenes.size() << " genes are reachable? " << lattice.isReachable(c.getTags(), &reachSeqTags) << endl;

            for(int j = 0; j < transcripts.size(); j++) {
                Transcript &t = *transcripts[j];
                printer.displayTranscript(format, (std::ofstream &)cout, c.id(), t);

                cout << "\t" << "transcript " << t.getId() << "(" << t.transcriptLen() << ")" << " is reachable? " << reachSeqTags[j] << endl;
            }

            fe.deleteSeqFilters();
            fe.releaseSequence();
            re.deleteResourceSeqContents();
        }

        // Storing parameters
        reStream.close();
        flStream.close();
        // releasing tags
        SeqTags::releaseMem();
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
