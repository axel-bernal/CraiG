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
  fd << "CRAIG v. " << CRAIG_VERSION << " tool for listing all possible orfs found in FASTA_FILE.\n";
  fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
  cerr << "Usage : " << pname << " [Options] FASTA_FILE > ORF_LIST\n\n";
  fd << "positional arguments:\n";
  ArgParseUtils::displayOption(fd, "FASTA_FILE", "Name of the file containing the input query sequence(s) in fasta format");
  fd << "optional arguments:\n";
  ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
  ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
  ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
  ArgParseUtils::displayOption(fd, "--min-len=MIN_LEN", "Minimum ORF length to extract [1]");
  ArgParseUtils::displayOption(fd, "-t --report-truncated", "Report orfs without stop codon");  
  ArgParseUtils::displayOption(fd, "--resources=RES_FILE", "Resources file RES_FILE will be used instead of default one in model directory");
  ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "prefix for fetching files used as input resources");
  ArgParseUtils::displayOption(fd, "--species=SPECIES", "Species name, given by the common suffix of file present in the models subdirectory. Possible values are: human, celegans, human-est and so on [human]");
  ArgParseUtils::displayOption(fd, "ORF_LIST", "The list of tags which represent the genes as they appear in the sequence");
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
  int minOrfLen = 1;
  bool reportTruncated = false;
  int i;

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
      else if(!strncmp(argv[i], "--min-len=", 10)) {
	sscanf(argv[i] + 10, "%d", &minOrfLen);
      }
      else if(!strncmp(argv[i], "--report-truncated", 18) 
	      || !strncmp(argv[i], "-t", 2)) {  
	reportTruncated = true;
      }
      else if(!strncmp(argv[i], "--resources=", 12)) {
	resourcesFile = ::string(argv[i] + 12);
      }
      else if(!strncmp(argv[i], "--prefix-evidence=", 18)) {
	prefix_evidence = string(argv[i] + 18);
      }
      else if(!strncmp(argv[i], "--help", 6) 
         || !strncmp(argv[i], "-h", 2)) {
        printHelp("extractORFs", (std::ofstream &)cout); 
        exit(0);  
      }
      else if(!strncmp(argv[i], "--version", 9)) { 
        PRINT_VERSION(cerr, "extractORFs", "program to extract ORFs from sequences");
        PRINT_DISCLAIMER(cerr, "extractORFs"); 
        exit(0);  
      }
      else if(!strncmp(argv[i], "--species=", 10)) {
        species = std::string(argv[i] + 10);
      }
      else
        break;
    }

    if(argc - i < 1)
      throw EXCEPTION(BAD_USAGE, "insufficient arguments");
        
    seqFileName = std::string(argv[i++]);

    modelPath = species.find('/', 0) == std::string::npos ?
      modelPath + species + "." :
      species + ".";

    if(!resourcesFile.length())
      resourcesFile = modelPath + "resources";
    
    std::ifstream reStream(resourcesFile.c_str());
    std::ifstream flStream((modelPath + "filters").c_str());

    map<std::string, std::string> resource_subs;
    resource_subs["PREFIX_EVIDENCE"] = prefix_evidence;
    resource_subs["~"] = modelPath;
    ResourceEngine re(reStream, &resource_subs);

    Sigma *sigma = (Sigma *)re.getResource(std::string("dna-alpha"));
    FSM & fsm = *(FSM *)re.getResource(topology);
    fsm.setParseStrand((TStrand)BOTH_STRANDS);
    GeneEvaluator evaluator(fsm);
    FilterEngine fe(re, flStream);
    GeneTagPrinter printer;

    //loading annotSeqs and gene annotations
    InpFile seqFile("default", seqFileName, FASTA, sigma, true);
    list<Sequence *> & annotSeqs = (list<Sequence *> &)seqFile.sequences();

    // Initializing special variables
    TypedFilter<EdgeInst> **signals = new TypedFilter<EdgeInst> * [NUM_EDGE_INSTS];
    signals[START] = (TypedFilter<EdgeInst> *)fe.getFilter("Start-Signal");
    signals[STOP] = (TypedFilter<EdgeInst> *)fe.getFilter("Stop-Signal");
    signals[DONOR] = (TypedFilter<EdgeInst> *)fe.getFilter("Donor-Signal");
    signals[ACCEPTOR] = (TypedFilter<EdgeInst> *)fe.getFilter("Acceptor-Signal");
          
    TypedFilter<UCHAR> *contexts = (TypedFilter<UCHAR> *)fe.findFilter("GC-Content");    

    list<Sequence *>::iterator cit = annotSeqs.begin(); 
    for( ; cit != annotSeqs.end(); cit++) {
      Sequence &c = *(*cit);
      fe.setSequence(c);

      for(int s = 0; s < NUM_STRANDS; s++) {
	TStrand strand = (TStrand)s;
	signals[START]->allocValArrays(&c, strand);
	signals[START]->computeVals(&c, strand);
	signals[STOP]->allocValArrays(&c, strand);
	signals[STOP]->computeVals(&c, strand);
	
	for(int i = 1; i <= 3; i++) {
	  vector<pair<int, int> > orfs;
	  GeneUtils::extractFrameOrfs(orfs, c, signals, i, strand, minOrfLen, reportTruncated, true);
	  for(int j = 0; j < orfs.size(); j++) {
	    pair<int, int> & orf  = orfs[j];
	    std::ostringstream geneId;
	    geneId << c.id() << "." << strand << i << "." << orf.first;
	    cout << ">" << geneId.str() << "\t" << c.id() << "_" << orf.first << "_" << orf.second << "\n";
	  }
	}
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
