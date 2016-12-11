/*****************************************************************************
 * This program reports an analysis of the UTR regions in gene annotations 
 * provided as input
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
  fd << "CRAIG v. " << CRAIG_VERSION << " Tool for computing gene coverages using input filter scores containing coverage information.\n";
  fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
  fd << "Usage : " << pname << " [options] FASTA_FILE ANNOT_FILE COVERAGE_FILE > GENE_COVERAGE_VALS\n\n";
  ArgParseUtils::displayOption(fd, "FASTA_FILE", "Name of the file containing the input query sequence(s) in fasta format");
  ArgParseUtils::displayOption(fd,  "ANNOT_FILE", "Name of the file containing annotated genes in locs format");
  ArgParseUtils::displayOption(fd,  "COVERAGE_FILE", "Name of the file containing sequence coverage in score filter format");
  fd << "optional arguments:\n";
  ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
  ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
  ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
  ArgParseUtils::displayOption(fd, "-s --stranded", "RNA-Seq data is stranded");
  fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";
}

bool verbose = false;

int main(int argc, char *argv[]) {
  ::string fastaFileName(""), annotFileName(""), coverageFileName("");
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
      else if(!strncmp(argv[i], "--help", 6) 
              || !strncmp(argv[i], "-h", 2)) {
        printHelp("analyze_annot_vs_rnaseq", (std::ofstream &)cout);
        exit(0);  
      }
      else if(!strncmp(argv[i], "--version", 9)) { 
        PRINT_VERSION(cerr, "get_coverage_genes", "Tool for computing gene coverages using input filter scores containing coverage information");
        PRINT_DISCLAIMER(cerr, "get_coverage_genes"); 
        exit(0);  
      }
      else if(!strncmp(argv[i], "--stranded", 10) 
              || !strncmp(argv[i], "-s", 2)) {
	stranded = true;
      }
      else
        break;
    }

    if(argc - i < 2)
      throw EXCEPTION(BAD_USAGE, "insufficient arguments");
        
    fastaFileName = std::string(argv[i++]);
    annotFileName = std::string(argv[i++]);
    coverageFileName = std::string(argv[i++]);

    DNASigma sigma;
    FilterEngine fe;

    /* Creating filters */    
    InpFile covFile("default", coverageFileName, EXTENDED_LDSCORE, NULL, true, false);
    
    FL_ScoreFile<ScoreSeq<double>, ScoreSeq<double>, double, double> coverage(0,  "RNASeq-Cov", &covFile, 1, 0, FT_DOUBLE);
    fe.setFilter(coverage.getName(), &coverage, 0, 1, 0, false, false, false);
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
    list<Sequence *> &fastaSeqs = (list<Sequence *> &)fastaFile.sequences();
    list<Gene> geneSet;
    GeneUtils::loadGenes(geneSet, annotFileName.c_str(), false, false, 
			 fastaSeqs, &sigma);

    list<Sequence *> *fastaSequences = &fastaSeqs;
    list<Sequence *>::iterator cit = fastaSequences->begin();
    
    cerr << "Processing sequences" << endl;
    
    for(; cit != fastaSequences->end(); cit++) {
      Sequence &c = *(*cit);
      list<Gene *> &listGenes = (list<Gene *> &)c.getAnnotBioFeats();
      
      if(!listGenes.size())
	continue;

      //      cout << "cid " << c.id() << endl;
      fe.setSequence(c);

      for(int s = 0; s < NUM_STRANDS; s++) 
        fe.computeSeqFilters((TStrand)s);
      
      int numStrands = stranded ? NUM_STRANDS : 1;
      vector<float> *contig_cov = new vector<float> [numStrands];

      RSqUtils::computeSeqCoverage(c, &coverage,
				   contig_cov,
				   stranded);
      
      for(int s = 0; s < numStrands; s++)  {
	for(i = 1; i < contig_cov[s].size(); i++)
	  contig_cov[s][i] += contig_cov[s][i - 1];
      }
    
      list<Transcript *> transcripts;
      list<Gene *>::iterator git = listGenes.begin();
      for( ; git != listGenes.end(); git++) {
	Gene &gene = *(*git);
	for(int j = 0; j < gene.transcripts().size(); j++) {
	  Transcript &t = gene.transcripts()[j];
	  transcripts.push_back(&t);
	}
      }
      
      cerr.precision(4);
      list<Transcript *>::iterator it = transcripts.begin();

      for( ; it != transcripts.end(); it++) {
	Transcript &t = *(*it);
	TStrand cstrand = (stranded ? t.getStrand() : STRAND_FWD);
	vector<Exon> t_exons = t.FPexons();
	t_exons.insert(t_exons.end(), t.exons().begin(), t.exons().end());
	t_exons.insert(t_exons.end(), t.TPexons().begin(), t.TPexons().end());
	
	double gene_cov = 0;
	int gene_len = 0;
	for(int i = 0; i < t_exons.size(); i++) {
	  Exon &e = t_exons[i];
	  int beg = e.begin(), end = e.end();
	  /*	  if(t.getStrand() == STRAND_COMP) {
	    cout << "compl " << beg << " " << end << "\n";
	    beg = end;
	    end = e.begin();
	    }*/

	  gene_len += end - beg + 1;
	  //	  cout << "exon " << i << "\t" << beg << "\t" << end << "\t" << contig_cov[cstrand][end] << "\t" << contig_cov[cstrand][beg - 1] << endl;
	  gene_cov += contig_cov[cstrand][end] - contig_cov[cstrand][beg - 1];
	}

	if(gene_len > 0)
	  cout << t.getId() << "\t" << gene_cov/gene_len << endl;
      }
      
      delete [] contig_cov;
      
      fe.deleteSeqFilters();
      fe.releaseSequence();
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
  
