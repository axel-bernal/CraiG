/*****************************************************************************
 * This program computes change points in rnaseq coverage data and assign them
 * a signal occurrence. Either a PAS or a TSS
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
#include "FSM.h"
#include "TagUtils.h"
#include "InpFile.h"
#include "RSqUtils.h"
#include "FL_ScoreFile.h"
#include <boost/regex.hpp>
#include <numeric>
#define NUM_PHASES 3

using namespace craig;

void printHelp(const char * pname, ::ofstream &fd) {
  cerr << "CRAIG v. " << CRAIG_VERSION <<  "tool to compute the signal and sigscore filters in FASTA_FILE.xsignal\nand FASTA_FILE.xsigscore respectively\n";
  fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
  cerr << "  usage :" << pname << " [options] FASTA_FILE\n\n";
  ArgParseUtils::displayOption(fd, "FASTA_FILE", "Name of the file containing the input query sequence(s) in fasta format");
  fd << "optional arguments:\n";
  ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
  ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
  ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
  ArgParseUtils::displayOption(fd, "--num-permutations=NUM_PERMUTATIONS", "Number of permutations to perform to detect change points in coverage");
  ArgParseUtils::displayOption(fd, "--chunk-size=CHUNK_SIZE", "Contigs larger than CHUNK_SIZE will be partition into fragments of CHUNK_SIZE size [10000]");
  ArgParseUtils::displayOption(fd, "--min-cov=MIN_COVERAGE", "Minimum coverage required. If given, regions below that threshold will be ignored. This option speeds up the program dramatically and it is highly recommended to be set when working on entire chromosomes/long contigs [0]");
  ArgParseUtils::displayOption(fd, "--block-len=BLOCK_LEN", "Block length for finding change points[20]");
  ArgParseUtils::displayOption(fd, "--prefix-output=PREFIX_EVIDENCE.OUTPUT", "Prefix for output. The output will be reported in PREFIX_EVIDENCE.xsignal and PREFIX_EVIDENCE.xsigscore");
  ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "Prefix for associated RNAseq files. According to naming convention, junction and small reads must be located in files PREFIX_EVIDENCE.junctions.signal and PREFIX_EVIDENCE.cov respectively. The output will be reported in PREFIX_EVIDENCE.xsignal and PREFIX_EVIDENCE.xsigscore unless PREFIX_EVIDENCE.OUTPUT has been specified with option --prefix-output");
  ArgParseUtils::displayOption(fd, "-s --stranded", "RNA-Seq data is stranded");
  fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";
}

bool verbose = false;

int main(int argc, char *argv[]) {
  ::string fastaFileName("");
  int i;
  ::string prefix_evidence = "", prefix_output = "";
  int min_cov = 0;
  int num_perms = 1000;
  int chunk_size = 10000;
  bool stranded = false;
  int block_len = 20;
  string symbols("XZZX");

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
        printHelp("get_covxsignals", (std::ofstream &)cout);
        exit(0);  
      }
      else if(!strncmp(argv[i], "--version", 9)) { 
        PRINT_VERSION(cerr, "get_covxsignals", "tool for reporting AS events in RNA-seq data");
        PRINT_DISCLAIMER(cerr, "get_covxsignals"); 
        exit(0);  
      }
      else if(!strncmp(argv[i], "--prefix-evidence=", 18)) {
        prefix_evidence = std::string(argv[i] + 18);
      }
      else if(!strncmp(argv[i], "--prefix-output=", 16)) {
        prefix_output = std::string(argv[i] + 16);
      }
      else if(!strncmp(argv[i], "--num-permutations=", 19)) {
        sscanf(argv[i] + 19, "%d", &num_perms);
      }
      else if(!strncmp(argv[i], "--chunk-size=", 13)) {
        sscanf(argv[i] + 13, "%d", &chunk_size);
      }      

      else if(!strncmp(argv[i], "--stranded", 10) 
              || !strncmp(argv[i], "-s", 2)) {
	stranded = true;
      }
      else if(!strncmp(argv[i], "--block-len=", 12)) {
        sscanf(argv[i] + 12, "%d", &block_len);
      }
      else if(!strncmp(argv[i], "--min-cov=", 10)) {
        sscanf(argv[i] + 10, "%d", &min_cov);
      }
      else
        break;
    }

    if(argc - i < 1)
      throw EXCEPTION(BAD_USAGE, "insufficient arguments");
        
    fastaFileName = std::string(argv[i++]);
    
    if(!prefix_output.length())
      prefix_output = prefix_evidence;
    
    std::ofstream sigd((prefix_output + ".xsignal").c_str());
    std::ofstream sigscored((prefix_output + ".xsigscores").c_str());

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
    list<Sequence *>::iterator cit = fastaSequences.begin();

    cerr << "Processing sequences"<< endl;
    
    std::tr1::hash<std::string> hash_fn;
    
    for(; cit != fastaSequences.end(); cit++) {
      Sequence &superC = *(*cit);
      cerr << "CONTIG " <<  superC.id() << "\t" << superC.length() << endl;
      srand ( hash_fn(superC.id()) );
      int nsync_beg, nsync_end;
      bool fragmented = false;
      int this_chunk_size;
      int numStrands = stranded ? NUM_STRANDS : 1;      
      int begin = 0, end;
      string cid = superC.id();
      int clen = superC.length();
      int offset = 1; 
      boost::RegEx rExscont("^(\\S+)\\.SPLIT_CONTIG\\.(\\d+)\\.(\\d+)$");
      bool match = rExscont.Match(superC.id());

      if(match) {
	cid = rExscont[1];
	if(!sscanf(rExscont[2].c_str(), "%d", &clen)) 
	  assert(0);
	if(!sscanf(rExscont[3].c_str(), "%d", &offset)) 
	  assert(0);

	std::ostringstream ostr;      
	chunk_size = superC.length();
	ostr << cid << "_" << offset << "_" << offset + superC.length() - 1;
	std::string new_cid = ostr.str();
	superC.setId(new_cid);
      }

      sigscored << ">" << cid << "\t" << clen << "\t3" << endl;
      sigd << ">" << cid << "\t" << clen << endl;
      
      vector<RSqChangePoint> pchanges;
      vector<vector<char> > 
	sig_filter(NUM_STRANDS, vector<char>(clen + 1, '-'));
      vector< vector<vector<int> > >
	score_filter(NUM_STRANDS, vector<vector<int> >(clen + 1));
      bool empty[] = {true, true};

      while(begin <= superC.length()) {
	// calculating end.
	this_chunk_size = Utils::min(superC.length() - begin, chunk_size);
	end = begin + this_chunk_size + 1;
	nsync_beg = begin + 1;
	nsync_end = end - 1;
	
	cerr << "working on [" << nsync_beg << ", " << nsync_end << "] out of " << superC.length() << endl;
	
	fragmented = (this_chunk_size < superC.length());
	
	Sequence *c = fragmented ?
	  (Sequence *)superC.getSubSequence(nsync_beg, nsync_end) :
	  &superC;
	
	fe.setSequence(*c);
	
	// allocating memory
	for(int s = 0; s < NUM_STRANDS; s++)
	  fe.computeSeqFilters((TStrand)s);
	
	vector<float> *contig_cov = new vector<float> [numStrands];
	vector<float> *junction_cov = new vector<float> [numStrands];
	
	RSqUtils::computeSeqCoverage(*c, &coverage,
				     contig_cov,
				     stranded);
	
	RSqUtils::computeJunctionCoverage(*c, &junctions,
					  junction_cov,
					  stranded);
	
	double intron_mean, exon_mean;
		
	for(int s = 0; s < numStrands; s++)  {
	  TStrand strand = (TStrand)s;
	  pchanges.clear();
	  for(i = 1; i < contig_cov[s].size(); i++) {
	    contig_cov[s][i] += junction_cov[s][i];
	    float cov = contig_cov[s][i];
	    contig_cov[s][i] = (cov > 1 ? log(cov) : 1e-10);
	  }
	  int length = end - begin;
	  
	  pair<float, float> covmm = Utils::findMinMax(contig_cov[s],
						       0, length);
	  intron_mean = covmm.first;
	  exon_mean = 2*covmm.second/3.0;
	  	    
	  int blen = block_len;
	  if(length/block_len < 25)
	    blen = length/25 + 1;

	  vector<float> tmp_cov(length + 40, 0.1);
	  for(i = 1; i < length; i++)
	    tmp_cov[i + 20] = contig_cov[s][i];
	  
	  RSqUtils::computeRSqChangePoints(tmp_cov, pchanges, num_perms,
					   0, tmp_cov.size(),
					   blen, 80, NULL);//, &spl_posmap);
	  
	  RSqUtils::refineRSqChangePoints(tmp_cov, pchanges, num_perms,
					  0, tmp_cov.size(),
					  blen, 80, NULL);
	  
	  for(i = 0; i < pchanges.size(); i++)
	    pchanges[i].pos -= 20;

	  if(pchanges.size()) {
	    if(pchanges[0].pos == 0)
	      pchanges[0].pos = 1;
	    if(pchanges[pchanges.size() - 1].pos == c->length() + 1)
	      pchanges[pchanges.size() - 1].pos = c->length();
	    empty[s] = false;
	  }

	  int c_offset = (strand == STRAND_FWD ? 
			  begin + offset - 1:
			  clen - (end + offset - 1) + 1);
	  
	  RSqUtils::changePoints2Filter(c_offset, 
					symbols, pchanges,
					contig_cov[s], 80.0,
					c->length(),
					sig_filter[s], score_filter[s],
					strand);
	}
	
	if(!stranded) {
	  RSqUtils::changePoints2Filter(clen - (end + offset - 1) + 1,
					symbols, pchanges,
					contig_cov[0], 80.0,
					c->length(),
					sig_filter[1], score_filter[1],
					STRAND_COMP);
	  empty[1] = empty[0];
	}	

	begin = end;

	delete [] contig_cov;
	delete [] junction_cov;
	
	fe.deleteSeqFilters();
	fe.releaseSequence();
	junctionsFile.releaseSeqContents();
	covFile.releaseSeqContents();
	
	if(fragmented)
	  delete c;
	
      }
    
      for(int s = 0; s < NUM_STRANDS; s++)  {
	if(!empty[s]) {
	  RSqUtils::filter2XFastaFile(sig_filter[s], sigd);   
	  RSqUtils::filter2MultiScoreFile(score_filter[s], 3, sigscored);
	}
	else 
	  sigd << "- x " << clen << " 1\n";
	
	if(!s) {
	  sigd << "//" << endl;
	  sigscored << "//" << endl;
	}
      }
    }

    cerr << "Done!\n";
    
    sigd.close();
    sigscored.close();
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
  
