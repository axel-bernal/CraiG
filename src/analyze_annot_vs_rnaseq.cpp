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
  fd << "CRAIG v. " << CRAIG_VERSION <<  "tool to estimate the UTRs of transcripts based on RNA-Seq evidence.\nThe tool uses some heuristics whose objective is to reduce the number of false positives.\nThe output can be used as training set for craigTrain\n\n";
  fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
  fd << "  usage :" << pname << " [options] FASTA_FILE ANNOT_FILE PRED_FILE\n\n";
  ArgParseUtils::displayOption(fd, "FASTA_FILE", "Name of the file containing the input query sequence(s) in fasta format");
  ArgParseUtils::displayOption(fd,  "ANNOT_FILE", "Name of the file containing annotated genes in locs format");
  ArgParseUtils::displayOption(fd,  "PRED_FILE", "Name of the file containing predicted genes in locs format");
  fd << "optional arguments:\n";
  ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
  ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
  ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
  ArgParseUtils::displayOption(fd, "--utr-char=UTR_CHAR", "Character used for plotting UTR information[+]");
  ArgParseUtils::displayOption(fd, "-a --abs-coords", "Report is made using absolute chromosome coordinates");
  ArgParseUtils::displayOption(fd, "--fasta-ids=FASTA_ID_LIST", "Comma separated list of ids in FASTA file to be processed");
  ArgParseUtils::displayOption(fd, "-s --stranded", "RNA-Seq data is stranded");
  ArgParseUtils::displayOption(fd, "--prefix-evidence=PREFIX_EVIDENCE", "prefix for fetching files used as input resources. According to naming convention, PREFIX_EVIDENCE.xsignal and PREFIX_EVIDENCE.xsigscores are the files containing the TSS and PSS signals and scores, respectively");
  fd << "Report bugs to <" << AUTHOR_EMAIL << ">\n";
}

bool verbose = false;

int main(int argc, char *argv[]) {
  ::string fastaFileName(""), annotFileName(""), predFileName("");
  ::string prefix_evidence = "";
  ::vector<string> fastaIds;
  bool abs_coords = false;
  char utr_char = '+';
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
      else if(!strncmp(argv[i], "--fasta-ids=", 12)) {
	::string fidsString = string(argv[i] + 12);
	fastaIds.push_back(fidsString);
      }
      else if(!strncmp(argv[i], "--help", 6) 
              || !strncmp(argv[i], "-h", 2)) {
        printHelp("analyze_annot_vs_rnaseq", (std::ofstream &)cout);
        exit(0);  
      }
      else if(!strncmp(argv[i], "--version", 9)) { 
        PRINT_VERSION(cerr, "analyze_annot_vs_rnaseq", "tool for reporting AS events in RNA-seq data");
        PRINT_DISCLAIMER(cerr, "analyze_annot_vs_rnaseq"); 
        exit(0);  
      }
      else if(!strncmp(argv[i], "--stranded", 10) 
              || !strncmp(argv[i], "-s", 2)) {
	stranded = true;
      }
      else if(!strncmp(argv[i], "--utr-char=", 11)) {
        sscanf(argv[i] + 11, "%c", &utr_char);
      }
      else if(!strncmp(argv[i], "-a", 2) || !strncmp(argv[i], "--abs-coords=", 13)) {
	abs_coords = true;
      }
      else
        break;
    }

    if(argc - i < 2)
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
    
    InpFile xsignalsFile("default", prefix_evidence+".xsignal", XFASTA, &evsigma, true, false);
    InpFile xsigscoresFile("default", prefix_evidence+".xsigscores", EXTENDED_MULTI_INTSCORE, NULL, true, false);

    InpFile fastaFile("default", fastaFileName, FASTA, &sigma, true);     
    list<Sequence *> &oFastaSeqs = (list<Sequence *> &)fastaFile.sequences();
    list<Gene> geneSet, pgeneSet;
    GeneUtils::loadGenes(geneSet, annotFileName.c_str(), false, false, 
			 oFastaSeqs, &sigma);
    GeneUtils::loadGenes(pgeneSet, predFileName.c_str(), false, false, 
			 oFastaSeqs, &sigma, PRED_SET);

    list<Sequence *> *fastaSequences = &oFastaSeqs;
    list<Sequence *> qFastaSeqs;

    if(fastaIds.size()) {
      for(i = 0; i < fastaIds.size(); i++) {
	BasicSeq *c = fastaFile.findSeq(fastaIds[i]);
	if(!c) 
	  throw EXCEPTION(CONTIG_UNAVAILABLE, fastaIds[i]);
	qFastaSeqs.push_back((Sequence *)c);
      }
      fastaSequences = &qFastaSeqs;
    }
    
    list<Sequence *>::iterator cit = fastaSequences->begin();
    
    cerr << "Processing sequences" << endl;
    
    std::tr1::hash<std::string> hash_fn;

    for(; cit != fastaSequences->end(); cit++) {
      Sequence &c = *(*cit);
      list<Gene *> &listGenes = (list<Gene *> &)c.getAnnotBioFeats();
      list<Gene *> &plistGenes = (list<Gene *> &)c.getAnnotBioFeats(PRED_SET);
      fe.setSequence(c);
      srand ( hash_fn(c.id()) );//unsigned ( time (NULL) ) );

      int offset = 1; 
      string cid = c.id();
      if(abs_coords) {
	boost::RegEx rExscont("^(\\S+)\\.([^\\.]+)$");
	bool match = rExscont.Match(c.id());
	
	if(match) {
	  cid = rExscont[1];
	  if(!sscanf(rExscont[2].c_str(), "%d", &offset)) 
	    assert(0);
	}
	else cerr << "Absolute coordinates not available\nProceed with relative ones";
      }

      cout << ">" << c.id() << " " << c.length() << endl;
      cout << "Position= " <<  cid << ":" << offset << ".." << offset + c.length() - 1 << endl;
      
      for(int s = 0; s < NUM_STRANDS; s++) 
        fe.computeSeqFilters((TStrand)s);
      
      int numStrands = stranded ? NUM_STRANDS : 1;
      vector<float> *contig_cov = new vector<float> [numStrands];
      vector<float> *junction_cov = new vector<float> [numStrands];

      RSqUtils::computeSeqCoverage(c, &coverage,
				   contig_cov,
				   stranded);
      
      RSqUtils::computeJunctionCoverage(c, &junctions,
					junction_cov,
					stranded);
      
      vector<double> intron_mean(numStrands), exon_mean(numStrands);
      vector<RSqChangePoint> *pchanges = new vector<RSqChangePoint> [numStrands];

      // find ChangePoints
      MultiScoreSeq<int> *xsigscore_c = (MultiScoreSeq<int> *)xsigscoresFile.findSeq(c.id());
      Sequence *xsignal_c = (Sequence *)xsignalsFile.findSeq(c.id());
      if(!xsigscore_c || !xsignal_c) {
	assert(0);
	throw EXCEPTION(CONTIG_UNAVAILABLE, c.id());      
      }
      
      /*      for(int s = 0; s < numStrands; s++)  {
	cout << "xsignal strand " << s << "\n";
	for(i = 0; i < c.length(); i++) {	
	  char symbol = (*xsignal_c)((TStrand)s, i);
	  if(symbol != '-')
	    cout << s << "," << i << " " << symbol<< endl;
	}
	}*/

      cerr.precision(4);
      for(int s = 0; s < numStrands; s++)  {
	TStrand strand = (TStrand)s;
	pair<float, float> covmm = Utils::findMinMax(contig_cov[s], 0, 
						     contig_cov[s].size());
	intron_mean[s] = (covmm.first > 1 ? log(covmm.first) : 1e-10);
	exon_mean[s] = 2*(covmm.second > 1 ? log(covmm.second) : 1e-10)/3.0;
	
	cerr << "intron mean " << intron_mean[s] << endl;
	cerr << "exon mean " << exon_mean[s] << endl;
	
	for(i = 1; i < contig_cov[s].size(); i++) {
	  contig_cov[s][i] += junction_cov[s][i];
	  float cov = contig_cov[s][i];
	  contig_cov[s][i] = (cov > 1 ? log(cov) : 1e-10);
	}

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
	  //	  cerr << pos << " " << mean_l << " " << mean_r << endl;
	  if(strand == STRAND_COMP) 
	    pchanges[s].insert(pchanges[s].begin(), RSqChangePoint(pos, 100, mean_l, mean_r));
	  else
	    pchanges[s].push_back(RSqChangePoint(pos, 100, mean_l, mean_r));
	}
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
      
      transcripts.sort(greater<Transcript *>());

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
	next_begin = (nit == transcripts.end() ? 
		      contig_cov[gstrand].size() + 40 : 
		      (*nit)->codingBegin() - 1);


	cout << "GeneID= " << transcript->getId() << endl << "Strand= " << (transcript->getStrand() == STRAND_FWD ? "Sense" : "Antisense") << endl;
	
	cout << (transcript->getStrand() == STRAND_FWD ? "5UTR=" : "3UTR=") << transcript->begin() + offset - 1 << endl << (transcript->getStrand() == STRAND_FWD ? "5CDS=" : "3CDS=") << tcod_beg + offset - 1 << endl;
	cout << (transcript->getStrand() == STRAND_FWD ? "3UTR=" : "5UTR=") << tcod_end + offset - 1 << endl << (transcript->getStrand() == STRAND_FWD ? "3CDS=" : "5CDS=") << transcript->end() + offset - 1 << endl << endl;
	cout << (transcript->getStrand() == STRAND_FWD ? "5UTR2START" : "3UTR2STOP") << endl;

	RSqUtils::changePoints2Histogram(pchanges[gstrand], last_end + offset - 1, tcod_beg + offset, transcript->begin() + offset, utr_char, offset);

	cout << endl << (transcript->getStrand() == STRAND_FWD ? "STOP23UTR" : "START25UTR") << endl;
	
	RSqUtils::changePoints2Histogram(pchanges[gstrand], tcod_end + offset, next_begin + offset, transcript->end() + offset + 1, utr_char, offset);
	last_end = transcript->end() + 1;
	cout << endl;
      }
      
      delete [] contig_cov;
      delete [] junction_cov;
      
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
  
