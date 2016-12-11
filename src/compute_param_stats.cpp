/*****************************************************************************
 * This is the main of the compute_param_stats program. It computes the 
 * feature vector associated to both training and validation tag annotation
 * present int he configuration file.
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

void printHelp(const char * pname) {
  cerr << "CRAIG v. " << CRAIG_VERSION << " discriminative gene prediction tool for eukarya.\nWritten by Axel E. Bernal (abernal@seas.upenn.edu)\n\n";
  cerr << "  usage : " << pname << " [options] PARAMS_FILE > avg_vals\n\n";
  cerr << "  PARAMS_FILE is the name of the file containing the gene model parameters;\n             if file parameters is not found, craig assumes the filename\n             to be $(CRAIG_HOME)/models/PARAMS_FILE\n";
  cerr << "  options:\n" ;
  cerr << "    --version\t\tPrint version name and license information\n";
  cerr << "    --logFile=<logFile>\tCompute stats for log file logFile\n";
  cerr << "    -h --help\t\tPrint this message\n";
  cerr << "    -v --verbose\tTurns on output of debugging information\n";
  cerr << " output :";
  cerr << " The scores for each seq tag present in the input sequence\n\n";
  cerr << "  Report bugs to <abernal@seas.upenn.edu>\n\n";
}

bool verbose = false;

void reportMaxValues(GlobalVector &gv, FeatureEngine &fte,
		     ofstream &fd, int numVals = 10) {

  list<pair<FeatVectorInd, double> > maxv;
  gv.computeMaxParamVals(maxv, numVals);
  list<pair<FeatVectorInd, double> >::iterator it2 = maxv.begin();
  
  for(; it2 != maxv.end(); it2++)
    fd << fte.getFeatures()[it2->first.ind]->getName() << "\t" 
	 << fte.getFeatDescriptions()[it2->first.ind] << "\t"
	 << it2->first.first << " " << it2->first.second << "\t"
	 << exp(it2->second) << endl;
}


int main(int argc, char *argv[]) {
  std::string paramsFile(""), logFile("");
  std::string topology("partial");
  int origStrand = BOTH_STRANDS;
  int i;

  try {

    for(i = 1; i < argc; i++) {

      if(!strncmp(argv[i], "--verbose", 9) 
         || !strncmp(argv[i], "-v", 2)) {  
        verbose = true;  
      }
      else if(!strncmp(argv[i], "--help", 6) 
	      || !strncmp(argv[i], "-h", 2)) {
        printHelp("compute_param_stats"); 
        exit(0);  
      }
      else if(!strncmp(argv[i], "--logFile=", 10)) {
	logFile = ::string(argv[i] + 10);
      }
      else if(!strncmp(argv[i], "--version", 9)) { 
        PRINT_VERSION(cerr, "compute_param_stats", "discriminative gene prediction tool for eukarya");
        PRINT_DISCLAIMER(cerr, "compute_param_stats"); 
        exit(0);  
      }
      else
        break;
    }

    if(argc - i < 1)
      throw EXCEPTION(BAD_USAGE, "insufficient arguments");

    paramsFile = std::string(argv[i++]);
    ::ifstream *paramStream, ifstr1(paramsFile.c_str()), ifstr2;
    paramStream = &ifstr1;

    if(!ifstr1.is_open()) {
      if(!getenv("CRAIG_HOME"))
        throw EXCEPTION(NOT_ANNOTATED,
                            "CRAIG_HOME must be initialized!. See README for details");
      
      paramsFile = std::string(getenv("CRAIG_HOME")) + "/models/" + paramsFile;
      ifstr2.open(paramsFile.c_str());
      
      if(!ifstr2.is_open()) 
        throw EXCEPTION(FILE_UNAVAILABLE, paramsFile);

      paramStream = &ifstr2;
    }
    
    // Retrieving resources
    ResourceEngine re(*paramStream);
    FSM & fsm = *(FSM *)re.getResource(topology);
    fsm.setParseStrand((TStrand)origStrand);
    FilterEngine fe(re, *paramStream);
    FeatureEngine fte(fsm, fe, *paramStream);

    GlobalVector params(fte.getFeatures(), paramStream);
    fte.setParamVector(&params);

    /*    cout << "Avg Values\n";    
    list<pair<int, double> > avgs;
    params.computeAvgParamVals(avgs);
    list<pair<int, double> >::iterator it = avgs.begin();    
    for(; it != avgs.end(); it++)
      cout << fte.getFeatures()[it->first]->getName() << "\t"
	   << fte.getFeatDescriptions()[it->first] << "\t"
	   << exp(it->second) << endl;
    */

    if(logFile.length() > 0) {
      std::ifstream logStream(logFile.c_str());
      int dummy;
      logStream >> dummy >> dummy;
      GlobalVector tmp1(fte.getFeatures(), &logStream, 0);
      GlobalVector tmp2(fte.getFeatures(), &logStream, 0);
      GlobalVector sigma_gv(fte.getFeatures(), &logStream, 0, true, 1.0);
      GlobalVector accsigma_gv(fte.getFeatures(), &logStream, 0);
      cout << "Covariance Values\n";
      reportMaxValues(sigma_gv, fte, (ofstream &)cout);
      logStream.close();
    }
    
    cout << "Max Values\n";
    reportMaxValues(params, fte, (ofstream &)cout);
    
  } catch(exception *e) { 
    perror(e->what());

    GenLibExcept *gle = (GenLibExcept *)e;
    if(gle->error() == BAD_USAGE)
      cerr << "use --help for more information\n";

    delete e;
    exit(1);
  }
}
