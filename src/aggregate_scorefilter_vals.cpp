#include <iostream>
#include <streambuf>
#include <fstream>
#include <string>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <vector>
#include <cfloat>
#include <map>
#include <iomanip>
#include "ArgParseUtils.h"

using namespace std;
using namespace craig;
using namespace lless;

void printHelp(const char * pname, ::ofstream &fd) {
  fd << "CRAIG v. " << CRAIG_VERSION << " tool for aggregating filter scores for one sequence spread over multiple entries contained in the STDIN.\n";
  fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
  fd << "Usage : " << pname << " [options] < SCOREFILTER_VALS > AGGREGATED_SCOREFILTER_VALS\n\n";
  fd << "optional arguments:\n";
  ArgParseUtils::displayOption(fd, "--stranded", "Scoring Filter has reverse strand information. If not specified, the natural complement is assumed");
  ArgParseUtils::displayOption(fd, "--num-cols=NUM_COLS", "Scoring Filter has value that spawns more than one column [1]");
  ArgParseUtils::displayOption(fd, "--static-vals", "Scoring Filter has value that is static, i.e. it won't accumulate and if it appears in two different entries it should be the same(it will be tested) [1]");
  ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
  ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
  ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
}


// only one contig is allowed
void createFilterValArrays(vector<float>*** v, int clen, int num_cols) {
  
  (*v) = new vector<float> * [2];
  for(int i = 0; i < 2; i++) {
    (*v)[i] = new vector<float> [clen + 1];
    for(int j = 0; j <= clen; j++)
      (*v)[i][j] = vector<float>(num_cols, 0);
  }
}

void freeFilterValArrays(vector<float> *** v) {
  for(int i = 0; i < 2; i++)
    delete [] (*v)[i];
  delete [] *v;
}

void printScoreFilterVals(vector<float> **v, char* cid, int num_cols, 
			  int clen, ostream &os, bool stranded) {
			  
  bool empty = true;
  for(int t = 0; t < 2; t++) {
    if(!t) {
      os << ">" << cid << "\t" << clen;
      if(num_cols > 1) os << "\t" << num_cols;
      os << "\n";
    }
    if(t && (stranded == true || empty))
      os << "//\n";
    
    int i = 2;
    vector<float> last_fval = v[t][1];
    int last_pos = 1;
    
    for( ; i <= clen; i++) {
      if(v[t][i] != last_fval) {
	bool iszero = true;
	for(int nc = 0; nc < num_cols; nc++)
	  if(last_fval[nc]) {
	    iszero = false;
	    break;
	  }
	if(!iszero) {
	  os << last_pos;
	  if(last_pos < i - 1) 
	    os <<  ".." << i - 1;
	  for(int nc = 0; nc < num_cols; nc++)
	    os << "\t" << setprecision(10) << last_fval[nc];
	  os << "\n";
	  empty = false;
	}

	last_fval = v[t][i];
	last_pos = i;
      }
    }
    //       print STDERR "\"$t $cid \"$last_fval\" $last_pos \"\n";
    bool iszero = true;
    for(int nc = 0; nc < num_cols; nc++)
      if(last_fval[nc]) {
	iszero = false;
	break;
      }
    if(!iszero) {
      os << last_pos;
      if(last_pos < clen)
	os << ".." << clen;
      
      for(int nc = 0; nc < num_cols; nc++)
	os << "\t" << setprecision(10) <<  last_fval[nc];
      os << "\n";
      empty = false;
    }
  }
}


bool verbose = false;

int main(int argc, char *argv[]) {
  string line;
  char cid[1000] = "";
  int clen; 
  bool stranded = false, static_vals = false;;
  bool arrays_created = false;
  int i;
  int num_cols = 1;

  try {
    for(i = 1; i < argc; i++) {
      if(!strncmp(argv[i], "--stranded", 10))
	stranded = true;
      else if(!strncmp(argv[i], "--static-vals", 13))
	static_vals = true;
      else if(!strncmp(argv[i], "--num-cols=", 11))
	sscanf(argv[i] + 11, "%d", &num_cols);
      else if(!strncmp(argv[i], "--verbose", 9) 
	      || !strncmp(argv[i], "-v", 2)) {  
        verbose = true;  
      }
      else if(!strncmp(argv[i], "--help", 6)
	      || !strncmp(argv[i], "-h", 2)) {
        printHelp("craig", (std::ofstream &)cout); 
        exit(0);  
      }
      else if(!strncmp(argv[i], "--version", 9)) { 
        PRINT_VERSION(cerr, "aggregate_scorefilter_vals", "a tool for aggregating filter values from multiple entries");
        PRINT_DISCLAIMER(cerr, "craig"); 
        exit(0);  
      }
      else  
	throw EXCEPTION(BAD_USAGE, string("Unrecognized argument ") + argv[i]);
      
    }    

    bool fileOk = std::getline(cin, line);

    vector<float> ** v;
    while(fileOk && sscanf(line.c_str(), ">%s %d", cid, &clen) == 2) {
      if(!arrays_created) {
	createFilterValArrays(&v, clen, num_cols);
	arrays_created = true;
      }
      int strand = 0;
      vector<float> sigVals(num_cols);
      
      while((fileOk = std::getline(cin, line)) &&
	  sscanf(line.c_str(), ">%s %d", cid, &clen) != 2) {
	int sigPos1 = 1, sigPos2 = 0;
        std::istringstream isline(line);
	std::string sigPos;
	isline >> sigPos;
	for(int nc = 0; nc < num_cols; nc++)
	  isline >> sigVals[nc];

	if(sscanf(sigPos.c_str(), "%d..%d", &sigPos1, &sigPos2) == 2)
	  ;
	else if(sscanf(sigPos.c_str(), "%d", &sigPos1) == 1)
 	  sigPos2 = sigPos1;
	else if(!line.compare("//")) {
	  sigPos1 = 1; 
	  strand = 1;
	  if(!stranded)
	    throw EXCEPTION(BAD_USAGE, "Non stranded filter contains rev strand info");
	  continue;
	}
	else
	  throw EXCEPTION(BAD_USAGE, "Unrecognized scoring filter format "+line);
	
	for(int p = sigPos1; p <= sigPos2; p++)
	  if(static_vals)
	    for(int nc = 0; nc < num_cols; nc++)
	      v[strand][p][nc] = sigVals[nc];
	  else
	    for(int nc = 0; nc < num_cols; nc++)
	      v[strand][p][nc] += sigVals[nc];
      }
      
      // check whether the filter did not contain a comp strand info and so 
      // the comp strand has reversed info from the fwd.
      //    if(strand == 0)
      //      for(int p = 1; p <= clen; p++)
      //    	v[1][p] = v[strand][clen - p + 1];
    }
    printScoreFilterVals(v, cid, num_cols, clen, cout, stranded);
    freeFilterValArrays(&v);
    cerr << "Done\n";

  } catch(exception *e) { 
    perror(e->what());
    
    GenLibExcept *gle = (GenLibExcept *)e;
    if(gle->error() == BAD_USAGE)
      cerr << "use --help for more information\n";
    
    delete e;
    exit(1);
  }
}



