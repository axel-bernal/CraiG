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
  fd << "CRAIG v. " << CRAIG_VERSION << " tool for aggregating signal scores for one sequence spread over multiple entries contained in the STDIN.\n";
  fd << "Written by Axel E. Bernal (" << AUTHOR_EMAIL << ")\n\n";
  fd << "Usage : " << pname << " [options] < SIGNALSCORE_VALS > AGGREGATED_SIGNALSCORE_VALS\n\n";
  fd << "optional arguments:\n";
  ArgParseUtils::displayOption(fd, "--stranded", "Filter Signal has reverse strand information. If not specified, the natural complement is assumed");
  ArgParseUtils::displayOption(fd, "--nosignal-char", "Character used to denote the lack of signal information[-]");
  ArgParseUtils::displayOption(fd, "--version", "Print version name and license information");
  ArgParseUtils::displayOption(fd, "-h --help", "Show this message and exit");
  ArgParseUtils::displayOption(fd, "-v --verbose", "Turns on output of debugging information");
}


// only one contig is allowed
void createFilterValArrays(char*** v, int clen, char fillchr) {
  (*v) = new char * [2];
  for(int i = 0; i < 2; i++) {
    (*v)[i] = new char [clen + 1];
    for(int j = 0; j <= clen; j++)
      (*v)[i][j] = fillchr;
  }
}


void freeFilterValArrays(char*** v) {
  for(int i = 0; i < 2; i++) 
    delete [] (*v)[i];
  delete [] *v;
}

void printSignalFilterVals(char** v, char* cid, int clen, 
			  ostream &os, bool stranded) {
  bool empty = true;
  for(int t = 0; t < 2; t++) {
    if(!t) 
      os << ">" << cid << "\t" << clen << "\n";
    
    if(t && (stranded == true || empty))
      os << "//\n";
    
    int i = 2;
    char last_fval = v[t][1];
    int last_pos = 1;
    
    for( ; i <= clen; i++) {
      if(v[t][i] != last_fval) {
	os << last_fval;
	if(last_pos < i - 1)
	  os  << " x " << i - last_pos << " " << last_pos;
	os << "\n";

	empty = false;
	last_fval = v[t][i];
	last_pos = i;
      }
    }
    //       print STDERR "\"$t $cid \"$last_fval\" $last_pos \"\n";
    os << last_fval;
    if(last_pos < clen)
      os  << " x " << clen + 1 - last_pos << " " << last_pos;
    os << "\n";
    empty = false;
  }
}


bool verbose = false;

int main(int argc, char *argv[]) {
  string line;
  char** v;
  char cid[1000] = "";
  int clen; 
  bool stranded = false;
  bool arraysCreated = false;
  char fillchr = '-';
  int i;

  try {
    for(i = 1; i < argc; i++) {
      if(!strncmp(argv[i], "--stranded", 10))
	stranded = true;
      else if(!strncmp(argv[i], "--nosignal-char", 16))
	sscanf(argv[i] + 16, "%c", &fillchr);
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
        PRINT_VERSION(cerr, "aggregate_signalfilter_vals", "a tool for aggregating signal filter values from multiple entries");
        PRINT_DISCLAIMER(cerr, "craig"); 
        exit(0);  
      }
      else  
	throw EXCEPTION(BAD_USAGE, string("Unrecognized argument") + argv[i]);
      
    }    

    bool fileOk = std::getline(cin, line);
    while(fileOk && sscanf(line.c_str(), ">%s %d", cid, &clen) == 2) {
      if(!arraysCreated) {
	createFilterValArrays(&v, clen, fillchr);
	arraysCreated = true;
      }
      int strand = 0;
      int currPos = 1;
      while((fileOk = std::getline(cin, line)) &&
	    sscanf(line.c_str(), ">%s %d", cid, &clen) != 2) {
	int sigPos1 = 1, len = 1, pos;  
	char sigVal;

	if(!line.compare("//")) {
	  currPos = 1; 
	  strand = 1;
	  if(!stranded)
	    throw EXCEPTION(BAD_USAGE, "Non stranded signal filter contains rev strand info");
	}
	else if(sscanf(line.c_str(), "%c x %d %d", &sigVal, &len, &pos) == 3) {
	  if(currPos != pos)
	    throw EXCEPTION(BAD_USAGE, "third column position is inconsistent");
	  if(sigVal != fillchr)
	    for(int p = 1; p <= len; p++) {
	      if(v[strand][currPos] != fillchr && sigVal != fillchr)
		throw EXCEPTION(BAD_USAGE, "character should be nosignal");
	      v[strand][currPos++] = sigVal;
	    }
	  else
	    currPos += len;
	}
	else if(sscanf(line.c_str(), "%c", &sigVal) == 1) {
	  if(sigVal != fillchr) {
	    if(v[strand][currPos] != fillchr && sigVal != fillchr) 
	      throw EXCEPTION(BAD_USAGE, "character should be nosignal");
	    
	    v[strand][currPos++] = sigVal;
	  }
	  else currPos++;
	}
	else
	  throw EXCEPTION(BAD_USAGE, "Unrecognized signal filter format "+line);
      }
      
      // check whether the filter did not contain a comp strand info and so 
      // the comp strand has reversed info from the fwd.
      //    if(strand == 0)
      //      for(int p = 1; p <= clen; p++)
      //    	v[1][p] = v[strand][clen - p + 1];
    }
    printSignalFilterVals(v, cid, clen, cout, stranded);
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



