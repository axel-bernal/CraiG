/*****************************************************************************
 *   print rnaseq_intron filters
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

using namespace lless;
using namespace std;
using namespace google;

typedef dense_hash_map<pair<int, int>, pair<int, bool>, IntPairHash, IntPairEqual > TJunctionMap;
typedef dense_hash_map<string, dense_hash_map<string, TJunctionMap> > TSeqJunctions;

void printHelp(const char * pname) {
    cerr << "CRAIG v. " << CRAIG_VERSION << " tool (in development)  for processing RNAseq data.\nWritten by Axel E. Bernal (abernal@seas.upenn.edu)\n\n";
    cerr << "  usage :" << pname << " [options] coverage junctions\n\n";
    cerr << "  coverage contains the coverage and junctions contain the splice junctions\n";
    cerr << "  options:\n" ;
    cerr << "    --version\t\tPrint version name and license information\n";
    cerr << "    -h --help\t\tPrint this message\n";
    cerr << "    -v --verbose\t\tTurns on output of debugging information\n";
    cerr << "  Report bugs to <abernal@seas.upenn.edu>\n\n";
}
/*
  template<class T> void printIntervalMap(T map) {
  TSeqReads::iterator it = map.begin();
  dense_hash_map<string, dense_hash_map<string, T> >::
  iterator it = map.begin();
  while(it != map.end()) {
  cout << "Sequence " << it->first << endl;
  dense_hash_map<string, T>::iterator sit = it->second.begin();

  while(sit != it->second.end()) {
  typename T::iterator iit = sit->second.begin();
  cout << "Experiment " << sit->first << endl;

  while(iit != sit->second.end()) {
  discrete_interval<int> boundaries = iit->first;
  cout << "\t[" << boundaries.lower() << " - " << boundaries.upper() << "]"
  << ": " << (*iit++).second << endl;
  }
  sit++;
  }
  it++;
  }
  }
*/

void printFilter(TSeqJunctions & introns, list<Sequence *> & annotSeqs) {
    list<Sequence *>::iterator cit = annotSeqs.begin();
    for( ; cit != annotSeqs.end(); cit++) {
        Sequence &c = *(*cit);

        TSeqJunctions::iterator it = introns.find(c.id());
        if(it == introns.end())
            throw std::runtime_error("no id found\n");

        cout << ">" << c.id() << "\t" << c.length() << "\t3\n";
        dense_hash_map<string, TJunctionMap>::iterator sit = it->second.begin();

        while(sit != it->second.end()) {
            TJunctionMap::iterator jit = sit->second.begin();
            //	cout << "Experiment " << sit->first << endl;

            while(jit != sit->second.end()) {
                pair<int, int> boundaries = jit->first;
                cout << "\t[" << boundaries.first << " " << boundaries.second << "]"
                     << ": " << (*jit++).second.first << " " << jit->second.second << endl;
            }
            sit++;
        }
        it++;
    }
}

void readJunctions(::ifstream &jstream, TSeqJunctions &introns)  {
    // Now we can use class RnSqIntron with interval containers:

    //  boost::regex delim("\\s+");
    boost::RegEx rExDelim("\\s+");
    std::string line;
    std::getline(jstream, line); // header
    introns.set_empty_key("");

    while(std::getline(jstream, line), !jstream.eof()) {
        ::vector<std::string> cols;
        rExDelim.Split(cols, line, 0, 100);
        //    boost::algorithm::split_regex(cols, line, delim);

        int start, end, score;
        bool known;
        sscanf(cols[6].c_str(),"%d", &score);

        if(score < 2) continue;

        sscanf(cols[4].c_str(),"%d", &start);
        sscanf(cols[5].c_str(),"%d", &end);
        sscanf(cols[7].c_str(),"%d", &known);
        //    cout << cols[0] << " " << cols[3] << " " << start << " " << end << " " << cols[1] << cols[2] << endl;

        TSeqJunctions::iterator it = introns.find(cols[3]);

        if(it == introns.end()) {
            introns[cols[3]] = dense_hash_map<string, TJunctionMap>();
            introns[cols[3]].set_empty_key("");
        }

        dense_hash_map<string, TJunctionMap>::iterator sit = introns[cols[3]].find(cols[1]+cols[2]);
        if(sit == introns[cols[3]].end()) {
            introns[cols[3]][cols[1]+cols[2]] = TJunctionMap();
            introns[cols[3]][cols[1]+cols[2]].set_empty_key(pair<int,int>(0,0));
        }

        TJunctionMap::iterator jit = introns[cols[3]][cols[1]+cols[2]].find(pair<int,int>(start, end));
        if(jit != introns[cols[3]][cols[1]+cols[2]].end())  {
            jit->second = pair<double, bool>(jit->second.first + score, known);
            if(jit->second.second != known)
                throw std::runtime_error("known intron?");
        }
        else
            introns[cols[3]][cols[1]+cols[2]][pair<int, int>(start, end)] = pair<double,bool>(score, known);
    }
}

bool verbose = false;

int main(int argc, char *argv[]) {

    int i = 1;
    std::string coverageFile, junctionsFile;
    std::string seqFileName("");
    try {
        if(!strncmp(argv[i], "--help", 6)
           || !strncmp(argv[i], "-h", 2)) {
            printHelp("train_cimm");
            exit(0);
        }

        if(argc - i < 3)
            throw EXCEPTION(BAD_USAGE, "insufficient arguments");

        coverageFile = std::string(argv[i++]);
        junctionsFile = std::string(argv[i++]);
        seqFileName = std::string(argv[i++]);

        ::ifstream cstream(coverageFile.c_str()), jstream(junctionsFile.c_str());
        DNASigma sigma;

        InpFile annotSeqFile("default", seqFileName, FASTA, &sigma, true);
        //loading annotSeqs
        list<Sequence *> & annotSeqs = (list<Sequence *> &)annotSeqFile.sequences();

        TSeqJunctions introns;
        readJunctions(jstream, introns);
        //    printIntervalMap<TJunctionMap>(introns);
        printFilter(introns, annotSeqs);

        cstream.close();
        jstream.close();

    } catch(exception *e) {
        perror(e->what());

        GenLibExcept *gle = (GenLibExcept *)e;
        if(gle->error() == BAD_USAGE)
            cerr << "use --help for more information\n";

        delete e;
        exit(1);
    }

    return 0;
}
