/****************************************************************************
* FilterEngine.h - part of the lless namespace, a general purpose
*                  linear semi-markov structure prediction library
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
*
****************************************************************************/

#ifndef _FILTER_ENGINE_H_
#define _FILTER_ENGINE_H_
#include "FSM.h"
#include "NodeInst.h"
#include "Filter.h"
#include "EdgeInst.h"
#include "ObjectFactory.h"
#include "Sequence.h"
#include <vector>

namespace lless {

  class ResourceEngine;

  /**
   * The class FilterEngine read Filter Header definitions and creates Filter
   * objects dynamically.
   * In the same way as Resource objects, Filter objects can also be defined
   * in an input filter file. The FeatureEngine class in this case, performs a
   * similar task to its Resource counterpart, the ResourceEngine class.
   * The general format for Filter definitions is the following:\n\n
   * Filter    filter-name     FilterClass     space-separated-parameters \n\n
   * FilterClass is either a Filter or any Class derived from Filter, whereas
   * filter-name is a unique identifier.
   **************************************************************************/

  class FilterEngine {
   protected:
    void **paddedFilterValArrays[NUM_STRANDS];  // padded for buffer overruns
    int numArrayEntries;
    
    std::string filterHeaders;

    std::vector<Filter *> filters;
    std::map<std::string, Filter *> namedFilters;
    std::map<std::string, vector<std::string> > multiFilters;
    std::map<std::string, pair<string, vector<std::string> > > filterSets;
    
    Sequence *C;        // associated annotSeq
    int seqLen;         // sequence length

    ResourceEngine *re;

    ObjectFactory5< Filter, int, vector<std::string> &, int &, ResourceEngine *, FilterEngine *, std::string > filterFactory;  
    
   public:
    FilterEngine() {
      this->re = NULL;
      numArrayEntries = 0;
      resetFilterValArrays();
    }

    FilterEngine(ResourceEngine &re, std::ifstream & fd) {
      vector<std::string> lines;
      std::string line;
      boost::RegEx rExEndofFile("^\\s*\\/\\/\\s*$");
      boost::RegEx rExComment("^\\s*\\#.*");
      
      /*
       * Load filters
       */
      while(std::getline(fd, line), !fd.eof()) {      
	lines.push_back(line);
	if(rExEndofFile.Match(line))
	  break;
      }
      
      initialize(re, lines);

    }

    FilterEngine(ResourceEngine &re, std::ifstream & fd, 
		 std::string **sig_info, int num_sigfilters);
    
    FilterEngine(ResourceEngine &re, vector<std::string> & lines) {
      initialize(re, lines);
    }
    
    void initialize(ResourceEngine &re, vector<std::string> & lines) {
      this->re = &re;
      numArrayEntries = 0;
      resetFilterValArrays();
      registerDefaultFilters();
      readHeaderDefinitions(lines);
      allocateFilterValArrays();
    }

    inline void resetFilterValArrays() {
      for(int strand = 0; strand < NUM_STRANDS; strand++)
	this->paddedFilterValArrays[strand] = NULL;
    }

    /*
     * allocating memory for filter values
     */
    inline void allocateFilterValArrays() {
      for(int strand = 0; strand < NUM_STRANDS; strand++) {
	this->paddedFilterValArrays[strand] = new void *[numArrayEntries];
	
	for(int i = 0; i < numArrayEntries; i++)
	  this->paddedFilterValArrays[strand][i] = NULL;
      }
    }

    inline void setSequence(Sequence &c) {
      this->C = &c;
      this->seqLen = c.length();
    }  

    inline Sequence *getSequence() { return C; }
    inline int seqLength() {    return seqLen; }
    inline vector<Filter *> &getFilters() {  return filters; }
    inline char *getSeq(TStrand strand) { return C->getSeq(strand); }
    inline void saveFilters(std::ofstream &fd) {
      fd << filterHeaders;
    }

    inline ResourceEngine * getResourceEngine() {
      return re;
    }

    inline Filter * findFilter(std::string filterName) {
      std::map<std::string, Filter *>::iterator it;
      it = namedFilters.find(filterName);
      
      if(it == namedFilters.end())
        return NULL;
      
      return it->second;
    }
    
    inline vector<std::string> *multiFilterEntries(string multiFilterName) {
      std::map<string, vector<std::string> >::iterator it;
      it = multiFilters.find(multiFilterName);
      
      if(it == multiFilters.end())
        return NULL;
      
      return &it->second;
    }

    inline pair<string, vector<string> > * findFilterSet(string filterSetName) {
      std::map<string, pair<string, vector<string> > >::iterator it;
      it = filterSets.find(filterSetName);
      
      if(it == filterSets.end())
        return NULL;
      
      return &it->second;
    }
    
    inline Filter * getFilter(std::string filterName) {
      Filter *f = findFilter(filterName);
      
      if(!f) {
        assert(0);
        throw EXCEPTION( PARSE_ERROR, std::string("undefined Filter ") + filterName);
      }
      return f;
    } 

    void registerDefaultFilters();
    bool readFilterHeader(int & id, std::string & line);
    void readHeaderDefinitions(vector<std::string> & filter_lines);
    void setFilter(std::string &name, Filter *f,
		   int id, double power, int period,
		   bool lazy_match, 
		   bool sparse_match,
		   bool log_match);
    void computeSeqFilters(TStrand strand);
    void deleteSeqFilters();
    void releaseSequence();

    virtual ~FilterEngine();
  };

}

#endif
