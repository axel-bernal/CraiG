/****************************************************************************
* FT_BinnedSegment.h - part of the lless namespace, a general purpose
*                       linear semi-markov structure prediction library
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

#ifndef FT_BINNED_SEGMENT_FEAT_H
#define FT_BINNED_SEGMENT_FEAT_H

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Sequence.h"
#include "Feature.h"
#include "Filter.h"
#include "FilterEngine.h"
#include "FeatureEngine.h"
#include "FT_Segment.h"
#include "FT_Bin.h"

namespace lless {

  /**
   * FT_BinnedSegment is a subclass of TypedFeature whose value of this 
   * feature (and its derivate classes) is the score obtained by accumulating
   * values of the contained filter over a tag which must be a segment
   * (node with variable length).
   * This feature works well with bimodal filter values, in which scores above
   * a certain threshold (the mean) should be classified differently.
   * Invariants: the contained filter must accumulate.
   *
   ***************************************************************************/
  template <class TClass1, class TClass2> 
  class FT_BinnedSegment : public FT_Segment<TClass2> {
   protected:
    int numBinObjects;
    ULONG domain, range;
    TypedFeature<TClass1> **binObjects;
    double *_cachedBinZeroValues;
    vector<int> accVal4fvals;
    vector<int> binObj4accVals;  //!< which bin object to use for each filter val
    bool countZeroLenBlocks;
    vector<int> *hblockLengths;
    int _numPhases;
    double _featValues[2];
    bool updated;

   public:  
    FT_BinnedSegment(
		     int fInd, 
		     int paramInd, 
		     int parsingFrames,
		     char *name, 
		     FilterEngine *fe, 
		     TypedFilter<TClass2> *filter, 
		     int period,
		     int numBinObjects,
		     bool countZeroLenBlocks,
		     int phaseDependent)
      : FT_Segment<TClass2>(fInd, paramInd, 
			    parsingFrames, name, 
			    fe, filter, period, 0) {
      
      initialize(numBinObjects);
      this->countZeroLenBlocks = countZeroLenBlocks;
      //\todo domain = 0 here, needs to be updated
      _numPhases = phaseDependent ? period : 1;

      _cachedBinZeroValues = new double [numPhases()*domain];
      for(int i = 0; i < numPhases()*domain; i++)
	_cachedBinZeroValues[i] = 0;

      hblockLengths = new vector<int> [numPhases()];
      for(int i = 0; i < numPhases(); i++)
	hblockLengths[i] = vector<int>(domain, 0);

      this->setMaxNumParams(numPhases()*domain, range);      

    }

    FT_BinnedSegment(
		     int fInd, 
		     int paramInd, 
		     vector<std::string> & fargs, 
		     int & offset,
		     bool phaseDependent,
		     FilterEngine *fe, 
		     FeatureEngine *fte
		     )
      : FT_Segment<TClass2>(fInd, paramInd,
			    fargs, offset, 
			    fe, 0) {
      
      int numBinObjects;
      updated = false;

      if(!sscanf(fargs[offset++].c_str(), "%d", &numBinObjects))
        assert(0);    

      initialize(numBinObjects);
      
      TypedFeature<TClass1> *binObj;
      vector<bool> used_fv(this->filter->maxNumFilterValues(), false);
      
      for(int i = 0; i < numBinObjects; i++) {
        ULONG numfvals, fval; 	

        binObj = (TypedFeature<TClass1> *)fte->getFeature(fargs[offset++]);
	binObjects[i] = binObj;
	bool accumulate = Utils::stringToBoolean(fargs[offset++]);
	
	if(!sscanf(fargs[offset++].c_str(), "%ld", &numfvals))
          assert(0);    

	if(numfvals == this->filter->maxNumFilterValues()) {
	  for(int j = 0; j < numfvals; j++) {
	    accVal4fvals[j] = domain;
	    binObj4accVals[domain]= i;
	    used_fv[j] = true;

	    if(!accumulate)
	      domain++;

	  }
	}
	else if(!numfvals) { // assign all non used filter vals to bin i
	  for(int j = 0; j < this->filter->maxNumFilterValues(); j++) {
	    if(used_fv[j])
	      continue;

	    accVal4fvals[j] = domain;
	    binObj4accVals[domain]= i;
	    used_fv[j] = true;

	    if(!accumulate)
	      domain++;

	  }	  
	}
	else {
	  for(int j = 0; j < numfvals; j++) {
	    if(!sscanf(fargs[offset++].c_str(), "%ld", &fval))
	      assert(0);    

	    accVal4fvals[fval] = domain;
	    binObj4accVals[domain]= i;
	    used_fv[fval] = true;
	    
	    if(!accumulate)
	      domain++;
	    
	  }
	}
	
	if(accumulate)
	  domain++;

	if(binObj->maxNumFeatValues() > range)
	  range = binObj->maxNumFeatValues();
	
      }
      
      _numPhases = phaseDependent ? this->period() : 1;


      _cachedBinZeroValues = new double [numPhases()*domain];
      for(int i = 0; i < numPhases()*domain; i++)
	_cachedBinZeroValues[i] = 0;

      hblockLengths = new vector<int> [numPhases()];
      for(int i = 0; i < numPhases(); i++)
	hblockLengths[i] = vector<int>(domain, 0);

      this->setMaxNumParams(numPhases()*domain, range);      
      this->countZeroLenBlocks = true;
      
      if(offset < fargs.size())
	this->countZeroLenBlocks = Utils::stringToBoolean(fargs[offset++]);
      
    }
    
    void initialize(int numBinObjects) {
      domain = range = 0;
      this->numBinObjects = numBinObjects;
      binObjects = new TypedFeature<TClass1> * [numBinObjects];     
      accVal4fvals = vector<int>(this->filter->maxNumFilterValues());
      binObj4accVals = vector<int>(this->filter->maxNumFilterValues());
    }

    inline int numPhases() {
      return _numPhases;
    }
    
    virtual void computeHBlockLengths(int beg, int end, TStrand strand) = 0;

    // cache all those values in the feature vector that give a zero score
    // when computing the dot product with the parameter vector
    inline void updCacheParams(const FeatureVector ***params, int fConjInd = 0) {
      int i = 0, phase = 0;
      _featValues[0] = _featValues[1] = 0;
      fConjInd = fConjInd*domain*numPhases();
      const FeatureVector **nparams = params[this->paramInd()];
      
      for(phase = 0; phase < numPhases(); phase++) {
	int offset  = fConjInd + domain*phase;
	
	for(i = 0; i < domain; i++) {
	  TypedFeature<TClass1> *binObj = binObjects[binObj4accVals[i]];

	  if(binObj->isFeatureSet())
	    _cachedBinZeroValues[offset + i] = binObj->dotParamV(_featValues, nparams, offset + i);
	  else 
	    _cachedBinZeroValues[offset + i] = binObj->dotParamV(0.0, nparams, offset + i);
	}
      }     
    }

    inline double dotParamV(Tag *ge, int frame, 
                            const FeatureVector **params,
			    int fConjInd = 0,
			    TypedFilter<UCHAR> *filter = NULL) {
      
      assert(!filter);
      int beg = this->begTag(ge, frame);
      int end = this->endTag(ge, frame);
      int i = 0, phase = 0;

      computeHBlockLengths(beg, end, ge->getStrand());
      
      double result = 0;
      fConjInd = fConjInd*domain*numPhases();

      for( ; phase < numPhases(); phase++) {
	int offset  = fConjInd + domain*phase;
	
	for(i = 0; i < domain; i++) {

	  if(!hblockLengths[phase][i]) {	    
	    if(countZeroLenBlocks)
 	      result += _cachedBinZeroValues[offset + i];
	    continue;
	  }

	  TypedFeature<TClass1> *binObj = binObjects[binObj4accVals[i]];

	  if(binObj->isFeatureSet()) {
	    _featValues[1] = hblockLengths[phase][i];
	    _featValues[0] = _featValues[1]*100.0/ge->getLen();
	    result += binObj->dotParamV(_featValues, params, offset + i); 
	  }
	  else
	    result += binObj->dotParamV(hblockLengths[phase][i], params, offset + i);
	  
	}
      }
 
      return result;
    }
  
    inline void updParamV(double updVal, Tag *ge, int frame, 
			  FeatureVector **params, int fConjInd = 0,
			  TypedFilter<UCHAR> *filter = NULL) {
      
      updated = true;
      assert(!filter);
      int beg = this->begTag(ge, frame);
      int end = this->endTag(ge, frame);
      int i, phase = 0;
            
      computeHBlockLengths(beg, end, ge->getStrand());
      
      fConjInd = fConjInd*domain*numPhases();
      //      cerr << this->getName() << " " << ge->getPos() << " " << ge->getPos() + ge->getLen() - 1 << " " << numPhases() << endl;
      
      for( ; phase < numPhases(); phase++) {
	int offset  = fConjInd + domain*phase;
	for(i = 0; i < domain; i++) {
	  
	  if(!hblockLengths[phase][i] && !countZeroLenBlocks)
	    continue;
	  
	  TypedFeature<TClass1> *binObj = binObjects[binObj4accVals[i]];
	  
	  //	  cerr << "\t" << offset + i << " " << hblockLengths[phase][i] << endl;
	  if(binObj->isFeatureSet()) {
	    _featValues[1] = hblockLengths[phase][i];
	    _featValues[0] = _featValues[1]*100/ge->getLen();	    
	    binObj->updParamV(updVal, _featValues, params, offset + i);
	  }
	  else
	    binObj->updParamV(updVal, hblockLengths[phase][i], params, offset + i);
	}
      }
    }
    
    virtual ~FT_BinnedSegment() {
      delete [] binObjects;
      binObjects = NULL;
      delete [] hblockLengths;
      hblockLengths = NULL;
    }
    
  };
    
  
  template <class TClass1, class TClass2> 
  class FT_DenseBinnedSegment : public FT_BinnedSegment<TClass1, TClass2> {
   public:  
    FT_DenseBinnedSegment(
			  int fInd, 
			  int paramInd, 
			  int parsingFrames,
			  char *name, 
			  FilterEngine *fe, 
			  TypedFilter<TClass2> *filter, 
			  int period,
			  int numBinObjects,
			  bool countZeroLenBlocks)
      : FT_BinnedSegment<TClass1, TClass2>(fInd, paramInd, 
					   parsingFrames, name, 
					   fe, filter, period, 
					   numBinObjects,
					   countZeroLenBlocks, false) {
      
    }

    FT_DenseBinnedSegment(
			  int fInd, 
			  int paramInd, 
			  vector<std::string> & fargs, 
			  int & offset, 
			  FilterEngine *fe, 
			  FeatureEngine *fte
			  )
      : FT_BinnedSegment<TClass1, TClass2>(fInd, paramInd,
					   fargs, offset,
					   false, fe, 0) {
      
    }
      
    void computeHBlockLengths(int beg, int end, TStrand strand) {
      int phase = 0, i;
      
      for( ; phase < this->numPhases(); phase++) 
	for(i = 0; i < this->hblockLengths[phase].size(); i++)
	  this->hblockLengths[phase][i] = 0;
      
      for(i = beg; i <= end; i++) {
	TClass2 fval = this->filter->value(i, strand);
	int accval = this->accVal4fvals[fval];

	if(!this->subFeatIsOn(accval))
	  continue;

	phase = (i - beg) % this->numPhases();
	this->hblockLengths[phase][accval]++;
	
      }
    }
    
  };


  template <class TClass1, class TClass2, class TClass3> 
    class FT_SparseBinnedSegment : 
  public FT_BinnedSegment<TClass2, TClass3> {
   public:  
    FT_SparseBinnedSegment(
			   int fInd, 
			   int paramInd, 
			   int parsingFrames,
			   char *name, 
			   FilterEngine *fe, 
			   TypedFilter<TClass3> *filter, 
			   int period,
			   int numBinObjects,
			   bool countZeroLenBlocks)
      : FT_BinnedSegment<TClass2, TClass3>(fInd, paramInd, 
					   parsingFrames, name, 
					   fe, filter, period, 
					   numBinObjects,
					   countZeroLenBlocks, false) {
      
    }
    
    FT_SparseBinnedSegment(
			   int fInd, 
			   int paramInd, 
			   vector<std::string> & fargs, 
			   int & offset, 
			   FilterEngine *fe, 
			   FeatureEngine *fte)
      : FT_BinnedSegment<TClass2, TClass3>(fInd, paramInd, 
					   fargs, offset, 
					   false, fe, fte) {
      
    }

    void computeHBlockLengths(int beg, int end, TStrand strand) {
      FL_SparseGram<TClass1, TClass3> *cgram =
	(FL_SparseGram<TClass1, TClass3> *)this->filter;
      typename vector<HBlock<TClass3> >::iterator block = cgram->findHBlock(beg, strand);

      int phase = 0;
 
      for( ; phase < this->numPhases(); phase++) 
	for(int i = 0; i < this->hblockLengths[phase].size(); i++)
	  this->hblockLengths[phase][i] = 0;

      for( ; block != cgram->end(strand); ++block) {

	if(end < block->min)
	  break;

	int accval = this->accVal4fvals[block->val];

	if(!this->subFeatIsOn(accval))
	  continue;

	int bl = block->length, rl = block->max - block->min + 1;
	int ibeg = (beg > block->min ? beg : block->min);
	int iend = (end < block->max ? end : block->max);
	
	if(iend >= ibeg) {
	  int il = iend - ibeg + 1;

	  assert(il <= rl);

	  if(il < rl) 
	    bl = bl*il/rl;

	  phase = (ibeg - beg) % this->numPhases();
	  this->hblockLengths[phase][accval] += bl;

	}
      }
    }

    ~FT_SparseBinnedSegment() {}
    
  };


  template <class TClass1, class TClass2, class TClass3> 
    class FT_PhasedSparseBinnedSegment : 
  public FT_BinnedSegment<TClass2, TClass3> {
   protected:

   public:  
    FT_PhasedSparseBinnedSegment(
				 int fInd, 
				 int paramInd, 
				 int parsingFrames,
				 char *name, 
				 FilterEngine *fe, 
				 TypedFilter<TClass3> *filter, 
				 int period,
				 int numBinObjects,
				 bool countZeroLenBlocks)
      : FT_BinnedSegment<TClass2, TClass3>(fInd, paramInd, 
					   parsingFrames, name, 
					   fe, filter, period, 
					   numBinObjects,
					   countZeroLenBlocks, true) {
      

    }
    
    FT_PhasedSparseBinnedSegment(
				 int fInd, 
				 int paramInd, 
				 vector<std::string> & fargs, 
				 int & offset, 
				 FilterEngine *fe, 
				 FeatureEngine *fte)
      : FT_BinnedSegment<TClass2, TClass3>(fInd, paramInd, 
					   fargs, offset, 
					   true, fe, fte) {
      
    }

    void computeHBlockLengths(int beg, int end, TStrand strand) {
      FL_SparseGram<TClass1, TClass3> *cgram =
	(FL_SparseGram<TClass1, TClass3> *)this->filter;
      typename vector<HBlock<TClass3> >::iterator block = cgram->findHBlock(beg, strand);
      
      int phase = 0;

      for( ; phase < this->numPhases(); phase++) 
	for(int i = 0; i < this->hblockLengths[phase].size(); i++)
	  this->hblockLengths[phase][i] = 0;

      for( ; block != cgram->end(strand); ++block) {

	if(end < block->min)
	  break;

	int accval = this->accVal4fvals[block->val];

	if(!this->subFeatIsOn(accval))
	  continue;

	int bl = block->length, rl = block->max - block->min + 1;
	int ibeg = (beg > block->min ? beg : block->min);
	int iend = (end < block->max ? end : block->max);
	
	if(iend >= ibeg) {
	  int il = iend - ibeg + 1;

	  assert(il <= rl);

	  if(il < rl) 
	    bl = bl*il/rl;

	  phase = (ibeg - beg) % this->numPhases();
	  this->hblockLengths[phase][accval] += bl;
	}
      }
    }

    ~FT_PhasedSparseBinnedSegment() {}
    
  };
}

#endif
