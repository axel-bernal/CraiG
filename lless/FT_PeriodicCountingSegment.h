/****************************************************************************
* FT_PeriodicCountingSegment.h - part of the lless namespace, a general purpose
*                                linear semi-markov structure prediction 
*                                library
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

#ifndef F_PERIODICCOUNTING_SEGMENT_H
#define F_PERIODICCOUNTING_SEGMENT_H
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "FT_CountingSegment.h"
#include "Sequence.h"
#include "NGram.h"


namespace lless {  

  /**
   * FT_PeriodicCountingSegment is a subclass of FT_CountingSegment and 
   * computes its feature value as sets of filter values over segments.
   * The filter values are periodic.
   * This feature works very efficiently during decoding when the provided 
   * filter has previously accumulated its values along the sequence
   ***************************************************************************/
  template <class TClass>
  class FT_PeriodicCountingSegment : public FT_CountingSegment<TClass> {
  
   protected:
    int _phase, _period;
    bool normalize;
    double mean;

   public:
    FT_PeriodicCountingSegment(
                               int fInd, 
                               int paramInd,
                               int parsingFrames,
                               char *name,
                               FilterEngine *fe, 
                               TypedFilter<TClass> *filter, 
                               int period, 
                               int phase = 0,
                               int normalize = false,
                               double mean = 1.0
                               )
      : FT_CountingSegment<TClass>(fInd,
				   paramInd,
				   parsingFrames, 
				   name, 
				   fe,
				   filter) {
      
      this->_period = period;
      this->_phase = phase;
      this->normalize = normalize;
      this->mean = mean;

    }
  
    /**
     * Constructor of a PeriodicCountingSegment in the feature file
     *
     * Feature  id  PeriodicCountingSegment 3 3-gram 0 3 1 false 1
     *
     * The definition above creates a PeriodicCountingSegment feature, which
     * extends CountingSegment to deal with filters whose values are periodic.
     * This feature has 3 phases, the filter 3-gram has values which are 
     * 3-periodic and the starting phase of the feature is 1. The last two
     * parameters indicate that the feature should not normalize the countings
     * (false) and instead should use 1 as default value for the standard 
     * deviation of the counts distribution.
     *
     */

    FT_PeriodicCountingSegment(
                               int fInd, 
                               int paramInd, 
                               vector<std::string> & fargs, 
                               int & offset,
                               FilterEngine *fe,
                               FeatureEngine *fte
                               )
      : FT_CountingSegment<TClass>(fInd, 
				   paramInd,
				   fargs,
				   offset, 
				   fe) {
      
      if(!sscanf(fargs[offset++].c_str(), "%d", &_period))
        assert(0);
      if(!sscanf(fargs[offset++].c_str(), "%d", &_phase))
        assert(0);
      
      normalize = Utils::stringToBoolean(fargs[offset++]);
      if(!normalize) {
        if(!sscanf(fargs[offset++].c_str(), "%lf", &mean))
          assert(0);
      }

    }

    inline int phase() {
      return _phase;
    }
  
    inline int period() {
      return _period;
    }

    virtual inline int begTag(Tag *ge, int frame) {
      int beg = ge->getPos();
      return beg + period()*(frame != 0) - frame;
    }
    
    virtual inline int lenTag(Tag *ge, int frame) {
      return period()*(ge->getLen()/period());
    }
    
    virtual inline int endTag(Tag *ge, int frame) {
      int end = begTag(ge, frame) + lenTag(ge, frame);
      end -= period()*(ge->getLen() % period() == 0);

      if(end > this->fe->seqLength())  
	end -= period();

      return end;
    }

    virtual inline double scaleDotParamV(Tag *ge, double result) {
      return (normalize ? result/ge->getLen() :  result/mean);
    }    

    ~FT_PeriodicCountingSegment() {

    }
  
  };

  
  /**
   * FT_PeriodicCountingSegmentWUpperLimit is a subclass of
   * FT_PeriodicCountingSegment and 
   * computes its feature value as sets of filter values over segments.
   * The filter values are periodic and the segment lengths are limited 
   * to be no longer than limit.
   * This feature works very efficiently during decoding when the provided 
   * filter has previously accumulated its values along the sequence
   ***************************************************************************/

  template <class TClass>
  class FT_PeriodicCountingSegmentWUpperLimit : public FT_PeriodicCountingSegment<TClass> {
  
   protected:
    int limit;
  
   public:
    FT_PeriodicCountingSegmentWUpperLimit(
                                          int fInd, 
                                          int paramInd, 
                                          int dotInd, 
                                          char *name, 
                                          FilterEngine *fe, 
                                          TypedFilter<TClass> *filter, 
                                          int period, 
                                          int phase = 0, 
                                          bool normalize = false,
                                          double mean = 1.0,
                                          int limit = 1
                                          ) 
      : FT_PeriodicCountingSegment<TClass>(fInd, paramInd, 
					   dotInd, name,
					   fe, filter, 
					   period, phase, 
					   normalize, mean) {
      
      this->limit = limit;  
    }
    
    FT_PeriodicCountingSegmentWUpperLimit(
                                          int fInd, 
                                          int paramInd, 
                                          vector<std::string> & fargs, 
                                          int & offset, 
                                          FilterEngine *fe, 
                                          FeatureEngine *fte
                                          )
      : FT_PeriodicCountingSegment<TClass>(fInd, 
					   paramInd, 
					   fargs,
					   offset, 
					   fe,
					   fte) {
      
      if(!sscanf(fargs[offset++].c_str(), "%d", &limit))
        assert(0);
  
    }
  
    inline int lenTag(Tag *ge, int frame) {
      int len = (ge->getLen() > limit)
	? limit
        : ge->getLen();

      return this->period()*(len/this->period());
    }

    inline int endTag(Tag *ge, int frame) {
      int len = (ge->getLen() > limit)
	? limit
        : ge->getLen();

      int end = this->begTag(ge, frame) + this->period()*(len/this->period());
      end -= this->period()*(len % this->period() == 0);
      
      //      int end = ge->getPos() + ((ge->getLen() > (limit + 1))
      //        ? (limit + 1)
      //        : ge->getLen()) - 1;

      //      end = beg + this->period()*((end - beg)/this->period());
      
      if(end > this->fe->seqLength())  
	end -= this->period();
      //      if(end < beg) end = beg;

      return end;

    }

    ~FT_PeriodicCountingSegmentWUpperLimit() {
  
    }
    
  };


  /**
   * The FT_AggregatedPeriodicCountingSegment is a subclass of
   * FT_PeriodicCountingSegment which
   * computes its feature value as sets of filter values over segments.
   * The main difference with FT_PeriodicCountingSegment is that in this class
   * there are only two filter value types: nonAggrFVal which is any 
   * enabled  value from the original filter and aggrFVal, to which 
   * any other value of the original filter is assigned to.
   * This feature works very efficiently during decoding when the provided 
   * filter has previously accumulated its values along the sequence
   ***************************************************************************/

  template <class TClass>
  class FT_AggregatedPeriodicCountingSegment : public FT_PeriodicCountingSegment<TClass> {
    
   protected:
    TClass aggrFVal, nonAggrFVal;

   public:
    FT_AggregatedPeriodicCountingSegment(
                                    int fInd, 
                                    int paramInd,
                                    int parsingFrames,
                                    char *name,
                                    FilterEngine *fe, 
                                    TypedFilter<TClass> *filter, 
                                    int period, 
                                    int phase,
                                    int normalize,
                                    double mean,
                                    TClass nonAggrFVal
                                    ) 
      : FT_PeriodicCountingSegment<TClass>(fInd,
					   paramInd, 
					   parsingFrames, 
					   name, 
					   fe,
					   filter,
					   period,
					   phase,
					   normalize,
					   mean
					   ) {
      
      this->nonAggrFVal = nonAggrFVal;
      findAggrIndex();
      
    }
    
    FT_AggregatedPeriodicCountingSegment(
                                    int fInd, 
                                    int paramInd, 
                                    vector<std::string> & fargs, 
                                    int & offset,
                                    FilterEngine *fe, 
                                    FeatureEngine *fte
                                    )
      : FT_PeriodicCountingSegment<TClass>(fInd, 
					   paramInd,
					   fargs,
					   offset, 
					   fe, 
					   fte) {
      
      if(!sscanf(fargs[offset++].c_str(), "%d", &nonAggrFVal))
        assert(0);
      findAggrIndex();
      
    }

    /**
     * Finds a valid index which is an enabled filter value
     * to which all other filter values except the one given as nonAggregated,
     * will aggregate, i.e. will be assigned to.
     */
    inline void findAggrIndex() {
      aggrFVal = 0;
      for(; aggrFVal < this->maxNumFeatValues(); aggrFVal++) {
        if(aggrFVal == nonAggrFVal)
          continue;
        break;
      }
    }

    /**
     * @return either aggrFVal or nonAggrFVal or the original filter value
     * if it has not been enabled.
     */
    inline double posfValue(int pos, TStrand strand) {
      TClass fVal = this->filter->value(pos, strand);
      if(fVal == nonAggrFVal)
        return fVal;

      return aggrFVal;
    }

    ~FT_AggregatedPeriodicCountingSegment() {

    }

  };

}

#endif
