/****************************************************************************
* FT_ScoringSegment.h - part of the lless namespace, a general purpose
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

#ifndef FT_SCORING_SEGMENT_FEAT_H
#define FT_SCORING_SEGMENT_FEAT_H

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Sequence.h"
#include "Feature.h"
#include "Filter.h"
#include "FilterEngine.h"
#include "FeatureEngine.h"
#include "FT_Segment.h"

namespace lless {

  /**
   * FT_ScoringSegment is a subclass of TypedFeature whose value of this 
   * feature (and its derivate classes) is the score obtained by accumulating
   * values of the contained filter over a tag which must be a segment
   * (node with variable length).
   * This feature works well with bimodal filter values, in which scores above
   * a certain threshold (the mean) should be classified differently.
   * Invariants: the contained filter must accumulate.
   *
   ***************************************************************************/
  
  template <class TClass> class FT_ScoringSegment : public FT_Segment<TClass> {
    
   protected:
    double mean;
    bool normalize;

   public:  
    FT_ScoringSegment(
                      int fInd, 
                      int paramInd, 
                      int parsingFrames,
                      char *name, 
                      FilterEngine *fe, 
                      TypedFilter<TClass> *filter, 
                      int period,
                      bool normalize,
                      double mean
                      )
      : FT_Segment<TClass>(fInd, paramInd, 
                           parsingFrames, name, 
                           fe, filter, period, 2) {
      
      this->normalize = normalize;
      this->mean = mean;
      assert(this->filter->period());
    }

    FT_ScoringSegment(
                      int fInd, 
                      int paramInd, 
                      vector<std::string> & fargs, 
                      int & offset, 
                      FilterEngine *fe,
                      bool normalize,
                      double mean)
      : FT_Segment<TClass>(fInd, paramInd,
                        fargs, offset, 
                        fe, 2) {

      this->normalize = normalize;
      this->mean = mean;
      assert(this->filter->period());
    }
    
    FT_ScoringSegment(
                      int fInd, 
                      int paramInd, 
                      int parsingFrames,
                      vector<std::string> & fargs, 
                      int & offset, 
                      FilterEngine *fe,
                      int period,
                      bool normalize,
                      double mean)
      : FT_Segment<TClass>(fInd, paramInd,
                           parsingFrames,
                           fargs, offset, 
                           fe, period, 2) {
      
      this->normalize = normalize;
      this->mean = mean;
      assert(this->filter->period());
    }

    FT_ScoringSegment(
                      int fInd, 
                      int paramInd, 
                      vector<std::string> & fargs, 
                      int & offset, 
                      FilterEngine *fe, 
                      FeatureEngine *fte
                      )
      : FT_Segment<TClass>(fInd, paramInd,
                        fargs, offset, 
                        fe, 2) {
      
      if(!sscanf(fargs[offset++].c_str(), "%lf", &mean))
        assert(0);

      normalize = Utils::stringToBoolean(fargs[offset++]);
      assert(this->filter->period());
    }

    virtual inline double featValue(Tag *ge, int frame) {
      int beg = ge->getPos(), end = beg + ge->getLen() - 1;
      Utils::fixWithinBounds(beg, end, this->fe->seqLength());
      
      if(end < beg)
        return DOUBLE_INFINITY;

      int len = end - beg + 1;

      beg = this->begTag(ge, frame);
      end = this->endTag(ge, frame);
      
      double segVal = this->filter->accValue(end, 0, ge->getStrand())
        - this->filter->accValue(beg - this->period(), 0, ge->getStrand());

      //      cerr << "(" << ge->getPos() << " " << ge->getLen() << " " << frame << ") ";
      //      cerr << this->getName() << " " << this->period() << " " << this->filter->accValue(beg - this->period(), 0, ge->getStrand()) << " " << end << " " << this->filter->accValue(end, 0, ge->getStrand()) << " " << len << " " << segVal << endl;
      return normalize ? segVal/len : segVal;

    }
    
    inline double dotParamV(double featVal, 
                            const FeatureVector **params, 
			    int fConjInd = 0) {

      if(featVal == DOUBLE_INFINITY)
        return 0;
      //      cerr << "return " << featVal << " " << (*params[fConjInd])[featVal > mean ? 1 : 0] << " " << mean << endl;
      return featVal*(*params[fConjInd])[featVal > mean ? 1 : 0];
    }
    
    inline void updParamV(double updVal, double featVal, 
                         FeatureVector **params, int fConjInd = 0) {
      if(featVal == DOUBLE_INFINITY)
        return;
      //      cerr << this->getName() << " " << featVal << endl;
      (*params[fConjInd])[featVal > mean ? 1 : 0] += featVal*updVal;
    }  

    ~FT_ScoringSegment() {
      
    }
    
  };

  /**
   * FT_WeightedScoringSegment is a subclass of TypedFeature whose value of this 
   * feature (and its derivate classes) is the score obtained by accumulating
   * values of the contained filter over a tag which must be a segment
   * (node with variable length).
   * This feature works well with bimodal filter values, in which scores above
   * a certain threshold (the mean) should be classified differently.
   * Invariants: the contained filter must accumulate.
   *
   ***************************************************************************/
  
  template <class TClass> class FT_WeightedScoringSegment : public FT_Segment<TClass> {
    
   protected:
    double weight;

   public:  
    FT_WeightedScoringSegment(
                      int fInd, 
                      int paramInd, 
                      int parsingFrames,
                      char *name, 
                      FilterEngine *fe, 
                      TypedFilter<TClass> *filter, 
                      int period,
                      double weight
                      )
      : FT_Segment<TClass>(fInd, paramInd, 
                           parsingFrames, name, 
                           fe, filter, period, 1) {
      
      this->weight = weight;
      assert(this->filter->period() && weight);
    }
    

    FT_WeightedScoringSegment(
                      int fInd, 
                      int paramInd, 
                      vector<std::string> & fargs, 
                      int & offset, 
                      FilterEngine *fe, 
                      FeatureEngine *fte
                      )
      : FT_Segment<TClass>(fInd, paramInd,
                        fargs, offset, 
                        fe, 1) {
      
      if(!sscanf(fargs[offset++].c_str(), "%lf", &weight))
        assert(0);

      assert(this->filter->period() && weight);
    }

    virtual inline double featValue(Tag *ge, int frame) {
      int beg = ge->getPos(), end = beg + ge->getLen() - 1;
      Utils::fixWithinBounds(beg, end, this->fe->seqLength());
      
      if(end < beg)
        return DOUBLE_INFINITY;

      beg = this->begTag(ge, frame);
      end = this->endTag(ge, frame);
      
      double segVal = this->filter->accValue(end, 0, ge->getStrand())
        - this->filter->accValue(beg - this->period(), 0, ge->getStrand());

      //      cerr << "(" << ge->getPos() << " " << ge->getLen() << " " << frame << ") ";
      //      cerr << this->getName() << " " << this->period() << " " << this->filter->accValue(beg - this->period(), 0, ge->getStrand()) << " " << end << " " << this->filter->accValue(end, 0, ge->getStrand()) << " " << len << " " << segVal << endl;
      return weight*segVal;
      
    }
    
    inline double dotParamV(double featVal, 
                            const FeatureVector **params, 
			    int fConjInd = 0) {

      if(featVal == DOUBLE_INFINITY)
        return 0;
      //      cerr << "return " << featVal << " " << (*params[fConjInd])[featVal > mean ? 1 : 0] << " " << mean << endl;
      return featVal*(*params[fConjInd])[0];
    }
    
    inline void updParamV(double updVal, double featVal, 
                         FeatureVector **params, int fConjInd = 0) {
      if(featVal == DOUBLE_INFINITY)
        return;
      //      cerr << this->getName() << " " << featVal << endl;
      (*params[fConjInd])[0] += featVal*updVal;
    }  

    ~FT_WeightedScoringSegment() {
      
    }
    
  };

  /**
   * FT_RatioSegment is a subclass of TypedFeature<double> whose value is
   * the difference (in exponential space) between the values of two
   * FT_SegmentScore features with parameters filter1, _period1 and filter2,
   * _period2 respectively.
   */
  class FT_RatioSegment : public FT_ScoringSegment<double> {
    
  protected:
    TypedFilter<double> *filter2;
    int _period2;
    
  public:    
    FT_RatioSegment(
                    int fInd, 
                    int paramInd, 
                    int parsingFrames, 
                    char *name, 
                    FilterEngine *fe, 
                    TypedFilter<double> *filter, 
                    int _period, 
                    TypedFilter<double> *filter2, 
                    int _period2
                    )
      : FT_ScoringSegment<double>(fInd, paramInd, 
                                  parsingFrames, name, fe,
                                  filter, _period, true, 0) { 
      
      this->filter2 = filter2;
      this->_period2 = _period2;
    }

    FT_RatioSegment(
                    int fInd, 
                    int paramInd, 
                    vector<std::string> & fargs, 
                    int & offset, 
                    FilterEngine *fe, 
                    FeatureEngine *fte
                    )
      : FT_ScoringSegment<double>(fInd, paramInd,
                                  fargs, offset,
                                  fe, true, 0) {
      
      this->filter2 = (TypedFilter<double> *)fe->getFilter(fargs[offset++]);
      if(!sscanf(fargs[offset++].c_str(), "%d", &_period2))
        assert(0);

    }

    inline double featValue(Tag *ge, int frame) {

      int beg = ge->getPos(), end = beg + ge->getLen() - 1;
      Utils::fixWithinBounds(beg, end, this->fe->seqLength());
      
      if(end < beg)
        return DOUBLE_INFINITY;

      int len = end - beg + 1;

      beg = this->begTag(ge, frame);
      end = this->endTag(ge, frame);

      double res1 = this->filter->accValue(end, 0, ge->getStrand())
        - this->filter->accValue(beg - this->period(), 0, ge->getStrand());
      double res2 = filter2->accValue(end, 0, ge->getStrand())
        - filter2->accValue(beg - this->period(), 0, ge->getStrand());
      //      cerr << this->getName() << " " << this->period() << " " << beg - this->period() << " " << res1 << " " << end << " " << res2 << " " << len << endl;

      return (res1 - res2)/len;
    }
    
    ~FT_RatioSegment() {
    }  
    
  };

  /**
   * FT_CodingDifferential is a subclass of TypedFeature<double> whose value
   * is the coding differential as defined in the literature, i.e. the 
   * in-frame score given by a FT_ScoringSegment feature minus the average of
   * the out-of-frame scores given by similar features.
   */
  class FT_CodingDifferential : public FT_ScoringSegment<double> {
   protected:
    int _normalizer;
   public:
    FT_CodingDifferential(
                          int fInd, 
                          int paramInd, 
                          char *name, 
                          FilterEngine *fe, 
                          TypedFilter<double> *filter,
                          int normalizer = 3
                          )
      : FT_ScoringSegment<double>(fInd, paramInd, 
                                  3, name, fe, 
                                  filter, 3, true, 0) {
      
      // Coding differential has _period 3 always
      _normalizer = normalizer;
    }
    
    FT_CodingDifferential(
                          int fInd, 
                          int paramInd, 
                          vector<std::string> & fargs, 
                          int & offset, 
                          FilterEngine *fe, 
                          FeatureEngine *fte
                          )
      : FT_ScoringSegment<double>(fInd, paramInd, 3,
                                  fargs, offset,
                                  fe, 3, true, 0) {

      _normalizer = 2;

      if(offset < fargs.size())
        if(!sscanf(fargs[offset++].c_str(), "%d", &_normalizer))
          assert(0);
      
    }
    
    inline double featValue(Tag *ge, int frame) {
      
      double results[] = {0, 0};

      int beg = ge->getPos(), end = beg + ge->getLen() - 1;
      Utils::fixWithinBounds(beg, end, this->fe->seqLength());
      
      if(end < beg)
        return DOUBLE_INFINITY;

      int len = end - beg + 1;
      TStrand strand = ge->getStrand();

      for(int k = 0; k < this->period(); k++) {
        beg = this->begTag(ge, k);
        end = this->endTag(ge, k);

        results[k == frame] += this->filter->accValue(end, 0, strand) -
          this->filter->accValue(beg - this->period(), 0, strand);

        //        cerr << this->getName() << " " << this->period() << " " << k << " " << beg - this->period() << " " << results[k == frame] << " " << beg << " " << end << " " << len << endl;      
        
      }
      return (results[1] - results[0]/_normalizer)/len; 
    }    

    ~FT_CodingDifferential() {
    }

  };
  
  /**
   * FT_ScoringAtPos computes its value as the score of the contained filter 
   */
  template <class TClass> class FT_ScoringAtPos : public TypedFeature<TClass> {
  protected:
    TypedFilter<TClass> *filter;
    double weight;
    double mean;
  public:  
    FT_ScoringAtPos(
                    int fInd,
                    int paramInd, 
                    int parsingFrames,
                    char *name,
                    FilterEngine *fe, 
                    TypedFilter<double> *filter, 
                    double weight,
                    double mean
                    )
      : TypedFeature<TClass>(fInd, paramInd, parsingFrames,
                             name, fe, 2) {
      
      this->filter = filter;
      this->weight = weight;
      this->mean = mean;
    }
    
    FT_ScoringAtPos(
             int fInd, 
             int paramInd, 
             vector<std::string> & fargs, 
             int & offset, 
             FilterEngine *fe, 
             FeatureEngine *fte
             )
      : TypedFeature<TClass>(fInd, paramInd,
                             fargs, offset, 
                             fe, 2) {
      
      this->filter = (TypedFilter<TClass> *)fe->getFilter(fargs[offset++]);
      
      if(!sscanf(fargs[offset++].c_str(), "%lf", &weight))
        assert(0);
      if(!sscanf(fargs[offset++].c_str(), "%lf", &mean))
        assert(0);

    }

    inline double featValue(Tag *ge, int frame) {
      
      int beg = ge->getPos();
      return weight*filter->value(beg, ge->getStrand());
    }
    
    inline double dotParamV(double featVal, 
                            const FeatureVector **params,
			    int fConjInd = 0) {
      
      return featVal*(*params[fConjInd])[featVal > mean ? 1 : 0];
    }

    inline void updParamV(double updVal, double featVal, 
                         FeatureVector **params, int fConjInd = 0) {
      
      //      cerr << this->getName() << " " << featVal << endl;
      (*params[fConjInd])[featVal > mean ? 1 : 0] += featVal*updVal;
    }  

    ~FT_ScoringAtPos() {
    }
    
  };
}

#endif
