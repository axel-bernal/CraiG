/****************************************************************************
* FT_ZScore.h - part of the lless namespace, a general purpose
*               linear semi-markov structure prediction library
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

#ifndef _FT_ZSCORE_H_
#define _FT_ZSCORE_H_
  
#include "Feature.h"
#include "FeatureEngine.h"


namespace lless {

  /**
   * FT_ZScore is a subclass of TypedFeature<double> whose value is the Z 
   * score of a distribution induced by the scores of the contained feature 
   * ft scores over Tag objects.
   **************************************************************************/

  class FT_ZScore : public TypedFeature<double> {
  
   private:
    double mean;
    double std;
    TypedFeature<double> *ft;
  
   public:
    FT_ZScore(
           int fInd, 
           int paramInd,
           char *name,
           FilterEngine *fe,
           TypedFeature<double> *ft,
           double mean,
           double std
           ) 
      : TypedFeature<double>(fInd, paramInd, 
                             ft->parsingFrames(), name,
                             fe, 1, FT_DOUBLE) {
  
    }
    
    FT_ZScore(
           int fInd, 
           int paramInd, 
           vector<std::string> & fargs, 
           int & offset, 
           FilterEngine *fe, 
           FeatureEngine *fte
           ) 
      : TypedFeature<double>(fInd, paramInd,
                             fte->getFeature(fargs[offset + 2])->parsingFrames(),
                             fargs, offset, 
                             fe, 1, FT_DOUBLE) {
      
      this->ft = (TypedFeature<double> *)fte->getFeature(fargs[offset++]);
      if(!sscanf(fargs[offset++].c_str(), "%lf", &mean))
        assert(0);
      if(!sscanf(fargs[offset++].c_str(), "%lf", &std))
        assert(0);
    }

    inline double featValue(Tag *ge, int frame) {
      return (ft->featValue(ge, frame) - mean)/std;
    }
  
    inline double dotParamV(double featVal, 
			    const FeatureVector **params,
			    int fConjInd = 0) { 
  
      return featVal*(*params[fConjInd])[0];
    }
    
    inline void updParamV(double updVal, double featVal, 
			  FeatureVector **params, int fConjInd = 0) {
  
      (*params[fConjInd])[0] += featVal*updVal;
    }  

    ~FT_ZScore() {
  
    }  

  };
    
}
  
#endif
