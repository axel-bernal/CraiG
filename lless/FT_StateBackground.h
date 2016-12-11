/****************************************************************************
* FT_StateBackground.h - part of the lless namespace, a general purpose
*                        linear semi-markov structure prediction library
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

#ifndef _FT_STATE_BACKGROUND_H_
#define _FT_STATE_BACKGROUND_H_

#include "Feature.h"

namespace lless {

  class FeatureEngine;
  class FeatureEngine;


  /**
   * The FT_StateBackground class is a feature which value is proportional
   * to the length of the state. Generally scaled by the mean length of the 
   * state.
   ***************************************************************************/
  
  class FT_StateBackground: public TypedFeature<double> {
   private:
    double meanLength;
   public:
    FT_StateBackground(
                       int fInd, 
                       int paramInd,
                       int numFrames, 
                       char *name, 
                       double meanLength,
                       FilterEngine *fe
                       ) 
      : TypedFeature<double>(fInd, paramInd,
                             numFrames, name, 
                             fe, 1, FT_DOUBLE) {

      this->meanLength = meanLength;

    }
    
    FT_StateBackground(
                       int fInd, 
                       int paramInd, 
                       vector<std::string> & fargs, 
                       int & offset, 
                       FilterEngine *fe, 
                       FeatureEngine *fte
                       ) 
      : TypedFeature<double>(fInd, paramInd,
                             fargs, offset, 
                             fe, 1, FT_DOUBLE) {

      if(!sscanf(fargs[offset++].c_str(), "%lf", &meanLength))
        assert(0);      
      
    }
    
    inline double meanLen() {
      return meanLength;
    }

    inline double featValue(Tag *ge, int frame) {
      return ge->getLen()/meanLength;
    }
    
    inline double dotParamV(double featVal, 
			    const FeatureVector **params,
			    int fConjInd = 0) {
 
      double res = featVal*(*params[fConjInd])[0];
      //      cerr << featVal << " (" << res << ")" << endl;
      return res;
    }
    
    inline void updParamV(double updVal, double featVal, 
                         FeatureVector **params, int fConjInd = 0) {
      //    cerr << updVal << " " << this->getName() << "[" << fConjInd << "] = " << (*params[fConjInd])[0] << " -> "; 
      (*params[fConjInd])[0] += featVal*updVal;
      //    cerr << (*params[fConjInd])[0] << endl;
    }  

    ~FT_StateBackground() {
    }
    
  };
  
}

#endif
