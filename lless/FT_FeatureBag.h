/****************************************************************************
* FT_FeatureBag.h - part of the lless namespace, a general purpose
*                   linear semi-markov structure prediction library
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

#ifndef _FT_FEATURE_BAG_H_
#define _FT_FEATURE_BAG_H_
  
#include "Utils.h"
#include "Feature.h"
#include "FeatureEngine.h"
#include <stdio.h>
#include "ParamModel.h"

namespace lless {
  
  /**
   * The FT_FeatureBag class is a bag of features, previously trained, and
   * which represent an alternative view of the data. The class resources and
   * filters are assumed to have been loaded previously in main.
   ***************************************************************************/

  class FT_FeatureBag : public TypedFeature<double> {
  
   private:
    ParamModel *bag;
  
   public:
    FT_FeatureBag(
                  int fInd, 
                  int paramInd, 
                  int parsingFrames,
                  char *name, 
                  ParamModel *bag,
                  FilterEngine *fe,
                  FeatureEngine *fte
                  )
      : TypedFeature<double>(
                             fInd, paramInd,
                             parsingFrames, name, 
                             fe, 1, FT_DOUBLE
                             ) {

      this->bag = bag;
      bag->readParams(fte->getFSM(), fe);
      
    }
  
    FT_FeatureBag(
                   int fInd,
                   int paramInd, 
                   vector<std::string> & fargs, 
                   int & offset, 
                   FilterEngine *fe,
                   FeatureEngine *fte
                   )
      : TypedFeature<double>(
                             fInd, paramInd,
                             fargs, offset, 
                             fe, 1, FT_DOUBLE
                             ) {
  
      ResourceEngine *re = fe->getResourceEngine();
      bag = (ParamModel *)re->getResource(fargs[offset++]);
      bag->readParams(fte->getFSM(), fe);
      
    }
    

    inline  bool isFeatureBag() {
      return true;
    }


    inline double featValue(Tag *ge, int frame) {
      FeatureEngine *fte = bag->featureEngine();

      if(ge->getGEClass() == EDGE_INST)
        return fte->dotParamV4Edges((EdgeInst *)ge, frame);
      
      return fte->dotParamV4Nodes((NodeInst *)ge, frame);

    }

    inline void doPreComputation(double **array,
                                 int preCompArrayIndex, 
                                 int beg, int end, TStrand strand) {
      
      bag->preComputeFeatures(beg, end, strand);

    }

    inline void undoPreComputation() {
      
      bag->deletePreComputedFeatures();

    }

    inline double dotParamV(double featVal, 
                            const FeatureVector **params,
			    int fConjInd = 0) {
      
      return featVal*(*params[fConjInd])[0];
    }
    
    inline void updParamV(double updVal, double featVal, 
                          FeatureVector **params, int fConjInd = 0) {
      
      //      cerr << this->getName() << " " << featVal <<  endl;
      (*params[fConjInd])[0] += featVal*updVal;
    }  

    ~FT_FeatureBag() {
    }

  };

}
  
#endif
