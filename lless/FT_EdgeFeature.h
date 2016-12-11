/****************************************************************************
* FT_EdgeFeature.h - part of the lless namespace, a general purpose
*                    linear semi-markov structure prediction library
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

#ifndef _FT_EDGE_FEATURE_H_
#define _FT_EDGE_FEATURE_H_
#include "Feature.h"
#include "Filter.h"
#include "FT_Edge.h"
#include "FilterEngine.h"
#include "FeatureEngine.h"

namespace lless {

  /**
   * FT_EdgeFeature is a subclass of FT_Edge whose value is the phase
   * of the Tag object received as parameter.
   */
  template <class TClass> 
  class FT_EdgeFeature : public FT_Edge<TClass> {
   protected:
    TypedFeature<TClass> *f;
   public:
    /**
     * Default constructor
     * @param fInd A feature index.     
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * @param parsingFrames The maximum number of phases(frames) of any Tag
     * object to which this feature is tied to.
     * @param name A unique feature identifier.     
     * @param fe A pointer to a FilterEngine object.
     * @param maxNumFeatVals The number of possible feature values. 
     * If greater than zero then the Feature object's possible values are
     * countable.
     * @param signal the signal this feature is associated with.
     */
    FT_EdgeFeature(
		   int fInd,
		   int paramInd,
		   int parsingFrames,
		   char *name,
		   FilterEngine *fe,
		   FL_Signal<EdgeInst> *signal,
		   TypedFeature<TClass> *f
                 ) 
      : FT_Edge<TClass>(fInd,
		       paramInd, 
		       parsingFrames, 
		       name, fe, 
		       f->maxNumFeatValues(),
		       signal,
		       f->valType()) {
      this->f = f;
    }

    /**
     * Constructor from a Header string definition. 
     * @param fInd A feature index.     
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * @param fargs The Header string definition loaded as a vector of strings.
     * The Header has the following form:\n\n
     * Feature name EdgeFeature parsingFrames signal maxNumFeatVals\n\n
     * The description of the fields could be found in 
     * the other constructor(s)     
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param fte A pointer to a FeatureEngine object.
     */      
    FT_EdgeFeature(
		   int fInd, 
		   int paramInd, 
		   vector<std::string> & fargs, 
		   int & offset, 
		   FilterEngine *fe, 
		   FeatureEngine *fte
		   ) 
      : FT_Edge<TClass>(fInd, 
			paramInd,
			fargs, 
			offset, 
			fe,
			0) { 
      
      f = (TypedFeature<TClass> *)fte->getFeature(fargs[offset++]);
      this->_numParams = f->maxNumParams();
    }
    
    inline double dotParamV(Tag *ge, int frame, 
			    const FeatureVector **params,
			    int fConjInd = 0,
			    TypedFilter<UCHAR> *filter = NULL) {
      
      int offset = fConjInd;

      if(filter && !f->isFeatureSet()) {
	offset += filter->value(ge->getPos(), ge->getStrand());
	filter = NULL;
      }
      
      return f->dotParamV(ge, frame, params, offset, filter);
    }

    virtual inline void updParamV(double updVal, Tag *ge, int frame, 
				  FeatureVector **params,
				  int fConjInd = 0,
				  TypedFilter<UCHAR> *filter = NULL) {      

      int offset = fConjInd;

      if(filter && !f->isFeatureSet()) {
	offset += filter->value(ge->getPos(), ge->getStrand());
	filter = NULL;
      }
      //      cerr << "updating parameters for " << this->getName() << " " << updVal << endl;
      f->updParamV(updVal, ge, frame, params, offset, filter);
    }

    ~FT_EdgeFeature() {}

  };

}  

#endif
