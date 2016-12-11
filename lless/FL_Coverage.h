/****************************************************************************
* FL_Coverage.h - part of the lless namespace, a general purpose
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

#ifndef _FILTER_COVERAGE_H_
#define _FILTER_COVERAGE_H_

#include "Utils.h"
#include "FL_Score.h"
#include "ContextIMM.h"
#include "FilterEngine.h"
#include "ResourceEngine.h"
#include "IMM.h"
#include "Motif.h"
  
namespace lless {
  
  /**
   *
   * FL_Coverage is a subclass of FL_Score. The scoring function matches
   * a given motif profile against each position of the sequence. The score
   * reflects the strength of the match.
   */
  template<class TClass>
    class FL_Coverage : public FL_Score<TClass> {
    
   private:
    TypedFilter<TClass> *cov;
    FL_EvdEdgeAligner *edges;

   public:
    //! Default constructor
    FL_Coverage(int fInd,             //!< A unique identifier.
		std::string & name,   //!< A unique name.
		int period,               //!< Filter values period.
		int smoothWindow,    //!< Window size for smoothing scores 
		
		TypedFilter<TClass> *cov,    //!< Coverage information
		FL_EvdEdgeAligner *edges    //!< Edge information
		)
      :  FL_Score<TClass>(fInd,
                          name, 
			  period,
			  smoothWindow,
                          cov->valType()) {

      this->cov = cov;
      this->edges = edges;
      
    }
  
    /**
     * Constructor from a Header string definition
     */
    FL_Coverage(int fInd,            //!< A unique identifier.
                  /*! The Header string definition, loaded as a vector of
                   * strings.
                   * The Header has the following form:\n\n
                   * Filter name Coverage arrInd period smoothenWSize
                   * threshold motif\n\n
                   * The description of the fields could be found in the 
                   * other constructor(s)
                   */
                  vector<std::string> & params, 
                  int & offset,       //!< The index for vector params.
                  ResourceEngine *re, /*!< A pointer to the ResourceEngine 
                                        object. */
                  FilterEngine *fe    /*!< A pointer to the FilterEngine 
                                        object. */
                  )
      : FL_Score<TClass>(fInd, 
                         params,
                         offset, 
                         fe) {

      this->cov = (TypedFilter<TClass> *)fe->getFilter(params[offset++]);
      assert(this->cov);
      this->edges = (FL_EvdEdgeAligner *)fe->getFilter(params[offset++]);
      assert(this->edges);
    }
    
    /**
     * A member function that extends the functionality of the original
     * function definedin FL_Score to change the contained gram object's
     * default strand to strand
     * @param strand the new default strand
     */
    inline void setDefaultStrand(TStrand strand) {
      this->defaultStrand = strand;
      cov->setDefaultStrand(strand);
      edges->setDefaultStrand(strand);
    }

    /**
     * The implemented scoring function. At each position p the score is the
     * match strength of the fragmeng seq[p..p+length(motif)]  against the 
     * motif profile, if the score is above threshold, or zero otherwise.
     * @param seq the input sequence.
     * @param filterVals the array of filter values.
     * @param len the length to consider. It could be shorter than the
     * actual length of the input sequence seq.
     */
    void scoreInputSeq(char *seq, TClass *filterVals, int len) {
      int leftB = 1;
      int rightB;
      filterVals[0] = 0;
      for( ; leftB <= len; leftB++)
	filterVals[leftB] = cov->value(leftB);
      
      for(leftB = 1; leftB <= len; leftB++) {
	EvdEdges &eds = edges->value(leftB);
	vector<EvdEdge> &fwdEdges = eds.allFwdEdges();
	
	if(!fwdEdges.size())  continue;
	
	for(int j = 0; j < fwdEdges.size(); j++) {
	  int k;
	  rightB = fwdEdges[j].boundary;
	  int weight = fwdEdges[j].weight;
	  for(k = leftB; k <= rightB; k++) 
	    filterVals[k] += weight;
	}	  
      }

    }  

    ~FL_Coverage() {
      
    }

  };

}

#endif
