/****************************************************************************
* FL_Cluter.h - part of the lless namespace, a general purpose
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
****************************************************************************/

#ifndef _FILTER_CLUSTER_H_
#define _FILTER_CLUSTER_H_

#include "FL_Score.h"
#include <sstream>
#include "FilterEngine.h"

namespace lless {
  class ResourceEngine;
  
  /**************************************************************************
   * FL_Cluster
   * is a subclass of TypedFilter<int> which computes a cluster, using 
   * the centroids as anchor points
   *
   **************************************************************************/
  template <class TClass> 
  class FL_Cluster: public TypedFilter<UCHAR> {
  private:
    int sldwSize;
    FL_Score<TClass> *filter;
    vector<TClass> centroids;
  public:
    //! Default constructor
    FL_Cluster(int fInd,                //!< A unique identifier.
               std::string & name,      //!< A unique name.
               int sldwSize          /*!< A sliding window's size.
                                          If equal to zero, the window size
                                          is set automatically to the length
                                          of the current input sequence.
                                         */
	       FL_Score<TClass> *filter,
               vector<TClass> centroids, //!< Context level ranges.
               ) : 
    TypedFilter<UCHAR>(fInd,name,
		       centroids.size(), FT_UCHAR) {

      this->sldwSize = sldwSize;
      this->filter = filter;
      this->centroids = centroids;

    }

    /**
     * Constructor from a Header string definition
     */
    FL_Cluster(int fInd,           //!< A unique identifier.
               /*! The Header string definition, loaded as a vector of strings.
                * The Header has the following form:\n\n
                * Filter name Cluster num_centroids 
                * centroid_1 [centroid_2 .. ] filter sldwindow log_scale\n\n
                */
               vector<std::string> &params, 
               int & offset,       //!< The index for vector params.
               ResourceEngine *re, //!< A pointer to the ResourceEngine object.
               FilterEngine *fe    //!< A pointer to the FilterEngine object.
               ) : 
      TypedFilter<UCHAR>(fInd, params, offset, 1, FT_UCHAR) {
	
      if(!sscanf(params[offset++].c_str(), "%d", &sldwSize))
        assert(0);
      
      this->filter = (FL_Score<TClass> *)fe->getFilter(params[offset++]);
      
      if(!sscanf(params[offset++].c_str(), "%d", &maxFilterVals))
        assert(0);
      
      double centroid;
      for(int i = 0; i < maxFilterVals; i++) {
        if(!sscanf(params[offset++].c_str(), "%lf", &centroid))
          assert(0);
        centroids.push_back(centroid);
      }
    }

    /**
     * A member function that extends the functionality of the original
     * function definedin FL_Score to change the contained gram object's
     * default strand to strand
     * @param strand the new default strand
     */
    inline void setDefaultStrand(TStrand strand) {
      this->defaultStrand = strand;
      filter->setDefaultStrand(strand);
    }
    
    inline void setSlidingWindow(int window) {
      sldwSize = window;
    }
    
    /**
     * A member function for computing the context level.
     * @param ctxContent the relative content(abundance).
     * @return The computed context level.
     */
    UCHAR closestCentroid(double score) {
      int last = 0;
      
      for(unsigned int j = 1; j < centroids.size(); j++) {
	TClass threshold = (centroids[j] + centroids[last])/2.0;
        if(score <= threshold)
          break;
	
        last = j;
      }
      return (UCHAR)last;
    }
    
    /**
     * A member function for computing filter values.
     * 
     * The value at each position of the input sequence is the centroid 
     * closest to the smoothed value in a sliding window centered at said
     * position.
     *
     * @param seq the input sequence.
     * @param filterVals the array of filter values.
     * @param len the length to consider. It could be shorter than the
     * actual length of the input sequence seq.
     */
    inline void computeVals(char *seq, UCHAR *filterVals, int len) {
      assert(len);
      int i = 0;
      double score = 0;
      int wSize = sldwSize;
      
      if(!wSize || len < wSize)
        wSize = len;
      
      // 5 prime
      for(i = 0; i < wSize; i++)
	score += filter[i + 1];
      
      UCHAR centroid = closestCentroid(score/wSize);
      
      for(i = -5; i < wSize/2; i++)
        filterVals[i] = centroid;
      
      // middle
      for( ; i < len - wSize/2; i++) {
	score -= filter[i - wSize/2 + 1];
	score += filter[i + wSize/2 + 1];
        centroid = closestCentroid(score/wSize);
        filterVals[i] = centroid;
      }
      
      // 3 prime
      for( ; i <= len + 5; i++)
        filterVals[i] = centroid;
    }
    
    ~FL_Cluster() {
    }
  };
}

#endif
