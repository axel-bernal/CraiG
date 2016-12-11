/****************************************************************************
* FL_Score.h - part of the lless namespace, a general purpose
*              linear semi-markov structure prediction library
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

#ifndef _FILTER_SCORE_H_
#define _FILTER_SCORE_H_
  
#include "Utils.h"
#include "Filter.h"
#include "ContextIMM.h"
#include "FilterEngine.h"
#include "ResourceEngine.h"
#include "IMM.h"
#include "Motif.h"
  
namespace lless {

  /**
   *
   * FL_Score is a subclass of Filter which uses a scoring function to
   * compute a score at each position of the input sequence. The
   * details of the implementation of the scoring function are left to 
   * each derived class. This is an abstract base class for other 
   * FL_Score-derived classes.
   *
   * Filter values can accumulate along the sequence. This makes the
   * computation of segment scores very efficient.
   *
   ***************************************************************************/

  template <class TClass> class FL_Score : public TypedFilter<TClass> {
   protected:
    int smoothWindow;

   public:
    //! Default constructor
    FL_Score(
             int fInd,       //!< A unique identifier.
             std::string & name,       //!< A unique name.
             int period,               //!< Filter values period.
             int smoothWindow = 0,    //!< Window size for smoothing scores.
             TValType type = FT_DOUBLE //!< Filter value type.
            ) : 
      TypedFilter<TClass>(fInd, 
                          name,
                          1, FT_DOUBLE) {
  
      this->_period = period;
      this->smoothWindow = smoothWindow;

    }

    /**
     * First constructor from a Header string definition
     */    
    FL_Score(
             int fInd,         //!< A unique identifier.
             vector<std::string> &params, //!< Header string definition
             int & offset,       //!< The index for vector params.
             FilterEngine *fe    //!< A pointer to the FilterEngine object.
             ) 
      : TypedFilter<TClass>(fInd, 
                            params, 
                            offset, 
                            1) {
  
      if(!sscanf(params[offset++].c_str(), "%d", &this->_period))
        assert(0);
  
      if(!sscanf(params[offset++].c_str(), "%d", &smoothWindow))
        assert(0);

    }

    /**
     * Second constructor from a Header string definition.
     * Parameters defined as above.
     */      
    FL_Score(int fInd, 
             vector<std::string> &params, 
             int & offset, 
             FilterEngine *fe,
             TValType type
             ) 
      : TypedFilter<TClass>(fInd, 
                            params, 
                            offset, 
                            1, type) {
      
      if(!sscanf(params[offset++].c_str(), "%d", &this->_period))
        assert(0);
  
      if(!sscanf(params[offset++].c_str(), "%d", &smoothWindow))
        assert(0);

    }

    /**
     * A member function for computing filter values.
     *
     * It calls the implementation-dependent scoring function for computing
     * filter values, smoothens these values if smoothenWSize is not zero and
     * finally, it accumulates the values along the input sequence.
     * @param seq the input sequence.
     * @param filterVals the array of filter values.
     * @param len the length to consider. It could be shorter than the
     * actual length of the input sequence seq.
     */
    inline void computeVals(char *seq, TClass *filterVals, int len) {
      //      cerr << this->getName() << endl;
      scoreInputSeq(seq, filterVals, len);

      if(smoothWindow)
	Utils::smoothScores(filterVals, smoothWindow, 1, len);
    }
  
    virtual void scoreInputSeq(char *seq, TClass *filterVals, int len) = 0;
  
    virtual ~FL_Score() {
  
    }
  
  };

}
  
#endif
