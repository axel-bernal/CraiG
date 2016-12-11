/****************************************************************************
* FL_ChangePtSignal.h - part of the lless namespace, a general purpose
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

#ifndef _FL_CHANGE_PTSIGNAL_H_
#define _FL_CHANGE_PTSIGNAL_H_

#include "Utils.h"
#include "Filter.h"
#include "ChangePtUtils.h"
  
namespace lless {

  /**
   *
   * FL_ChangePoint is a subclass of Filter which computes the 0th order
   * statistic dynamically and very efficiently over intervals. Other order
   * statistics are not implemented. This filter accumulates
   *
   ***************************************************************************/

  template <class TClass>
    class FL_ChangePtSignal : public TypedFilter<UCHAR> {
   protected:
    string symbols;
    FL_ChangePoint<TClass> *filter;
   public:
    FL_ChangePtSignal(int fInd,       //!< A unique identifier.
		      std::string & name,       //!< A unique name.
		      std::string &symbols,      //!< symbols 4 representing cps
		      FL_ChangePoint<TClass> *filter  //!< filter to detect cps
		      ) : 
    TypedFilter<UCHAR>(fInd, 
		       name,
		       1, FT_UCHAR) {
      
      this->symbols = symbols;
      this->filter = filter;
    }
    
    /**
     * First constructor from a Header string definition
     */    
    FL_ChangePtSignal(
		      int fInd,         //!< A unique identifier.
		      vector<std::string> &params, //!< Header string definition
		      int & offset,       //!< The index for vector params.
		      ResourceEngine *re, /*!< A pointer to the ResourceEngine
					    object. */
		      FilterEngine *fe    /*!< A pointer to the FilterEngine 
					    object. */
		      ) 
      : TypedFilter<UCHAR>(fInd, 
			   params, 
			   offset, 
			   1, FT_UCHAR) {
      
      this->symbols = params[offset++];
      this->filter = (FL_ChangePoint<TClass> *)fe->getFilter(params[offset++]);

    }
    
    inline void computeVals(char *seq, char *vals, int len) {
      
      pair<TClass, TClass> &range = filter->range();
      vector<ChangePoint> & pchanges = 
	filter->getChangePoints(this->defaultStrand);

      ChangePtUtils::onsets2Filter(pchanges,
				   symbols[this->defaultStrand*2],
				   vals,
				   range.first, 2*range.second/3,
				   len, this->defaultStrand);
      ChangePtUtils::offsets2Filter(pchanges, 
				    symbols[this->defaultStrand*2 + 1],
				    vals,
				    range.first, 2*range.second/3,
				    len, this->defaultStrand);
      
    }

    ~FL_ChangePtSignal() {

    }
  };
}

#endif
