/****************************************************************************
* FT_EvdAlignment.h - part of the lless namespace, a general purpose
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

#ifndef FT_EVD_ALIGNMENT_FEAT_H
#define FT_EVD_ALIGNMENT_FEAT_H

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Sequence.h"
#include "Feature.h"
#include "Filter.h"
#include "FilterEngine.h"
#include "FeatureEngine.h"
#include "FT_Segment.h"
#include "FL_EvdEdgeAligner.h"
#include "FSM.h"

namespace lless {


  /**
   * FT_EvdEvdAlignmentType computes a feature for the current tag whose value 
   * is 0 if the current tag is not aligned with any evidence source's eddes,
   * 1 if either 5 or 3 prime are aligned and 2 if both boundaries are 
   * aligned
   */
  class FT_EvdAlignmentType : public TypedFeature<UCHAR> {
   protected:
    FL_EvdEdgeAligner *filter;
    FSM *fsm;
   public:

    FT_EvdAlignmentType(
			int fInd, 
			int paramInd, 
			vector<std::string> & fargs, 
			int & offset, 
			FilterEngine *fe, 
			FeatureEngine *fte
			)
      : TypedFeature<UCHAR>(fInd, paramInd,
			    fargs, offset, 
			    fe, 4, FT_UCHAR) {
      
      this->filter = (FL_EvdEdgeAligner *)fe->getFilter(fargs[offset++]);
      this->fsm = fte->getFSM();
    }
    
    inline double featValue(Tag *tag, int frame) {
      Node *node = fsm->node(tag->getParseType());
      int type = filter->alignmentType(node, tag, frame);
      //      cerr << "type = " << type << "\n";
      return type;
    }
  };  

  /**
   * FT_EvdEvdAlignmentWeight computes the alignment weight for the
   * current tag whose value.
   */
  class FT_EvdAlignmentWeight : public TypedFeature<double> {
   protected:
    FL_EvdEdgeAligner *filter;
    FSM *fsm;
   public:

    FT_EvdAlignmentWeight(
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
      
      this->filter = (FL_EvdEdgeAligner *)fe->getFilter(fargs[offset++]);
      this->fsm = fte->getFSM();
    }
    
    inline double featValue(Tag *tag, int frame) {
      Node *node = fsm->node(tag->getParseType());
      double weight = filter->alignmentWeight(node, tag, frame);
      //      cerr << "weight = " << weight << "\n";
      return weight;
    }
  };  
  
}

#endif
