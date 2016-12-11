/****************************************************************************
* LatticeEdge.h - part of the lless namespace, a general purpose
*             linear semi-markov structure prediction library
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

#ifndef _LATTICE_EDGE_H_
#define _LATTICE_EDGE_H_
#include "EdgeInst.h"
#include "NodeInst.h"
#include "FeatureEngine.h"
#include <algorithm>
#include "FSMQueue.h"
#include "FT_HiddenSeq.h"

namespace lless { 

  class LatticeVertex;

  /**
   * LatticeEdge is class which implements a back edge within a lattice.
   * It provides all the needed information to backtrack a path from the
   * the ending to the beginning node within the Lattice object.
   */
  class LatticeEdge {

   public:
    Tag st;                             //!< current state
    UCHAR _backpType : BITS4P_NODES;    //!< last state's parse type
    UCHAR backPhase : BITS4P_NODEPHASE; //!< last state's phase
    UCHAR backKth;                      //!< last state's k-best entry
    UCHAR wordLength;                   //!< length of word in #nodes.
                                        //!< from the last sync state
    double _score;                      //!< score up to this entry
  
    LatticeEdge(Tag &st, 
                int wordLength, 
                TParseNode backpType, 
                int backPhase, 
                int backKth, 
                double score);
    
    LatticeEdge(const LatticeEdge &ledge) {
      (*this) = ledge;
    }

    LatticeEdge() {}
    
    inline TParseNode backpType() {
      return (TParseNode)_backpType;
    }

    inline double score() {
      return _score;
    }
    
    LatticeEdge *backtrack(LatticeVertex ***vptr);
    
    inline void update(TParseNode backpType, int backPhase, int backKth) {
      this->_backpType = backpType;
      this->backPhase = backPhase;
      this->backKth = backKth;
    }
    
    /**
     * This function is needed to implement the k-best algorithm using a
     * a priority queue based on a min-heap
     */
    inline bool operator<(LatticeEdge &ledge) {
      return this->_score > ledge.score(); 
    }
    
    ~LatticeEdge() {};
  }; 
  
} 

#endif
