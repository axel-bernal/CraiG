/****************************************************************************
* LatticeVertex.h - part of the lless namespace, a general purpose
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

#ifndef _LATTICEVERTEX_H_
#define _LATTICEVERTEX_H_
#include "EdgeInst.h"
#include "NodeInst.h"
#include "FeatureEngine.h"
#include <algorithm>
#include "FSMQueue.h"
#include "FT_HiddenSeq.h"
#include "LatticeEdge.h"

namespace lless { 
 
  /**
   * LatticeNode implements a node in a Lattice object. 
   */
  class LatticeVertex { 
  public:
    vector<LatticeEdge> *_topKEdges;  //!< stores the k-best entries
    UCHAR capacity;                   //!< total capacity
    
    LatticeVertex() {
      capacity = 0;
      _topKEdges = NULL;
    }
    
    LatticeVertex(int capacity) {
      this->capacity = capacity;
      this->_topKEdges = new vector<LatticeEdge>;
      assert(this->_topKEdges);
    }
  
    inline int size() {
      return _topKEdges->size();
    }
    
    inline LatticeEdge * operator[](int kth) {
      return &(*_topKEdges)[kth];
    }

    /**
     * In adding a back edge, we are cautious on maintaining a priority
     * queue, so that we can retrieve the back edge with the best score
     * in constant time while we are decoding.
     */
    inline bool addLatticeEdge(LatticeEdge &ledge) {
      assert(_topKEdges);
      if(size() < capacity) {
        _topKEdges->push_back(ledge);
        if(size() == capacity)
          make_heap(_topKEdges->begin(), _topKEdges->end());
        return true;
      }
      if((*_topKEdges)[0].score() < ledge.score()) {
        pop_heap(_topKEdges->begin(), _topKEdges->end());
        _topKEdges->pop_back();
        _topKEdges->push_back(ledge);
        push_heap(_topKEdges->begin(), _topKEdges->end());
        return true;
      }

      /*Tag *b = &sbt[indexMinVal]->st;
        cerr << "    rejected " << " " << st.getPos() << "-" << st.getPos() + st.getLen() - 1 << " " << (backPhase + st.getLen()) % 3 << " by " << " " << b->getPos() << "-" << b->getPos() + b->getLen() - 1 << "  " << (sbt[indexMinVal]->backPhase + b->getLen()) % 3 << " " << sbt[indexMinVal]->score << endl;
        } */
      return false;

    }

    ~LatticeVertex() {
      if(_topKEdges) {
        delete _topKEdges;
        _topKEdges = NULL;
      }
    }
  };

}

#endif
