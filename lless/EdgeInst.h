/****************************************************************************
* EdgeInst.h - part of the lless namespace, a general purpose
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

#ifndef _EDGE_INST_H_
#define _EDGE_INST_H_

#include "Utils.h"
#include "Tag.h"


namespace lless {

  /**
   * The EdgeInst class represents instances of lattice edges in the sequence.
   * it derives from the Tag class the parseType and position among other
   * things. Although the _eclipsed field has been implemented with DNA 
   * sequences
   * in mind, it only is used effectively when the FSM specification contains
   * disruptor signals whose TParseType is different from INVALIG_EDGE.
   ***************************************************************************/

  class EdgeInst : public Tag {
   protected:
    UCHAR type : BITS4P_ID2EDGES;
    bool _eclipsed[3];
    int len_prev, len_next;
   public:
    /**
     * Default constructor
     */
    EdgeInst() : Tag(INVALID_EDGE, 0, NO_STRAND) {
      type = NO_EDGE_INST;
      this->geClass = EDGE_INST;
      setLen(0);
      uneclipse();
      len_prev = len_next = 0;
    }

    /*
     * Constructor with no parsing Type
     * @param type the id2 of the edge, one of TEdgeId2 enumerate.
     * @param pos the starting position of the node in the sequence.
     * @param strand the strand of the node.
     */
     
    EdgeInst(int type, int pos, TStrand strand) :
      Tag(INVALID_EDGE, pos, strand) {

      this->type = type;
      this->geClass = EDGE_INST;
      setLen(0);
      uneclipse();
      len_prev = len_next = 0;
    }

    /**
     * Copy constructor
     */
    EdgeInst(const EdgeInst &bs) {
      *this = bs;
    }

    /**
     * Main constructor
     * @param type the id2 of the edge, one of TEdgeId2 enumerate.
     * @param pType the parsing type of the edge, one of TParseEdge enumerate.
     * @param pos the starting position of the node in the sequence.
     * @param strand the strand of the node.
     */
    EdgeInst(int type, 
             TParseEdge parseType, 
             int pos, 
             TStrand strand) : Tag(parseType, pos, strand) {

      this->type = type;
      this->geClass = EDGE_INST;
      setLen(0);
      uneclipse();
      len_prev = len_next = 0;
    }

    inline bool operator==(Tag &ge) {
      EdgeInst &bs = (EdgeInst &)ge;
      return pos == (unsigned int)bs.getPos() &&
        type == bs.getType() && 
        (TStrand)strand == bs.getStrand();
    }
    
    inline void setAssocNodeLens(int lprev, int lnext) {
      this->len_prev = lprev;
      this->len_next = lnext;
    }

    inline int lenPrevNode() {
      return this->len_prev;
    }
    
    inline int lenNextNode() {
      return this->len_next;
    }
    
    inline int getLen() {  return 0; }

    inline int getType() {  return type; }

    /**
     * Eclipse this edge's phase. 
     * \todo the maximum number of phases should come as parameter, we need
     * to change this. Right now is set to 3.
     */
    inline void eclipse(int phase) {
      _eclipsed[phase] = true;
    }

    /*
     * @return true is edge has been eclipse in phase, false otherwise
     */
    inline bool eclipsed(int phase) {
      return _eclipsed[phase];
    }

    inline void uneclipse() {
      for(int i = 0; i < 3; i++) 
        _eclipsed[i] = false;
    }

    inline void uneclipse(int phase) {
      _eclipsed[phase] = false;
    }

    inline bool operator==(int j) {
      return type == (UCHAR)j;
    }

    ~EdgeInst() { }
  };

}

#endif
