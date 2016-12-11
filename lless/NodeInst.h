/****************************************************************************
* NodeInst.h - part of the lless namespace, a general purpose
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

#ifndef _NODE_INST_H_
#define _NODE_INST_H_

#include "Utils.h"
#include "Tag.h"


namespace lless {

  /**
   * The NodeInst class represents instances of lattice nodes in the sequence.
   * it derives from the Tag class the parseType and position among other
   * things. It adds up to this, a length and an internal type.
   ***************************************************************************/
  class NodeInst : public Tag {
   protected:
    UCHAR type : BITS4P_ID3NODES;

   public:
    /**
     * Main constructor
     * @param type the id3 of the node, one of TNodeId3 enumerate.
     * @param pType the parsing type of the node, one of TParseNode enumerate.
     * @param pos the starting position of the node in the sequence.
     * @param strand the strand of the node.
     * @param len the length of the node.
     */
    NodeInst(TNodeId3 type, 
             TParseNode pType, 
             int pos, 
             TStrand strand, 
             int len) : Tag(pType, pos, strand) {

      this->type = type;
      this->geClass = NODE_INST;
      setLen(len);

    }
    
    /**
     * Copy constructor
     */
    NodeInst(const NodeInst & st) {
      *this = st;
    }

    /**
     * Constructor from a Tag object
     * @param ge the Tag object.
     * @param type the id3 of the node, one of TNodeId3 enumerate.
     * @param len the length of the node.
     */
    NodeInst(Tag & ge, TNodeId3 type, int len) : Tag(ge) {
      this->type = type;
      this->len = len;
    }

    NodeInst() :  Tag(INVALID_NODE, 0, NO_STRAND) {
      this->type = NO_NODE_INST;
      this->len = 0;
      this->geClass = NODE_INST;
    }
    
    inline int getType() {
      return type;
    }
        
    ~NodeInst() { 
    }
  };

}

#endif

