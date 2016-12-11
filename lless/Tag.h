/****************************************************************************
* Tag.h - part of the lless namespace, a general purpose
*               linear semi-markov structure prediction library
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

#ifndef __TAG_H__
#define __TAG_H__
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <TypeDefs.h>



namespace lless {
  /**
   * An enumerate to distinguish between Tag objects that are nodes and tag
   * objects that are edges.
   */
  typedef enum {NODE_INST, EDGE_INST, WORD_INST} TTagClass;
  
  /**
   * The Tag class is either an instance of a node or an edge as they appear in
   * the sequence.
   ***************************************************************************/
  class Tag {
  protected:
    UINT pos;
    UINT len;
    UCHAR pType : BITS4P_TAG;
    UCHAR geClass : BITS4P_TAG_CLASS;
    UCHAR strand : BITS4P_STRAND; 

  public:
    //! Default constructor
    Tag() {}

    /**
     * Main constructor
     * @param pType the parsing Type of the tag object.
     * @param pos the position in the sequence.
     * @param strand the strand of the Tag object.
     */
    Tag(int pType, int pos = 0, TStrand strand = NO_STRAND) {
      this->pType = pType;
      this->pos = pos;
      this->strand = strand;
    }

    Tag(const Tag &tag) {
      *this = tag;
    }

    inline int getLen() {
      return len;
    }

    inline void setLen(int len) {
      this->len = len;
    }
    
    inline void reInitialize(int pos, int len) {
      setPos(pos);
      setLen(len);
    }
    
    // getters
    inline int getPos() {
      return pos;
    }

    inline int getParseType() {  return pType; }

    inline bool operator<(const Tag &tag) {
      Tag &t = (Tag &)tag;

      if(this->strand != t.getStrand()) 
        return this->strand < t.getStrand();

      if(this->pos < t.getPos())
        return this->pos < t.getPos();
      
      if(this->geClass != t.getGEClass()) 
        return this->geClass < t.getGEClass();
      
      return this->pType < t.getParseType();

    }

    inline bool operator!=(const Tag &tag) {
      return !((*this) == tag);
    }

    inline bool operator==(const Tag &tag) {
      Tag &t = (Tag &)tag;
      
      if(this->strand != t.getStrand()) 
        return false;

      if(this->pos != t.getPos() || this->len != t.getLen())
        return false;
      
      if(this->geClass != t.getGEClass() || this->pType != t.getParseType())
        return false;
      
      return true;      
    }

    inline TStrand getStrand() {  return (TStrand)strand; }
    inline TTagClass getGEClass() {  return (TTagClass)geClass; }
    // setters
    inline void setPos(int pos) {  this->pos = pos; }
    inline void setStrand(TStrand strand) {  this->strand = strand; } 
    inline void setParseType(int pType) {  this->pType = pType; }

    inline void reInitialize(int pType, int pos, TStrand strand) {
      setParseType(pType);
      setPos(pos);
      setStrand(strand);
    }

    ~Tag() {
    }
    
  };

}

namespace std {
  /**                                                                            
   * Routine for sorting STL containers instantiated with Tag* types
   */
  template <>
    struct greater<lless::Tag *> {
      bool operator()(const lless::Tag* t1, const lless::Tag* t2)
      {
        //Defined for ascending sorting
        if(!t1)
          return true;
        if(!t2)
          return false;
        return ((lless::Tag &)(*t1)) < (*t2);
      }
  };
}

#endif

