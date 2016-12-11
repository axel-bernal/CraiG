/****************************************************************************
* PrefixTree.h - part of the lless namespace, a general purpose
*                linear semi-markov structure prediction library
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

#ifndef _PREFIX_TREE_H_
#define _PREFIX_TREE_H_

namespace lless {

  /**
   * The PrefixTree class implements a prefix tree. It is used for storing
   * strings and efficiently do dictionary queries to it.
   ***************************************************************************/


  class PrefixTree {
   private:
    int _order;
    PrefixTree **_children;
   public:
    PrefixTree(int order) {
      _order = order;
      _children = new PrefixTree *[_order]; 
      for(int i = 0; i < _order; i++)
        _children[i] = NULL;
    }

    inline PrefixTree *operator[](int index) {  return _children[index];}
    inline PrefixTree **children() {  return _children; }  
    inline PrefixTree *addChild(int index) {
      assert(index < _order);
      if(!_children[index])
        _children[index] = new PrefixTree(_order);

      return _children[index];
    }

    ~PrefixTree() {
      for(int i = 0; i < _order; i++) {
        if(_children[i] != NULL) 
          delete _children[i];
      }

      delete [] _children;

    }
  };

}

#endif
