/****************************************************************************
 * IntervalSet.h - part of the lless namespace, a general purpose
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
 *
 ****************************************************************************/

#ifndef _INTERVAL_SET_H_
#define _INTERVAL_SET_H_
#include <list>
#include <vector>
#include "VectorUtils.h"

namespace lless {

    /***************************************************************************
     * The Vector abstract base class defines the interface of a vector, which
     * can be implemented as a sparse vector or dense vector.
     ***************************************************************************/

    struct IntPairHash {
    private:
        std::hash<int> ah;
        std::hash<int> bh;
    public:
    IntPairHash() : ah(), bh() {}
        size_t operator()(const std::pair<int, int> &p) const {
            return ah(p.first) ^ bh(p.second);
        }
    };

    struct IntPairEqual {
    public:
        bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const {
            return a.first == b.first && a.second == b.second;
        }
    };

    template <class TClass> class IntervalSet : public DENSE_HASH< pair<int, int>,TClass, IntPairHash, IntPairEqual >  {
    protected:

    public:
    };

}

#endif
