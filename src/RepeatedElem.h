/****************************************************************************
 * RepeatedElem.h - part of the craig namespace, a genomics library
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

#ifndef _REPEATED_ELEM_H_
#define _REPEATED_ELEM_H_

#include "BioFeature.h"

namespace craig {

    /**
     * The RepeatedElem class represents a previously-masked low density region
     * in the DNA
     ***************************************************************************/

    class RepeatedElem : public BioFeature {
    private:
        int _begin, _end;
    public:
    RepeatedElem(std::string & iden,
                 int begin, int end,
                 Sequence *c,
                 bool inOneStrand,
                 TStrand strand
        ) : BioFeature(iden, c, inOneStrand, strand) {
            _begin = begin;
            _end = end;
        }

        RepeatedElem(const RepeatedElem &re) {
            *this = (RepeatedElem &)re;
        }

        RepeatedElem & operator=(const RepeatedElem & repeatedElem) {
            RepeatedElem &re = (RepeatedElem &)repeatedElem;
            (BioFeature &)*this = (BioFeature &)re;
            _begin = re.begin();
            _end = re.end();
            return *this;
        }

        bool operator==(const RepeatedElem & repeatedElem) {
            RepeatedElem &re = (RepeatedElem &)repeatedElem;
            bool equal = (_begin == re.begin() && _end == re.end());
            return equal &&  (BioFeature &)*this == (BioFeature &)re;
        }

        bool operator<(const RepeatedElem & repeatedElem) {
            RepeatedElem &re = (RepeatedElem &)repeatedElem;
            int b = getSequence()->id().compare(re.getSequence()->id());
            if(b)
                return b < 0;

            assert(inOneStrand == re.isInOneStrand());
            return _begin > re.begin();
        }

        inline void toOneStrand() {
            _end = annotSeq->length() - _end + 1;
            _begin = annotSeq->length() - _begin + 1;
        }

        inline void toTwoStrand() {
            _end = annotSeq->length() - _end + 1;
            _begin = annotSeq->length() - _begin + 1;
        }

        inline int begin() {  return _begin; }
        inline int end() {  return _end; }

        ~RepeatedElem() {  }
    };
}

#endif
