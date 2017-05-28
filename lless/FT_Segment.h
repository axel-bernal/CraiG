/****************************************************************************
 * FT_Segment.h - part of the lless namespace, a general purpose
 *                        linear semi-markov structure prediction library
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

#ifndef FT_SEGMENT_H
#define FT_SEGMENT_H
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Feature.h"
#include "Sequence.h"
#include "NGram.h"


namespace lless {

    /**
     * The FT_Segment class computes its feature value as sets of filter
     * values over segments. It is the base class implementation for more
     * sophisticated ways to count the filter values, such as
     * FT_PeriodicSegment, which works with filters whose values are
     * periodic.
     * This feature works very efficiently during decoding when the provided
     * filter has previously accumulated its values along the sequence
     ***************************************************************************/

    template <class TClass> class FT_Segment : public TypedFeature<TClass> {
    protected:
        TypedFilter<TClass> *filter;
        int _period;

    public:
    FT_Segment(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass> *filter,
        int period
        )
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               parsingFrames,
                               name,
                               fe,
                               filter->maxNumFilterValues(),
                               filter->valType()) {

            this->filter = filter;
            this->_period = period;
        }

    FT_Segment(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass> *filter,
        int period,
        int numFeatVals
        )
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               parsingFrames,
                               name,
                               fe,
                               numFeatVals,
                               filter->valType()) {

            this->filter = filter;
            this->_period = period;

        }


    FT_Segment(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe)
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               fargs,
                               offset,
                               fe,
                               fe->getFilter(fargs[offset + 3])
                               ->maxNumFilterValues(),
                               fe->getFilter(fargs[offset + 3])
                               ->valType()) {

            this->filter = (TypedFilter<TClass> *)fe->getFilter(fargs[offset++]);
            if(!sscanf(fargs[offset++].c_str(), "%d", &_period))
                assert(0);

        }

    FT_Segment(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        int numFeatVals)
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               fargs,
                               offset,
                               fe,
                               numFeatVals,
                               fe->getFilter(fargs[offset + 3])
                               ->valType()) {

            this->filter = (TypedFilter<TClass> *)fe->getFilter(fargs[offset++]);
            if(!sscanf(fargs[offset++].c_str(), "%d", &_period))
                assert(0);

        }

    FT_Segment(
        int fInd,
        int paramInd,
        int parsingFrames,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        int period,
        int numFeatVals)
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               parsingFrames,
                               fargs,
                               offset,
                               fe,
                               numFeatVals,
                               fe->getFilter(fargs[offset + 2])
                               ->valType()) {
            this->_period = period;
            this->filter = (TypedFilter<TClass> *)fe->getFilter(fargs[offset++]);

        }

        int period() {
            return _period;
        }

        virtual inline int begTag(Tag *ge, int frame) {
            int beg = ge->getPos() + period()*(frame != 0) - frame;
            int excess = -(beg - 1);
            if(excess > 0)
                beg += period()*(excess/period());

            return beg;
        }

        virtual inline int lenTag(Tag *ge, int frame) {
            return period()*(ge->getLen()/period());
        }

        virtual inline int endTag(Tag *ge, int frame) {
            int beg = ge->getPos() + period()*(frame != 0) - frame;
            int end = beg + lenTag(ge, frame) - period()*(ge->getLen() % period() == 0);
            int excess = end - this->fe->seqLength();
            if(excess > 0)
                end -= period()*(excess/period());
            if(end < beg) end = beg;

            return end;
        }

        inline TypedFilter<TClass> *getFilter() {
            return filter;
        }

        virtual inline double featValue(Tag *ge, int frame) {
            assert(0);
            throw EXCEPTION(BAD_USAGE, this->getName() + " is a virtual feature.\n");
            return 0;
        }

        virtual ~FT_Segment() {
        }

    };

}

#endif
