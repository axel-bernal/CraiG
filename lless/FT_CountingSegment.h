/****************************************************************************
 * FT_CountingSegment.h - part of the lless namespace, a general purpose
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

#ifndef FT_COUNTING_SEGMENT_H
#define FT_COUNTING_SEGMENT_H
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Sequence.h"
#include "NGram.h"
#include "Feature.h"

namespace lless {

    /**
     * The FT_BaseCountingSegment class computes its feature value as sets of filter
     * values over segments. It is the base class implementation for more
     * sophisticated ways to count the filter values, such as
     * FT_PeriodicBaseCountingSegment, which works with filters whose values are
     * periodic.
     * This feature works very efficiently during decoding when the provided
     * filter has previously accumulated its values along the sequence
     ***************************************************************************/

    template <class TClass>
        class FT_BaseCountingSegment : public TypedFeature<TClass> {
    protected:

        TypedFilter<TClass> *filter;
        bool accumulated;   //!< true if the values can be computed by accumulation

    public:
    FT_BaseCountingSegment(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass> *filter)
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               parsingFrames,
                               name,
                               fe,
                               filter->maxNumFilterValues()) {

            this->filter = filter;
            this->accumulated = false;

        }

    FT_BaseCountingSegment(
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
                               fe->getFilter(fargs[offset + 3])->maxNumFilterValues()) {

            this->filter = (TypedFilter<TClass> *)fe->getFilter(fargs[offset++]);
            this->accumulated = false;

            int numDisabledVals;
            if(!sscanf(fargs[offset++].c_str(), "%d", &numDisabledVals))
                assert(0);

            for(int i = 0; i < numDisabledVals; i++) {
                int fVal;
                if(!sscanf(fargs[offset++].c_str(), "%d", &fVal))
                    assert(0);
            }

        }

        virtual int period() = 0;
        virtual int phase() = 0;
        virtual int begTag(Tag *ge, int frame) = 0;
        virtual int endTag(Tag *ge, int frame) = 0;
        virtual int lenTag(Tag *ge, int frame) = 0;
        virtual double scaleDotParamV(Tag *ge, double result) = 0;

        virtual double posDotParamV(TClass & fValue, int pos, TStrand strand,
                                    const FeatureVector **params, int fConjInd = 0) = 0;

        virtual void posUpdParamV(double updVal, TClass & fValue, int pos, TStrand strand,
                                  FeatureVector **params, int fConjInd = 0) = 0;

        virtual inline TClass & posfValue(int pos, TStrand strand) {
            return this->filter->value(pos, strand);
        }

        inline TypedFilter<TClass> *getFilter() {
            return filter;
        }

        inline double featValue(Tag *ge, int frame) {
            assert(0);
            throw EXCEPTION(BAD_USAGE, this->getName() + " is a multiple-valued feature.\n");
            return 0;
        }

        virtual inline int begPreComp(int beg, int frame) {
            return beg + frame;
        }

        inline int endPreComp(int end, int frame) {
            return end + this->parsingFrames() - frame;
        }

        inline int stepPreComp() {
            return this->parsingFrames();
        }

        inline bool isFeatureSet() {
            return true;
        }

        /**
         * A member function that performs dot product for feature precomputation.
         * The precomputation occurs at position pos and frame, the  word found
         * at such position is used to choose the right entry in array *params.
         */
        inline double dotParamVPreComp(int pos, int frame,
                                       TStrand strand,
                                       const FeatureVector **params,
                                       int fConjInd = 0,
                                       TypedFilter<UCHAR> *filter = NULL) {

            pos = pos + phase();
            int offset = fConjInd;

            if(filter)
                offset += filter->value(pos, strand);

            TClass & fValue = posfValue(pos, strand);
            return posDotParamV(fValue, pos, strand, params, offset);

        }


        inline void accumulatePreCompEntries(double **array,
                                             int preCompArrayIndex,
                                             int arrayEntryPeriod,
                                             int beg, int end) {

            assert(arrayEntryPeriod == period());

            for(int phase = 0; phase < period(); phase++) {
                int arrayInd = preCompArrayIndex + this->_preCompArrayIndex[phase];
                for(int i = beg; i < end + period() + 3; i += period())
                    array[arrayInd][i + phase] += array[arrayInd][i + phase - period()];
            }

            /*      int k;
                    cerr << getName() << " " << preCompArrayIndex + _preCompArrayIndex[0] << endl;
                    for(k = 1; k < 10; k++)
                    cerr << array[preCompArrayIndex + _preCompArrayIndex[0]][k] << "\t";
                    cerr << endl;
                    for(k = fe->seqLength() - 10; k <= fe->seqLength(); k++)
                    cerr << array[preCompArrayIndex + _preCompArrayIndex[0]][k] << "\t";
                    cerr << endl;
            */
            this->accumulated = true;
        }

        inline double preComputedDotParamV(double **array,
                                           int preCompArrayIndex,
                                           Tag *ge, int frame) {

            int arrayInd = preCompArrayIndex + this->_preCompArrayIndex[frame];
            double result;

            if(!this->accumulated)
                result =  array[arrayInd][ge->getPos()];
            else {
                int beg = begTag(ge, frame) - period();
                int end = endTag(ge, frame);
                //        int end = endTag(ge, frame);

                //        if(end < beg)
                //          return 0;

                //        result = array[arrayInd][end]
                //	cerr << "\n" << this->getName() << " " << end << " " << beg << " " <<  array[arrayInd][end]  << " " << array[arrayInd][beg] << "\n";
                result = array[arrayInd][end] - array[arrayInd][beg];

            }

            return scaleDotParamV(ge, result);

        }

        virtual inline double dotParamV(Tag *ge, int frame,
                                        const FeatureVector **params,
                                        int fConjInd = 0,
                                        TypedFilter<UCHAR> *filter = NULL) {

            int beg = begTag(ge, frame) + phase();
            int end = endTag(ge, frame) + phase();

            double result = 0;
            //      cerr << this->getName() << " " << beg << " " << end << endl;
            for(int i = beg; i <= end; i += period()) {
                int offset = fConjInd;

                if(filter)
                    offset += filter->value(i, ge->getStrand());

                TClass &fValue = posfValue(i, ge->getStrand());
                result += posDotParamV(fValue, i, ge->getStrand(), params, offset);

                //	cerr << result << " ";
            }

            return scaleDotParamV(ge, result);
        }

        virtual inline void updParamV(double updVal,
                                      Tag *ge, int frame,
                                      FeatureVector **params, int fConjInd = 0,
                                      TypedFilter<UCHAR> *filter = NULL) {

            int beg = begTag(ge, frame) + phase();
            int end = endTag(ge, frame) + phase();

            updVal = scaleDotParamV(ge, updVal);

            for(int i = beg; i <= end; i += period()) {
                int offset = fConjInd;

                if(filter)
                    offset += filter->value(i, ge->getStrand());

                TClass & fValue = posfValue(i, ge->getStrand());
                posUpdParamV(updVal, fValue, i, ge->getStrand(), params, offset);

            }
        }

        virtual ~FT_BaseCountingSegment() {

        }

    };


    /**
     * The FT_BaseCountingSegment class computes its feature value as sets of filter
     * values over segments. It is the base class implementation for more
     * sophisticated ways to count the filter values, such as
     * FT_PeriodicBaseCountingSegment, which works with filters whose values are
     * periodic.
     * This feature works very efficiently during decoding when the provided
     * filter has previously accumulated its values along the sequence
     ***************************************************************************/

    template <class TClass>
        class FT_CountingSegment : public FT_BaseCountingSegment<TClass> {
    protected:

    public:
    FT_CountingSegment(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass> *filter)
        : FT_BaseCountingSegment<TClass>(fInd,
                                         paramInd,
                                         parsingFrames,
                                         name,
                                         fe,
                                         filter) {

        }

    FT_CountingSegment(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe)
        : FT_BaseCountingSegment<TClass>(fInd,
                                         paramInd,
                                         fargs,
                                         offset,
                                         fe) {

        }

        virtual inline double posDotParamV(TClass & fValue, int pos,
                                           TStrand strand,
                                           const FeatureVector **params,
                                           int fConjInd = 0) {

            assert(fValue >= 0 && fValue < this->maxNumFeatValues());
            return (*params[fConjInd])[fValue];

        }

        virtual inline void posUpdParamV(double updVal,
                                         TClass & fValue,
                                         int pos, TStrand strand,
                                         FeatureVector **params,
                                         int fConjInd = 0) {

            assert(fValue >= 0 && fValue < this->maxNumFeatValues());
            (*params[fConjInd])[fValue] += updVal;
        }

        virtual ~FT_CountingSegment() {

        }

    };

    /**
     * The FT_CountingAtPos is a subclass of FT_CountingSegment which
     * computes its feature value as the filter value occurring at the Tag's
     * position given as parameter
     ***************************************************************************/

    template <class TClass>
        class FT_CountingAtPos : FT_CountingSegment<TClass> {
    protected:
        bool normalize;
        double mean;

    public:
    FT_CountingAtPos(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass> *filter
        )
        : FT_CountingSegment<TClass>(fInd,
                                     paramInd,
                                     parsingFrames,
                                     name,
                                     fe,
                                     filter) {

            this->normalize = normalize;
            this->mean = mean;

        }

    FT_CountingAtPos(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_CountingSegment<TClass>(fInd,
                                     paramInd,
                                     fargs,
                                     offset,
                                     fe) {

            normalize = Utils::stringToBoolean(fargs[offset++]);
            if(!normalize) {
                if(!sscanf(fargs[offset++].c_str(), "%lf", &mean))
                    assert(0);
            }

        }

        inline int period() {
            return 1;
        }

        inline int phase() {
            return 0;
        }

        inline int begTag(Tag *ge, int frame) {
            return ge->getPos();
        }

        inline int lenTag(Tag *ge, int frame) {
            return 1;
        }

        inline int endTag(Tag *ge, int frame) {
            return ge->getPos();
        }

        inline double scaleDotParamV(Tag *ge, double result) {
            return (normalize ? result/ge->getLen() :  result/mean);
        }

        ~FT_CountingAtPos() {

        }

    };

    /**
     * The FT_TrimmedCountingSegment is a subclass of FT_CountingSegment which
     * computes its feature value as sets of filter values over segments.
     * The main difference with FT_CountingSegment is that the segment in this
     * case has been trimmed from both sides to avoid counting spurious
     * words that appear on the segment boundaries.
     * This feature works very efficiently during decoding when the provided
     * filter has previously accumulated its values along the sequence
     ***************************************************************************/

    template <class TClass>
        class FT_TrimmedCountingSegment : public FT_CountingSegment<TClass> {

    protected:
        int fromLeft, toRight;

    public:
    FT_TrimmedCountingSegment(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass> *filter,
        int fromLeft,
        int toRight
        )
        : FT_CountingSegment<TClass>(fInd,
                                     paramInd,
                                     parsingFrames,
                                     name,
                                     fe,
                                     filter) {

            this->fromLeft = fromLeft;
            this->toRight = toRight;
        }

    FT_TrimmedCountingSegment(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_CountingSegment<TClass>(fInd,
                                     paramInd,
                                     fargs,
                                     offset,
                                     fe) {

            if(!sscanf(fargs[offset++].c_str(), "%d", &fromLeft))
                assert(0);
            if(!sscanf(fargs[offset++].c_str(), "%d", &toRight))
                assert(0);

        }

        inline int period() {
            return 1;
        }

        inline int phase() {
            return 0;
        }

        inline int begTag(Tag *ge, int frame) {
            return ge->getPos() + fromLeft;
        }

        inline int lenTag(Tag *ge, int frame) {
            return ge->getLen() - (fromLeft + toRight);
        }

        inline int endTag(Tag *ge, int frame) {
            int end = ge->getPos() + ge->getLen() - toRight;

            if(end > this->fe->seqLength())
                end = this->fe->seqLength();

            return end;
        }

        virtual inline double scaleDotParamV(Tag *ge, double result) {
            return result;
        }

        ~FT_TrimmedCountingSegment() {

        }

    };

    /**
     * The FT_AggregatedTrimmedCountingSegment is a subclass of
     * FT_TrimmedCountingSegment which
     * computes its feature value as sets of filter values over segments.
     * The main difference with FT_TrimmedCountingSegment is that in this class
     * there are only two filter value types: nonAggrFVal which is any
     * enabled  value from the original filter and aggrFVal, to which
     * any other value of the original filter is assigned to.
     * This feature works very efficiently during decoding when the provided
     * filter has previously accumulated its values along the sequence
     ***************************************************************************/

    template <class TClass>
        class FT_AggregatedTrimmedCountingSegment : public FT_TrimmedCountingSegment<TClass> {
    protected:
        TClass aggrFVal, nonAggrFVal;

    public:
    FT_AggregatedTrimmedCountingSegment(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass> *filter,
        int fromLeft,
        int toRight,
        TClass nonAggrFVal
        )
        : FT_TrimmedCountingSegment<TClass>(fInd,
                                            paramInd,
                                            parsingFrames,
                                            name,
                                            fe,
                                            filter,
                                            fromLeft,
                                            toRight) {

            this->nonAggrFVal = nonAggrFVal;
            findAggrIndex();

        }

    FT_AggregatedTrimmedCountingSegment(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_TrimmedCountingSegment<TClass>(fInd,
                                            paramInd,
                                            fargs,
                                            offset,
                                            fe,
                                            fte) {

            if(!sscanf(fargs[offset++].c_str(), "%d", &nonAggrFVal))
                assert(0);
            findAggrIndex();
        }

        /**
         * Finds a valid index which is an enabled filter value
         * to which all other filter values except the one given as nonAggregated,
         * will aggregate, i.e. will be assigned to.
         */
        inline void findAggrIndex() {
            aggrFVal = 0;
            for(; aggrFVal < this->maxNumFeatValues(); aggrFVal++) {
                if(aggrFVal == nonAggrFVal)
                    continue;
                break;
            }
        }

        /**
         * @return either aggrFVal or nonAggrFVal or the original filter value
         * if it has not been enabled.
         */
        inline TClass & posfValue(int pos, TStrand strand) {
            TClass fVal = this->filter->value(pos, strand);
            assert(fVal < this->maxNumFeatValues());
            if(fVal == nonAggrFVal)
                return fVal;

            return aggrFVal;
        }

        ~FT_AggregatedTrimmedCountingSegment() {

        }

    };

}

#endif
