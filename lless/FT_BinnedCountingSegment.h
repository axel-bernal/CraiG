/****************************************************************************
 * FT_BinnedCountingSegment.h - part of the lless namespace, a general purpose
 *                       linear semi-markov structure prediction library
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

#ifndef FT_BINNEDCOUNT_SEGMENT_FEAT_H
#define FT_BINNEDCOUNT_SEGMENT_FEAT_H

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Sequence.h"
#include "Feature.h"
#include "Filter.h"
#include "FilterEngine.h"
#include "FeatureEngine.h"
#include "FT_CountingSegment.h"
#include "FT_Bin.h"

namespace lless {

    /**
     * FT_BinnedCountingSegment is a subclass of TypedFeature whose value of this
     * feature (and its derivate classes) is the score obtained by accumulating
     * values of the contained filter over a tag which must be a segment
     * (node with variable length).
     * This feature works well with bimodal filter values, in which scores above
     * a certain threshold (the mean) should be classified differently.
     * Invariants: the contained filter must accumulate.
     *
     ***************************************************************************/
    template <class TClass1, class TClass2>
        class FT_BinnedCountingSegment : public FT_BaseCountingSegment<TClass2> {

    protected:
        int numBinObjects;
        ULONG domain, range;
        vector<FT_BaseBin<TClass1> *> binObjects;
        vector<TypedFilter<TClass1> *> binValues;
        vector<int> binObj4fvals;  //!< which bin object to use for each filter val

    public:
    FT_BinnedCountingSegment(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass2> *filter,
        int numBinObjects
        )
        : FT_BaseCountingSegment<TClass2>(fInd,
                                          paramInd,
                                          parsingFrames,
                                          name,
                                          fe,
                                          filter) {

            assert(this->filter->period());
            initialize(numBinObjects);

        }

    FT_BinnedCountingSegment(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_BaseCountingSegment<TClass2>(fInd, paramInd,
                                          fargs, offset, fe) {

            int numBinObjects;

            if(!sscanf(fargs[offset++].c_str(), "%d", &numBinObjects))
                assert(0);

            initialize(numBinObjects);

            FT_BaseBin<TClass1> *binObj;
            TypedFilter<TClass1> *binVal;
            vector<bool> used_fv(this->filter->maxNumFilterValues(), false);

            for(int i = 0; i < numBinObjects; i++) {
                int numfvals, fval;

                binObj = (FT_BaseBin<TClass1> *)fte->getFeature(fargs[offset++]);
                binVal = (TypedFilter<TClass1> *)fe->findFilter(fargs[offset++]);
                addBinObj(i, binObj, binVal);

                if(!sscanf(fargs[offset++].c_str(), "%d", &numfvals))
                    assert(0);

                if(numfvals == this->filter->maxNumFilterValues()) {
                    for(int j = 0; j < numfvals; j++) {
                        binObj4fvals[j] = i;
                        used_fv[j] = true;
                    }
                    continue;
                }
                // assign all non used filter values to this bin.
                if(!numfvals) {
                    for(int j = 0; j < this->filter->maxNumFilterValues(); j++) {
                        if(used_fv[j])
                            continue;
                        binObj4fvals[j] = i;
                        used_fv[j] = true;
                    }
                    continue;
                }

                for(int j = 0; j < numfvals; j++) {
                    if(!sscanf(fargs[offset++].c_str(), "%d", &fval))
                        assert(0);
                    binObj4fvals[fval] = i;
                    used_fv[fval] = true;
                }
            }
        }

        void initialize(int numBinObjects) {
            range = 0;
            domain = this->filter->maxNumFilterValues();
            this->setMaxNumParams(domain, range);
            this->numBinObjects = numBinObjects;
            this->binObjects = vector<FT_BaseBin<TClass1> *>(numBinObjects, NULL);
            this->binValues = vector<TypedFilter<TClass1> *>(numBinObjects, NULL);
            this->binObj4fvals = vector<int>(domain, -1);
        }

        void addBinObj(int bin, FT_BaseBin<TClass1> *binObj,
                       TypedFilter<TClass1> *binVal) {

            assert(bin <= numBinObjects);
            binObjects[bin] = binObj;
            binValues[bin] = binVal;

            if(binObj->maxNumFeatValues() > range) {
                range = binObj->maxNumFeatValues();
                this->_numParams.s = range;
            }
        }


        inline double posDotParamV(UCHAR &fValue,
                                   int pos, TStrand strand,
                                   const FeatureVector **params,
                                   int fConjInd = 0) {
            int nv = (int)fValue;
            return posDotParamV(nv, pos, strand, params, fConjInd);

        }


        inline double posDotParamV(int &fValue,
                                   int pos, TStrand strand,
                                   const FeatureVector **params,
                                   int fConjInd = 0) {

            int offset = fConjInd*domain;

            if(binObj4fvals[fValue] < 0 || !this->subFeatIsOn(fValue))
                return 0;

            int objNo = binObj4fvals[fValue];
            TClass1 binVal = binValues[objNo]->value(pos, strand);

            return binObjects[objNo]->dotParamV((double)binVal,
                                                params, offset + fValue);


        }


        inline double posDotParamV(SPARSE_HASH<UCHAR, int> & fValue,
                                   int pos, TStrand strand,
                                   const FeatureVector **params,
                                   int fConjInd = 0) {

            int offset = fConjInd*domain;
            typename SPARSE_HASH<UCHAR, int>::iterator it = fValue.begin();

            double pdotp = 0;
            for( ; it != fValue.end(); it++) {

                UCHAR fv = it->first;
                int binVal = it->second;

                if(binObj4fvals[fv] < 0 || !this->subFeatIsOn(fv))
                    continue;

                int objNo = binObj4fvals[fv];
                pdotp += binObjects[objNo]->dotParamV((double)binVal,
                                                      params, offset + fv);
            }

            return pdotp;
        }

        inline void posUpdParamV(double updVal,
                                 UCHAR & fValue,
                                 int pos, TStrand strand,
                                 FeatureVector **params,
                                 int fConjInd = 0) {
            int nv = (int)fValue;
            posUpdParamV(updVal, nv, pos, strand, params, fConjInd);
        }

        inline void posUpdParamV(double updVal,
                                 int & fValue,
                                 int pos, TStrand strand,
                                 FeatureVector **params,
                                 int fConjInd = 0) {

            int offset = fConjInd*domain;

            if(binObj4fvals[fValue] < 0 || !this->subFeatIsOn(fValue))
                return;

            int objNo = binObj4fvals[fValue];
            TClass1 binVal = binValues[objNo]->value(pos, strand);

            //      cerr << this->getName() << " " << strand << " " << pos << " " << fValue << " " << (double)binVal << endl;;

            binObjects[objNo]->updParamV(updVal, (double)binVal,
                                         params, offset + fValue);

        }


        inline void posUpdParamV(double updVal,
                                 SPARSE_HASH<UCHAR, int> &fValue,
                                 int pos, TStrand strand,
                                 FeatureVector **params,
                                 int fConjInd = 0) {

            int offset = fConjInd*domain;
            typename SPARSE_HASH<UCHAR, int>::iterator it = fValue.begin();

            for( ; it != fValue.end(); it++) {
                UCHAR fv = it->first;
                int binVal = it->second;

                if(binObj4fvals[fv] < 0 || !this->subFeatIsOn(fv))
                    continue;

                int objNo = binObj4fvals[fv];
                binObjects[objNo]->updParamV(updVal, (double)binVal,
                                             params, offset + fv);
            }
        }

        ~FT_BinnedCountingSegment() { }

    };


    /**
     * FT_BinnedPeriodicCountingSegment is a subclass of FT_BinnedCountingSegment and
     * computes its feature value as sets of filter values over segments.
     * The filter values are periodic.
     * This feature works very efficiently during decoding when the provided
     * filter has previously accumulated its values along the sequence
     ***************************************************************************/
    template <class TClass1, class TClass2>
        class FT_BinnedPeriodicCountingSegment : public FT_BinnedCountingSegment<TClass1, TClass2> {

    protected:
        int _phase, _period;
        bool normalize;
        double mean;

    public:
    FT_BinnedPeriodicCountingSegment(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass2> *filter,
        int numBinObjects,
        int period,
        int phase = 0,
        int normalize = false,
        double mean = 1.0
        )
        : FT_BinnedCountingSegment<TClass1, TClass2>(fInd,
                                                     paramInd,
                                                     parsingFrames,
                                                     name,
                                                     fe,
                                                     filter,
                                                     numBinObjects) {

            this->_period = period;
            this->_phase = phase;
            this->normalize = normalize;
            this->mean = mean;

        }

        /**
         * Constructor of a BinnedPeriodicCountingSegment in the feature file
         *
         * Feature  id  BinnedPeriodicCountingSegment 3 3-gram 0 3 1 false 1
         *
         * The definition above creates a BinnedPeriodicCountingSegment feature,
         * which extends BinnedCountingSegment to deal with filters whose values
         * are periodic. This feature has 3 phases, the filter 3-gram has values
         * which are 3-periodic and the starting phase of the feature is 1. The
         * last two parameters indicate that the feature should not normalize
         * the countings (false) and instead should use 1 as default value for
         * the standard deviation of the counts distribution.
         *
         */

    FT_BinnedPeriodicCountingSegment(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_BinnedCountingSegment<TClass1, TClass2>(fInd,
                                                     paramInd,
                                                     fargs,
                                                     offset,
                                                     fe,
                                                     fte) {

            if(!sscanf(fargs[offset++].c_str(), "%d", &_period))
                assert(0);
            if(!sscanf(fargs[offset++].c_str(), "%d", &_phase))
                assert(0);

            normalize = Utils::stringToBoolean(fargs[offset++]);
            if(!normalize) {
                if(!sscanf(fargs[offset++].c_str(), "%lf", &mean))
                    assert(0);
            }

        }

        inline int phase() {
            return _phase;
        }

        inline int period() {
            return _period;
        }

        virtual inline int begTag(Tag *ge, int frame) {
            int beg = ge->getPos();
            return beg + period()*(frame != 0) - frame;
        }

        virtual inline int lenTag(Tag *ge, int frame) {
            return period()*(ge->getLen()/period());
        }

        virtual inline int endTag(Tag *ge, int frame) {
            int end = begTag(ge, frame) + lenTag(ge, frame);
            end -= period()*(ge->getLen() % period() == 0);

            if(end > this->fe->seqLength())
                end -= period();

            return end;
        }

        virtual inline double scaleDotParamV(Tag *ge, double result) {
            return (normalize ? result/ge->getLen() :  result/mean);
        }

        ~FT_BinnedPeriodicCountingSegment() {

        }
    };


    /**
     * FT_BinnedPeriodicCountingSegmentWUpperLimit is a subclass of
     * FT_BinnedCountingSegment and computes its feature value as sets of
     * filter values over segments.
     * The filter values are periodic.
     * This feature works very efficiently during decoding when the provided
     * filter has previously accumulated its values along the sequence
     ***************************************************************************/
    template <class TClass1, class TClass2>
        class FT_BinnedPeriodicCountingSegmentWUpperLimit : public FT_BinnedPeriodicCountingSegment<TClass1, TClass2> {

    protected:
        int limit;

    public:
    FT_BinnedPeriodicCountingSegmentWUpperLimit(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass2> *filter,
        int numBinObjects,
        int period,
        int phase = 0,
        int normalize = false,
        double mean = 1.0,
        int limit = 15
        )
        : FT_BinnedPeriodicCountingSegment<TClass1, TClass2>(fInd,
                                                             paramInd,
                                                             parsingFrames,
                                                             name,
                                                             fe,
                                                             filter,
                                                             numBinObjects,
                                                             period, phase,
                                                             normalize, mean) {

            this->limit = limit;

        }

        /**
         * Constructor of a BinnedPeriodicCountingSegmentWUpperLimit in the feature file
         *
         * Feature  id  BinnedPeriodicCountingSegmentWUpperLimit 3 3-gram 0 3 1 false 1
         *
         * The definition above creates a BinnedPeriodicCountingSegmentWUpperLimit feature,
         * which extends BinnedCountingSegment to deal with filters whose values
         * are periodic. This feature has 3 phases, the filter 3-gram has values
         * which are 3-periodic and the starting phase of the feature is 1. The
         * last two parameters indicate that the feature should not normalize
         * the countings (false) and instead should use 1 as default value for
         * the standard deviation of the counts distribution.
         *
         */

    FT_BinnedPeriodicCountingSegmentWUpperLimit(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_BinnedPeriodicCountingSegment<TClass1, TClass2>(fInd,
                                                             paramInd,
                                                             fargs,
                                                             offset,
                                                             fe,
                                                             fte) {

            if(!sscanf(fargs[offset++].c_str(), "%d", &limit))
                assert(0);

        }

        inline int lenTag(Tag *ge, int frame) {
            int len = (ge->getLen() > limit)
                ? limit
                : ge->getLen();

            return this->period()*(len/this->period());
        }

        inline int endTag(Tag *ge, int frame) {
            int len = (ge->getLen() > limit)
                ? limit
                : ge->getLen();

            int end = this->begTag(ge, frame) + this->period()*len/this->period();
            end -= this->period()*(len % this->period() == 0);

            if(end > this->fe->seqLength())
                end -= this->period();

            return end;

        }

        ~FT_BinnedPeriodicCountingSegmentWUpperLimit() {

        }
    };

}

#endif
