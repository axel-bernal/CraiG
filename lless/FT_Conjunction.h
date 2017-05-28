/****************************************************************************
 * FT_Conjunction.h - part of the lless namespace, a general purpose
 *                    linear semi-markov structure prediction library
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

#ifndef _FT_CONJUNCTION_H_
#define _FT_CONJUNCTION_H_

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Feature.h"

using namespace lless;

namespace lless {

    /**
     * The FT_Conjunction class
     ***************************************************************************/

    template <class TClass> class FT_Conjunction : public TypedFeature<TClass> {
    protected:
        TypedFeature<TClass> *f2;
        vector<int> mapf1Vals2ArrayIndexes;

    public:
    FT_Conjunction(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFeature<TClass>* f2)
        : TypedFeature<TClass>(fInd, paramInd,
                               parsingFrames,
                               name,
                               fe, 0,
                               f2->valType()) {

            this->f2 = f2;

        }

    FT_Conjunction(
        int fInd,
        int paramInd,
        int parsingFrames,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        TypedFeature<TClass>* f2)
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               parsingFrames,
                               fargs,
                               offset,
                               fe, 0,
                               f2->valType()) {

            this->f2 = f2;

        }

        virtual int numPreCompEntries() = 0;
        virtual int numf1Params() = 0;
        virtual int f1Value(Tag *ge, int frame) = 0;


        inline FL_Signal<EdgeInst> *signal() {
            return f2->signal();
        }

        inline TypedFeature<TClass> *compFeature() {
            return f2;
        }

        inline bool isPreComputable() {
            return f2->isPreComputable();
        }

        /*
         * Returns the index of the last entry used in array for precomputation
         */
        inline int preComputeTo(bool collapseFrames, int initArrayIndex = 0) {
            this->_collapseFrames = collapseFrames;
            int rangeUsedIndexes = f2->preComputeTo(collapseFrames, initArrayIndex) -
                initArrayIndex + 1;

            for(int i = 0; i < numPreCompEntries(); i++) {
                mapf1Vals2ArrayIndexes.push_back(initArrayIndex);
                initArrayIndex += rangeUsedIndexes;
            }

            return initArrayIndex - 1;
        }


        virtual inline void doPreComputation(double **array,
                                             int preCompArrayIndex,
                                             int beg, int end, TStrand strand,
                                             const FeatureVector **params,
                                             int frame,
                                             int fConjInd = 0,
                                             TypedFilter<UCHAR> *filter = NULL) {

            if(!this->isPreComputable()) {
                assert(0);
                throw EXCEPTION(BAD_USAGE, this->getName() + " cannot be precomputed.\n");
            }

            int offset = fConjInd*numPreCompEntries();
            //      cerr << "offset for " << this->getName() << " " << offset << endl;

            for(int i = 0; i < numPreCompEntries(); i++) {
                int arrayIndex = preCompArrayIndex + mapf1Vals2ArrayIndexes[i];
                //        cerr << i << " arrayIndex " << arrayIndex << endl;

                f2->doPreComputation(array, arrayIndex, beg, end, strand,
                                     params, frame, offset + i, filter);
            }
        }

        inline void accumulatePreCompEntries(double **array,
                                             int preCompArrayIndex,
                                             int period, int beg, int end) {

            for(int i = 0; i < numPreCompEntries(); i++) {
                int arrayIndex = preCompArrayIndex + mapf1Vals2ArrayIndexes[i];
                f2->accumulatePreCompEntries(array, arrayIndex, period, beg, end);
            }
        }

        inline double preComputedDotParamV(double **array,
                                           int preCompArrayIndex,
                                           Tag *ge, int frame) {

            int arrayIndex = preCompArrayIndex +
                mapf1Vals2ArrayIndexes[f1Value(ge, frame)];

            return f2->preComputedDotParamV(array, arrayIndex, ge, frame);
        }


        virtual inline double dotParamV(Tag *ge, int frame,
                                        const FeatureVector **params,
                                        int fConjInd = 0,
                                        TypedFilter<UCHAR> *filter = NULL) {

            int offset = fConjInd*numf1Params();
            return f2->dotParamV(ge, frame, params,
                                 offset + f1Value(ge, frame), filter);
        }

        virtual inline void updParamV(double updVal,
                                      Tag *ge, int frame,
                                      FeatureVector **params,
                                      int fConjInd = 0,
                                      TypedFilter<UCHAR> *filter = NULL) {

            int offset = fConjInd*numf1Params();
            f2->updParamV(updVal, ge, frame, params,
                          offset + f1Value(ge, frame), filter);
        }

        ~FT_Conjunction() { }
    };



    /**
     * The FT_FeatureXFeature class computes the conjunction of two features. The
     * restriction on them is that the first feature must be an integer-valued
     * one, so that its values can be used as indexes which will define different
     * sets of values for the second.
     * this class only works if f1 is a feature returning integer values
     ***************************************************************************/

    template<class TClass1, class TClass2>
        class FT_FeatureXFeature : public FT_Conjunction<TClass2> {

    protected:
        TypedFeature<TClass1> *f1;

    public:
    FT_FeatureXFeature(
        int fInd,
        int paramInd,
        char *name,
        FilterEngine *fe,
        TypedFeature<TClass1> *f1,
        TypedFeature<TClass2> *f2
        )
        : FT_Conjunction<TClass2>(fInd, paramInd,
                                  f2->parsingFrames(),
                                  name, fe, f2) {

            this->f1 = f1;

            initNumParams();

        }

    FT_FeatureXFeature(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_Conjunction<TClass2>(
            fInd,
            paramInd,
            fte->getFeature(fargs[offset+3])->parsingFrames(),
            fargs, offset,
            fe,
            (TypedFeature<TClass2>*)fte->getFeature(fargs[offset+3])
            ) {

            this->f1 = (TypedFeature<UCHAR> *)fte->getFeature(fargs[offset++]);
            offset++;

            initNumParams();
        }

        inline int numPreCompEntries() {
            return f1->maxNumFeatValues();
        }

        inline int numf1Params() {
            return f1->maxNumFeatValues();
        }

        inline void initNumParams() {
            Pair<ULONG> & p1 = f1->maxNumParams();
            if(p1.f != 1)
                throw EXCEPTION(BAD_USAGE, f1->getName() + " cannot be used as f1 feature in conjuntion.\n");
            Pair<ULONG> p2 = this->f2->maxNumParams();

            p2.f = p1.s*p2.f;
            this->setMaxNumParams(p2);

        }

        inline int f1Value(Tag *ge, int frame) {
            return (int)f1->featValue(ge, frame);
        }

        inline double featValue(Tag *ge, int frame) {
            assert(0);
            throw EXCEPTION(BAD_USAGE, this->getName() + " is a multiple-valued feature.\n");
            return 0;
        }

        ~FT_FeatureXFeature() {

        }

    };

    /**
     * The FT_FilterXFeature class computes the conjunction between a filter
     * and a feature. The filter must be integer-valued, in similar fashion to
     * FT_FeaturexFeature
     ***************************************************************************/

    template<class TClass1, class TClass2>
        class FT_FilterXFeature : public FT_Conjunction<TClass2> {
    protected:
        TypedFilter<TClass1> *f1;

    public:
    FT_FilterXFeature(
        int fInd,
        int paramInd,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass1> *f1,
        TypedFeature<TClass2> *f2
        )
        : FT_Conjunction<TClass2>(
            fInd, paramInd, f2->parsingFrames(),
            name, fe,
            f2
            ) {
            this->f1 = f1;

            if(this->f1->valType() != FT_UCHAR) {
                assert(0);
                throw EXCEPTION(BAD_USAGE, this->getName() + " must be FT_UCHAR");
            }

            initNumParams();

        }

    FT_FilterXFeature(
        int fInd,
        int paramInd,
        std::vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_Conjunction<TClass2>(fInd, paramInd,
                                  fte->getFeature(fargs[offset + 3])->parsingFrames(),
                                  fargs, offset,
                                  fe,
                                  (TypedFeature<TClass2>*)fte->getFeature(fargs[offset+3])
            ) {

            this->f1 = (TypedFilter<UCHAR> *)fe->getFilter(fargs[offset++]);

            if(this->f1->valType() != FT_UCHAR) {
                assert(0);
                throw EXCEPTION(BAD_USAGE, this->getName() + " must be FT_UCHAR");
            }

            offset++;
            initNumParams();

        }

        inline void initNumParams() {
            Pair<ULONG> p = this->f2->maxNumParams();
            p.f = f1->maxNumFilterValues()*p.f;
            this->setMaxNumParams(p);
        }

        virtual inline int numPreCompEntries() {
            return f1->maxNumFilterValues();
        }

        inline int numf1Params() {
            return f1->maxNumFilterValues();
        }

        virtual inline int f1Value(Tag *ge, int frame) {
            return f1->value(ge->getPos(), ge->getStrand());
        }

        inline TypedFilter<UCHAR> *compFilter() {
            return f1;
        }

        inline double featValue(Tag *ge, int frame) {
            assert(0);
            throw EXCEPTION(BAD_USAGE, this->getName() + " is a multiple-valued feature.\n");
            return 0;
        }

        ~FT_FilterXFeature() {

        }

    };


    /**
     * FT_LazyFilterXFeature is a subclass of FT_FilterXFeature which does not
     * precompute f2 for all possible filter f1's values and for each position.
     * In this respect is more efficient FT_FilterXFeature as it delays f2
     * precomputation,  until f2 knows about the current value of f1.
     * f2 on the other hand uses f1's value as an offset to call dotParamV.
     *
     */

    template<class TClass1, class TClass2>
        class FT_LazyFilterXFeature : public FT_FilterXFeature<TClass1, TClass2> {

    protected:

    public:
    FT_LazyFilterXFeature(
        int fInd,
        int paramInd,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass1> *f1,
        TypedFeature<TClass2> *f2
        ) :
        FT_FilterXFeature<TClass1, TClass2>(fInd, paramInd,
                                            name, fe,
                                            f1, f2) {

        }

    FT_LazyFilterXFeature(
        int fInd,
        int paramInd,
        std::vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        ) :
        FT_FilterXFeature<TClass1, TClass2>(fInd, paramInd,
                                            fargs, offset,
                                            fe, fte) {

        }

        inline int numPreCompEntries() {
            return 1;
        }

        inline int f1Value(Tag *ge, int frame) {
            return 0;
        }

        inline void doPreComputation(double **array,
                                     int preCompArrayIndex,
                                     int beg, int end, TStrand strand,
                                     const FeatureVector **params,
                                     int frame,
                                     int fConjInd = 0,
                                     TypedFilter<UCHAR> *filter = NULL) {

            if(!this->isPreComputable() || !this->f2->isFeatureSet() || filter) {
                assert(0);
                throw EXCEPTION(BAD_USAGE, this->getName() + " cannot be lazy-precomputed.\n");
            }

            int offset = fConjInd*this->f1->maxNumFilterValues();
            int arrayIndex = preCompArrayIndex + this->mapf1Vals2ArrayIndexes[0];

            this->f2->doPreComputation(array, arrayIndex,
                                       beg, end, strand, params, frame,
                                       offset, this->f1);

        }

        inline double dotParamV(Tag *ge, int frame,
                                const FeatureVector **params,
                                int fConjInd = 0,
                                TypedFilter<UCHAR> *filter = NULL) {

            assert(!filter && this->f2->isFeatureSet());
            int offset = fConjInd*this->numf1Params();

            return this->f2->dotParamV(ge, frame, params, offset, this->f1);

        }

        inline void updParamV(double updVal,
                              Tag *ge, int frame,
                              FeatureVector **params, int fConjInd = 0,
                              TypedFilter<UCHAR> *filter = NULL) {

            assert(!filter && this->f2->isFeatureSet());
            int offset = fConjInd*this->numf1Params();

            this->f2->updParamV(updVal, ge, frame, params, offset, this->f1);

        }

        ~FT_LazyFilterXFeature() {}

    };

}

#endif
