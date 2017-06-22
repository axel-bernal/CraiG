/****************************************************************************
 * FT_Composition.h - part of the lless namespace, a general purpose
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

#ifndef _FT_COMPOSITION_H_
#define _FT_COMPOSITION_H_

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Feature.h"

namespace lless {

    /***************************************************************************
     * The FT_Composition class
     ***************************************************************************/
    /*
     * \todo add precomputation, then features like BinnedCountingSegment, WWAMs and
     * EdgeFeature can be implemented in terms of features compositions
     */

    template <class TClass> class FT_Composition : public TypedFeature<TClass> {
    protected:
        TypedFeature<TClass> *f2;

    public:
    FT_Composition(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        Feature* f2)
        : TypedFeature<TClass>(fInd, paramInd,
                               parsingFrames,
                               name, fe,
                               f2->maxNumFeatValues(), f2->valType()) {

            this->f2 = (TypedFeature<TClass> *)f2;

        }

    FT_Composition(
        int fInd,
        int paramInd,
        int parsingFrames,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        Feature* f2)
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               parsingFrames,
                               fargs,
                               offset,
                               fe,
                               f2->maxNumFeatValues(), f2->valType()) {

            this->f2 = (TypedFeature<TClass> *)f2;

        }

        /**
         * A member function that performs dot product for feature precomputation.
         * The precomputation occurs at position pos, the  word found
         * at such position is used to choose the right entry in array *params.
         */
        inline double dotParamVPreComp(int pos, int frame,
                                       TStrand strand,
                                       const FeatureVector **params,
                                       int fConjInd = 0) {

            throw EXCEPTION( NOT_SUPPORTED,
                             this->getName() + std::string(" cannot precompute"));


        }

        inline int fConjIndPreComp(int pos,
                                   TStrand strand,
                                   int fConjInd = 0,
                                   TypedFilter<UCHAR> *filter = NULL) {

            throw EXCEPTION( NOT_SUPPORTED,
                             this->getName() + std::string(" cannot precompute"));

        }


        virtual inline double dotParamV(Tag *ge, int frame,
                                        const FeatureVector **params,
                                        int fConjInd = 0,
                                        TypedFilter<UCHAR> *filter = NULL) {

            double featVal = this->featValue(ge, frame);
            //      cerr << this->getName() << " " << ge->getPos() << "_" << ge->getLen() << "(" << featVal << ")" << endl;

            return f2->dotParamV(featVal, params, fConjInd);

        }

        virtual inline void updParamV(double updVal, double featVal,
                                      FeatureVector **params,
                                      int fConjInd = 0) {
            //      cerr << this->getName() << "(" << updVal << " " << featVal << ")" << endl;
            f2->updParamV(updVal, featVal, params, fConjInd);

        }

        ~FT_Composition() { }
    };

    /*
     * \todo what is going on with the precomputation routines?
     */
    template <class TClass> class FT_FeatureOFeature : public FT_Composition<TClass> {
    protected:
        Feature *f1;

    public:
    FT_FeatureOFeature(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        Feature* f1,
        Feature* f2)
        : FT_Composition<TClass>(fInd, paramInd,
                                 f1->parsingFrames(),
                                 name, fe, f2) {

            this->f1 = f1;

        }

    FT_FeatureOFeature(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte)
        : FT_Composition<TClass>(fInd,
                                 paramInd,
                                 fte->getFeature(fargs[offset + 2])->parsingFrames(),
                                 fargs,
                                 offset,
                                 fe,
                                 fte->getFeature(fargs[offset + 3])) {

            this->f1 = fte->getFeature(fargs[offset++]);
            offset++;

        }

        inline double featValue(Tag *ge, int frame) {
            return f1->featValue(ge, frame);
        }

        /**
         * @param array the precomputation array.
         * @param preCompArrayIndex index in the precomputation.
         * array in which the precomputed dot product values will be stored.
         * @param ge Tag object to which this feature is tied.
         * @param frame Tag object's current phase(or frame).
         * @return the dot product of feature and its corresponding parameter(s)
         *
         */
        inline double preComputedDotParamV(double **array,
                                           int preCompArrayIndex,
                                           Tag *ge, int frame) {

            throw EXCEPTION( NOT_SUPPORTED,
                             this->getName() + std::string(" cannot precompute"));

        }

        inline void accumulatePreCompEntries(double **array,
                                             int preCompArrayIndex,
                                             int period, int beg,
                                             int end) {

            throw EXCEPTION( NOT_SUPPORTED,
                             this->getName() + std::string(" cannot precompute"));

        }

        ~FT_FeatureOFeature() {}

    };


    template <class TClass1, class TClass2>
        class FT_FilterOFeature : public FT_Composition<TClass2> {
    protected:
        TypedFilter<TClass1> *f1;

    public:
    FT_FilterOFeature(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass1>* f1,
        Feature* f2)
        : FT_Composition<TClass2>(fInd, paramInd, 1,
                                  name, fe, f2) {

            this->f1 = f1;

        }

    FT_FilterOFeature(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte)
        : FT_Composition<TClass2>(fInd,
                                  paramInd, 1,
                                  fargs,
                                  offset,
                                  fe,
                                  fte->getFeature(fargs[offset + 3])) {

            this->f1 = (TypedFilter<TClass1> *)fe->getFilter(fargs[offset++]);
            offset++;

        }

        inline double featValue(Tag *ge, int frame) {
            return (double)f1->value(ge->getPos(), ge->getStrand());
        }


        ~FT_FilterOFeature() {}

    };

    /**
     * This class describes a composition between features at the update level
     * the updVal parameter which is used by f2 corresponds to the value of
     * feature f1
     */
    template <class TClass> class FT_FeatureO2Feature : public FT_Composition<TClass> {
    protected:
        Feature *f1;

    public:
    FT_FeatureO2Feature(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        Feature* f1,
        Feature* f2)
        : FT_Composition<TClass>(fInd, paramInd, 1,
                                 name, fe, f2) {

            this->f1 = f1;

        }

    FT_FeatureO2Feature(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte)
        : FT_Composition<TClass>(fInd,
                                 paramInd, 1,
                                 fargs,
                                 offset,
                                 fe,
                                 fte->getFeature(fargs[offset + 3])) {

            this->f1 = fte->getFeature(fargs[offset++]);
            offset++;

        }

        inline double featValue(Tag *ge, int frame) {
            return this->f2->featValue(ge, frame);
        }

        inline double dotParamV(Tag *ge, int frame,
                                const FeatureVector **params,
                                int fConjInd = 0,
                                TypedFilter<UCHAR> *filter = NULL) {

            double featVal = this->featValue(ge, frame);
            double res = f1->featValue(ge,frame)*this->f2->dotParamV(featVal, params, fConjInd);
            //      cerr << this->getName() << " " << ge->getPos() << "_" << ge->getLen() << " " << featVal << " (" << res << ")" << endl;
            return res;
        }

        inline void updParamV(double updVal, Tag *ge, int frame,
                              FeatureVector **params,
                              int fConjInd = 0,
                              TypedFilter<UCHAR> *filter = NULL) {

            updVal *= f1->featValue(ge, frame);
            double featVal = this->featValue(ge, frame);
            //      cerr << this->getName() << " " << updVal << " " << featVal << endl;
            this->f2->updParamV(updVal, featVal, params, fConjInd);

        }

        ~FT_FeatureO2Feature() {}

    };

}

#endif
