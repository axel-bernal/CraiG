/****************************************************************************
 * FT_PWMxPWM.h - part of the lless namespace, a general purpose
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

#ifndef _FT_PWMxPWM_FEAT_H_
#define _FT_PWMxPWM_FEAT_H_
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "FT_PWM.h"


namespace lless {


    /**
     * FT_PWMxFT_PWM is a subclass of FT_Edge<int> for computing the
     * feature conjunction between two FT_PWM features.
     * It can prove to be very useful for finding correlations between signal
     * positions, which might not be able to be captured by a WWAM model
     * alone.
     **************************************************************************/

    template <class TClass>
        class FT_PWMxPWM : public FT_Edge<TClass> {

    protected:
        FT_WWAM<TClass> *pwm1, *pwm2;

    public:
    FT_PWMxPWM(
        int fInd,
        int paramInd,
        char *name,
        FilterEngine *fe,
        FT_WWAM<TClass> *pwm1,
        FT_WWAM<TClass> *pwm2)
        : FT_Edge<TClass>(fInd, paramInd,
                          Utils::max(pwm1->parsingFrames(),
                                     pwm2->parsingFrames()),
                          name,
                          fe,
                          pwm2->maxNumFeatValues(),
                          pwm1->signal(),
                          pwm1->valType()) {

            this->pwm1 = pwm1;
            this->pwm2 = pwm2;
            assert(pwm1->step() == pwm2->step() == 1);
            this->setMaxNumParams(this->pwm1->maxNumFeatValues(),
                                  this->pwm2->maxNumFeatValues());

        }

    FT_PWMxPWM(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        ) :
        FT_Edge<TClass>(
            fInd,
            paramInd,
            Utils::max(fte->getFeature(fargs[offset + 2])->parsingFrames(),
                       fte->getFeature(fargs[offset + 3])->parsingFrames()),
            fargs, offset, fe, 0,
            fte->getFeature(fargs[offset + 3])->valType(),
            ((FT_Edge<TClass> *)fte->getFeature(fargs[offset + 3]))->signal()
            ) {

            pwm1 = (FT_WWAM<TClass> *)fte->getFeature(fargs[offset++]);
            pwm2 = (FT_WWAM<TClass> *)fte->getFeature(fargs[offset++]);
            assert(pwm1->step() == pwm2->step() == 1);

            this->setMaxNumParams(pwm1->maxNumFeatValues(),
                                  pwm2->maxNumFeatValues());

        }

        ~FT_PWMxPWM() {

        }

        inline double featValue(Tag *ge, int frame) {
            assert(0);
            throw EXCEPTION(BAD_USAGE, this->getName() + " is a multiple-valued feature.\n");
            return 0;
        }

        inline double dotParamV(Tag *ge, int frame,
                                const FeatureVector **params,
                                int fConjInd = 0,
                                TypedFilter<UCHAR> *filter = NULL) {

            assert(!filter);

            double result = 0;
            register int i, j;
            int end1 = pwm1->numUpdates() - pwm1->order();
            int end2 = pwm2->numUpdates();
            int offset = fConjInd*pwm1->maxNumFeatValues();
            /*
             * Pointer comparison. If objects are the same
             * then no need to correlate all positions
             */
            if(pwm1 == pwm2) {
                int *values = pwm1->computeValues(ge);

                for(i = 0; i < end1; i++) {
                    assert(offset + values[i] < pwm1->maxNumFeatValues());
                    const FeatureVector &mpar = *params[offset + values[i]];

                    for(j =  i + 1; j < end2; j++) {
                        assert(values[j] < pwm2->maxNumFeatValues());
                        result += mpar[values[j]];
                    }
                }
                /*        TClass *values = pwm1->computeValues(ge);

                          for(i = 0; i <= end1; i += pwm1->step()) {
                          int offset_i = offset + i*pwm1->domainSize() + values[i];
                          j =  i + pwm2->step();

                          for( ; j <= end2; j += pwm2->step())
                          result += (*params[offset_i])[j*pwm2->domainSize() + values[j]];
                          }
                */
            }
            else {
                int *values1 = pwm1->computeValues(ge);
                int *values2 = pwm2->computeValues(ge);

                for(i = 0; i < end1; i++) {
                    assert(offset + values1[i] < pwm1->maxNumFeatValues());
                    const FeatureVector &mpar = *params[offset + values1[i]];

                    for(j = 0; j < end2; j++) {
                        assert(values2[j] < pwm2->maxNumFeatValues());
                        result += mpar[values2[j]];
                    }
                }
            }

            return result;
        }

        inline void updParamV(double updVal, Tag *ge, int frame,
                              FeatureVector **params, int fConjInd = 0,
                              TypedFilter<UCHAR> *filter = NULL) {

            assert(!filter);

            int i, j;
            int end1 = pwm1->numUpdates() - pwm1->order();
            int end2 = pwm2->numUpdates();
            int offset = fConjInd*pwm1->maxNumFeatValues();

            /*
             * Pointer comparisson. If objects are the same
             * then no need to correlate all positions
             */
            if(pwm1 == pwm2) {
                int *values = pwm1->computeValues(ge);

                for(i = 0; i < end1; i++) {
                    assert(offset + values[i] < pwm1->maxNumFeatValues());
                    FeatureVector &mpar = *params[offset + values[i]];

                    for(j =  i + 1; j < end2; j++) {
                        assert(values[j] < pwm2->maxNumFeatValues());
                        mpar[values[j]] += updVal;
                    }
                }
                /*        for(i = 0; i <= end1; i += step_1) {
                          int offset_i = offset + i*pwm1->domainSize() + values[i];
                          j = i + pwm2->step();

                          for( ; j <= end2; j += step_2)
                          (*params[offset_i])[j*pwm2->domainSize() + values[j]] += updVal;
                          }*/
            }
            else {
                int *values1 = pwm1->computeValues(ge);
                int *values2 = pwm2->computeValues(ge);

                for(i = 0; i < end1; i++) {
                    assert(offset + values1[i] < pwm1->maxNumFeatValues());
                    FeatureVector &mpar = *params[offset + values1[i]];

                    for(j = 0; j < end2; j++) {
                        assert(values2[j] < pwm2->maxNumFeatValues());
                        mpar[values2[j]] += updVal;
                    }
                }
            }
        }
    };

    /**
     * FT_PWMxFeature is a subclass of FT_Edge<int> for computing the
     * feature conjunction between a FT_PWM and any other feature.
     * It can prove to be very useful for finding correlations between signal
     * positions and other feature values.
     **************************************************************************/

    template <class TClass1, class TClass2>
        class FT_PWMxFeature : public FT_Edge<TClass2> {

    protected:
        FT_WWAM<TClass1> *pwm1;
        TypedFeature<TClass2> *f2;

    public:
    FT_PWMxFeature(
        int fInd,
        int paramInd,
        char *name,
        FilterEngine *fe,
        FT_WWAM<TClass1> *pwm1,
        TypedFeature<TClass2> *f2)
        : FT_Edge<TClass2>(fInd, paramInd,
                           Utils::max(pwm1->parsingFrames(),
                                      f2->parsingFrames()),
                           name,
                           fe,
                           f2->maxNumFeatValues(),
                           pwm1->signal(),
                           f2->valType()) {

            this->pwm1 = pwm1;
            this->f2 = f2;
            this->setMaxNumParams(this->pwm1->maxNumFeatValues(),
                                  this->f2->maxNumFeatValues());

        }

    FT_PWMxFeature(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        ) :
        FT_Edge<TClass2>(
            fInd,
            paramInd,
            Utils::max(fte->getFeature(fargs[offset + 2])->parsingFrames(),
                       fte->getFeature(fargs[offset + 3])->parsingFrames()),
            fargs, offset, fe, 0,
            fte->getFeature(fargs[offset + 3])->valType(),
            ((FT_Edge<TClass1> *)fte->getFeature(fargs[offset + 2]))->signal()
            ) {

            pwm1 = (FT_WWAM<TClass1> *)fte->getFeature(fargs[offset++]);
            f2 = (TypedFeature<TClass2> *)fte->getFeature(fargs[offset++]);

            this->setMaxNumParams(pwm1->maxNumFeatValues(),
                                  f2->maxNumFeatValues());

        }

        ~FT_PWMxFeature() {

        }

        inline double featValue(Tag *ge, int frame) {
            assert(0);
            throw EXCEPTION(BAD_USAGE, this->getName() + " is a multiple-valued feature.\n");
            return 0;
        }

        inline double dotParamV(Tag *ge, int frame,
                                const FeatureVector **params,
                                int fConjInd = 0,
                                TypedFilter<UCHAR> *filter = NULL) {

            assert(!filter);

            double result = 0;
            register int i;
            int end1 = pwm1->numUpdates() - pwm1->order();
            int offset = fConjInd*pwm1->maxNumFeatValues();


            int *values = pwm1->computeValues(ge);

            for(i = 0; i < end1; i++) {
                assert(offset + values[i] < pwm1->maxNumFeatValues());
                assert(values[i] < pwm1->maxNumFeatValues());
                result += f2->dotParamV(ge, frame, params,
                                        offset + values[i]);
            }

            return result;
        }

        inline void updParamV(double updVal, Tag *ge, int frame,
                              FeatureVector **params, int fConjInd = 0,
                              TypedFilter<UCHAR> *filter = NULL) {

            assert(!filter);

            int i;
            int end1 = pwm1->numUpdates() - pwm1->order();
            int offset = fConjInd*pwm1->maxNumFeatValues();


            int *values = pwm1->computeValues(ge);

            for(i = 0; i < end1; i++) {
                assert(offset + values[i] < pwm1->maxNumFeatValues());
                assert(values[i] < pwm1->maxNumFeatValues());
                f2->updParamV(updVal, ge, frame, params,
                              offset + values[i]);
            }
        }
    };

}

#endif
