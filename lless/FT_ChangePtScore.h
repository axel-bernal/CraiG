/****************************************************************************
 * FT_ChangePtScore.h - part of the lless namespace, a general purpose
 *              linear semi-markov structure prediction library
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

#ifndef _FT_CHANGE_PTSCORE_H_
#define _FT_CHANGE_PTSCORE_H_

#include "Utils.h"
#include "Filter.h"
#include "ChangePtUtils.h"

namespace lless {

    /**
     *
     * FT_ChangePtScore
     *
     ***************************************************************************/

    class FT_ChangePtScore : public TypedFeature<double> {
    protected:
        TypedFeature<double> *leftMean;
        TypedFilter<double> *leftOvMean;
        TypedFeature<double> *rightMean;
        TypedFilter<double> *rightOvMean;

    public:
    FT_ChangePtScore(
        int fInd,
        int paramInd,
        char *name,
        FilterEngine *fe,
        TypedFeature<double> *leftMean,
        TypedFilter<double> *leftOvMean,
        TypedFeature<double> *rightMean,
        TypedFilter<double> *rightOvMean
        )
        : TypedFeature<double>(fInd, paramInd,
                               leftMean->parsingFrames(), name,
                               fe, 1, FT_DOUBLE) {

            this->leftMean = leftMean;
            this->leftOvMean = leftOvMean;
            this->rightMean = rightMean;
            this->rightOvMean = rightOvMean;

        }

    FT_ChangePtScore(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : TypedFeature<double>(fInd, paramInd,
                               fte->getFeature(fargs[offset + 2])->parsingFrames(),
                               fargs, offset,
                               fe, 1, FT_DOUBLE) {

            this->leftMean = (TypedFeature<double> *)fte->getFeature(fargs[offset++]);
            this->leftOvMean = (TypedFilter<double> *)fe->getFilter(fargs[offset++]);
            this->rightMean = (TypedFeature<double> *)fte->getFeature(fargs[offset++]);
            this->rightOvMean = (TypedFilter<double> *)fe->getFilter(fargs[offset++]);
        }

        inline double featValue(Tag *ge, int frame) {
            float mean_diff =
                pow(leftMean->featValue(ge, frame) -
                    leftOvMean->value(ge->getPos(), ge->getStrand()), 2.0) +
                pow(rightMean->featValue(ge, frame) -
                    rightOvMean->value(ge->getPos(), ge->getStrand()), 2.0);

            return mean_diff;
        }

        inline double dotParamV(double featVal,
                                const FeatureVector **params,
                                int fConjInd = 0) {

            return featVal*(*params[fConjInd])[0];
        }

        inline void updParamV(double updVal, double featVal,
                              FeatureVector **params, int fConjInd = 0) {

            (*params[fConjInd])[0] += featVal*updVal;
        }

        ~FT_ChangePtScore() {

        }
    };
}

#endif
