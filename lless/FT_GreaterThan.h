/****************************************************************************
 * FT_GreaterThan.h - part of the lless namespace, a general purpose
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

#ifndef _FT_GREATER_THAN_H_
#define _FT_GREATER_THAN_H_

#include "Feature.h"
#include "FeatureEngine.h"


namespace lless {

    /**
     * The FT_GreaterThan class is a binary feature which is equivalent the
     * predicate [f->featValue(..) > threshold]
     ***************************************************************************/

    template <class TClass> class FT_GreaterThan : public TypedFeature<TClass> {

    private:
        double threshold;
        double weight;
        TypedFeature<TClass> *f;

    public:
    FT_GreaterThan(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFeature<TClass> *f,
        double threshold,
        double weight
        )
        : TypedFeature<TClass>(
            fInd, paramInd,
            parsingFrames, name,
            fe, 1
            ) {

            this->threshold = threshold;
            this->weight = weight;
            this->f = f;
        }

    FT_GreaterThan(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : TypedFeature<TClass>(
            fInd, paramInd,
            fargs, offset,
            fe, 1
            ) {

            if(!sscanf(fargs[offset++].c_str(), "%lf", &threshold))
                assert(0);
            if(!sscanf(fargs[offset++].c_str(), "%lf", &weight))
                assert(0);

            this->f = (TypedFeature<TClass> *)fte->getFeature(fargs[offset++]);
        }

        inline void turnOff() {
            this->threshold = DBL_MAX;
        }

        inline double featValue(Tag *ge, int frame) {

            return weight*(f->featValue(ge, frame) > threshold);
        }

        inline double dotParamV(double featVal,
                                const FeatureVector **params,
                                int fConjInd = 0) {


            return featVal*(*params[fConjInd])[0];
        }

        inline void updParamV(double updVal, double featVal,
                              FeatureVector **params, int fConjInd = 0) {

            //      cerr << this->getName() << " " << featVal <<  endl;
            (*params[fConjInd])[0] += featVal*updVal;
        }

        ~FT_GreaterThan() {

        }

    };

}

#endif
