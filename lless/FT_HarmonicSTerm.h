/****************************************************************************
 * FT_HarmonicSTerm.h - part of the lless namespace, a general purpose
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

#ifndef _FT_HARMONIC_STERM_H_
#define _FT_HARMONIC_STERM_H_

#include "Feature.h"
#include "FeatureEngine.h"


namespace lless {

    /**
     * The FT_HarmonicSTerm class is a feature whose value is the meanLength of
     * the tag object (given as parameter) divided by the current's tag object's
     * length.
     ***************************************************************************/

    class FT_HarmonicSTerm : public TypedFeature<double> {
        double meanLength;
    private:
        TypedFeature<int> *ft;

    public:
    FT_HarmonicSTerm(
        int fInd,
        int paramInd,
        char *name,
        double meanLength,
        FilterEngine *fe,
        TypedFeature<int> *ft)
        : TypedFeature<double>(fInd, paramInd,
                               ft->parsingFrames(), name,
                               fe, 1, FT_DOUBLE) {

            this->ft = ft;
            this->meanLength = meanLength;
        }

    FT_HarmonicSTerm(
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

            this->ft = (TypedFeature<int> *)fte->getFeature(fargs[offset++]);
            if(!sscanf(fargs[offset++].c_str(), "%lf", &meanLength))
                assert(0);

        }

        inline double featValue(Tag *ge, int frame) {
            return meanLength/ft->featValue(ge, frame);
            //   return exp(-ft->featValue(ge, frame));
        }

        inline double dotParamV(double featVal,
                                const FeatureVector **params, int fConjInd = 0) {

            return featVal*(*params[fConjInd])[0];
        }

        inline void updParamV(double updVal, double featVal,
                              FeatureVector **params, int fConjInd = 0) {

            (*params[fConjInd])[0] += featVal*updVal;
        }

        ~FT_HarmonicSTerm() {

        }

    };

}

#endif
