/****************************************************************************
 * FT_PWM.h - part of the lless namespace, a general purpose
 *            linear semi-markov structure prediction library
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

#ifndef _FT_PWM_FEAT_H_
#define _FT_PWM_FEAT_H_
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "FT_WWAM.h"


namespace lless {


    /**
     * The FT_PWM class is a subtype of FT_WWAM with window size equal to 1.
     * They have been implemented as separate classes
     * given that they are commonly used in the literature.
     * @see FT_WWAM
     ***************************************************************************/
    template <class TClass>
        class FT_PWM : public FT_WWAM<TClass> {

    public:
    FT_PWM(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        FL_Signal<EdgeInst> *signal,
        TValType type,
        int domSize,
        int left,
        int right,
        int step)
        : FT_WWAM<TClass>(fInd,
                          paramInd, parsingFrames,
                          name, fe,
                          signal, type,
                          1,
                          domSize,
                          left, right,
                          step) {

            ;
        }

    FT_PWM(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte,
        int domSize
        )
        : FT_WWAM<TClass>(fInd,
                          paramInd,
                          fargs, offset,
                          fe,
                          1, domSize) {

        }

        ~FT_PWM() {

        }
    };


    template <class TClass1, class TClass2>
        class FT_BinnedPWM : public FT_BinnedWWAM<TClass1, TClass2> {

    public:
    FT_BinnedPWM(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        FL_Signal<EdgeInst> *signal,
        TValType type,
        int domSize,
        int left,
        int right,
        int step,
        FT_BaseBin<TClass1> *binObj)
        : FT_BinnedWWAM<TClass1, TClass2>(fInd,
                                          paramInd, parsingFrames,
                                          name, fe,
                                          signal, type,
                                          1,
                                          domSize,
                                          left, right,
                                          step, binObj) {

            ;
        }

    FT_BinnedPWM(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte,
        int domSize
        )
        : FT_BinnedWWAM<TClass1, TClass2>(fInd,
                                          paramInd,
                                          fargs, offset,
                                          fe, fte,
                                          1, domSize) {

        }

        ~FT_BinnedPWM() {

        }
    };


    /**
     * FT_GramPWM is a FT_GramWWAM with window size equal to 1.
     * @see FT_GramWWAM
     **************************************************************************/
    template <class TClass1, class TClass2>
        class FT_GramPWM : public FT_GramWWAM<TClass1, TClass2> {

    public:
    FT_GramPWM(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        FL_Signal<EdgeInst> *signal,
        FL_BaseGram<TClass1,TClass2> *gram,
        int left,
        int right,
        int step)
        : FT_GramWWAM<TClass1, TClass2>(fInd, paramInd,
                                        parsingFrames, name,
                                        fe, signal, 1, gram,
                                        left, right, step) {
            ;
        }

    FT_GramPWM(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_GramWWAM<TClass1, TClass2>(fInd,
                                        paramInd,
                                        fargs, offset,
                                        fe,
                                        1) {

        }

        ~FT_GramPWM() {

        }
    };

    /**
     * FT_BinnedGramPWM is a FT_BinnedGramWWAM with window size equal to 1.
     * @see FT_BinnedGramWWAM
     **************************************************************************/
    template <class TClass1, class TClass2, class TClass3>
        class FT_BinnedGramPWM : public FT_BinnedGramWWAM<TClass1, TClass2, TClass3> {

    public:
    FT_BinnedGramPWM(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        FL_Signal<EdgeInst> *signal,
        FL_FastGram<TClass2, SPARSE_HASH<TClass3,int> > *gram,
        int left,
        int right,
        int step,
        FT_BaseBin<TClass1> *binObj)
        : FT_BinnedGramWWAM<TClass1, TClass2, TClass3>(fInd, paramInd,
                                                       parsingFrames, name,
                                                       fe, signal, 1, gram,
                                                       left, right, step,
                                                       binObj) {
            ;
        }

    FT_BinnedGramPWM(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_BinnedGramWWAM<TClass1, TClass2, TClass3>(fInd,
                                                       paramInd,
                                                       fargs, offset,
                                                       fe, fte,
                                                       1) {

        }

        ~FT_BinnedGramPWM() {

        }
    };


    /**
     * FT_ConsensusPWM is a F_ConsensusWWAM with window size equal to 1.
     * @see FT_ConsensusWWAM
     **************************************************************************/

    template <class TClass>
        class FT_ConsensusPWM : public FT_ConsensusWWAM<TClass> {

    public:
    FT_ConsensusPWM(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        FL_Signal<EdgeInst> *signal,
        TypedFilter<TClass> *gram,
        char *consensus,
        int left,
        int right,
        int step)
        : FT_ConsensusWWAM<TClass>(fInd, paramInd,
                                   parsingFrames, name,
                                   fe, signal, 1,
                                   gram, consensus,
                                   left, right, step) {

            ;
        }

    FT_ConsensusPWM(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_ConsensusWWAM<TClass>(fInd,
                                   paramInd,
                                   fargs, offset,
                                   fe,
                                   1) {

        }

        ~FT_ConsensusPWM() {

        }
    };

    template <class TClass>
        class FT_PWMUnion : public FT_WWAMUnion<TClass> {

    public:
    FT_PWMUnion(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        FL_Signal<EdgeInst> *signal,
        int step,
        vector<FT_WWAM<TClass> *> & wwams)
        : FT_WWAMUnion<TClass>(fInd,
                               paramInd, parsingFrames,
                               name, fe,
                               signal,
                               1, step,
                               wwams) {
            ;
        }

    FT_PWMUnion(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_WWAMUnion<TClass>(fInd,
                               paramInd,
                               fargs, offset,
                               fe, fte,
                               1) {

        }

        ~FT_PWMUnion() {

        }
    };

}

#endif
