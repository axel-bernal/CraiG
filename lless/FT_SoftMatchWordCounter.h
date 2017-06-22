/****************************************************************************
 * FT_SoftMatchWordCounter.h - part of the lless namespace, a general purpose
 *                             linear semi-markov structure prediction library
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

#ifndef _FT_SOFT_MATCH_COUNTS_H_
#define _FT_SOFT_MATCH_COUNTS_H_

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Feature.h"
#include "FT_Edge.h"
#include "FT_Bin.h"
#include "Sequence.h"
#include "NGram.h"
#include "Motif.h"


using namespace lless;

namespace lless {

    /**
     * The FT_SoftMatchWordCounter class is an array of bin features which use
     * soft matches of the Tag object against a contained DBSoftMatchMotif
     * resource database.
     * The binned match counts is what is used as feature values.
     ***************************************************************************/

    template <class TClass> class FT_SoftMatchWordCounter : public FT_Edge<TClass> {

    protected:
        int maxHD;
        int _left, _right;
        int len;
        FL_Gram<UCHAR>* gram;
        FT_Bin<int> **hDists;
        DBSoftMatchMotif *dbMatches;

    public:
    FT_SoftMatchWordCounter(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        int left,
        int right,
        FL_Signal<EdgeInst>* signal,
        TypedFilter<UCHAR>* gram,
        DBSoftMatchMotif *dbMatches,
        TValType type = FT_INTEGER)
        : FT_Edge<TClass>(fInd,
                          paramInd,
                          parsingFrames,
                          name, fe, 0,
                          signal, type) {

            this->dbMatches = dbMatches;
            _left = left;
            _right = right;
            this->gram = (FL_Gram<UCHAR> *)gram;

            initialize();

        }

    FT_SoftMatchWordCounter(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_Edge<TClass>(fInd,
                          paramInd,
                          fargs, offset,
                          fe, 0) {

            if(!sscanf(fargs[offset++].c_str(), "%d", &_left))
                assert(0);
            if(!sscanf(fargs[offset++].c_str(), "%d", &_right))
                assert(0);

            ResourceEngine *re = fe->getResourceEngine();
            dbMatches = (DBSoftMatchMotif *)re->getResource(fargs[offset++]);
            gram = (FL_Gram<UCHAR> *)fe->getFilter(fargs[offset++]);

            initialize();

            int numFeatsToAdd;

            if(!sscanf(fargs[offset++].c_str(), "%d", &numFeatsToAdd))
                assert(0);

            for(int i = 0; i < numFeatsToAdd; i++) {
                int hd;
                FT_Bin<int> *feature;

                if(!sscanf(fargs[offset++].c_str(), "%d", &hd))
                    assert(0);

                feature = (FT_Bin<int> *)fte->getFeature(fargs[offset++]);
                addHammingDistBin(hd, feature);
            }
        }

        void initialize() {
            maxHD = dbMatches->maxHammingDistance();
            hDists = new FT_Bin<int> * [maxHD + 1];
            len = _left + _right;

            assert(len == dbMatches->motifLength() &&
                   gram->alphabetSize() == dbMatches->alphabet()->alphabetSize());

            this->setMaxNumParams(maxHD + 1, 0);
        }

        void addHammingDistBin(int hd, FT_Bin<int> *dist) {
            assert(hd <= maxHD);
            hDists[hd] = ((FT_Bin<int> *)dist);

            if(dist->maxNumFeatValues() > this->_numParams.s)
                this->_numParams.s = dist->maxNumFeatValues();

        }

        inline double featValue(Tag *ge, int frame) {
            assert(0);
            return 0;
        }

        inline int length() {
            return len;
        }

        inline int left() {
            return _left;
        }

        inline int right() {
            return _right;
        }

        virtual inline int valueAtPos(Tag *ge, int i) {
            throw EXCEPTION( FORB_VIRTUAL_FUNCTION, string(" FT_SoftMatchWordCounter.h"));
        }
        virtual inline int valueAtPos(Tag *ge, int i, int *vals) {
            throw EXCEPTION( FORB_VIRTUAL_FUNCTION, string(" FT_SoftMatchWordCounter.h"));
        }

        inline double dotParamV(Tag *ge, int frame,
                                const FeatureVector **params,
                                int fConjInd = 0,
                                TypedFilter<UCHAR> *filter = NULL) {

            assert(!filter);

            int i, beg = ge->getPos() - left();
            int baseCode = 0;
            double result = 0;
            fConjInd = fConjInd*(maxHD + 1);
            int alphSize = this->gram->alphabetSize();

            for(i = 0; i < len; i++)
                baseCode = baseCode*alphSize + this->gram->value(beg + i, ge->getStrand());

            //      cerr << " soft counts for " << this->getName() << "@" << ge->getPos() << endl;

            for(i = 0; i <= maxHD; i++) {
                double hdRes = hDists[i]->dotParamV(dbMatches->softFreq(i, baseCode),
                                                    params, fConjInd + i);

                //        cerr << "hdist(" << i << ") = " << dbMatches->softFreq(i, baseCode) << " " << hdRes << endl;

                result += hdRes;

            }
            return result;
        }

        virtual inline void updParamV(double updVal, Tag *ge, int frame,
                                      FeatureVector **params, int fConjInd = 0,
                                      TypedFilter<UCHAR> *filter = NULL) {

            assert(!filter);

            int i, beg = ge->getPos() - left();
            int baseCode = 0;
            fConjInd = fConjInd*(maxHD + 1);
            int alphSize = this->gram->alphabetSize();
            for(i = 0; i < len; i++)
                baseCode = baseCode*alphSize + this->gram->value(beg + i, ge->getStrand());

            //    cerr << "fvals = " << dbMatches->softFreq(0, baseCode) << " " << dbMatches->(1, baseCode) << endl;

            for(i = 0; i <= maxHD; i++) {
                //      cerr << this->getName() << " HD = " << i << " " << ge->getPos() << " " << dbMatches->softFreq(i, baseCode) << endl;
                hDists[i]->updParamV(updVal, dbMatches->softFreq(i, baseCode), params, fConjInd + i);
            }
            //    for(i = 0; i < len; i++) {
            //      int ind = this->gram->value(beg + i, ge->getStrand());
            //      cerr << this->gram->alphabet()->contC(ind);
            //    }
            //    cerr << endl;
        }

        ~FT_SoftMatchWordCounter() {

            if(hDists)
                delete [] hDists;

            hDists = NULL;
        }

    };

}

#endif
