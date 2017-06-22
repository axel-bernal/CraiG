/****************************************************************************
 * FL_MovingQuantile.h - part of the lless namespace, a general purpose
 *                   linear semi-markov structure prediction library
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

#ifndef _FILTER_MOVINGQUANTILE_H_
#define _FILTER_MOVINGQUANTILE_H_

#include "Utils.h"
#include "Filter.h"
#include "FilterEngine.h"
#include <set>

namespace lless {

    /**
     *
     * FL_MovingQuantile is a subclass of FL_Score. The scoring function computes
     * the quantile computed around a window at each position
     */
    template<class TClass>
        class FL_MovingQuantile: public TypedFilter<TClass> {
    private:
        int sldwSize;
        int quantile;
        TypedFilter<TClass> *cov;
        multiset<double> lwHalf[NUM_STRANDS], upHalf[NUM_STRANDS];
    public:
        //! Default constructor
    FL_MovingQuantile(int fInd,                //!< A unique identifier.
                      std::string & name,      //!< A unique name.
                      int sldwSize,          /*!< A sliding window's size.
                                               If equal to zero, the window size
                                               is set automatically to the length
                                               of the current input sequence.
                                             */
                      int quantile,       //!< The quantile
                      TypedFilter<TClass> *cov
        ) :
        TypedFilter<TClass>(fInd, name, 1) {
            this->cov = cov;
            this->sldwSize = sldwSize;
            this->quantile = quantile;
        }

        /**
         * Constructor from a Header string definition
         */
    FL_MovingQuantile(int fInd,           //!< A unique identifier.
                      /*! The Header string definition, loaded as a vector of strings.
                       * The Header has the following form:\n\n
                       * Filter name MovinQuantile sldwSize quantile\n\n
                       */
		      vector<std::string> &params,
		      int & offset,       //!< The index for vector params.
		      ResourceEngine *re, //!< A pointer to the ResourceEngine
		                          // object.
		      FilterEngine *fe    //!< A pointer to the FilterEngine
		                          // object.
        ) :
        TypedFilter<TClass>(fInd, params, offset, 1) {

            this->cov = (TypedFilter<TClass> *)fe->getFilter(params[offset++]);
            assert(this->cov);
            if(!sscanf(params[offset++].c_str(), "%d", &sldwSize))
                assert(0);
            if(!sscanf(params[offset++].c_str(), "%d", &quantile))
                assert(0);

        }

        void freeValArrays(TStrand strand) {
            lwHalf[strand].clear();
            upHalf[strand].clear();
            TypedFilter<TClass>::freeValArrays(strand);
        }

        void setDefaultStrand(TStrand strand) {
            this->defaultStrand = strand;
            cov->setDefaultStrand(strand);
        }

        double popFirstElem(multiset<double> &S) {
            double result = *(S.begin());
            S.erase( S.begin() );
            return result;
        }

        double popLastElem(multiset<double> &S) {
            multiset<double>::iterator it = S.end();
            it--;
            double result = *it;
            S.erase( it );
            return result;
        }

        double getLastElem(multiset<double> &S) {
            multiset<double>::iterator it = S.end();
            it--;
            return *it;
        }

        /**
         * A member function for computing filter values.
         *
         * The value at each position of the input sequence is the context level
         * obtained through computation of the relative abundance of ctxSymbols
         * in a sliding window centered at said position.
         *
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        inline void computeVals(char *seq, TClass *filterVals, int len) {
            assert(len);
            int i = 0, j = -5;
            int wSize = sldwSize;
            TStrand strand = this->defaultStrand;

            if(!wSize || len < wSize)
                wSize = len;

            for( ; i <= len; i++) {
                // add the new value
                int K = (wSize+1)*quantile/100;
                int L = wSize*(100 - quantile)/100;
                double cval = cov->value(i);
                if (lwHalf[strand].size() == K && getLastElem(lwHalf[strand]) <= cval)
                    upHalf[strand].insert(cval);
                else
                    lwHalf[strand].insert(cval);

                // remove the old value
                if(i >= wSize) {
                    double ocval = cov->value(i - wSize);
                    if(ocval <= getLastElem(lwHalf[strand])) {
                        multiset<double>::iterator it = lwHalf[strand].find(ocval);
                        lwHalf[strand].erase(it);
                    }
                    else {
                        multiset<double>::iterator it = upHalf[strand].find(ocval);
                        upHalf[strand].erase(it);
                    }
                }

                // update the counts of elements in our multisets

                while(lwHalf[strand].size() > K)
                    upHalf[strand].insert(popLastElem(lwHalf[strand]));
                while (upHalf[strand].size() > L)
                    lwHalf[strand].insert(popFirstElem(upHalf[strand]));

                // compute the quantile
                if(i == wSize - 1) {
                    double elem = getLastElem(lwHalf[strand]);
                    for( ; j < wSize/2; j++) {
                        //	    cerr << "first " << j << " " << len << " " << i << endl;
                        filterVals[j] = elem;
                    }
                }
                else if(i == len) {
                    double elem = getLastElem(lwHalf[strand]);
                    for( ; j <= len + 5; j++)  {
                        //	    cerr << "last " << j << " " << len << " " << i << endl;
                        filterVals[j] = elem;
                    }
                }
                else if(i >= wSize) {
                    //	  cerr << "middle "<< j << " " << len << " " << i << endl;
                    filterVals[j++] = getLastElem(lwHalf[strand]);
                }
            }
            assert(j == len + 6);

        }

        ~FL_MovingQuantile() {
        }
    };

}

#endif
