/****************************************************************************
 * NGram.h - part of the lless namespace, a general purpose
 *           linear semi-markov structure prediction library
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

#ifndef _NGRAM_H
#define _NGRAM_H

#include <streambuf>
#include <string>
#include <string.h>
#include "Utils.h"
#include "Sigma.h"
#include "FL_Gram.h"
#include <math.h>
#include <assert.h>

#define MAX_ORDER 8
#define MIN_FREQ 1e-50
#define MIN_PROB exp((double)-DBL_MAX_EXP)
#define MAX_GRAM_ORDER 9

namespace lless {

    /**
     * NGram groups FL_Gram objects together and perform operations which need
     * all FL_Gram objects, such as interpolation and computation of backoff
     * probabilities.
     */
    class NGram {
    private:
        Sigma *alphabet;
        FL_Gram<int> *grams[MAX_GRAM_ORDER];
        int _maxOrder;

        void _initialize() {
            for(int i = 0; i < MAX_GRAM_ORDER; i++)
                grams[i] = NULL;
        }

    public:
        NGram(NGram &);

        //! Minimal constructor
        NGram(int maxOrder) {
            _initialize();
            _maxOrder = maxOrder;
        }
        NGram(Sigma *alphabet, int maxOrder = 3);
        NGram(::ifstream *, Sigma *alphabet);

        ~NGram();
        inline int maxOrder() {
            return _maxOrder;
        }

        inline int alphabetSize() {
            return grams[0]->alphabetSize();
        }

        FL_Gram<int>* gram(int ord);

        /**
         * A member function that deletes a contained FL_Gram object.
         * @param order the order of the FL_Gram to be removed.
         */
        inline void removeFL_Gram(int order) {
            if(grams[order - 1]) {
                delete grams[order - 1];
                grams[order - 1] = NULL;
            }
        }

        void generateSeqKmers(char *seq, int **seqFL_Grams, int len);

        void storeParams(::ofstream **);
        NGram & operator+=(NGram & ng);

        /**
         * @see FL_Gram::buildCp()
         */
        inline void buildCp() {
            for(int i = 0; i < _maxOrder; i++)
                if(grams[i])
                    grams[i]->buildCp();
        }
    };

}

#endif
