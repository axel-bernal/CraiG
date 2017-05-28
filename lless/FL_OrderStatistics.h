/****************************************************************************
 * FL_OrderStatistics.h - part of the lless namespace, a general purpose
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

#ifndef _FILTER_ORDER_STATISTICS_H_
#define _FILTER_ORDER_STATISTICS_H_

#include "Utils.h"
#include "Filter.h"
#include "ContextIMM.h"
#include "FilterEngine.h"
#include "ResourceEngine.h"
#include "IMM.h"
#include "Motif.h"

namespace lless {

    /**
     *
     * FL_OrderStatistics is a subclass of Filter which computes the 0th order
     * statistic dynamically and very efficiently over intervals. Other order
     * statistics are not implemented. This filter accumulates
     *
     ***************************************************************************/

    class FL_OrderStatistics : public TypedFilter<vector<int> > {
    protected:
        vector<int> orders;
        double *currA[NUM_STRANDS];
    public:
        //! Default constructor
    FL_OrderStatistics(
        int fInd,       //!< A unique identifier.
        std::string & name,       //!< A unique name.
        int order, //!< The n^th statistic -- 1 allowed only
        TValType type = FT_INTVECTOR //!< Filter value type.
        ) :
        TypedFilter<vector<int> >(fInd,
                                  name,
                                  1, FT_INTVECTOR) {

            assert(order == 1);
            for(int i = 0; i < NUM_STRANDS; currA[i++] = NULL);
            orders.push_back(order);
        }

        /**
         * First constructor from a Header string definition
         */
    FL_OrderStatistics(
        int fInd,         //!< A unique identifier.
        vector<std::string> &params, //!< Header string definition
        int & offset,       //!< The index for vector params.
        ResourceEngine *re, /*!< A pointer to the ResourceEngine
                              object. */
        FilterEngine *fe    /*!< A pointer to the FilterEngine
                              object. */
        )
        : TypedFilter<vector<int> >(fInd,
                                    params,
                                    offset,
                                    1, FT_INTVECTOR) {

            int order;
            if(!sscanf(params[offset++].c_str(), "%d", &order))
                assert(0);
            assert(order == 1);
            for(int i = 0; i < NUM_STRANDS; currA[i++] = NULL);
            orders.push_back(order);

        }

        /**
         * A member function for computing filter values.
         *
         * It calls the implementation-dependent scoring function for computing
         * filter values, smoothens these values if smoothenWSize is not zero and
         * finally, it accumulates the values along the input sequence.
         * @param seq the input sequence.
         * @param M the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        inline void computeVals(char *seq, vector<int> *M, int len) {
            double *A = currA[this->defaultStrand] = (double *)seq;

            int logN = log((double)len)/log(2.0), i,j;

            for(i = 0; i < len; i++)
                M[i] = vector<int>(logN + 1);

            // Generating Table of M[n][log(n)]

            for(i = 0; i < len; i++)
                M[i][0] = i;
            for(j = 1; (1<<j) <= len; j++)
                for (i = 0; i+(1<<j)-1 < len; i++)
                    M[i][j] = A[M[i][j-1]] <= A[M[i+(1<<(j-1))][j-1]]?
                        M[i][j-1] :
                        M[i+(1<<(j-1))][j-1];
        }

        double RMQ(Tag *tag) {
            TStrand strand = tag->getStrand();
            int i = tag->getPos() - 1, j = tag->getPos() + tag->getLen() - 2;
            double *A = currA[strand];
            vector<int> *M = this->values(strand);

            int diff=j-i;
            diff = 31 - __builtin_clz(diff+1);
            return A[M[i][diff]] <=A[M[j-(1<<diff)+1][diff]] ?
                M[i][diff] :
                M[j-(1<<diff)+1][diff];

        }
        ~FL_OrderStatistics() {

        }

    };

}

#endif
