/****************************************************************************
 * FL_ChangePoint.h - part of the lless namespace, a general purpose
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

#ifndef _FL_CHANGE_POINT_H_
#define _FL_CHANGE_POINT_H_

#include "Utils.h"
#include "Filter.h"
#include "ChangePtUtils.h"

namespace lless {

    /**
     *
     * FL_ChangePoint is a subclass of Filter which computes the 0th order
     * statistic dynamically and very efficiently over intervals. Other order
     * statistics are not implemented. This filter accumulates
     *
     ***************************************************************************/

    template <class TClass>
        class FL_ChangePoint : public TypedFilter<ChangePoint> {
    protected:
        int block_len;
        double min_conf;
        int num_bootstraps;
        TypedFilter<TClass> *filter;
        pair<TClass, TClass> _range;
        vector<ChangePoint> pchanges[NUM_STRANDS];
    public:
        //! Default constructor
    FL_ChangePoint(
        int fInd,       //!< A unique identifier.
        std::string & name,       //!< A unique name.
        int block_len,            //!< block length
        double min_conf,
        int num_bootstraps,
        TypedFilter<TClass> *filter      //!< filter to detect cps
        ) :
        TypedFilter<ChangePoint>(fInd,
                                 name,
                                 1, FT_CHANGEPT) {

            this->block_len = block_len;
            this->min_conf = min_conf;
            this->num_bootstraps = num_bootstraps;
            this->filter = filter;
        }

        /**
         * First constructor from a Header string definition
         */
    FL_ChangePoint(
        int fInd,         //!< A unique identifier.
        vector<std::string> &params, //!< Header string definition
        int & offset,       //!< The index for vector params.
        ResourceEngine *re, /*!< A pointer to the ResourceEngine
                              object. */
        FilterEngine *fe    /*!< A pointer to the FilterEngine
                              object. */
        )
        : TypedFilter<ChangePoint>(fInd,
                                   params,
                                   offset,
                                   1, FT_CHANGEPT) {

            if(!sscanf(params[offset++].c_str(), "%d", &this->block_len))
                assert(0);
            if(!sscanf(params[offset++].c_str(), "%lf", &this->min_conf))
                assert(0);
            if(!sscanf(params[offset++].c_str(), "%d", &this->num_bootstraps))
                assert(0);
            this->filter = (TypedFilter<TClass> *)fe->getFilter(params[offset++]);

        }

        vector<ChangePoint> & getChangePoints(TStrand strand) {
            return pchanges[strand];
        }

        pair<TClass, TClass> &range() {
            return _range;
        }

        void freeValArrays(TStrand strand) {
            pchanges[strand].clear();
        }

        /**
         * A member function for computing filter values.
         *
         * It calls the implementation-dependent scoring function for computing
         * filter values, smoothens these values if smoothenWSize is not zero and
         * finally, it accumulates the values along the input sequence.
         * @param seq the input sequence.
         * @param vals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        inline void computeVals(char *seq, ChangePoint *vals, int len) {
            vector<double> tmp_cov(len + 40, 0.1);
            TClass *fvals = filter->values(this->defaultStrand);
            for(int i = 1; i < len; i++)
                tmp_cov[i + 20] = fvals[i];

            _range = Utils::findMinMax(tmp_cov, 20, len + 20);

            int blen = block_len;

            if(len/block_len < 25)
                blen = len/25 + 1;

            ChangePtUtils::computeChangePoints(tmp_cov,
                                               pchanges[this->defaultStrand],
                                               num_bootstraps,
                                               0, len,
                                               blen, min_conf);

            ChangePtUtils::refineChangePoints(tmp_cov,
                                              pchanges[this->defaultStrand],
                                              num_bootstraps,
                                              0, len,
                                              blen, min_conf);


            for(int i = 0; i < pchanges[this->defaultStrand].size(); i++)
                pchanges[this->defaultStrand][i].pos -= 20;

        }

        ~FL_ChangePoint() {

        }

    };

}


#endif
