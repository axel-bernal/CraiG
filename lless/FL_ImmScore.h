/****************************************************************************
 * FL_ImmScore.h - part of the lless namespace, a general purpose
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

#ifndef _FILTER_IMMSCORE_H_
#define _FILTER_IMMSCORE_H_

#include "Utils.h"
#include "FL_Score.h"
#include "ContextIMM.h"
#include "FilterEngine.h"
#include "ResourceEngine.h"
#include "IMM.h"
#include "Motif.h"

namespace lless {

    /**
     *
     * FL_ImmScore is a subclass of FL_Score whose scoring function depends on
     * cimm scores on the input sequence.
     */
    class FL_ImmScore: public FL_Score<double> {
    private:
        TypedFilter<UCHAR> *contextFilter;
        ContextIMM *cimm;

    public:
    FL_ImmScore(int fInd,                //!< A unique identifier.
                std::string & name,      //!< A unique name.
                int period,              //!< Filter values period.
                int smoothenWSize,       //!< Window size for smoothing scores.
                TypedFilter<UCHAR> *contextFilter, /*!< A pointer to a context
                                                     object */
                ContextIMM *cimm         //!< A pointer to a cimm object
        )
        :  FL_Score<double>(fInd,
                            name,
                            period,
                            smoothenWSize,
                            FT_DOUBLE
            ) {

            this->contextFilter = contextFilter;
            this->cimm = cimm;
        }

        /**
         * Constructor from a Header string definition
         */
    FL_ImmScore(
        int fInd,        //!< A unique identifier.
        /*! The Header string definition, loaded as a vector
         * of strings.
         * The Header has the following form:\n\n
         * Filter name ImmScore  arrInd period smoothenWSize
         * contextFilter cimmFile\n\n
         * If cimmFile equals "null", then the cimm Contents are to
         * be read right after the Header. The description of the
         * fields could be found in the other constructor(s)
         */
        vector<std::string> & params,
        int & offset,       //!< The index for vector params.
        ResourceEngine *re, /*!< A pointer to the
                              ResourceEngine object. */
        FilterEngine *fe    //!< A pointer to the FilterEngine object.
        )
        : FL_Score<double>(fInd,
                           params,
                           offset,
                           fe,
                           FT_DOUBLE) {

            this->contextFilter = (TypedFilter<UCHAR> *)fe->getFilter(params[offset++]);
            this->cimm = (ContextIMM *)re->getResource(params[offset++]);
        }

        inline void setDefaultStrand(TStrand strand) {
            this->defaultStrand = strand;
            contextFilter->setDefaultStrand(strand);
        }

        /**
         * The implemented scoring function. At each position p the score is the
         * output of the cimm object evaluated at p. The context filter decides
         * which imm to use.
         * For that, the function divides the sequence in fragments of size length
         * each of them having the same context filter value.
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        void scoreInputSeq(char *seq, double *filterVals, int len) {
            IMM **imm = cimm->getIMMModels();
            int i, j, phase;
            int length;

            for(int pos = 0; pos < len; pos += length) {
                length = 2;

                while(1) {

                    if(pos + length >= len)
                        break;

                    if(contextFilter->value(pos + length) != contextFilter->value(pos + length - 1))
                        break;

                    length++;
                }

                double *scores = new double [length + 100];
                assert(scores);

                for(j = 0; j < length + 100; j++)
                    scores[j] = 0;

                for(phase = 1; phase <= 3; phase++) {
                    imm[contextFilter->value(pos + 1)]->seqScore(scores, seq + pos + phase - 1, 1, 1 + pos % period(), length - phase + 1);

                    for(i = pos; i < pos + length - phase + 1;  i += 3) {
                        double score = 0;

                        for(int ph = 0; ph < 3; ph++)
                            score += scores[i - pos + ph]; // scores are taken as logs

                        filterVals[i + pos % 3 + phase] = score;
                    }
                }

                delete [] scores;
            }
        }

        ~FL_ImmScore() {}
    };

}

#endif
