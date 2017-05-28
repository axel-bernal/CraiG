/****************************************************************************
 * FL_SpacerMotifScore.h - part of the lless namespace, a general purpose
 *                         linear semi-markov structure prediction library
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

#ifndef _FILTER_SPACERMOTIFSCORE_H_
#define _FILTER_SPACERMOTIFSCORE_H_

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
     * FL_SpacerMotifScore is a subclass of FL_Score. It is similar to
     * FL_MotifScore, but here, the score at each sequence position depends
     * on the match strength against two different Motif objects (motif profiles)
     * and the separating distance between them.
     */
    class FL_SpacerMotifScore : public FL_Score<double> {
    protected:
        Motif *motif1, *motif2;
        double threshold1, threshold2;
        int rangeMin, rangeMax;

    public:
    FL_SpacerMotifScore(int fInd,              //!< A unique identifier.
                        std::string & name,    //!< A unique name.
                        Motif *motif1,         //!< First motif profile
                        Motif *motif2,         //!< Second motif profile
                        double threshold1,     //!< Threshold motif profile1
                        double threshold2,     //!< Threshold motif profile2
                        int rangeMin,          //!< min distance between motifs
                        int rangeMax           //!< max distance between motifs
        )
        :  FL_Score<double>(fInd,
                            name,
                            1,
                            0, FT_DOUBLE
            ) {

            this->motif1 = motif1;
            this->threshold1 = threshold1;
            this->motif2 = motif2;
            this->threshold2 = threshold2;
            this->rangeMin = rangeMin;
            this->rangeMax = rangeMax;

        }

        /**
         * Constructor from a Header string definition
         */
    FL_SpacerMotifScore(int fInd,        //!< A unique identifier.
                        /*! The Header string definition, loaded as a vector
                         * of strings.
                         * The Header has the following form:\n\n
                         * Filter name SpacerMotifScore arrInd period
                         * smoothenWSize threshold1 motif1 threshold2 motif2
                         * minRange maxRange \n\n
                         * The description of the fields could be found in
                         * the other constructor(s)
                         */
                        vector<std::string> & params,
                        int & offset,       //!< The index for vector params.
                        ResourceEngine *re, /*!< A pointer to the
                                              ResourceEngine object. */
                        FilterEngine *fe    /*!< A pointer to the FilterEngine
                                              object. */
        )
        : FL_Score<double>(fInd,
                           params,
                           offset,
                           fe,
                           FT_DOUBLE) {

            if(!sscanf(params[offset++].c_str(), "%lf", &threshold1))
                assert(0);
            this->motif1 = (Motif *)re->getResource(params[offset++]);

            if(!sscanf(params[offset++].c_str(), "%lf", &threshold2))
                assert(0);
            this->motif2 = (Motif *)re->getResource(params[offset++]);

            assert(this->motif1 && this->motif2);

            if(!sscanf(params[offset++].c_str(), "%d", &rangeMin))
                assert(0);
            if(!sscanf(params[offset++].c_str(), "%d", &rangeMax))
                assert(0);

        }

        /**
         * The implemented scoring function. At each position p the score for
         * motif1 is computed, then the function checks if motif2 appears
         * between positions p + rangeMin and p + rangeMax, if it does the
         * final score is the sum of the two motif scores, of zero otherwise
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        void scoreInputSeq(char *seq, double *filterVals, int len) {
            //      cerr << this->getName() << endl;
            int i, j;
            double m2Score = 0, res;

            for(i = 0; i < len - rangeMax; i++)
                filterVals[i + 1] = motif1->score(threshold1, seq + i);

            for(i = 0; i < len - rangeMax; i++) {
                if(!filterVals[i + 1])
                    continue;

                for(j = i + rangeMin; j < i + rangeMax; j++) {
                    res =  motif2->score(threshold2, seq + j);

                    if(res > m2Score)
                        m2Score = res;
                }

                filterVals[i + 1] += m2Score;

                if(m2Score != 0)
                    filterVals[i + 1] *= 2;
                //        cerr << i << " " << filterVals[i + 1] << " " << m2Score << "\n";
            }
        }

        ~FL_SpacerMotifScore() {

        }

    };

}

#endif
