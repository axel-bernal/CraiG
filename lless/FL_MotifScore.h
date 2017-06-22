/****************************************************************************
 * FL_MotifScore.h - part of the lless namespace, a general purpose
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

#ifndef _FILTER_MOTIFSCORE_H_
#define _FILTER_MOTIFSCORE_H_

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
     * FL_MotifScore is a subclass of FL_Score. The scoring function matches
     * a given motif profile against each position of the sequence. The score
     * reflects the strength of the match.
     */
    class FL_MotifScore : public FL_Score<double> {

    private:
        Motif *motif;
        double threshold;

    public:
        //! Default constructor
    FL_MotifScore(int fInd,             //!< A unique identifier.
                  std::string & name,   //!< A unique name.
                  double threshold,     /*!< A motif match score below this
                                          threshold will be replaced by zero
                                        */
                  Motif *motif          //!< A motif profile
        )
        :  FL_Score<double>(fInd,
                            name,
                            1,
                            0, FT_DOUBLE
            ) {

            this->threshold = threshold;
            this->motif = motif;

        }

        /**
         * Constructor from a Header string definition
         */
    FL_MotifScore(int fInd,            //!< A unique identifier.
                  /*! The Header string definition, loaded as a vector of
                   * strings.
                   * The Header has the following form:\n\n
                   * Filter name MotifScore arrInd period smoothenWSize
                   * threshold motif\n\n
                   * The description of the fields could be found in the
                   * other constructor(s)
                   */
                  vector<std::string> & params,
                  int & offset,       //!< The index for vector params.
                  ResourceEngine *re, /*!< A pointer to the ResourceEngine
                                        object. */
                  FilterEngine *fe    /*!< A pointer to the FilterEngine
                                        object. */
        )
        : FL_Score<double>(fInd,
                           params,
                           offset,
                           fe,
                           FT_DOUBLE) {

            if(!sscanf(params[offset++].c_str(), "%lf", &threshold))
                assert(0);

            this->motif = (Motif *)re->getResource(params[offset++]);
            assert(this->motif);
        }

        /**
         * The implemented scoring function. At each position p the score is the
         * match strength of the fragmeng seq[p..p+length(motif)]  against the
         * motif profile, if the score is above threshold, or zero otherwise.
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        void scoreInputSeq(char *seq, double *filterVals, int len) {
            //      cerr << this->getName() << endl;
            for(int i = 0; i < len; i++) {
                filterVals[i + 1] = motif->score(threshold, seq + i);
                //        if(filterVals[i + 1])
                //          cerr << i << " " << filterVals[i + 1] << "\n";
            }
        }

        ~FL_MotifScore() {

        }

    };

}

#endif
