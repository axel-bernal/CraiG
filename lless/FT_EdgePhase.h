/****************************************************************************
 * FT_EdgePhase.h - part of the lless namespace, a general purpose
 *                  linear semi-markov structure prediction library
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

#ifndef _FT_EDGE_PHASE_H_
#define _FT_EDGE_PHASE_H_
#include "FT_Edge.h"
#include "Filter.h"
#include "FilterEngine.h"
#include "FeatureEngine.h"

namespace lless {

    /**
     * FT_EdgePhase is a subclass of FT_Edge whose value is the phase
     * of the Tag object received as parameter.
     */

    class FT_EdgePhase : public FT_Edge<UCHAR> {
    public:
        /**
         * Default constructor
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param parsingFrames The maximum number of phases(frames) of any Tag
         * object to which this feature is tied to.
         * @param name A unique feature identifier.
         * @param fe A pointer to a FilterEngine object.
         * @param maxNumFeatVals The number of possible feature values.
         * If greater than zero then the Feature object's possible values are
         * countable.
         * @param signal the signal this feature is associated with.
         */
    FT_EdgePhase(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        int maxNumFeatVals,
        FL_Signal<EdgeInst> *signal
        )
        : FT_Edge<UCHAR>(fInd,
                         paramInd,
                         parsingFrames,
                         name, fe,
                         maxNumFeatVals,
                         signal,
                         FT_UCHAR) {
        }

        /**
         * Constructor from a Header string definition.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param fargs The Header string definition loaded as a vector of strings.
         * The Header has the following form:\n\n
         * Feature name EdgePhase parsingFrames signal maxNumFeatVals\n\n
         * The description of the fields could be found in
         * the other constructor(s)
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param fte A pointer to a FeatureEngine object.
         */
    FT_EdgePhase(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_Edge<UCHAR>(fInd,
                         paramInd,
                         fargs,
                         offset,
                         fe,
                         FT_UCHAR) {

        }

        /**
         * @param ge Tag object to which this feature is tied.
         * @param frame Tag object's current phase(or frame).
         * @return ge's phase
         */
        inline double featValue(Tag *ge, int frame) {
            return frame;
        }

        ~FT_EdgePhase() {}

    };

    class FT_NodePhase : public TypedFeature<UCHAR> {
    public:
    FT_NodePhase(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        int maxNumFeatVals
        )
        : TypedFeature<UCHAR>(fInd,
                              paramInd,
                              parsingFrames,
                              name, fe,
                              maxNumFeatVals,
                              FT_UCHAR) {
        }

        /**
         * Constructor from a Header string definition.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param fargs The Header string definition loaded as a vector of strings.
         * The Header has the following form:\n\n
         * Feature name NodePhase parsingFrames maxNumFeatVals\n\n
         * The description of the fields could be found in
         * the other constructor(s)
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param fte A pointer to a FeatureEngine object.
         */
    FT_NodePhase(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : TypedFeature<UCHAR>(fInd,
                              paramInd,
                              fargs,
                              offset,
                              fe,
                              0,
                              FT_UCHAR) {

            if(!sscanf(fargs[offset++].c_str(), "%ld", &this->_numParams.s))
                assert(0);

        }

        /**
         * @param ge Tag object to which this feature is tied.
         * @param frame Tag object's current phase(or frame).
         * @return ge's phase
         */
        inline double featValue(Tag *ge, int frame) {
            return frame;
        }

        ~FT_NodePhase() {}

    };

}

#endif
