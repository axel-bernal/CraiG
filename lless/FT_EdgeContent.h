/****************************************************************************
 * FT_EdgeContent.h - part of the lless namespace, a general purpose
 *                    linear semi-markov structure prediction library
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

#ifndef _FT_EDGE_CONTENT_H_
#define _FT_EDGE_CONTENT_H_
#include "FT_Edge.h"
#include "Filter.h"
#include "FilterEngine.h"
#include "FeatureEngine.h"

namespace lless {

    /**
     * FT_EdgeContent is subclass of FT_Edge which measures the content
     * before(to the right) and after(to the left) the occurrence of a signal
     * and computes its difference as feature value.
     * The content measure is given by a content feature which is a parameter
     * of the class.
     */
    template<class TClass> class FT_EdgeContent : public FT_Edge<TClass> {

    protected:
        TypedFeature<TClass> *contentFeature;
        int upsRegion, dwsRegion;
        int nextNodePos;

    public:
        /**
         * Default constructor.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param name A unique feature identifier.
         * @param fe A pointer to a FilterEngine object.
         * @param signal the signal this feature is associated with.
         * @param nextNodePos the relative position of the next node to the
         * location of the signal.
         * @param contentFeature the feature used to measure content around
         * the signal.
         * @param upsRegion the size of the region to the right of the signal.
         * @param dwsRegion the size of the region to the left of the signal.
         * @param type the feature value type, one of TValType enumerate type
         */
    FT_EdgeContent(
        int fInd,
        int paramInd,
        char *name,
        FilterEngine *fe,
        FL_Signal<EdgeInst> *signal,
        int nextNodePos,
        TypedFeature<TClass> *contentFeature,
        int upsRegion,
        int dwsRegion,
        TValType type = FT_INTEGER)
        : FT_Edge<TClass>(fInd,
                          paramInd,
                          contentFeature->parsingFrames(),
                          name,
                          fe, 2,
                          signal, type) {

            this->upsRegion = upsRegion;
            this->dwsRegion = dwsRegion;
            this->contentFeature = contentFeature;
            this->nextNodePos = nextNodePos;

        }

        /**
         * Constructor from a Header string definition.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * object to which this feature is tied to.
         * @param fargs The Header string definition loaded as a vector of strings.
         * The Header has the following form:\n\n
         * Feature name EdgeContent<TClass> signal nextNodePos contentFeature
         * upsRegion dwsRegion\n\n
         * The description of the fields could be found in
         * the other constructor(s)
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param fte A pointer to a FeatureEngine object.
         */
    FT_EdgeContent(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_Edge<TClass>(fInd,
                          paramInd,
                          fte->getFeature(fargs[offset+4])->parsingFrames(),
                          fargs,
                          offset,
                          fe,
                          2) {

            if(!sscanf(fargs[offset++].c_str(), "%d", &nextNodePos))
                assert(0);

            this->contentFeature = (TypedFeature<TClass> *)fte->getFeature(fargs[offset++]);

            if(!sscanf(fargs[offset++].c_str(), "%d", &upsRegion))
                assert(0);
            if(!sscanf(fargs[offset++].c_str(), "%d", &dwsRegion))
                assert(0);

        }

        /**
         * @param ge Tag object to which this feature is tied.
         * @param frame Tag object's current phase(or frame).
         * @return the difference in content, as measured by contentFeature,
         * in the region before and the region after ge->getPos()
         */
        inline double featValue(Tag *ge, int frame) {
            int bPos = ge->getPos() + nextNodePos;
            /*
             * creating auxiliary states
             */
            NodeInst rightBS(NO_NODE_INST,
                             INVALID_NODE,
                             bPos,
                             ge->getStrand(),
                             dwsRegion + 1);
            NodeInst leftBS(NO_NODE_INST,
                            INVALID_NODE,
                            bPos - upsRegion,
                            ge->getStrand(),
                            upsRegion + 1);

            double upsRes = contentFeature->featValue(&leftBS, frame);
            double dwsRes = contentFeature->featValue(&rightBS, frame);

            //      cerr << this->getName() << " ups " << upsRes << " ds " << dwsRes << endl;

            if(!upsRegion)
                return dwsRes;

            if(!dwsRegion)
                return upsRes;

            return dwsRes - upsRes;
        }

        /**
         * @param featVal the feature value, computed with
         * featValue(Tag *, int). It is upcasted to double at this
         * point since the result of the dot product is also double.
         * @param params parameter vector to be used for computing the dot
         * product.
         * @param fConjInd is an offset for the array **params, which is
         * different from zero when this feature has been used in a Feature
         * conjunction or Feature disjunction operation (defined as the class
         * FeaturexFeature, in which case fConjInd usually denotes the value of
         * the first feature; in the latter case, the values of the first feature
         * must countable.
         * @return The dot product between feature value and the feature's
         * corresponding parameter(s).
         */
        inline double dotParamV(double featVal,
                                const FeatureVector **params,
                                int fConjInd = 0) {

            return featVal*(*params[fConjInd])[featVal > 0 ? 1 : 0];
        }

        /**
         * @param updVal the value of the update to be made to this feature's
         * corresponding parameter(s)
         * @param featVal the feature value, computed with
         * featValue(Tag *, int). It is upcasted to double at this
         * point since the result of the dot product is also double.
         * @param params parameter vector to be used for computing the dot
         * product.
         * @param fConjInd is an offset for the array **params, which is
         * different from zero when this feature has been used in a Feature
         * conjunction or Feature disjunction operation (defined as the class
         * FeaturexFeature, in which case fConjInd usually denotes the value of
         * the first feature; in the latter case, the values of the first feature
         * must countable.
         */
        inline void updParamV(double updVal, double featVal,
                              FeatureVector **params, int fConjInd = 0) {

            (*params[fConjInd])[featVal > 0 ? 1 : 0] += featVal*updVal;
        }

        ~FT_EdgeContent() {

        }

    };

}

#endif
