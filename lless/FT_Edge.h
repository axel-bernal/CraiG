/****************************************************************************
 * FT_Edge.h - part of the lless namespace, a general purpose
 *             linear semi-markov structure prediction library
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

#ifndef _FT_EDGE__H_
#define _FT_EDGE__H_
#include "Feature.h"
#include "Filter.h"
#include "FilterEngine.h"
#include "FeatureEngine.h"

namespace lless {

    /**
     * The FT_Edge class is the base class for all features that are tied
     * to edges which have a signal object associated with them.
     * There are many constructors from Header definitions available for this
     * class, depending on the presence or absence of the following parameters:
     * parsingFrames, numFeatVals, signal and Type. Derived classes need this
     * many constructors.
     *
     ***************************************************************************/

    template <class TClass> class FT_Edge: public TypedFeature<TClass> {

    protected:
        //! signal type the edge feature is associated with
        FL_Signal<EdgeInst> *_signal;

    public:
        /**
         * Default constructor.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param name A unique feature identifier.
         * @param fe A pointer to a FilterEngine object.
         * @param maxNumFeatVals The number of possible feature values.
         * If greater than zero then the Feature object's possible values are
         * countable.
         * @param signal the signal this feature is associated with.
         * @param type the feature value type, one of TValType enumerate type
         */
    FT_Edge(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        int maxNumFeatVals,
        FL_Signal<EdgeInst> *signal,
        TValType type = FT_INTEGER
        )
        : TypedFeature<TClass>(fInd, paramInd,
                               parsingFrames,
                               name, fe,
                               maxNumFeatVals, type) {

            this->_signal = signal;
            assert(signal);
        }

        /**
         * Constructor from a Header string definition. The feature's value type,
         * the parsingFrames and maxNumFeatVals need to be determined from
         * information in the Header.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param fargs The Header string definition loaded as a vector of strings.
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param fte A pointer to a FeatureEngine object.
         */
    FT_Edge(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               fargs,
                               offset,
                               fe, 0) {

            extractSignal(fargs, offset);

            if(!sscanf(fargs[offset++].c_str(), "%d", &this->_numParams.s))
                assert(0);

        }

        /**
         * Constructor from a Header string definition. The parsingFrames and
         * maxNumFeatVals need to be determined from information in the Header.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param fargs The Header string definition loaded as a vector of strings.
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param type the feature value type, one of TValType enumerate type
         */
    FT_Edge(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        TValType type
        )
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               fargs,
                               offset,
                               fe,
                               0, type) {

            extractSignal(fargs, offset);

            if(!sscanf(fargs[offset++].c_str(), "%d", &this->_numParams.s))
                assert(0);

        }

        /**
         * Constructor from a Header string definition. The parsingFrames needs
         * to be determined from information in the Header.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param fargs The Header string definition loaded as a vector of strings.
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param maxNumFeatVals The number of possible feature values.
         * If greater than zero then the Feature object's possible values are
         * countable.
         * @param type the feature value type, one of TValType enumerate type
         */
    FT_Edge(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        int maxNumFeatValues,
        TValType type
        )
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               fargs,
                               offset,
                               fe,
                               maxNumFeatValues,
                               type) {

            extractSignal(fargs, offset);

        }


        /**
         * Constructor from a Header string definition. The feature's value type
         * and parsingFrames need to be determined from information in the Header.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param fargs The Header string definition loaded as a vector of strings.
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param maxNumFeatVals The number of possible feature values.
         * If greater than zero then the Feature object's possible values are
         * countable.
         */
    FT_Edge(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        int maxNumFeatVals
        )
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               fargs,
                               offset,
                               fe,
                               maxNumFeatVals) {

            extractSignal(fargs, offset);

        }

        /**
         * Constructor from a Header string definition. The feature's value type
         * needs to be determined from information in the Header.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param parsingFrames The maximum number of phases(frames) of any Tag
         * object to which this feature is tied to.
         * @param fargs The Header string definition loaded as a vector of strings.
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param maxNumFeatVals The number of possible feature values.
         * If greater than zero then the Feature object's possible values are
         * countable.
         */
    FT_Edge(
        int fInd,
        int paramInd,
        int parsingFrames,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        int maxNumFeatVals
        )
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               parsingFrames,
                               fargs,
                               offset,
                               fe,
                               maxNumFeatVals) {

            extractSignal(fargs, offset);

        }

        /**
         * Constructor from a Header string definition. The feature's value type
         * needs to be determined from information in the Header.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param parsingFrames The maximum number of phases(frames) of any Tag
         * object to which this feature is tied to.
         * @param fargs The Header string definition loaded as a vector of strings.
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param maxNumFeatVals The number of possible feature values.
         * If greater than zero then the Feature object's possible values are
         * countable.
         * @param type the feature value type, one of TValType enumerate type
         * @param signal the signal this feature is associated with.
         */
    FT_Edge(
        int fInd,
        int paramInd,
        int parsingFrames,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        int maxNumFeatVals,
        TValType type,
        FL_Signal<EdgeInst> *signal
        )
        : TypedFeature<TClass>(fInd,
                               paramInd,
                               parsingFrames,
                               fargs,
                               offset,
                               fe,
                               maxNumFeatVals,
                               type) {

            _signal = signal;

        }

        /**
         * A member function for extracting the type information from the
         * Header definition
         */
        void extractSignal(vector<std::string> & fargs, int & offset) {
            this->_signal = NULL;

            if(fargs[offset++].compare("null")) {
                _signal = (FL_Signal<EdgeInst> *)this->fe->getFilter(fargs[offset - 1]);
                assert(_signal);
            }
        }

        inline FL_Signal<EdgeInst> *signal() {
            return _signal;
        }

        /**
         * @param array the precomputation array.
         * @param preCompArrayIndex index in the precomputation.
         * array in which the precomputed dot product values will be stored.
         * @param ge Tag object to which this feature is tied.
         * @param frame Tag object's current phase(or frame).
         * @return the dot product of feature and its corresponding parameter(s)
         *
         */
        virtual inline double preComputedDotParamV(double **array,
                                                   int preCompArrayIndex,
                                                   Tag *ge,
                                                   int frame) {

            return array[preCompArrayIndex + this->_preCompArrayIndex[frame]][ge->getPos()];

        }


        /**
         * A member function that performs dot product for feature precomputation.
         * The precomputation occurs at position pos and frame, the signal found
         * at such position must be of right type. The precomputation is only is
         * allowed when TClass = EdgeInst
         */
        inline double dotParamVPreComp(int pos, int frame,
                                       TStrand strand,
                                       const FeatureVector **params,
                                       int fConjInd = 0,
                                       TypedFilter<UCHAR> *filter = NULL) {

            assert(_signal);
            double result = 0;

            if(strand == STRAND_COMP)
                pos = this->fe->seqLength() - pos + 1;

            EdgeInst *ge = &_signal->value(pos, strand);
            int offset = fConjInd;

            if(filter && !this->isFeatureSet()) {
                offset += filter->value(pos, strand);
                filter = NULL;
            }

            if(_signal->type() == (int)ge->getType()) {
                if(strand == STRAND_COMP)
                    ge->setPos(this->fe->seqLength() - ge->getPos() + 1);

                result = this->dotParamV(ge, frame, params, offset, filter);

                if(strand == STRAND_COMP)
                    ge->setPos(this->fe->seqLength() - ge->getPos() + 1);
            }

            return result;

        }

        ~FT_Edge() {

        }

    };

}

#endif
