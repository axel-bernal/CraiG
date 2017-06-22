/****************************************************************************
 * FT_GramWWAM.h - part of the lless namespace, a general purpose
 *                 linear semi-markov structure prediction library
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

#ifndef FT_GRAMWWAM_FEAT_H
#define FT_GRAMWWAM_FEAT_H
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "FT_WWAM.h"
#include "NGram.h"
#include "Consensus.h"
#include "FT_Edge.h"


namespace lless {

    /**
     * FT_GramWWAM is a subtype of FT_WWAM that uses gram values at each signal
     * position.
     **************************************************************************/
    template <class TClass1, class TClass2>
        class FT_GramWWAM : public FT_WWAM<TClass2> {

    protected:
        FL_BaseGram<TClass1,TClass2> *gram;

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
         * @param signal the signal this feature is associated with.
         * @param windowSize the size of the window around any position in the
         * WWAM which is used for counting symbols or words; windowSize/size
         * must be an odd number.
         * @param gram a pointer to a FL_Gram object which defines the words that
         * will be counted at position.
         * @param offset the length of the WWAM, offset of the signal occurrence.
         * @param length  the length of the WWAM, length of the signal occurrence.
         * The total length of the WWAM would be then offset + length.
         * @param step the number of positions to jump within the WWAM, before
         * counting again.
         */
    FT_GramWWAM(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        FL_Signal<EdgeInst> *signal,
        TValType type,
        int windowSize,
        FL_BaseGram<TClass1, TClass2> *gram,
        int offset,
        int length,
        int step
        )
        : FT_WWAM<TClass2>(fInd,
                           paramInd,
                           parsingFrames,
                           name,
                           fe,
                           signal,
                           type,
                           windowSize,
                           gram->maxNumFilterValues(),
                           offset,
                           length,
                           step) {

            this->gram = gram;
            initialize();
        }

        /**
         * Constructor from a Header string definition. windowSize should be read
         * from the Header definition.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * object to which this feature is tied to.
         * @param fargs The Header string definition loaded as a vector of strings.
         * The Header has the following form:\n\n
         * Feature name FT_GramWWAM parsingFrames signal offset length step windowSize
         * gram \n\n
         * The description of the fields could be found in
         * the other constructor(s)
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param fte A pointer to a FeatureEngine object.
         */
    FT_GramWWAM(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_WWAM<TClass2>(fInd,
                           paramInd,
                           fargs, offset,
                           fe,
                           fe->getFilter(fargs[offset + 8])->maxNumFilterValues()) {


            this->gram = (FL_BaseGram<TClass1,TClass2> *)fe->getFilter(fargs[offset++]);
            initialize();

        }

    FT_GramWWAM(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        int windowSize
        )
        : FT_WWAM<TClass2>(fInd,
                           paramInd,
                           fargs, offset,
                           fe,
                           windowSize,
                           fe->getFilter(fargs[offset + 7])->maxNumFilterValues()) {

            this->gram = (FL_BaseGram<TClass1,TClass2> *)fe->getFilter(fargs[offset++]);
            initialize();

        }

        void initialize() {
            this->_order = (int)gram->order();
            this->_numPosVals = this->windowSize() - this->order() + 1;
            this->_numUpdates = ((this->length() - 1)/this->step() + 1)*this->numPosVals();
            this->_numParams.s = ((this->length() - 1)/this->step() + 1)*this->domainSize();
            //      this->_numParams.s = (_offset + _length + 1)*domSize;

            this->vals = new int [this->numUpdates()];
        }


        /**
         * A member function that computes protected member vals. Any entry
         * appearing at some position within array vals is assigned the
         * gram appearing at that position.
         * The positions in the sequence that are covered by vals are between
         * ge->getPos() + offset and ge->getPos() + offset + length
         * @return vals
         */
        inline int* computeValues(Tag *ge) {
            TClass2 *gramVals = gram->values(ge->getStrand()) + ge->getPos() + this->offset();
            for(register int i = 0; i < this->length(); i += this->step()) {
                int ni = i/this->step();
                int offset_v = ni*this->domainSize();
                int offset_i = ni*this->numPosVals();

                for(register int j = 0; j < this->numPosVals(); j++)
                    //          cerr << "\t" << i << " " << j << " " << vals[j] << " " << offset + vals[j] << endl;
                    this->vals[offset_i + j] = offset_v + gramVals[i+j];
            }
            return this->vals;
        }
        /*    inline TClass* computeValues(Tag *ge) {
              this->vals = gram->values(ge->getStrand()) + ge->getPos() - this->offset() - this->windowSize()/2;
              return this->vals + this->windowSize()/2;
              }*/


        ~FT_GramWWAM() { }

    };


    /**
     * FT_BinnedGramWWAM is a subtype of FT_GramWWAM that uses gram values at each signal
     * position.
     **************************************************************************/
    template <class TClass1, class TClass2, class TClass3>
        class FT_BinnedGramWWAM : public FT_BinnedWWAM<TClass1, TClass3> {
    protected:
        FL_BaseGram<TClass2, SPARSE_HASH<TClass3, int> > *gram;
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
         * @param signal the signal this feature is associated with.
         * @param windowSize the size of the window around any position in the
         * WWAM which is used for counting symbols or words; windowSize/size
         * must be an odd number.
         * @param gram a pointer to a FL_Gram object which defines the words that
         * will be counted at position.
         * @param offset the length of the WWAM, offset of the signal occurrence.
         * @param length  the length of the WWAM, length of the signal occurrence.
         * The total length of the WWAM would be then offset + length.
         * @param step the number of positions to jump within the WWAM, before
         * counting again.
         */
    FT_BinnedGramWWAM(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        FL_Signal<EdgeInst> *signal,
        TValType type,
        int windowSize,
        FL_FastGram<TClass2, SPARSE_HASH<TClass3, int> > *gram,
        int offset,
        int length,
        int step,
        FT_BaseBin<TClass1> *binObj
        )
        : FT_BinnedWWAM<TClass1, TClass3>(fInd,
                                          paramInd,
                                          parsingFrames,
                                          name,
                                          fe,
                                          signal,
                                          type,
                                          windowSize,
                                          gram->maxNumFilterValues(),
                                          offset,
                                          length,
                                          step,
                                          binObj) {

            this->gram = gram;
            initialize();

        }

        /**
         * Constructor from a Header string definition. windowSize should be read
         * from the Header definition.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * object to which this feature is tied to.
         * @param fargs The Header string definition loaded as a vector of strings.
         * The Header has the following form:\n\n
         * Feature name FT_BinnedGramWWAM parsingFrames signal offset length step windowSize
         * gram \n\n
         * The description of the fields could be found in
         * the other constructor(s)
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param fte A pointer to a FeatureEngine object.
         */
    FT_BinnedGramWWAM(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_BinnedWWAM<TClass1, TClass3>(fInd,
                                          paramInd,
                                          fargs, offset,
                                          fe,
                                          fte,
                                          fe->getFilter(fargs[offset + 9])->maxNumFilterValues()) {

            this->gram = (FL_FastGram<TClass2, SPARSE_HASH<TClass3, int> > *)fe->getFilter(fargs[offset++]);
            initialize();

        }

    FT_BinnedGramWWAM(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte,
        int windowSize
        )
        : FT_BinnedWWAM<TClass1, TClass3>(fInd,
                                          paramInd,
                                          fargs, offset,
                                          fe,
                                          fte,
                                          windowSize,
                                          fe->getFilter(fargs[offset + 8])->maxNumFilterValues()) {

            this->gram = (FL_FastGram<TClass2, SPARSE_HASH<TClass3, int> > *)fe->getFilter(fargs[offset++]);
            initialize();

        }

        void initialize() {
            this->_order = (int)gram->order();
            this->_numPosVals = this->windowSize() - this->order() + 1;
            this->_numUpdates = ((this->length() - 1)/this->step() + 1)*this->numPosVals();
            this->_numParams.f = ((this->length() - 1)/this->step() + 1)*this->domainSize();
            this->_numParams.s = this->binObj->maxNumFeatValues();
            //      this->_numParams.s = (_offset + _length + 1)*domSize;

            this->vals = new int [this->numUpdates()];
            this->counts = new SPARSE_HASH<TClass3, int> * [this->numUpdates()];

        }


        /**
         * A member function that computes protected member vals. Any entry
         * appearing at some position within array vals is assigned the
         * gram appearing at that position.
         * The positions in the sequence that are covered by vals are between
         * ge->getPos() + offset and ge->getPos() + offset + length
         * @return vals
         */
        inline int* computeValues(Tag *ge) {
            SPARSE_HASH<TClass3, int> *gramVals = gram->values(ge->getStrand()) +
                ge->getPos() + this->offset();

            for(register int i = 0; i < this->length(); i += this->step()) {
                int ni = i/this->step();
                int offset_v = ni*this->domainSize();
                int offset_i = ni*this->numPosVals();

                for(register int j = 0; j < this->numPosVals(); j++) {
                    //          cerr << "\t" << i << " " << j << " " << vals[j] << " " << offset + vals[j] << endl;
                    this->vals[offset_i + j] = offset_v;
                    this->counts[offset_i + j] = &gramVals[i+j];
                }
            }

            return this->vals;

        }

        ~FT_BinnedGramWWAM() {

        }

    };

}

#endif
