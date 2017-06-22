/****************************************************************************
 * FT_CSplineBin.h - part of the lless namespace, a general purpose
 *            linear semi-markov structure prediction library
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

#ifndef FT_CSPLINEBIN_FEATURE_H
#define FT_CSPLINEBIN_FEATURE_H

#include "FT_Bin.h"

namespace  lless {

    /**
     * FT_CSplineBin derives FT_Bin. It creates a set of equidistant ordered
     * disjoint control points which cover the range of values of the contained
     * feature and uses them to fit a c-spline function to the data between each
     * of them. The contained feature's value is used to find out between which
     * control points the value falls in.
     * This class is useful to model contained features with sparse or
     * multimodal distributions, such as lengths and scores.
     *
     ***************************************************************************/

    template <class TClass> class FT_CSplineBin : public FT_Bin<TClass> {

    public:
        /**
         * Default constructor
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param name A unique feature identifier.
         * @param fe A pointer to a FilterEngine object.
         * @param ft A pointer to a TypedFeature<TClass>, the contained feature.
         * @param maxfeatVals The number of bins that are created.
         * @param _rangeStart Where binning starts
         * @param _binWidth The width of each bin
         * values become more sparse, the larger they get.
         * @param type the feature value type, one of TValType enumerate type
         */
    FT_CSplineBin(
        int fInd,
        int paramInd,
        char *name,
        FilterEngine *fe,
        TypedFeature<TClass> *ft,
        int maxNumFeatVals,
        double rangeStart,
        double binWidth,
        TValType type = FT_INTEGER
        )
        : FT_Bin<TClass>(fInd, paramInd,
                         name, fe, ft,
                         maxNumFeatVals + 4,
                         rangeStart - 2*binWidth,
                         binWidth, type) {


        }

        /**
         * Constructor from a Header string definition. The feature value type
         * is provided as parameter.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * object to which this feature is tied to.
         * @param fargs The Header string definition loaded as a vector of strings.
         * The Header has the following form:\n\n
         * Feature name Bin<TClass> integerValuedFilter numBins _rangeStart
         * _binWidth monotone\n\n
         * The Header above creates a FT_Bin feature whose contained feature is
         * of type TClass and bins of size _binWidth. The feature can handle values
         * of up to _rangeStart + numBinsx_binWidth, but any value above
         * numBinsx_binWidth + _rangeStart is extrapolated, depending on the
         * model
         *
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param fte A pointer to a FeatureEngine object.
         */

    FT_CSplineBin(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_Bin<TClass>(fInd, paramInd,
                         fargs, offset,
                         fe, fte) {

            this->_rangeStart -= 2*this->_binWidth;
            /*
             * Modifying the number of parameters to account for boundary
             * control point interpolation
             */
            this->_numParams.s += 4;
        }

        int lastFloor() {
            return this->maxNumFeatValues() - 2;
        }

        int featValCs(double featVal, int *C) {
            C[1] = this->floor(featVal);

            if(this->outOfRange) {
                C[0] = C[1];
                return 1;
            }

            C[2] = this->ceiling(featVal, C[1]);

            assert(C[1] != C[2]);

            if(!C[1]  || C[2] == this->maxNumFeatValues() - 1) {
                C[0] = C[1]; C[1] = C[2];
                return 2;
            }

            C[0] = this->floor(featVal - this->binWidth());
            C[3] = this->ceiling(featVal + this->binWidth(), C[2]);

            return 4;

        }

        void featValVs(double featVal, const int *C,
                       double *V, int numCPUpdates) {

            if(this->outOfRange) {
                assert(numCPUpdates == 1);
                V[0] = 1.0;
                this->outOfRange = false;
                return;
            }

            if(numCPUpdates == 2) {
                V[0] = (this->cpValue(C[1]) - featVal)/this->binWidth();
                V[1] = (featVal - this->cpValue(C[0]))/this->binWidth();
                return;
            }

            assert(numCPUpdates == 4);

            double t = (featVal - this->cpValue(C[1]))/this->binWidth();
            double t2 = t*t;
            double t3 = t2*t;

            V[0] = (-t3 + 2*t2 - t)/2;
            V[1] = (3*t3 - 5*t2 + 2)/2;
            V[2] = (-3*t3 + 4*t2 + t)/2;
            V[3] = (t3 - t2)/2;

        }

        ~FT_CSplineBin() {

        }

    };


    /**
     * FT_CSplineCustomBin inherits from FT_CustomBin and creates a set of
     * arbitrarily spaced and ordered disjoint control points which cover the
     * range of values of the contained feature  and uses them to fit a c-spline
     * function of the data between each of them. The contained feature's value
     * is used to find out between which. The contained feature's value is used
     * to find out between which control points the value falls in.
     * This class is useful to model contained features with sparse or
     * multimodal distributions, such as lengths and scores.
     * \todo make the interpolation work for non-equidistant points
     *
     ***************************************************************************/

    template <class TClass> class FT_CSplineCustomBin : public FT_CustomBin<TClass> {

    public:
        /**
         * Default constructor
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param name A unique feature identifier.
         * @param fe A pointer to a FilterEngine object.
         * @param ft A pointer to a TypedFeature<TClass>, the contained feature.
         * @param cps the custom-made control points that define the bins
         * @param type the feature value type, one of TValType enumerate type
         */
    FT_CSplineCustomBin(
        int fInd,
        int paramInd,
        char *name,
        FilterEngine *fe,
        TypedFeature<TClass> *ft,
        vector<double> &cps,
        TValType type = FT_INTEGER
        )
        : FT_CustomBin<TClass>(fInd, paramInd,
                               name, fe, ft,
                               cps, type) {

            this->_numParams.s += 5;

            this->cps.insert(this->cps.begin(), 2*this->cps[0] - this->cps[1]);
            this->cps.insert(this->cps.begin(), 2*this->cps[0] - this->cps[1]);
            this->cps.push_back(2*this->cps[this->cps.size() - 1] -
                                this->cps[this->cps.size() - 2]);
            this->cps.push_back(2*this->cps[this->cps.size() - 1] -
                                this->cps[this->cps.size() - 2]);

        }

        /**
         * Constructor from a Header string definition. The feature value type
         * is provided as parameter.
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * object to which this feature is tied to.
         * @param fargs The Header string definition loaded as a vector of strings.
         * The Header has the following form:\n\n
         * Feature name CustomRBin<TClass> integerValuedFilter numBins [bins ..]
         *
         * The Header above creates a FT_CustomRBin feature whose contained feature
         * is of type TClass and numBins bins.
         *
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param fte A pointer to a FeatureEngine object.
         */

    FT_CSplineCustomBin(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_CustomBin<TClass>(fInd, paramInd,
                               fargs, offset,
                               fe, fte) {

            /*
             * Modifying the number of parameters to reflect the actual number of
             * control points (and not bins)
             */
            this->_numParams.s += 5;

            this->cps.insert(this->cps.begin(), 2*this->cps[0] - this->cps[1]);
            this->cps.insert(this->cps.begin(), 2*this->cps[0] - this->cps[1]);
            this->cps.push_back(2*this->cps[this->cps.size() - 1] -
                                this->cps[this->cps.size() - 2]);
            this->cps.push_back(2*this->cps[this->cps.size() - 1] -
                                this->cps[this->cps.size() - 2]);

        }

        int lastFloor() {
            return this->maxNumFeatValues() - 2;
        }

        int featValCs(double featVal, int *C) {
            C[1] = this->floor(featVal);
            C[2] = this->ceiling(featVal, C[1]);

            assert(C[1] != C[2]);

            if(!C[1]  || C[2] == this->maxNumFeatValues() - 1) {
                C[0] = C[1]; C[1] = C[2];
                return 2;
            }

            C[0] = (C[1] < 1 ? C[1] : C[1] - 1);
            C[3] = (C[2] >= this->maxNumFeatValues() ? C[2] : C[2] + 1);

            return 4;

        }

        void featValVs(double featVal, const int *C,
                       double *V, int numCPUpdates) {

            if(numCPUpdates == 2) {
                V[0] = (this->cpValue(C[1]) - featVal)/this->binWidth(C[0]);
                V[1] = (featVal - this->cpValue(C[0]))/this->binWidth(C[0]);
                return;
            }

            assert(numCPUpdates == 4);

            double t = (featVal - this->cpValue(C[1]))/this->binWidth(C[0]);
            double t2 = t*t;
            double t3 = t2*t;

            V[0] = (-t3 + 2*t2 - t)/2;
            V[1] = (3*t3 - 5*t2 + 2)/2;
            V[2] = (-3*t3 + 4*t2 + t)/2;
            V[3] = (t3 - t2)/2;

        }

        ~FT_CSplineCustomBin() {

        }

    };
}

#endif
