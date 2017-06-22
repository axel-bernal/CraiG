/****************************************************************************
 * FT_MultiBin.h - part of the lless namespace, a general purpose
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

#ifndef FT_MULTIBIN_FEATURE_H
#define FT_MULTIBIN_FEATURE_H

#include "FT_Bin.h"
#include "Filter.h"
#include "FilterEngine.h"
#include "FeatureEngine.h"

namespace  lless {


    /**
     * FT_MultiBin creates a set of 2 FT_Bin objects and
     * performs multilinear interpolation betweeen them.
     * /todo make it work for N bin Objects, needs to use recursive formula
     * for this
     *
     ***************************************************************************/

    template <class TClass> class FT_MultiBin : public TypedFeature<TClass> {
    private:
        double *featValues;

    protected:
        int numBinObjs;
        FT_BaseBin<TClass> **binObjs;
        Feature **featObjs;
        int *C;
        double *V;
        bool deleteObjs;

    public:
        /**
         * Default constructor
         * @param fInd A feature index.
         * @param paramInd feature's parameter vector index. If negative, the
         * feature does not have any ties to Tag objects.
         * @param name A unique feature identifier.
         * @param fe A pointer to a FilterEngine object.
         * @param type the feature value type, one of TValType enumerate type
         */
    FT_MultiBin(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        int numBinObjs,
        Feature **featObjs,
        FT_BaseBin<TClass> **binObjs,
        TValType type = FT_INTEGER
        )
        : TypedFeature<TClass>(fInd, paramInd,
                               parsingFrames,
                               name, fe,
                               1, type) {

            this->numBinObjs = numBinObjs;
            this->featObjs = featObjs;
            this->binObjs = binObjs;
            this->deleteObjs = false;
            this->_numParams.s = 1;
            featValues = new double [numBinObjs];
            C = new int[MAX_FEATVAL_UPDATES*numBinObjs];
            V = new double[MAX_FEATVAL_UPDATES*numBinObjs];

            for(int i = 0; i < numBinObjs; i++)
                this->_numParams.s *= binObjs[i]->maxNumFeatValues();
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
         * Feature name MultiLinearBin<TClass> integerValuedFilter numBins
         * binWidth \n\n
         * The Header above creates a FT_Bin feature whose contained feature is
         * of type TClass and bins of size binWidth. The feature can handle values
         * of up to rangeStart + numBinsxbinWidth, but any value above
         * numBinsxbinWidth + rangeStart is extrapolated, depending on the
         * model
         *
         * @param offset The index for vector fargs.
         * @param fe A pointer to a FilterEngine object.
         * @param fte A pointer to a FeatureEngine object.
         */

    FT_MultiBin(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte)
        : TypedFeature<TClass>(fInd, paramInd,
                               fargs, offset,
                               fe, 0) {

            if(!sscanf(fargs[offset++].c_str(), "%d", &numBinObjs))
                assert(0);

            binObjs = new FT_BaseBin<TClass> * [numBinObjs];
            featObjs = new Feature * [numBinObjs];
            this->deleteObjs = true;
            this->_numParams.s = 1;
            featValues = new double [numBinObjs];
            C = new int[MAX_FEATVAL_UPDATES*numBinObjs];
            V = new double[MAX_FEATVAL_UPDATES*numBinObjs];

            FT_BaseBin<TClass> *bin;

            for(int i = 0; i < numBinObjs; i++) {
                featObjs[i] = fte->getFeature(fargs[offset++]);
                bin = (FT_BaseBin<TClass> *)fte->getFeature(fargs[offset++]);
                binObjs[i] = bin;
                this->_numParams.s *= bin->maxNumFeatValues();
            }
        }

        inline bool isFeatureSet() {
            return true;
        }

        inline double _dotParamV(double *featValues, const FeatureVector *params,
                                 int numFeatVals, int N, int fConjInd = 0) {

            if(!N) throw EXCEPTION(NOT_SUPPORTED, "At least one bin object must be given");

            FT_BaseBin<TClass> *bin = binObjs[N - 1];
            double featValue = featValues[N - 1];
            double result = 0;
            int i;
            int *C = this->C + MAX_FEATVAL_UPDATES*(N-1);
            double *V = this->V + MAX_FEATVAL_UPDATES*(N-1);

            int numCPUpdates = bin->featValCs(featValue, C);
            bin->featValVs(featValue, C, V, numCPUpdates);

            if(N == 1) {
                for(i = 0; i < numCPUpdates; i++) {
                    assert(fConjInd + C[i] >= 0 && fConjInd + C[i] < this->maxNumFeatValues());
                    result += V[i]*(*params)[fConjInd + C[i]];
                }
                return result;
            }

            numFeatVals = numFeatVals/bin->maxNumFeatValues();

            for(i = 0; i < numCPUpdates; i++) {
                double r = _dotParamV(featValues, params, numFeatVals,
                                      N - 1, fConjInd + C[i]*numFeatVals);
                result += V[i]*r;
            }

            return result;

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
         * parameter which corresponds to the bin that is 'on'.
         * @see double TypedFeature<int>::dotParamV(double, double **, int)
         */

        inline double dotParamV(Tag *ge, int frame,
                                const FeatureVector **params,
                                int fConjInd = 0,
                                TypedFilter<UCHAR> *filter = NULL) {

            assert(!filter);

            for(int i = 0; i < numBinObjs; i++)
                featValues[i] = featObjs[i]->featValue(ge, frame);

            return _dotParamV(featValues, params[fConjInd],
                              this->_numParams.s, numBinObjs, 0);

        }

        inline double dotParamV(double *featValues,
                                const FeatureVector **params,
                                int fConjInd = 0) {

            return _dotParamV(featValues, params[fConjInd],
                              this->_numParams.s, numBinObjs, 0);

        }

        inline void _updParamV(double updVal, double *featValues,
                               FeatureVector *params,
                               int numFeatVals, int N, int fConjInd = 0) {

            if(!N) throw EXCEPTION(NOT_SUPPORTED, "At least one bin object must be given");

            FT_BaseBin<TClass> *bin = binObjs[N - 1];
            double featValue = featValues[N - 1];
            int i;
            int *C = this->C + MAX_FEATVAL_UPDATES*(N-1);
            double *V = this->V + MAX_FEATVAL_UPDATES*(N-1);

            int numCPUpdates = bin->featValCs(featValue, C);
            bin->featValVs(featValue, C, V, numCPUpdates);

            if(N == 1) {
                for(i = 0; i < numCPUpdates; i++) {
                    //	  cerr << "updating " << V[i] << " " << updVal << " to " << fConjInd << "(" << C[i] << ")" <<endl;
                    (*params)[fConjInd + C[i]] += V[i]*updVal;
                }
                return;
            }

            numFeatVals = numFeatVals/bin->maxNumFeatValues();

            for(i = 0; i < numCPUpdates; i++) {
                //	cerr << "level " << N << " " << featValue << " " << C[i] << " " << V[i] << endl;
                _updParamV(V[i]*updVal, featValues, params, numFeatVals,
                           N - 1, fConjInd + C[i]*numFeatVals);
            }
        }

        /**
         * A member function to update the feature's parameter which corresponds
         * to the bin that is 'on'.
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
        inline void updParamV(double updVal, Tag *ge, int frame,
                              FeatureVector **params, int fConjInd = 0,
                              TypedFilter<UCHAR> *filter = NULL) {

            assert(!filter);
            //      cerr << this->getName() << " " << ge->getPos() << " " << ge->getPos() +
            //	ge->getLen() - 1 << " " << ge->getParseType() << " " << frame << " " << updVal << endl;
            for(int i = 0; i < numBinObjs; i++) {
                featValues[i] = featObjs[i]->featValue(ge, frame);
                //	cerr << i << " " << featValues[i] << endl;
            }

            _updParamV(updVal, featValues, params[fConjInd],
                       this->_numParams.s, numBinObjs, 0);

        }

        inline void updParamV(double updVal, double *featValues,
                              FeatureVector **params, int fConjInd = 0) {

            return _updParamV(updVal, featValues, params[fConjInd],
                              this->_numParams.s, numBinObjs, 0);

        }


        ~FT_MultiBin() {
            delete [] featValues;
            delete [] C;
            delete [] V;
            if(this->deleteObjs) {
                delete [] binObjs;
                delete [] featObjs;
            }
        }

    };

}

#endif
