#include "GlobalVector.h"

/****************************************************************************
 * GlobalVector.cpp - part of the lless namespace, a general purpose
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

using namespace lless;

namespace lless {

    /**
     * Allocates memory for the parameter vector
     */
    void GlobalVector::allocate() {
        int i, j;
        unsigned int ind;
        this->featVals = new FeatureVector ** [features->size()];

        for(ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0)
                continue;

            Pair<ULONG> p = f->maxNumParams();

            this->featVals[f->paramInd()] = new FeatureVector *[p.f];
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                if(f->isSparse() || this->isSparse())
                    fv[i] = new SparseFeatureVector(p.s, defaultVal);
                else
                    fv[i] = new DenseFeatureVector(p.s);
        }
    }

    /**
     * Computes the average values of each feature. Stores information in
     * vector v.
     */
    void GlobalVector::computeAvgParamVals(list<pair<int, double> > &v) const {
        int i, j;
        int ind;
        double result = 0;

        for(ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];
            if(f->paramInd() < 0 || f->frozen())
                continue;

            Pair<ULONG> p = f->maxNumParams();
            const FeatureVector **fv =
                (const FeatureVector **)this->featVals[f->paramInd()];

            double avg = 0;
            for(int i = 0; i < p.f; i++)
                avg += fv[i]->averageValue();

            if(avg) avg /= p.f;

            list<pair<int, double> >::iterator it = v.begin();
            while(it != v.end() && avg < it->second)
                it++;
            v.insert(it, pair<int, double>(ind, avg));
        }
    }


    /**
     * Computes the maximal values of each feature. Stores information in
     * vector v.
     */
    void GlobalVector::computeMaxParamVals(list<pair<FeatVectorInd, double> > &v,
                                           int numMaxVals) const {

        double result = 0;

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];
            if(f->paramInd() < 0 || f->frozen())
                continue;

            Pair<ULONG> p = f->maxNumParams();
            const FeatureVector **fv =
                (const FeatureVector **)this->featVals[f->paramInd()];

            list<pair<FeatVectorInd, double> >::iterator it = v.begin();
            list<pair<FeatVectorInd, double> > maxV;

            for(int i = 0; i < p.f; i++)
                fv[i]->maxValues(maxV, ind, i, numMaxVals*p.f);

            list<pair<FeatVectorInd, double> >::iterator it2 = maxV.begin();
            for( ; it2 != maxV.end(); it2++) {
                //	  cout << f->getName() << " " << it2->first << " " << it2->second << endl;
                while(it != v.end() && it2->second < it->second)
                    it++;

                v.insert(it, *it2);
                //	  cout << " inserted " << " " << it->second << endl;
            }
        }
    }

    /**
     * Reinitializes the parameter values to zero
     */
    GlobalVector & GlobalVector::operator=(double val) {
        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0 || f->frozen())
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                (*fv[i]) = val;

        }

        return *this;
    }

    /**
     * Computes the dot product operation of *this and parameter vector gv
     * @return the result of the dot product.
     */
    double GlobalVector::operator*(const GlobalVector & gv) const {

        assert(features == gv.features && !defaultVal);
        double result = 0;

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];
            if(f->paramInd() < 0 || f->frozen())
                continue;

            Pair<ULONG> p = f->maxNumParams();
            const FeatureVector **fv =
                (const FeatureVector **)this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                result += (*fv[i])*(*gv[f->paramInd()][i]);

        }

        return result;
    }

    /**
     * Computes this[i] = this[i]*gv[i], for each entry i.
     * @return the resulting parameter vector (*this)
     */
    GlobalVector & GlobalVector::product(const GlobalVector & gv) {

        assert(features == gv.features && !defaultVal);

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0 || f->frozen())
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                fv[i]->product(*gv[f->paramInd()][i]);

        }

        return *this;
    }

    GlobalVector & GlobalVector::productInverse(const GlobalVector & gv) {

        assert(features == gv.features && !defaultVal);

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0 || f->frozen())
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                fv[i]->productInverse(*gv[f->paramInd()][i]);

        }

        return *this;
    }

    /**
     * Adds gv to *this.
     * @return the resulting parameter vector (*this)
     */
    GlobalVector & GlobalVector::operator+=(const GlobalVector & gv) {

        assert(features == gv.features);

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0 || f->frozen())
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                (*fv[i])+=(*gv[f->paramInd()][i]);

        }

        return *this;
    }

    /**
     * Set *this equal to gv.
     * @return *this
     */
    GlobalVector & GlobalVector::operator=(const GlobalVector & gv) {

        assert(features == gv.features);

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0)
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                (*fv[i])=(*gv[f->paramInd()][i]);

        }

        return *this;
    }

    /**
     * Substracts gv from *this
     * @return the resulting parameter vector (*this)
     */
    GlobalVector & GlobalVector::operator-=(const GlobalVector & gv) {

        assert(features == gv.features);

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0 || f->frozen())
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                (*fv[i])-=(*gv[f->paramInd()][i]);

        }

        return *this;

    }

    GlobalVector & GlobalVector::substractInverse(const GlobalVector & gv) {

        assert(features == gv.features);

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0 || f->frozen())
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                fv[i]->substractInverse(*gv[f->paramInd()][i]);

        }

        return *this;

    }

    /**
     * Checks whether *this is equal to gv.
     * @return true if *this equals gv, false otherwise
     */
    bool GlobalVector::operator==(const GlobalVector & gv) const {

        assert(features == gv.features && !defaultVal);

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0)
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++) {
                bool eq = ((*fv[i])==(*gv[f->paramInd()][i]));
                if(!eq)
                    return false;
            }
        }

        return true;
    }

    /**
     * Divides every element in *this by denom.
     * @return the resulting parameter vector (*this)
     */
    GlobalVector & GlobalVector::operator/(double denom) {

        assert(!defaultVal && denom != 0);

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0 || f->frozen())
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                (*fv[i])/denom;

        }

        return *this;
    }

    /**
     * Multiplies every element in *this by denom.
     * @return the resulting parameter vector (*this)
     */
    GlobalVector & GlobalVector::operator*(double factor) {

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0 || f->frozen())
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                (*fv[i])*factor;

        }
        return *this;
    }

    /**
     * Divides every element of gv by accumIterations; the result is
     * stored in *this.
     * @return the resulting parameter vector (*this)
     */
    GlobalVector & GlobalVector::average(int accumIterations,
                                         const GlobalVector & gv,
                                         bool forceFrozen) {

        assert(accumIterations && features == gv.features && !defaultVal);

        // Average state this->vals

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0 || (!forceFrozen && f->frozen()))
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                fv[i]->average(accumIterations, *gv[f->paramInd()][i]);

        }

        return *this;
    }

    /**
     * Read parameter vector from file input stream params
     */
    void GlobalVector::retrieve(::ifstream & params, int numFeatures) {
        int ind;
        std::string header;
        double val;
        bool isBinary = true;

        if(!numFeatures)
            numFeatures = features->size();

        params >> header; // string "model="

        if(!header.compare("model="))
            isBinary = false;

        if(isBinary && header.compare("binary"))
            throw EXCEPTION(BAD_USAGE,
                            header + string(" as header instead of \"model=\" or \"binary\""));

        for(ind = 0; ind < numFeatures; ind++) {
            Feature *f = (*features)[ind];
            if(f->paramInd() < 0)
                continue;

            std::string featName;
            int numfvs, i = 0;
            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            params >> featName >> numfvs;

            //      cerr << "Name" << " " << featName << " " << numfvs << endl;
            //      cerr << "Coords " << p.f << " " << p.s << endl;

            if(featName.compare(f->getName()))
                throw EXCEPTION(BAD_USAGE,
                                featName + string(" differs from ") + f->getName());

            if(params.eof())
                return;

            if(isBinary)
                params.get(); // new line

            while(i < numfvs) {
                //	cerr << i << endl;
                fv[i]->retrieve(params, isBinary);
                i++;
            }
        }

    }

    /**
     * Store parameter vector int file output stream params
     */
    void GlobalVector::store(::ofstream & params, bool isBinary) {
        if(isBinary)
            params << "binary" << endl;
        else
            params << "model=" << endl;

        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0)
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            params << f->getName() << "\t" << p.f << endl;

            for(int i = 0; i < p.f; i++)
                fv[i]->store(params, isBinary);

        }
    }

    /**
     * A member function that checks whether *this is an element of array
     * gvs of size size
     */
    bool GlobalVector::isMember(GlobalVector **gvs, int size) {
        for(int i = 0; i < size; i++) {
            if((*this) == (*gvs[i]))
                return true;
        }
        return false;
    }


    /**
     * Free memory used to store parameter vector
     */
    void GlobalVector::deallocate() {
        // Freeing this->vals
        for(UINT ind = 0; ind < features->size(); ind++) {
            Feature *f = (*features)[ind];

            if(f->paramInd() < 0)
                continue;

            Pair<ULONG> p = f->maxNumParams();
            FeatureVector **fv = this->featVals[f->paramInd()];

            for(int i = 0; i < p.f; i++)
                delete fv[i];

            delete [] fv;
        }

        delete [] this->featVals;
    }

} // end less namespace
