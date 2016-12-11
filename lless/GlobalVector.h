/****************************************************************************
* GlobalVector.h - part of the lless namespace, a general purpose
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

#ifndef _GLOBAL_VECTOR_H_
#define _GLOBAL_VECTOR_H_
#include <vector>
#include "Feature.h"
#include "FeatureVector.h"

namespace lless {

  /**
   * The GlobalVector class implements a global feature vector associated with
   * features passed as arguments in the constructor
   ***************************************************************************/
  
  class GlobalVector {
   private:
    vector<Feature *> *features;
    FeatureVector ***featVals;
    bool sparse;
    double defaultVal;
    
   public: 
    //! Default contructor
    GlobalVector(bool sparse = false, double defaultVal = 0.0) {
      this->sparse = sparse;
      this->defaultVal = defaultVal;
    }
    
    /**
     * Main constructor
     * @param features the list of features. Features which are tied to Tag
     * elements need to access the parameter vector to retreive their own
     * parameters.
     * @param fd the file containing the parameter values.
     */
    GlobalVector(vector<Feature *> & features, bool sparse = false, 
		 double defaultVal = 0.0) {
      this->sparse = sparse;
      this->defaultVal = defaultVal;
      setFeatures(features);
    }

    GlobalVector(vector<Feature *> & features, ::ifstream *fd, 
		 bool sparse = false, double defaultVal = 0.0) {

      this->sparse = sparse;
      this->defaultVal = defaultVal;
      setFeatures(features);
      retrieve(*fd);

    }		 

    GlobalVector(vector<Feature *> & features, ::ifstream *fd, 
		 int addedFeatures, bool sparse = false, double defaultVal = 0.0) {

      this->sparse = sparse;
      this->defaultVal = defaultVal;
      setFeatures(features);
      retrieve(*fd, features.size() - addedFeatures);

    }		 
    
    bool isSparse() {
      return sparse;
    }

    inline FeatureVector ***values() {
      return featVals;
    }

    inline void setFeatures(vector<Feature *> &features) {
      this->features = &features;
      allocate();
    }

    FeatureVector **operator[](int pInd) {
      return featVals[pInd];
    }

    const FeatureVector **operator[](int pInd) const {
      return (const FeatureVector **)featVals[pInd];
    }

    inline GlobalVector & operator*=(double f) {
      return (*this)*f;
    }
    
    inline GlobalVector & operator/=(double f) {
      return (*this)/f;
    }

    inline void print(std::ostream &ost) {
      //      for(unsigned int ind = 0; ind < features->size(); ind++) {
      for(unsigned int ind = 107; ind <= 107; ind++) {
	Feature *f = (*features)[ind];
	f->printParams(ost, (const FeatureVector ***)featVals);
      }
    }

    inline GlobalVector & highFilter(const GlobalVector &gv, double v) {
      for(unsigned int ind = 0; ind < features->size(); ind++) {
	Feature *f = (*features)[ind];
	if(f->paramInd() < 0 || f->frozen())
	  continue;
	
	Pair<ULONG> p = f->maxNumParams();
	FeatureVector **fv = this->featVals[f->paramInd()];

	for(int i = 0; i < p.f; i++)
	  fv[i]->highFilter(*gv[f->paramInd()][i], v);
	
      }
      return *this;
    }

    void allocate();
    double operator*(const GlobalVector& gv) const;
    GlobalVector & operator=(double);
    GlobalVector & product(const GlobalVector & gv);
    GlobalVector & operator+=(const GlobalVector & gv);
    GlobalVector & operator=(const GlobalVector & gv);
    GlobalVector & operator-=(const GlobalVector & gv);
    GlobalVector & productInverse(const GlobalVector & gv);
    GlobalVector & addInverse(const GlobalVector & gv);
    GlobalVector & assignInverse(const GlobalVector & gv);
    GlobalVector & substractInverse(const GlobalVector & gv);
    
    void computeAvgParamVals(list<pair<int, double> > &v) const;    
    void computeMaxParamVals(list<pair<FeatVectorInd, double> > &,
			     int maxNumVals) const;
    
    bool operator==(const GlobalVector & gv) const;
    GlobalVector & operator/(double);
    GlobalVector & operator*(double);
    GlobalVector & average(int accumIterations, const GlobalVector & gv, bool = false);

    void retrieve(::ifstream & params, int numFeatures = 0);
    void store(::ofstream &params, bool isBinary = false);
    bool isMember(GlobalVector **gvs, int size);
    void deallocate();

    ~GlobalVector() {
      deallocate();
    }    

  };

}

#endif
