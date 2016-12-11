/****************************************************************************
* Feature.h - part of the lless namespace, a general purpose
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

#ifndef _FEATURE_H_
#define _FEATURE_H_

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <set>
#include "Filter.h"
#include "FL_Signal.h"
#include "Utils.h"
#include "FilterEngine.h"
#include "FeatureVector.h"
#include <sstream>

#define INVALID_FEATVAL exp((double)DBL_MAX_EXP - 1)
#define NUM_FEAT_TYPES 3

namespace lless {

  //! enumerate for types of features, depending on the Tag object they tie.
  typedef enum {ST, TR, WO} TFeatureTypeEnum;
 
  class FeatureEngine;

  /**
   * A Feature is a property defined at each position of the input sequence and
   * for all valid taggings at that position. This means that a Feature is 
   * usually tied to a subset of all possible Tags. Each Feature object knows
   * about its associated parameter(s) and as such, it can compute the 
   * "dot product" between its value and its parameter(s), or simply update 
   * its parameter(s).
   *
   * The "dot product" can be precomputed for certain features which have a
   * small number of valid Tags tied at each position which are known before
   * decoding. In thoses cases, there could be a significant speed up of the
   *
   ***************************************************************************/

  class Feature {
   protected:
    int fInd;                //!< feature index
    std::string name;        //!< name
    int pInd;                //!< parameter index -1, if no params
    FeatureVector ***params; // feature parameters after phase is resolved
    const FeatureVector ***cparams; // const version of feature parameters
    int _parsingFrames;      //!< number of frames of this feature at decoding time
    int *_preCompArrayIndex; //!< the array index per frame. != -1 if this 
                             //!< feature's precomputed its values

    double noiseFactor;
    TValType type;
    FilterEngine *fe;
    bool _collapseFrames;    /*!< For precomputable features, if false
                               all frames of this feature are
                               pre-computed in separated array entries.
                             */                               
    bool _phaseDependent;    /*!< true if a set of parameters is needed
                               for each frame. */
    bool _sparse;            //!< uses a sparse parameter vector if true     
    Pair<ULONG> _numParams;    //!< number of parameters (in two dimensions)

    double mean_fval;        //! mean of feature values
    double stdev_fval;       //!standard deviation of feature values
      
    bool _frozen;
    bool _off;
    int subfVal;             /*!< different from -1 when this is a Feature
                              * which we only want to update on subfVal 
                              */
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
     */
    Feature(
            int fInd, 
            int paramInd, 
            char *name, 
            FilterEngine *fe, 
            ULONG maxNumFeatVals = 0
            ) { 

      this->name = std::string(name);
      initFields(fInd, paramInd, fe, maxNumFeatVals, false, false);
    }
    
    /**
     * Constructor from a Header string definition.
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
    Feature(
            int fInd, 
            int paramInd, 
            std::vector<std::string> &fargs, 
            int & offset, 
            FilterEngine *fe, 
            ULONG maxNumFeatVals = 0
            ) { 

      this->name = fargs[offset++];
      initFields(fInd, paramInd, fe, maxNumFeatVals, false, false);
    }

    /**
     * A member function that initiliazes the Feature object's members
     */
    void initFields(int fInd, int paramInd,
                    FilterEngine *fe, ULONG maxNumFeatVals, 
                    bool frozen, bool off) { 

      this->fInd = fInd;
      this->pInd = paramInd;
      this->fe = fe;
      this->_frozen = frozen;
      this->_off = off;
      this->noiseFactor = 0;
      this->_phaseDependent = false;
      this->_sparse = false;
      setMaxNumParams(1, maxNumFeatVals);
      turnAllSubFeatsOn();
      params = NULL;
      cparams = (const FeatureVector ***)params;
      mean_fval = 0;
      stdev_fval = 1;
    }

    /**
     * A member function that the initialized the arrays for precomputation. 
     * By default all precomputation array indexes are initialized to -1.
     * @param parsingFrames the number of parsing frames (or phases) of this
     * feature
     */
    inline void initPreCompArrays(int parsingFrames) {

      this->_parsingFrames = parsingFrames;
      _preCompArrayIndex = new int [parsingFrames];

      for(int fr = 0; fr < _parsingFrames; fr++)
        _preCompArrayIndex[fr] = -1;
    }

    /**
     * @param frame the frame(or phase) which maps the 
     * precomputation array index.
     * @return the precomputation index 
     */
    inline int preCompArrayIndex(int frame) {
      return _preCompArrayIndex[frame];
    }

    inline TValType valType() {
      return type;
    }

    inline void turnSubFeatOn(int subfVal) {
      this->subfVal = subfVal;
    }

    inline void turnAllSubFeatsOn() {
      this->subfVal = -1;
    }

    inline bool subFeatIsOn(int subfVal) {
      return (this->subfVal == -1 || this->subfVal == subfVal);
              
    }

    inline std::string & getName() {
      return name;
    }
    
    inline int ind() {
      return fInd;
    }

    inline int paramInd() {
      return pInd;
    }

    inline ULONG maxNumFeatValues() {
      return this->_numParams.s;
    }

    inline Pair<ULONG> & maxNumParams() {
      return this->_numParams;
    }

    inline void setMaxNumParams(Pair<ULONG> &p) {
      this->_numParams.f = p.f;
      this->_numParams.s = p.s;
    }

    inline void setMaxNumParams(ULONG f, ULONG s) {
      this->_numParams.f = f;
      this->_numParams.s = s;
    }

    inline bool isOff() {
      return _off;
    }

    virtual inline bool isFeatureBag() {
      return false;
    }

    /**
     * @param ge Tag object to which this feature is tied.
     * @param frame Tag object's current phase(or frame).
     * @return the feature value, which depends tightly on the Tag object
     * to which the feature ties and its phase(frame), both of them given
     * as parameters.
     */
    virtual inline double featValue(Tag *ge, int frame) {
      assert(!isFeatureSet());
      return 0;
    }

    virtual inline FL_Signal<EdgeInst> * signal() { return NULL;}

    /**
     * A member function that freezes any update to the feature's 
     * corresponding parameter(s)
     */
    inline void freeze() {
      _frozen = true;
    }

    inline bool frozen() {
      return _frozen;
    }

    inline void unfreeze() {
      _frozen = false;
    }

    inline int parsingFrames() {
      return _parsingFrames;
    }  

    inline void setParams(FeatureVector ***params) {
      this->params = params;
      this->cparams = (const FeatureVector ***)params;
    }

    inline void setParamInd(int paramInd) {
      this->pInd = paramInd;
    }

    inline void makePhaseDependent() {
      this->_phaseDependent = true;
      _numParams.f = this->parsingFrames()*_numParams.f;

    }

    inline void makeSparse() {
      this->_sparse = true;
    }

    inline bool isSparse() {      
      return this->_sparse;
    }

    inline bool isPhaseDependent() {
      return this->_phaseDependent;
    }

    virtual inline bool isFeatureSet() {
      return false;
    }

    inline bool operator<(const Feature &f) {
      return ind() < ((Feature &)f).ind();
    }

    inline bool operator==(const Feature &f) {
      return ind() == ((Feature &)f).ind();
    }

    virtual inline bool isPreComputable() {
      return _preCompArrayIndex[0] >= 0;
    }

    /**
     * A member function that makes the feature return zero as the value
     * of the dot product between feature value and its corresponding
     *  parameter(s)
     */
    virtual inline void turnOff() {
      _off = true;
    }

    virtual inline void turnOn() {
      _off = false;
    }

    /**
     * @return the index of the last entry used in the array for precomputation
     */
    virtual inline int preComputeTo(bool collapseFrames, int initArrayIndex = 0) {
      int lastUsedArrayIndex = initArrayIndex;
      _collapseFrames = collapseFrames;
      //      cerr << "pcIndexes for " << this->getName() << endl;

      for(int frame = 0; frame < parsingFrames(); frame++) {
        if(!collapseFrames)
          lastUsedArrayIndex = initArrayIndex + frame;

        _preCompArrayIndex[frame] = lastUsedArrayIndex;
        //        cerr << "frame " << frame << " " << _preCompArrayIndex[frame];

      }

      return lastUsedArrayIndex;
    }

    virtual inline int begPreComp(int beg, int frame) {
      return beg;
    }

    virtual inline int endPreComp(int end, int frame) {
      return end;
    }

    virtual inline int stepPreComp() {
      return 1;
    }
    
    /**
     * A member function that performs dot product for feature precomputation.
     * The precomputation occurs at position pos and frame.
     */
    virtual inline double dotParamVPreComp(int pos,
					   int frame,
					   TStrand strand,
					   const FeatureVector **params,
					   int fConjInd = 0,
					   TypedFilter<UCHAR> *filter = NULL) {

      /* 
       * This is the best we can do with this information. Derived feature
       * classes should re-implement this function accodingly
       */
      Tag tag(INVALID_EDGE, pos, strand);
      int offset = fConjInd;

      if(filter && !this->isFeatureSet()) {
	offset += filter->value(pos, strand);
	filter = NULL;
      }

      return this->dotParamV(&tag, frame, params, offset, filter);

    }
    
    /**
     * A member function that performs dot product precomputation.
     * @param array the precomputation array
     * @param preCompArrayIndex index in the precomputation
     * array in which the precomputed dot product values will be stored.
     * @param beg position in the sequence at which precomputation starts
     * @param end position in the sequence at which precomputation stops
     * @param strand strand to precompute.
     * @param params parameter vector to be used for computing the dot 
     * product.
     * @param filter a pointer to a filter whose values along the sequence
     * will be used to choose the right entry in params.
     * @param fConjInd is an offset for the feature vector params, which is
     * different from zero when this feature has been used in a Feature 
     * conjunction or Feature disjunction operation (defined as the class
     * FeaturexFeature, in which case fConjInd usually denotes the value of 
     * the first feature; in the latter case, the values of the first feature 
     * must be countable and not too large.
     */
    virtual inline void doPreComputation(double **array,
                                         int preCompArrayIndex, 
                                         int beg, int end, TStrand strand, 
                                         const FeatureVector **params, 
                                         int frame,
                                         int fConjInd = 0,
                                         TypedFilter<UCHAR> *filter = NULL) {

      if(!this->isPreComputable()) {
        assert(0);
        throw EXCEPTION(BAD_USAGE, this->getName() + " cannot be precomputed.\n");
      }

      int arrayInd = preCompArrayIndex + _preCompArrayIndex[frame];
      int pos = begPreComp(beg, frame);
      //      cerr << this->getName() << " " <<  strand << " " << frame << " ";
      //      double fres = 0;
      for( ; pos <= endPreComp(end, frame); pos += stepPreComp()) {
	double res =  dotParamVPreComp(pos, frame,
				       strand, params, 
				       fConjInd, filter);

        array[arrayInd][pos] += res;
	//        fres += res;
      }
      //      cerr << fres << endl;
    }


    /**
     * @see void Feature::doPreComputation(double **, int, int, int, TStrand,
     * double **, TypedFilter<UCHAR> *, int = 0)
     */
    virtual inline void doPreComputation(double **array,
                                         int preCompArrayIndex, 
                                         int beg, int end, TStrand strand) {
      
      for(int frame = 0; frame < this->parsingFrames(); frame++) {
        int fConjInd = 0;

        if(this->isPhaseDependent())
          fConjInd = frame;

        doPreComputation(array, preCompArrayIndex, beg, end, 
                         strand, this->cparams[paramInd()], 
                         frame, fConjInd);
                         
      }

    }
    
    /**
     * Routine for caching  precomputed values of features in special 
     * feature-dependent structures
     */
    virtual inline void updCacheParams(const FeatureVector ***params, int fConjInd = 0) {
      
      return;
      
    }

    /**
     * A member function to accumulate precomputed dot products. Features
     * that need to compute a 'global' dot product value along a segment which
     * can be computed by summing over 'local' dot product values at each
     * segment position need to implement this function.
     * @param array the precomputation array
     * @param preCompArrayIndex index in the precomputation
     * array in which the precomputed dot product values will be stored.
     * @param period period of the precomputation array.
     * @param beg position in the sequence at which precomputation starts.
     * @param end position in the sequence at which precomputation stops.
     */
    virtual inline void accumulatePreCompEntries(double **array, 
                                                 int preCompArrayIndex, 
                                                 int period, int beg, 
                                                 int end) {

      ;
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
                                               Tag *ge, int frame) {
      
      return 0;
    }

    /**
     * @see double Feature::dotParamV(Tag *, int, double **, int)
     */
    virtual inline double dotParamV(Tag *ge, int frame) {
      if(isOff())
        return 0;

      if(this->isPhaseDependent())
        return dotParamV(ge, frame, this->cparams[paramInd()], frame); 
       
      return dotParamV(ge, frame, this->cparams[paramInd()]);

    }

    /**
     * @see double Feature::dotParamV(double, double **, int)
     */
    inline double dotParamV(double featVal, int fConjInd = 0) {
      return dotParamV(featVal, this->cparams[paramInd()], fConjInd);
    }


    /**
     * @param ge Tag object to which this feature is tied.
     * @param frame Tag object's current phase(or frame).
     * @param fConjInd is an offset for the feature vector params, which is
     * different from zero when this feature has been used in a Feature 
     * conjunction or Feature disjunction operation (defined as the class
     * FeaturexFeature, in which case fConjInd usually denotes the value of 
     * the first feature; in the latter case, the values of the first feature 
     * must countable.
     * @see double Feature::dotParamV(double, double **, int)
     */

    virtual inline double dotParamV(Tag *ge, int frame, 
                                    const FeatureVector **params,
				    int fConjInd = 0,
				    TypedFilter<UCHAR> *filter = NULL) {

      assert(!filter);

      double featVal = featValue(ge, frame);
      return dotParamV(featVal, params, fConjInd);

    }

    /**
     * Default dotParamV function implementation. The behaviour of this 
     * function works only for integer-valued features.
     * Features with other value types must overload this function.
     * @param featVal the feature value, computed with 
     * featValue(Tag *, int). It is upcasted to double at this
     * point since the result of the dot product is also double.
     * @param params parameter vector to be used for computing the dot 
     * product.
     * @param fConjInd is an offset for the feature vector params, which is
     * different from zero when this feature has been used in a Feature 
     * conjunction or Feature disjunction operation (defined as the class
     * FeaturexFeature, in which case fConjInd usually denotes the value of 
     * the first feature; in the latter case, the values of the first feature 
     * must countable.
     * @return The dot product between feature value and the feature's 
     * corresponding parameter(s).
     * @see double Feature::dotParamV(Tag *, int, double **, int, TypedFilter<UCHAR> *)
     */
    virtual inline double dotParamV(double featVal, 
				    const FeatureVector **params,
                                    int fConjInd = 0) {

      return (*params[fConjInd])[(int)featVal];

    }

    /**
     * For feature sets
     * @see double Feature::dotParamV(Tag *, int, double **, int, TypedFilter<UCHAR> *)
     */
    virtual inline double dotParamV(double *featValues, 
				    const FeatureVector **params,
                                    int fConjInd = 0) {

      return (*params[fConjInd])[(int)featValues[0]];

    }


    /**
     * @see void Feature::updParamV(double updVal, Tag *, int, vector<FeatureVector *>, int)
     */
    virtual inline void updParamV(double updVal, Tag *ge, int frame) {
      if(this->frozen())
        return;

      if(noiseFactor != 0)
        updVal *= noiseFactor;

      if(this->isPhaseDependent())
        updParamV(updVal, ge, frame, this->params[paramInd()], frame); 
      else
	updParamV(updVal, ge, frame, this->params[paramInd()]);
      
    }


    /**
     * @see void Feature::updParamV(double, double, double **, int)
     */
    inline void updParamV(double updVal, double featVal, int fConjInd = 0) {
      updParamV(updVal, featVal, this->params[paramInd()], fConjInd);
    }

    /**
     * A member function that updates the feature's corresponding parameter(s).
     * @param updVal the value of the update to be made to this feature's
     * corresponding parameter(s)
     * @param ge Tag object to which this feature is tied.
     * @param frame Tag object's current phase(or frame).
     * @param params parameter vector to be used for computing the dot 
     * product.
     * @param fConjInd is an offset for the feature vector params, which is
     * different from zero when this feature has been used in a Feature 
     * conjunction or Feature disjunction operation (defined as the class
     * FeaturexFeature, in which case fConjInd usually denotes the value of 
     * the first feature; in the latter case, the values of the first feature 
     * must countable.
     * @see void Feature::updParamV(double, double, double **, int)    
     */
    virtual inline void updParamV(double updVal, Tag *ge, int frame, 
				  FeatureVector **params,
				  int fConjInd = 0,
				  TypedFilter<UCHAR> *filter = NULL) {

      assert(!filter);

      double result = featValue(ge, frame);
      //      cerr << this->getName() << " " << ge->getGEClass() << " " << ge->getPos() << " " << result << " " << updVal << endl;
      updParamV(updVal, result, params, fConjInd);
    }


    /**
     * Default updParamV function implementation, the function that updates 
     * the feature's parameter(s). The behaviour of this function works only
     * for integer-valued features.
     * Features with other value types must overload this function.
     * @param updVal the value of the update to be made to this feature's
     * corresponding parameter(s)
     * @param featVal the feature value, computed with 
     * featValue(Tag *, int). It is upcasted to double at this
     * point since the result of the dot product is also double.
     * @param params parameter vector to be used for computing the dot 
     * product.
     * @param fConjInd is an offset for the feature vector params, which is
     * different from zero when this feature has been used in a Feature 
     * conjunction or Feature disjunction operation (defined as the class
     * FeaturexFeature, in which case fConjInd usually denotes the value of 
     * the first feature; in the latter case, the values of the first feature 
     * must countable.
     */
    virtual inline void updParamV(double updVal, double featVal, 
				  FeatureVector **params,
				  int fConjInd = 0) {

      (*params[fConjInd])[(int)featVal] += updVal;
    }  

    /**
     * For feature sets
     * @see void Feature::updParamV(double, double, double **, int)
     */
    virtual inline void updParamV(double updVal, double *featValues,
				  FeatureVector **params,
				  int fConjInd = 0) {

      (*params[fConjInd])[(int)featValues[0]] += updVal;
    }  

    /**
     * A member function that prints the parameters associated with this
     * Feature object.
     */
    virtual inline void printParams(std::ostream &ost, const FeatureVector ***par) {
      int i,j;

      if(paramInd() < 0) 
        return;

      Pair<ULONG> p = maxNumParams();
      ostringstream osst;

      for(i = 0; i < p.f; i++) 
	par[paramInd()][i]->print(osst, i);

      if(!osst.str().length())
        return;

      ost << "\tFeature " << ind() << " " << getName() << "\n";
      ost << "\tnum[" << p.f << "," << p.s << "]" << endl;
      ost  << osst.str();
      ost << endl;
    }


    virtual ~Feature() {
      delete [] _preCompArrayIndex;
    }    

  };

  
  /**
   * The class TypedFeature<TClass> is a subtype of Feature and is used as a 
   * wrapper of the latter to introduce type variables in the values of the
   * Feature objects.
   */
  template <class TClass> class TypedFeature : public Feature {

   public:
    /**
     * Default constructor
     * @param fInd A feature index.
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * @param name A unique feature identifier.
     * @param fe A pointer to a FilterEngine object.
     * @param maxNumFeatVals The number of possible feature values. 
     * If greater than zero then the Feature object's possible values are
     * countable.
     * @param type the feature value type, one of TValType enumerate type
     */
    TypedFeature(
                 int fInd, 
                 int paramInd, 
                 int parsingFrames, 
                 char *name, 
                 FilterEngine *fe, 
                 ULONG maxNumFeatVals = 0, 
                 TValType type = FT_INTEGER
                 ) 
      :   Feature(fInd, paramInd,
                  name, fe, maxNumFeatVals 
                  ) {

      this->type = type;
      initPreCompArrays(parsingFrames);
    }

    /**
     * Constructor from a Header string definition. The feature value type
     * is provided as parameter.
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
     */
    TypedFeature(
                 int fInd, 
                 int paramInd, 
                 int parsingFrames, 
                 std::vector<std::string> &fargs, 
                 int & offset, 
                 FilterEngine *fe, 
                 ULONG maxNumFeatVals, 
                 TValType type
                 ) : 
      Feature(fInd, paramInd, fargs, 
              offset, fe, maxNumFeatVals) {

      offset++; //uniqueClsId
      this->type = type;
      initPreCompArrays(parsingFrames);

    }

    /**
     * Constructor from a Header string definition. The parsingFrames is 
     * part of the Header and not provided as parameter.
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
    TypedFeature(
                 int fInd, 
                 int paramInd, 
                 std::vector<std::string> &fargs, 
                 int & offset, 
                 FilterEngine *fe, 
                 ULONG maxNumFeatVals, 
                 TValType type
                 )
      :   Feature(fInd, paramInd, fargs, 
                  offset, fe, maxNumFeatVals
                  ) {

      int pframes;

      offset++; //uniqueClsId
      this->type = type;
      if(!sscanf(fargs[offset++].c_str(), "%d", &pframes))
        assert(0);
      initPreCompArrays(pframes);
    }

    /**
     * Constructor from a Header string definition. The feature's value type 
     * needs to be determined from information in the Header.
     * is provided as parameter.
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
    TypedFeature(
                 int fInd, 
                 int paramInd, 
                 int parsingFrames, 
                 std::vector<std::string> &fargs, 
                 int & offset, 
                 FilterEngine *fe, 
                 ULONG maxNumFeatVals
                 ) : 
      Feature(fInd, paramInd, fargs, 
              offset, fe, maxNumFeatVals) {

      std::string & uClsId = fargs[offset++];
      string::size_type index =  uClsId.find('<', 0);
      assert(index != std::string::npos);
      std::string t = uClsId.substr(index + 1, uClsId.length() - index - 2);
      type = Utils::stringToTValType(t);

      initPreCompArrays(parsingFrames);
    }

    /**
     * Constructor from a Header string definition. Both, the parsingFrames 
     * and feature's value type need to be determined from information in the
     * Header.
     * is provided as parameter.
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
    TypedFeature(
                 int fInd, 
                 int paramInd, 
                 std::vector<std::string> &fargs, 
                 int & offset, 
                 FilterEngine *fe, 
                 ULONG maxNumFeatVals
                 ) : 
      Feature(fInd, paramInd, fargs,
              offset, fe, maxNumFeatVals) {

      std::string & uClsId = fargs[offset++];
      type = Utils::extractValTypeFromClassName(uClsId);

      int pframes;

      if(!sscanf(fargs[offset++].c_str(), "%d", &pframes))
        assert(0);
      initPreCompArrays(pframes);

    }


    ~TypedFeature() {
      ;
    }

  };  

  template <class TClass> class FT_FilterWrapper : public TypedFeature<TClass> {
   protected:
    TypedFilter<TClass> *f1;
   public:
    FT_FilterWrapper(
		  int fInd, 
		  int paramInd, 
		  char *name, 
		  FilterEngine *fe,
		  TypedFilter<TClass>* f1)
      : TypedFeature<TClass>(fInd, paramInd, 
                             1,
                             name, fe,
                             f1->maxNumFilterValues()) {
      
      this->f1 = f1;
      
    }
    
    FT_FilterWrapper(
		  int fInd, 
		  int paramInd, 
		  vector<std::string> & fargs, 
		  int & offset, 
		  FilterEngine *fe, 
		  FeatureEngine *fte)
      : TypedFeature<TClass>(fInd, 
                             paramInd,                             
			     1,
                             fargs,
                             offset, 
                             fe,
			     ((TypedFilter<TClass> *)
			      fe->getFilter(fargs[offset + 2]))->maxNumFilterValues()) {

      this->f1 = (TypedFilter<TClass> *)fe->getFilter(fargs[offset++]);
      
    }
    
    inline double featValue(Tag *ge, int frame) {
      double ret_val = (double)f1->value(ge->getPos(), ge->getStrand());
      //      cerr << f1->getName() << " " << ge->getPos() << " " << ret_val << endl;
      return ret_val;
    }

  };

}

namespace std {

  /**
   * Routine for sorting STL containers instantiated with Feature* types
   */
  template <>
  struct less<Feature *> {
    bool operator()(const Feature* f1, const Feature* f2)
    {
      //Defined for ascending sorting
      if(!f1)
        return true;
      if(!f2)
        return false;
      return ((Feature &)(*f1)) < (*f2);
    }
  };

  /**
   * Routine for sorting STL containers instantiated with Feature* types
   */  
  template <>
  struct equal_to<Feature *> {
    bool operator()(const Feature* f1, const Feature* f2)
    {
      //Defined for ascending sorting
      if(!f1 || !f2)
        return false;
      return ((Feature &)(*f1)) == (*f2);
    }
  };
}

#endif

