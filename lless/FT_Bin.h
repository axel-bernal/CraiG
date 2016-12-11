/****************************************************************************
* FT_Bin.h - part of the lless namespace, a general purpose
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

#ifndef FT_BIN_FEATURE_H
#define FT_BIN_FEATURE_H
  
#include "Feature.h"
#include "Filter.h"
#include "FilterEngine.h"
#include "FeatureEngine.h"

#define MAX_FEATVAL_UPDATES 5

namespace  lless {
  
  /**
   * FT_BaseBin is the abstract base class for all types of bin features.
   * Bin features create a set of ordered disjoint bins which cover the range
   * of values of the contained feature. The contained feature's value is used
   * to find out the bin which is 'on', i.e. the bin in which the value of 
   * the contained feature falls in.
   * This class is useful to model contained features with sparse or 
   * multimodal distributions, such as lengths and scores.
   *
   ***************************************************************************/

  template <class TClass> class FT_BaseBin : public TypedFeature<TClass> {
  
   protected:
    TypedFeature<TClass> *ft;
    bool outOfRange;
    int C[MAX_FEATVAL_UPDATES];
    double V[MAX_FEATVAL_UPDATES];
    
   public:
    /**
     * Default constructor
     * @param fInd A feature index.
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * @param name A unique feature identifier.
     * @param fe A pointer to a FilterEngine object.
     * @param ft A pointer to a TypedFeature<TClass>, the contained feature.
     * @param maxNumFeatVals The number of bins that are created. 
     * @param type the feature value type, one of TValType enumerate type
     */
    FT_BaseBin(
	       int fInd,
	       int paramInd, 
	       int parsingFrames,
	       char *name, 
	       FilterEngine *fe,
	       TypedFeature<TClass> *ft,
	       int maxNumFeatVals,
	       TValType type = FT_INTEGER
	       ) 
      : TypedFeature<TClass>(fInd, paramInd, 
                             parsingFrames,
                             name, fe,
                             maxNumFeatVals, type) {
      
      this->ft = ft;
      this->outOfRange = false;
    }
  
    /**
     * Constructor from a Header string definition. The feature value type
     * is provided as parameter. 
     * @see FT_Bin for Header details
     *
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param fte A pointer to a FeatureEngine object.
     */

    FT_BaseBin(
           int fInd, 
           int paramInd, 
           int parsingFrames,
           vector<std::string> & fargs, 
           int & offset, 
           FilterEngine *fe, 
           FeatureEngine *fte
	       )
      : TypedFeature<TClass>(fInd, paramInd,
                             parsingFrames,
                             fargs, offset,
                             fe, 0) {
      
      this->ft = NULL;
      this->outOfRange = false;

      if(fargs[offset++].compare("null")) {
        ft = (TypedFeature<TClass> *)fte->getFeature(fargs[offset - 1]);
        assert(this->type == ft->valType());
      }
      
      if(!sscanf(fargs[offset++].c_str(), "%d", &this->_numParams.s))
        assert(0);            

    }

    virtual int featValCs(double featVal, int *C) {
      C[0] = floor(featVal);
      return 1;
    }
    
    virtual void featValVs(double featVal, const int *C,
			   double *V, int numCPUpdates) {
      
      assert(numCPUpdates == 1);
      V[0] = 1.0;
      //      cerr << this->getName() << " " << featVal << " " << this->cpValue(C[0])
      //	   << " " << this->binWidth(C[0]) << endl;
      outOfRange = false;
    }

    /**
     * A virtual member function that returns the value of the last valid
     * floor for this bin arrangement
     */
    virtual int lastFloor() {
      return this->maxNumFeatValues() - 1;
    }

    /**
     * A virtual member function that computes the floor control point closest
     * to featVal
     */
    virtual int floor(double featVal) = 0;

    /**
     * A virtual member function that computes the value corresponding to cp
     * the control point passed as parameter
     */
    virtual double cpValue(int cp) = 0;
 
    /**
     * A virtual member function that returns the binWidth starting at control
     * point cp
     */
    virtual double binWidth(int cp = 0) = 0;

    /**
     * A member function that computes the ceiling control point closest to featVal
     */
    inline int ceiling(double featVal, int & floor) {

      if(outOfRange || floor > lastFloor()) {
	assert(0);
	throw EXCEPTION(BAD_USAGE, "can't call ceiling if point out of range");
      }
      return floor + 1;
    }    

    
    /**
     * @param ge Tag object to which this feature is tied.
     * @param frame Tag object's current phase(or frame).
     * @return the contained feature value 
     * @see TypedFeature<int>::featValue(Tag *, int)
     */
    inline double featValue(Tag *ge, int frame) {
      assert(ft);
      return ft->featValue(ge, frame);
    }

    /**
     * @see FT_BaseBin::dotParamV(double, double **, int)
     */
    inline double dotParamV(double featVal, int fConjInd = 0) {
      return dotParamV(featVal, this->cparams[this->paramInd()], fConjInd);
    }
  
    /**
     * @see FT_BaseBin::dotParamV(double, double **, int)
     */
    inline double dotParamV(Tag *ge, int frame, 
                            const FeatureVector **params, int fConjInd = 0,
			    TypedFilter<UCHAR> *filter = NULL) {

      assert(!filter);
      return dotParamV(ft->featValue(ge, frame), params, fConjInd);
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
    inline double dotParamV(double featVal, 
			    const FeatureVector **params,
			    int fConjInd = 0) {
      
      int numCPUpdates = featValCs(featVal, C);
      featValVs(featVal, C, V, numCPUpdates);
      double result = 0;
      //      cerr << this->getName() << endl;
      for(int i = 0; i < numCPUpdates; i++) {
	//	cerr << V[i] << " " << C[i] << " " << (*params[fConjInd])[C[i]] << endl;
	result += V[i]*(*params[fConjInd])[C[i]];
      }
      return result;

    }
    
    /**
     * @see FT_BaseBin::updParamV(double, double, double **, int)
     */    
    inline void updParamV(double updVal, double featVal, int fConjInd = 0) {
      updParamV(updVal, featVal, this->params[this->paramInd()], fConjInd);
    }
  
    /**
     * @see FT_BaseBin::updParamV(double, double, double **, int)
     */    
    inline void updParamV(double updVal, Tag *ge, int frame,
			  FeatureVector **params, int fConjInd = 0,
			  TypedFilter<UCHAR> *filter = NULL) {

      assert(!filter);
      updParamV(updVal, ft->featValue(ge, frame), params, fConjInd);
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
    inline void updParamV(double updVal, double featVal,
			  FeatureVector **params, int fConjInd = 0) {
      
      int numCPUpdates = featValCs(featVal, C);
      featValVs(featVal, C, V, numCPUpdates);

      for(int i = 0; i < numCPUpdates; i++)
	(*params[fConjInd])[C[i]] += V[i]*updVal;
      
    }

    virtual ~FT_BaseBin() {
  
    }

  };


  /**
   * FT_LogScaledBin inherits from FT_Bin and creates a set of ordered
   * disjoint  bins which cover the range of values of the contained feature
   * in logarithmic scale. The logarithm of the contained feature's value is
   * used to find out the bin which is 'on', i.e. the bin in which the log 
   * value of the contained feature falls in.
   * This class is useful to model contained features with sparse or 
   * multimodal distributions, such as lengths and scores.

  /**
   * FT_Bin derives FT_BaseBin. It creates a set of equidistant ordered 
   * disjoint bins which cover the range of values of the contained feature.
   * The contained feature's value is used to find out the bin which is 'on',
   * i.e. the bin in which the value of the contained feature falls in.
   * This class is useful to model contained features with sparse or 
   * multimodal distributions, such as lengths and scores.
   *
   ***************************************************************************/

  template <class TClass> class FT_Bin : public FT_BaseBin<TClass> {
  
   protected:
    double _rangeStart;
    double _binWidth;
    double _logBase;    

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
     * @param rangeStart Where binning starts
     * @param binWidth The width of each bin
     * @param logBase the base for computing the logarithms
     * @param type the feature value type, one of TValType enumerate type
     */
    FT_Bin(
        int fInd,
        int paramInd, 
        char *name, 
        FilterEngine *fe,
        TypedFeature<TClass> *ft,
        int maxNumFeatVals,
        double rangeStart,
        double binWidth, 
	double logBase = 0,
        TValType type = FT_INTEGER
        ) 
      : FT_BaseBin<TClass>(fInd, paramInd, 
                           ft ? ft->parsingFrames() : 1,
                           name, fe, ft,
                           maxNumFeatVals, type) {
      
      this->_rangeStart = rangeStart;
      this->_binWidth = binWidth;
      this->_logBase = logBase;
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
     * _binWidth \n\n
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

    FT_Bin(
           int fInd, 
           int paramInd, 
           vector<std::string> & fargs, 
           int & offset, 
           FilterEngine *fe, 
           FeatureEngine *fte
           )
      : FT_BaseBin<TClass>(fInd, paramInd,
                           (fargs[offset + 2].compare("null") == 0) ? 1 : 
                           fte->getFeature(fargs[offset + 2])->parsingFrames(),
                           fargs, offset,
                           fe, fte) {

      if(!sscanf(fargs[offset++].c_str(), "%lf", &_logBase))
        assert(0);      
      if(!sscanf(fargs[offset++].c_str(), "%lf", &_rangeStart))
        assert(0);
      if(!sscanf(fargs[offset++].c_str(), "%lf", &_binWidth))
        assert(0);    

    }

    inline double cpValue(int cp) { 
      double times = cp;
      if(cp && logBase() > 1)
	times = exp(log(logBase())*(cp - 1));
      
      return times*_binWidth + this->rangeStart();

    }

    inline double rangeStart() {
      return _rangeStart;
    }

    inline double binWidth(int cp = 0) {
      double times = 1;

      if(cp && logBase() > 1)
	times = exp(log(logBase())*cp) - exp(log(logBase())*(cp - 1));
      
      return times*_binWidth;
      
    }
      
    inline double logBase() {
      return _logBase;
    }

    /**
     * A member function that computes the lowest control point closest
     * to featVal/_binWidth.  If the bin is larger than
     * highestRangeLimit, then highest possible floor is returned
     * @see FT_BaseBin::floor(double)
     */  
    virtual inline int floor(double featVal) {
      int bin = (int)((featVal - this->rangeStart())/_binWidth);

      if(logBase() > 1)
	bin = Utils::logScaleIt(bin, logBase());

      if(bin > this->lastFloor() || featVal < this->rangeStart()) {
	this->outOfRange = true;
	if(bin <= 0)  bin = 0;
	else  bin = this->maxNumFeatValues() - 1;
      }

      return bin;

    }       


    virtual ~FT_Bin() {
  
    }

  };


  /**
   * FT_RBin derives FT_BaseBin. It creates a set of equidistant ordered 
   * disjoint bins which cover the range of values of the contained feature.
   * The contained feature's value is used to find out the bin which is 'on',
   * i.e. the bin in which the value of the contained feature falls in.
   * This class is useful to model contained features with sparse or 
   * multimodal distributions, such as lengths and scores. It differs from 
   * FT_Bin in that, the weights corresponding to each of the bins are updated
   * using the feature values instead of predicate values, i.e. 1's
   *
   ***************************************************************************/
  
  template <class TClass> class FT_RBin : public FT_BaseBin<TClass> {
  
   protected:
    double _rangeStart;
    double _binWidth;
    double _logBase;
    
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
     * @param logBase the base for computing the logarithms
     * @param type the feature value type, one of TValType enumerate type
     */
    FT_RBin(
        int fInd,
        int paramInd, 
        char *name, 
        FilterEngine *fe,
        TypedFeature<TClass> *ft,
        int maxNumFeatVals,
        double rangeStart,
        double binWidth,
	double logBase = 0,
        TValType type = FT_INTEGER
        ) 
      : FT_BaseBin<TClass>(fInd, paramInd, 
                           ft ? ft->parsingFrames() : 1,
                           name, fe, ft,
                           maxNumFeatVals, type) {
      
      this->rangeStart = rangeStart;
      this->binWidth = binWidth;
      this->_logBase = logBase;
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
     * Feature name RBin<TClass> integerValuedFilter numBins _rangeStart 
     * _binWidth mean\n\n
     * The Header above creates a FT_RBin feature whose contained feature is 
     * of type TClass and bins of size _binWidth. The feature can handle values
     * of up to _rangeStart + numBinsx_binWidth, but any value above 
     * numBinsx_binWidth + _rangeStart is extrapolated, depending on the 
     * model.
     *
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param fte A pointer to a FeatureEngine object.
     */

    FT_RBin(
           int fInd, 
           int paramInd, 
           vector<std::string> & fargs, 
           int & offset, 
           FilterEngine *fe, 
           FeatureEngine *fte
           )
      : FT_BaseBin<TClass>(fInd, paramInd,
                           (fargs[offset + 2].compare("null") == 0) ? 1 : 
                           fte->getFeature(fargs[offset + 2])->parsingFrames(),
                           fargs, offset,
                           fe, fte) {

      if(!sscanf(fargs[offset++].c_str(), "%lf", &_rangeStart))
        assert(0);      
      if(!sscanf(fargs[offset++].c_str(), "%lf", &_binWidth))
        assert(0);    

    }

    virtual void featValVs(double featVal, const vector<int> & C,
			   vector<double> & V) {

      assert(C.size() == 0);
      V.push_back(sqrt(featVal));
    }

    inline double cpValue(int cp) { 

      if(cp && logBase() > 1)
	cp = exp(log(logBase())*(cp - 1));

      return cp*_binWidth + this->rangeStart();
    }

    double rangeStart() {
      return _rangeStart;
    }

    double binWidth(int cp = 0) {
      int times = 1;

      if(cp && logBase() > 1)
	times = exp(log(logBase())*cp) - exp(log(logBase())*(cp - 1));

      return times*_binWidth;
    }
    
    inline double logBase() {
      return _logBase;
    }

    /**
     * A member function that computes the lowest control point closest
     * to featVal/_binWidth.  If the bin is larger than
     * highestRangeLimit, then highest possible floor is returned
     * @see FT_BaseBin::floor(double)
     */  

    inline int floor(double featVal) {
      int bin = (int)((featVal - this->rangeStart())/_binWidth);

      if(logBase() > 1) 
	bin = Utils::logScaleIt(bin, logBase());
      
      if(bin > this->lastFloor() || featVal < this->rangeStart()) {
	this->outOfRange = true;
	if(bin <= 0)  bin = 0;
	else  bin = this->maxNumFeatValues() - 1;
      }

      return bin;

    }

    virtual ~FT_RBin() {
  
    }

  };


  /**
   * FT_CustomBin inherits from FT_BaseBin and creates a set of arbitrarily
   * spaced and ordered disjoint bins which cover the range of
   * values of the contained feature. The contained feature's value is used
   * to find out the bin which is 'on', i.e. the bin in which the value of 
   * the contained feature falls in.
   * This class is useful to model contained features with sparse or 
   * multimodal distributions, such as lengths and scores.
   *
   ***************************************************************************/

  template <class TClass> class FT_CustomBin : public FT_BaseBin<TClass> {
  
   protected:
    vector<double> cps;

   public:
    /**
     * Default constructor
     * @param fInd A feature index.
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * @param name A unique feature identifier.
     * @param fe A pointer to a FilterEngine object.
     * @param cps the custom-made control points that define the bins
     * @param type the feature value type, one of TValType enumerate type
     */
    FT_CustomBin(
                 int fInd,
                 int paramInd, 
                 char *name, 
                 FilterEngine *fe,
                 TypedFeature<TClass> *ft,
                 vector<double> &cps,
                 bool binary,
                 TValType type = FT_INTEGER
                 ) 
      : FT_BaseBin<TClass>(fInd, paramInd, 
                           ft ? ft->parsingFrames() : 1, 
                           name, fe, ft,
                           cps.size() - 1, type) {

      if(cps.size() < 2)
	throw EXCEPTION(BAD_USAGE, "# of control points should be at least two");

      this->cps = cps;
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
     * Feature name CustomBin<TClass> integerValuedFilter numBins [bins ..]
     *
     * The Header above creates a FT_CustomBin feature whose contained feature
     * is of type TClass and numBins bins.
     *
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param fte A pointer to a FeatureEngine object.
     */

    FT_CustomBin(
                 int fInd, 
                 int paramInd, 
                 vector<std::string> & fargs, 
                 int & offset, 
                 FilterEngine *fe, 
                 FeatureEngine *fte
                 )
      : FT_BaseBin<TClass>(fInd, paramInd,
                           (fargs[offset + 2].compare("null") == 0) ? 1 : 
                           fte->getFeature(fargs[offset + 2])->parsingFrames(),
                           fargs, offset,
                           fe, fte) {      

      for(int i = 0; i < this->maxNumFeatValues(); i++) {
        double bin;
        if(!sscanf(fargs[offset++].c_str(), "%lf", &bin))
          assert(0);
        cps.push_back(bin);
      }

      if(cps.size() < 2)
	throw EXCEPTION(BAD_USAGE, "# of control points should be at least two");

      this->_numParams.s--;

    }
    
    inline double cpValue(int cp) { 
      return cps[cp];
    }

    inline double binWidth(int cp = 0) {
      return cps[cp+1] - cps[cp];
    }

    /**
     * A member function that computes the current bin 'ON' from
     * the feature value taken as parameter. 
     * @see FT_BaseBin::floor(double)
     */  
    inline int floor(double featVal) {
      int bin = 0;

      for( ; bin < cps.size(); bin++) {
        if(featVal < cps[bin])
          break;
      }

      if(bin == 0 || bin == cps.size()) {
	this->outOfRange = true;

	if(bin == cps.size())
	  return this->maxNumFeatValues() - 1;

      }

      return (bin == 0) ? bin : bin - 1;

    }
    
    virtual ~FT_CustomBin() {
      
    }

  };



  /**
   * FT_CustomRBin inherits from FT_BaseBin and creates a set of arbitrarily
   * spaced and ordered disjoint bins which cover the range of
   * values of the contained feature. The contained feature's value is used
   * to find out the bin which is 'on', i.e. the bin in which the value of 
   * the contained feature falls in.
   * This class is useful to model contained features with sparse or 
   * multimodal distributions, such as lengths and scores.
   *
   ***************************************************************************/

  template <class TClass> class FT_CustomRBin : public FT_BaseBin<TClass> {
  
   protected:
    vector<double> cps;

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
    FT_CustomRBin(
                 int fInd,
                 int paramInd,
                 char *name, 
                 FilterEngine *fe,
                 TypedFeature<TClass> *ft,
                 vector<double> &cps,
                 TValType type = FT_INTEGER
                 ) 
      : FT_BaseBin<TClass>(fInd, paramInd, 
                           ft ? ft->parsingFrames() : 1, 
                           name, fe, ft,
                           cps.size() - 1, type) {

      if(cps.size() < 2)
	throw EXCEPTION(BAD_USAGE, "# of control points should be at least two");
      
      this->cps = cps;
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

    FT_CustomRBin(
                 int fInd, 
                 int paramInd, 
                 vector<std::string> & fargs, 
                 int & offset, 
                 FilterEngine *fe, 
                 FeatureEngine *fte
                 )
      : FT_BaseBin<TClass>(fInd, paramInd,
                           (fargs[offset + 2].compare("null") == 0) ? 1 : 
                           fte->getFeature(fargs[offset + 2])->parsingFrames(),
                           fargs, offset,
                           fe, fte) {      
      
      for(int i = 0; i < this->maxNumFeatValues(); i++) {
        double bin;
        if(!sscanf(fargs[offset++].c_str(), "%lf", &bin))
          assert(0);
        cps.push_back(bin);
      }
      
      if(cps.size() < 2)
	throw EXCEPTION(BAD_USAGE, "# of control points should be at least two");

      this->_numParams.s--;

    }
    
    void featValVs(double featVal, const int *C,
		   double *V, int numCPUpdates) {
      
      V[0] = sqrt(featVal);
    }
    
    inline double cpValue(int cp) { 
      return cps[cp];
    }

    inline double binWidth(int cp = 0) {
      return cps[cp+1] - cps[cp];
    }

    /**
     * A member function that computes the current bin 'ON' from
     * the feature value taken as parameter. 
     * @see FT_BaseBin::floor(double)
     */  
    inline int floor(double featVal) {
      int bin = 0;

      for( ; bin < cps.size(); bin++) {
        if(featVal < cps[bin])
          break;
      }

      if(bin == 0 || bin == cps.size()) {
	this->outOfRange = true;

	if(bin == cps.size())
	  return this->maxNumFeatValues() - 1;
      }

      return (bin - 1 < 0) ? 0 : bin - 1;
      
    }
    
    virtual ~FT_CustomRBin() {
      
    }
  };

}  

#endif
