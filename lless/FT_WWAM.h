/****************************************************************************
* FT_WWAM.h - part of the lless namespace, a general purpose
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

#ifndef FT_WWAM_FEAT_H
#define FT_WWAM_FEAT_H
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Feature.h"
#include "NGram.h"
#include "Consensus.h"
#include "FT_Edge.h"


namespace lless {
  
  /**
   * The FT_WWAM class is a subclass of FT_Edge and provides an 
   * interface for all the WWAM-derived models that are often used for 
   * signal characterization. The class provides general enough 
   * implementations of updParamV and dotParamV and most of the time the
   * only function that derived classes must implement is computeValues(Tag *).
   * The added bonus is that because the WWAM has
   * been implemented as a Feature object, we can train the WWAM parameters
   * as part of the global structure.
   ***************************************************************************/

  template <class TClass>
  class FT_WWAM : public FT_Edge<TClass> {
    
   protected:
    int _step;
    int _offset, _length;
    int _windowSize;
    int domSize;
    int _order;
    int _numPosVals;
    int _numUpdates;
    int *vals;
    
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
     * @param windowSize the size of the window around any position in the 
     * WWAM which is used for counting symbols or words; windowSize/size 
     * must be an odd number.
     * @param domSize the total number of symbols which the sequence may have.
     * @param offset the length of the WWAM, offset of the signal occurrence.
     * @param length  the length of the WWAM, length of the signal occurrence. 
     * The total length of the WWAM would be then offset + length.
     * @param step the number of positions to jump within the WWAM, before
     * counting again.
     */
    FT_WWAM(
            int fInd,
            int paramInd, 
            int parsingFrames, 
            char *name, 
            FilterEngine *fe, 
            FL_Signal<EdgeInst> *signal, 
	    TValType type,
            int windowSize, 
            int domSize, 
            int offset, 
            int length, 
            int step
            )
      : FT_Edge<TClass>(fInd, 
			paramInd, 
			parsingFrames, 
			name, 
			fe, 
			0,
			signal,
			type) {
      
      this->_length = length;
      this->_offset = offset;
      this->_step = step;
      this->_windowSize = windowSize;
      this->domSize = domSize;
      // _windowSize/size must be an odd number
      assert((_windowSize/_step) % 2);
      assert(!(_length%_step));

    }

    /**
     * Constructor from a Header string definition. windowSize should be read
     * from the Header definition.
     * @param fInd A feature index.     
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * object to which this feature is tied to.
     * @param fargs The Header string definition loaded as a vector of strings.
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param domSize the total number of symbols which the sequence may have.
     */
    
    FT_WWAM(
            int fInd, 
            int paramInd, 
            vector<std::string> & fargs, 
            int & offset, 
            FilterEngine *fe, 
            int domSize
            )
      : FT_Edge<TClass>(fInd, 
			paramInd, 
			fargs, offset, 
			fe, 0) {
      
      readCoordinates(fargs, offset);
      
      if(!sscanf(fargs[offset++].c_str(), "%d", &_windowSize))
        assert(0);
      
      this->domSize = domSize;
      // _windowSize/size must be an odd number
      assert((_windowSize/_step) % 2);
      assert(!(_length%_step));

    }
    
    /**
     * Constructor from a Header string definition. windowSize is provided
     * as parameter.
     * @param fInd A feature index.     
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * object to which this feature is tied to.
     * @param fargs The Header string definition loaded as a vector of strings.
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param windowSize the size of the window around any position in the 
     * WWAM which is used for counting symbols or words; windowSize/size 
     * must be an odd number.
     * @param domSize the total number of symbols which the sequence may have.
     */

    FT_WWAM(
            int fInd, 
            int paramInd, 
            vector<std::string> & fargs, 
            int & offset, 
            FilterEngine *fe, 
            int windowSize,
            int domSize
            )
      : FT_Edge<TClass>(fInd, 
			paramInd, 
			fargs, offset, 
			fe, 0) {
      
      readCoordinates(fargs, offset);
      this->_windowSize = windowSize;
      this->domSize = domSize;
      // _windowSize/size must be an odd number
      assert((_windowSize/_step) % 2);
      assert(!(_length%_step));

    }
    
    /**
     * A member function for reading the WWAM coordinate parameters
     * @param fargs The Header string definition loaded as a vector of strings.
     * @param offset The index for vector fargs.
     */
    inline void readCoordinates(vector<std::string> &fargs, int & offset) {
      if(!sscanf(fargs[offset++].c_str(), "%d", &_offset))
        assert(0);
      if(!sscanf(fargs[offset++].c_str(), "%d", &_length))
        assert(0);
      if(!sscanf(fargs[offset++].c_str(), "%d", &_step))
        assert(0);
    }
    
    inline bool isFeatureSet() {
      return true;
    }

    inline int windowSize() {
      return _windowSize;
    }
    
    inline int length() {
      return _length;
    }
    
    inline int offset() {
      return _offset;
    }
    
    inline int domainSize() {
      return domSize;
    }
    
    inline int step() {
      return  _step;
    }
    
    inline int order() { 
      return _order;
    }
    
    /**                                                                                  
     * @return the size array vals                                                    
     */
    inline int numPosVals() {
      return _numPosVals;
    }

    inline int numUpdates() {
      return _numUpdates;
    }

    virtual inline int* computeValues(Tag *ge) {
      throw EXCEPTION( FORB_VIRTUAL_FUNCTION, this->getName());
      return NULL;
    }

    /**
     * This function is unaltered in any derived class from FT_WWAM.
     * (except for BinnedWWAM)
     * @param ge Tag object to which this feature is tied.
     * @param frame Tag object's current phase(or frame).
     * @param params parameter vector to be used for computing the dot 
     * product.
     * @param fConjInd is an offset for the array **params, which is
     * different from zero when this feature has been used in a Feature 
     * conjunction or Feature disjunction operation (defined as the class
     * FeaturexFeature, in which case fConjInd usually denotes the value of 
     * the first feature; in the latter case, the values of the first feature 
     * must countable.
     * @return the dot product of WWAM and its corresponding parameter(s)
     * \todo function is not working with filter parameter for Conjunction
     * features as of now. This behaviour is mimicked by PWMxPWM and PWMUnion 
     * features so it is not needed.. maybe
     */    
    virtual inline double dotParamV(Tag *ge, int frame, 
                                    const FeatureVector **params,
				    int fConjInd = 0,
				    TypedFilter<UCHAR> *filter = NULL) {      
      
      assert(!filter);

      double result = 0;
      computeValues(ge);
      const FeatureVector &mpar = *params[fConjInd];

      //      cerr << this->getName() << " " << ge->getPos() << endl;
      for(register int i = 0; i < this->numUpdates(); i++) {
	assert(this->vals[i] < this->_numParams.s);	       
	result += mpar[this->vals[i]];
	//	cerr << "i " << i << " " << this->vals[i] << " " << mpar[this->vals[i]] << " " << result << endl;
      }
      
      /*
      double result = 0;
      computeValues(ge);
      int offset;
      TClass *posVals;
      //     cerr << this->getName();
      for(register int i = 0; i < length(); i += step()) {
vals[j] << endl;                                                                                 posVals = posValues(ge, i);
        offset = i*domainSize();
	
        for(register int j = 0; j < this->numPosVals(); j++) {
          //          cerr << "\t" << i << " " << j << " " << vals[j] << " " << offset + vals[j] << endl;
          result += (*params[fConjInd])[offset + posVals[j]];
          
        }
	}*/
      return result;
    }
    
    /**
     * This function is unaltered in all derived classes from FT_WWAM.
     * (except for BinnedWWAM)
     * @param updVal the value of the update to be made to this feature's
     * corresponding parameter(s)
     * @param ge Tag object to which this feature is tied.
     * @param frame Tag object's current phase(or frame).
     * @param params parameter vector to be used for computing the dot 
     * product.
     * @param fConjInd is an offset for the array **params, which is
     * different from zero when this feature has been used in a Feature 
     * conjunction or Feature disjunction operation (defined as the class
     * FeaturexFeature, in which case fConjInd usually denotes the value of 
     * the first feature; in the latter case, the values of the first feature 
     * must countable.
     * \todo function is not working with filter parameter for Conjunction
     * features as of now. This behaviour is mimicked by PWMxPWM and PWMUnion 
     * features so it is not needed.. maybe
     */
    virtual inline void updParamV(double updVal, Tag *ge, 
				  int frame, FeatureVector **params, int fConjInd = 0,
				  TypedFilter<UCHAR> *filter = NULL) {      
      
      assert(!filter);      

      computeValues(ge);
      FeatureVector &mpar = *params[fConjInd];

      //      cerr << this->getName() << " " << fConjInd << " " << ge->getStrand() << " " << ge->getPos() << endl;
      for(register int i = 0; i < this->numUpdates(); i++) {
	assert(this->vals[i] < this->_numParams.s);
	//	if(!this->getName().compare("Start-Peptide-WWAM")) 
	//	  cerr << "updating " << this->vals[i] << " " << mpar[this->vals[i]] << endl;
        mpar[this->vals[i]] += updVal;
      }

      /*      int offset;
      computeValues(ge);
      TClass *posVals;
      
      for(int i = 0; i < length(); i += step()) {
        posVals = posValues(ge, i);
        offset = i*domainSize();
        
        for(int j = 0; j < this->numPosVals(); j++) 
          (*params[fConjInd])[offset + posVals[j]] += updVal;
	  } */
    }

    
    /**
     * A member function for displaying the  WWAM's computed parameters.
     * @see Feature::printParams()
     */

    /*    virtual inline void printParams(std::ostream &ost, const FeatureVector ***par) {
      if(this->paramInd() < 0)
        return;

      printParams(ost, par[this->paramInd()]);
      }*/

    inline void printParams(std::ostream &ost,
			    const FeatureVector **par, int offset = 0) {

      int i, j, k, limit;   
      Pair<ULONG> params = this->maxNumParams();
      
      ost << "\tFeature " << this->ind() << " " << this->getName() << "\n";
      ost << "\tnum[" << domainSize() << "," << length() << "]" << endl;
      
      double *cols = new double[length()];
      
      ost.precision(10);
      for(int phase = 0; phase < params.f; phase++) {
        ost << "\tPhase " << phase << endl;
        for(i = 0; i < length(); i++)
          cols[i] = 0;
        
        for(i = 0; i < length() - order() + 1; i += 10) {
	  ost << "\n*";
          limit = Utils::min(i + 10, length() - order() + 1);
          
          for(j = i; j < limit; j += step())
            ost << "\t" << j;
          
          for(j = 0; j < domainSize(); j++) {
            for(k = i; k < limit; k += step()) 
              cols[k] += (*par[phase])[offset + k*domainSize() + j];
          }
          
          for(j = 0; j < domainSize(); j++) {
            ost << "\n" << j;
            
            for(k = i; k < limit; k += step()) {
              //            ost << "\t" << par[paramInd()][phase][i*domainSize() + j]/cols[i]*100.0;
              ost << "\t" << (*par[phase])[offset + k*domainSize() + j];
            }
          }
          ost << endl;
        }
        ost << endl;
      }
      delete [] cols;
    } 
    
    ~FT_WWAM() {
      delete [] this->vals;
      this->vals = NULL;      
    }

  };
  

  /**
   * The FT_BinnedWWAM class is a subclass of FT_WWAM and providesa
   * binning implementation over filter value counts at each position
   * in the WWAM.
   ***************************************************************************/

  template <class TClass1, class TClass2>
  class FT_BinnedWWAM : public FT_WWAM<TClass2> { 
   protected:
    FT_BaseBin<TClass1> *binObj;
    SPARSE_HASH<TClass2, int> **counts;
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
     * @param windowSize the size of the window around any position in the 
     * WWAM which is used for counting symbols or words; windowSize/size 
     * must be an odd number.
     * @param domSize the total number of symbols which the sequence may have.
     * @param offset the length of the WWAM, offset of the signal occurrence.
     * @param length  the length of the WWAM, length of the signal occurrence. 
     * The total length of the WWAM would be then offset + length.
     * @param step the number of positions to jump within the WWAM, before
     * counting again.
     */
    FT_BinnedWWAM(
		  int fInd,
		  int paramInd, 
		  int parsingFrames, 
		  char *name, 
		  FilterEngine *fe, 
		  FL_Signal<EdgeInst> *signal, 
		  TValType type,
		  int windowSize, 
		  int domSize, 
		  int offset, 
		  int length, 
		  int step,
		  FT_BaseBin<TClass1> *binObj)
      : FT_WWAM<UCHAR>(fInd, 
		       paramInd,
		       parsingFrames,
		       name,
		       fe,
		       signal,
		       type,
		       windowSize,
		       domSize,
		       offset,
		       length,
		       step) {
      
      this->binObj = binObj;
    }

    /**
     * Constructor from a Header string definition. windowSize should be read
     * from the Header definition.
     * @param fInd A feature index.     
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * object to which this feature is tied to.
     * @param fargs The Header string definition loaded as a vector of strings.
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param domSize the total number of symbols which the sequence may have.
     */
    
    FT_BinnedWWAM(
		  int fInd, 
		  int paramInd, 
		  vector<std::string> & fargs, 
		  int & offset, 
		  FilterEngine *fe, 
		  FeatureEngine *fte, 
		  int domSize
		  )
      : FT_WWAM<TClass2>(fInd, 
			 paramInd,
			 fargs, offset,
			 fe,
			 domSize) {
      
      this->binObj = (FT_BaseBin<TClass1> *)fte->getFeature(fargs[offset++]);

    }
    
    /**
     * Constructor from a Header string definition. windowSize is provided
     * as parameter.
     * @param fInd A feature index.     
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * object to which this feature is tied to.
     * @param fargs The Header string definition loaded as a vector of strings.
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param windowSize the size of the window around any position in the 
     * WWAM which is used for counting symbols or words; windowSize/size 
     * must be an odd number.
     * @param domSize the total number of symbols which the sequence may have.
     */

    FT_BinnedWWAM(
            int fInd, 
            int paramInd, 
            vector<std::string> & fargs, 
            int & offset, 
            FilterEngine *fe, 
	    FeatureEngine *fte,
            int windowSize,
            int domSize)
      : FT_WWAM<TClass2>(fInd, 
			 paramInd, 
			 fargs, offset, 
			 fe, windowSize, 
			 domSize) {
      
      this->binObj = (FT_BaseBin<TClass1> *)fte->getFeature(fargs[offset++]);

    }
    

    /**
     * @param ge Tag object to which this feature is tied.
     * @param frame Tag object's current phase(or frame).
     * @param params parameter vector to be used for computing the dot 
     * product.
     * @param fConjInd is an offset for the array **params, which is
     * different from zero when this feature has been used in a Feature 
     * conjunction or Feature disjunction operation (defined as the class
     * FeaturexFeature, in which case fConjInd usually denotes the value of 
     * the first feature; in the latter case, the values of the first feature 
     * must countable.
     * @return the dot product of WWAM and its corresponding parameter(s)
     * \todo function is not working with filter parameter for Conjunction
     * features as of now. This behaviour is mimicked by PWMxPWM and PWMUnion 
     * features so it is not needed.. maybe
     */    
    inline double dotParamV(Tag *ge, int frame, 
			    const FeatureVector **params,
			    int fConjInd = 0,
			    TypedFilter<UCHAR> *filter = NULL) {      
      
      assert(!filter);
      
      double result = 0;
      this->computeValues(ge);
      typename SPARSE_HASH<TClass2, int>::iterator it;

      fConjInd = fConjInd*this->numUpdates()*this->domainSize()/this->numPosVals();
      
      for(register int i = 0; i < this->numUpdates(); i++) {
	it = counts[i]->begin();

	for( ; it != counts[i]->end(); it++) 
	  result += binObj->dotParamV((double)it->second, 
				      params, fConjInd + this->vals[i] + it->first);
      }
      
      return result;
    }
    
    /**
     * This function is unaltered in all derived classes from FT_WWAM.
     * @param updVal the value of the update to be made to this feature's
     * corresponding parameter(s)
     * @param ge Tag object to which this feature is tied.
     * @param frame Tag object's current phase(or frame).
     * @param params parameter vector to be used for computing the dot 
     * product.
     * @param fConjInd is an offset for the array **params, which is
     * different from zero when this feature has been used in a Feature 
     * conjunction or Feature disjunction operation (defined as the class
     * FeaturexFeature, in which case fConjInd usually denotes the value of 
     * the first feature; in the latter case, the values of the first feature 
     * must countable.
     * \todo function is not working with filter parameter for Conjunction
     * features as of now. This behaviour is mimicked by PWMxPWM and PWMUnion 
     * features so it is not needed.. maybe
     */
    virtual inline void updParamV(double updVal, Tag *ge, 
				  int frame, FeatureVector **params, int fConjInd = 0,
				  TypedFilter<UCHAR> *filter = NULL) {      
      
      assert(!filter);      
      
      this->computeValues(ge);
      typename SPARSE_HASH<TClass2, int>::iterator it;

      fConjInd = fConjInd*this->numUpdates()*this->domainSize()/this->numPosVals();
      
      //      cerr << this->getName() << " " << fConjInd << " " << ge->getStrand() << " " << ge->getPos() << endl;
      for(register int i = 0; i < this->numUpdates(); i++) {
	it = counts[i]->begin();
	//	cerr << "\t" << i << " ";
	for( ; it != counts[i]->end(); it++)  {
	  //	  cerr << "\t\t" << fConjInd << " " << this->vals[i] << " " << (int)it->first << " " << it->second;
	  binObj->updParamV(updVal, (double)it->second, 
			    params, fConjInd + this->vals[i] + it->first);
	}
	//	cerr << endl;
      }
    }
    
    /**
     * A member function for displaying the  WWAM's computed parameters.
     * @see Feature::printParams()
     */
    /*    virtual inline void printParams(std::ostream &ost, const FeatureVector ***par) {
      Feature::printParams(ost, par);
      }*/

    ~FT_BinnedWWAM() {
      delete [] this->counts;
      this->counts = NULL;      
    }

  };
  
  
  /**
   * FT_WWAMUnion is a subtype of FT_WWAM that concatenates the values of
   * FT_WWAM objects at each signal position.
   **************************************************************************/
  template <class TClass>
  class FT_WWAMUnion : public FT_WWAM<TClass> {
    
   protected:
    vector<FT_WWAM<TClass> *> _wwams;

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
    FT_WWAMUnion(
		 int fInd, 
		 int paramInd, 
		 int parsingFrames,
		 char *name,
		 FilterEngine *fe,
		 FL_Signal<EdgeInst> *signal, 
		 int windowSize, 
		 int step,
		 vector<FT_WWAM<TClass> *> & wwams)
      : FT_WWAM<TClass>(fInd,
			paramInd,
			parsingFrames,
			name,
			fe,
			signal,
			windowSize,
			wwams[0]->domainSize(),
			0, 0,
			step) {

      this->_wwams = wwams;
      initMembers();
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
     * Feature name FT_WWAMUnion parsingFrames signal offset length step windowSize 
     * gram \n\n
     * The description of the fields could be found in 
     * the other constructor(s)
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param fte A pointer to a FeatureEngine object.
     */
    FT_WWAMUnion(  
		 int fInd, 
		 int paramInd, 
		 vector<std::string> & fargs, 
		 int & offset, 
		 FilterEngine *fe, 
		 FeatureEngine *fte)
      : FT_WWAM<TClass>(fInd,
			paramInd,
			fargs, offset,
			fe,
	       ((FT_WWAM<TClass> *)fte->getFeature(fargs[offset + 9]))->domainSize()) {
      
      FT_WWAM<TClass> *w;
      while(offset < fargs.size()) {
	w = (FT_WWAM<TClass> *)fte->getFeature(fargs[offset++]);
	_wwams.push_back(w);
      }
      initMembers();
    }
    
    FT_WWAMUnion(  
		 int fInd, 
		 int paramInd, 
		 vector<std::string> & fargs, 
		 int & offset, 
		 FilterEngine *fe, 
		 FeatureEngine *fte,
		 int windowSize)
      : FT_WWAM<TClass>(fInd,
			paramInd,
			fargs, offset,
			fe, windowSize,
		((FT_WWAM<TClass> *)fte->getFeature(fargs[offset + 8]))->domainSize()) {
	
      FT_WWAM<TClass> *w;
      while(offset < fargs.size()) {
	w = (FT_WWAM<TClass> *)fte->getFeature(fargs[offset++]);
	_wwams.push_back(w);
      }
      initMembers();
    }

    inline void initMembers() {
      assert(_wwams.size());
      this->_order = _wwams[0]->order();
      this->_length = -1;
      this->_offset = -1;
      this->_numUpdates = 0;
      this->_numParams.s = 0;

      for(int i = 0; i < _wwams.size(); i++) {
	this->_numUpdates += _wwams[i]->numUpdates();
	this->_numParams.s += _wwams[i]->numUpdates()*_wwams[i]->domainSize();
      }

      this->_numPosVals = 1;
      this->vals = new int [this->numUpdates()];
    }
      
    /**
     * A member function that computes protected member vals, using the vals
     * from contained WWAM *objects;
     * @return vals
     */
    inline int* computeValues(Tag *ge) {
      int base = 0;
      int offset = 0;
      int size = _wwams.size();

      for(register int i = 0; i < size; i++) {
	FT_WWAM<TClass> *w = _wwams[i];
	int *wvals = w->computeValues(ge);

	for(register int j = base; j < base + w->numUpdates(); j++)
	  this->vals[j] = offset + wvals[j - base]; 

	offset += w->numUpdates()*w->domainSize();
	base += w->numUpdates();
      }

      assert(base == this->numUpdates());
      return this->vals;
    }

    /*    virtual inline void printParams(std::ostream &ost, const FeatureVector ***par) {
      int offset = 0;

      for(int i = 0; i < _wwams.size(); i++) {
    	_wwams[i]->printParams(ost, par[this->paramInd()], offset);
	offset += _wwams[i]->numUpdates()*_wwams[i]->domainSize();
      }
      }*/

    ~FT_WWAMUnion() {
    }
    
  };

}

#endif
