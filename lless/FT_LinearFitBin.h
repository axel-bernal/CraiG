/****************************************************************************
* FT_LinearFitBin.h - part of the lless namespace, a general purpose
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

#ifndef FT_LINEARFITBIN_FEATURE_H
#define FT_LINEARFITBIN_FEATURE_H
  
#include "FT_Bin.h"

namespace  lless {

  /**
   * FT_LinearFitBin derives FT_Bin. It creates a set of equidistant ordered 
   * disjoint control points which cover the range of values of the contained
   * feature and uses them to fit a linear function to the data between each 
   * of them. The contained feature's value is used to find out between which
   * control points the value falls in.
   * This class is useful to model contained features with sparse or 
   * multimodal distributions, such as lengths and scores.
   *
   ***************************************************************************/

  template <class TClass> class FT_LinearFitBin : public FT_Bin<TClass> {
  
   protected:

    
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
    FT_LinearFitBin(
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
                       maxNumFeatVals, 
                       rangeStart,
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

    FT_LinearFitBin(
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
      
    }

    int lastFloor() {
      return this->maxNumFeatValues() - 2;
    }

    int featValCs(double featVal, int *C) {
      C[0] = this->floor(featVal);
      
      if(this->outOfRange)
	return 1;

      C[1] = this->ceiling(featVal, C[0]);
      
      assert(C[0] != C[1]);
      return 2;
    }
    
    void featValVs(double featVal, const int *C,
		   double *V, int numCPUpdates) {
      
      if(this->outOfRange) {
	assert(numCPUpdates == 1);
	V[0] = 1.0;
	this->outOfRange = false;
	return;
      }
      
      //      cerr << this->getName() << " " << featVal << " " << this->cpValue(C[0])
      //      	   << " " << this->cpValue(C[1]) << " " << this->binWidth(C[0]) << endl;

      assert(numCPUpdates == 2);

      V[0] = (this->cpValue(C[1]) - featVal)/this->binWidth(C[0]);
      V[1] = (featVal - this->cpValue(C[0]))/this->binWidth(C[0]);
      
    }

    ~FT_LinearFitBin() {
  
    }

  };


  /**
   * FT_LinearFitCustomBin inherits from FT_CustomBin and creates a set of 
   * arbitrarily spaced and ordered disjoint control points which cover the
   * range of values of the contained feature  and uses them to fit a linear
   * function of the data between each of them. The contained feature's value
   * is used to find out between which. The contained feature's value is used
   * to find out between which control points the value falls in.
   * This class is useful to model contained features with sparse or 
   * multimodal distributions, such as lengths and scores.
   * \todo make the interpolation work for non-equidistant points
   *
   ***************************************************************************/

  template <class TClass> class FT_LinearFitCustomBin : public FT_CustomBin<TClass> {

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
    FT_LinearFitCustomBin(
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

    FT_LinearFitCustomBin(
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
      this->_numParams.s++;

    }
    
    int lastFloor() {
      return this->maxNumFeatValues() - 2;
    }

    int featValCs(double featVal, int *C) {

      C[0] = this->floor(featVal);

      if(this->outOfRange)
	return 1;

      C[1] = this->ceiling(featVal, C[0]);

      assert(C[0] != C[1]);
      return 2;
    }
    
    void featValVs(double featVal, const int *C,
		   double *V, int numCPUpdates) {

      if(this->outOfRange) {
	assert(numCPUpdates == 1);
	V[0] = 1.0;
	this->outOfRange = false;
	return;
      }
      
      assert(numCPUpdates == 2);
      
      V[0] = (this->cpValue(C[1]) - featVal)/this->binWidth(C[0]);
      V[1] = (featVal - this->cpValue(C[0]))/this->binWidth(C[0]);
      
    }

    ~FT_LinearFitCustomBin() {
  
    }

  };
}

#endif
