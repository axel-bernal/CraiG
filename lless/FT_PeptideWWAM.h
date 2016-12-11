/****************************************************************************
* FT_PeptideWWAM.h - part of the lless namespace, a general purpose
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

#ifndef FT_PEPTIDEWWAM_FEAT_H
#define FT_PEPTIDEWWAM_FEAT_H
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "FT_WWAM.h"
#include "NGram.h"
#include "Consensus.h"
#include "FT_Edge.h"


namespace lless {
  
  /**
   * FT_PeptideWWAM a subtype of FT_WWAM for peptide sequences, it behaves
   * like a FT_GramWWAM but each position is a codon and its value is the 
   * codon's aminoacid index which is computed by an aminoacid alphabet.
   *
   **************************************************************************/
  template <class TClass1, class TClass2>
    class FT_PeptideWWAM : public FT_WWAM<UCHAR> {
    
   protected:
    FL_BaseGram<TClass1, TClass2> *gram;;
    AASigma *aa;
    
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
     * @param gram a pointer to a FL_Gram object which defines the words that
     * will be counted at position.
     * @param offset the length of the WWAM, offset of the signal occurrence.
     * @param length  the length of the WWAM, length of the signal occurrence. 
     * The total length of the WWAM would be then offset + length.
     * @param step the number of positions to jump within the WWAM, before
     * counting again.
     */
    FT_PeptideWWAM(
                   int fInd,
                   int paramInd,
                   int parsingFrames,
                   char *name,
                   FilterEngine *fe,
                   FL_Signal<EdgeInst> *signal,
                   int windowSize, 
		   FL_FastGram<TClass1, TClass2> *gram,
                   int offset, 
                   int length, 
                   int step) 
      : FT_WWAM<UCHAR>(fInd, 
		       paramInd,
		       parsingFrames,
		       name,
		       fe,
		       signal, 
		       windowSize, 
		       AA_ALPHABET_SIZE, 
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
     * Feature name PeptideWWAM parsingFrames signal offset length step windowSize
     * gram \n\n
     * The description of the fields could be found in 
     * the other constructor(s)
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param fte A pointer to a FeatureEngine object.
     */    
    FT_PeptideWWAM(
                   int fInd, 
                   int paramInd, 
                   vector<std::string> & fargs, 
                   int & offset, 
                   FilterEngine *fe, 
                   FeatureEngine *fte
                   )
      : FT_WWAM<UCHAR>(fInd,
		       paramInd,
		       fargs, offset,
		       fe,
		       AA_ALPHABET_SIZE
		       ) {
      
      this->gram = (FL_BaseGram<TClass1,TClass2> *)fe->getFilter(fargs[offset++]);
      initialize();
      
    }
    
    void initialize() {
      assert(!(_offset % 3) && !(_length % 3) && !(_step % 3) && !(_windowSize % 3));
      aa = new AASigma();
      //      this->vals = new int[(_offset + _length + _windowSize) / 3];

      this->_order = (int)gram->order();
      this->_numPosVals = (windowSize() - _order)/3 + 1;
      this->_numUpdates = ((this->length() - 1)/this->step() + 1)*this->numPosVals();
      this->_numParams.s = ((this->length() - 1)/this->step() + 1)*this->domainSize();
      this->vals = new int[this->numUpdates()];

      assert(_order == 3);
      
    }
    
    
    /**
     * A member function that computes protected member vals. Any entry
     * appearing at some position within array vals is assigned the 
     * amino acid appearing at that position.
     * The positions in the sequence that are covered by vals are between
     * ge->getPos() + offset and ge->getPos() + offset + length
     * @return vals
     * @see int* FT_GramWWAM::computeValues(Tag *)
     */    
    inline int* computeValues(Tag *ge) {
      TClass2 *gramVals = gram->values(ge->getStrand()) + ge->getPos() + this->offset();
      for(register int i = 0; i < this->length(); i += this->step()) {
	int ni = i/this->step();
	int offset_v = ni*this->domainSize();
	int offset_i = ni*this->numPosVals();
	
	for(register int j = 0; j < this->numPosVals(); j++) {
	  assert(offset_i + j < this->numUpdates());
	  this->vals[offset_i + j] = offset_v + aa->AAIndex((int)gramVals[i+3*j]);
	  //	  cerr << "\t" << i << " " << j << " " << vals[offset_i + j] << endl;
	}
      }
      return this->vals;
    }
    /*
    inline UCHAR* computeValues(Tag *ge) {
      int numVals = offset() + length() + windowSize() - 1;
      int index = ge->getPos() - offset() - windowSize()/2;
      
      for(int j = 0; j < numVals; j += 3)
        this->vals[j/3] = (UCHAR)aa->AAIndex((int)(*_gramVals)[ge->getStrand()][j + index]);
      
      return this->vals + windowSize()/6;
      } */
    
    ~FT_PeptideWWAM() {
      delete aa;
    }
    
  };

}

#endif
