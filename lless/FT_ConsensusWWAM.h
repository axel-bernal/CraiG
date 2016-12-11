/****************************************************************************
* FT_ConsensusWWAM.h - part of the lless namespace, a general purpose
*                      linear semi-markov structure prediction library
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

#ifndef FT_CONSENSUSWWAM_FEAT_H
#define FT_CONSENSUSWWAM_FEAT_H
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "FT_WWAM.h"
#include "NGram.h"
#include "Consensus.h"
#include "FT_Edge.h"


namespace lless {
    
  /**
   * FT_ConsensusWWAM is a subtype of FT_WWAM that uses indicator variables
   * derived from a Consensus sequence given as parameter. The indicator
   * variables are defined at each position of the signal, so the consensus
   * length must match the signal length.
   *
   **************************************************************************/
  template <class TClass>
  class FT_ConsensusWWAM : public FT_WWAM<UCHAR> {
    
   protected:
    FL_Gram<TClass> *gram;
    Consensus *consensus;
    
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
     * @param consensus consensus pattern.
     * @param offset the length of the WWAM, offset of the signal occurrence.
     * @param length  the length of the WWAM, length of the signal occurrence. 
     * The total length of the WWAM would be then offset + length.
     * @param step the number of positions to jump within the WWAM, before
     * counting again.
     */
    FT_ConsensusWWAM(
                     int fInd, 
                     int paramInd, 
                     int parsingFrames,
                     char *name,
                     FilterEngine *fe,
                     FL_Signal<EdgeInst> *signal,
                     int windowSize, 
		     FL_Gram<TClass> *gram,
                     char *consensus,
                     int offset, 
                     int length,
                     int step
                     ) 
      : FT_WWAM<UCHAR>(fInd, 
		       paramInd, 
		       parsingFrames,
		       name,fe,
		       signal, 
		       FT_UCHAR,
		       windowSize, 2, 
		       offset, length,
		       step) {
      
      this->gram = gram;
      initialize(consensus);

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
     * Feature name ConsensusWWAM parsingFrames signal offset length step 
     * windowSize gram consensus \n\n
     * The description of the fields could be found in 
     * the other constructor(s)
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param fte A pointer to a FeatureEngine object.
     */
    
    FT_ConsensusWWAM(
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
			2) {
      
      this->gram = (FL_Gram<TClass> *)fe->getFilter(fargs[offset++]);
      initialize((char *)(fargs[offset++].c_str()));

    }

    /**
     * Constructor from a Header string definition. windowSize should be passed
     * as parameter
     * @param fInd A feature index.     
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * object to which this feature is tied to.
     * @param fargs The Header string definition loaded as a vector of strings.
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param windowSize the window size
     */
    FT_ConsensusWWAM(
                     int fInd, 
                     int paramInd, 
                     vector<std::string> & fargs, 
                     int & offset, 
                     FilterEngine *fe, 
                     int windowSize
                     )
      : FT_WWAM<UCHAR>(fInd,
			paramInd,
			fargs, offset,
			fe, windowSize,
			2) {
      
      this->gram = (FL_Gram<TClass> *)fe->getFilter(fargs[offset++]);      
      initialize((char *)(fargs[offset++].c_str()));

    }
    
    inline void initialize(char *cons) {
      this->_order = (int)gram->order();
      /*
       * consensus must have the sides of the signal window covered up to
       * _windowSize/2, each side
       */
      assert(_order == 1);

      this->_numPosVals = this->windowSize() - this->order() + 1;
      //      this->vals = new UCHAR[_offset + _length + _windowSize];
      this->_numUpdates = ((this->length() - 1)/this->step() + 1)*this->numPosVals();
      this->_numParams.s = ((this->length() - 1)/this->step() + 1)*this->domainSize();
      this->vals = new int[this->numUpdates()];

      char *final_cons = interpretConsensus(cons);
      this->consensus = new Consensus(final_cons, gram->alphabet(),
				      gram->order(), 
				      gram->maxNumFilterValues());
      delete [] final_cons;
      
    }

    char *interpretConsensus (char *consensus) {
      char *cons = new char[this->numUpdates() + 1];
      
      if(consensus[1] == '*') {
        int num;
        
        if(!sscanf(consensus + 2, "%d", &num))
          assert(0);
        assert(num == this->numUpdates());
        for(int i = 0; i < num; i++)
          cons[i] = consensus[0];
        
        cons[num] = '\0';
      }
      else  {
	assert(strlen(consensus) == this->numUpdates());
        strncpy(cons, consensus, this->numUpdates());
	cons[this->numUpdates()] = 0;
      }
      
      return cons;
    }    

    /**
     * A member function that computes protected member vals. Any entry
     * appearing at some position within array vals is assigned 1 if the
     * gram at that position matches the consensus sequence or 0 otherwise.
     * The positions in the sequence that are covered by vals are between
     * ge->getPos() + offset and ge->getPos() + offset + length
     * @return vals
     * @see int* FT_GramWWAM::computeValues(Tag *)
     */
    inline int* computeValues(Tag *ge) {
      TClass *gramVals = gram->values(ge->getStrand()) + ge->getPos() + this->offset();
      for(register int i = 0; i < this->length(); i += this->step()) {
	int ni = i/this->step();
	int offset_v = ni*this->domainSize();
	int offset_i = ni*this->numPosVals();
	
	for(register int j = 0; j < this->numPosVals(); j++)
	  //          cerr << "\t" << i << " " << j << " " << vals[j] << " " << offset + vals[j] << endl;
	  this->vals[offset_i + j] = offset_v + consensus->matches((int)gramVals[i+j], i+j, _order);
      }
      return this->vals;
    }

    /*    inline UCHAR* computeValues(Tag *ge) {
      
      int numVals = offset() + length() + windowSize() - 1;
      int index = ge->getPos() - offset() - windowSize()/2;
      
      for(int j = 0; j < numVals; j++) {
        TClass nt = (*_gramVals)[ge->getStrand()][j + index];
        this->vals[j] = consensus->matches((int)nt, j, _order);
      }
      
      return this->vals + windowSize()/2;
      }*/


    ~FT_ConsensusWWAM() {
      delete consensus;
    }
    
  };
  
}

#endif
