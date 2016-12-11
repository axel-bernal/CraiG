/****************************************************************************
* FL_ScoreFile.h - part of the lless namespace, a general purpose
*                     linear semi-markov structure prediction library
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

#ifndef _FL_SCORE_FILE_H_
#define _FL_SCORE_FILE_H_
  
#include "FL_FastaxFilter.h"

/**
 * A namespace created for classes belonging to lless, a general purpose linear
 * semi-markov structure prediction library.
 */
namespace lless {
  
  /**
   * Class to load scores from a file and access them directly as filter
   * values. It is derived from FL_InpFile
   */

  template <class TClass1, class TClass2, class TClass3, class TClass4> 
    class FL_ScoreFile : public FL_InpFile<TClass1, TClass3> {
   protected:
    int smoothWindow;
   public:
    //! Default constructor
    FL_ScoreFile(
                 int fInd,          //!< A unique identifier.
                 std::string name, //!< A unique name.
                 IOSource *file,     //!< A Score file
                 int maxFilterVals,  //!< number of filter values
                 int seqNo,          //!< sequence number inside file
                 TValType type = FT_DOUBLE, //!< Filter value type.
                 int period = 0,            //!< period for accumulation
		 int smoothWindow = 0 //!<smoothing window size
                 ) : FL_InpFile<TClass1, TClass3>(fInd, 
                                                  name,
                                                  file,
                                                  maxFilterVals,
                                                  seqNo,
                                                  type
                                                  ) {
      
      if(this->maxFilterVals > 1)
        assert(this->file->type() == EXTENDED_INTSCORE || 
	       this->file->type() == EXTENDED_MULTI_INTSCORE);

      assert(this->file->type() != FASTA &&
             this->file->type() != XFASTA &&
             this->file->type() != MULTIZ &&
             this->file->type() != EDGEANNOT_FASTA);
      
      this->_period = period;
      this->smoothWindow = smoothWindow;
    }
    
    /**
     * Constructor from Header string definition
     */ 
    FL_ScoreFile(
                 int fInd,       //!< A unique identifier.
                 /*! The Header string definition, loaded as a vector of 
                  * strings.
                  * The Header has the following form:\n\n
                  * Filter name ScoreFile arrInd file\n\n
                  * The meaning of these fields are as in the default
                  * constructor
                  * @see FL_FastaxFilter()
                  */
                 vector<std::string> & params, 
                 int & offset,       //!< The index for vector params.
                 ResourceEngine *re, /*!< A pointer to the ResourceEngine 
                                       object. */
                 FilterEngine *fe    /*!< A pointer to the FilterEngine
                                       object. */
                 ) : FL_InpFile<TClass1, TClass3>(fInd, 
                                                  params, 
                                                  offset, 
                                                  0,
                                                  re) {

      this->maxFilterVals = this->file->alphabetSize();
      if(this->maxFilterVals > 1)
        assert(this->file->type() == EXTENDED_INTSCORE || 
	       this->file->type() == EXTENDED_MULTI_INTSCORE);
      
      assert(this->file->type() != FASTA &&
             this->file->type() != XFASTA &&
             this->file->type() != MULTIZ &&
             this->file->type() != EDGEANNOT_FASTA);

      this->_period = 0;
      this->smoothWindow = 0;
      if(offset < params.size()) {
        if(!sscanf(params[offset++].c_str(), "%d", &this->_period))
          assert(0);
      }
      if(offset < params.size()) {
        if(!sscanf(params[offset++].c_str(), "%d", &this->smoothWindow))
          assert(0);
      }
    }

    /**
     * A member function that will copy each character from Sequence
     * obejct c into the filter values. In case of sparse sequences, 
     * the function first will make a dense clone and extract the
     * values from there.
     */
    int copySequence(Sequence *c, TStrand strand) {
      assert(!this->isSparse());

      if(this->arrayInd() < 0)
        return 0;

      if(!this->file->loaded())
        return 0;
      
      this->setDefaultStrand(strand);

      bool isSubSeq = false;
      TClass2 *fileSeq
        = (TClass2 *)this->findSeqInFile(c->id(), isSubSeq);
                                                  
      /*
       * Set the filter values variable to point to the fileSeq's sequence
       * of values
       */
      if(!fileSeq) {
        assert(0);
        throw EXCEPTION(BAD_USAGE, c->id() + string("not found in file ") + this->getName());
      }

      int zero = (fileSeq->isZeroBased() ? 0 : 1);
      int length = fileSeq->length();
      TClass4 *clone = (TClass4 *)fileSeq->cloneSeq(strand);

      this->filterVals[strand][0] = TClass3();
      for(int i = 0; i < length; i++)
	this->filterVals[strand][i + 1] = (TClass3)clone[i + zero];
      
      fileSeq->deleteClone(clone);

      if(isSubSeq)  delete fileSeq;      

      return length;
    }

    /**
     * A member function for computing filter values.
     * 
     * The value at each position of c is whatever value the private member 
     * filter computes at said position. 
     *
     * @param c a pointer to a Sequence object which has been read from the
     * fasta file.
     * @param strand the strand to use for accessing the c's string sequence.
     * @see FL_InpFile::computeVals(Sequence *, TStrand)
     */  
    inline void computeVals(Sequence *c, TStrand strand) {
      int length = this->copySequence(c, strand);
      // smooth out sequence
      if(smoothWindow)
	Utils::smoothScores(this->filterVals[strand], smoothWindow,
			    this->period(), length);
      
      if(this->isSparse())
	return;

      // log scale 
      if(this->isLogScaled())
	this->logScaleVals(this->filterVals[strand], length);
      // power to this->power()
      if(this->power() != 1)
	raiseVals2Pow(this->filterVals[strand], length, this->power());
      // accumulate
      if(this->period() && length)
	this->accumVals(this->filterVals[strand],
			this->accFilterVals[strand], 
			length);   
    }
  };

  /**
   * FL_NormalizedScoreFile class derives from class FL_ScoreFile. adds
   * normalization procedure before computing the values
   *
   */
  template <class TClass1, class TClass2, class TClass3, class TClass4> 
    class FL_NormalizedScoreFile : public 
    FL_ScoreFile<TClass1, TClass2, TClass3, TClass4> {
   public:
    //! Default constructor
    FL_NormalizedScoreFile(
                 int fInd,            //!< A unique identifier.
                 std::string &name,   //!< A unique name.
                 IOSource *file,      //!< A Score file
                 int maxFilterVals,   //!< number of filter values
                 int seqNo,           //!< sequence number inside file
                 TValType type = FT_DOUBLE, //!< Filter value type.
                 int period = 0,         //!< period for accumulation
		 int smoothWindow = 0 //!< smooth window size
			   ) :
      FL_ScoreFile<TClass1, TClass2, TClass3, TClass4>(fInd, 
						       name,
						       file,
						       maxFilterVals,
						       seqNo,
						       type,
						       period,
						       smoothWindow
						       ) {
	
    }
    
    /**
     * Constructor from Header string definition
     */ 
    FL_NormalizedScoreFile(
                 int fInd,       //!< A unique identifier.
                 /*! The Header string definition, loaded as a vector of 
                  * strings.
                  * The Header has the following form:\n\n
                  * Filter name NormalizedScoreFile arrInd file\n\n
                  * The meaning of these fields are as in the default
                  * constructor
                  * @see FL_FastaxFilter()
                  */
                 vector<std::string> & params, 
                 int & offset,       //!< The index for vector params.
                 ResourceEngine *re, /*!< A pointer to the ResourceEngine 
                                       object. */
                 FilterEngine *fe    /*!< A pointer to the FilterEngine
                                       object. */

                 ) : 
      FL_ScoreFile<TClass1, TClass2, TClass3, TClass4>(fInd, 
						       params, 
						       offset, 
						       re,
						       fe) {
      
    }

    /**
     * A member function for computing filter values.
     * 
     * The value at each position of c is whatever value the private member 
     * filter computes at said position. 
     *
     * @param c a pointer to a Sequence object which has been read from the
     * fasta file.
     * @param strand the strand to use for accessing the c's string sequence.
     * @see FL_ScoreFile::computeVals(Sequence *, TStrand)
     */  
    inline void computeVals(Sequence *c, TStrand strand) {
      int length = this->copySequence(c, strand);
      // log scale it
      if(this->isLogScaled())
	this->logScaleVals(this->filterVals[strand], length);

      // smooth out sequence
      if(this->smoothWindow)
	Utils::smoothScores(this->filterVals[strand], this->smoothWindow,
			    this->period(), length);
      // normalize sequence
      TClass3 maxValue = TClass3();
      for(int i = 1; i <= length; i++)
	if(this->filterVals[strand][i] > maxValue)
	  maxValue = this->filterVals[strand][i];
      
      if(maxValue)
	for(int i = 1; i <= length; i++)
	  this->filterVals[strand][i] /= maxValue;
     
      // accumulate
      if(this->period() && !this->isSparse() && length)
	this->accumVals(this->filterVals[strand],
			this->accFilterVals[strand], 
			length);
    }
  };
}

#endif
