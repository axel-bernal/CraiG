/****************************************************************************
* FL_DiffCountingScore.h - part of the lless namespace, a general purpose
*                          linear semi-markov structure prediction library
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

#ifndef _FILTER_DIFFCOUNTINGSCORE_H_
#define _FILTER_DIFFCOUNTINGSCORE_H_
  
#include "Utils.h"
#include "FL_Score.h"
#include "ContextIMM.h"
#include "FilterEngine.h"
#include "ResourceEngine.h"
#include "IMM.h"
#include "Motif.h"
  
namespace lless {
  
  /**
   *
   * FL_DiffCountingScore is a subclass of FL_Score whose scoring function
   * depends on the difference of countings of two words which are part of 
   * gram given as parameter.
   */
  class FL_DiffCountingScore : public FL_Score<int> {
  
   private:
    FL_Gram<UCHAR> *gram; // gram to use as the source of countings
    UCHAR indPos; // count gram values equal to indPos
    UCHAR indNeg; // and substract gram values equal to indNeg 
    
   public:
    //! Default constructor.
    FL_DiffCountingScore(int fInd,             //!< A unique identifier.
                         std::string & name,   //!< A unique name.
                         int period,           //!< The filter values period
                         FL_Gram<UCHAR> *gram,           //!< A pointer to a gram object
                         UCHAR indPos,           /*!< gram index of word to 
                                                 count as positive */
                         UCHAR indNeg            /*!< gram index of word to 
                                                 count as negative */
                         )
      :  FL_Score<int>(fInd,
                       name, 
                       period,
                       0, FT_INTEGER
                       ) {
      
      this->gram = gram;
      this->indPos = indPos;
      this->indNeg = indNeg;
    }

    /**
     * Constructor from a Header string definition
     */  
    FL_DiffCountingScore(int fInd,        //!< A unique identifier.
                         /*! The Header string definition, loaded as a vector 
                          * of strings.
                          * The Header has the following form:\n\n
                          * Filter name DiffCountingScore period smoothenWSize 
                          * gram indPos indNeg\n\n
                          * The description of the fields could be found
                          * in the other constructor(s) 
                          */
                         vector<std::string> & params, 
                         int & offset,       //!< The index for vector params.
                         ResourceEngine *re, /*!< A pointer to the 
                                               ResourceEngine object. */
                         FilterEngine *fe    /*!< A pointer to the 
                                               FilterEngine object. */
                         )
      : FL_Score<int>(fInd, 
                      params,
                      offset, 
                      fe,
                      FT_INTEGER) {
      
      this->gram  = (FL_Gram<UCHAR> *)fe->getFilter(params[offset++]);
      assert(gram);
      
      if(!sscanf(params[offset++].c_str(), "%d", &indPos)) 
        assert(0);
      
      if(!sscanf(params[offset++].c_str(), "%d", &indNeg)) 
        assert(0);
      
    }

    /**
     * A member function that extends the functionality of the original
     * function definedin FL_Score to change the contained gram object's
     * default strand to strand
     * @param strand the new default strand
     */
    inline void setDefaultStrand(TStrand strand) {
      this->defaultStrand = strand;
      gram->setDefaultStrand(strand);
    }

    /**
     * The implemented scoring function. At each position p the score is 1(-1)
     * if sigma->ind(seq[p..p+order]) matches indPos(indNeg) or zero otherwise.
     * @param seq the input sequence.
     * @param filterVals the array of filter values.
     * @param len the length to consider. It could be shorter than the
     * actual length of the input sequence seq.
     */  
    void scoreInputSeq(char *seq, int *filterVals, int len) {
      //      cerr << this->getName() << " " << ind << "\n";
      for(int i = 0; i < len; i++) {
        UCHAR v = gram->value(i + 1);
        filterVals[i + 1] = ((v == indPos) ? 1 : ((v == indNeg) ? -1 : 0));
        //        cerr << seq[i] << " " << filterVals[i + 1] << endl;
      }
    }
  
    ~FL_DiffCountingScore() {
  
    }

  };

  /**
   * FL_RelPeptideScore
   */
  class FL_RelPeptideScore : public FL_Score<double> {
   protected:
    FL_Score<double> *scores;
    int maxRange;

   public:
    FL_RelPeptideScore(int fInd,              //!< A unique identifier.
		       std::string & name,    //!< A unique name.
		       FL_Score<double> *scores,
		       int maxRange)
      :  FL_Score<double>(fInd,
			  name, 
			  0,
			  0, FT_DOUBLE) {

      this->scores = scores;
      this->maxRange = maxRange;

    }

    /**
     * Constructor from a Header string definition
     */
    FL_RelPeptideScore(int fInd,        //!< A unique identifier. 
		       /*! The Header string definition, loaded as a vector 
			* of strings.
			* The Header has the following form:\n\n
			* Filter name SpacerMotifScore arrInd period 
			* smoothenWSize threshold1 motif1 threshold2 motif2
			* minRange maxRange \n\n
			* The description of the fields could be found in 
			* the other constructor(s)
			*/
		       vector<std::string> & params, 
		       int & offset,       //!< The index for vector params.
		       ResourceEngine *re, /*!< A pointer to the 
					     ResourceEngine object. */
		       FilterEngine *fe    /*!< A pointer to the FilterEngine 
					     object. */
		       )
      : FL_Score<double>(fInd, 
			 params,
			 offset, 
			 fe,
			 FT_DOUBLE) {
      
      scores = (FL_Score<double> *)fe->getFilter(params[offset++]);
      if(!sscanf(params[offset++].c_str(), "%d", &maxRange))
	assert(0);
    }
      
    inline void setDefaultStrand(TStrand strand) {
      this->defaultStrand = strand;
      scores->setDefaultStrand(strand);
    }

    /**
     * The implemented scoring function. At each position p the score is 1(-1)
     * if sigma->ind(seq[p..p+order]) matches indPos(indNeg) or zero otherwise.
     * @param seq the input sequence.
     * @param filterVals the array of filter values.
     * @param len the length to consider. It could be shorter than the
     * actual length of the input sequence seq.
     */  
    void scoreInputSeq(char *seq, double *filterVals, int len) {
      //      cerr << this->getName() << "\n";
      vector< pair<int, double> > trail(3), lead(3);
      int i;

      for(i = 1; i <= len; i++) {
        double v = scores->value(i);

	if(trail[i%3].first && i - trail[i%3].first > maxRange) 
	  trail[i%3].second = 0;

	filterVals[i] = trail[i%3].second;
	  
	if(v)
	  trail[i%3] = pair<int, double>(i, v);
      }

      for(i = len; i >= 1; i--) {
        double v = scores->value(i);

	if(lead[i%3].first && lead[i%3].first - i > maxRange) 
	  lead[i%3].second = 0;

	if(lead[i%3].second > filterVals[i])
	  filterVals[i] = lead[i%3].second;
	
	if(v)
	  lead[i%3] = pair<int, double>(i, v);

      }

      /*      for(i = 1; i <= len; i++) 
      	cerr << i << "|" << scores->value(i) << "|" << filterVals[i] << " ";
	cerr << endl; */
    }

    ~FL_RelPeptideScore() {
    }
    
  };


  /**
   * FL_RelPeptideOffset
   */
  class FL_RelPeptideOffset : public FL_Score<double> {
   protected:
    FL_Score<double> *scores;
    int maxRange;

   public:
    FL_RelPeptideOffset(int fInd,              //!< A unique identifier.
			std::string & name,    //!< A unique name.
			FL_Score<double> *scores,
			int maxRange)
      :  FL_Score<double>(fInd,
			  name, 
			  0,
			  0, FT_DOUBLE) {

      this->scores = scores;
      this->maxRange = maxRange;

    }

    /**
     * Constructor from a Header string definition
     */
    FL_RelPeptideOffset(int fInd,        //!< A unique identifier. 
			/*! The Header string definition, loaded as a vector 
			 * of strings.
			 * The Header has the following form:\n\n
			 * Filter name SpacerMotifScore arrInd period 
			 * smoothenWSize threshold1 motif1 threshold2 motif2
			 * minRange maxRange \n\n
			 * The description of the fields could be found in 
			 * the other constructor(s)
			 */
			vector<std::string> & params, 
			int & offset,       //!< The index for vector params.
			ResourceEngine *re, /*!< A pointer to the 
					      ResourceEngine object. */
			FilterEngine *fe    /*!< A pointer to the FilterEngine 
					      object. */
			)
      : FL_Score<double>(fInd, 
			 params,
			 offset, 
			 fe,
			 FT_DOUBLE) {
      
      scores = (FL_Score<double> *)fe->getFilter(params[offset++]);
      if(!sscanf(params[offset++].c_str(), "%d", &maxRange))
	assert(0);
    }
      
    inline void setDefaultStrand(TStrand strand) {
      this->defaultStrand = strand;
      scores->setDefaultStrand(strand);
    }

    /**
     * The implemented scoring function. At each position p the score is 1(-1)
     * if sigma->ind(seq[p..p+order]) matches indPos(indNeg) or zero otherwise.
     * @param seq the input sequence.
     * @param filterVals the array of filter values.
     * @param len the length to consider. It could be shorter than the
     * actual length of the input sequence seq.
     */  
    void scoreInputSeq(char *seq, double *filterVals, int len) {
      //      cerr << this->getName() << "\n";
      vector< pair<int, double> > trail(3), lead(3);
      vector<double> auxVals(len + 1);
      int i;

      for(i = 1; i <= len; i++) {
        double v = scores->value(i);
	
	if(trail[i%3].first)
	  if(i - trail[i%3].first > maxRange)
	    trail[i%3].second = 0;
	  else
	    filterVals[i] = i - trail[i%3].first;

	auxVals[i] = trail[i%3].second;

	if(v)
	  trail[i%3] = pair<int, double>(i, v);
      }

      for(i = len; i >= 1; i--) {
        double v = scores->value(i);

	if(lead[i%3].first && lead[i%3].first - i > maxRange) 
	  lead[i%3].second = 0;

	if(lead[i%3].second > auxVals[i]) {
	  auxVals[i] = lead[i%3].second;
	  filterVals[i] = lead[i%3].first - i;
	}
 
	if(v)
	  lead[i%3] = pair<int, double>(i, v);

      }

      /*      for(i = 1; i <= len; i++) 
      	cerr << i << "|" << scores->value(i) << "|" << filterVals[i] << " ";
      cerr << endl;
      */
    }

    ~FL_RelPeptideOffset() {
    }
    
  };

}

#endif
