/****************************************************************************
*   FL_Reranker.h - part of the lless namespace, a general purpose
*                   linear semi-markov structure prediction library
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

#ifndef _SEQTAG_SCORE_H_
#define _SEQTAG_SCORE_H_
  
#include "Utils.h"
#include "Filter.h"
#include "ContextIMM.h"
#include "FilterEngine.h"
#include "ResourceEngine.h"
#include "SequenceUtils.h"
#include "TagUtils.h"
#include "IMM.h"
#include "Motif.h"
  
namespace lless {

  /**
   *
   * FL_Reranker is a base class for filter that compute filter values 
   * to be used in a reranking setting
   *
   ***************************************************************************/

  template <class TClass> class FL_Reranker : public TypedFilter<TClass> {
   protected:
    TSetType setType;
    TTagClass tagClass;
    int tagType;
    TypedFilter<TClass> *filter;

   public:
    //! Default constructor
    FL_Reranker(
                int fInd,       //!< A unique identifier.
                std::string & name,    //!< A unique name.
                TSetType setType,       //!< The SeqTags are under this setType
                TTagClass tagClass,
                int tagType,
                TypedFilter<TClass> *filter, //!< Contained filter
                TValType type = FT_DOUBLE //!< Filter value type.
                ) : 
    TypedFilter<TClass>(fInd, 
			name,
			1, type) {
      
      this->setType = setType;
      this->tagClass = tagClass;
      this->tagType = tagType;
      this->filter = filter;

    }
    
    /**
     * First constructor from a Header string definition
     */    
    FL_Reranker(
                int fInd,         //!< A unique identifier.
                vector<std::string> &params, //!< Header string definition
                int & offset,       //!< The index for vector params.
                FilterEngine *fe    //!< A pointer to the FilterEngine object.
                ) 
      : TypedFilter<TClass>(fInd, 
                            params, 
                            offset, 
                            1) {
      
      this->setType = SequenceUtils::stringToTSetType(params[offset++]);
      this->tagClass = TagUtils::stringToTTagClass(params[offset++]);

      if(this->tagClass == EDGE_INST) 
	this->tagType = TypeDefs::stringToTEdgeId2(params[offset++]);
      else 
	this->tagType = TypeDefs::stringToTNodeId3(params[offset++]);

      this->filter = (TypedFilter<TClass> *)fe->getFilter(params[offset++]);
      
    }
    
    int findLowestRank(vector<SeqTags> & seqTags) {

      int i = 0;
      int index = -1;
      int minRank = INT_MAX;
      for(i = 0; i < seqTags.size(); i++)
        if(seqTags[i].rank() < minRank) {
          index = i;
          minRank = seqTags[i].rank();
        }

      assert(index >= 0);
      return index;
    }

    Tag * findTargetTag(SeqTags &stObj) {
      double bstScore = -DOUBLE_INFINITY, score;
      SeqTags::iterator it = stObj.begin();
      Tag *tag = NULL;
      for( ; it != stObj.end(); it++) {
	Tag *t = (*it);
        if(t->getGEClass() != tagClass ||
           TagUtils::tagType(t) != tagType ||
	   t->getStrand() != this->defaultStrand)
          continue;
        score = filter->value(t->getPos(), t->getStrand());

        if(score > bstScore) {
          bstScore = score;
          tag = t;
        }
      }
      
      return tag;
    }
    
    virtual void computeRerankVals(Tag *tag, vector<SeqTags> &seqTags,
                                   int lrank_ind, TStrand strand) = 0;
    
    
    /**
     * A member function for computing filter values.
     * 
     */
    inline void computeVals(Sequence *c, TStrand strand) {

      assert(!this->isSparse());

      if(this->arrayInd() < 0)
        return;

      this->setDefaultStrand(strand);
      //      cerr << c->id() << " " << this->getName() << " " << strand\n";

      vector<SeqTags> & seqTags = c->getTags(setType);
      int lrank_ind = seqTags[0].rank();
      Tag *targetTag = findTargetTag(seqTags[seqTags[0].rank()]);

      computeRerankVals(targetTag, seqTags, lrank_ind, strand);
      
    }
    
    virtual ~FL_Reranker() { }
    
  };
  
  /**
   *
   * FL_RerankScore computes the difference in score between the top 
   * scoring SeqTag object and any other and assigns it to the position
   * at which the start of the SeqTag objects begins. The position of 
   * the Top scoring SeqTag in the list is assumed to be at the top
   * This filter works in reranking settings
   * /todo make it work for topK instead of top1.
   *
   ***************************************************************************/
  
  template <class TClass> class FL_RerankScore : public FL_Reranker<TClass> {
   public:
    //! Default constructor
    FL_RerankScore(
                    int fInd,       //!< A unique identifier.
                    std::string & name,   //!< A unique name.
                    TSetType setType,      //!< The SeqTags are under setType
                    TTagClass tagClass,
                    int tagType,
                    TypedFilter<TClass> *filter, //!< Contained filter
                    TValType type = FT_DOUBLE    //!< Filter value type.
		   ) : 
    FL_Reranker<TClass>(fInd, 
			name,
			setType,
			tagClass,
			tagType,
			filter,
			type) {
      
    }
    
    /**
     * First constructor from a Header string definition
     */    
    FL_RerankScore(
		   int fInd,         //!< A unique identifier.
		   vector<std::string> &params, //!< Header string definition
		   int & offset,       //!< The index for vector params.
		   ResourceEngine *re,
		   FilterEngine *fe    //!< A pointer to the FilterEngine obj.
		   ) 
      : FL_Reranker<TClass>(fInd, 
			    params, 
			    offset, 
			    fe) {
      
    }
    
    /**
     * A member function for computing filter values using reranking
     * information.
     * 
     */
    
    inline void computeRerankVals(Tag *tag, vector<SeqTags> &seqTags,
                                  int lrank_ind, TStrand strand) {
      
      TClass lrank_score = TClass();
      if(tag)
        lrank_score = this->filter->value(tag->getPos(), tag->getStrand());
      
      for(int i = 0; i < seqTags.size(); i++) {
        if(i == lrank_ind)
          continue;
	
        SeqTags & stObj = seqTags[i];
        SeqTags::iterator it = stObj.begin();
        
        for( ; it != stObj.end(); it++) {
          Tag *t = (*it);
          if(t->getGEClass() != this->tagClass ||
	     TagUtils::tagType(t) != this->tagType ||
             t->getStrand() != strand)
            continue;

          TClass score = this->filter->value(t->getPos(), strand);
          this->filterVals[strand][t->getPos()] = score - lrank_score;

        }
      }
    }
    
    ~FL_RerankScore() { }
  
  };

  /**
   *
   * FL_RerankOffset computes the offset between the top scoring SeqTag
   * and any other and assigns it to the position at which the SeqTag
   * object begins. The position of the Top scoring SeqTag in the list
   * is assumed to be at the top. 
   * This filter works in reranking settings
   * /todo make it work for topK instead of top1.
   *
   ***************************************************************************/

  template <class TClass> class FL_RerankOffset : public FL_Reranker<TClass> {
   public:
    //! Default constructor
    FL_RerankOffset(
                    int fInd,       //!< A unique identifier.
                    std::string & name,  //!< A unique name.
                    TSetType setType,     //!< set label the SeqTags are under 
                    TTagClass tagClass,
                    int tagType,
                    TypedFilter<TClass> *filter, //!< Contained filter
                    TValType type = FT_DOUBLE //!< Filter value type.
                    ) : 
      FL_Reranker<TClass>(fInd, 
                          name,
                          setType,
                          tagClass,
                          tagType,
                          filter,
                          type) {

    }

    /**
     * First constructor from a Header string definition
     */    
    FL_RerankOffset(
                    int fInd,         //!< A unique identifier.
                    vector<std::string> &params, //!< Header string definition
                    int & offset,       //!< The index for vector params.
		    ResourceEngine *re,
                    FilterEngine *fe    //!< A pointer to the FilterEngine obj.
                    ) 
      : FL_Reranker<TClass>(fInd, 
                            params, 
                            offset, 
                            fe) {
  
    }

    /**
     * A member function for computing filter values.
     * 
     */
    inline void computeRerankVals(Tag *tag, vector<SeqTags> &seqTags,
                                  int lrank_ind, TStrand strand) {

      for(int i = 0; i < seqTags.size(); i++) {
        if(i == lrank_ind)
          continue;
        SeqTags & stObj = seqTags[i];
        SeqTags::iterator it = stObj.begin();
        
        for( ; it != stObj.end(); it++) {
          Tag *t = (*it);
          if(t->getGEClass() != this->tagClass ||
	     TagUtils::tagType(t) != this->tagType ||
             t->getStrand() != strand)
            continue;

          int offset = 0;
          if(tag)
            offset = abs(t->getPos() - tag->getPos());
	  
          this->filterVals[strand][t->getPos()] = offset;
          
        }
      }
    }

    ~FL_RerankOffset() { }

  };
  
}

#endif
