/****************************************************************************
* BioFeature.h - part of the lless namespace, a general purpose
*                linear semi-markov structure prediction library
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

#ifndef _BIO_FEATURE_H_
#define _BIO_FEATURE_H_
#include "Sequence.h"


namespace lless {

  /**
   * The BioFeature class is the base class for bio-type features such as 
   * Genes, Transcripts and so on. The most significant characteristic of this
   * class is the fact that the internal representation can be stated in 
   * OneStrand or TwoStrand format which are explained elsewhere in detail.
   ***************************************************************************/

  class BioFeature {
   protected:
    Sequence *annotSeq;
    string id;
    int _begin, _end;
    TStrand strand;
    bool inOneStrand; //!< either OneStrand if true or TwoStrand otherwise

   public:
    /**
     * Constructor
     * @param iden a unique identifier for this biofeature
     * @param annotSeq a pointer to the Sequence object.
     * @param inOneStrand biofeature in OneStrand format if true or in
     * TwoStrand format otherwise.
     * @param strand the biofeature's strand.
     */
    BioFeature(std::string &iden, 
               Sequence *annotSeq, 
               bool inOneStrand, 
               TStrand strand) {

      id = iden;
      this->annotSeq = annotSeq;
      this->_begin = INT_MAX;
      this->_end = -1;
      this->strand = strand;
      this->inOneStrand = inOneStrand;
    }

    BioFeature(std::string &iden, 
               Sequence *annotSeq, 
	       int begin, int end,
               bool inOneStrand, 
               TStrand strand) {

      id = iden;
      this->annotSeq = annotSeq;
      this->_begin = begin;
      this->_end = end;
      this->strand = strand;
      this->inOneStrand = inOneStrand;
    }


    /**
     * Default constructor, no parameters. annotSeq is made equal to NULL.
     */
    BioFeature() { annotSeq = NULL; }

    virtual void toOneStrand() {;}
    virtual void toTwoStrand() {;}

    inline BioFeature & operator=(const BioFeature &biofeature) {
      BioFeature &bf = (BioFeature &)biofeature;
      id = bf.getId();
      annotSeq = bf.getSequence();
      _begin = bf.begin();      
      _end = bf.end();
      strand = bf.getStrand();
      inOneStrand = bf.isInOneStrand();
      return *this;
    }
    
    inline bool operator<(const BioFeature &biofeat) {
      BioFeature &bf = (BioFeature &)biofeat;

      if(!inOneStrand)
	throw EXCEPTION(BAD_USAGE, " biofeatures must be inOneStrand for sorting");
      if(this->strand != bf.getStrand()) 
        return this->strand < bf.getStrand();
     
      return (this->begin() < bf.begin() || 
	      (this->begin() == bf.begin() && this->end() < bf.end()));
    }

    inline bool operator==(const BioFeature &biofeature) {
      BioFeature &bf = (BioFeature &)biofeature;
      return id == bf.getId() &&
        annotSeq == bf.getSequence() &&
        strand == bf.getStrand();
    }
    
    inline Sequence * getSequence() {  return annotSeq; }
    inline string & getSequenceId() {  return annotSeq->id(); }
    inline std::string & getId() {  return id; }
    inline int begin() {  return _begin; }
    inline int end() {  return _end; }
    inline bool isInOneStrand() {  return inOneStrand; }
    inline TStrand getStrand() {  return strand; }

    inline void setSequence(Sequence *annotSeq) {  this->annotSeq = annotSeq; }
    inline void setStrand(TStrand strand) {  this->strand = strand;}
    inline void setBegin(int begin) {  _begin = begin; }
    inline void setEnd(int end) {  _end = end; }

    virtual ~BioFeature() {

    }
  };

}

namespace std {
  /**                                                                            
   * Routine for sorting STL containers instantiated with BioFeature* types
   */
  template <>
    struct greater<lless::BioFeature *> {
      bool operator()(const lless::BioFeature* b1, const lless::BioFeature* b2)
      {
        //Defined for ascending sorting
        if(!b1)
          return true;
        if(!b2)
          return false;
        return ((lless::BioFeature &)(*b1)) < (*b2);
      }
  };
}

#endif
