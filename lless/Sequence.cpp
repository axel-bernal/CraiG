#include <algorithm>
#include "Utils.h"
#include <ctype.h>
#include "IMM.h"
#include "Gene.h"
#include "Sequence.h"
#include <sstream>

/****************************************************************************
* Sequence.cpp - part of the lless namespace, a general purpose
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


using namespace std;

namespace lless {

  list<Tag *> SeqTags::_allTags;

  void SeqTags::reverseComplement(FSM &fsm, int numPhases, int seqLength) {
    int swapPhase = _initPhase;
    _initPhase = (numPhases - _endPhase) % numPhases;
    _endPhase = (numPhases - swapPhase) % numPhases;
    this->reverse();
    SeqTags::iterator it = this->begin();

    while(it != this->end()) {
      Tag &tag = *(*it);

      std::map<Tag *, int, std::less<Tag *> >::iterator bit = phaseBreaks.find(*it);
      
      int begin = tag.getPos();
      tag.setStrand(tag.getStrand() == STRAND_FWD ? STRAND_COMP : STRAND_FWD);
      
      if(tag.getGEClass() == NODE_INST) {
        int end = begin + ((NodeInst &)tag).getLen() - 1;       
        TParseNode node = (TParseNode)tag.getParseType();
        node = fsm.node(node)->complementNode()->id();
        tag.setParseType(node);
        tag.setPos(seqLength - end + 1);
      }
      else {
        TParseEdge edge = (TParseEdge)tag.getParseType();
        edge = fsm.edge(edge)->complementEdge()->id();
        tag.setParseType(edge);
        tag.setPos(seqLength - begin + 1);
      }

      if(bit != phaseBreaks.end()) {
        int pBreak = bit->second;
        phaseBreaks.erase(bit);
        /* 
         * bug here as the inserted phase break can be different from zero
         * but I cannot recover it as this point
         */
        insertPhaseBreak(*it, 0); 
      }
      it++;
    }
  }     

  void SeqTags::releaseMem() {
    list<Tag *>::iterator it = _allTags.begin(); 
    for( ; it != _allTags.end(); it++) {
      if(*it)
        delete *it;
      *it = NULL;
    }
    _allTags.clear();
  }

  BasicSeq::BasicSeq(std::string &cId, bool circ, bool zero_b) {
    seqId = cId;
    circular = circ;
    zero_based = zero_b;
    seqLen[STRAND_FWD] = seqLen[STRAND_COMP] = 0;
  }

  /**
   * Copy constructor
   */
  BasicSeq::BasicSeq(const BasicSeq & seq) {
    *this = seq;
  }

  BasicSeq& BasicSeq::operator= (const BasicSeq & seq) {
    BasicSeq & s = (BasicSeq &)seq;
    seqId = s.id();
    circular = s.isCircular();
    zero_based = s.isZeroBased();

    for(int i = 0; i < NUM_STRANDS; i++)
      seqLen[i] = s.length((TStrand)i);
  }


  BasicSeq* BasicSeq::findSequence(std::string & id, 
                                   list<BasicSeq *> & seqs) {

    for(list<BasicSeq *>::iterator cit = seqs.begin(); cit != seqs.end(); cit++) {
      BasicSeq *c = *cit;
      if(!c->id().compare(id)) {
        return c;
        break;
      }
    }
    return NULL;
  }


  /**
   * Main constructor.
   * @param cId unique sequence identifier.
   * @param cSeq the sequence as a stream of TClass objects.
   * @param the alphabet that recognizes all symbols present in cSeq.
   * @param circ true if the sequence is circular, false otherwise.
   */
  Sequence::Sequence(std::string & cId, 
                     char *cSeq,
                     Sigma *alphabet, 
                     bool circ) : BasicSeq(cId, circ) {

    this->sigma = alphabet;

    for(int strand = 0; strand < NUM_STRANDS; strand++)
      seq[strand] = "";
    
    if(cSeq && strlen(cSeq))
      setSeq(cSeq);
  }

  Sequence& Sequence::operator= (const Sequence & seq) {
    int strand;
    (BasicSeq &)*this = (BasicSeq &)seq;
    Sequence & s = (Sequence &)seq;
    sigma = s.alphabet();    

    for(strand = 0; strand < NUM_STRANDS; strand++) {

      if(s.getSeq((TStrand)strand))
        this->seq[strand] = s.getSeq((TStrand)strand);
      else
        this->seq[strand] = "";

    }

    return *this;
  }

  /**                                                                           
   * Sets the sequence on strand strand to cSeq. 
   */
  void Sequence::setSeq(char *cSeq, TStrand strand) {
    if(!cSeq)
      return;
    
    this->seqLen[strand] = strlen(cSeq);
    seq[strand] = cSeq;
    
    if(this->seqLen[STRAND_FWD] && this->seqLen[STRAND_COMP] &&
       this->seqLen[STRAND_FWD] != this->seqLen[STRAND_COMP]) {
      assert(0);
      throw EXCEPTION( MISSING_ARGUMENT, string("different seq lengths"));
    }
  }
  

  /**
   * Complements the sequence. If typeSet is different from NO_SET it also
   * complements the Tag objects contained in the sequences
   */
  void Sequence::reverseComplement(FSM &fsm,
                                   int numPhases, 
                                   TSetType typeSet) {
    
    std::string swap = seq[STRAND_FWD];
    seq[STRAND_FWD] = seq[STRAND_COMP];
    seq[STRAND_COMP] = swap;
    
    if(typeSet != NO_SET && tags[typeSet].size()) {
      for(int i = 0; i < tags[typeSet].size(); i++) 
        tags[typeSet][i].reverseComplement(fsm, numPhases, length());
    }
  }

  /**
   * Computes the complementary sequence using the alphabet's complement
   * function applied on the forward sequence, i.e. the natural
   * complement.
   */
  void Sequence::setNatComplementSeq() {
    if(seq[STRAND_FWD].length() != 0) {
      
      for(int i = 0; i < seqLen[STRAND_FWD]; i++)
        seq[STRAND_COMP].push_back(sigma->toComplement(seq[STRAND_FWD][i]));
      
      seqLen[STRAND_COMP] = seqLen[STRAND_FWD];
      Utils::reverse(seq[STRAND_COMP]);
    }
    else {
      assert(0);
      throw EXCEPTION( MISSING_ARGUMENT, string("no sequence in fwd strand"));
    }
  }


  /*                                                                            
   * @return a Sequence object defined from beg to end coordinates.             
   * \todo No single annotated element is copied over, future dev.
   */
  BasicSeq *Sequence::getSubSequence(int beg, int end) {
    assert(end > beg && beg >= 1 && end <= length());

    std::ostringstream ostr;

    ostr << id() << "_" << beg << "_" << end;
    std::string seqId = ostr.str();
    int len = end - beg + 1;

    Sequence *c = new Sequence(seqId, NULL, 
                               alphabet(), false);

    for(int i = 0; i < NUM_STRANDS; i++) {
      TStrand strand = (TStrand)i;
      int offset = beg;

      if(strand == STRAND_COMP)
        offset = seqLen[strand] - end + 1;
      
      // set the character sequence
      char *region = NULL;
      if(this->seq[strand].length()) {
        region = new char[len + 1];
        strncpy(region, getSeq(strand) + offset - 1, len);
        region[len] = '\0';
        c->setSeq(region, strand);
        delete [] region;
      }        
    }

    return c;
  }

  /**                                                                           
   * @return a subsequence of this sequence defined by begin and end (and       
   * possibly upStream and downStream, if they are different from zero) as      
   * a char * object in region. It takes the strand into account.               
   */
  void Sequence::getLocation(char *region,
                             int upStream,
                             int downStream,
                             int begin, int end,
                             TStrand strand) {

    int regBeg, regEnd;
    assert(!circular); // not implemented yet                                   
    int len = end - begin + 1;

    if(strand == STRAND_COMP) {
      begin = seqLen[strand] - begin + 1;
      end = seqLen[strand] - end + 1;
      len = begin - end + 1;
    }

    if(!len) {
      region[0] = '\0';
      return;
    }

    if(end >= begin && strand != STRAND_COMP)
      strand = STRAND_FWD;
    else if(end <= begin && strand != STRAND_FWD)
      strand = STRAND_COMP;
    else {
      assert(0);
      ostringstream ostr;
      ostr << this->seqId << "[" << this->length() << "] : " << " one of " << begin << " or " << end << " are out of bounds";
      throw EXCEPTION(INVALID_CONTIG_COORDS, ostr.str());		      
    }

    if(begin <= 0 || end <= 0 || begin > seqLen[strand] || end > seqLen[strand]) {
      assert(0);
      ostringstream ostr;
      ostr << this->seqId << "[" << this->length() << "] : " << " one of " << begin << " or " << end << " are out of bounds";
      throw EXCEPTION(INVALID_CONTIG_COORDS, ostr.str());
    }

    if(strand == STRAND_FWD) {
      regEnd = end + downStream;
      regBeg = begin - upStream;

      if(regEnd > seqLen[strand])
        regEnd = seqLen[strand];
      else if(regEnd < 1)
        regEnd = 1;

      if(regBeg < 1)
        regBeg = 1;
      else if(regBeg > seqLen[strand])
        regBeg = seqLen[strand];

      if(regEnd < regBeg)  // truncated, zero length                            
        region[0] = '\0';
      else {
        strncpy(region, 
                seq[strand].c_str() + regBeg - 1, 
                regEnd - regBeg + 1);
        region[regEnd - regBeg + 1] = '\0';
      }
    }
    else {
      regBeg = begin + upStream;
      regEnd = end - downStream;

      if(regEnd > seqLen[strand])
        regEnd = seqLen[strand];
      else if(regEnd < 1)
        regEnd = 1;

      if(regBeg < 1)
        regBeg = 1;
      else if(regBeg > seqLen[strand])
        regBeg = seqLen[strand];

      if(regBeg < regEnd)  // truncated                                         
        region[0] = '\0';
      else {
        int lenR = regBeg - regEnd + 1;
        region[lenR] = '\0';
        //        if(seq[STRAND_COMP].length())
        strncpy(region, seq[strand].c_str() + seqLen[strand] - regBeg + 2, lenR);
        /*        else {
          strncpy(region, seq[STRAND_FWD].c_str() + regEnd - 1, lenR);
          Utils::revComp(region, sigma);
          } */
      }
    }
  }
  
  /** 
   * Appends SeqTags from Sequence c into *this. Tags from c are considered
   * to appear starting at position offset in *this. All SeqTag objects from
   * *this are supposed to end in a unique edge, so that SeqTag objects from
   * c are independent of *this's SeqTag objects.
   * The method requires that both c and *this have the same number of 
   * SeqTag objects
   * \todo Make it work for Seqtags with incompatible edge instances at the 
   * boundaries. Which edge to keep? nTag.front() or sTag.back()?
   * This happens when a longer word is being split into sTag and nTag. 
   * Return 0 if everything is ok, or the position at which the word starts
   * - it is the position at which the last syncNode3 appeared in each sTag
   * object
   */
  void Sequence::appendTags(TParseNode &bottlNode,
                            FSM &fsm,
                            Sequence &c, 
                            TSetType set,
                            int offset,
                            int limit) {

    unsigned int i;
    int phaseBreak = 0;

    if(!tags[set].size())
      tags[set].resize(c.getTags(set).size());

    assert(c.getTags(set).size() == tags[set].size());
    
    for(i = 0; i < c.getTags(set).size(); i++) {
      SeqTags &sTags = tags[set][i];
      SeqTags &nTags = c.getTags(set)[i];
      
      assert(nTags.size() >= 2);
      
      if(sTags.rank() < 0)
	sTags.setRank(nTags.rank());
      
      if(sTags.size() >= 2) {
        // remove edges at the boundaries if they do not have a signal
        
        Tag *tag = nTags.front();

        Edge *next_e = fsm.edge((TParseEdge)tag->getParseType());
        
        if(!next_e->hasSignal())
          nTags.pop_front();

        phaseBreak = nTags.initPhase();
        
      }
      else // only initialize phase if this sequence has no tags associated
        sTags.setInitPhase(nTags.initPhase());

      // append Tag objects from nTags to sTags, one at a time
      SeqTags::iterator bit = nTags.begin();
      //      cerr << "init phase " << nTags.initPhase() << endl;

      for( ; bit != nTags.end(); bit++) {
        Tag &tag = *(*bit);
        int begin = tag.getPos(), end = begin - 1;
        //                cerr << "orig " << tag.getGEClass() << " " << tag.getParseType() << " " << begin << " " << end << " " << limit << " clen " << c.length() << " " << tag.getStrand() << endl;

        if(tag.getGEClass() == NODE_INST &&
           tag.getParseType() == bottlNode) {

          NodeInst &node = (NodeInst &)tag;
          end = begin + node.getLen() - 1; 
          
          if(end > limit) {
            end = limit;
            node.setLen(end - begin + 1);
          }
        }

        if(begin > limit || end > limit)
          continue;
        
        int nTagBreak = nTags.phaseBreak(*bit);

        tag.setPos(tag.getPos() + offset);

        //        cerr << "Inserted from " << c.id() << " " << tag.getParseType() << " " << tag.getPos() << endl;
        sTags.push_back(&tag);

        if(tag.getGEClass() != NODE_INST)
          continue;
        
        if(phaseBreak > 0) {
          sTags.insertPhaseBreak(sTags.back(), phaseBreak);
          //          cerr << " inserted break @ " << sTags.back()->getPos() << " " << sTags.back()->getParseType() << " " << phaseBreak << endl;
          phaseBreak = 0;
        }        
        else if(nTagBreak > 0)
          sTags.insertPhaseBreak(sTags.back(), nTagBreak);
        
      }
      //      cerr << "end phase " << nTags.endPhase() << endl;
      sTags.setEndPhase(nTags.endPhase());

    }
  }

  pair<int, int> EdgeAnnotSeq::findEdge(pair<int, int> &coords, 
					TStrand strand) {
					
    vector<int> *sseq = seq[strand];
    pair<int, int> p(-1, -1);
    if(!sseq)
      return p;
        
    vector<int> &ibackedges = sseq[coords.first];
    vector<int> &jbackedges = sseq[coords.second];
    for(int j = 0; j < jbackedges.size(); j += 3) {
      if(jbackedges[j + 1] == coords.first)
	for(int i = 0; i < ibackedges.size(); i += 3) {   
	  if(ibackedges[i] == jbackedges[j] &&
	     ibackedges[i + 2] == jbackedges[j + 2])
	    return pair<int, int>(i, j);
	}
    }
    return p;
  }

  EdgeAnnotSeq& EdgeAnnotSeq::operator= (const EdgeAnnotSeq & c) {
    EdgeAnnotSeq & myc = (EdgeAnnotSeq &)c;
    _numPhases = myc.numPhases();    
    
    for(int t = 0; t < NUM_STRANDS; t++) {
      TStrand strand = (TStrand)t;
      vector<int>* e = NULL;
      
      if(myc.getSeq(strand)) {
        e = new vector<int> [myc.length(strand) + 1];
        
        for(int pos = 0; pos <= myc.length(strand); pos++)
	  e[pos] = myc.getSeq(strand)[pos];

      }

      this->setSeq(e, myc.length(strand), strand);
    }
  }


  BasicSeq *EdgeAnnotSeq::getSubSequence(int beg, int end) {
    assert(end > beg && beg >= 1 && end <= length());


    std::ostringstream ostr;  
    ostr << id() << "_" << beg << "_" << end;
    std::string seqId = ostr.str();
    int len = end - beg + 1;
 
    EdgeAnnotSeq *c = new EdgeAnnotSeq(seqId, _numPhases,
                                       alphabet(), false);
    
    for(int t = 0; t < NUM_STRANDS; t++) {
      TStrand strand = (TStrand)t;
      int offset = beg;

      if(strand == STRAND_COMP)
        offset = seqLen[strand] - end + 1;
            
      // set the edge annotation information      
      vector<int>* subseq = NULL;

      if(this->getSeq(strand)) {
        subseq = new vector<int> [len + 1];
        for(int i = 1; i <= len; i++) {
          subseq[i] = getSeq(strand)[offset + i - 1];
	  // changing relative positions
	  for(int j = 2; j < subseq[i].size(); j += 3) {
	    subseq[i][j] = subseq[i][j] - offset + 1;
	    if(subseq[i][j] <= 0)
	      subseq[i][j] = -1;
	  }
        }
      }
      
      c->setSeq(subseq, len, strand);

    }
    
    return c;

  }
  
}
