
#include "IMM.h"
#include "Utils.h"
#include "Gene.h"

/****************************************************************************
* Gene.cpp - part of the craig namespace, a genomics library
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


                                                                              
namespace craig {
  
  Gene::Gene(std::string &iden,
             Sequence *annotSeq,
             bool inOneStrand,
             TStrand strand) 
    : BioFeature(iden, 
                 annotSeq, 
                 inOneStrand, 
                 strand) {

    _state = ORF;
  }                        
  
  Gene::Gene(const Gene& o) {
    *this = (Gene &)o;
  }
  
  
  std::string & Gene::setId(std::string & iden) {
    id = iden;
    return id;
  }                           
  
  Gene *Gene::findGene(std::string &id, list<Gene> & genes) {
    list<Gene>::iterator it;

    for(it = genes.begin(); it != genes.end(); it++)
      if(!id.compare(it->getId()))
        return &(*it);

    return NULL;
  }
  
  Gene& Gene::operator=(const Gene& gene) {
    Gene &g = (Gene &)gene;
    (BioFeature &)*this = (BioFeature &)g;
    this->_begin = g.begin();
    this->_end = g.end();
    _state = g.state();

    return *this;
  }                
  
  bool Gene::operator<(const Gene& gene) {
    Gene &g = (Gene &)gene;
    int b = getSequence()->id().compare(g.getSequence()->id());

    if(b)
      return b < 0;
    
    int begin1 = begin();
    int begin2 = g.begin();
    
    if(strand == g.getStrand()) {
      assert(isInOneStrand() == g.isInOneStrand());

      if(strand == STRAND_FWD) 
        return (begin1 <= begin2);
      else
        if(isInOneStrand())
          return (begin1 <= begin2);        
        else
          return (begin1 > begin2);

    }
    else {

      if(strand == STRAND_FWD) {
        if(!g.isInOneStrand())
          begin2 = getSequence()->length() - g.end() + 1;
      }
      else
        if(!isInOneStrand())
          begin1 = getSequence()->length() - end() + 1;

      return (begin1 <= begin2);

    }
  }
  
  bool Gene::operator==(const Gene& gene) {
    Gene &g = (Gene &)gene;
    bool equal = true;

    if(getSequenceId().compare(g.getSequenceId()))
      return 0;

    unsigned int i = 0;

    while(i < _transcripts.size()) {
      equal = equal && (_transcripts[i] == g.transcript(i));
      i++;
    }

    return equal && (_transcripts.size() == g.transcripts().size());
  }                                       
  
  
  void Gene::makeSuspect() {
    assert(_state & IS_GENE);
    _state = (TGene) (_state | TURN_SUSPECT_ON);
  }                    
  
  Transcript &Gene::addTranscript(Transcript & t) {
    assert(t.getStrand() == getStrand());
    
    if(_transcripts.empty())
      inOneStrand = t.isInOneStrand();
    else {
      if(inOneStrand != t.isInOneStrand()) {
        if(inOneStrand) t.toOneStrand();
        else t.toTwoStrand();
      }
    }
    
    vector<Transcript>::iterator it = _transcripts.size() ? _transcripts.begin() : _transcripts.end();
    while(it != _transcripts.end() && *it < t)
      it++;
    
    it = _transcripts.insert(it, (const Transcript &)t);
    it->setGene(this);
    updateLocation(t);
    
    return *it;
  }
  
  /**
   * Computes the gap between g1 and g2. If there is a gap in between g1
   * and g2 it will return a positive amount, otherwise return value will
   * be negative or zero if no exons are present in either gene 
   */
  int gap(Gene & g1, Gene & g2) {
    unsigned int i;
    if(!g1.isInOneStrand() || !g2.isInOneStrand()) {
      assert(0);
      throw EXCEPTION(INCOMPAT_TRANSCRIPT, string(g1.getId()) + string(g2.getId()));
    }

    Transcript &t1 = g1.transcripts()[0];
    Transcript &t2 = g2.transcripts()[0];

    if(t1.begin() < t2.begin())  {
      int end = INT_MAX;

      for(i = 0; i < g1.transcripts().size(); i++) 
        if(end < g1.transcripts()[i].end())
          end = g1.transcripts()[i].end();

      return t2.begin() - end;

    }
    else {

      int end = INT_MAX;

      for(i = 0; i < g2.transcripts().size(); i++) 
        if(end < g2.transcripts()[i].end())
          end = g2.transcripts()[i].end();

      return t1.begin() - end;

    }
  }
  
  void Gene::toOneStrand() {
    if(inOneStrand || strand == STRAND_FWD)
      return;
    
    this->_begin = INT_MAX;
    this->_end = -1;
    
    for(unsigned int i = 0; i < _transcripts.size(); i++)
      _transcripts[i].toOneStrand();
    
    std::reverse(_transcripts.begin(), _transcripts.end());
    
    for(unsigned int i = 0; i < _transcripts.size(); i++) {
      if(this->_begin > _transcripts[i].begin()) 
        this->_begin = _transcripts[i].begin();
      
      if(this->_end < _transcripts[i].end()) 
        this->_end = _transcripts[i].end(); 
    }
    
    inOneStrand = true;
  }
  
  void Gene::toTwoStrand() {
    if(!inOneStrand || strand == STRAND_FWD)
      return;
    
    this->_begin = INT_MAX;
    this->_end = -1;
    
    for(unsigned int i = 0; i < _transcripts.size(); i++)
      _transcripts[i].toTwoStrand();
    
    std::reverse(_transcripts.begin(), _transcripts.end());
    
    for(unsigned int i = 0; i < _transcripts.size(); i++) {
      if(this->_begin > _transcripts[i].begin()) 
        this->_begin = _transcripts[i].begin();

      if(this->_end < _transcripts[i].end())
        this->_end = _transcripts[i].end(); 
    }

    inOneStrand = false;
  }
  
  void Gene::addSimFact(Sim *s) {  
    for(unsigned int i = 0; i < _transcripts.size(); i++)
      ;//      _transcripts[i].addSimFact(s);
  }
  
  bool Gene::isAlternativeTranscript(Transcript &transcript) {
    if(transcript.getStrand() != getStrand())
      return false;
    
    for(unsigned int i = 0; i < _transcripts.size(); i++) {
      Transcript & t = _transcripts[i];
      
      if(t.begin() < transcript.end() && t.end() > transcript.begin())
        return true;
      if(transcript.begin() < t.end() && transcript.end() > t.begin())
        return true;
    }
    
  return false;
  }
}
