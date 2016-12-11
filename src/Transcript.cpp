#include "Transcript.h"
#include "Gene.h"

/****************************************************************************
 * Trancript.cpp - part of the craig namespace, a genomics library
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
 ****************************************************************************/


namespace craig {

  Transcript::Transcript(std::string & iden, 
                         Gene *gene, 
                         bool hasStart, 
                         bool hasStop,
                         bool hasPeptide, 
                         int phase) 
    : BioFeature(iden, 
                 gene->getSequence(), 
                 gene->isInOneStrand(),
                 gene->getStrand()) {
    
    this->_gene = gene;
    this->_hasStart = hasStart;
    this->_hasStop = hasStop;
    this->_hasPeptide = hasPeptide;
    this->_phase = phase;
    this->_transcriptLen = 0;
  }

  Transcript::Transcript(const Transcript & t) {
    *this = (Transcript &)t;
  }             


  Exon & Transcript::addExon(Exon ex, vector<Exon> & exons) {
    vector<Exon>::iterator it = exons.size() ? exons.begin() : exons.end();

    while(it != exons.end() && *it < ex)
      it++;
    
    _transcriptLen += abs(ex.end() - ex.begin()) + 1;
      
    it = exons.insert(it, (const Exon &)ex);
    it->setPhase(_phase + it->phase() % 3);
    it->setTranscript(this);

    if(this->_begin > it->begin()) 
      this->_begin = it->begin();

    if(this->_end < it->end())
      this->_end = it->end();  

    _gene->updateLocation(*this);

    return *it;
  //   cout << getId() << " " << ex.beg << " " << ex.end << " " << strand <<"\n";
  }

  void Transcript::setExons(vector<Exon> &theseExons, vector<Exon> & exons) {
    theseExons = (const vector<Exon>)exons;
    vector<Exon>::iterator it = theseExons.begin();

    while(it != theseExons.end()) { // changing the owner transcript...
      it->setTranscript(this);
      it++;
    }
  }

  Transcript & Transcript::operator=(const Transcript & transcript) {
    Transcript &t = (Transcript &)transcript;
    this->_gene = t.gene();
    (BioFeature &)*this = (BioFeature &)t;
    this->_gene = t.gene();
    
    setExons(_5exons, t.FPexons());
    setExons(_exons, t.exons());
    setExons(_3exons, t.TPexons());
    
    _transcriptLen = t.transcriptLen();
    _5utrSeq = t.fiveUTRSeq();
    _codingSeq = t.codingSeq();
    _3utrSeq = t.threeUTRSeq();
    _hasStart = t.hasStart();
    _hasStop = t.hasStop();
    _hasPeptide = t.hasPeptide();
    _phase = t.phase();

    return *this;
  }

  bool Transcript::operator==(const Transcript &transcript) {
    Transcript &t = (Transcript &)transcript;
    unsigned int i = 0;
    bool equal = (_phase == t.phase());

    if(!equal)
      return equal;

    while(i++ < _5exons.size()) 
      equal = equal && (_5exons[i] == t.FPexons()[i]);

    i = 0;

    while(i++ < _exons.size()) 
      equal = equal && (_exons[i] == t.exon(i));

    i = 0;

    while(i++ < _3exons.size()) 
      equal = equal && (_3exons[i] == t.TPexons()[i]);
    
    return equal && (BioFeature &)*this == (BioFeature &)t;
  }

  /**
   * Operator<. Coordinates for Transcripts should be in same format
   */
  bool Transcript::operator<(const Transcript & transcript) {
    Transcript &t = (Transcript &)transcript;

    if(inOneStrand != t.isInOneStrand()) {
      assert(0);
      throw EXCEPTION(INCOMPAT_TRANSCRIPT, std::string(t.getId()));
    }

    if(getSequence() == t.getSequence()) {
      return begin() < t.begin();
    }
    else 
      return getSequence() < t.getSequence();
  }

  const char * Transcript::fiveUTRSeq() {
    if(_5utrSeq.length() == 0) 
      getSplicedSeq(_5utrSeq, _5exons);
    
    return _5utrSeq.c_str();
  }

  const char * Transcript::threeUTRSeq() {
    if(_3utrSeq.length() == 0) 
      getSplicedSeq(_3utrSeq, _3exons);
    
    return _3utrSeq.c_str();

  }

  const char * Transcript::codingSeq() {
    if(_codingSeq.length() == 0) 
      getSplicedSeq(_codingSeq, _exons);

    return _codingSeq.c_str();
  }        

  /**
   * This function will transform the coordinates of the exons of a gene 
   * in complementary strand to the forward strand, meaning all exon 
   * coordinates will be read from the start of the forward annotSeq
   */
  void Transcript::toOneStrand() {
    if(inOneStrand || strand == STRAND_FWD)
      return;

    this->_begin = INT_MAX;
    this->_end = -1;

    exonsToOneStrand(_5exons);
    exonsToOneStrand(_exons);
    exonsToOneStrand(_3exons);

    inOneStrand = true;
  }

  void Transcript::exonsToOneStrand(vector<Exon> & exons) {
    unsigned int size = exons.size();
    if(!size)
      return;

    for(unsigned int i = 0; i < size; i++) {
      int begin = exons[i].begin();
      exons[i].setBegin(getSequence()->length() - exons[i].end() + 1);
      exons[i].setEnd(getSequence()->length() - begin + 1);
      updateLocation(exons[i]);
    }

    // reversing exons
    std::reverse(exons.begin(), exons.end());
  //  cout << exons[0].beg << " " << exons[0].end << endl;

  }

  void Transcript::toTwoStrand() {
    if(!inOneStrand || strand == STRAND_FWD)
      return;

    this->_begin = INT_MAX;
    this->_end = -1;

    exonsToTwoStrand(_5exons);
    exonsToTwoStrand(_exons);
    exonsToTwoStrand(_3exons);

    inOneStrand = false;
  }

  void Transcript::exonsToTwoStrand(vector<Exon> & exons) {
    unsigned int size = exons.size();
    if(!size) 
      return;

    // reversing exons
    std::reverse(exons.begin(), exons.end());

    for(unsigned int i = 0; i < size; i++) {
      int begin = exons[i].begin();
      exons[i].setBegin(getSequence()->length() - exons[i].end() + 1);
      exons[i].setEnd(getSequence()->length() - begin + 1);
      updateLocation(exons[i]);
    }
  }


  //  void Transcript::addSimFact(Sim *s) {
  //    // TODO: determine which exon overlaps the similarity entry 
  //    int exonNumber = 0;
  //    _exons[exonNumber].addSimFact(s);
  //  }

  void Transcript::getSplicedSeq(std::string &seq, vector<Exon> & exons) {
    Sequence *c = _gene->getSequence();
    unsigned int i;

    if(c == NULL) {
      assert(0);
      throw EXCEPTION(CONTIG_UNAVAILABLE, " ");
    }

    for(i = 0; i < exons.size(); i++) {
      assert(seq.length() + abs(exons[i].end() - exons[i].begin()) + 1 <= transcriptLen());
      const char *exonSeq = exons[i].exonSeq();
      seq = seq + exonSeq;
    }
    if(inOneStrand && strand == STRAND_COMP) {
      for(i = 0; i < seq.length(); i++)
        seq[i] = c->alphabet()->toComplement(seq[i]);
      Utils::reverse(seq);
    }
  }

  bool Transcript::getNeighboorExons4Intron(pair<int, int> &intron, 
					    Exon* &lexon, Exon* &rexon) {
    bool found = getNeighboorExons4Intron(_5exons, intron, lexon, rexon);
    if(!found)
      found =  getNeighboorExons4Intron(_exons, intron, lexon, rexon);
    if(!found)
      found =  getNeighboorExons4Intron(_3exons, intron, lexon, rexon);
    return found;
  }

  bool Transcript::getNeighboorExons4Intron(vector<Exon> &exons, 
					    pair<int, int> &intron, 
					    Exon* &lexon, Exon* &rexon) {
    int i;
    if(_gene->getStrand() == STRAND_FWD) {
      for(i = 1; i < exons.size() ; i++)
	if(intron.first == exons[i - 1].end() + 1 && 
	   intron.second == exons[i].begin() - 2) {
	  lexon = &exons[i - 1];
	  rexon = &exons[i];
	  return true;
	}
    }
    else {
      for(i = exons.size() - 2; i >= 0; i--)
	if(intron.first == exons[i].end() + 2 &&
	   intron.second == exons[i + 1].begin() - 1) {
	  lexon = &exons[i + 1];
	  rexon = &exons[i];
	  return true;
	}
    }
    return false;
  }
  
  void Transcript::extractIntrons(vector<pair<int, int> > & introns) {
    extractIntrons(_5exons, introns);
    extractIntrons(_exons, introns);
    extractIntrons(_3exons, introns);
  }

  void Transcript::extractIntrons(vector<Exon> &exons,
				 vector<pair<int, int> > & introns) {
    int i;
    if(_gene->getStrand() == STRAND_FWD)
      for(i = 1; i < exons.size() ; i++) {
	pair<int, int> p(exons[i - 1].end() + 1, exons[i].begin() - 2);
	introns.push_back(p);
      }
    else 
      for(i = exons.size() - 2; i >= 0; i--) {
	pair<int, int> p(exons[i].end() + 2, exons[i + 1].begin() - 1);
	introns.push_back(p);
      }
  }
}
