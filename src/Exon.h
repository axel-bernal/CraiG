/****************************************************************************
 * Exon.h - part of the craig namespace, a genomics library
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

#ifndef _EXON_H_
#define _EXON_H_

#include "Overlap.h"
#include <list>
#include "Sim.h"


namespace craig {
  //forward declarations
  class Transcript;

  /**
   * The Exon class represents an exonic fragment in a transcript. It has a
   * phase,  two coordinates and a pointer to the owner Transcript object 
   * (declared forward) here. 
   ***************************************************************************/
  
  class Exon {
   private:
    Transcript *_transcript; //!< owner object of type Transcript;
    int _begin;
    int _end;
    int _phase;
    /*    list<Overlap *> _overlapRegions;
    list<Exon *> _overlapExons;
    list<Sim *> _simFacts;
    double _scores[7];  */
    double _codingScore;
    std::string _exonSeq;

   public:
    Exon(Transcript *transcript,
         int beg,
         int end,
         int phase,
         char *exonSeq = NULL,
         double codScore = -1
       );

    Exon(const Exon &);

    Exon & operator=(const Exon & e);
    const char * exonSeq();
    void filterOlapExons();  
    int getMaxFrameScore(int frame = 0);
    
    /* Exons should belong to the same transcript */
    
    inline Transcript *transcript() {  return _transcript; }
    inline int begin() {  return _begin; }
    inline int end() { return _end; }
    inline  int phase() { return _phase; }
    //    inline int numSims() {  return _simFacts.size(); }
    //    inline void addSimFact(Sim * s) {  _simFacts.push_back(s); }
    inline bool operator<(Exon & e) {  
      return _begin == e.begin() ? _end < e.end() : _begin < e.begin(); 
    }  
    inline bool operator==(Exon & e) {
      return(_begin == e.begin() && _end == e.end());
    }
    /*    inline void setFrameScore(int frame, double val) {
      _scores[abs(frame) + 3*(frame <0)] = val;
    }
    inline double * scores() {
      return _scores;
    }
    */
    inline void setTranscript(Transcript *transcript) {  _transcript = transcript; }
    inline void setBegin(int begin) {  _begin = begin; }
    inline void setEnd(int end) {  _end = end; }
    inline void setPhase(int phase) {  _phase = phase; }
/*    inline list<Overlap *>& overlapRegions() {
      return _overlapRegions;
    }
    inline list<Exon *>& overlapExons() {
      return _overlapExons;
    } 
    inline list<Sim *>& simFacts() {
      return _simFacts;
    }
    inline friend void markOverlap(Exon *exonOit, Overlap *ov, Exon *exonIt) {
      exonOit->overlapRegions().push_back(ov);
      exonOit->overlapExons().push_back(exonIt);
      exonIt->overlapRegions().push_back(ov);
      exonIt->overlapExons().push_back(exonOit);
    }
    inline bool hasOverlap() {
      return _overlapRegions.size() == 0 ? false : true;
      } */
    inline double codingScore() {
      return _codingScore;
    }
    inline void setCodingScore(double cScore) {
      _codingScore = cScore;
    }

    ~Exon();
  };
}

#endif
