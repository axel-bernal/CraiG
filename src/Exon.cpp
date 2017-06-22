#include "Exon.h"
#include "Transcript.h"
#include "Gene.h"

/****************************************************************************
 * Exon.cpp - part of the craig namespace, a genomics library
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

    Exon::Exon(Transcript * ownerTranscript,
               int b, int e, int phase,
               char * seq, double codScore) {

        assert(b <= e && phase >= 0 && phase < 3);
        _end = e;
        _begin = b;
        _phase = phase;
        _transcript = ownerTranscript;

        if(seq)
            this->_exonSeq = seq;

        //    for(int i = 0; i < 7; i++)
        //      _scores[i] = -1;

        _codingScore = codScore;

    }

    Exon::Exon(const Exon &e) {
        *this = (const Exon &)e;
    }


    Exon & Exon::operator=(const Exon & exon) {
        Exon &e  = (Exon &)exon;
        this->_begin = e.begin();
        this->_end = e.end();
        this->_phase = e.phase();
        this->_transcript = e.transcript();
        //    this->_overlapRegions = e.overlapRegions();
        //    this->_overlapExons = e.overlapExons();
        //    this->_simFacts = e.simFacts();
        this->_exonSeq = std::string(e.exonSeq());
        this->_codingScore = e.codingScore();

        //    for(int i = 0; i < 7; i++)
        //      _scores[i] = e.scores()[i];

        return *this;
    }

    const char * Exon::exonSeq() {
        Transcript *t = _transcript;

        if(!_exonSeq.length()) {
            char *seq;
            seq = new char[abs(_end - _begin) + 2];

            if(t->isInOneStrand())
                (t->getSequence())->getLocation(seq, 0, 0, _begin, _end);
            else
                (t->getSequence())->getLocation(seq, 0, 0, _begin, _end, t->getStrand());

            _exonSeq = std::string(seq);

            delete [] seq;
        }
        assert(_exonSeq.length() == (unsigned)abs(_end - _begin) + 1);

        return _exonSeq.c_str();
    }

    /*
      void Exon::filterOlapExons() {
      list<Exon *>::iterator oit = _overlapExons.begin();
      list<Overlap *>::iterator it = _overlapRegions.begin();

      while(it != _overlapRegions.end() && oit != _overlapExons.end()) {
      Gene *g = ((*oit)->transcript())->gene();

      if(g->state() == OVLAPD_ORF || g->state() & IS_GENE) {
      it++;
      oit++;
      }
      else {
      oit = _overlapExons.erase(oit);
      it = _overlapRegions.erase(it);
      }
      }
      assert( it == _overlapRegions.end() && oit == _overlapExons.end());
      }

      int Exon::getMaxFrameScore(int frame) {
      int max = abs(frame) + 3*(frame < 0);

      for(int i = 0; i < 7; i++) {
      if(_scores[i] > _scores[max]*1.002)
      max = i;
      }

      assert(_scores[max] > 0);

      if(max > 3)  return 3-max;
      else  return max;
      } */

    Exon::~Exon() {

        //    if(!_overlapRegions.empty())
        //      _overlapRegions.erase(_overlapRegions.begin(), _overlapRegions.end());

        //    if(!_overlapExons.empty())
        //      _overlapExons.erase(_overlapExons.begin(), _overlapExons.end());

        //    if(!_simFacts.empty())
        //      _simFacts.erase(_simFacts.begin(), _simFacts.begin());

    }
}
