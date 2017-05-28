/****************************************************************************
 * Transcript.h - part of the craig namespace, a genomics library
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


#ifndef _TRANSCRIPT_H_
#define _TRANSCRIPT_H_

#include "Exon.h"
#include "BioFeature.h"

using namespace craig;

namespace craig {

    // forward declarations
    class Gene;

    /**
     * The Transcript class derives from the class BioFeature and is represented
     * as a vector of exons (including noncoding ones).
     * It also contains information about whether the transcript has a
     * translation start, stop and signal peptide.
     ***************************************************************************/

    class Transcript : public BioFeature {
    private:
        Gene *_gene; //!< owner object of type Gene
        vector<Exon> _5exons;
        vector<Exon> _exons;
        vector<Exon> _3exons;
        bool _hasStart;
        bool _hasStop;
        bool _hasPeptide;
        int _phase;
        unsigned int _transcriptLen;
        ::string _5utrSeq, _codingSeq, _3utrSeq;

    public:
        Transcript(std::string & iden,
                   Gene *gene,
                   bool hasStart = true,
                   bool hasStop = true,
                   bool hasPeptide = false,
                   int phase = 0
            );
        Transcript(const Transcript &t);
        Transcript() { }

        inline int phase() {  return _phase; }
        inline Gene *gene() {  return _gene; }
        inline bool hasStart() {  return _hasStart; }
        inline bool hasStop() {  return _hasStop; }
        inline bool hasPeptide() {  return _hasPeptide; }
        inline vector<Exon>& FPexons() {  return _5exons; }
        inline vector<Exon>& exons() {  return _exons; }
        inline vector<Exon>& TPexons() {  return _3exons; }
        inline Exon & exon(int numExon = 0) {  return _exons[numExon]; }
        inline unsigned int numExons() {  return _exons.size(); }
        inline unsigned int transcriptLen() {
            return _transcriptLen;
        }
        inline void setGene(Gene *gene) {  _gene = gene; }

        inline void setPeptide(bool peptidePresent) {
            _hasPeptide = peptidePresent;
        }
        inline void setStart(bool startPresent) {
            _hasStart = startPresent;
        }
        inline void setStop(bool stopPresent) {
            _hasStop = stopPresent;
        }
        inline void setPhase(int phase) {  _phase = phase; }
        inline Exon & add5Exon(int b, int e, char *exonSeq = NULL) {
            _5utrSeq.clear();
            return addExon(Exon(this, b, e, 0, exonSeq), _5exons);
        }
        inline Exon & addExon(int b, int e, int phase, char *exonSeq = NULL) {
            _codingSeq.clear();
            return addExon(Exon(this, b, e, phase, exonSeq), _exons);
        }
        inline Exon & add3Exon(int b, int e, char *exonSeq = NULL) {
            _3utrSeq.clear();
            return addExon(Exon(this, b, e, 0, exonSeq), _3exons);
        }
        inline void updateLocation(Exon &e) {
            if(this->_begin > e.begin()) this->_begin = e.begin();
            if(this->_end < e.end()) this->_end = e.end();
        }
        inline int codingBegin() {
            return (_exons.front().begin() > _exons.back().end()) ?
                _exons.back().end() :
                _exons.front().begin();
        }
        inline int codingEnd() {
            return (_exons.front().begin() > _exons.back().end()) ?
                _exons.front().begin() :
                _exons.back().end();
        }

        Exon & addExon(Exon e, vector<Exon> &exons);
        void setExons(vector<Exon> &theseExons, vector<Exon> & exons);
        Transcript & operator=(const Transcript & t);
        bool operator==(const Transcript & t);
        bool operator<(const Transcript & t);
        const char *fiveUTRSeq();
        const char *codingSeq();
        const char *threeUTRSeq();
        void toOneStrand();
        void exonsToOneStrand(vector<Exon> & exons);
        void toTwoStrand();
        void exonsToTwoStrand(vector<Exon> & exons);
        //    void addSimFact(Sim *s);
        void getSplicedSeq(std::string &seq, vector<Exon> & exons);
        bool getNeighboorExons4Intron(pair<int, int> &intron,
                                      Exon* &lexon, Exon* &rexon);
        bool getNeighboorExons4Intron(vector<Exon> &exons,
                                      pair<int, int> &intron,
                                      Exon* &lexon, Exon* &rexon);
        void extractIntrons(vector<Exon> &exons, vector<pair<int, int> > &introns);
        void extractIntrons(vector<pair<int, int> > &introns);

        ~Transcript() {  }
    };
}

namespace std {
    /**
     * Routine for sorting STL containers instantiated with Transcript* types
     */
    template <>
        struct greater<Transcript *> {
        bool operator()(const Transcript* t1, const Transcript* t2)
        {
            //Defined for ascending sorting
            if(!t1)
                return true;
            if(!t2)
                return false;
            return ((Transcript &)(*t1)) < (*t2);
        }
    };
}

#endif
