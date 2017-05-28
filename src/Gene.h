/****************************************************************************
 * Gene.h - part of the craig namespace, a genomics library
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


#ifndef _GENE_H_
#define _GENE_H_

#include "IMM.h"
#include <assert.h>
#include "Sequence.h"
#include "Sim.h"
#include "BioFeature.h"
#include "Transcript.h"
#include <ctype.h>
#include <string.h>
#include <list>
#include <vector>
#include <algorithm>
#include <functional>

using namespace craig;

namespace craig {

    //! enumerate for the different states a gene can be in
    typedef enum {  TRUNCATED_ORF      =  0x01000000,
                    ORF                =  0x01000001,
                    SHDWED_ORF         =  0x01000010,
                    OVLAPD_ORF         =  0x01000100,
                    HOMOLOG_ORF        =  0x01001000,
                    OTHER_ORF          =  0x01010000,
                    SUSPECT_GENE       =  0x10100001,
                    SUSPECT_SHDWEDGENE =  0x10100010,
                    SUSPECT_OVLAPDGENE =  0x10100100,
                    SUSPECT_HOMOLOGENE =  0x10101000,
                    GENE               =  0x10000001,
                    SHDWED_GENE        =  0x10000010,
                    OVLAPD_GENE        =  0x10000100,
                    HOMOLOG_GENE       =  0x10001000,
                    OTHER_GENE         =  0x10010000,
                    TURN_SUSPECT_ON    =  0x00100000,
                    IS_GENE            =  0x10000000,
                    IS_ORF             =  0x01000000
    } TGene;


    /**
     * The Gene class derives from the class BioFeature and is represented as a
     * vector of transcripts. It also mantains two variables for specifying the
     * coordinates of the gene  which are automatically updated everytime
     * transcripts and removed or added.
     ***************************************************************************/

    class Gene : public BioFeature {
    private:
        TGene _state;
        vector<Transcript> _transcripts;
    public:
        Gene(std::string &,
             Sequence *annotSeq,
             bool inOneStrand,
             TStrand strand = STRAND_FWD
            );
        Gene(const Gene& o);
        Gene() { }

        inline TGene state() {  return _state; }
        inline void setState(TGene state) {  _state = state; }
        inline Transcript & transcript(int i) {  return _transcripts[i]; }
        inline vector<Transcript> & transcripts() {  return _transcripts; }
        inline bool hasTranscripts() {  return _transcripts.size() != 0; }
        inline void updateLocation(Transcript &t) {
            if(_begin > t.begin()) _begin = t.begin();
            if(_end < t.end()) _end = t.end();
        }

        static Gene *findGene(std::string &, list<Gene> &);
        std::string & setId(std::string &);
        Gene& operator= (const Gene& o);
        bool operator< (const Gene& o);
        bool operator< (const Gene& o) const;
        bool operator== (const Gene& o);
        void makeSuspect();
        Transcript &addTranscript(Transcript &t);
        friend int gap(Gene &, Gene &);
        void toOneStrand();
        void toTwoStrand();
        void addSimFact(Sim *s);
        bool isAlternativeTranscript(Transcript &t);

        ~Gene() {
            if(!_transcripts.empty())
                _transcripts.erase(_transcripts.begin(),_transcripts.end());
        }
    };
}

namespace std {
    /**
     * Routine for sorting STL containers instantiated with Gene* types
     */
    template <>
        struct greater<Gene *> {
        bool operator()(const Gene* o1, const Gene* o2)
        {
            //Defined for ascending sorting
            if(!o1)
                return true;
            if(!o2)
                return false;
            return ((Gene &)(*o1)) < (*o2);
        }
    };
}
#endif
