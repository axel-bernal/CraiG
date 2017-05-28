/****************************************************************************
 * GeneUtils.h - part of the craig namespace, a genomics library
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

#ifndef _GEN_UTIL_GF_H
#define _GEN_UTIL_GF_H

#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <map>
#include <streambuf>
#include "GenLibExcept.h"
#include <sys/stat.h>
#include <vector>
#include "Utils.h"
#include "Tag.h"
#include "FSM.h"
#include "Gene.h"
#include "FeatureEngine.h"
#include "Sequence.h"
#include "RepeatedElem.h"

#define MAX_STOPS_PER_FRAME 1000000

namespace craig {


    /**
     * The GeneUtils class implements routines for converting Tag to Gene
     * objects and visceversa. It also contains utility functions for
     * Gene, Transcripts and Exon objects
     ***************************************************************************/

    class GeneUtils {

    protected:
        static EdgeInst *_getSignal(Edge *fsmEgde, int sigPos,
                                    TypedFilter<EdgeInst> **signals,
                                    string msg);

    public:
        GeneUtils() {

        }

        static void loadIGenicRegions(map<string, list<BioFeature> > &igenics,
                                      string &igenics_fn);

        static void getGenicRegions(list<BioFeature> &igenics,
                                    list<BioFeature> &genics,
                                    int percent_igenic = 50);

        static void genes2transcripts(Sequence &c, TSetType typeSet,
                                      list<Transcript *> &transcripts);

        static void loadGenes(list<Gene> & listGenes, const char *geneFile,
                              bool buildstartInfo, bool readfastaSequence,
                              list<Sequence *> &annotSeqs, Sigma *sigma,
                              TSetType typeSet = DEFAULT_SET);

        static void extractFrameOrfs(vector<pair<int,int> > &orfs,
                                     Sequence &c,
                                     TypedFilter<EdgeInst> **signals,
                                     int frame, TStrand strand,
                                     int mingeneLen = 0,
                                     bool reportTruncated = false,
                                     bool includeAllStarts = false);

        static void findNextGene(Sequence &c, SeqTags &st, SeqTags::iterator it,
                                 int &beg, int &end, TStrand &strand);

        static void updateTranscript(Transcript *t, bool hasStart, bool hasStop,
                                     bool hasPeptide, int phase);
        static void tags2Genes(Sequence &c,
                               SeqTags & listTags,
                               const char *idPrefix,
                               FSM &,
                               list<Gene> & listGenes,
                               int offset = 0,
                               bool reportNonCodingGenes = false);

        static void updateAltTranscript(list<Gene> &listGenes, Gene &gene,
                                        Transcript *t, bool hasStart,
                                        bool hasStop, bool hasPeptide,
                                        int phase);

        static void tags2Genes(Sequence &c,
                               vector<SeqTags> & listTags,
                               const char * idPrefix,
                               FSM &,
                               list<Gene> & listGenes,
                               int offset = 0,
                               bool reportNonCodingGenes = false);

        static void exonStructure2Tags(vector<Exon> & v,
                                       FSM *fsm,
                                       TypedFilter<UCHAR> *contexts,
                                       TypedFilter<EdgeInst> **signals,
                                       Sequence &c,
                                       Edge **fsmEdge,
                                       Node **fsmNode,
                                       std::string &intronNode,
                                       std::string &lintronNode,
                                       const char *singleExonNode,
                                       const char *initExonNode,
                                       const char *innerExonNode,
                                       const char *lastExonNode,
                                       SeqTags &listGEs,
                                       int lExonEnd,
                                       int nExonBeg,
                                       TStrand strand
            );
        static void annotSeq2Tags(vector<SeqTags> &,
                                  FSM *fsm,
                                  TypedFilter<UCHAR> *contexts,
                                  TypedFilter<EdgeInst> **signals,
                                  Sequence & c,
                                  TSetType typeSet = DEFAULT_SET);

        static void parseSequenceLocation(std::string & location,
                                          std::string & annotSeqId,
                                          vector<Pair<int> > &,
                                          vector<Pair<int> > &,
                                          vector<Pair<int> > &,
                                          bool &, bool &, int &);

        static void readExonLocations(std::string location,
                                      std::string & annotSeqId,
                                      vector<Pair<int> > & v);

        static void computeRepeatedElems(list<RepeatedElem> & reps, Sequence &c,
                                         int minLen, TStrand strand);

        static void computeCodingDifferential(FeatureEngine *params,
                                              TSetType geneSet, TStrand strand);

        static void computeSignalDifferential(FSM *fsm, FilterEngine *,
                                              TypedFilter<UCHAR> *contexts,
                                              TypedFilter<EdgeInst> **signals,
                                              FeatureEngine *params,
                                              TSetType geneSet, TStrand strand);

        ~GeneUtils() {

        }

    };
}

#endif
