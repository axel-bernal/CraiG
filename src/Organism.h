/****************************************************************************
 * Organism.h - part of the craig namespace, a genomics library
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

#ifndef _ORGANISM_H
#define _ORGANISM_H

#include "IMM.h"
#include "Utils.h"
#include "EdgeInst.h"
#include "FeatureEngine.h"
#include <list>
#include <vector>
#include "Sequence.h"
#include "Sim.h"
#include "Gene.h"
#include "FilterEngine.h"
#include "Lattice.h"
#include "FSM.h"
#include <assert.h>
#include "NGram.h"

#define MAX_STARTCODONS 4

namespace craig {

    typedef list<Gene> TGeneList;


    /**
     * The Organism class is a legacy class. It contains routines that should be
     * moved to either SequenceUtils.cpp or GeneUtils.cpp. I'll clean them up
     * eventually
     ***************************************************************************/

    class Organism {
    private:
        Sigma *sigma;
        TGeneList genes[NUM_SET_TYPES]; //!< array of lists of genes
        list<Overlap> listOverlaps; //!< track overlaps between gene
        list<Sequence *> annotSeqs; //!< annotSeqs loaded for the organism
        vector<Sim> listSims; //!< list of similarities, for each gene .. hmmm need to change this
        double countstartCodons[MAX_STARTCODONS]; //!< start codon statistics
        vector<unsigned int> countgenomicBases; //!< base frequency statistics
        int lenRBS;
        int windowRBS;
        double *weightRBS[DNA_ALPHABET_SIZE];
        double *vectorRBS;
    protected:
        Sequence *findSequence(char *annotSeq);
    public:
        Organism(Sigma &sigma);

        void setDefaults();

        inline list<Sequence *>& getSequences() {
            return annotSeqs;
        }

        inline void setSequences(list<Sequence *> &annotSeqs) {
            this->annotSeqs = annotSeqs;
        }

        inline void setGenesList(TGeneList& list, TSetType typeSet = DEFAULT_SET) {
            genes[typeSet] = list;
        }

        void getcodonTables(double ***im, int numPos, int minPlets, int maxPlets, int bC, TGene cond, TSetType typeSet = DEFAULT_SET);

        void genesInOneStrand(TSetType = DEFAULT_SET);
        void attachGenes2Sequences(TSetType typeSet = DEFAULT_SET);
        void filterGenesByLength(int minGeneLen = 0, TSetType typeSet = DEFAULT_SET);
        void getannotSeqGaps(int mingapLen = 0, TSetType = DEFAULT_SET, TSetType gapSet = PRED_SET);
        //    void processSIMS(::ifstream *, TSetType typeSet = DEFAULT_SET);
        void buildRBSVector(int numGenes, int *offsetRBS);
        void buildRBSweightMatrix(char **upsRegs, int *offsetRBS, int numGenes, double ** weightFreq);
        void refineRBSweightMatrix(char **upsRegs, int *offsetRBS, int numGenes, double **weightFreq, int percentage);

        void computestartFRQ(TGene, TSetType typeSet = DEFAULT_SET);
        void computegenomFreqs();
        void computeRBSprofile(int lenSD, TGene cond, int windowLen, int percentage, TSetType typeSet = DEFAULT_SET);

        void storegenomFreqs(::ofstream *);
        void storestartFRQ(::ofstream *);
        void storeRBSprofile(::ofstream *);
        void storeStats(::ofstream *);
        void printStats();
        void retrievestartFRQ(::ifstream *);
        void retrievegenomFreqs(::ifstream *);
        void retrieveRBSprofile(::ifstream *);
        void retrieveStats(::ifstream *);

        void pickcorrectStart(int *typeStart, int minGeneLen, int numPos = 3, TSetType typeSet = DEFAULT_SET);

        void changeStateGenes(TGene state, TSetType typeSet = DEFAULT_SET);
        void freeStartInformation(TGene state, TSetType typeSet = DEFAULT_SET);
        void allocateStartInformation(TGene state, TSetType typeSet = DEFAULT_SET);
        void findProkGenes(IMM *, int, double, int, int, int, int, int, int, int, bool);

        inline TGeneList & getGenes(TSetType typeSet) {
            return genes[typeSet];
        }

        ~Organism();

    };
}
#endif
