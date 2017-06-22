/****************************************************************************
 * Evaluator.h - part of the lless namespace, a general purpose
 *               linear semi-markov structure prediction library
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

#ifndef _EVALUATOR_H_
#define _EVALUATOR_H_

#include "Sequence.h"
#include <list>
#include "FSM.h"


namespace lless {

    /**
     * The Evaluator class is an abstract class which is to be implemented in
     * order to compute the loss function between expected and predicted parses,
     * which are  used in the optimization algorithm.
     ***************************************************************************/

    class Evaluator {
    protected:
        FSM *_fsm;
        int **annot4Loss[2];
        double *seqLoss;
        int noiseLevel;
        int numTags4Loss;
        int seqLength;
        double *predAccuracy;
        int numEntries4Pred;
        int segmentWait;
        TLossType currLoss, currEdgeLoss;

    public:

        Evaluator(FSM &fsm,
                  int numTags4Loss,
                  int numEntries4Pred,
                  TLossType currLoss = LF_HAMMING,
                  TLossType currEdgeLoss = LF_NONE,
                  int noiseLevel = 0);

        inline void resetPredAccuracy() {
            for(int i = 0; i < numEntries4Pred; i++)
                predAccuracy[i] = 0;
        }

        virtual inline void reportAccuracy(ofstream &ost) {
            for(int i = 0; i < numEntries4Pred; i++)
                ost << predAccuracy[i] << " ";
            ost << endl;
        }

        inline void setLossType(TLossType lossType) {
            currLoss = lossType;
        }

        void setAnnot4Loss(Sequence &c, SeqTags &annotList);

        void computePredAccuracy(Sequence &c,
                                 SeqTags & expTags,
                                 SeqTags & predTags,
                                 bool releaseAnnot = true);

        void computePredAccuracy(list<Sequence *> &seqs,
                                 TSetType trainSet,
                                 TSetType predSet);

        void releaseAnnot4Loss();

        virtual  double nodeLoss(NodeInst &node) = 0;
        virtual void computeNodeAccuracy(NodeInst &node) = 0;
        virtual void computeEdgeAccuracy(Sequence &, SeqTags &, SeqTags &) = 0;
        virtual double loss(Sequence &c, SeqTags &expTags, SeqTags &predTags, int t = 0) = 0;

        virtual int edge4Loss(EdgeInst &edge) {
            return edge.getParseType();
        }

        virtual int node4Loss(NodeInst &node) {
            return node.getParseType();
        }

        virtual ~Evaluator() {
            releaseAnnot4Loss();

            for(int l = 0; l < 2; l++) {
                delete [] annot4Loss[l];
                annot4Loss[l] = NULL;
            }

            delete [] predAccuracy;
            predAccuracy = NULL;

        }

    };

}

#endif
