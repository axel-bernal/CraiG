#include "Evaluator.h"

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
 ****************************************************************************/

namespace lless {


    Evaluator::Evaluator(FSM &fsm,
                         int numTags4Loss,
                         int numEntries4Pred,
                         TLossType currLoss,
                         TLossType currEdgeLoss,
                         int noiseLevel) {



        this->_fsm = &fsm;
        this->noiseLevel = noiseLevel;
        this->seqLoss = NULL;
        this->currLoss = currLoss;
        this->currEdgeLoss = currEdgeLoss;
        if(!numTags4Loss)
            numTags4Loss = fsm.numParseNodes();
        this->numTags4Loss = numTags4Loss;

        for(int l = 0; l < 2; l++) {
            annot4Loss[l] = new int* [numTags4Loss];

            for(int i = 0; i < numTags4Loss; i++)
                annot4Loss[l][i] = NULL;
        }

        this->numEntries4Pred = numEntries4Pred;
        predAccuracy = new double [numEntries4Pred];
        resetPredAccuracy();

        seqLength = 0;
    }


    /**
     * \todo make it work for alternative splicing and also
     * specify that synchroniation state has always Tag4Loss
     * value equal to zero
     */
    void Evaluator::setAnnot4Loss(Sequence &c,
                                  SeqTags &annotList) {

        assert(annotList.size());

        int i, j;
        seqLength = c.length();

        for(int l = 0; l < 2; l++) {
            for(j = 0; j < numTags4Loss; j++) {
                annot4Loss[l][j] = new int [seqLength + 1];
                for(i = 0; i <= seqLength; i++)
                    annot4Loss[l][j][i] = 0;
            }
        }

        seqLoss = new double [seqLength + 1];
        for(i = 0; i <= seqLength; i++)
            seqLoss[i] = 1;

        SeqTags::iterator it = annotList.begin();
        EdgeInst *edge = NULL;

        for( ; it != annotList.end(); it++) {

            if((*it)->getGEClass() != NODE_INST) {
                edge = (EdgeInst *)(*it);
                continue;
            }

            NodeInst &node = (NodeInst &)*(*it);

            int beg = node.getPos(), end = node.getPos() + node.getLen() - 1;

            assert(edge && beg > 0 && end <= seqLength);

            int pt = node4Loss(node);

            if(pt < 0)  continue;

            if(noiseLevel && !pt) {
                int mybeg = (beg == 1) ? beg : beg + 200;
                int myend = (end == seqLength) ? end : end - 200;
                for(i = mybeg; i < myend; i++)
                    seqLoss[i] = 1 - noiseLevel/10;
            }

            for(j = 0; j < numTags4Loss; j++) {
                int *hamming = annot4Loss[0][j];
                int *segment = annot4Loss[1][j];
                segment[beg] = segment[beg - 1] + (j == 0) + (j == pt)*(j != 0);
                hamming[beg] = hamming[beg - 1] + (j == pt);

                for(i = beg + 1; i <= end; i++) {
                    hamming[i] = hamming[i - 1] + (j == pt);
                    segment[i] = segment[i - 1];
                }
            }
        }

        for(i = 1; i <= seqLength; i++)
            seqLoss[i] += seqLoss[i - 1];

    }


    void Evaluator::computePredAccuracy(Sequence &c,
                                        SeqTags & expTags,
                                        SeqTags & predTags,
                                        bool releaseAnnot) {

        SeqTags::iterator it = predTags.begin();
        setAnnot4Loss(c, expTags);
        resetPredAccuracy();
        EdgeInst *edge = NULL;

        for( ; it != predTags.end(); it++) {

            if((*it)->getGEClass() != NODE_INST) {
                edge = (EdgeInst *)(*it);
                continue;
            }

            NodeInst &node = (NodeInst &)*(*it);
            assert(edge);
            computeNodeAccuracy(node);
        }

        if(releaseAnnot)
            releaseAnnot4Loss();

        computeEdgeAccuracy(c, expTags, predTags);

    }


    void Evaluator::computePredAccuracy(list<Sequence *> &seqs,
                                        TSetType trainSet,
                                        TSetType predSet) {

        int k;
        double *accTable1 = new double [numEntries4Pred];
        double *accTable2 = new double [numEntries4Pred];

        for(k = 0; k < numEntries4Pred; k++)
            accTable1[k] = 0;

        list<Sequence *>::iterator cit = seqs.begin();
        for( ; cit != seqs.end(); cit++) {
            Sequence &c = *(*cit);
            double minLoss = DBL_MAX, myLoss;

            // check if prediction was performed in this sequence
            if(!c.getTags(predSet).size())
                continue;

            for(unsigned int i = 0; i < c.getTags(trainSet).size(); i++) {
                for(unsigned int j = 0; j < c.getTags(predSet).size(); j++) {
                    myLoss = loss(c, c.getTags(trainSet)[i], c.getTags(predSet)[j]);

                    if(myLoss < minLoss) {
                        minLoss = myLoss;
                        for(k = 0; k < numEntries4Pred; k++)
                            accTable2[k] = predAccuracy[k];

                    }
                }
            }

            for(k = 0; k < numEntries4Pred; k++)
                accTable1[k] += accTable2[k];
        }

        for(k = 0; k < numEntries4Pred; k++)
            predAccuracy[k] = accTable1[k];

        delete [] accTable1;
        delete [] accTable2;
    }

    void Evaluator::releaseAnnot4Loss() {
        for(int l = 0; l < 2; l++) {
            if(annot4Loss[l])
                for(int i = 0; i < numTags4Loss; i++) {
                    if(annot4Loss[l][i])
                        delete [] annot4Loss[l][i];
                    annot4Loss[l][i] = NULL;
                }
        }

        if(seqLoss) {
            delete [] seqLoss;
            seqLoss = NULL;
        }

        seqLength = 0;

    }

}
