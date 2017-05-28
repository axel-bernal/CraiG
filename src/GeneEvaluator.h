/****************************************************************************
 * GeneEvaluator.h - part of the craig namespace, a genomics library
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


#ifndef _GENE_EVALUATOR_H_
#define _GENE_EVALUATOR_H_
#include "Evaluator.h"
#include "Sequence.h"
#include "Gene.h"

namespace craig {


    /**
     * These class computes the loss function between gene structures which is to
     * be used during training. It inherits from class less::Evaluator all its
     * functionality and provides additional accuracy prediction functions which
     * receive lists of genes instead of the generic lists of node and edge
     * instances. See class Evaluator for more details
     ***************************************************************************/

    class GeneEvaluator : public Evaluator {
    protected:
        double addStartLoss; // to be added as a percentage of the current loss
        double edgeLossFactor;
    public:
    GeneEvaluator(FSM &fsm,
                  TLossType lossType = LF_HAMMING,
		  TLossType edgeLossType = LF_NONE,
                  bool noisyAnnotation = false,
		  double addStartLoss = 0)
        : Evaluator(fsm,
                    5, 12, lossType,
                    edgeLossType, noisyAnnotation) {

            this->addStartLoss = addStartLoss;
            edgeLossFactor = 0;

            if(edgeLossType == LF_EDGE)
                edgeLossFactor = 10;
            else if(edgeLossType == LF_SEGMENT || edgeLossType == LF_SOFT_EDGE)
                edgeLossFactor = 20;
        }

        inline void setEdgeLossFactor(double edgeLossFactor) {
            if(currEdgeLoss != LF_EDGE && currEdgeLoss != LF_SEGMENT &&
               currEdgeLoss != LF_SOFT_EDGE)
                throw EXCEPTION(NOT_SUPPORTED, "lossfactor not supported with loss function choice");

            this->edgeLossFactor = edgeLossFactor;
        }

        inline int edge4Loss(EdgeInst &edge) {
            return 8 + (edge.getType() != START)*2;
        }

        inline int node4Loss(NodeInst &node) {
            if(node.getType() != EXON && node.getType() != UTR)
                return 0;

            int offset = (node.getType() == EXON) ? 1 : 3;
            return offset + (node.getStrand() == STRAND_COMP);
        }

        /**
         * \todo make it work for alternative splicing
         */

        inline double nodeLoss(NodeInst &node) {
            if(!seqLength)
                return 0;

            int pt = node4Loss(node);

            if(pt < 0)  return 0;

            int fneg = 0;
            int fpos = 0;
            int beg = node.getPos() - 1, end = node.getPos() + node.getLen() - 1;
            int *hamming = annot4Loss[0][pt];
            int hmatches = hamming[end] - hamming[beg];
            int base_result = 0, edge_result = 0;
            double lossFactor = (seqLoss[end] - seqLoss[beg] + 1)/node.getLen();

            if(currLoss == LF_HAMMING) {
                if(pt) {
                    fpos = node.getLen() - hmatches;

                    int *thamming = annot4Loss[0][pt - 1];
                    if(pt % 2)
                        thamming = annot4Loss[0][pt + 1];

                    fneg = thamming[end] - thamming[beg];
                }
                else
                    fneg = node.getLen() - hmatches;
                base_result = fpos + fneg;
            }
            else if(currLoss != LF_NONE)
                throw EXCEPTION(NOT_SUPPORTED, "loss function is not supported");

            if(currEdgeLoss == LF_SEGMENT) {
                int *segment = annot4Loss[1][pt];
                int smatches = segment[end] - segment[beg];

                if(pt) {
                    int *tsegment = annot4Loss[1][pt - 1];

                    if(hmatches < node.getLen()) {
                        fpos++;

                        if(pt % 2)
                            tsegment = annot4Loss[1][pt + 1];

                        fneg = (tsegment[end] - tsegment[beg])/2;
                    }
                    else {
                        tsegment = annot4Loss[1][0];

                        if(smatches == 1 &&
                           (end == seqLength || tsegment[end + 1] == tsegment[end] + 1) &&
                           segment[beg + 1] == segment[beg] + 1)
                            ;
                        else
                            fpos++;

                    }
                }
                else
                    fneg = smatches/2;

                edge_result = fpos + fneg;
            }
            else if(currEdgeLoss != LF_NONE)
                throw EXCEPTION(NOT_SUPPORTED, "edge loss function is not supported");

            return lossFactor*(base_result + edgeLossFactor*edge_result);

        }

        void computeNodeAccuracy(NodeInst &node);
        void computeEdgeAccuracy(Sequence &, SeqTags &, SeqTags &);
        double corrCoeff();
        double loss(Sequence &c, SeqTags &expTags, SeqTags &predTags, int t = 0);


        ~GeneEvaluator() {  }
    };

}


#endif
