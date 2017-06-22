#include "GeneEvaluator.h"
#include "GeneTagPrinter.h"
#include "GeneUtils.h"

/****************************************************************************
 * GeneEvaluator.cpp - part of the craig namespace, a genomics library
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


extern bool verbose;

namespace craig {

    /**
     * Computes the accuracy of predTags taking expTags as the correct set. It
     * converts the Tag objects to Gene objects and then proceeds.
     * @see computePredAccuracy(double *, Sequence &, SeqTags &, SeqTags &)
     */

    /* Legacy
       void GeneEvaluator::computePredAccuracy(Sequence &c,
       SeqTags & expTags,
       SeqTags & predTags) {

       list<Gene> predGenes, expGenes;
       resetPredAccuracy();
       GeneUtils::tags2Genes(c, expTags, "EXP", *_fsm, expGenes);
       GeneUtils::tags2Genes(c, predTags, "PRED", *_fsm, predGenes);
       computePredAccuracy(expGenes, predGenes, 99.9);
       } */


    void GeneEvaluator::computeNodeAccuracy(NodeInst &node) {
        assert(seqLength && node.getLen());

        int pt = node4Loss(node);
        if(pt < 0)  return;

        int fneg = 0;
        int beg = node.getPos() - 1, end = node.getPos() + node.getLen() - 1;
        int *hamming = annot4Loss[0][pt];
        int hmatches = hamming[end] - hamming[beg];

        if(pt) {
            predAccuracy[0] += hmatches;
            predAccuracy[1] += node.getLen() - hmatches;

            int *thamming = annot4Loss[0][pt - 1];
            if(pt % 2)
                thamming = annot4Loss[0][pt + 1];

            fneg = thamming[end] - thamming[beg];
        }
        else {
            predAccuracy[3] += hmatches;
            predAccuracy[2] += node.getLen() - hmatches;
        }

        predAccuracy[3] += node.getLen() - fneg;
        predAccuracy[2] += fneg;

        int *segment = annot4Loss[1][pt];
        int smatches = segment[end] - segment[beg];
        int wrong = 0;

        if(pt) {
            int *tsegment = annot4Loss[1][pt - 1];
            if(hmatches < node.getLen()) {
                predAccuracy[5]++;

                if(pt % 2)
                    tsegment = annot4Loss[1][pt + 1];

                predAccuracy[6] += (tsegment[end] - tsegment[beg])/2;
                if(tsegment[end] < tsegment[beg])
                    cerr << predAccuracy[6] << " " << tsegment[end] << " " << tsegment[beg] << endl;
            }
            else {
                tsegment = annot4Loss[1][0];

                if(smatches == 1 &&
                   (end == seqLength || tsegment[end + 1] == tsegment[end] + 1) &&
                   segment[beg + 1] == segment[beg] + 1)
                    predAccuracy[4]++;
                else
                    predAccuracy[5]++;

            }
            predAccuracy[7]++;

        }
        else {
            predAccuracy[6] += smatches/2;
            if(smatches < 0) {
                cerr << predAccuracy[6] << " " << pt << " " << smatches << " " << end << " " <<  segment[end] << " " <<  beg << " " << segment[beg] << " node " << node.getPos() << " " << node.getLen() << endl;
                for(int i = 0; i <= seqLength; i++)
                    cerr << "[" << i << "]=" << segment[i] << " ";
                cerr << endl;
            }
        }
        //    cerr << "report after " << pt << " " << beg << " " << end << " ";
        //    reportAccuracy((ofstream &)cerr);

    }

    // predAccuracy[8],[9] FP and FN for STARTS, [10],[11] FP and FN for other signals

    void GeneEvaluator::computeEdgeAccuracy(Sequence &c,
                                            SeqTags & expTags,
                                            SeqTags & predTags) {

        SeqTags::iterator it = expTags.begin();
        SeqTags::iterator pit = predTags.begin();
        int diffs = 0;
        int ed, ped;
        //    cerr << "sizes = " << expTags.size() << " " << predTags.size() << endl;

        while(it != expTags.end() && pit != predTags.end()) {
            if((*it)->getGEClass() != EDGE_INST ||
               (*it)->getPos() == 0 || (*it)->getPos() > c.length()) {
                it++;
                continue;
            }
            if((*pit)->getGEClass() != EDGE_INST ||
               (*pit)->getPos() == 0 || (*pit)->getPos() > c.length()) {
                pit++;
                continue;
            }

            //      cerr << (*it)->getPos() << " " << (*it)->getParseType() << "  " << (*pit)->getPos() << " "
            // << (*pit)->getParseType() << " " << diffs << endl;

            ed = edge4Loss((EdgeInst &)*(*it));
            ped = edge4Loss((EdgeInst &)*(*pit));

            if((*it)->getPos() == (*pit)->getPos()) {
                if((*it)->getParseType() != (*pit)->getParseType()) {
                    predAccuracy[ed] += 0.5;
                    predAccuracy[ped + 1] += 0.5;
                    //	  predAccuracy[ed]++;
                    //	  predAccuracy[ped + 1]++;

                }
                it++; pit++;
            }
            else if((*it)->getPos() > (*pit)->getPos()) {
                predAccuracy[ped]++;
                pit++;
            }
            else {
                predAccuracy[ed + 1]++;
                it++;
            }
        }

        while(it != expTags.end()) {
            if((*it)->getGEClass() == EDGE_INST &&
               (*it)->getPos() != 0 && (*it)->getPos() <= c.length()) {
                ed = edge4Loss((EdgeInst &)*(*it));
                predAccuracy[ed + 1]++;
                //	cerr << (*it)->getPos() << " " << (*it)->getParseType() << " " << diffs << endl;
            }
            it++;
        }

        while(pit != predTags.end()) {
            if((*pit)->getGEClass() == EDGE_INST &&
               (*pit)->getPos() != 0 && (*pit)->getPos() <= c.length()) {
                ped = edge4Loss((EdgeInst &)*(*pit));
                predAccuracy[ped]++;
                //	cerr << (*pit)->getPos() << " " << (*pit)->getParseType() << " " << diffs << endl;
            }
            pit++;
        }

        //    cerr << "loss " << diffs << endl;

    }

    /**
     * Loss functions for base level. It is defined as
     * the percentage of correctly predicted nucleotides.
     */
    double GeneEvaluator::corrCoeff() {
        double *pa = predAccuracy;
        double d = (pa[0] + pa[2])*(pa[3] + pa[1])*
            (pa[0] + pa[1])*(pa[3] + pa[2]);

        if(!d)
            return 0;

        return (pa[0]*pa[3] - pa[2]*pa[1])/pow(d, 0.5);

    }


    double GeneEvaluator::loss(Sequence &c,
                               SeqTags &expTags,
                               SeqTags &predTags, int t) {

        double base_result = 0, edge_result = 0;
        computePredAccuracy(c, expTags, predTags, false);
        double lossFactor = (seqLoss[c.length()] - seqLoss[1] + 1)/c.length();
        releaseAnnot4Loss();
        int numDiffEdges = 0;
        double hamming = predAccuracy[1] + predAccuracy[2], lambda;

        switch(currLoss) {
        case LF_HAMMING:
            base_result = hamming; break;
        case LF_CORR_COEF:
            base_result = 1 - corrCoeff(); break;
        case LF_ZERO_ONE:
            base_result = (hamming > 0) ? 1 : 0; break;
        case LF_FSCORE:
            base_result = hamming/(2*predAccuracy[0] + hamming);
        default:break;
        }

        switch(currEdgeLoss) {
        case LF_EDGE:
            for(int i = 8; i < 12; i++) edge_result += predAccuracy[i];
            break;
        case LF_SOFT_EDGE:
            lambda = 1.0/(t^(1/4) + 1);
            edge_result += (1 - lambda)*(predAccuracy[8] + predAccuracy[9]);
            edge_result += lambda*(predAccuracy[10] + predAccuracy[11]);
            break;
        case LF_SEGMENT:
            edge_result = (predAccuracy[5] + predAccuracy[6]);
        default:break;
        }

        edge_result += addStartLoss*(predAccuracy[8] + predAccuracy[9])/100.0;
        //    cerr << result << " " << addStartLoss << "  " << predAccuracy[8] << " " << predAccuracy[9] << " " << lossFactor << " " << edgeLossFactor << endl;

        return lossFactor*(base_result + edgeLossFactor*edge_result);
    }

    /**
     * Computes the accuracy of predGenes taking expGenes as the correct set.
     * TP FP FN TN are obvious. CE WE ME EE PE are shorts for exon accuracy
     * ct is a table [TP FP FN TN CE WE ME EE PE]
     * @param minPercCov is the minimum percentage that expected exon must be
     * covered in order to be declared correct considered match. Default is
     * 99.9.
     */

    /*  Legacy
        void GeneEvaluator::computePredAccuracy(list<Gene> & expGenes,
        list<Gene> & predGenes,
        float minPercCov) {

        int i, j, t;
        expGenes.sort();
        predGenes.sort();

        list<Gene>::iterator itP = predGenes.begin(), itE = expGenes.begin();
        list<Gene *> ppList, peList;
        std::string fmt = "locs";

        char *cP = NULL, *cE = NULL;
        char *lastcP, *lastcE;
        int *annotSeqReg[2];
        int annotSeqLen = 0;
        int toE = 0, toP = 0;
        int toCE = 0, toME = 0, toWE = 0;
        int strand;

        // gene lists are empty?
        if(itP != predGenes.end())
        cP = (char *)itP->getSequence()->id().c_str();

        if(itE != expGenes.end())
        cE = (char *)itE->getSequence()->id().c_str();

        vector<Transcript> *transcripts = NULL;
        vector<Exon> *exons = NULL;

        while(1) {
        annotSeqReg[STRAND_FWD] = NULL; annotSeqReg[STRAND_COMP] = NULL;
        lastcP = cP;
        lastcE = cE;

        // we are in the same annotSeq in both lists

        if(cE && cP && !strcmp(cP, cE)) {
        annotSeqLen = itP->getSequence()->length();

        for(strand = 0; strand < NUM_STRANDS; strand++)
        annotSeqReg[strand] = new int [annotSeqLen + 2];

        //  all of them are true negatives at beginning

        for(i = 0; i <= annotSeqLen + 1; i++)
        for(strand = 0; strand < NUM_STRANDS; strand++)
        annotSeqReg[strand][i] = 3;

        // iterate over the prediction list until we move to the next annotSeq

        while(cP && !strcmp(cP, lastcP)) {
        strand = itP->getStrand();
        transcripts = &itP->transcripts();

        for(t = 0; (unsigned)t < transcripts->size(); t++) {
        exons = &(*transcripts)[t].exons();

        for(i = 0; (unsigned)i < exons->size(); i++) {

        for(j = (*exons)[i].begin(); j <= (*exons)[i].end(); j++)
        annotSeqReg[strand][j] -= 2; // false positive
        }
        toP += exons->size();
        }
        ppList.push_back(&(*itP));
        itP++;
        cP = (itP == predGenes.end()) ? NULL: (char *)itP->getSequence()->id().c_str();
        }

        while(cE && !strcmp(cE, lastcE)) {
        strand = itE->getStrand();
        transcripts = &itE->transcripts();

        for(t = 0; (unsigned)t < transcripts->size(); t++) {
        exons = &(*transcripts)[t].exons();

        for(i = 0; (unsigned)i < exons->size(); i++) {

        for(j = (*exons)[i].begin(); j <= (*exons)[i].end(); j++)
        annotSeqReg[strand][j] -= 1; // true positive or false neg

        }
        }

        peList.push_back(&(*itE));
        itE++;
        cE =  (itE == expGenes.end()) ? NULL: (char *)itE->getSequence()->id().c_str();
        }
        }
        // prediction list's gene is smaller iterate over pred list then
        else if((cE && cP && strcmp(cP, cE) < 0) || (itE == expGenes.end() && cP)) {
        annotSeqLen = itP->getSequence()->length();

        for(strand = 0; strand < NUM_STRANDS; strand++)
        annotSeqReg[strand] = new int [annotSeqLen + 2];

        for(i = 0; i <= annotSeqLen + 1; i++)

        for(strand = 0; strand  <NUM_STRANDS; strand++)
        annotSeqReg[strand][i] = 3; //  all of them are true negatives at beginning

        while(cP && !strcmp(cP, lastcP)) {
        strand = itP->getStrand();
        transcripts = &itP->transcripts();

        for(t = 0; (unsigned)t < transcripts->size(); t++) {
        exons = &(*transcripts)[t].exons();

        for(i = 0; (unsigned)i < exons->size(); i++) {

        for(j = (*exons)[i].begin(); j <= (*exons)[i].end(); j++)
        annotSeqReg[strand][j] -= 2; // false positive
        }
        toP += exons->size();
        }

        ppList.push_back(&(*itP));
        itP++;
        cP = (itP == predGenes.end()) ? NULL: (char *)itP->getSequence()->id().c_str();
        }
        }
        // expected list's gene is smaller. Iterate over expected list then
        else if((cE && cP && strcmp(cP, cE) > 0) || (itP == predGenes.end() && cE)) {
        // cP > cE
        annotSeqLen = itE->getSequence()->length();

        for(strand = 0; strand < NUM_STRANDS; strand++)
        annotSeqReg[strand] = new int [annotSeqLen + 2];

        for(i = 0; i <= annotSeqLen + 1; i++)

        for(strand = 0; strand < NUM_STRANDS; strand++)
        annotSeqReg[strand][i] = 3; //  all of them are true negatives at beginning

        while(cE && !strcmp(cE, lastcE)) {
        strand = itE->getStrand();
        transcripts = &itE->transcripts();

        for(t = 0; (unsigned)t < transcripts->size(); t++) {
        exons = &(*transcripts)[t].exons();

        for(i = 0; (unsigned)i < exons->size(); i++) {

        for(j = (*exons)[i].begin(); j <= (*exons)[i].end(); j++)
        annotSeqReg[strand][j] -= 1; // true positive or false neg
        }
        }
        peList.push_back(&(*itE));
        itE++;
        cE =  (itE == expGenes.end()) ? NULL: (char *)itE->getSequence()->id().c_str();
        }
        }
        else
        break;

        list<Gene *>::iterator pit;

        for(pit = peList.begin(); pit != peList.end(); pit++) {
        strand = (*pit)->getStrand();
        transcripts = &(*pit)->transcripts();

        for(t = 0; (unsigned)t < transcripts->size(); t++) {
        exons = &(*transcripts)[t].exons();

        for(i = 0; (unsigned)i < exons->size(); i++) {
        int numTP = 0;

        for(j = (*exons)[i].begin(); j <= (*exons)[i].end(); j++)

        if(!annotSeqReg[strand][j])
        numTP++;

        if((100.0*numTP)/abs((*exons)[i].end() - (*exons)[i].begin() + 1) >= minPercCov)
        toCE++;
        else
        toME++;
        }
        }
        }

        for(pit = ppList.begin(); pit != ppList.end(); pit++) {
        strand = (*pit)->getStrand();
        transcripts = &(*pit)->transcripts();

        for(t = 0; (unsigned)t < transcripts->size(); t++) {
        exons = &(*transcripts)[t].exons();

        for(i = 0; (unsigned)i < exons->size(); i++) {
        int numTP = 0;

        for(j = (*exons)[i].begin(); j <= (*exons)[i].end(); j++)

        if(!annotSeqReg[strand][j])
        numTP++;

        if((100.0*numTP)/abs((*exons)[i].end() - (*exons)[i].begin() + 1) < minPercCov)
        toWE++;
        }
        }
        }

        for(strand = 0; strand < NUM_STRANDS; strand++) {

        for(i = 1; i <= annotSeqLen; i++) {
        predAccuracy[annotSeqReg[strand][i]]++;
        }

        delete [] annotSeqReg[strand];

        }
        ppList.clear();
        peList.clear();
        }

        predAccuracy[4] += toCE; predAccuracy[5] += toWE;
        predAccuracy[6] += toME;
        predAccuracy[7] += toP;
        return;
        } */

}
