/****************************************************************************
 * FT_MeanScoringSegment.h - part of the lless namespace, a general purpose
 *                           linear semi-markov structure prediction library
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

#ifndef FT_MEANSCORING_SEGMENT_FEAT_H
#define FT_MEANSCORING_SEGMENT_FEAT_H

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Sequence.h"
#include "Feature.h"
#include "Filter.h"
#include "FilterEngine.h"
#include "FeatureEngine.h"
#include "FT_Segment.h"

namespace lless {

    /**
     * FT_MeanScoringSegment is a subclass of FT_Segment whose feature value
     * (and its derivate classes) is the score obtained by accumulating
     * values of the contained filter over a tag which must be a segment
     * (node with variable length).
     *
     ***************************************************************************/

    template <class TClass>
        class FT_MeanScoringSegment : public FT_Segment<TClass> {

    protected:
        TypedFilter<int> *starts, *ends;
    public:
    FT_MeanScoringSegment(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass> *filter,
        TypedFilter<int> *starts,
        TypedFilter<int> *ends,
        int period
        )
        : FT_Segment<TClass>(fInd, paramInd,
                             parsingFrames, name,
                             fe, filter, period, 1) {

            this->starts = starts;
            this->ends = ends;
            assert(this->filter->period());
        }

    FT_MeanScoringSegment(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_Segment<TClass>(fInd, paramInd,
                             fargs, offset,
                             fe, 1) {

            this->starts = (TypedFilter<int> *)fe->getFilter(fargs[offset++]);
            this->ends = (TypedFilter<int> *)fe->getFilter(fargs[offset++]);
            assert(this->filter->period());

        }

        virtual inline double featValue(Tag *ge, int frame) {
            TStrand strand = ge->getStrand();
            int lpos = starts->value(ge->getPos(), strand);
            int rpos = ends->value(ge->getPos(), strand);
            int len = rpos - lpos;

            double segVal = this->filter->accValue(rpos - 1, 0, strand)
                - this->filter->accValue(lpos - 1, 0, strand);

            //      cerr << this->getName() << " (" << lpos << "," << rpos << ") "
            // << strand <<" " << (len ? segVal/len : 0) << endl;

            return len ?  segVal/len : 0;

        }

        inline double dotParamV(double featVal,
                                const FeatureVector **params,
                                int fConjInd = 0) {

            if(featVal == DOUBLE_INFINITY)
                return 0;
            //      cerr << "return " << featVal << " " << (*params[fConjInd])[featVal > mean ? 1 : 0] << " " << mean << endl;
            return featVal*(*params[fConjInd])[0];
        }

        inline void updParamV(double updVal, double featVal,
                              FeatureVector **params, int fConjInd = 0) {
            //      cerr << this->getName() << "'s val " << featVal << endl;
            if(featVal == DOUBLE_INFINITY)
                return;

            (*params[fConjInd])[0] += featVal*updVal;
        }

        ~FT_MeanScoringSegment() {

        }

    };


    template <class TClass>
        class FT_MeanScoringSignal : public FT_Segment<TClass> {
    protected:
        FSM *fsm;
        bool left;
    public:
    FT_MeanScoringSignal(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass> *filter,
        int period,
        FSM *fsm,
        bool left
        )
        : FT_Segment<TClass>(fInd, paramInd,
                             parsingFrames, name,
                             fe, filter, period, 1) {

            assert(this->filter->period());
            this->fsm = fsm;
            this->left = left;

        }

    FT_MeanScoringSignal(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_Segment<TClass>(fInd, paramInd,
                             fargs, offset,
                             fe, 1) {

            assert(this->filter->period());
            this->fsm = fte->getFSM();
            left = false;
            if(!fargs[offset++].compare("left"))
                left = true;

        }

        virtual inline double featValue(Tag *tag, int frame) {
            Edge *edge = fsm->edge((TParseEdge)tag->getParseType());
            EdgeInst *ei = (EdgeInst *)tag;
            TStrand strand = tag->getStrand();

            if(strand == STRAND_COMP)
                edge = edge->complementEdge();

            int len = ei->lenNextNode();
            int lpos = ei->getPos() + edge->nextNodePos();
            int rpos = lpos + ei->lenNextNode() - 1;

            if(left) {
                len = ei->lenPrevNode();
                rpos = ei->getPos() + edge->nextNodePos() - 1;
                lpos = rpos - ei->lenPrevNode() + 1;
            }

            if(lpos < 1) lpos = 1;
            if(rpos > this->fe->seqLength() + 1) rpos = this->fe->seqLength() + 1;

            double segVal = this->filter->accValue(rpos - 1, 0, strand)
                - this->filter->accValue(lpos - 1, 0, strand);

            /*      cerr << "MSSig " << ei->lenPrevNode() << " " << ei->lenNextNode() <<
                    " " << this->getName() << " " << ei->getPos() << " " <<
                    ei->getLen() << " " << ei->getParseType() << " (" << lpos <<
                    "," << rpos << ") " << strand <<" " <<
                    this->filter->accValue(lpos - 1, 0, strand) << " " <<
                    this->filter->accValue(rpos - 1, 0, strand) <<  " " <<
                    (len ? segVal/len : 0) << endl; */

            return len ?  segVal/len : 0;

        }

        inline double dotParamV(double featVal,
                                const FeatureVector **params,
                                int fConjInd = 0) {

            if(featVal == DOUBLE_INFINITY)
                return 0;
            //      cerr << "return " << featVal << " " << (*params[fConjInd])[0] << endl;
            return featVal*(*params[fConjInd])[0];
        }

        inline void updParamV(double updVal, double featVal,
                              FeatureVector **params, int fConjInd = 0) {
            //      cerr << this->getName() << "'s val " << featVal << endl;
            if(featVal == DOUBLE_INFINITY)
                return;
            //      cerr << "updating " << this->getName() << endl;
            (*params[fConjInd])[0] += featVal*updVal;
        }

        ~FT_MeanScoringSignal() {

        }

    };

    template <class TClass>
        class FT_RampScoringSignal : public FT_Segment<TClass> {
    protected:
        FSM *fsm;
        bool onset;
    public:
    FT_RampScoringSignal(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass> *filter,
        int period,
        FSM *fsm,
        bool onset
        )
        : FT_Segment<TClass>(fInd, paramInd,
                             parsingFrames, name,
                             fe, filter, period, 1) {

            assert(this->filter->period());
            this->fsm = fsm;
            this->onset = onset;

        }

    FT_RampScoringSignal(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_Segment<TClass>(fInd, paramInd,
                             fargs, offset,
                             fe, 1) {

            assert(this->filter->period());
            this->fsm = fte->getFSM();
            this->onset = Utils::stringToBoolean(fargs[offset++]);

        }

        virtual inline double featValue(Tag *tag, int frame) {
            Edge *edge = fsm->edge((TParseEdge)tag->getParseType());
            EdgeInst *ei = (EdgeInst *)tag;
            TStrand strand = tag->getStrand();

            if(strand == STRAND_COMP)
                edge = edge->complementEdge();

            int rlen = ei->lenNextNode();
            int rposl = ei->getPos() + edge->nextNodePos();
            int rposr = rposl + rlen - 1;

            if(rposl < 1) rposl = 1;
            if(rposr > this->fe->seqLength() + 1) rposr = this->fe->seqLength() + 1;

            double rsegVal = this->filter->accValue(rposr - 1, 0, strand)
                - this->filter->accValue(rposl - 1, 0, strand);

            int llen = ei->lenPrevNode();
            int lposr = ei->getPos() + edge->nextNodePos() - 1;
            int lposl = lposr - llen + 1;

            if(lposl < 1) lposl = 1;
            if(lposr > this->fe->seqLength() + 1) lposr = this->fe->seqLength() + 1;

            double lsegVal = this->filter->accValue(lposr - 1, 0, strand)
                - this->filter->accValue(lposl - 1, 0, strand);

            /*      cerr << "MSSig " << llen << " " << rlen <<
                    " " << this->getName() << " " << ei->getPos() << " " <<
                    ei->getLen() << " " << ei->getParseType() << " (" << lposl <<
                    "-" << lposr << "," << rposl << "-" << rposr << ") " << strand <<
                    " " << lsegVal/llen << " " << rsegVal/rlen  <<  " " << endl; */

            return onset ? ((rlen ? rsegVal/rlen : 0) - (llen ? lsegVal/llen : 0)) :
                ((llen ? lsegVal/llen : 0) - (rlen ? rsegVal/rlen : 0));

        }

        inline double dotParamV(double featVal,
                                const FeatureVector **params,
                                int fConjInd = 0) {

            if(featVal == DOUBLE_INFINITY)
                return 0;
            //      cerr << "return " << featVal << " " << (*params[fConjInd])[0] << endl;
            return featVal*(*params[fConjInd])[0];
        }

        inline void updParamV(double updVal, double featVal,
                              FeatureVector **params, int fConjInd = 0) {
            //      cerr << this->getName() << "'s val " << featVal << endl;
            if(featVal == DOUBLE_INFINITY)
                return;
            //      cerr << "updating " << this->getName() << endl;
            (*params[fConjInd])[0] += featVal*updVal;
        }

        ~FT_RampScoringSignal() {

        }

    };


    /**
     * FT_StDevScoringSegment is a subclass of FT_Segment whose value of this
     * feature (and its derivate classes) is the score obtained by accumulating
     * values of the contained filter over a tag which must be a segment
     * (node with variable length).
     *
     ***************************************************************************/

    template <class TClass>
        class FT_StDevScoringSegment : public FT_Segment<TClass> {

    protected:
        TypedFilter<TClass> *filter2;
    public:
    FT_StDevScoringSegment(
        int fInd,
        int paramInd,
        int parsingFrames,
        char *name,
        FilterEngine *fe,
        TypedFilter<TClass> *filter,
        TypedFilter<TClass> *filter2,
        int period
        )
        : FT_Segment<TClass>(fInd, paramInd,
                             parsingFrames, name,
                             fe, filter, period, 1) {

            this->filter2 = filter2;
            assert(this->filter->period() == this->filter2->period());
        }

    FT_StDevScoringSegment(
        int fInd,
        int paramInd,
        vector<std::string> & fargs,
        int & offset,
        FilterEngine *fe,
        FeatureEngine *fte
        )
        : FT_Segment<TClass>(fInd, paramInd,
                             fargs, offset,
                             fe, 1) {

            this->filter2 = (TypedFilter<TClass> *)fe->getFilter(fargs[offset++]);
            assert(this->filter->period() == this->filter2->period());
        }

        virtual inline double featValue(Tag *ge, int frame) {
            int beg = ge->getPos(), end = beg + ge->getLen() - 1;
            Utils::fixWithinBounds(beg, end, this->fe->seqLength());

            if(end < beg)  return DOUBLE_INFINITY;
            if(end == beg)  return 0;

            double s0 = end - beg + 1;
            beg = this->begTag(ge, frame);
            end = this->endTag(ge, frame);

            double s1 = this->filter->accValue(end, 0, ge->getStrand())
                - this->filter->accValue(beg - this->period(), 0, ge->getStrand());
            double s2 = this->filter2->accValue(end, 0, ge->getStrand())
                - this->filter2->accValue(beg - this->period(), 0, ge->getStrand());

            double fVal =  sqrt(fabs((1/(s0 - 1))*(s2 - (s1*s1)/s0)));
            //      cerr << this->getName() << " (" << beg << " " << end << ") " << s0
            //      	   << " " << s1 << " " << s2 << " " << fVal << endl;

            return fVal;

        }

        inline double dotParamV(double featVal,
                                const FeatureVector **params,
                                int fConjInd = 0) {

            if(featVal == DOUBLE_INFINITY)
                return 0;
            //      cerr << "return " << featVal << " " << (*params[fConjInd])[0] << endl;
            return featVal*(*params[fConjInd])[0];
        }

        inline void updParamV(double updVal, double featVal,
                              FeatureVector **params, int fConjInd = 0) {
            //      cerr << this->getName() << "'s val " << featVal << endl;
            if(featVal == DOUBLE_INFINITY)
                return;

            (*params[fConjInd])[0] += featVal*updVal;
        }

        ~FT_StDevScoringSegment() {

        }

    };

}

#endif
