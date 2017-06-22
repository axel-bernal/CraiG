/****************************************************************************
 * Sequence.h - part of the lless namespace, a general purpose
 *              linear semi-markov structure prediction library
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

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <ctype.h>
#include <assert.h>
#include <iostream>
#include <cfloat>
#include "EdgeInst.h"
#include "NodeInst.h"
#include "Sigma.h"
#include "FSM.h"
#include <list>
#include <vector>
#include <string>
#include <map>
#include "FeatureVector.h"

#define NUM_SET_TYPES 5

using namespace lless;

namespace lless {

    typedef enum {DEFAULT_SET, PRED_SET, TRAIN_SET, EXTRA1_SET, EXTRA2_SET, NO_SET} TSetType;

    // forward declarations
    class BioFeature;
    class Feature;

    /**
     * SeqTags is a list of Tag objects with an initial and end phase and a
     * score.
     * The initial phase corresponds to the phase of the firt Tag object in the
     * list, whereas the end phase corresponds to the last Tag object in the
     * list
     */
    class SeqTags : public list<Tag *> {
    private:
        static std::list<Tag *> _allTags;

    protected:
        double _score;
        int _initPhase;
        int _endPhase;
        std::map<Tag *, int, std::less<Tag *> > phaseBreaks;
        int _rank;

    public:
        /**
         * Default constructor. All fields are initialized to zero
         * @param score the score of this sequence tagging
         * @param initPhase phase at the beginning of SeqTag
         * @param endPhase phase at the end of SeqTag
         */
        SeqTags(double score = 0.0, int initPhase = 0,
                int endPhase = 0, int rank = -1) {
            initialize(score, initPhase, endPhase, rank);
        }

        void initialize(double score, int initPhase, int endPhase, int rank) {
            _score = score;
            _initPhase = initPhase;
            _endPhase = endPhase;
            _rank = rank;
        }

        void clone(SeqTags & st) {
            _score = st.score();
            _initPhase = st.initPhase();
            _endPhase = st.endPhase();
            _rank = st.rank();

            SeqTags::iterator it = st.begin();
            for( ; it != st.end(); it++) {
                this->push_back(*it);
                if((*it)->getGEClass() == NODE_INST) {
                    int phaseBreak = st.phaseBreak(*it);
                    if(phaseBreak > 0) // phaseBreak only could happen in between genes
                        this->insertPhaseBreak(this->back(), phaseBreak);
                }
            }

        }

        inline void insertPhaseBreak(Tag *t, int initPhase) {
            phaseBreaks[t] = initPhase;
            //      cerr << "phase break @ " << t->getPos() << " " << initPhase << endl;
        }

        int phaseBreak(Tag *t) {
            std::map<Tag *, int, std::less<Tag *> >::iterator it = phaseBreaks.find(t);
            //      for(; it != phaseBreaks.end(); it++)
            //        cerr << "break "  << it->first->getPos() << " " << it->first->getParseType() << " " << it->first->getStrand() << " " << t->getPos() << " " << t->getParseType()  << " " << t->getStrand() << endl;

            if(it == phaseBreaks.end())
                return -1;

            //      cerr << " found break @ " << it->first->getPos() << " " << it->first->getParseType() << " " << t->getParseType() << " " << t->getPos() << " " << it->second << endl;

            return it->second;
        }

        inline void tpush_back(Tag *tag) {
            _allTags.push_back(tag);
            this->push_back(tag);
        }
        inline void tpush_front(Tag *tag) {
            _allTags.push_back(tag);
            this->push_front(tag);
        }

        bool operator==(const SeqTags & stags) {
            SeqTags &st = (SeqTags &)stags;
            SeqTags::iterator mit = this->begin(), it = st.begin();
            for( ; mit != this->end() && it != st.end(); mit++, it++) {
                if((**it) != (**mit))
                    return false;
            }

            if(mit != this->end() || it != st.end())
                return false;

            return true;
        }

        bool equals(const SeqTags & stags, TNodeId3 syncNode3) {
            SeqTags &st = (SeqTags &)stags;
            SeqTags::iterator mit = this->begin(), it = st.begin();
            for( ; mit != this->end() && it != st.end(); mit++, it++) {
                if((*mit)->getGEClass() == EDGE_INST)
                    continue;

                NodeInst & b = (NodeInst &)*(*it);
                if(b.getType() == syncNode3)
                    continue;

                if(b != (**mit))
                    return false;
            }

            if(mit != this->end() || it != st.end())
                return false;

            return true;
        }

        inline int rank() {  return _rank; }
        inline double score() {  return _score; }
        inline int initPhase() {  return _initPhase; }
        inline int endPhase() {  return _endPhase; }
        inline void setRank(int rank) {  _rank = rank; }
        inline void setScore(double score) {   _score = score;  }
        inline void setInitPhase(int phase) {   _initPhase = phase;  }
        inline void setEndPhase(int phase) {   _endPhase = phase;  }

        NodeInst *buildTag(TNodeId3, TParseNode, int from, int to, TStrand strand);
        EdgeInst *buildTag(TEdgeId2, TParseEdge, int pos, TStrand strand);
        void reverseComplement(FSM &fsm, int numPhases, int seqLength);
        static void releaseMem();

        ~SeqTags() {}

    };

    class BasicSeq {
    protected:
        std::string seqId;
        int seqLen[NUM_STRANDS];
        bool circular;
        bool zero_based;
    public:
        BasicSeq(bool zero_based = true) {
            seqLen[STRAND_FWD] = seqLen[STRAND_COMP] = 0;
            this->zero_based = zero_based;
        }

        BasicSeq(std::string &, bool, bool = true);
        BasicSeq(const BasicSeq &);
        BasicSeq& operator= (const BasicSeq & c);

        inline bool operator==(const BasicSeq & c) {
            return id().compare(((BasicSeq &)c).id()) == 0;
        }

        inline bool operator<(const BasicSeq & c) {
            int b = id().compare(((BasicSeq &)c).id());
            return b <= 0;
        }

        virtual inline BasicSeq *operator[](int seqNo) {
            return this;
        }

        /*
         * Makes sure that beg is greater than zero and end is not greater
         * than the sequence's length. If any of them do not comply with these
         * restrictions, it is 'fixed'.
         */
        inline void fixLocation(int &beg, int &end, TStrand strand = STRAND_FWD) {
            Utils::fixWithinBounds(beg, end, seqLen[strand]);
        }


        inline std::string & id() {
            return seqId;
        }

        bool isCircular() {
            return circular;
        }

        bool isZeroBased() {
            return zero_based;
        }

        inline int length(TStrand strand = STRAND_FWD) {
            return seqLen[strand];
        }

        inline int setLength(int length, TStrand strand = STRAND_FWD) {
            return seqLen[strand] = length;
        }

        inline void setId(std::string &id) {
            seqId = id;
        }

        static BasicSeq* findSequence(std::string &seq, list<BasicSeq *> &);
        virtual BasicSeq *getSubSequence(int beg, int end) = 0;

        virtual ~BasicSeq() {}
    };

    /**
     * The ScoreSeq class represents a input sequence of TClass scores used
     * often to store filter values read from an external file, in conjunction
     * with objects of type FL_InpFile.
     ***************************************************************************/
    template <class TClass>
        class ScoreSeq : public BasicSeq {
    private:
        TClass *seq[NUM_STRANDS];
    public:

        /**
         * Main constructor. Constructs a 1-based ScoreSeq object
         * @param cId unique sequence identifier.
         * @param cSeq the sequence as a stream of TClass objects.
         */
    ScoreSeq(std::string & cId,
             TClass *cSeq = NULL,
             int seqLen = 0) : BasicSeq(cId, false, false) {

            for(int strand = 0; strand < NUM_STRANDS; strand++)
                seq[strand] = NULL;

            if(cSeq && seqLen)
                setSeq(cSeq, seqLen);
        }

    ScoreSeq(std::string & cId,
             int seqLen) : BasicSeq(cId, false, false) {

            for(int strand = 0; strand < NUM_STRANDS; strand++) {
                seq[strand] = new TClass [seqLen + 1];
                for(int i = 0; i <= seqLen; i++)
                    seq[strand][i] = TClass();
                this->seqLen[strand] = seqLen;
            }
        }

        ScoreSeq(const ScoreSeq<TClass> & seq) {
            *this = seq;
        }

        ScoreSeq<TClass> & operator= (const ScoreSeq<TClass> & c) {
            (BasicSeq &)*this = (BasicSeq &)c;
            ScoreSeq<TClass> & s = (ScoreSeq<TClass> &)c;

            for(int i = 0; i < NUM_STRANDS; i++) {
                TStrand strand = (TStrand)i;
                if(s.getSeq(strand)) {
                    this->seq[strand] = new TClass[s.length(strand) + 1];

                    for(int j = 0; j <= s.length(strand); j++)
                        this->seq[strand][j] = s.getSeq(strand)[j];

                }
                else
                    this->seq[strand] = NULL;

            }

            return *this;
        }

        TClass *cloneSeq(TStrand strand) {
            return seq[strand];
        }

        void deleteClone(TClass *clone) {
            ;
        }

        inline TClass & operator()(TStrand strand, int i) {
            return seq[strand][i];
        }

        inline const TClass operator()(TStrand strand, int i) const {
            return seq[strand][i];
        }

        /**
         * Sets the sequence's strand to be equal to cSeq.
         * @param cSeq the sequence
         * @param strand the strand
         */
        void setSeq(TClass *cSeq, int seqLen, TStrand strand = STRAND_FWD) {
            if(!cSeq || seqLen == 0)
                return;

            this->seqLen[strand] = seqLen;
            seq[strand] = cSeq;

            if(this->seqLen[STRAND_FWD] && this->seqLen[STRAND_COMP] &&
               this->seqLen[STRAND_FWD] != this->seqLen[STRAND_COMP]) {
                assert(0);
                throw EXCEPTION( MISSING_ARGUMENT, string("different seq lengths"));
            }
        }

        void setNatComplementSeq() {

            if(this->seqLen[STRAND_FWD] != 0) {
                if(!seq[STRAND_COMP])
                    seq[STRAND_COMP] = new TClass[this->seqLen[STRAND_FWD] + 1];

                for(int i = 1; i <= this->seqLen[STRAND_FWD]; i++)
                    seq[STRAND_COMP][i] = seq[STRAND_FWD][this->seqLen[STRAND_FWD] - i + 1];

                this->seqLen[STRAND_COMP] = this->seqLen[STRAND_FWD];
            }
            else {
                assert(0);
                throw EXCEPTION( MISSING_ARGUMENT, string("no sequence in fwd strand"));
            }
        }


        inline TClass *getSeq(TStrand strand) {
            return seq[strand];
        }

        /*
         * @return a ScoreSeq object defined from beg to end coordinates.
         * No single annotated element is copied over
         */
        BasicSeq *getSubSequence(int beg, int end) {
            assert(end > beg && beg >= 1 && end <= length());

            std::ostringstream ostr;
            ostr << id() << "_" << beg << "_" << end;
            std::string seqId = ostr.str();
            ScoreSeq<TClass> *c = new ScoreSeq<TClass>(seqId);
            int len = end - beg + 1;
            TClass *region[NUM_STRANDS];


            for(int i = 0; i < NUM_STRANDS; i++) {
                TStrand strand = (TStrand)i;

                int offset = beg;
                if(strand == STRAND_COMP)
                    offset = seqLen[strand] - end + 1;

                region[strand] = NULL;

                if(getSeq(strand)) {
                    region[strand] = new TClass[len + 1];
                    region[strand][0] = TClass(0);
                    for(int j = 1; j <= len; j++)
                        region[strand][j] = getSeq(strand)[offset + j - 1];

                    c->setSeq(region[strand], len, strand);
                }
            }

            return c;
        }

        ~ScoreSeq() {
            for(int i = 0; i < NUM_STRANDS; i++)
                delete [] seq[i];
        }
    };

    template <class TClass>
        class MultiScoreSeq : public BasicSeq {
    private:
        int _numSequences;
        ScoreSeq<TClass> **seqObjs;
    public:

        /**
         * Main constructor.
         * @param cId unique sequence identifier.
         * @param numSequences the number of sequences
         */
    MultiScoreSeq(std::string & cId,
		  int numSequences)
        : BasicSeq(cId, false, false) {

            seqObjs = new ScoreSeq<TClass> * [numSequences];
            for(int i = 0; i < numSequences; i++)
                seqObjs[i] = NULL;

            this->_numSequences = numSequences;

        }

        MultiScoreSeq(const MultiScoreSeq<TClass> & seq) {
            *this = seq;
        }

        inline int numSequences() {
            return _numSequences;
        }

        MultiScoreSeq<TClass> & operator= (const MultiScoreSeq<TClass> & c) {
            (BasicSeq &)*this = (BasicSeq &)c;
            MultiScoreSeq<TClass> & s = (MultiScoreSeq<TClass> &)c;

            _numSequences = s.numSequences();
            this->seqObjs = new ScoreSeq<TClass> *[s.numSequences()];

            for(int i = 0; i < s.numSequences(); i++)
                this->seqObjs[i] = new ScoreSeq<TClass>((ScoreSeq<TClass> &)*s[i]);

            return *this;
        }

        void setSeq(int seqNo, ScoreSeq<TClass> *c) {
            seqObjs[seqNo] = c;

            if(!this->seqLen[STRAND_FWD]) {
                this->seqLen[STRAND_FWD] = c->length(STRAND_FWD);
                this->seqLen[STRAND_COMP] = c->length(STRAND_COMP);
            }
            else if(this->seqLen[STRAND_FWD] != c->length(STRAND_FWD) ||
                    this->seqLen[STRAND_COMP] != c->length(STRAND_COMP)) {
                assert(0);
                throw EXCEPTION( MISSING_ARGUMENT, string("different seq lengths"));
            }

        }

        inline void setNatComplementSeq() {
            for(int i = 0; i < numSequences(); i++)
                setNatComplementSeq(i);
        }

        inline void setNatComplementSeq(int seqNo) {
            seqObjs[seqNo]->setNatComplementSeq();
            this->seqLen[STRAND_COMP] = this->seqLen[STRAND_FWD];

        }

        inline BasicSeq *operator[](int seqNo) {
            return seqObjs[seqNo];
        }

        inline TClass & operator()(TStrand strand, int seqNo, int i) {
            return (*seqObjs[seqNo])(strand, i);
        }

        inline const TClass operator()(TStrand strand, int seqNo, int i) const {
            return (*seqObjs[seqNo])(strand, i);
        }

        inline BasicSeq *getSubSequence(int beg, int end) {
            assert(end > beg && beg >= 1 && end <= length());
            std::ostringstream ostr;
            ostr << id() << "_" << beg << "_" << end;
            std::string seqId = ostr.str();
            int len = end - beg + 1;
            MultiScoreSeq<TClass> *c = new MultiScoreSeq<TClass>(seqId, numSequences());

            for(int i = 0; i < numSequences(); i++) {
                BasicSeq *seqObj = getSubSequence(i, beg, end);
                c->setSeq(i, (ScoreSeq<TClass> *)seqObj);
            }
            return c;
        }

        /*
         * @return a ScoreSeq<TClass> object defined from beg to end coordinates.
         * No single annotated element is copied over
         */
        inline BasicSeq *getSubSequence(int seqNo, int beg, int end) {
            return seqObjs[seqNo]->getSubSequence(beg, end);
        }

        ~MultiScoreSeq() {

            int i = 0;
            for( ; i < numSequences(); i++)
                if(seqObjs[i])
                    delete seqObjs[i];

            delete [] seqObjs;

        }

    };


    /**
     * The Sequence class represents a generic input sequence of characters,
     * with a fixed alphabet and annotated list of Tag objects (or SeqTags
     * objects). It also considers the possibility that the sequence may have
     * two strands
     ***************************************************************************/

    class Sequence : public BasicSeq {
    protected:
        std::string seq[NUM_STRANDS];
        Sigma *sigma;
        list<BioFeature *> bioFeats[NUM_SET_TYPES];
        vector< SeqTags >  tags[NUM_SET_TYPES];

    public:
        Sequence(std::string &,
                 char *cSeq = NULL,
                 Sigma * = NULL,
                 bool circ = false);

        /**
         * Copy constructor
         */
        Sequence(const Sequence & seq) {
            *this = seq;
        }

        Sequence() { }

        inline Sigma *alphabet() {
            return sigma;
        }

        inline char *getSeq(TStrand strand) {
            return (char *)seq[strand].c_str();
        }

        char *cloneSeq(TStrand strand) {
            return (char *)seq[strand].c_str();
        }

        void deleteClone(char *clone) {
            ;
        }


        /**
         * @return true if there are annotated BioFeature objects in the list
         * specified by entry typeSet.
         */
        inline bool hasAnnotBioFeats(TSetType typeSet = DEFAULT_SET) {
            return (bioFeats[typeSet].size() != 0);
        }

        inline list<BioFeature *> & getAnnotBioFeats(TSetType typeSet = DEFAULT_SET) {
            return bioFeats[typeSet];
        }

        /**
         * adds BioFeature object bioFeat to the list specified by typeSet.
         */
        inline void addBioFeature(BioFeature *bioFeat,
                                  TSetType typeSet = DEFAULT_SET) {

            bioFeats[typeSet].push_back(bioFeat);
        }

        /**
         * @return true if there are annotated Tag objects in the SeqTags
         * object specified by entry  typeSet.
         */
        inline bool hasTags(TSetType typeSet = DEFAULT_SET) {
            return (tags[typeSet].size() != 0);
        }

        inline vector<SeqTags> & getTags(TSetType typeSet = DEFAULT_SET) {
            return tags[typeSet];
        }

        /**
         * adds SeqTags object tags to the SeqTags object specified by typeSet.
         */
        inline SeqTags & addTags(SeqTags & tags,
                                 TSetType typeSet = DEFAULT_SET) {

            this->tags[typeSet].push_back(tags);
            return this->tags[typeSet].back();
        }

        /**
         * Resets and clears all the existing SeqTag objects
         */
        inline vector<SeqTags> & resetTags(TSetType set) {
            vector<SeqTags> & ges = getTags(set);

            for(unsigned kth = 0; kth < ges.size(); kth++)
                ges[kth].clear();

            ges.erase(ges.begin(), ges.end());
            return ges;
        }

        inline SeqTags & addTags(int labelSet, TSetType set = DEFAULT_SET) {
            if(tags[set].size() < (unsigned)labelSet + 1)
                tags[set].resize(labelSet + 1);
            return this->tags[set][labelSet];
        }

        inline char & operator()(TStrand strand, int i) {
            return seq[strand][i];
        }

        inline const char operator()(TStrand strand, int i) const {
            return seq[strand][i];
        }

        Sequence& operator= (const Sequence & c);
        void setSeq(char *cSeq, TStrand strand = STRAND_FWD);
        void reverseComplement(FSM &, int numPhases, TSetType typeSet = NO_SET);
        void setNatComplementSeq();
        void getLocation(char *region, int upStream,
                         int downStream,
                         int begin, int end,
                         TStrand = NO_STRAND);
        BasicSeq *getSubSequence(int beg, int end);

        void appendTags(TParseNode &bottlNode, FSM &fsm,
                        Sequence &c, TSetType set, int offset, int limit);

        ~Sequence() {

        }
    };

    /**
     * The EdgeAnnotSequence class represents a input sequence with signal
     * annotations. The signals can be  predictions made by an external program
     * for example.
     ***************************************************************************/

    class EdgeAnnotSeq : public BasicSeq {
    private:
        int _numPhases;
        Sigma *_sigma;
        vector<int> *seq[NUM_STRANDS];

    public:
    EdgeAnnotSeq(std::string &id,
                 int numPhases,
                 Sigma *sigma = NULL,
                 bool circ = false)
        : BasicSeq(id, circ, false) {

            _numPhases = numPhases;
            _sigma = sigma;

            for(int t = 0; t < NUM_STRANDS; t++)
                seq[t] = NULL;

        }

        /**
         * Copy constructor
         */
        EdgeAnnotSeq(const EdgeAnnotSeq & seq) {
            *this = seq;
        }

        vector<int>* cloneSeq(TStrand strand) {
            return seq[strand];
        }

        void deleteClone(vector<int> *clone) {
            ;
        }

        inline Sigma *alphabet() {
            return _sigma;
        }

        inline vector<int>* getSeq(TStrand strand) {
            return seq[strand];
        }

        void setNatComplementSeq() {
            return ;
        }

        inline int activePhase(int pos, TStrand strand) {
            return seq[strand][pos].size() > 0 ? seq[strand][pos][0] : 0;
        }

        inline int numPhases() {
            return _numPhases;
        }

        inline void setSeq(vector<int>* seq, int seqLen, TStrand strand) {
            if(!seq)  return;
            this->seqLen[strand] = seqLen;
            this->seq[strand] = seq;

            if(this->seqLen[STRAND_FWD] && this->seqLen[STRAND_COMP] &&
               this->seqLen[STRAND_FWD] != this->seqLen[STRAND_COMP]) {
                assert(0);
                throw EXCEPTION( MISSING_ARGUMENT, string("different seq lengths"));
            }

        }

        pair<int, int> findEdge(pair<int, int> &coords, TStrand);
        EdgeAnnotSeq& operator= (const EdgeAnnotSeq & c);
        BasicSeq *getSubSequence(int beg, int end);

        ~EdgeAnnotSeq() {
            for(int t = 0; t < NUM_STRANDS; t++)
                if(seq[t]) {
                    delete [] seq[t];
                    seq[t] = NULL;
                }
        }
    };


    /**
     * The SparseSeq<TClass> class represents a input sequence of TClass scores
     * used often to store filter values read from an external file, in conjunction
     * with objects of type FL_InpFile.
     ***************************************************************************/
    template <class TClass>
        class SparseSeq : public BasicSeq {
    private:
        SparseVector<TClass> seq[NUM_STRANDS];

    public:

        /**
         * Main constructor.
         * @param cId unique sequence identifier.
         * @param cSeq the sequence as a stream of TClass objects.
         */

        SparseSeq() {}

    SparseSeq(std::string & cId, int length)
        : BasicSeq(cId, false, false) {

            this->seqLen[STRAND_FWD] = this->seqLen[STRAND_COMP] = length;
        }

        SparseSeq(const SparseSeq<TClass> &c) {
            *this = c;
        }

        SparseSeq<TClass> & operator= (const SparseSeq<TClass> & c) {
            (BasicSeq &)*this = (BasicSeq &)c;

            for(int i = 0; i < NUM_STRANDS; i++) {
                TStrand strand = (TStrand)i;
                this->seq[i] = c.getSeq(strand);
            }
        }

        TClass *cloneSeq(TStrand strand) {
            TClass *clone = new TClass [length() + 1];

            for(int i = 0; i <= length(); i++)
                clone[i] = TClass();

            typename DENSE_HASH<ULONG, TClass>::iterator it = seq[strand].begin();

            for( ; it != seq[strand].end(); it++)
                clone[it->first] = it->second;

            return clone;
        }

        void deleteClone(TClass *clone) {
            delete [] clone;
        }

        void setNatComplementSeq() {

            if(this->seqLen[STRAND_FWD] != 0) {
                typename DENSE_HASH<ULONG, TClass>::iterator it = seq[STRAND_FWD].begin();

                for( ; it != seq[STRAND_FWD].end(); it++)
                    seq[STRAND_COMP][this->seqLen[STRAND_FWD] - it->first + 1] = it->second;

                this->seqLen[STRAND_COMP] = this->seqLen[STRAND_FWD];
            }
            else {
                assert(0);
                throw EXCEPTION( MISSING_ARGUMENT, string("no sequence in fwd strand"));
            }
        }

        inline TClass & operator()(TStrand strand, int i) {
            return seq[strand][i];
        }

        inline const TClass operator()(TStrand strand, int i) const {
            return seq[strand][i];
        }

        inline const SparseVector<TClass> &getSeq(TStrand strand) const {
            return seq[strand];
        }

        /*
         * @return a SparseSeq object defined from beg to end coordinates.
         * No single annotated element is copied over
         */
        BasicSeq *getSubSequence(int beg, int end) {
            assert(end > beg && beg >= 1 && end <= length());

            std::ostringstream ostr;
            ostr << id() << "_" << beg << "_" << end;
            std::string seqId = ostr.str();
            int len = end - beg + 1;
            SparseSeq<TClass> *c = new SparseSeq<TClass>(seqId, len);

            const SparseSeq<TClass> &siv = (const SparseSeq<TClass> &)(*this);

            for(int i = 0; i < NUM_STRANDS; i++) {
                TStrand strand = (TStrand)i;
                int offset = beg;

                if(strand == STRAND_COMP)
                    offset = seqLen[strand] - end + 1;

                for(int j = 1; j <= len; j++) {
                    TClass cval = siv(strand, offset + j - 1);
                    if(cval != TClass())
                        (*c)(strand, j) = cval;
                }
            }

            return c;
        }

        ~SparseSeq() {}

    };

    template <class TClass>
        class MultiSparseSeq : public BasicSeq {
    private:
        int _numSequences;
        SparseSeq<TClass> **seqObjs;
    public:

        /**
         * Main constructor.
         * @param cId unique sequence identifier.
         * @param numSequences the number of sequences
         */
    MultiSparseSeq(std::string & cId,
                   int numSequences)
        : BasicSeq(cId, false, false) {

            seqObjs = new SparseSeq<TClass>* [numSequences];
            for(int i = 0; i < numSequences; i++)
                seqObjs[i] = NULL;

            this->_numSequences = numSequences;

        }

        MultiSparseSeq(const MultiSparseSeq<TClass> &c) {
            *this = c;
        }

        MultiSparseSeq & operator= (const MultiSparseSeq<TClass> & c) {
            (BasicSeq &)*this = (BasicSeq &)c;
            MultiSparseSeq<TClass> & s = (MultiSparseSeq<TClass> &)c;

            _numSequences = s.numSequences();
            this->seqObjs = new SparseSeq<TClass> *[s.numSequences()];

            for(int i = 0; i < s.numSequences(); i++)
                this->seqObjs[i] = new SparseSeq<TClass>((SparseSeq<TClass> &)*s[i]);

            return *this;
        }

        inline int numSequences() {
            return _numSequences;
        }

        void setSeq(int seqNo, SparseSeq<TClass> *c) {
            seqObjs[seqNo] = c;

            if(!this->seqLen[STRAND_FWD]) {
                this->seqLen[STRAND_FWD] = c->length(STRAND_FWD);
                this->seqLen[STRAND_COMP] = c->length(STRAND_COMP);
            }
            else if(this->seqLen[STRAND_FWD] != c->length(STRAND_FWD) ||
                    this->seqLen[STRAND_COMP] != c->length(STRAND_COMP)) {
                assert(0);
                throw EXCEPTION( MISSING_ARGUMENT, string("different seq lengths"));
            }


        }

        inline void setNatComplementSeq() {
            for(int i = 0; i < numSequences(); i++)
                setNatComplementSeq(i);
        }

        inline void setNatComplementSeq(int seqNo) {
            seqObjs[seqNo]->setNatComplementSeq();
            this->seqLen[STRAND_COMP] = this->seqLen[STRAND_FWD];

        }

        inline BasicSeq *operator[](int seqNo) {
            return seqObjs[seqNo];
        }

        inline TClass & operator()(TStrand strand, int seqNo, int i) {
            return (*seqObjs[seqNo])(strand, i);
        }

        inline const TClass operator()(TStrand strand, int seqNo, int i) const {
            const SparseSeq<TClass> &siv = (const SparseSeq<TClass> &)(*seqObjs[seqNo]);
            return siv(strand, i);
        }

        inline BasicSeq *getSubSequence(int beg, int end) {
            assert(end > beg && beg >= 1 && end <= length());
            std::ostringstream ostr;
            ostr << id() << "_" << beg << "_" << end;
            std::string seqId = ostr.str();
            int len = end - beg + 1;
            MultiSparseSeq<TClass> *c = new MultiSparseSeq<TClass>(seqId, numSequences());

            for(int i = 0; i < numSequences(); i++) {
                BasicSeq *seqObj = getSubSequence(i, beg, end);
                c->setSeq(i, (SparseSeq<TClass> *)seqObj);
            }

            return c;
        }

        /*
         * @return a SparseSeq object defined from beg to end coordinates.
         * No single annotated element is copied over
         */
        inline BasicSeq *getSubSequence(int seqNo, int beg, int end) {
            return seqObjs[seqNo]->getSubSequence(beg, end);
        }

        ~MultiSparseSeq() {
            int i = 0;
            for( ; i < numSequences(); i++)
                if(seqObjs[i])
                    delete seqObjs[i];

            delete [] seqObjs;
        }

    };
}

namespace std {
    /**
     * Routine for sorting STL containers instantiated with BasicSeq* types
     */
    template <>
        struct greater<BasicSeq *> {
        bool operator()(const BasicSeq* c1, const BasicSeq* c2)      {
            //Defined for ascending sorting
            if(!c1)
                return true;
            if(!c2)
                return false;
            return ((BasicSeq &)(*c1)) < (*c2);
        }
    };


    /**
     * Routine for sorting STL containers instantiated with Sequence* types
     */
    template <>
        struct greater<Sequence *> {
        bool operator()(const Sequence* c1, const Sequence* c2)
        {
            //Defined for ascending sorting
            if(!c1)
                return true;
            if(!c2)
                return false;
            return ((Sequence &)(*c1)) < (*c2);
        }
    };

}

#endif
