#include "FSMQueue.h"

/****************************************************************************
 * FSMQueue.cpp - part of the lless namespace, a general purpose
 *                linear semi-markov structure prediction library
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


namespace lless {

    /**
     * Main constructor
     * @param fe a pointer to a FilterEngine object
     * @param fsm a pointer to a FSM object
     * @param signals the filter containing the EdgeInst occurrences in the
     * input sequence.
     */
    FSMQueue::FSMQueue(FilterEngine &fe, FSM &fsm,
                       TypedFilter<EdgeInst> **signals) {

        int i, j, k;
        this->_fe = &fe;
        this->_fsm = &fsm;
        this->_signals = signals;
        Node *node, *nextNode;
        Edge *edge;
        _maxQSize = -1;

        loStrand = 0;
        hiStrand = NUM_STRANDS;

        if(fsm.parsingStrand() != BOTH_STRANDS) {
            loStrand = fsm.parsingStrand();
            hiStrand = loStrand + 1;
        }

        /*
         * compute node2edgeQs
         */
        _node2edgeQs = new vector<TEdgeQ> [fsm.numParseNodes()];

        std::map<std::string, TEdgeQ> cachedNodes;
        _phaseDisruptorQs = new EdgeQueue ** [fsm.numParseEdges()];

        for(i = 0; i < fsm.numParseEdges(); i++)
            _phaseDisruptorQs[i] = NULL;

        for(i = 0; i < fsm.numParseNodes(); i++) {
            node = fsm.node((TParseNode)i);

            if(!hasEdgeQ(node))
                continue;

            vector<Node *> &nextNodes = node->nextNodes();

            for(j = 0 ; (unsigned)j < nextNodes.size() ; j++) {
                nextNode = nextNodes[j];

                if(nextNode->id() == SYNC_END_STATE)
                    continue;

                edge = node->nextEdge(nextNode->id());
                std::ostringstream key;
                key << edge->name() << node->minLength() <<
                    node->maxLength() << node->phaseDisruptor();

                //	cerr << "edge queues for " << node->name() << " " << edge->name() << " " << nextNode->name() << " " << key.str() << endl;

                std::map<std::string, TEdgeQ>::iterator it = cachedNodes.find(key.str());

                if(it != cachedNodes.end()) {
                    TEdgeQ & e = it->second;
                    bool inserted = false;
                    //	  cerr << "key found\n";

                    for(k = 0; k < _node2edgeQs[node->id()].size(); k++)
                        if(e.first == _node2edgeQs[node->id()][k].first) {
                            //	      cerr << "edge parsetype " << e.first << " found\n";
                            inserted = true;
                        }

                    if(!inserted)
                        _node2edgeQs[node->id()].push_back(e);

                    continue;
                }

                //	cerr << "inserting new queue for node " << node->name() << endl;
                EdgeQueue *newQ = new EdgeQueue(node->maxLength());
                newQ->setTailPos(0);
                _edgeQs.push_back(newQ);

                _nodes.push_back(node);
                _node2edgeQs[node->id()].push_back(TEdgeQ(edge->id(), newQ));

                cachedNodes[key.str()] = _node2edgeQs[node->id()].back();

                _nextEdge.push_back(edge->id2());

                TParseEdge disruptor = node->phaseDisruptor();

                if(disruptor != INVALID_EDGE) {
                    int numPhases = node->numPhases();

                    if(!_phaseDisruptorQs[disruptor]) {
                        _phaseDisruptorQs[disruptor] = new EdgeQueue * [numPhases];

                        for(k = 0; k < numPhases; k++)
                            _phaseDisruptorQs[disruptor][k] = new EdgeQueue(node->maxLength());
                    }
                }

                EdgeInst *synch_edge = NULL;

                if(node->isSyncEnd()) {
                    edge = fsm.edge(node->id(), SYNC_END_STATE);
                    synch_edge = new EdgeInst(_nextEdge.back(),
                                              edge->id(),
                                              -1,
                                              edge->strand());
                }

                //	cerr << " " << edge->strand() << "-" << _nextEdge.back() << "<-" << node->id() << endl;
                _syncEdge.push_back(synch_edge);

                if(_maxQSize < node->maxLength())
                    _maxQSize = node->maxLength();

            }
        }
    }

    /**
     * Resets the edge queues and queues the first signals(EdgeInst objects)
     * - up to the length of the node the precedes them -
     * which are read from a signal filter into queues. This function is called
     * when starting decoding of a new sequence.
     */

    void FSMQueue::resetQueues() {
        int i, j;
        int numPhases;
        Node *node;
        TParseEdge disruptor;
        EdgeQueue *sigQ;

        for(i = 0; i < _edgeQs.size(); i++) {
            sigQ = _edgeQs[i];
            node = _nodes[i];
            sigQ->makeEmpty();
            sigQ->setTailPos(0);
            sigQ->setHeadPos(0);
            disruptor = node->phaseDisruptor();

            if(disruptor != INVALID_EDGE) {
                numPhases = node->numPhases();

                for(j = 0; j < numPhases; j++)
                    _phaseDisruptorQs[disruptor][j]->makeEmpty();

            }

            if(_syncEdge[i])
                _syncEdge[i]->setPos(_fe->seqLength() + 1);

        }

        int end = Utils::min(_maxQSize, _fe->seqLength() + 1);

        for(i = 1; i <= end; i++) {
            //      cerr << "position " << i << endl;
            queue();
        }
    }

    /**
     * Queues the signals present at currPos for each edge queue. currPos is
     * incremented before queueing takes place.
     */
    void FSMQueue::queue() {
        int currPos;
        int i, type;
        Node *node;
        EdgeInst *sig = NULL;
        EdgeQueue *sigQ;

        for(i = 0; i < _edgeQs.size(); i++) {
            sigQ = _edgeQs[i];
            node = _nodes[i];
            currPos = sigQ->tailPos() + 1;

            if(currPos - sigQ->headPos() > node->maxLength())
                continue;

            if(currPos > _fe->seqLength() + 1)
                continue;

            sig = &_signals[_nextEdge[i]]->value(currPos, node->strand());

            //      cerr << " signal " << currPos << " " << _nextEdge[i] << " " << sig->getType() << " " << node->name() << endl;

            if(currPos <= _fe->seqLength()) {

                TParseEdge disruptor = node->phaseDisruptor();

                if(disruptor != INVALID_EDGE) {
                    EdgeInst *dsig = &_signals[_fsm->edge(disruptor)->id2()]->value(currPos, node->strand());
                    if(dsig->getType() != NO_EDGE_INST) {
                        int phase = dsig->getPos() % node->numPhases();
                        EdgeInst *last_disruptor = _phaseDisruptorQs[disruptor][phase]->peekTail();

                        if(!last_disruptor || !(*last_disruptor == *dsig)) {
                            _phaseDisruptorQs[disruptor][phase]->queue(dsig);
                            //	      cerr << "disruptor queued " << dsig->getPos() << " " << dsig->getStrand() << " " << phase << " for " << dsig->getType() << endl;
                        }
                    }
                }
            }

            sigQ->setTailPos(currPos);

            if(currPos - sigQ->headPos() < node->minLength())
                continue;

            if(currPos > _fe->seqLength() && _syncEdge[i]) {
                queueSignal(i, _syncEdge[i]);
                continue;
            }

            if(sig->getType() == _nextEdge[i])
                queueSignal(i, sig);
        }
    }

    /**
     * Queues signal sig into the queue identified by idQ.
     * It also eclipses sig if it appears in phase with any disruptor signal
     * which has been queued already.
     */

    void FSMQueue::queueSignal(int idQ, EdgeInst *sig) {
        int numPhases;
        EdgeInst *disruptSig;
        EdgeQueue *sigQ = _edgeQs[idQ];
        TStrand strand = sig->getStrand();
        Node *node = _nodes[idQ];

        sigQ->queue(sig);

        TParseEdge disruptor = node->phaseDisruptor();

        if(disruptor == INVALID_EDGE)
            return;

        sig->uneclipse();
        //    cerr << "uneclipsing and queueing " << sig->getType() << "@" << sig->getStrand() << ":" << sig->getPos() << " for node " << node->name() << endl;

        /*
         * Eclipse signal if it appears in phase with any disruptor signal
         * which has been queued already
         */

        numPhases = node->numPhases();

        for(int i = 0; i < numPhases; i++) {
            disruptSig = _phaseDisruptorQs[disruptor][i]->peek();

            if(!disruptSig)
                continue;

            if(sig->getPos()  > disruptSig->getPos() + 2*(strand == STRAND_FWD)) {
                sig->eclipse(i);
                //	cerr << "eclipsing for node " << node->name() << " " << sig->getPos() << " " << disruptSig->getPos() << " " << disruptSig->getStrand() << " " << i << endl;
            }
        }
    }

    /**
     * Dequeues all signals that appear up to position decodingPos, including
     * disruptor signals.
     */
    void FSMQueue::dequeue(int decodingPos) {
        EdgeInst *sig;
        std::map<TParseEdge, EdgeInst*> disruptSigs;
        EdgeQueue *sigQ;

        int i, j, numPhases;

        for(i = 0; i < _edgeQs.size(); i++) {
            sigQ = _edgeQs[i];
            sig = sigQ->peek();
            Node *node = _nodes[i];

            int pos = decodingPos + node->minLength();

            /*
             * Dequeue variable length/phase node signals
             */
            while(sig && sig->getPos() < pos) {
                sigQ->dequeue();
                //	cerr <<  "dequeueing  " << sig->getType() << " " << sig->getPos() << " " << sig->getStrand() << endl;
                sig = sigQ->peek();
            }

            sigQ->setHeadPos(pos - 1);

            /*
             * Dequeue disruptor signals
             */
            numPhases = node->numPhases();
            TParseEdge disruptor = node->phaseDisruptor();

            if(disruptor == INVALID_EDGE)
                continue;

            for(j = 0; j < numPhases; j++) {
                sigQ = _phaseDisruptorQs[disruptor][j];

                if(!sigQ)
                    continue;

                sig = sigQ->peek();

                while(sig && sig->getPos() <= decodingPos) {
                    disruptSigs[disruptor] = sigQ->dequeue();
                    //	  cerr <<  "dequeueing disruptor " << sig->getType() << " " << sig->getPos() << " " << sig->getStrand() << endl;
                    sig = sigQ->peek();
                }

                sigQ->setHeadPos(decodingPos);
            }
        }

        map<TParseEdge, EdgeInst *>::iterator it = disruptSigs.begin();

        for( ; it != disruptSigs.end(); it++)
            uneclipseQs(it->first, it->second);

    }

    /**
     * Uneclipses all signals which were previously eclipsed when disruptSig
     * was queued.
     * @param disruptor the parsing Type of disruptSig
     * @param disruptSig the disruptor signal which is about to be dequeued.
     */
    void FSMQueue::uneclipseQs(TParseEdge disruptor, EdgeInst *disruptSig) {
        EdgeInst *sig;
        EdgeQueue *sigQ;

        for(int i = 0; i < _edgeQs.size(); i++) {
            Node *node = _nodes[i];
            if(node->phaseDisruptor() != disruptor)
                continue;

            int eclStart = disruptSig->getPos();
            int eclEnd = eclStart + node->maxLength();
            int phase = eclStart % node->numPhases();

            TStrand strand = disruptSig->getStrand();
            /*
             * Shrink eclipsing scope until the occurrence of the next disruptor
             * in the same phase
             */
            if((sig = _phaseDisruptorQs[disruptor][phase]->peek()) != NULL)
                eclEnd = sig->getPos() + 2*(strand == STRAND_FWD);

            sigQ = _edgeQs[i];

            QITERATOR head = sigQ->head();
            sig = sigQ->peek();

            while(sig && (sig->getPos() <= eclEnd)) {
                sig->uneclipse(phase);
                //	cerr << "uneclipsing for node " << node->name() << " " << sig->getType() << " " << sig->getPos() << " " << sig->getStrand() << " " << disruptSig->getPos() << " " << phase << endl;
                sigQ->dequeue();
                sig = sigQ->peek();
            }

            sigQ->restoreHead(head);
        }
    }

    FSMQueue::~FSMQueue() {
        int i, j;
        int numPhases;
        Node *node;

        for(i = 0; i < _edgeQs.size(); i++) {
            if(_edgeQs[i])
                delete _edgeQs[i];
            _edgeQs[i] = NULL;

            if(_syncEdge[i])
                delete _syncEdge[i];
            _syncEdge[i] = NULL;

            node = _nodes[i];
            TParseEdge disruptor = node->phaseDisruptor();

            if(disruptor != INVALID_EDGE && _phaseDisruptorQs[disruptor]) {
                numPhases = node->numPhases();

                for(j = 0; j < numPhases; j++) {
                    delete _phaseDisruptorQs[disruptor][j];
                    _phaseDisruptorQs[disruptor][j] = NULL;
                }

                delete [] _phaseDisruptorQs[disruptor];
                _phaseDisruptorQs[disruptor] = NULL;
            }
        }

        delete [] _phaseDisruptorQs;
        _phaseDisruptorQs = NULL;

        delete [] _node2edgeQs;
        _node2edgeQs = NULL;
    }

}
