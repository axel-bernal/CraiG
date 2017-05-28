/****************************************************************************
 * FSMQueue.h - p part of the lless namespace, a general purpose
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
 ****************************************************************************/

#ifndef _GENE_FSM_QUEUE_
#define _GENE_FSM_QUEUE_

#include "FilterEngine.h"
#include "EdgeQueue.h"
#include "FSM.h"
#include "TypeDefs.h"

#define QITERATOR int

namespace lless {

    /**
     * FSMQueue creates and manages all the EdgeInst queues that are needed for
     * viterbi decoding. It speeds up the decoding time considerably.
     *
     ***************************************************************************/

    typedef pair<TParseEdge, EdgeQueue *> TEdgeQ;

    class FSMQueue {

    protected:
        int _maxQSize;
        FilterEngine *_fe;
        FSM *_fsm;
        TypedFilter<EdgeInst> **_signals;

        /*
         * The next queues are defined per each node->nextEdge(), for all
         * node \in fsm.
         * Notice that if two nodes have the same nextEdge(), naturally they will
         * share one single edge queue, but if they have different length
         * restrictions, then two edge queues will be created, one for each node.
         */
        vector<EdgeQueue *> _edgeQs;

        EdgeQueue ***_phaseDisruptorQs;

        vector<EdgeInst *> _syncEdge;
        vector<TEdgeId2> _nextEdge;

        //! Mapping from nodes to edge queue indexes
        vector<TEdgeQ> *_node2edgeQs;
        vector<Node *> _nodes;

        //! temporary variable to work only the strand used for decoding
        int loStrand;
        //! temporary variable to work only the strand used for decoding
        int hiStrand;

    public:

        FSMQueue(FilterEngine &fe, FSM &fsm, TypedFilter<EdgeInst> **);

        inline vector<TEdgeQ> & edgeQueues(Node *node) {
            return _node2edgeQs[node->id()];
        }

        /**
         * @return true if node has an associated EdgeQueue object, false
         * otherwise. A node has
         * one if it is not a FIXEDLEN_NODE and it is neither a SYNC_BEG_STATE
         * or a SYNC_END_STATE.
         */
        inline bool hasEdgeQ(Node *node) {
            if(node->strand() < loStrand || node->strand() >= hiStrand)
                return false;
            if(node->nodeType() == FIXEDLEN_NODE)
                return false;
            if(node->id() == SYNC_BEG_STATE || node->id() == SYNC_END_STATE)
                return false;

            return true;
        }

        void resetQueues();
        void queue();
        void queueSignal(int idQ, EdgeInst *);
        void dequeue(int pos);
        void uneclipseQs(TParseEdge, EdgeInst *);

        ~FSMQueue();

    };

}

#endif
