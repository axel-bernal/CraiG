/****************************************************************************
 * SelfLattice.h - part of the lless namespace, a general purpose
 *                 linear semi-markov structure prediction library
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

#ifndef _SELFLATTICE_H_
#define _SELFLATTICE_H_

#include "Lattice.h"

namespace lless {

    /**
     * The Lattice class represents the graph composed of LatticeVertice nodes,
     * connected by LatticeEdge backedges. It is used by viterbi or other
     * similar algorithms for decoding an input sequence.
     ***************************************************************************/

    class SelfLattice : public Lattice {
    public:
    SelfLattice(FilterEngine &fe, FSM &fsm,
                int numPhases,
                Evaluator &ev,
                int maxWordLength = INT_MAX,
                TNodeId3 syncNode3 = NO_NODE_INST,
                TypedFilter<EdgeInst> **signals = NULL,
                TypedFilter<UCHAR> *contexts = NULL) :
        Lattice(fe,
                fsm,
                numPhases,
                ev,
                maxWordLength,
                syncNode3,
                signals,
                contexts, NULL) {

        }

        void viterbiEdge();
        void updModelParams(SeqTags &,
                            FeatureEngine &p,
                            double updVal);

        double dotModelParams(SeqTags &,
                              FeatureEngine &p);

        ~SelfLattice() {}

    };
}

#endif
