/****************************************************************************
* Lattice.h - part of the lless namespace, a general purpose
*             linear semi-markov structure prediction library
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

#ifndef _LATTICE_H_
#define _LATTICE_H_
#include "EdgeInst.h"
#include "NodeInst.h"
#include "FeatureEngine.h"
#include <algorithm>
#include "FSMQueue.h"
#include "FT_HiddenSeq.h"
#include "LatticeVertex.h"
#include "Evaluator.h"
#include "IntervalSet.h"

namespace lless { 
  
  /**
   * The Lattice class represents the graph composed of LatticeVertice nodes, 
   * connected by LatticeEdge backedges. It is used by viterbi or other 
   * similar algorithms for decoding an input sequence.
   ***************************************************************************/
  
  class Lattice {

   protected:
    // cache variables for decoding
    UCHAR _vcontext;
    double _vcarriedOverScore;
    bool _vresetWordLength;
    int _vkBest;
    FeatureEngine *_vp;
    
    // current state
    double _vscore;
    int _vbeg, _vlen, _vend;
    TParseNode _vpType;
    TNodeId3 _vbType;
    TStrand _vstrand;
    int _vnextPhase;
    Node * _vfsmNode, *_vfsmLastNode;
    Edge *_vfsmEdge, *_vfsmNextEdge;
    LatticeVertex *_vlattVertex;
    EdgeInst *_vsig;  
    
    //! last state
    TParseNode _vlstPType;
    int _vphase;

    // main members
    LatticeVertex ****lvMat; //!< matrix with lattice vertices 
    int numPhases;           //!< vertex number of phases 
    int maxWordLength;       //!< maximum word length allowed, in #nodes
    FSM *fsm;                //!< variable with information about the FSM
    FilterEngine *fe;        //!< filter information
    FSMQueue *fsmQ;          //!< for managing the signal queues

    /*
     * Internal variables 
     */
    Tag initState, endState;
    EdgeInst noSignal;
    double **initNodeScores;

    /*
     * The following are special variables needed by the decoder:
     * en evaluator, two filters,
     * one label and one feature
     */ 
    Evaluator *ev;           /*!<pointer to evaluator object , for loss
                               augmented decoding */
    TNodeId3 syncNode3;      //!< default node label which must be present  between words

    TypedFilter<EdgeInst> **signals; /*! Through this filter the decoder can
                                       access the sequence signals */
    TypedFilter<UCHAR> *contexts;       //!< Filter to access the context values
    FT_HiddenSeq *_vftHidden;    /*! Feature which perform two functions:
                                        compute dotParamV and check if a
                                        disruptor signal is present in the 
                                        hidden sequence    */
  
   public:
    Lattice(FilterEngine &fe, FSM &fsm, 
            int numPhases,
            Evaluator &ev,
            int maxWordLength = INT_MAX,
            TNodeId3 syncNode3 = NO_NODE_INST,
            TypedFilter<EdgeInst> **signals = NULL,
            TypedFilter<UCHAR> *contexts = NULL,
            FT_HiddenSeq *vftHidden = NULL);
    
    inline FSM *getFSM() {
      return fsm;
    }
  
    inline FilterEngine *filterEngine() {
      return fe;
    }

    /**
     * Makes a LatticeVertex at position pos of type pType
     */
    inline void makeLatticeVertex(int pos, 
                                  TParseNode pType, 
                                  int phases, int k) {

      lvMat[pos][pType] = new LatticeVertex * [phases];
      assert(lvMat[pos][pType]);

      for(int phase = 0; phase < phases; phase++) {
        lvMat[pos][pType][phase] = new LatticeVertex(k);
        assert(lvMat[pos][pType][phase]);
      }
    }

    /**
     * Builds a Lattice (an array of LatticeVertex objects) in the current 
     * position pos
     */
    inline void pos2Lattice(int pos) {
      lvMat[pos] = new LatticeVertex ** [fsm->numParseNodes()];
      assert(lvMat[pos]);
      for(int lst = 0; lst < fsm->numParseNodes(); lst++)
        lvMat[pos][lst] = NULL;
    }

    inline void allocateInitScores() {
      initNodeScores = new double * [fsm->numParseNodes()];
      for(int i = 0; i < fsm->numParseNodes(); i++) {
        initNodeScores[i] = new double [numPhases];
        for(int j = 0; j < numPhases; j++) 
          initNodeScores[i][j] = 0;
      }
    }
    
    inline void deleteInitScores() {
      for(int i = 0; i < fsm->numParseNodes(); i++) {
        delete [] initNodeScores[i];
        initNodeScores[i] = NULL;
      }
      
      delete [] initNodeScores;
      initNodeScores = NULL;
    }

    inline void setInitScore(double score, TParseNode node, 
                             int phase = -1) { 
      if(node == INVALID_NODE)
        return;

      if(phase < 0) {
        for(int j = 0; j < numPhases; j++) 
          initNodeScores[node][j] = score;
      }
      else
        initNodeScores[node][phase] = score;

    }

    inline void resetInitScores() {
      // Resetting any initial node score to ZERO
      for(int i = 0; i < fsm->numParseNodes(); i++)
        for(int j = 0; j < numPhases; j++) 
          initNodeScores[i][j] = 0;
    }

    virtual void viterbiEdge();

    void viterbiNode();
    double viterbi(FeatureEngine & p, int kBest);
    void backtrackPath(vector<SeqTags> &, int kBest);
    int findNextBottleneck(TParseNode & bottlNode,
                           double &bottlScore, 
                           int minLenSyncSt);
    int findNextSafeCut(Sequence &,
                        int begin, int chunkSize, 
                        int syncStSize, TSetType set);

    void updWordInstParams(SeqTags &seqTags,
                           FeatureEngine &p, 
                           double updVal);      

    virtual void updModelParams(SeqTags &, 
                                FeatureEngine &p, 
                                double updVal);    

    double dotWordInstParams(SeqTags &seqTags,
			     FeatureEngine &p);
    
    virtual double dotModelParams(SeqTags &, 
				  FeatureEngine &p,
				  IntervalSet<double> * = NULL);
    
    bool isReachable(vector<SeqTags> &, vector<bool> * = NULL);

    void sortById3TagLen(vector<SeqTags> &, 
			 vector<int> &,
			 vector<TNodeId3> &);

    void deleteVars(); 

    virtual ~Lattice() { 
      deleteVars(); 
      deleteInitScores();
      delete fsmQ;
    }
  
  };

}

#endif
