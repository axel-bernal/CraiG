#include "Utils.h"
#include "Lattice.h"

/****************************************************************************
* Lattice.cpp - part of the lless namespace, a general purpose
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

extern bool verbose;

namespace lless {
  

  /**
   * Main constructor
   * @param fe a reference to a FilterEngine object
   * @param fsm a reference to a FSM object
   * @param numPhases maximum number of phases of any Tag object present in 
   * the sequence.
   * @param ev Evaluator object, for computing the loss-augmented decoding
   * is computed.
   * @param maxWordLength maximum number of consecutive nodes between
   * syncNode3 occurrences
   * @param syncNode3 id of the synchronization node. If not provided
   * NO_NODE_INST will be assigned to it.
   * @param signals the filter for accessing EdgeInst objects at each
   * position of the sequence.
   * @param contexts the filter for accesing the Context levels at each 
   * position of the sequence.
   * @param vftHidden the feature that turns on if a hidden sequence is 
   * found during decoding.
   */
  Lattice::Lattice(FilterEngine &fe, FSM &fsm, 
                   int numPhases,
                   Evaluator &ev,
                   int maxWordLength,
                   TNodeId3 syncNode3,
                   TypedFilter<EdgeInst> **signals,
                   TypedFilter<UCHAR> *contexts,
                   FT_HiddenSeq *vftHidden) {
    
    this->fe = &fe;
    this->fsm = &fsm;
    this->numPhases = numPhases;
    this->ev = &ev;
    this->fsmQ = new FSMQueue(fe, fsm, signals);
    
    this->maxWordLength = maxWordLength;
    this->syncNode3 = syncNode3;
    this->signals = signals;
    this->contexts = contexts;
    this->_vftHidden = vftHidden;
    
    lvMat = NULL;
    initState = Tag(SYNC_BEG_STATE, 0, STRAND_FWD);
    endState = Tag(SYNC_END_STATE, 0, STRAND_FWD);
    
    allocateInitScores();
  }
  
  /**
   * Performs viterbi at the edge level
   */
  void Lattice::viterbiEdge() {
                             
    int seqLen = fe->getSequence()->length();
    double edgeScore, nodeScore;
    int wordLength = 0;
    Tag *ge;
    Edge *fsmEdge2;
    LatticeEdge *lattEdge;
    LatticeVertex & lattVertex = (*_vlattVertex);
    int sigPos;
    //    cerr << "size of backedges struct " << lattVertex.size() << endl;
    for(int kth = 0; kth < lattVertex.size(); kth++) {    
      lattEdge = lattVertex[kth];
      /*
       * ge is the node associated with this kth back edge 
       */
      ge = &lattEdge->st;
      _vscore = _vcarriedOverScore + lattEdge->score();
      wordLength = lattEdge->wordLength;
      fsmEdge2 = _vfsmEdge;
      NodeInst nb(_vbType, _vpType, _vbeg, _vstrand, _vlen); 
      
      if(_vfsmEdge->toSelf()) {
        nb.setPos(ge->getPos());
        nb.setLen(_vlen + _vbeg - ge->getPos());
	
        if(nb.getLen() > _vfsmNode->maxLength(_vcontext)) {
	  //	  cerr << "self too long " << nb.getLen() << endl;
          continue;
        }
        lattEdge = lattEdge->backtrack(lvMat[nb.getPos() - 1]);
        ge = &lattEdge->st;

	//	cerr << _vfsmLastNode->id() << " " << ge->getPos()  << " " << _vfsmNode->id() << " " <<nb.getPos() << " " << nb.getLen() << endl;
        _vscore = lattEdge->score();

        fsmEdge2 = fsm->edge((TParseNode)ge->getParseType(), 
                             (TParseNode)nb.getParseType());

        sigPos = nb.getPos() - fsmEdge2->nextNodePos();
      }
      else {
        if(_vfsmNode->phaseDisruptor() != INVALID_EDGE) { 
          /*
           * Dealing with disruptor signals which appear in hidden sequences
           */
          if(_vphase && _vftHidden) {
            if(ge->getStrand() == STRAND_COMP)
              _vftHidden->computeHiddenSeq(ge, 
                                           seqLen - nb.getPos() + 2,
                                           seqLen - ge->getPos() + 1,
                                           (numPhases - _vphase) % numPhases);
            else
              _vftHidden->computeHiddenSeq(ge, 
                                           ge->getPos(), 
                                           nb.getPos() - 1,
                                           _vphase);
            
            if(_vftHidden->hiddenSeqHasStop()) {
	      //	      cerr << "has hidden stop\n";
              continue;
	    }
          }

          wordLength++;
        }

        if(wordLength > maxWordLength)
          continue;

        /*
         * updating and checking the current word length
         */
        if(_vresetWordLength)
          wordLength = 0;

        sigPos = nb.getPos() + _vfsmLastNode->nextEdgePos(_vpType);
      }

      /*
      * getting signal info
      */
      if(fsmEdge2->hasSignal()) {
        _vsig = &signals[fsmEdge2->id2()]->value(sigPos, fsmEdge2->strand());
        _vsig->setParseType(fsmEdge2->id());
      }
      else {
        noSignal.reInitialize(fsmEdge2->id(), sigPos, fsmEdge2->strand());
        _vsig = &noSignal;
      }

      //      cerr << "Sig\t" << sigPos << " " << _vpType << " " << nb.getPos() << " " << nb.getLen() << " " << nb.getStrand() << endl;
      
      /*
       * computing edge score
       */
      //      if(_vbeg == 2343 && nb.getLen() == 199 && _vsig->getStrand() == STRAND_FWD)
      //      	cerr << "hey\n";

      if(fsmEdge2->strand() == STRAND_COMP) {
	_vsig->setPos(seqLen - _vsig->getPos() + 1);
	_vsig->setAssocNodeLens(nb.getLen(), ge->getLen());
	edgeScore = _vp->dotParamV4Edges(_vsig, (numPhases - _vphase) % numPhases);
	_vsig->setPos(seqLen - _vsig->getPos() + 1);
      }
      else {
	_vsig->setAssocNodeLens(ge->getLen(), nb.getLen());
	edgeScore = _vp->dotParamV4Edges(_vsig, _vphase);
      }

      _vscore += edgeScore;
      
      /*
       * computing node score
       */
      if(nb.getStrand() == STRAND_COMP) { 
        nb.setPos(seqLen - nb.getPos() - nb.getLen() + 2);
        nodeScore = _vp->dotParamV4Nodes(&nb,
                                         (numPhases - _vnextPhase) % numPhases);
        nb.setPos(seqLen - nb.getPos() - nb.getLen() + 2);
      }
      else
        nodeScore = _vp->dotParamV4Nodes(&nb, _vphase);

      _vscore += nodeScore;

      /*
       * Add loss-augmented term to the node score
       */
      _vscore += ev->nodeLoss(nb);
      
      if(lvMat[_vend] == NULL) 
        pos2Lattice(_vend);

      if(lvMat[_vend][_vpType] == NULL)
        makeLatticeVertex(_vend, _vpType, numPhases, _vkBest);

      LatticeVertex & lattVertex2 = *lvMat[_vend][_vpType][_vnextPhase];
      LatticeEdge ledge(nb, wordLength, _vlstPType, _vphase, kth, _vscore);

      if(_vfsmEdge->toSelf())
        ledge.update(lattVertex[kth]->backpType(), 
                     lattVertex[kth]->backPhase, 
                     lattVertex[kth]->backKth);
      
      bool inserted = lattVertex2.addLatticeEdge(ledge);
      
      //            if(inserted) {
      //	cerr << "inserted";
	
        //        if(!fe->getSequence()->id().compare("HSENAGENO"))
      //      }
      //        cerr << "  iteration " <<  seqLen - _vbeg + 1 << " " << TypeDefs::nodeName((TParseNode)ge->getParseType()) << " " << TypeDefs::edgeName((TParseEdge)fsmEdge2->id()) << " " << TypeDefs::nodeName(_vpType) << " " << lattEdge->score() << " " << nodeScore << " " << edgeScore << " " << _vscore << " coords " << seqLen - nb.getPos() + 2 << "-" << seqLen - ge->getPos() + 1 << " " << seqLen - nb.getPos() - nb.getLen() + 2 << "-" << seqLen - nb.getPos() + 1 << " " << nb.getStrand() <<  "  " << _vphase << " " << wordLength << endl;
      //      cerr << inserted << " iteration_" <<  _vbeg << " " << fsm->node((TParseNode)ge->getParseType())->name() << " " << fsmEdge2->name() << " " << _vfsmNode->name()  << " " << lattEdge->score() << " " << nodeScore << " " << edgeScore << " " << _vscore << " coords " << ge->getPos() << "-" << nb.getPos() - 1 << " " << nb.getPos() << "-" << _vend << " " << _vphase << " " << _vnextPhase << " " << wordLength << endl;
	      //	    }
    }
  }
  
  /**
   * Performs viterbi at the node level
   */  
  void Lattice::viterbiNode() {

    int seqLen = fe->getSequence()->length();
    vector<Node *> &nextNodes = _vfsmLastNode->nextNodes();
    //    cerr << "# next nodes for last node " << _vfsmLastNode->name() << " " << nextNodes.size() <<endl;
    for(int i = 0 ; (unsigned)i < nextNodes.size() ; i++) {
      _vfsmNode = nextNodes[i];
      //      cerr << i << "th node " << _vfsmNode->name() << " " << _vphase << endl; 
      if(_vphase >= _vfsmNode->maxInPhases())
        continue;
      
      _vcontext = contexts ? contexts->value(_vbeg, _vfsmNode->strand()) : 0;
      _vpType = _vfsmNode->id();
      _vbType = _vfsmNode->id3();
      _vfsmEdge = _vfsmLastNode->nextEdge(_vpType);
      _vstrand = _vfsmNode->strand();
      _vcarriedOverScore = 0;
      _vresetWordLength = false;
      _vnextPhase = _vphase;
      _vsig = NULL;

      // check if we need to carry an initial score assoc with bottlNode
      if(_vlstPType == SYNC_BEG_STATE)
        _vcarriedOverScore = initNodeScores[_vpType][_vphase];

      if(_vbType == syncNode3)       // if we are in a synchronization node
        _vresetWordLength = true;    // then reset the word length

      if(_vfsmEdge->hasSignal()) {
        int fivePSigPos = _vbeg + _vfsmLastNode->nextEdgePos(_vpType);
        _vsig = &signals[_vfsmEdge->id2()]->value(fivePSigPos, _vfsmEdge->strand());
        
        if(_vsig->getType() != _vfsmEdge->id2()) {
	  //	  cerr << "signal isn't right at " << fivePSigPos << " between nodes " << _vfsmLastNode->name() << " " << _vfsmNode->name() << " " << _vsig->getType() << " " << _vfsmEdge->id2() << " " << _vfsmEdge->name() << " " << _vfsmEdge->strand() << endl;
          continue;
	}        
      }
      
      /*
       * Handling self-jumping states
       */
      int minLen = _vfsmNode->minLength(_vcontext);
      
      if(_vfsmNode->nodeType() == FIXEDLEN_NODE) {
        if(_vfsmEdge->toSelf())
          _vlen = 1; // self-transition
        else
          _vlen = minLen;
        _vend = _vbeg + _vlen - 1;
        
        if(_vend > seqLen)
          continue;      
        //        if(!fe->getSequence()->id().compare("HUMHIS102"))
	//	cerr << "    self_"  << _vbeg << " " << _vend << "   state " << fsm->node(_vpType)->name() << " " << _vstrand << " " << _vphase << endl;
        
        _vnextPhase = _vfsmNode->nextPhase(_vphase, _vlen);

        if(_vnextPhase >= _vfsmNode->maxOutPhases())
          continue;

        viterbiEdge();
        continue;
      }
      
      /*
       * Handling variable length/phase states
       */
      vector<TEdgeQ> &sigQs = fsmQ->edgeQueues(_vfsmNode);
      //      cerr << "queues for node " << _vfsmNode->name() << " " << sigQs.size() << endl;
      for(int q = 0; q < sigQs.size(); q++) {
        _vend = _vbeg;
        int sigType = -1;

        _vfsmNextEdge = fsm->edge(sigQs[q].first);
        EdgeQueue *sigQ = sigQs[q].second;
	//	cerr << "queue # " << q << " " << _vfsmNextEdge->name() << endl;
        int expSigType = _vfsmNextEdge->id2();
        QITERATOR head = sigQ->head();
        EdgeInst *sig = NULL;
        
        int phaseStop = _vfsmNode->disruptingPhase(_vbeg, _vphase, _vstrand);

        while((sig = sigQ->dequeue())) {
          sigType = sig->getType();
          
          if(sig->getPos() == seqLen + 1) {
            if(!_vfsmNode->isSyncEnd())  continue;
            else  _vfsmNextEdge = fsm->edge(_vpType, SYNC_END_STATE);
          }

          _vend = sig->getPos() + _vfsmNextEdge->nextNodePos() - 1;
          _vlen = _vend - _vbeg + 1;
	  //	  cerr << "node " << _vfsmNode->name() << " " << _vbeg << " " << _vphase << " signal @ " << sig->getPos() << " " << phaseStop << " " << _vstrand << " " << (phaseStop >= 0 ? sig->eclipsed(phaseStop) : NAN) << endl;
          if(phaseStop >= 0 && sig->eclipsed(phaseStop)) {
	    //	    cerr << "there is a stop disrupting signal @ " << sig->getPos() << " " << phaseStop << " " << _vstrand << endl;
            break;
          }

          /*
           * Rejecting signals which don't correspond or if the node length 
           * is too small
           */
          if(sigType != expSigType || _vlen < minLen) {
	    //	    cerr << "sigs are different " << sigType << " " << expSigType << endl;
            continue;
          }
          _vnextPhase = _vfsmNode->nextPhase(_vphase, _vlen);
          
          if(_vnextPhase >= _vfsmNode->maxOutPhases()) {
	    //	    cerr << "phase is above max " << _vnextPhase << " " << _vfsmNode->maxOutPhases() << endl;
            continue;
          }
          //        if(!fe->getSequence()->id().compare("HUMHIS102"))
	  //	  cerr << "   state " << _vbeg << "-" << _vend << " " << _vfsmNode->name()  << " " << phaseStop << " " << _vfsmNextEdge->name() << " " << sig->eclipsed(0) << " " << sig->eclipsed(1) << " " <<  sig->eclipsed(2) << " " << sigType << " " << expSigType << " " << sig->getPos() << endl << "enter viterbi back edge\n";
          
          viterbiEdge();
        }
        
        sigQ->restoreHead(head);
	//	cerr << "   state-f " << _vbeg << " " << ((TParseNode)_vpType) << " " << " " << _vstrand << " " << _vlen << " " << _vphase << endl;
      }
    }
  }


  /**
   * General viterbi routine. 
   * It outputs a list of Tag* objects in OneStrand format coordinates.
   * Notice that while decoding, the Tag objects need to be converted to 
   * TwoStrand format in order to invoque dotParamV4Nodes and dotParamV4Edges
   * routines.
   * OneStrand format is when coordinates for objects located in STRAND_COMP
   * are with respect to the beginning of the forward Sequence.
   * If an object is in TwoStrand format then its coordinates, when located
   * in the STRAND_COMP are with respect to the beginning of the reversed 
   * and complemented Sequence object.
   */  
  
  double Lattice::viterbi(FeatureEngine & p, 
                          int kBest) {

    int seqLen = fe->getSequence()->length();
    this->_vkBest = kBest;
    this->_vp = &p;
    endState.setPos(seqLen + 1);
    int i = 0, j = 0, l = 0, kth = 0;

    /*
     * initializing lattice vertice matrix
     */
    lvMat = new LatticeVertex *** [seqLen + 2]; 
    assert(lvMat);

    for(i  = 0; i <= seqLen + 1; i++) 
      lvMat[i] = NULL;

    pos2Lattice(0);
    makeLatticeVertex(0, SYNC_BEG_STATE, numPhases, 1);
    
    /* 
     * Initializing with initState vertex at position 0.
     */
    for(_vphase = 0; _vphase < numPhases; _vphase++) {
      LatticeEdge ledge(initState, 0, SYNC_BEG_STATE, 0, 0, 0);
      lvMat[0][SYNC_BEG_STATE][_vphase]->addLatticeEdge(ledge);
    }

    fsmQ->resetQueues();
    for(_vbeg = 1; _vbeg <= seqLen; fsmQ->dequeue(_vbeg), _vbeg++) {

      fsmQ->queue();

      if(lvMat[_vbeg - 1] == NULL)
        continue;

      for(l = 0; l < fsm->numParseNodes(); l++) { 
        _vlstPType = (TParseNode)l;

        /*
         *  checking if states of type _vlstPType end at this position
         */
        if(lvMat[_vbeg - 1][_vlstPType] == NULL) 
	  continue;
	
        _vfsmLastNode = fsm->node(_vlstPType);
	
        for(_vphase = 0; _vphase < numPhases; _vphase++) {
          _vlattVertex = lvMat[_vbeg - 1][_vlstPType][_vphase];  
	  //	  cerr << "last node " << _vfsmLastNode->name() << " pos end " << _vbeg - 1 << " init phase " << _vphase << endl;
          if(_vlattVertex->size() == 0) {
	    //	    cerr << "empty\n";
            continue;
	  }          
          /*
           * these states cannot be in a phase different from zero
           */
          if(_vphase >= _vfsmLastNode->maxOutPhases()) {
	    //	    cerr << "num phases too high\n";
            continue;
	  }
	  //	  cerr << "last_node " << _vbeg << " " << _vfsmLastNode->name() << " " << _vphase << endl;
          viterbiNode();

        }
      }
    }

    int sigPos;
    double edgeScore;

    /*
     * finding the k-best state and final state scores
     */
    pos2Lattice(seqLen + 1);
    makeLatticeVertex(seqLen + 1, SYNC_BEG_STATE, numPhases, kBest);
    LatticeVertex & lattVertex2 = *lvMat[seqLen + 1][SYNC_BEG_STATE][0];

    for(l = 0; l < fsm->numParseNodes(); l++) {
      _vlstPType = (TParseNode)l;
      _vfsmNode = fsm->node(_vlstPType);

      if(!_vfsmNode->isSyncEnd())
        continue;

      _vfsmEdge = fsm->edge(_vlstPType, SYNC_END_STATE);
      TStrand strand = _vfsmEdge->strand();
      sigPos = endState.getPos() + _vfsmNode->nextEdgePos(SYNC_END_STATE);
      assert(sigPos <= seqLen + 1);

      if(_vfsmEdge->hasSignal()) {
        _vsig = &signals[_vfsmEdge->id2()]->value(sigPos, strand);
        _vsig->setParseType(_vfsmEdge->id());
      }
      else {
        noSignal.reInitialize(_vfsmEdge->id(), sigPos, strand);
        _vsig = &noSignal;        
      }

      for(_vphase = 0; _vphase < numPhases; _vphase++) {
        if(!lvMat[seqLen] || !lvMat[seqLen][_vlstPType])
          continue;
        /*
         * these states cannot be in a phase different from zero
         */
        if(_vphase >= _vfsmNode->maxOutPhases())
          continue;

        LatticeVertex &lattVertex = *lvMat[seqLen][_vlstPType][_vphase];

        for(kth = 0; kth < lattVertex.size(); kth++) {
	  _vsig->setAssocNodeLens((strand == STRAND_FWD ? lattVertex[kth]->st.getLen() : 0), (strand == STRAND_FWD ? 0 : lattVertex[kth]->st.getLen()));
	  edgeScore = _vp->dotParamV4Edges(_vsig, 
					   (strand == STRAND_COMP) ?
					   (numPhases - _vphase) % numPhases : 
					   _vphase);

          _vscore = lattVertex[kth]->score() + edgeScore;
          //      cerr << TypeDefs::nodeName((TParseNode)l) << " " << TypeDefs::edgeName((TParseEdge)_vfsmEdge->id()) << " " << lattVertex[kth]->score << " " << edgeScore << " " << score << endl;
          LatticeEdge ledge(endState, 
                            lattVertex[kth]->wordLength, 
                            _vlstPType, _vphase, kth, _vscore);
          lattVertex2.addLatticeEdge(ledge);

        }
      }
    }

    /*
     * ordering k-best parses, from higher to lower score
     */
    for(i = 0; i < lattVertex2.size() - 1; i++) {
      for(j = i + 1; j < lattVertex2.size(); j++)
        if(lattVertex2[i]->score() < lattVertex2[j]->score()) {
          // update
          LatticeEdge t = *lattVertex2[i];
          *lattVertex2[i] = *lattVertex2[j];
          *lattVertex2[j] = t;
        }
    }

    assert(lattVertex2[0]);
    return lattVertex2[0]->score();
    
  }


  /**
   * Backtracking routine
   * @return the kBest seqTags decoded in last viterbi run
   */
  void Lattice::backtrackPath(vector<SeqTags> & seqTags,
                              int kBest) {

    assert(lvMat != NULL);
    int seqLen = fe->getSequence()->length();
    int sigPos, kth;
    int backtrackPos;

    LatticeVertex & lattVertex = *lvMat[seqLen + 1][SYNC_BEG_STATE][0];
    NodeInst *last_b, *b;
    EdgeInst *sig;
    NodeInst endState(NO_NODE_INST, SYNC_END_STATE, seqLen + 1, STRAND_FWD, 1);

    for(kth = 0; kth < lattVertex.size(); kth++) {
      int initPhase = 0;
      LatticeEdge *lattEdge = lattVertex[kth];

      seqTags.push_back(SeqTags(lattEdge->score()));

      SeqTags & listTags = seqTags.back();
      list<Tag *>::iterator bit;

      lattEdge = lattEdge->backtrack(lvMat[seqLen]);
      
      /*
       * Inserting signal for synchronization@seqLen 
       */
      _vfsmNode = fsm->node((TParseNode)lattEdge->st.getParseType());
      last_b = new NodeInst(lattEdge->st, 
                            _vfsmNode->id3(),
                            seqLen + 1 - lattEdge->st.getPos());

      if(!_vfsmNode->isSyncEnd()) {
        delete last_b;
        throw EXCEPTION(PARSE_ERROR, "parse does not finish with sync node");
      }      

      _vfsmEdge = fsm->edge((TParseNode)last_b->getParseType(),
                            SYNC_END_STATE);
      TStrand strand = _vfsmEdge->strand();
      
      sigPos = endState.getPos() + _vfsmNode->nextEdgePos(SYNC_END_STATE);
      assert(sigPos <= seqLen + 1);
      
      if(_vfsmEdge->hasSignal()) {
        sig = new EdgeInst(signals[_vfsmEdge->id2()]->value(sigPos, strand));
        sig->setParseType(_vfsmEdge->id());
#ifndef NDEBUG
	cerr << " BT " << _vfsmEdge->name() << " "  << sigPos << endl;	     
#endif
      }
      else
        sig = new EdgeInst(NO_EDGE_INST,
                           _vfsmEdge->id(),
                           sigPos,
                           strand); 
      
      listTags.tpush_front((Tag *)sig);

      //      cerr << (listTags.front())->getPos() << " " << fsm->edge((TParseEdge)(listTags.front())->getParseType())->name() << endl;

      initPhase = lattEdge->backPhase;
      lattEdge = lattEdge->backtrack(lvMat[last_b->getPos() - 1]);
      _vfsmNode = fsm->node((TParseNode)lattEdge->st.getParseType());

      b = new NodeInst(lattEdge->st, _vfsmNode->id3(), 
                       last_b->getPos() - lattEdge->st.getPos());

      do {

        backtrackPos = b->getPos();
        _vfsmEdge = _vfsmNode->nextEdge(last_b->getParseType());

        if(_vfsmEdge->toSelf()) {
          last_b->setPos(b->getPos());
          // concatenate states with toSelf transitions
          last_b->setLen(last_b->getLen() + b->getLen());
          delete b;
        }
        else {
#ifndef NDEBUG
          cerr << " BT " << fsm->node((TParseNode)last_b->getParseType())->name() << " "  << last_b->getPos() << " " << last_b->getLen() << endl;
#endif
          listTags.tpush_front((Tag *)last_b);
	  sigPos = last_b->getPos() + _vfsmNode->nextEdgePos((TParseNode)last_b->getParseType());
          assert(sigPos <= seqLen + 1);
          EdgeInst *sig = NULL;

          if(_vfsmEdge->hasSignal()) {
            sig = new EdgeInst(signals[_vfsmEdge->id2()]->value(sigPos, _vfsmEdge->strand()));
            sig->setParseType(_vfsmEdge->id());
#ifndef NDEBUG
	    cerr << " BT " << _vfsmEdge->name() << " "  << sigPos << endl;
#endif
          }
          else
            sig = new EdgeInst(NO_EDGE_INST, _vfsmEdge->id(),
                               sigPos, _vfsmEdge->strand());

          listTags.tpush_front((Tag *)sig);

          if(b->getParseType() == SYNC_BEG_STATE) {
            delete b;
            break;
          }

          last_b = b;
          assert(last_b);
          assert(b->getPos() + b->getLen() != 1);
        }

        /*
         * Backtracking to the previous node in the lattice
         */
        initPhase = lattEdge->backPhase;
        lattEdge = lattEdge->backtrack(lvMat[backtrackPos - 1]);
        _vfsmNode = fsm->node((TParseNode)lattEdge->st.getParseType());

        b = new NodeInst(lattEdge->st, _vfsmNode->id3(), 
                         backtrackPos - lattEdge->st.getPos());

       } while(1);

      /*
       * set SeqTag's init and end phases
       */
      listTags.setInitPhase(initPhase);
      listTags.setRank(seqTags.size() - 1);
      listTags.setEndPhase(lattVertex[kth]->backPhase);
      
    }
  }


  /**
   * @return the closest point to the end of the parsing 
   * that is the end of a word, i.e., the beginning of a syncNode. The
   * point is a bottleneck in the sense that all the top k parses must go
   * throughout it.
   */
  int Lattice::findNextBottleneck(TParseNode & bottlNode, 
                                  double & bottlScore,
                                  int minLenSyncSt) {
    
    if(syncNode3 == NO_NODE_INST)
      throw EXCEPTION(NOT_ANNOTATED, "syncNode3 must be initialized for bottleneck computation");

    assert(lvMat != NULL);
    int seqLen = fe->getSequence()->length();
    int len2End = 0;
    double score;
    int backtrackPos;
    map<int, int, intDesCmp> bottlPos;
    map<int, double> bottlScores;
    map<int, TParseNode> bottlNodes;

    /*
     * Accessing the end synchronization vertex in decoding phase 0
     */
    LatticeVertex &lattVertex = *lvMat[seqLen + 1][SYNC_BEG_STATE][0];
    NodeInst *last_b, *b;
    int kth = 0;

    for( ; kth < lattVertex.size(); kth++) {
      // Find the back edge corresponding to the kth-best parse
      LatticeEdge *lattEdge = lattVertex[kth]; 

      lattEdge = lattEdge->backtrack(lvMat[seqLen]);
      
      _vfsmNode = fsm->node((TParseNode)lattEdge->st.getParseType());
      last_b = new NodeInst(lattEdge->st, 
                            _vfsmNode->id3(), 
                            seqLen + 1 - lattEdge->st.getPos());

      bottlNode = (TParseNode)last_b->getParseType();
      
      if(!_vfsmNode->isSyncEnd())
        throw EXCEPTION( PARSE_ERROR, 
                               std::string("Impossible parse: "));
      
      _vfsmEdge = fsm->edge(bottlNode, SYNC_END_STATE);
      
      lattEdge = lattEdge->backtrack(lvMat[last_b->getPos() - 1]);
      score = lattEdge->score();
      _vfsmNode = fsm->node((TParseNode)lattEdge->st.getParseType());
      
      b = new NodeInst(lattEdge->st, _vfsmNode->id3(), 
                       last_b->getPos() - lattEdge->st.getPos());

      do {

        backtrackPos = b->getPos();
        _vfsmEdge = _vfsmNode->nextEdge(bottlNode);
        
        if(_vfsmEdge->toSelf()) {
          last_b->setPos(b->getPos());
          last_b->setLen(last_b->getLen() + b->getLen());
          delete b;
        }
        else {
          len2End += last_b->getLen();
#ifndef NDEBUG
          cerr << " BT " << fsm->node((TParseNode)last_b->getParseType())->name() << " "  << last_b->getPos() << " " << last_b->getLen() << endl;
#endif

          // Have we reached a synchronization state?
          if(last_b->getType() == syncNode3 && len2End > minLenSyncSt) {
            last_b->setLen(last_b->getLen()/2 + 1);
            int pos = last_b->getPos() + last_b->getLen() - 1;
            //            cerr << "igenic beg,end @ " << last_b->getPos() << "," << pos << endl;
	    
            // computing node score
            if(last_b->getStrand() == STRAND_COMP)
              last_b->setPos(seqLen - last_b->getPos() - last_b->getLen() + 2);

            bottlScore = score + _vp->dotParamV4Nodes(last_b, 0);
	    
            if(bottlPos.find(pos) != bottlPos.end()) 
              bottlPos[pos]++;
            else {
              bottlPos[pos] = 1;              
              bottlNodes[pos] = bottlNode;
              bottlScores[pos] = bottlScore;
            }
          }
          
          delete last_b;
          last_b = NULL;
          
          if(b->getParseType() == SYNC_BEG_STATE) {
            delete b;
            break;
          }

          last_b = b;
          bottlNode = (TParseNode)last_b->getParseType();
          
          assert(last_b);
          assert(b->getPos() + b->getLen() != 1);
          
        }
        
        lattEdge = lattEdge->backtrack(lvMat[backtrackPos - 1]);
        score = lattEdge->score();
        _vfsmNode = fsm->node((TParseNode)lattEdge->st.getParseType());
        
        b = new NodeInst(lattEdge->st, _vfsmNode->id3(), 
                         backtrackPos - lattEdge->st.getPos());
        
      } while(1);
      
      if(last_b)
        delete last_b;
    }

    map<int, int, intDesCmp>::iterator it = bottlPos.begin();

    for( ;it != bottlPos.end(); it++) 
      if(it->second == kth) {
        bottlNode = bottlNodes[it->first];
        bottlScore = bottlScores[it->first];
        return it->first;        
      }

    /*
     * We couldn't find the start of any synchronization state in the middle
     * of the sequence. Most likely, all the sequence has been labeled as
     * a synchronization state; we need to take care of this before.
     */
    
    return 0;
      
  }

  /**
   * Find the next position to cut in Sequence's sequence after offset, 
   * which is 'safe', i.e., which does not truncate a word in the middle.
   */
  int Lattice::findNextSafeCut(Sequence &c,
                               int offset, int chunkSize, 
                               int syncStSize,
                               TSetType set) {
    
    if(syncNode3 == NO_NODE_INST)
      return c.length();

    int begin = INT_MAX, end = -1;

    if(chunkSize > c.length() - offset + 1)
      return c.length() - offset + 1;

    vector<SeqTags> & ges = c.getTags(set);
    
    for(unsigned i = 0; i < ges.size(); i++) {
      SeqTags::iterator git;
      
      for(git = ges[i].begin(); git != ges[i].end(); git++) {
        if((*git)->getGEClass() == EDGE_INST)
          continue;
        NodeInst & ge = (NodeInst &)*(*git);
        
        if(ge.getType() == syncNode3)
          continue;
        
        int b = ge.getPos();
        int e = ge.getPos() + ge.getLen() - 1;
        
        begin = (b < begin) ? b : begin;
        end = (e > end) ? e : end;
      }
    }

    if(chunkSize > end - offset + 1) {
      if(chunkSize + syncStSize > c.length())
        return c.length();
      return chunkSize;
    }

    int gap = begin - syncStSize - offset;

    if(gap > syncStSize)
      return gap > chunkSize ? chunkSize : gap;
    
    int cend = (end + syncStSize) > c.length() ? 
      c.length() : end + syncStSize;
    
    if(offset >= cend)
      return (offset + chunkSize) > c.length() ? 
        c.length() - offset + 1: chunkSize;
    else
      return cend - offset + 1;
  }


  void Lattice::updWordInstParams(SeqTags &seqTags,
                                  FeatureEngine &p, 
                                  double updVal) {
                               
    int i = 0;
    std::vector<WordInst *> wordInstances;

    for( ; i < fsm->numParseWords(); i++) {
      Word *word = fsm->words()[i];
      SeqTags::iterator it = seqTags.begin();

      while(it != seqTags.end()) {
        WordInst *wInst = new WordInst(*word, seqTags, it, fsm);
        if(wInst->getParseType() != INVALID_WORD)
          wordInstances.push_back(wInst);
        else
          delete wInst;
      }
    }

    for(i = 0; i < wordInstances.size(); i++) {
      p.updParamV4Words(wordInstances[i], 0, updVal);
      delete wordInstances[i];
    }
  }


  /**
   * Update model parameters, equivalently, it can be used to compute
   * the parameter vector associated with seqTags, if the initial value 
   * of the parameter vector contained in p is zero.
   */  
  void Lattice::updModelParams(SeqTags &seqTags,
                               FeatureEngine &p,
                               double updVal) {
    
    if(!updVal)
      return;
    
    int seqLen = fe->getSequence()->length();
    int phase = seqTags.initPhase();

    updWordInstParams(seqTags, p, updVal);
    
#ifndef NDEBUG
    cerr << fe->getSequence()->id() << " " << updVal << endl;
    cerr << "init phase " << phase << endl;
    int length = 0;
    int totalLength = 0;
#endif
    SeqTags::iterator it = seqTags.begin();
    SeqTags::iterator nit = seqTags.begin();
    NodeInst *b = NULL, *last_b = NULL;
    nit++;
    while(it != seqTags.end()) {
      if((*it)->getGEClass() == NODE_INST) {
        b = (NodeInst *)(*it);
        Node *node = fsm->node((TParseNode)b->getParseType());
	
        if(node->id3() == syncNode3)
          phase = 0;   // must reset phase @ syncNode occurrence       

        int phBreak = seqTags.phaseBreak(*it);
        if(phBreak > 0) // phaseBreak only could happen in between genes   
          phase = phBreak;
        
        if(b->getStrand() == STRAND_COMP) {
          b->setPos(seqLen - b->getPos() - b->getLen() + 2);
          phase = node->nextPhase(phase, b->getLen());
          p.updParamV4Nodes(b, (numPhases - phase) % numPhases, updVal);
          b->setPos(seqLen - b->getPos() - b->getLen() + 2);
        }
        else {
          p.updParamV4Nodes(b, phase, updVal);
          phase = node->nextPhase(phase, b->getLen());
        }

	last_b = b;

#ifndef NDEBUG

        if(node->nodeType() == VARPHASE_NODE) {
          length += b->getLen();
          totalLength += b->getLen();
          cerr << "variable length = " << length << endl;
        }
        else if(b->getType() == syncNode3)
          length = 0;

        cerr << "node " << node->name() << " " << b->getStrand() << " " << b->getPos() << " " << b->getLen() << " " << phase << endl;
#endif

      } 
      else { // EDGE_INST
        assert((*it)->getGEClass() == EDGE_INST);
        EdgeInst *sig = (EdgeInst *)(*it);
	NodeInst *next_b = (nit != seqTags.end() ? 
			    (NodeInst *)(*nit) : NULL);
	//	cerr << "updating edge " << sig->getPos() << sig->getStrand() << " " << updVal << endl;

        if(sig->getStrand() == STRAND_COMP) {
          sig->setPos(seqLen - sig->getPos() + 1);
	  sig->setAssocNodeLens(next_b ? next_b->getLen() : 0,
				last_b ? last_b->getLen() : 0);
          p.updParamV4Edges(sig, (numPhases - phase) % numPhases, updVal);
          sig->setPos(seqLen - sig->getPos() + 1) ;
        }
        else {
	  sig->setAssocNodeLens(last_b ? last_b->getLen() : 0,
				next_b ? next_b->getLen() : 0);
          p.updParamV4Edges(sig, phase, updVal);
	}          
#ifndef NDEBUG
        Edge *edge = fsm->edge((TParseEdge)(*it)->getParseType());
        cerr << "edge " << (*it)->getGEClass() << " " << edge->name() << " " << (*it)->getStrand() << " " << (*it)->getPos() << endl;
#endif
        
      }
      
      it++;
      nit++;
    }
#ifndef NDEBUG
    cerr << "total variable phase state length = " << totalLength << endl;
    cerr << "final phase = " << phase << " " << seqTags.endPhase() << endl;
#endif
    assert(phase == seqTags.endPhase());
    
  }
  
  double Lattice::dotWordInstParams(SeqTags &seqTags,
                                    FeatureEngine &p) {

    
    int i = 0;
    double dotParam = 0;
    std::vector<WordInst *> wordInstances;

    for( ; i < fsm->numParseWords(); i++) {
      Word *word = fsm->words()[i];
      SeqTags::iterator it = seqTags.begin();

      while(it != seqTags.end()) {
        WordInst *wInst = new WordInst(*word, seqTags, it, fsm);
        if(wInst->getParseType() != INVALID_WORD)
          wordInstances.push_back(wInst);
        else
          delete wInst;
      }
    }

    for(i = 0; i < wordInstances.size(); i++) {
      dotParam += p.dotParamV4Words(wordInstances[i], 0);
      delete wordInstances[i];
    }
    
    return dotParam;

  }
  
  double Lattice::dotModelParams(SeqTags &seqTags, FeatureEngine &p,
				 IntervalSet<double> *tag_scores) {

    int seqLen = fe->getSequence()->length();
    SeqTags::iterator it = seqTags.begin();
    SeqTags::iterator nit = seqTags.begin();
    NodeInst *b = NULL, *last_b = NULL;

    int phase = seqTags.initPhase();
    double result = dotWordInstParams(seqTags, p), dotParam, edgeDotParam = 0;
    
#ifndef NDEBUG
    cerr << fe->getSequence()->id() << endl;
    cerr << "init phase " << phase << " " << result << endl;
    int length = 0;
    int totalLength = 0;

#endif
    nit++;

    while(it != seqTags.end()) {
      
      if((*it)->getGEClass() == NODE_INST) {
        b = (NodeInst *)(*it);
        Node *node = fsm->node((TParseNode)b->getParseType());
        
        if(node->id3() == syncNode3)
          phase = 0;   // must reset phase @ syncNode occurrence       

        int phBreak = seqTags.phaseBreak(*it);
        if(phBreak > 0) { // phaseBreak only could happen in between genes 
	  phase = phBreak;
#ifndef NDEBUG
	  cerr << "phase break " << phase << endl;
#endif
	}
                
        if(b->getStrand() == STRAND_COMP) {
          b->setPos(seqLen - b->getPos() - b->getLen() + 2);
          phase = node->nextPhase(phase, b->getLen());
          dotParam = p.dotParamV4Nodes(b, (numPhases - phase) % numPhases);
          b->setPos(seqLen - b->getPos() - b->getLen() + 2);
        }
        else {
          dotParam = p.dotParamV4Nodes(b, phase);
          phase = node->nextPhase(phase, b->getLen());
        }
        
	result += dotParam;

	if(tag_scores) {
	  pair<int,int> node_coords(b->getPos(), b->getPos()+b->getLen()-1);
	  (*tag_scores)[node_coords] = dotParam + edgeDotParam;
	}

#ifndef NDEBUG

        if(node->nodeType() == VARPHASE_NODE) {
          length += b->getLen();
          totalLength += b->getLen();
          cerr << "variable phase state length = " << length << endl;
        }
        else if(b->getType() == syncNode3)
          length = 0;

        cerr << "node " << node->name() << " " << b->getStrand() << " " << b->getPos() << " " << b->getLen() << " " << phase << " " << dotParam << endl;
#endif

      } 
      else { // EDGE_INST
        assert((*it)->getGEClass() == EDGE_INST);
        EdgeInst *sig = (EdgeInst *)(*it);
	NodeInst *next_b = (nit != seqTags.end() ? 
			    (NodeInst *)(*nit) : NULL);
        if(sig->getStrand() == STRAND_COMP) {
          sig->setPos(seqLen - sig->getPos() + 1);
	  sig->setAssocNodeLens(next_b ? next_b->getLen() : 0,
				last_b ? last_b->getLen() : 0);
          edgeDotParam = p.dotParamV4Edges(sig, (numPhases - phase) % numPhases);
          sig->setPos(seqLen - sig->getPos() + 1) ;
        }
        else {
	  sig->setAssocNodeLens(last_b ? last_b->getLen() : 0,
				next_b ? next_b->getLen() : 0);	  
          edgeDotParam = p.dotParamV4Edges(sig, phase);
	}
	
	result += edgeDotParam; 
#ifndef NDEBUG
        Edge *edge = fsm->edge((TParseEdge)(*it)->getParseType());
        cerr << "edge " << (*it)->getGEClass() << " " << edge->name() << " " << (*it)->getStrand() << " " << (*it)->getPos() << " " << edgeDotParam << endl;
#endif
        
      }
      
      it++;
      nit++;
    }
#ifndef NDEBUG
    cerr << "total variable phase states length = " << totalLength << endl;
    cerr << "final phase = " << phase << " " << result << " " << seqTags.endPhase() << endl;
#endif
    assert(phase == seqTags.endPhase());

    return dotParam;
    
  }
  

  bool Lattice::isReachable(vector<SeqTags> & vseqTags, 
			    vector<bool> *reachSeqTags) {

    int seqLen = fe->getSequence()->length();
    bool reachable = false;

    for(int i = 0; i < vseqTags.size(); i++) {
      SeqTags::iterator it = vseqTags[i].begin();
      bool reachable_transcript = true;

      for( ;it != vseqTags[i].end(); it++) {

	if((*it)->getGEClass() == NODE_INST) {
	  NodeInst *b = (NodeInst *)(*it);
	  Node *node = fsm->node((TParseNode)b->getParseType());
	  
	  if(node->nodeType() != FIXEDLEN_NODE)
	    if(b->getLen() > node->maxLength() ||
	       b->getLen() < node->minLength()) {
	      if(verbose)
		cerr << "node " << node->name() << " " << fe->getSequence()->id() << "_" << b->getPos()
		     << "_" << b->getPos() + b->getLen() - 1 << " is unreachable \n";

	      reachable_transcript = false;
	    }
	} 
	else { // EDGE_INST
	  assert((*it)->getGEClass() == EDGE_INST);
	  EdgeInst *sig = (EdgeInst *)(*it);
	  Edge *edge = fsm->edge((TParseEdge)sig->getParseType());

	  if(!edge->hasSignal())
	    continue;

	  EdgeInst *esig = &signals[edge->id2()]->value(sig->getPos(), 
							edge->strand());
	  
	  if(esig->getType() != sig->getType()) {
	    if(verbose)
	      cerr << fe->getSequence()->id() << "_" << sig->getPos() << " " << 
		edge->name() << " is unreachable\n";

	    reachable_transcript = false;
	  }
	}
      }

      if(reachSeqTags)
	reachSeqTags->push_back(reachable_transcript);

      if(reachable_transcript)
	reachable = true;
    }
    
    return reachable;
    
  }

  void Lattice::sortById3TagLen(vector<SeqTags> & vseqTags, 
				vector<int> &lsorted,
				vector<TNodeId3> &id3Tags) {
    
    lsorted.clear();
    vector<pair<int, int> > lengths;
    int i = 0;

    for( ; i < vseqTags.size(); i++) {
      SeqTags::iterator it = vseqTags[i].begin();
      pair<int, int> labelLength;
    
      for( ;it != vseqTags[i].end(); it++) {

	if((*it)->getGEClass() != NODE_INST)
	  continue;

	NodeInst *b = (NodeInst *)(*it);
	
	bool found = false;
	int j = 0;
	for( ; j < id3Tags.size(); j++)
	  if(b->getType() == id3Tags[j]) {
	    found = true;
	    break;
	  }
	
	if(found) {
	  labelLength.first += b->getLen();
	  labelLength.second += (id3Tags.size() - j)*b->getLen();
	}
	
      }
      	
      vector<int>::iterator sit = lsorted.begin();
      vector<pair<int, int> >::iterator lit = lengths.begin();

      for( ; sit != lsorted.end(); sit++, lit++)
	if(labelLength.first > lit->first ||
	   (labelLength.first == lit->first && labelLength.second > lit->second))
	  break;
      
      lengths.insert(lit, labelLength);
      lsorted.insert(sit, i);
    }
    /*
    for(i = 0; i < lsorted.size(); i++) {
      cerr << "index " << i << " " << lsorted[i] << " " << lengths[i].first << " " << lengths[i].second << " " << vseqTags[lsorted[i]].size() << endl;
	
      }*/
  }

  void Lattice::deleteVars() {
    int i, phase;
    int lst;
    int seqLen;
    
    if(lvMat) {
      seqLen = fe->getSequence()->length();
      
      for(i  = 0; i <= seqLen + 1; i++) { 
        if(lvMat[i] == NULL) 
          continue;
        for(lst = 0; lst < fsm->numParseNodes(); lst++) {
          if(lvMat[i][lst] == NULL)
            continue;
          for(phase = 0; phase < numPhases; phase++) {
            if(lvMat[i][lst][phase])
              delete lvMat[i][lst][phase];
            lvMat[i][lst][phase] = NULL;
          }
          delete [] lvMat[i][lst];
          lvMat[i][lst] = NULL;
        }
        delete [] lvMat[i];
        lvMat[i] = NULL;
      }
      delete [] lvMat;
      lvMat = NULL;
    }
  }
}
