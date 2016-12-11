#include "SelfLattice.h"

/****************************************************************************
* SelfLattice.cpp - part of the lless namespace, a general purpose
*                   linear semi-markov structure prediction library
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
   * Performs viterbi at the edge level
   */
  void SelfLattice::viterbiEdge() {
                             
    int seqLen = fe->getSequence()->length();
    double edgeScore, nodeScore;
    int wordLength = 0;
    Tag *ge;
    Edge *fsmEdge2;
    LatticeEdge *lattEdge;
    LatticeVertex & lattVertex = (*_vlattVertex);
    int sigPos;

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
      sigPos = nb.getPos() + _vfsmLastNode->nextEdgePos(_vpType);
      
      /*
      * getting signal info
      */
      if(fsmEdge2->hasSignal()) {
        _vsig = &signals[_vfsmEdge->id2()]->value(sigPos, fsmEdge2->strand());
        _vsig->setParseType(fsmEdge2->id());
      }
      else {
        noSignal.reInitialize(fsmEdge2->id(), sigPos, fsmEdge2->strand());
        _vsig = &noSignal;
      }

      //    cerr << "Sig\t" << sigPos << " " << _vpType << " " << nb.getPos() << " " << nb.getLen() << " " << nb.getStrand() << endl;
      
      /*
       * computing edge score
       */
      edgeScore = _vp->dotParamV4Edges(_vsig, 
                                       (fsmEdge2->strand() == STRAND_COMP) ?
                                       (numPhases - _vphase) % numPhases :
                                       _vphase);
      _vscore += edgeScore;
      
      /*
       * computing node score
       */
      if(nb.getStrand() == STRAND_COMP) { 
        nb.setPos(seqLen - nb.getPos() - nb.getLen() + 2);
        nodeScore = _vp->dotParamV4Nodes(&nb,
                                         (numPhases - _vnextPhase)%numPhases);
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

      bool inserted = lattVertex2.addLatticeEdge(ledge);

      if(inserted) {
        //        cerr << "  iteration " <<  seqLen - _vbeg + 1 << " " << TypeDefs::nodeName((TParseNode)ge->getParseType()) << " " << TypeDefs::edgeName((TParseEdge)fsmEdge2->id()) << " " << TypeDefs::nodeName(_vpType) << " " << lattEdge->score() << " " << nodeScore << " " << edgeScore << " " << _vscore << " coords " << seqLen - nb.getPos() + 2 << "-" << seqLen - ge->getPos() + 1 << " " << seqLen - nb.getPos() - nb.getLen() + 2 << "-" << seqLen - nb.getPos() + 1 << " " << nb.getStrand() <<  "  " << _vphase << " " << wordLength << endl;
        //        if(!fe->getSequence()->id().compare("HSENAGENO"))
        //        cerr << "  iteration_" <<  _vbeg << " " << fsm->node((TParseNode)ge->getParseType())->name() << " " << fsmEdge2->name() << " " << _vfsmNode->name()  << " " << lattEdge->score() << " " << nodeScore << " " << edgeScore << " " << _vscore << " coords " << ge->getPos() << "-" << nb.getPos() - 1 << " " << nb.getPos() << "-" << _vend << " " << _vphase << " " << _vnextPhase << " " << wordLength << endl;
      }
    }
  }


  /**
   * Update model parameters, equivalently, it can be used to compute
   * the parameter vector associated with seqTags, if the initial value 
   * of the parameter vector contained in p is zero.
   * Additionally, if any node length is longer than 1, temporary 
   * self-transitions between nodes of length one will be created to 
   * represent the original nodes.
   * \todo make it work with phase breaks
   */  
  void SelfLattice::updModelParams(SeqTags &seqTags,
                                   FeatureEngine &p,
                                   double updVal) {
    
    if(!updVal)
      return;

    int seqLen = fe->getSequence()->length();
    SeqTags::iterator it= seqTags.begin();
    int phase = seqTags.initPhase();
    Edge *edge;

#ifndef NDEBUG
    cerr << fe->getSequence()->id() << " " << updVal << endl;
    cerr << "init phase " << phase << endl;
    int length = 0;
#endif

    while(it != seqTags.end()) {

      if((*it)->getGEClass() == NODE_INST) {
        NodeInst *this_b = (NodeInst *)(*it);
        Node *node = fsm->node((TParseNode)this_b->getParseType());

        if(node->id3() == syncNode3)
          phase = 0;   // must reset phase @ syncNode occurrence       
        
        edge = fsm->edge(node->id(), node->id());
        
        for(int i = 0; i < this_b->getLen(); i++) {

          NodeInst b(*this_b);
          b.reInitialize(b.getPos() + i, 1);
          
          if(i) {
            EdgeInst einst(edge->id2(), edge->id(), b.getPos(), edge->strand());
            
            if(einst.getStrand() == STRAND_COMP) {
              einst.setPos(seqLen - einst.getPos() + 1);
              p.updParamV4Edges(&einst, (numPhases - phase) % numPhases, updVal);
            }
            else
              p.updParamV4Edges(&einst, phase, updVal);
            
          }
          
          if(b.getStrand() == STRAND_COMP) {
            b.setPos(seqLen - b.getPos() + 1);
            phase = node->nextPhase(phase, 1);
            p.updParamV4Nodes(&b, (numPhases - phase) % numPhases, updVal);
          }
          else {
            p.updParamV4Nodes(&b, phase, updVal);
            phase = node->nextPhase(phase, 1);
          }
        }
        
#ifndef NDEBUG

        if(this_b->getType() == EXON) {
          length += this_b->getLen();
          cerr << "coding length = " << length << endl;
        }
        else if(this_b->getType() == INTERGENIC)
          length = 0;
        
        cerr << "node " << node->name() << " " << this_b->getStrand() << " " << this_b->getPos() << " " << this_b->getLen() << " " << phase << endl;
#endif

      } 
      else { // EDGE_INST
        assert((*it)->getGEClass() == EDGE_INST);
        EdgeInst *sig = (EdgeInst *)(*it);
        
        if(sig->getStrand() == STRAND_COMP) {
          sig->setPos(seqLen - sig->getPos() + 1);
          p.updParamV4Edges(sig, (numPhases - phase) % numPhases, updVal);
          sig->setPos(seqLen - sig->getPos() + 1) ;
        }
        else
          p.updParamV4Edges(sig, phase, updVal);
          
#ifndef NDEBUG
        edge = fsm->edge((TParseEdge)(*it)->getParseType());
        cerr << "edge " << (*it)->getGEClass() << " " << edge->name() << " " << (*it)->getStrand() << " " << (*it)->getPos() << endl;
#endif
        
      }
      
      it++;
    }
#ifndef NDEBUG
    cerr << "final phase = " << phase << " " << seqTags.endPhase() << endl;
#endif
    assert(phase == seqTags.endPhase());
    
  }

  double SelfLattice::dotModelParams(SeqTags &seqTags,
                                   FeatureEngine &p) {
    

    int seqLen = fe->getSequence()->length();
    SeqTags::iterator it= seqTags.begin();
    int phase = seqTags.initPhase();
    Edge *edge;
    double dotParam = 0;

#ifndef NDEBUG
    cerr << fe->getSequence()->id() << endl;
    cerr << "init phase " << phase << endl;
    int length = 0;
#endif

    while(it != seqTags.end()) {

      if((*it)->getGEClass() == NODE_INST) {
        NodeInst *this_b = (NodeInst *)(*it);
        Node *node = fsm->node((TParseNode)this_b->getParseType());

        if(node->id3() == syncNode3)
          phase = 0;   // must reset phase @ syncNode occurrence       
        
        edge = fsm->edge(node->id(), node->id());
        
        for(int i = 0; i < this_b->getLen(); i++) {

          NodeInst b(*this_b);
          b.reInitialize(b.getPos() + i, 1);
          
          if(i) {
            EdgeInst einst(edge->id2(), edge->id(), b.getPos(), edge->strand());
            
            if(einst.getStrand() == STRAND_COMP) {
              einst.setPos(seqLen - einst.getPos() + 1);
              dotParam += p.dotParamV4Edges(&einst, (numPhases - phase) % numPhases);
            }
            else
              dotParam += p.dotParamV4Edges(&einst, phase);
            
          }
          
          if(b.getStrand() == STRAND_COMP) {
            b.setPos(seqLen - b.getPos() + 1);
            phase = node->nextPhase(phase, 1);
            dotParam += p.dotParamV4Nodes(&b, (numPhases - phase) % numPhases);
          }
          else {
            dotParam += p.dotParamV4Nodes(&b, phase);
            phase = node->nextPhase(phase, 1);
          }
        }
        
#ifndef NDEBUG

        if(this_b->getType() == EXON) {
          length += this_b->getLen();
          cerr << "coding length = " << length << endl;
        }
        else if(this_b->getType() == INTERGENIC)
          length = 0;
        
        cerr << "node " << node->name() << " " << this_b->getStrand() << " " << this_b->getPos() << " " << this_b->getLen() << " " << phase << endl;
#endif

      } 
      else { // EDGE_INST
        assert((*it)->getGEClass() == EDGE_INST);
        EdgeInst *sig = (EdgeInst *)(*it);
        
        if(sig->getStrand() == STRAND_COMP) {
          sig->setPos(seqLen - sig->getPos() + 1);
          dotParam += p.dotParamV4Edges(sig, (numPhases - phase) % numPhases);
          sig->setPos(seqLen - sig->getPos() + 1) ;
        }
        else
          dotParam += p.dotParamV4Edges(sig, phase);
          
#ifndef NDEBUG
        edge = fsm->edge((TParseEdge)(*it)->getParseType());
        cerr << "edge " << (*it)->getGEClass() << " " << edge->name() << " " << (*it)->getStrand() << " " << (*it)->getPos() << endl;
#endif
        
      }
      
      it++;
    }
#ifndef NDEBUG
    cerr << "final phase = " << phase << " " << seqTags.endPhase() << endl;
#endif
    assert(phase == seqTags.endPhase());

    return dotParam;

  }

}
