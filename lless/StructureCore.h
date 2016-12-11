/****************************************************************************
* StructureCore.h - part of the lless namespace, a general purpose
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

#ifndef _STRUCTCORE_ALGORITHMS_H_
#define _STRUCTCORE_ALGORITHMS_H_
#include "Organism.h"
#include "FeatureEngine.h"
#include "FilterEngine.h"
#include "Lattice.h"
#include "FSM.h"
#include "ContextIMM.h"
#include "Evaluator.h"
#include "TagPrinter.h"
#include "GlobalVector.h"
#include "Core.h"

namespace lless {

  /**
   * The StructureCore class integrates resources, filters and features and 
   * implements all the core subroutines fo  training and predicting models 
   * for linear structure prediction. So far there are two training algorithms
   * which have been implemented: k-best perceptron and MIRA - including the 
   * hildreth subroutine to solve the quadratic program.
   * 
   ***************************************************************************/

  class StructureCore : public Core {

   protected:
    Lattice *_lattice;
    TMultiUpd _multiUpd;
    bool  _maxLoss;
    TOracleUpd _oracleUpd, _oracleUpdInWaiting;
    int _topK;
    int _oracleWait;
    list<Sequence *> *_validSequences; 
    double _sampleSizePerc;
    int _minLongerLimit;
    int _minLongerWait;
    vector<TNodeId3> _id34LenSorting;
    bool _addUnreachable;
    
   public:
    StructureCore(
                  FSM &fsm, 
                  Lattice &lattice, 
                  ResourceEngine &re,
                  FilterEngine &fe,
                  FeatureEngine &fte,
                  Evaluator &evaluator, 
                  TagPrinter &printer,
                  TAvgMethod avgMethod,
                  TCombMethod combMethod,
                  TStrand strand,
                  int topK,
                  bool maxLoss = false,
                  TMultiUpd = ML_ALL,
                  int oracleWait = -1,
                  TOracleUpd oracleUpd = OC_NONE,
		  bool addUnreachable = false
                  );
          
    inline void setValidationSet(list<Sequence *> & validSequences) {
      _validSequences = &validSequences;
    }

    inline void enableRandomSampling(double sampleSizePerc) {
      _sampleSizePerc = sampleSizePerc;
    }

    inline void setMinLongerWait(int minLongerWait) {
      _minLongerWait = minLongerWait;
    }

    inline void setId34LenSorting(vector<TNodeId3> & id34LenSorting) {
      _id34LenSorting = id34LenSorting;
    }

    int randomExpLabeling(int min, vector<double> &losses,
			  int numLosses);

    int prepareParamUpdate(Sequence &c,
                           vector<SeqTags> &expTags,  
                           vector<SeqTags> &predTags,
                           double learnRate);

    double argmax(Sequence &superC, TSetType, int topK);
		  
    
    int update(GlobalVector &gv, 
               Sequence &superC, TSetType trainSet,
               TSetType predSet, double learnRate, 
               TTrainMethod tm);
    
    void sortAnnotSeqs(list<Sequence *> &seqs, 
		       list<Sequence *> &annotSeqs,
		       TSetType trainSet, 
		       bool byNumSeqTags,
		       bool byLength = false);
    
    void activateInstances(vector<bool> & activeInstances,
                           int inputSize);
      
    int argmaxAndUpdate(GlobalVector &params, list<Sequence *> &annotSeqs, 
                        vector<bool> &,
                        TSetType trainSet, TSetType predSet, 
                        double learnRate, 
                        TTrainMethod = MIRA);

    void startTraining(GlobalVector & gv, list<Sequence *> &, 
                       TSetType, TSetType, 
                       double learnRate, int maxIterations, 
                       TTrainMethod tm = MIRA,
                       char *paramsFiles = "params",
                       int iteration = 0,
                       int accumIterations = 0);
    
    double predictTags(Sequence &c, TSetType, 
		       list<BioFeature> *predRegions = NULL);
    
    ~StructureCore() {
    }
    
  };

}

#endif
