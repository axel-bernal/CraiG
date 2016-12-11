/****************************************************************************
* RerankingCore.h - part of the lless namespace, a general purpose
*          linear semi-markov structure prediction library
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

#ifndef _RERANKINGCORE_ALGORITHMS_H_
#define _RERANKINGCORE_ALGORITHMS_H_
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
   * The Core class integrates resources, filters and features and implements
   * all the core subroutines fo  training and predicting models for linear 
   * structure prediction. So far there are two training algorithms which 
   * have been implemented: k-best perceptron and MIRA - including the 
   * hildreth subroutine to solve the quadratic program.
   * 
   ***************************************************************************/

  class RerankingCore : public Core {
   protected:
    TSetType _instSet;

   public:
  /**
   * RerankingCore object constructor.
   */
    RerankingCore(
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
                   TSetType instSet) :
      Core(fsm, lattice, re, fe, fte, 
           evaluator, printer, 
           avgMethod, combMethod, strand) {

      _instSet = instSet;

    }
    
    int prepareParamUpdate(Sequence &c,
                           vector<SeqTags> &expTags,  
                           vector<SeqTags> &predTags,
                           double learnRate);

    double argmax(Sequence &superC, TSetType, int topK);

    int update(GlobalVector &params,
	       Sequence &superC, TSetType trainSet,
	       TSetType predSet, double learnRate, 
	       TTrainMethod tm);
    
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

    double rerank(Sequence &superC, int topK);
    double predictTags(Sequence &superC, TSetType,
		       list<BioFeature> *predRegions = NULL);
        
    ~RerankingCore() {
    }

  };

}
#endif
