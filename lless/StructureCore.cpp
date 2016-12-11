#include "StructureCore.h"
#include "Utils.h"
#include <float.h>
#include <map>

/****************************************************************************
* StructureCore.cpp - part of the lless namespace, a general purpose
*                     linear semi-markov structure prediction library
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
   * Core object constructor.
   */
  StructureCore::StructureCore( 
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
                               bool maxLoss,
                               TMultiUpd multiUpd,
                               int oracleWait,
                               TOracleUpd oracleUpd,
			       bool addUnreachable
				) :
    
    Core(fsm, lattice, re, fe, fte, 
         evaluator, printer, avgMethod, 
         combMethod, strand) {
 
    _topK = topK;
    _maxLoss = maxLoss;
    _multiUpd = multiUpd;
    _sampleSizePerc = 0;
    _validSequences = NULL;
    _minLongerLimit = 1;
    _minLongerWait = 1;
    _oracleWait = oracleWait;
    _addUnreachable = addUnreachable;

    if(_oracleWait >= 0) {
      _oracleUpdInWaiting = oracleUpd;
      _oracleUpd = OC_NONE;
    }
    else
      _oracleUpd = oracleUpd;

  }

  int StructureCore::randomExpLabeling(int min, vector<double> &losses,
				       int numLosses) {

    double *probs = new double [numLosses];
    int i;
    double max_delta = 1;

    for(i = 0; i < numLosses; i++)
      if(losses[i] - losses[min] > max_delta)
        max_delta = losses[i] - losses[min] + 1;
    
    for(i = 0; i < numLosses; i++) {
      double delta_i = (losses[i] - losses[min])/max_delta + 1;
      probs[i] = exp(-delta_i*(t^(1/3)));
      if(i > 0)
        probs[i] += probs[i - 1];

      if(verbose)
        cerr << "range " << i << " " << losses[i] << " "  << probs[i]<< " " << delta_i << endl;
    }                             

    double sample = probs[numLosses - 1] * rand()/(RAND_MAX + 1.0);
    int pick = 0;
    for(i = 1; i < numLosses; i++)
      if(sample >= probs[i - 1] && sample <= probs[i]) {
        pick = i;
        break;
      }

    if(verbose)
      cerr << "sample " << sample << " pick @ " << t << " is " << pick << " with loss " << losses[pick] << endl;

    delete [] probs;

    return pick;
  }

  /**
   * A member function that performs the common tasks before performing
   * any of the parameter update methods. These tasks are computing for 
   * each y*, F(x, y') - F(x, y*), where y* is the predicted tagging and is
   * the expected tagging with minimum loss with respect to y*.
   * predictedTags must be in descending order with parsing score as key.
   * @return the number of restrictions F(x, y) - F(x, y*)
   */
  int StructureCore::prepareParamUpdate(Sequence &c,
                                        vector<SeqTags> &expectedTags,  
                                        vector<SeqTags> &predictedTags,
                                        double learnRate) {
    
    unsigned int kth1, kth2;
    int last_numRestr = 0, numRestr = 0, numPredTags;
    
    std::string fmt = "locs";

    loss_t min;

    vector<SeqTags> *expTags = &expectedTags;
    vector<SeqTags> *predTags = &predictedTags;
    vector<SeqTags> tmpTags1, tmpTags2, tmpTags3;
    vector<int> orac2Truth;
    vector<loss_t> oracLosses = Core::computeMinLoss(c, expectedTags, 
						     predictedTags); 

    if(_multiUpd == ML_LONGEST || _multiUpd == ML_MINLONGER) {
      vector<int> lsorted;
      lattice()->sortById3TagLen(expectedTags, lsorted, _id34LenSorting);
      
      for(kth1 = 0; kth1 < Utils::min(_minLongerLimit, lsorted.size()); kth1++) { 
	tmpTags3.resize(tmpTags3.size() + 1);
	tmpTags3[tmpTags3.size() - 1].clone(expectedTags[lsorted[kth1]]);
      }

      expTags = &tmpTags3;
    }

    if(_oracleUpd != OC_NONE && predictedTags.size() > 1) {

      numPredTags = 0;
      
      for(kth2 = 0; kth2 < predictedTags.size(); kth2++) { 
	if(oracLosses[kth2].first < 0) {
	  tmpTags2.resize(tmpTags2.size() + 1);
	  tmpTags2[tmpTags2.size() - 1].clone(predictedTags[kth2]);
	}
	else {
	  numPredTags = tmpTags2.size(); // for OC_BETTER update
	  orac2Truth.resize(orac2Truth.size() + 1);
	  orac2Truth[orac2Truth.size() - 1] = kth2;
	  tmpTags1.resize(tmpTags1.size() + 1);
	  tmpTags1[tmpTags1.size() - 1].clone(predictedTags[kth2]);
	}
      }

      expTags = &tmpTags1;
      predTags = &tmpTags2;

      if(_oracleUpd == OC_TOP)
        numPredTags = 1;
      else if(_oracleUpd == OC_ALL) 
        numPredTags = predTags->size();
      
      if(!numPredTags) { // the best prediction is the oracle, nothing to do
	cerr << "oracle is top prediction, nothing to do\n";
	return 0;
      }
    }
    else
      numPredTags = predTags->size();

    for(kth1 = 0; kth1 < expTags->size(); kth1++) { 
      gv1.push_back(computeGlobalVector((*expTags)[kth1], learnRate));
      
      // gv1 contains all corr_fv
      //      cerr << "expected " << gv1.size() - 1 << endl;
      //      gv1.back()->print(cerr);
      
    }

    vector<loss_t> predLosses = computeMinLoss(c, *expTags,
					       *predTags, numPredTags);
    
    for(kth2 = 0; kth2 < numPredTags; kth2++) { 
      
      min = predLosses[kth2];
      
      if(incorr_gv)
	delete incorr_gv;
      incorr_gv = computeGlobalVector((*predTags)[kth2], learnRate);

      // incorr_gv contains F(x, y*)
      //      cerr << "predicted " << kth2 << endl;
      //      incorr_gv->print(cerr);

      if(verbose) {
	cerr << "predicted[" << kth2 << "].size() = " << (*predTags)[kth2].size() << endl;
	_printer->displayTags(fmt, (std::ofstream &)cerr, c, (*predTags)[kth2], *_fsm, "PRED");

      }

      if(min.first < 0 || !min.second) {
	cerr << "correct!\n"; 
	continue;
      }

      // For ML_ALL, we update with all, ML_LEN has only one labeling so it
      // is also taken care of here.
      int firstLabeling = 0;
      int lastLabeling = expTags->size() - 1;

      vector<double> losses = Core::computeLosses(c, (*predTags)[kth2], 
						  *expTags);
      
      // If ML_MINLONGER is chosen, there may be more than one labeling, so we 
      // chose the one with the minimum loss.
      if(_multiUpd == ML_MIN || _multiUpd == ML_MINLONGER)
	firstLabeling = lastLabeling = min.first;
      else if(_multiUpd == ML_EXP) {
	firstLabeling = lastLabeling = randomExpLabeling(min.first,
							 losses, expTags->size());
      }
      
      for(kth1 = firstLabeling; kth1 <= lastLabeling; kth1++) {
	if(_oracleUpd != OC_NONE &&  predictedTags.size() > 1) { 
	  if(oracLosses[orac2Truth[kth1]].second > losses[kth1]) {
	    cerr << "oracle loss is too big, discarding update\n";
	    continue;
	  }
	}
	
	loss.push_back(losses[kth1]);
        
	assert(loss[numRestr] >= 0);
	// computing F(x, y) - F(x, y*)
	gv2.push_back(new GlobalVector(_fte->getFeatures(), true));
	(*gv2[numRestr]) = (*gv1[kth1]);
	(*gv2[numRestr]) -= (*incorr_gv);
	
	//      cerr << "restriction " << numRestr << endl;
	//      gv2[numRestr]->print(cerr);
        
	if(verbose) {
	  if(_oracleUpd != OC_NONE &&  predictedTags.size() > 1) { 
	    int truthInd = oracLosses[orac2Truth[kth1]].first;
	    cerr << "truth[" << truthInd << "].size() = " << expectedTags[truthInd].size() << endl;
	    _printer->displayTags(fmt, (std::ofstream &)cerr, c, expectedTags[truthInd], *_fsm, "TRUTH");
	  }
	  cerr << "expected[" << kth1 << "," << kth2 << "].size() = " << (*expTags)[kth1].size() << endl;
	  _printer->displayTags(fmt, (std::ofstream &)cerr, c, (*expTags)[kth1], *_fsm, "EXP");
	  cerr << "loss[" << numRestr << "] = " << loss[numRestr] << " " << endl;      
	}        
	
	numRestr++;
      }
      
      if(numRestr == last_numRestr)
	if(verbose)
	  cerr << kth2 << "th prediction is correct!\n";

      last_numRestr = numRestr;
    }
    
    return numRestr;

  }

  /* fragmented version is commented
  double StructureCore::argmax(Sequence &superC, TSetType trainSet,
			       TSetType predSet, int topK) {
  
    /*
     * Depending on the label annotation, divide annotSeq into variable
     * sized fragments 
     */
  /*
    int syncStSize = 10000;
    int chunkSize = INT_MAX;
    int begin = 1, begin2 = 1, end;
    bool fragmented = false;
    double score = 0;
    double syncScore = 0, lastSyncScore = 0;
    int syncPhase = -1;
    TParseNode syncNode = INVALID_NODE;
    
    while(begin < superC.length()) {
      end = begin - 1;
      int thisChunkSize = lattice()->findNextSafeCut(superC, 
                                                     begin,
                                                     chunkSize,
                                                     syncStSize,
                                                     trainSet);
      
      assert(thisChunkSize);
      
      fragmented = (thisChunkSize < superC.length());
      
#ifndef NDEBUG
      cerr << "subseq " << begin2 << " " << begin2 + thisChunkSize - 1 << endl;
#endif
      Sequence *c = fragmented ?
        (Sequence *)superC.getSubSequence(begin2, begin2 + thisChunkSize - 1) :
        &superC;          
      
      prePSequence(*c);
      
      lattice()->setInitScore(syncScore, syncNode, syncPhase);
      score = lattice()->viterbi(*_fte, topK);
      
      /*
       * Find the next lattice vertex in the current parse which is a
       * bottleneck for all the topK parses. Right now it only works
       * for best parse.
       */
  /*
      if(end + thisChunkSize < superC.length())
        end += lattice()->findNextBottleneck(syncNode, syncScore, 1);
      else  // we are at the end
        end += thisChunkSize - 1;
      //        cerr << "bottleneck " << end << endl;
      if(end <= begin2) {  // chunkSize is too small, increase it
        assert(fragmented);
        lattice()->deleteVars(); 
        lattice()->resetInitScores();
	postPSequence(!fragmented);
        
        begin += c->length() + 1;
        delete c;
        continue;
      }
      
      lattice()->backtrackPath(c->resetTags(predSet), topK);
      lattice()->deleteVars();
      lattice()->resetInitScores();
      
      if(fragmented) {
        vector<SeqTags> & tags = superC.getTags(predSet);
        
        for(unsigned int j = 0; j < tags.size(); j++)
          tags[j].setScore(tags[j].score() - lastSyncScore);
        
        lastSyncScore = syncScore;
        syncPhase = tags[0].endPhase();
        int limit = end - begin + ((end == superC.length() - 1) ? 3 : 1);          
        superC.appendTags(syncNode, *_fsm, *c, predSet, begin - 1, limit);
	postPSequence(!fragmented);

        delete c;
      }
      begin = begin2 = end + 1;
    }
    
    if(fragmented)
      prePSequence(superC);

    return score;
  }

  */
  
  double StructureCore::argmax(Sequence &superC, TSetType predSet, int topK) {

    lattice()->setInitScore(0, INVALID_NODE, -1);
    double score = lattice()->viterbi(*_fte, topK);
    lattice()->backtrackPath(superC.resetTags(predSet), topK);
    lattice()->deleteVars();
    lattice()->resetInitScores();
    return score;
  }
  

  int StructureCore::update(GlobalVector &params,
			    Sequence &superC, TSetType trainSet,
			    TSetType predSet, double learnRate, 
			    TTrainMethod tm) {
    
    int numUpdates = 0;
    /*
     * Update parameters if the learning rate is different from ZERO
     */
    if(!learnRate)
      return numUpdates;

    vector<SeqTags> &predTags = superC.getTags(predSet);
    vector<SeqTags> *expTags = &superC.getTags(trainSet);
    
    assert(expTags->size());

    vector<SeqTags> tmpTags;
    unsigned int numInstances = 1;

    if(_multiUpd == ML_SEP)
      numInstances = expTags->size();
    
    for(unsigned int j = 0; j < numInstances; j++) {
      
      if(_multiUpd == ML_SEP) {
        tmpTags.clear();
        tmpTags.resize(1);
        tmpTags[0].clone(superC.getTags(trainSet)[j]);
        expTags = &tmpTags;
      }

      freetmpFeatVectors();

      switch(tm) {
      case PERCEPTRON: // cannot train with multilabeled sequences
	updInstancePerceptron(superC, *expTags, predTags, learnRate);
	break;
      case MIRA:
	updInstanceMIRA(superC, *expTags, predTags, learnRate);
	break;
      case PEGASOS: // cannot train with multilabeled sequences
	updInstancePEGASOS(superC, *expTags, predTags, learnRate);
	break;
      case CWL: // cannot train with multilabeled sequences
	updInstanceCWL(superC, *expTags, predTags, learnRate);
	break;
      case ARROW: // cannot train with multilabeled sequences
	updInstanceARROW(superC, *expTags, predTags, learnRate);
	break;
      }
      
      if(_avgMethod != AVG_NONE)
	if(_combMethod == COMB_KL) {
	  GlobalVector params_x(_fte->getFeatures());
	  (*accsigma_gv) += (*sigma_gv);
	  params_x = params;
	  params_x.product(*sigma_gv);
	  (*accum_gv) += params_x;
	}
	else 
	  (*accum_gv) += params;
      
      numUpdates++;
      t++;
      //      cerr << "correct - incorrect\n";
      //      params.print(cerr);
    }

    freetmpFeatVectors();
    
    return numUpdates;

  }
  
  void StructureCore::sortAnnotSeqs(list<Sequence *> &seqs,
				    list<Sequence *> &annotSeqs,
				    TSetType set, 
				    bool byNumSeqTags,
				    bool byLength) {
    
    if(!byNumSeqTags && !byLength)
      return;
    
    list<Sequence *>::iterator cit = annotSeqs.begin();
      
    for( ; cit != annotSeqs.end(); cit++) {
      int len = (*cit)->length();
      int size = (*cit)->getTags(set).size();
      assert(size != 0);
      
      list<Sequence *>::iterator it = seqs.begin();      
      
      while(it != seqs.end()) {
	int t_len = (*it)->length();
	int t_size = (*it)->getTags(set).size();
	assert(t_size != 0);
	
	if(byNumSeqTags && t_size < size)
	  it++;
	else 
	  if((!byNumSeqTags || t_size == size) && byLength && t_len < len)
	    it++;
	  else
	    break;
       }
      
      seqs.insert(it, *cit);
    }
    
    assert(seqs.size() == annotSeqs.size());

  }

  void StructureCore::activateInstances(vector<bool> & activeInstances,
                                        int inputSize) { 
                        
    if(_sampleSizePerc) {
      while(activeInstances.size() < inputSize)
        activeInstances.push_back(false);

      int sampleSize = (int)(inputSize*_sampleSizePerc/100.0);
      
      int i = 0;
      while(i < sampleSize) {
        int sample = (int)(inputSize * rand()/(RAND_MAX + 1.0));
        
        while(activeInstances[sample])
          sample = (sample + 1) % inputSize;
        
        activeInstances[sample] = true;
        i++;
      }
    }
  }

  /**
   * A member function to viterbi-decode the input sequence and if learnRate is
   * different from zero, to update the parameter vector.
   * @param params a temporary parameter vector
   * @param annotSeqs the input training sequences
   * @param trainSet the set which contains the expected labeling of each 
   * Sequence. One of TSetType enumerate.
   * @param predSet the set which will contain the predicted labeling of each
   * Sequence object in annotSeqs after calling this function. One of 
   * TSetType enumerate.
   * @param learnRate a learning rate
   * @param tm the training method to use.
   * \todo make loss-augmented decoding work for fragmented sequences. 
   * \todo Make appendTags work for contiguous SeqTags which have incompatible
   * edges at the boundaries.
   */

  int StructureCore::argmaxAndUpdate(GlobalVector &params, 
                                     list<Sequence *> &annotSeqs, 
                                     vector<bool> &activeInstances,
                                     TSetType trainSet,
                                     TSetType predSet,
                                     double learnRate, 
                                     TTrainMethod tm) {
    
    unsigned int kth;
    int i = 0;
    int sampleSize = 0;
    int topK;

    _fte->setParamVector(&params);

    list<Sequence *>::iterator cit = annotSeqs.begin();

    for(cit = annotSeqs.begin(); cit != annotSeqs.end(); cit++, i++) {
      Sequence &superC = *(*cit);
      topK = _topK;
      // is it inside the sample?                                             
      if(activeInstances.size() && !activeInstances[i])
        continue;

      prePSequence(superC);

      /*
       * Set the annotation on the evaluator object.
       */

      bool isReachable = lattice()->isReachable(superC.getTags(trainSet));

      if(learnRate) {
	if(_maxLoss) // this does not work for multilabeled seqs
	  _evaluator->setAnnot4Loss(superC, superC.getTags(trainSet)[0]);

	if(_oracleUpd != OC_NONE && isReachable)
	  topK = 1;
      }
      else 
	topK = 1;

      if(verbose)
        cerr << ">" << superC.id() << " " << superC.length();

      double score = argmax(superC, predSet, topK);
      
      if(learnRate && _maxLoss)
        _evaluator->releaseAnnot4Loss();
      
      if(verbose)
        cerr << " result " << score << " " << superC.getTags(predSet).size() <<
          " " << superC.getTags(predSet)[0].size() << " reachable?" << 
	  isReachable  << endl;
      
      if(isReachable || _addUnreachable)
	sampleSize += update(params, superC, trainSet, predSet,
			     learnRate, tm);
      
      _fte->updFeatureCacheParams();

      postPSequence(true);

    }

    if(learnRate && _avgMethod == AVG_NONE)
      (*accum_gv) = params;

    return sampleSize;
  }
                             
  /**
   * The main function for training and validating the linear structure model.
   * @param finalParams the computed parameters
   * @param annotSeqs the input training sequences
   * @param trainSet the set which contains the expected labeling of each 
   * Sequence in annotSeqs. One of TSetType enumerate.
   * @param predSet the set which contains the predicted labeling of each
   * Sequence object in annotSeqs. One of TSetType enumerate.
   * @param learnRate a learning rate
   * @param tm the training method to use.
   * @param paramsFile the file name of the output parameters. It is used to
   * store parameters at the end of each iteration. The file names for the 
   * output parameters are paramFile1, paramFile2 and so on.
   */
  void StructureCore::startTraining(GlobalVector  &finalParams,
                                    list<Sequence *> &annotSeqs, 
                                    TSetType trainSet,
                                    TSetType predSet,
                                    double learnRate, 
                                    int maxIterations, 
                                    TTrainMethod tm, 
                                    char *paramsFile,
                                    int iteration,
                                    int accumIterations) {
    
    t = 1; 
    unsigned int k = 0;
    cerr.precision(10);
    
    if(!pivot_gv)
      pivot_gv = new GlobalVector(_fte->getFeatures());
    if(!accum_gv)
      accum_gv = new GlobalVector(_fte->getFeatures());

    list<Sequence *> tmpAnnotSeqs, *seqs = &annotSeqs;
    vector<bool> activeInstances;

    /*
     * Needs to reorder annotSeqs, putting those with a single
     * labeling at the front, to avoid local minima convergence 
     */
    if(_multiUpd == ML_MIN || _multiUpd == ML_MINLONGER) {
      sortAnnotSeqs(tmpAnnotSeqs, annotSeqs, trainSet, true);
      seqs = &tmpAnnotSeqs;
    }

    activateInstances(activeInstances, seqs->size());
    
    while(iteration < maxIterations) {
      
      list<Sequence *>::iterator cit;
      
      /*********** training *********/
      cerr << "Compute model params..." << endl;
      
      if(_avgMethod == AVG_LAST) {
        (*accum_gv) = 0;
	if(sigma_gv)
	  (*accsigma_gv) = 0;
      }

      if(_oracleWait == iteration)
        _oracleUpd = _oracleUpdInWaiting;
      
      int iterations = argmaxAndUpdate(*pivot_gv, *seqs,
                                       activeInstances,
                                       trainSet,
                                       predSet, 
                                       learnRate, 
                                       tm);

      switch(_avgMethod) {
      case AVG_ALL: accumIterations += iterations; break;
      case AVG_LAST: accumIterations = iterations; break;
      case AVG_NONE: accumIterations = 1; break;
      }

      _evaluator->computePredAccuracy(annotSeqs, trainSet, predSet);

      cit = annotSeqs.begin(); 
      for( ; cit != annotSeqs.end(); cit++)
        (*cit)->resetTags(predSet);

      iteration++;

      if(_multiUpd == ML_MINLONGER && !(iteration % _minLongerWait))
	_minLongerLimit++;

      if(_combMethod == COMB_KL) {
        finalParams = *accum_gv;
        finalParams.productInverse(*accsigma_gv);
      }
      else
        finalParams.average(accumIterations, *accum_gv);

      /******** validating *********
       * for validation only decode, don't update parameters (learnRate = 0)
       */

      cerr << "Validating model params..." << endl;
      vector<bool> dummy;
      argmaxAndUpdate(finalParams, *_validSequences, 
                      dummy,
                      trainSet, predSet, 0);
      
      if(verbose) {
	char validFile[500];
	sprintf(validFile, "%s_valid%d.locs", paramsFile, iteration);
	std::ofstream validStream(validFile);	
	string fmt("locs");
	list<Sequence *>::iterator cit = _validSequences->begin();
	for(cit = _validSequences->begin(); cit != _validSequences->end(); cit++) {
	  Sequence &c = *(*cit);
	  _printer->displayTags(fmt, validStream, c, c.getTags(predSet)[0], *_fsm, "VAL");
	}
	validStream.close();

	//        cerr << "\nOrig Params " << "\n";
	//        pivot_gv->print(cerr);
	//	cerr << "\nAvg Params " << "\n";
	//	finalParams.print(cerr);
      }
      
      cerr << "Iteration " << iteration << "\n";
      cerr << "Training Set\t";
      _evaluator->reportAccuracy((std::ofstream&)cerr);
      
      _evaluator->computePredAccuracy(*_validSequences, trainSet, predSet);           
      cit = _validSequences->begin();
      for( ; cit != _validSequences->end(); cit++)
        (*cit)->resetTags(predSet);

      cerr << "Validation Set\t";
      _evaluator->reportAccuracy((std::ofstream&)cerr);
      
      /* 
       * storing parameters for this iteration
       */
      storeIteration(paramsFile, iteration, accumIterations, &finalParams);

    }
    
  }

  /**
   * A member function to predict structures on an input sequence.
   * @param superC the input Sequence. If this sequence is too long, it is
   * fragmented within the function in regularly sized chunks.
   * @param predSet the set which contains the predicted labeling of
   * Sequence object superC. One of TSetType enumerate.
   * \todo Make appendTags work for contiguous SeqTags which have incompatible
   * edges at the boundaries.
   * @see Sequence::appendTags
   */
  double StructureCore::predictTags(Sequence &superC, TSetType predSet, 
				    list<BioFeature> *predRegions) {
    /* 
     * Divide annotSeq superC into regular size fragments
     */
    int nsync_beg, nsync_end;
    int begin = 1, end;
    int syncStSize = 10000;
    int chunkSize = 300000;

    bool splitAtUps = false, splitAtDws = false;
    bool fragmented = false;
    int myChunkSize = chunkSize;
    int thisChunkSize;
    int iteration = 1;
    double result = 0;
    TParseNode bottlNode = INVALID_NODE;
    int bottlPhase = -1;
    double bottlScore = 0;

    list<BioFeature> tmp_predRegs;
    list<BioFeature> &predRegs = (predRegions ? *predRegions : tmp_predRegs);
    list<BioFeature>::iterator prit = predRegs.begin();
    
    set<Feature *> syncBegFeats = _fte->syncBegFeatures();
    set<Feature *> syncEndFeats = _fte->syncEndFeatures();
    set<Feature *>::iterator it = syncEndFeats.begin();
    
    _fte->updFeatureCacheParams();

    for( ; it != syncEndFeats.end(); it++)  (*it)->turnOff();
    
    while(begin < superC.length()) {

      if(predRegs.size() > 0 && prit == predRegs.end()) 
	break;

      end = begin - 1;
      
      if(prit == predRegs.end()) {
	thisChunkSize = Utils::min(superC.length() - begin + 1, myChunkSize);
	nsync_beg = begin;
	nsync_end = begin + thisChunkSize - 1;
      }
      else {
	nsync_beg = prit->begin();
	nsync_end = prit->end();
	thisChunkSize = nsync_end - nsync_beg + 1;
      }

      splitAtDws = (nsync_end < superC.length());
      fragmented = (thisChunkSize < superC.length());
      
      Sequence *c = fragmented ?
        (Sequence *)superC.getSubSequence(nsync_beg, nsync_end) :
        &superC;

      prePSequence(*c);

      vector<SeqTags> & predTags = c->resetTags(predSet);

      if(!splitAtDws) {
        it = syncEndFeats.begin();

        for( ; it != syncEndFeats.end(); it++)   (*it)->turnOn();

      }
      
      if(splitAtUps) {      
        it = syncBegFeats.begin();        

        for( ;it != syncBegFeats.end(); it++)   (*it)->turnOff();

      }      

      //      cerr << "score " << bottlScore << " " << bottlNode << endl;
      
      lattice()->setInitScore(bottlScore, bottlNode, bottlPhase);
      double this_result = lattice()->viterbi(*_fte, _topK);

      if(prit == predRegs.end()) {
	/*
	 * Find the next lattice vertex in the current parse which is a
	 * bottleneck for all the topK parses. Right now it only works
	 * for best parse.
	 */
	if(end + thisChunkSize < superC.length())
	  end += lattice()->findNextBottleneck(bottlNode, bottlScore,
					       syncStSize);
	else // we are at the end
	  end += thisChunkSize - 1;
	
	//      cerr << "block " << c->id() << " " << begin << " " << end << " " << thisChunkSize << " " << myChunkSize << endl;
	
	if(end > begin + syncStSize) {
	  iteration = 1;
	  myChunkSize = chunkSize;
	}
	else { // chunkSize is too small, increase it
	  myChunkSize += iteration*syncStSize;
	  iteration = 7*iteration;
	}      
      }
      else {
	end = nsync_end - 1;
	prit++;
      }

      if(end > begin) {
	c->setId(superC.id());
	lattice()->backtrackPath(predTags, _topK);
	result += this_result;
	//	cerr << "predicted " << fragmented << " " << predTags[0].size() << " " << predTags[0].initPhase() << " " << predTags[0].endPhase() <<  endl;
	if(fragmented) {
	  // if at the end of the sequence, include the last EdgeInst object
	  int limit = end - begin + ((end == superC.length() - 1) ? 3 : 1);
	  superC.appendTags(bottlNode, *_fsm, *c, predSet,
			    begin - 1, limit);
	  bottlPhase = predTags[0].endPhase();
	  //          if(!seqTagOk) { // problems tying seqTags, need to rerun from checkpt
	}
	
	begin = end + 1;
	splitAtUps = true;
      }
      
      lattice()->deleteVars(); 
      lattice()->resetInitScores();      
      
      postPSequence(begin >= superC.length());
      
      if(fragmented)
        delete c;

    }

    it = syncBegFeats.begin();
          
    for( ; it != syncBegFeats.end(); it++)   (*it)->turnOn(); 

    return result;

  }
  
}
