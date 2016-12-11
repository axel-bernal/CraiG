#include "RerankingCore.h"
#include "Utils.h"
#include <float.h>
#include <map>

/****************************************************************************
* RerankingCore.cpp - part of the lless namespace, a general purpose
*            linear semi-markov structure prediction library
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
   * A member function that performs the common tasks before performing
   * any of the parameter update methods. These tasks are computing for 
   * each y*, F(x, y') - F(x, y*), where y* is the predicted tagging and is
   * the expected tagging with minimum loss with respect to y*.
   * @return the number of restrictions F(x, y) - F(x, y*)
   */
  int RerankingCore::prepareParamUpdate(Sequence &c,
                                        vector<SeqTags> &expTags,  
                                        vector<SeqTags> &predTags,
                                        double learnRate) {
    
    unsigned int kth1, kth2;
    std::string fmt = "locs";
    
    vector<SeqTags> & insTags = c.getTags(_instSet);
    vector<loss_t> instLosses(expTags.size());
    
    for(kth1 = 0; kth1 < expTags.size(); kth1++) {    

      instLosses[kth1] = Core::computeMinLoss(c, expTags[kth1],
					      insTags);

      gv1.push_back(computeGlobalVector(insTags[instLosses[kth1].first], learnRate));
      // gv1 contains all corr_fv
      cerr << "expected " << instLosses[kth1].first << endl;
      gv1.back()->print(cerr);
    }

    int numRestr = 0;
    
    for(kth2 = 0; kth2 < predTags.size(); kth2++) {

      loss_t min = computeMinLoss(c, predTags[kth2], 
				  expTags);
      assert(min.first >= 0);
      
      if(!min.second) {
        if(verbose) 
          cerr << "correct!\n";
        continue;
      }
      
      if(verbose) {
        cerr << "predicted[" << kth2 << "].size() = " << predTags[kth2].size() << endl;
        _printer->displayTags(fmt, (std::ofstream &)cerr, c, predTags[kth2], *_fsm, "PRED");
      }
      
      if(incorr_gv)
	delete incorr_gv;
      incorr_gv = computeGlobalVector(predTags[kth2], learnRate);
      
      // incorr_gv contains F(x, y*)
      cerr << "predicted " << kth2 << endl;
      incorr_gv->print(cerr);
      
      loss.push_back( min.second - instLosses[min.first].second); 
      
      assert(loss[numRestr] >= 0);
      // computing F(x, y) - F(x, y*)
      gv2.push_back(new GlobalVector(_fte->getFeatures(), true));
      (*gv2[numRestr]) = (*gv1[min.first]);
      (*gv2[numRestr]) -= (*incorr_gv);
      
      //      cerr << "restriction " << numRestr << endl;
      //      gv2[numRestr]->print(cerr);
      
      if(verbose) {
        cerr << "expected[" << min.first << "," << kth2 << "].size() = " << expTags[min.first].size() << endl;
        _printer->displayTags(fmt, (std::ofstream &)cerr, c, expTags[min.first], *_fsm, "EXP");
        cerr << "best[" << instLosses[min.first].first << "," << kth2 << "].size() = " << insTags[instLosses[min.first].first].size() << endl;
        _printer->displayTags(fmt, (std::ofstream &)cerr, c, insTags[instLosses[min.first].first], *_fsm, "BEST");
        cerr << "loss[" << numRestr << "] = " << loss[numRestr] << " " << endl;      
      }        
      numRestr++;
      
    }
    
    return numRestr;

  }

  /**
   * A member function that reranks the SeqTag lists according to the 
   * current set of parameters
   * @param superC the input Sequence. 
   */
  double RerankingCore::rerank(Sequence &superC, 
			       int topK) {
    
    double score;
    vector<SeqTags> & insTags = superC.getTags(_instSet);
    vector<double> scores;
    int i = 0;

    for( ; i < insTags.size(); i++) {
      score = lattice()->dotModelParams(insTags[i], *_fte);
      scores.push_back(score);
    }

    for(i = 0; i < scores.size() - 1; i++) {
      for(int j = i + 1; j < scores.size(); j++) 
        if(scores[i] < scores[j]) {
          score = scores[i];
          scores[i] = scores[j];
          scores[j] = score;
	  
          int rank = insTags[i].rank();
          insTags[i].setRank(insTags[j].rank());
          insTags[j].setRank(rank);
	  
        }
    }

    return scores[0];

  }

    // \todo make it work for topK != 1
  double RerankingCore::argmax(Sequence &superC, TSetType predSet, 
			       int topK) {
    
    vector<SeqTags> & insTags = superC.getTags(_instSet);
    double score = rerank(superC, topK);

    int size = superC.getTags(predSet).size();
    SeqTags &predTags = superC.addTags(size, predSet);
    predTags.clone(insTags[insTags[0].rank()]);

    return score;
  }
  
  int RerankingCore::update(GlobalVector &params,
			    Sequence &superC, TSetType trainSet,
			    TSetType predSet, double learnRate, 
			    TTrainMethod tm) {
    
    /*
     * Update parameters if the learning rate is different from ZERO
     */
    if(!learnRate)
      return 0;

    vector<SeqTags> & expTags = superC.getTags(trainSet);      
    vector<SeqTags> & predTags = superC.getTags(predSet);
    
    freetmpFeatVectors();

    switch(tm) {
    case PERCEPTRON: // cannot train with multilabeled sequences
      this->updInstancePerceptron(superC, expTags, predTags, learnRate);
      break;
    case MIRA:
      updInstanceMIRA(superC, expTags, predTags, learnRate);
      break;
    case PEGASOS: // cannot train with multilabeled sequences
      updInstancePEGASOS(superC, expTags, predTags, learnRate);
      break;
    case CWL: // cannot train with multilabeled sequences
      updInstanceCWL(superC, expTags, predTags, learnRate);
      break;
    case ARROW: // cannot train with multilabeled sequences
      updInstanceARROW(superC, expTags, predTags, learnRate);
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
    
    t++;
    //        cerr << "correct - incorrect\n";
    //        params.print(cerr);
    
    freetmpFeatVectors();
    
    return 1;
  }


  int RerankingCore::argmaxAndUpdate(GlobalVector &params, 
				     list<Sequence *> &annotSeqs, 
                                     vector<bool> &activeInstances,
                                     TSetType trainSet,
				     TSetType predSet, 
                                     double learnRate, 
                                     TTrainMethod tm) {
        
    unsigned int kth;
    int syncStSize = 10000;
    int chunkSize = INT_MAX;
    int sampleSize = 0;

    _fte->setParamVector(&params);

    list<Sequence *>::iterator cit = annotSeqs.begin(); 

    for(int i = 0; cit != annotSeqs.end(); cit++, i++) {
      Sequence &superC = *(*cit);

      // is it inside the sample?                                             
      if(activeInstances.size() && !activeInstances[i])
        continue;
      
      prePSequence(superC);

      if(verbose)
        cerr << ">" << superC.id() << " " << superC.length();

      double score = argmax(superC, predSet, 1);
      
      if(verbose)
        cerr << " result " << score << " " << superC.getTags(predSet).size() <<
          " " << superC.getTags(predSet)[0].size() << " " << endl;

      sampleSize += update(params, superC, trainSet, predSet,
                           learnRate, tm);
      
      _fte->updFeatureCacheParams(); // there might be a bug here if params are not zero initially

      postPSequence(true);
    }

    if(learnRate && _avgMethod == AVG_NONE)
      (*accum_gv) = params;

    return sampleSize;
  }


  /**
   * The main function for training and validating the linear structure model.
   * @param pAvg the computed parameters
   * @param annotSeqs the input training sequences
   * @param trainSet the set which contains the expected labeling of each 
   * Sequence in annotSeqs. One of TSetType enumerate.
   * @param learnRate a learning rate
   * @param tm the training method to use.
   * @param paramsFile the file name of the output parameters. It is used to
   * store parameters at the end of each iteration. The file names for the 
   * output parameters are paramFile1, paramFile2 and so on.
   */

  void RerankingCore::startTraining(
                                    GlobalVector  &pAvg,
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
    
    while(iteration < maxIterations) {
      
      list<Sequence *>::iterator cit;
      
      /*********** training *********/
      cerr << "Compute model params using training set..." << endl;

      if(_avgMethod == AVG_LAST) {
        (*accum_gv) = 0;
	if(sigma_gv)
	  (*accsigma_gv) = 0;
      }

      vector<bool> dummy;
      int iterations = argmaxAndUpdate(*pivot_gv, annotSeqs,
                                       dummy, trainSet, 
                                       predSet, learnRate, tm);

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

      if(_combMethod == COMB_KL) {
        pAvg = *accum_gv;
        pAvg.productInverse(*accsigma_gv);
      }
      else
        pAvg.average(accumIterations, *accum_gv);


      if(verbose) {
        cerr << "\nAvg Params " << "\n";
        pAvg.print(cerr);
      }
      
      cerr << "Iteration " << iteration << "\n";
      cerr << "Training Set\t";
      _evaluator->reportAccuracy((std::ofstream&)cerr);

      /* 
       * hack for storing parameters at each iteration
       */
      storeIteration(paramsFile, iteration, accumIterations, &pAvg);

    }
  }

  /**
   * A member function to predict structures on an input sequence.
   * @param superC the input Sequence. 
   * \todo make it work for topK != 1
   */
  double RerankingCore::predictTags(Sequence &superC, TSetType predSet,
				    list<BioFeature> *predRegions) {
    prePSequence(superC);
    double this_result = argmax(superC, predSet, 1);
    postPSequence(true);
    return this_result;
  }
  
  double scoreTags(Sequence &superC, TSetType tagSet,		      
		   list<BioFeature> *regions,
		   IntervalSet<double> &tag_scores) {
    cerr << "In development\n";
    return 0;
  }
  
}
