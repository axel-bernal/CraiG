#include "MPStructureCore.h"
#include "Utils.h"
#include <float.h>
#include <map>
#include "TagUtils.h"
#include "InpFile.h"

#ifndef WANT_MPI
#undef HAVE_MPI
#endif

/****************************************************************************
* MPStructureCore.cpp - part of the lless namespace, a general purpose
*                       linear semi-markov structure prediction library
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


  /************************************************************************
   *
   * MPStructureCore routines
   *
   ***********************************************************************/
  
#ifdef HAVE_MPI
  void MPStructureCore::MPIsend(MPI::Intracomm &comm, 
				std::string &msg, int to) {

    int len = msg.length();    
    comm.Send(&len, 1, MPI::INT, to, MSGTAG);
    comm.Send(msg.c_str(), len, MPI::CHAR, to, MSGTAG);
  }
  
  TMyMPITag MPStructureCore::MPIrecv(MPI::Intracomm &comm,
				     std::string &msg, int &from) {

    int len;
    MPI::Status status;
    
    comm.Recv(&len, 1, MPI::INT, MPI::ANY_SOURCE, 
                         MPI::ANY_TAG, status);
    from = ((int)status.Get_source());
    
    /* Check the tag of the received message. */
    if (status.Get_tag() == DIETAG && from == 0)
      return DIETAG;

    if(status.Get_tag() == SYNCTAG && !len)
      return SYNCTAG;

    assert(status.Get_tag() == MSGTAG);

    char *c_msg = new char [len + 1];
    comm.Recv(c_msg, len, MPI::CHAR, from, MPI::ANY_TAG, status);

    assert(status.Get_tag() == MSGTAG);
    c_msg[ len ] = '\0';
    msg = std::string(c_msg);
    delete [] c_msg;

    return MSGTAG;

  }
  
  void MPStructureCore::MPIbarrier(MPI::Intracomm &comm, int rangeIni, 
				   int rangeEnd, int exp_ack) {
    
    int ack, rank;
    
    MPI::Status status;
    
    for(int i = rangeIni; i < rangeEnd; i++) {
      //      cerr << "receiving in barrier (" << MPI::COMM_WORLD.Get_rank() << ")\n";
      comm.Recv(&ack, 1, MPI::INT, MPI::ANY_SOURCE, 
		MPI::ANY_TAG, status);
      rank = (( int )status.Get_source());
      //      cerr << "received " << status.Get_tag() << " from " << status.Get_source() 
      //	   <<  " (" << MPI::COMM_WORLD.Get_rank() << ")" << endl;
      if(ack != exp_ack || status.Get_tag() != SYNCTAG)
        throw EXCEPTION(BAD_USAGE, "something wrong with barrier");
      
    }    

  }

  void MPStructureCore::gatherNodeParams(GlobalVector &finalParams,
                                         char *paramsFile,			      
                                         int numProcs,
                                         int iteration,
                                         int &accumIterations) {
    
    char nodeFile[500];
    accumIterations = 0;
    finalParams = 0;
    (*accsigma_gv) = 0;
    GlobalVector sigma_x(_fte->getFeatures(), true);

    for(int rank = 0; rank < numProcs - 1; rank++) {
      sprintf(nodeFile, "node%d_%s%d", rank, paramsFile, iteration);
      std::ifstream paramStream(nodeFile);
      int this_iter, this_accumIter;
      paramStream >> this_accumIter >> this_iter;
      accumIterations += this_accumIter;
      pivot_gv->retrieve(paramStream);
      accum_gv->retrieve(paramStream);
      
      paramStream.get(); 

      if(!paramStream.eof()) {
        paramStream.unget();
        sigma_gv->retrieve(paramStream);
        sigma_x.retrieve(paramStream); // accsigma_gv for this node
      }
      
      if(_combMethod == COMB_KL) {
	(*accsigma_gv) += sigma_x;
        accum_gv->product(sigma_x);
	finalParams += (*accum_gv);
      }
      else
        finalParams += (*accum_gv);

      paramStream.close();

    }

    if(_combMethod == COMB_KL)
      finalParams.productInverse(*accsigma_gv);
    else
      finalParams.average(accumIterations, finalParams, true);

  }

  void MPStructureCore::master(GlobalVector  &finalParams,
                               list<Sequence *> &annotSeqs,
                               TSetType trainSet,
                               TSetType predSet,
                               double learnRate, 
                               int maxIterations, 
                               TTrainMethod tm, 
                               int blockSize,
                               bool shareModels,
                               bool balanceLoad,
                               char *paramsFile,
                               int iteration,
                               int accumIterations) {
    
    t = 1; 
    unsigned int k = 0;
    cerr.precision(10);
    int rank, i = 0,  ack = 0;
    list<Sequence *> tmpAnnotSeqs, *seqs = &annotSeqs;
    vector<bool> activeInstances;
    
    if(!accum_gv)
      accum_gv = new GlobalVector(_fte->getFeatures());
    if(!pivot_gv)
      pivot_gv = new GlobalVector(_fte->getFeatures());
    if(!sigma_gv)
      sigma_gv = new GlobalVector(_fte->getFeatures(), true, 1.0);
    if(!accsigma_gv)
      accsigma_gv = new GlobalVector(_fte->getFeatures());

    int numProcs = MPI::COMM_WORLD.Get_size();

    /*
     * Needs to reorder annotSeqs, putting those with a single
     * labeling at the front, to avoid local minima convergence 
     */

    if(_multiUpd == ML_MIN || balanceLoad) {
      sortAnnotSeqs(tmpAnnotSeqs, annotSeqs, trainSet,
		    _multiUpd == ML_MIN, balanceLoad);
      seqs = &tmpAnnotSeqs;
    }
    
    activateInstances(activeInstances, seqs->size());

    list<Sequence *>::iterator cit = seqs->end();

    while(iteration < maxIterations) {
      
      /*********** training *********/
      cerr << "Compute model params in parallel (#" << numProcs << ")" << endl;
      bool endReached = false;

      if(cit == seqs->end())  {
        cit = seqs->begin();
        i = 0;
      }
      
      while(!endReached && cit != seqs->end()) {
	
        for(rank = 1; rank < numProcs; rank++) {

	  std::string ids4Decoding = "";

	  for(int block = 0; block < blockSize; cit++, i++) {
	    if(cit == seqs->end()) {
	      endReached = true;
	      cit = seqs->begin();
	      i = 0;
	    }
	    
	    // is it inside the sample?                                             
	    if(activeInstances.size() && !activeInstances[i])
	      continue;

	    ids4Decoding = ids4Decoding + (*cit)->id() + string(" ");
	    block++;

	  }

	  //	  cerr << "sending " << (*cit)->id() << " to " << rank << endl;
	  MPIsend(MPI::COMM_WORLD, ids4Decoding, rank);
	  //	  cerr << "sent to " << rank << "\n";
	}

	//	cerr << "waiting barrier for child decoding (" << 0 << ")\n";

        MPIbarrier(MPI::COMM_WORLD, 1, numProcs, ack);
	//	cerr << "passed barrier for child decoding (0)\n";
      }
      
      iteration++;
      
      for(rank = 1; rank < numProcs; rank++)
	MPI::COMM_WORLD.Send(&ack, 1, MPI::INT, rank, SYNCTAG);
      
      //      msComm.Bcast(&ack, 1, MPI::INT, 0);

      MPIbarrier(MPI::COMM_WORLD, 1, numProcs, ack);

      gatherNodeParams(finalParams, paramsFile, numProcs, 
                       iteration, accumIterations);
      
      /******** validating *********
       * for validation only decode, don't update parameters (learnRate = 0)
       */
      
      cerr << "Validating model params..." << endl;
      vector<bool> dummy;
      argmaxAndUpdate(finalParams, *_validSequences, 
                      dummy,
                      trainSet, predSet, 0);
      
      
      if(verbose) {
        cerr << "\nAvg Params " << "\n";
        finalParams.print();
      }
      
      cerr << "Iteration " << iteration << "\n";
      
      _evaluator->computePredAccuracy(*_validSequences, trainSet, predSet);     
      list<Sequence *>::iterator cit = _validSequences->begin();
      
      for( ; cit != _validSequences->end(); cit++)
        (*cit)->resetTags(predSet);
      
      cerr << "Validation Set\t";
      _evaluator->reportAccuracy((std::ofstream&)cerr);
      
      /* 
       * storing parameters for this iteration
       */
      storeIteration(paramsFile, iteration, accumIterations, 
		     &finalParams, false);
    }
    
    /* Tell all the workers to exit by sending an empty message with the
       DIETAG. */
    for(rank = 1; rank < numProcs; rank++)
      MPI::COMM_WORLD.Send(&ack, 1, MPI::INT, rank, DIETAG);
    
      //    MPI::COMM_WORLD.Bcast(&ack, 1, MPI::INT, DIETAG);
    
  }
  
  void MPStructureCore::worker(InpDir &trainDir,
			       string & tagDir,
			       TSetType trainSet,
			       TSetType predSet,
			       double learnRate, 
			       int maxIterations, 
			       MPI::Intracomm &comm,
			       TTrainMethod tm, 
			       bool shareModels,
			       char *paramsFile,
			       int iteration,
			       int accumIterations) {
    
    t = 1; 
    unsigned int k = 0;
    cerr.precision(10);
    std::string ids4Decoding, seqId;
    int ack, p; 
    int iterations = 0;
    Sequence *c;
    int numProcs = comm.Get_size();
    int rank = comm.Get_rank();

    if(!pivot_gv)
      pivot_gv = new GlobalVector(_fte->getFeatures());
    if(!accum_gv)
      accum_gv = new GlobalVector(_fte->getFeatures());
    
    while (1) {
      
      /* Receive a message from the master */
      TMyMPITag myTag = MPIrecv(MPI::COMM_WORLD, ids4Decoding, p);

      if(myTag == DIETAG)
        break;

      if(myTag == SYNCTAG) { // received a end of iteration ack
        iteration++;

        /* 
         *  storing parameters for this iteration
         */
        switch(_avgMethod) {
        case AVG_ALL: accumIterations += iterations; break;
        case AVG_LAST: accumIterations = iterations; break;
        case AVG_NONE: accumIterations = 1; break;
        }      
        iterations = 0;
        
        char nodeFile[500];
	sprintf(nodeFile, "node%d_%s%d", rank, paramsFile, iteration);
	storeLog(nodeFile, iteration, accumIterations);
        
        if(_avgMethod == AVG_LAST)
          (*accum_gv) = 0;
        
        if(_oracleWait == iteration)
          _oracleUpd = _oracleUpdInWaiting;
        
	comm.Barrier();
        /* Send the ack that writing finished */        
	MPI::COMM_WORLD.Send(&ack, 1, MPI::INT, 0, SYNCTAG);

	continue;
      }
      
      if(!ids4Decoding.length() || myTag != MSGTAG || p)
	throw EXCEPTION(BAD_USAGE, "Something seriously wrong seqId");
	
      string::size_type base = 0, index;
      list<Sequence *> annotSeqs;
      vector<bool> dummy;
      std::ostringstream tagFile;

      while(index = ids4Decoding.find(" ", base),
	    index != std::string::npos) {
	
	seqId = ids4Decoding.substr(base, index - base);
	base = index + 1;
	
	c = (Sequence *)trainDir.findSeq(seqId);
	
	if(!c) {
	  assert(0);
	  throw EXCEPTION(CONTIG_UNAVAILABLE, seqId);
	}
	
	annotSeqs.push_back(c);
	
	tagFile.str("");
	tagFile << tagDir << "/" << seqId;
	TagUtils::loadSeqTags(tagFile.str(), annotSeqs, *_fsm, trainSet);

      }	
            
      iterations += argmaxAndUpdate(*pivot_gv, annotSeqs,
				    dummy,
				    trainSet,
				    predSet, 
				    learnRate, 
				    tm);
      

      /* synchronize with other ranks */
      tagFile.str("");
      tagFile << "hypothesis." << rank;
      ofstream tagOStream(tagFile.str().c_str());
      TagUtils::saveSeqTags(tagOStream, *c, *_fsm, predSet);
			 
      tagOStream.close();
      trainDir.uncacheSeqs();
      
      comm.Barrier();

      ack = 0;

      for(int i = 0; i < numProcs; i++) {
	if(i == rank)
	  continue;

	tagFile.str("");
	tagFile << "hypothesis." << i;
	ifstream tagIStream(tagFile.str().c_str());
	tagIStream >> seqId;
	seqId = seqId.substr(1);
	tagIStream.close();
	
        c = (Sequence *)trainDir.findSeq(seqId);
	//	cerr << "received " << seqId << " from " << p << " (" << rank << ")" << endl;
	
        if(!c) {
          assert(0);
          throw EXCEPTION(CONTIG_UNAVAILABLE, seqId);
        }

        annotSeqs.clear();	        
        annotSeqs.push_back(c);
        TagUtils::loadSeqTags(tagFile.str(), annotSeqs, *_fsm, predSet);

        tagFile.str("");
        tagFile << tagDir << "/" << seqId;
        TagUtils::loadSeqTags(tagFile.str(), annotSeqs, *_fsm, trainSet);
	
	/*	cerr << "received at " << rank << "\n";
	string fmt = "locs";
	cerr << "expected[0].size() = " << c->getTags(trainSet)[0].size() << endl;
	_printer->displayTags(fmt, (std::ofstream &)cerr, *c, c->getTags(trainSet)[0], *_fsm, "PRED");
	cerr << "predicted[0].size() = " << c->getTags(predSet)[0].size() << endl;
	_printer->displayTags(fmt, (std::ofstream &)cerr, *c, c->getTags(predSet)[0], *_fsm, "PRED");
	*/
	prePSequence(*c);
	iterations += update(*pivot_gv, *c, trainSet, predSet, learnRate, tm);
	postPSequence(true);
	
      }

      trainDir.uncacheSeqs();
      
      comm.Barrier();
      // tell master process we are done updating
      MPI::COMM_WORLD.Send(&ack, 1, MPI::INT, 0, SYNCTAG);
    }
  }
#endif  

}

