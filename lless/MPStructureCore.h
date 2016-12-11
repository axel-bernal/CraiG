/****************************************************************************
* MPStructureCore.h - part of the lless namespace, a general purpose
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

#ifndef _MPSTRUCTCORE_ALGORITHMS_H_
#define _MPSTRUCTCORE_ALGORITHMS_H_
#include "Organism.h"
#include "FeatureEngine.h"
#include "FilterEngine.h"
#include "Lattice.h"
#include "FSM.h"
#include "ContextIMM.h"
#include "Evaluator.h"
#include "TagPrinter.h"
#include "GlobalVector.h"
#include "StructureCore.h"
#include "InpFile.h"

#ifdef HAVE_CONFIG_H
#include "config/config.h"
#endif

#ifndef WANT_MPI
#undef HAVE_MPI
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

typedef enum {SYNCTAG=1, MSGTAG, DIETAG} TMyMPITag;

namespace lless {

  /**
   * The MPStructureCore is a StructureCore class which uses OpenMPI to
   * paralelize the decoding process
   * 
   ***************************************************************************/

#ifdef HAVE_MPI
  class MPStructureCore : public StructureCore {
  public:
    MPStructureCore(
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
                    TMultiUpd ml = ML_ALL,
                    int oracleWait = -1,
                    TOracleUpd oracleUpd = OC_NONE,
		    bool addUnreachable = false
                    ) : StructureCore(fsm, 
                                      lattice, re,
                                      fe, fte,
                                      evaluator, 
                                      printer,
                                      avgMethod,
                                      combMethod,
                                      strand,
                                      topK,
                                      maxLoss,
                                      ml, oracleWait,
                                      oracleUpd,
				      addUnreachable)  {
      
      
    }
    
    void MPIsend(MPI::Intracomm &comm, std::string & msg, int to);
    TMyMPITag MPIrecv(MPI::Intracomm &comm, std::string & msg, int &from);
    void MPIbarrier(MPI::Intracomm &comm, int rangeIni,
		    int rangeEnd, int exp_ack);

    void gatherNodeParams(GlobalVector &finalParams,
			  char *paramsFile,
                          int numProcs, 
			  int iteration,
			  int &accumIterations);
                          
    
    void master(GlobalVector  &finalParams,
                list<Sequence *> &annotSeqs,
                TSetType trainSet,
                TSetType predSet,
                double learnRate, 
                int maxIterations, 
                TTrainMethod tm = MIRA,
                int blockSize = 1,
                bool shareModels = false,
                bool balanceLoad = false,
                char *paramsFile = "params",
                int iteration = 0,
                int accumIterations = 0);
    
    void worker(InpDir &trainDir,
		string &tagDir,
		TSetType trainSet,
		TSetType predSet,
		double learnRate, 
		int maxIterations, 
		MPI::Intracomm &comm,
		TTrainMethod tm = MIRA, 
		bool shareModels = false,
		char *paramsFile = "params",
		int iteration = 0,
		int accumIterations = 0);
    
  };  
#endif
}

#endif
