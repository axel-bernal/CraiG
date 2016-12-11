#include "FeatureOptimizer.h"

/****************************************************************************
* FeatureOptimizer.cpp - part of the lless namespace, a general purpose
*                        linear semi-markov structure prediction library
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

  FeatureOptimizer::FeatureOptimizer(int numArrays) { 
    int ind;
    _numArrays = numArrays;

    for(int strand = 0; strand < NUM_STRANDS; strand++) {
      _tmpArray[strand] = new double * [_numArrays];
      _array[strand] = new double * [_numArrays];

      for(ind = 0; ind < _numArrays; ind++) {
        _array[strand][ind] = NULL;
        _tmpArray[strand][ind] = NULL;
      }

    }

    _entries = new PreComputedEntry * [_numArrays];
    _allocatedEntries = new bool [_numArrays];

    for(ind = 0; ind < _numArrays; ind++) {
      _entries[ind] = NULL;
      _allocatedEntries[ind] = false;
    }
  }

  int FeatureOptimizer::organize(Feature *f, 
                                 TFeatureTypeEnum type, 
                                 int elem, 
                                 int initArrayIndex,
                                 int subfVal) {
    
    int endArrayIndex = initArrayIndex;
    _allFeatures[type][elem].push_back(f);

    if(initArrayIndex >= 0) {
      assert(_entries[initArrayIndex] && subfVal == -1);
      _pcFeatures[type][elem].push_back(f);
      _entries[initArrayIndex]->addFeature(f);
      _indexes[type][elem].insert(initArrayIndex);
      endArrayIndex += f->preComputeTo(_entries[initArrayIndex]->collapseFrames(), 0);
      assert(endArrayIndex < _numArrays);

      for(int i = initArrayIndex; i <= endArrayIndex; i++)
        _allocatedEntries[i] = true;

    } 
    else {
      _npcFeatures[type][elem].push_back(f);
      _subfVals[type][elem].push_back(subfVal);
    }

    return endArrayIndex;
  }

  void FeatureOptimizer::allocatePreComputedFeatures(int beg, int end, 
                                                     TStrand strand) {

    this->seqLen = end - beg + 1;

    for(int ind = 0; ind < _numArrays; ind++) {

      if(!_allocatedEntries[ind])
        continue;

      _tmpArray[strand][ind] = new double [seqLen + PCARR_PADDING];
      _array[strand][ind] = _tmpArray[strand][ind] + PCARR_PADDING/2;

      for(int i = 0; i < seqLen + PCARR_PADDING; i++)
        _tmpArray[strand][ind][i] = 0.0;

    }
  }

  void FeatureOptimizer::deletePreComputedFeatures() {

    for(int strand = 0; strand < NUM_STRANDS; strand++) {

      if(!_tmpArray[strand])
        continue;

      for(int ind = 0; ind < _numArrays; ind++) {

        if(!_tmpArray[strand][ind])
          continue;

        assert(_allocatedEntries[ind]);
        delete [] _tmpArray[strand][ind];
        _tmpArray[strand][ind] = NULL;
        _array[strand][ind] = NULL;

      }
    }

  }


  void FeatureOptimizer::preComputeFeatures(int beg, int end, 
                                            TStrand strand) {
    /*
     * accumulate cumulative precomputed feature values
     * (those with period != 0)
     */
    for(int ind = 0; ind < _numArrays; ind++) {
      PreComputedEntry *e = _entries[ind];

      if(!e)
        continue;

      Feature *f;

      for(set<Feature *, less<Feature *> >::iterator it = e->features().begin(); it != e->features().end(); it++) {
        f = *it;
        assert(f->isPreComputable());
        f->doPreComputation(_array[strand], e->offset(), beg, end, strand);

	/*        cerr << "PRE " << f->getName() << endl;
	for(int fr = 0; fr < f->parsingFrames(); fr++) {
          for(int i = 0; i < 10; i++)
            cerr << _array[strand][e->offset() + (e->collapseFrames() ? 0 : fr)][i] << " ";
          cerr << endl;
        } */
        /*        if(!e->repFeature()->getName().compare("Donor-Background"))
          cerr << f->getName() << endl;
          for(int i = 61800; i <= 61820; i++)
            if(_array[strand][e->offset()][i])
            cerr << strand << " " << i << " " <<  _array[strand][e->offset()][i] << endl;
        */
      }
      
      if(e->period())
        e->accumulateEntryArrays(beg, end, strand);
      
    }
  }

}
