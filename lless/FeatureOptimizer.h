/****************************************************************************
* FeatureOptimizer.h - part of the lless namespace, a general purpose
*                      linear semi-markov structure prediction library
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

#ifndef _FEATURE_OPTIMIZER_H_
#define _FEATURE_OPTIMIZER_H_

#include "Feature.h"
#include <set>
#include <vector>
#include <TypeDefs.h>
#include "FSM.h"
#define PCARR_PADDING 100

namespace lless {

  /**
   * A class that is used exclusively for managing features that precompute 
   * their dot product value.
   */
  class PreComputedEntry {
   protected:
    int _offset; //!< this entry's offset in the array of precomputed values
    /**
     * period in the array. if != 0, then the array is cumulative
     */
    int _period;  

    bool _collapseFrames;  //!< collapse all frames into a single array entry 

    /**
     * Representative feature used to compute dotParamV and store it 
     * in the arrays with the precomputed values
     */
    Feature *_repFeature; 

    //! arrays with precomputed values, shared by all entries
    double ***_array;
    set<Feature *, less<Feature *> > _features;

   public:  
    PreComputedEntry(int offset = -1,
                     int period = 0, 
                     bool collapseFrames = false, 
                     Feature *f = NULL) {

      initialize(offset, period, collapseFrames, f);

    }

    inline int offset() {
      return _offset;
    }

    inline int period() {  
      return _period; 
    }

    inline void setArrayPreCompVals(double ***array) {
      _array = array;
    }

    inline void initialize(int offset, int period, 
                           bool collapseFrames, Feature *f) {  
      _offset = offset;
      _period = period;
      _collapseFrames = collapseFrames;
      _repFeature =f ;
    }

    inline void addFeature(Feature *f) {
      _features.insert(f);
    }

    inline set<Feature *, less<Feature *> > & features() {
      return _features;
    }

    inline bool collapseFrames() {
      return _collapseFrames;
    }  

    inline Feature* repFeature() {
      return _repFeature;
    }

    /**
     * Routine to be called only for periodic features which precompute
     * their values
     */
    inline void accumulateEntryArrays(int beg, int end, TStrand strand) {
      assert(period());
      _repFeature->accumulatePreCompEntries(_array[strand], offset(), 
                                            period(), beg, end);
    }

    inline double dotParamV(Tag *tag, int frame) {
      //      cerr << "name " << _repFeature->getName() << " " << _repFeature->preComputedDotParamV(_array[tag->getStrand()], offset(), tag, frame) << "\n";
      return _repFeature->preComputedDotParamV(_array[tag->getStrand()], 
                                               offset(), tag, frame);
    }

    ~PreComputedEntry() {
      _features.clear();
    }
  };

  /**
   * FeatureOptimizer
   *
   * The FeatureOptimizer class manages feature that precompute their values
   * and features that don't. The main objective of this class is to make 
   * this difference between features, transparent to the user while computing
   * dotParamV for either an edge or a node instance in a very efficient way.
   ***************************************************************************/

  class FeatureOptimizer {
   protected:
    int seqLen;
    int _numArrays;
    double **_tmpArray[NUM_STRANDS];
    double **_array[NUM_STRANDS];   //!< precomputed feature values
    bool *_allocatedEntries;
    PreComputedEntry **_entries;    //!< number of precomputed feature entries

    // All features
    vector<Feature *> _allFeatures[NUM_FEAT_TYPES][MAX_NUMPARSE_ELEMS];
    
    //! preComputable features
    vector<Feature *> _pcFeatures[NUM_FEAT_TYPES][MAX_NUMPARSE_ELEMS];
    vector<int> _subfVals[NUM_FEAT_TYPES][MAX_NUMPARSE_ELEMS];

    //! not preComputable features
    vector<Feature *> _npcFeatures[NUM_FEAT_TYPES][MAX_NUMPARSE_ELEMS]; 

    //! preComputeEntry indexes 
    set<int> _indexes[NUM_FEAT_TYPES][MAX_NUMPARSE_ELEMS];     

   public:
    FeatureOptimizer(int numArrays); 
    void preComputeFeatures(int beg, int end, TStrand);
    int organize(Feature *f, TFeatureTypeEnum type, int t, int preComputeIndex = -1, int subfVal = -1);
    void allocatePreComputedFeatures(int beg, int end, TStrand strand);
    void deletePreComputedFeatures();

    inline void setPreComputedEntry(int id, int period, bool collapseFrames, Feature *f) {
      assert(id >= 0 && id < _numArrays);
      _entries[id] = new PreComputedEntry(id, period, collapseFrames, f);
      _entries[id]->setArrayPreCompVals(_array);
    }

    inline bool entryInUse(int id) {
      assert(id >= 0 && id < _numArrays);    
      return _entries[id] != NULL;
    }

    inline int organize(Feature *f, Node *node, 
                        int preComputeIndex = -1, int subfVal = -1) {
      return organize(f, ST, (int)node->id(), preComputeIndex, subfVal);
    }

    inline int organize(Feature *f, Edge *edge, 
                        int preComputeIndex = -1, int subfVal = -1) {
      return organize(f, TR, (int)edge->id(), preComputeIndex, subfVal);
    }
    
    inline int organize(Feature *f, Word *word, 
                        int preComputeIndex = -1, int subfVal = -1) {
      return organize(f, WO, (int)word->id(), preComputeIndex, subfVal);
    }
    
    inline vector<Feature *> & features(TFeatureTypeEnum type, int elem) {
      return _allFeatures[type][elem];
    }

    inline double dotParamV(Tag *tag, TFeatureTypeEnum type, int frame) {
      Feature *f;
      int parseType = tag->getParseType();
      double score = 0;
      /*
       * dot product of precomputed features
       */
      set<int> & indexes = _indexes[type][parseType];
      set<int>::iterator it = indexes.begin(); 
      for( ; it != indexes.end(); it++) {
        score += _entries[*it]->dotParamV(tag, frame);
        //        if(seqLen == 2177)
	//	cerr << (*(_entries[*it]->features().begin()))->getName() << " " << _entries[*it]->dotParamV(tag, frame) << endl;
      } 
      /*
       * dot product of features that can't precompute
       */

      vector<Feature *> &feats = _npcFeatures[type][parseType];      
      vector<int> &subfVals = _subfVals[type][parseType];

      int numFeats = feats.size();
      for(int i = 0; i < numFeats; i++) {
        //        if(seqLen == 2177)
        feats[i]->turnSubFeatOn(subfVals[i]);
	//	cerr << feats[i]->getName() << " " << feats[i]->dotParamV(tag, frame) << endl;
	score += feats[i]->dotParamV(tag, frame);
        feats[i]->turnAllSubFeatsOn();
      }

      return score;
    }

    inline void updParamV(Tag *tag, TFeatureTypeEnum type, 
                          int frame, double updVal) {

      int parseType = tag->getParseType();
      vector<Feature *> &pcfeats = _pcFeatures[type][parseType];
      for(int i = 0; i < pcfeats.size(); i++)
        pcfeats[i]->updParamV(updVal, tag, frame);

      vector<Feature *> &npcfeats = _npcFeatures[type][parseType];
      vector<int> &subfVals = _subfVals[type][parseType];

      for(int i = 0; i < npcfeats.size(); i++) {
        npcfeats[i]->turnSubFeatOn(subfVals[i]);
        npcfeats[i]->updParamV(updVal, tag, frame);
        npcfeats[i]->turnAllSubFeatsOn();
      }
    }

    ~FeatureOptimizer() {
      deletePreComputedFeatures();

      for(int strand = 0; strand < NUM_STRANDS; strand++) {
        if(!_tmpArray[strand])
          continue;
        delete [] _tmpArray[strand];
        delete [] _array[strand];
      }      
      
      for(int ind = 0; ind < _numArrays; ind++)
        if(_entries[ind]) {
          delete _entries[ind];
          _entries[ind] =  NULL;
        }    
      
      if(_entries)
        delete [] _entries;
      _entries = NULL;      
      
      if(_allocatedEntries)
        delete [] _allocatedEntries;
      _allocatedEntries = NULL;

    }
  };

}

#endif
