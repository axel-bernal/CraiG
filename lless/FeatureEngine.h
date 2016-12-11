/****************************************************************************
* FeatureEngine.h - part of the lless namespace, a general purpose
*                linear semi-markov structure prediction library
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

#ifndef _FEAT_ENGINE_H_
#define _FEAT_ENGINE_H_

#include <vector>
#include "Feature.h"
#include "EdgeInst.h"
#include "NodeInst.h"
#include "WordInst.h"
#include "GlobalVector.h"
#include "FeatureOptimizer.h"
#include <map>

#define MAX_PRECOMP_ARRAYS 1024

using namespace lless;

namespace lless {

  class Feature;
  
  /**
   * The class FeatureEngine read Feature Header definitions and creates Feature
   * objects dynamically.
   * Features can also be read from an input feature file, in similar manner
   * as Filter and Resource objects. A typical input feature file can have 
   * three types of declarations: entries for feature precomputation, Feature
   * object and Tag set definitions.\n
   * The feature file can have three types of declarations: entries for 
   * feature precomputation Feature object definitions and tag set definitions.
   * Entries for feature precomputation have the following format:\n\n
   * PreComputeTo    precomp_entry   num_phases      collapse_phases \n\n
   * Each entry in the array for precomputation (precomp_array) will have 
   * length equal to the input sequence, so the precomputed values are defined
   * per input sequence position; Above, precomp_entry is just an offset 
   * precomp_array, num_phases and collapse_phases give information about how
   *  many subentries in precomp_array need to be created.
   * The rule of thumb is one additional subentry for each phase, unless
   * collapse_phases is set to false, which is the case for features which
   * produce phase-dependent values.\n
   * Each entry in the array for precomputation can have more than one feature
   * referencing to it, as long as all features share the same properties as 
   * the entry with respect to how the phases influence in the feature value.
   * A Feature object definition in the input feature file looks like this:\n\n
   * Feature   feature-name    FeatureClass feature-parameters
   * [-> num_tags tag1 [tag2 ..]] [@ precomp_entry]\n\n
   * The first four fields are similar to the description given above for
   * resources and filters. Optionally a feature can be tied to num_tags 
   * sequence Tag(s), in which case the feature is called tied feature. Only
   * tied features can update its parameters. Furthermore, a tied feature can
   * precompute its values to entry precomp_entry.\n
   * Tag set definitions in the input feature file have the following format (
   * for Node sets only):\n\n
   * NodeSet setName = num_tags node_tag1 [node_tag2 ..] 
   * [@ precomp_entry], .. \n\n
   * And similarly for EdgeSet. The objective of Tag set definitions is to
   * avoid cluttering the input feature file. For example the Node set \n\n
   * NodeSet CODING = 8 INIT_EXON_F LAST_EXON_F INTERNAL_EXON_F SINGLE_EXON_F 
   * INIT_EXON_B LAST_EXON_B INTERNAL_EXON_B SINGLE_EXON_B \n\n
   * would become very useful when tying features which are associated with all
   * types of coding regions.
   ***************************************************************************/

  class FeatureEngine {

    /**
     * Structure to store each entry used for feature precomputation
     * as is defined in the input feature file.
     */
    typedef struct {
      int arrIndex;
      int arrPeriod;
      bool collapseFrames;
    } TPCEntry;
    
   private:
    std::map<std::string, std::string> tagSets;
    std::map<int, TPCEntry> pcEntries;
    std::map<std::string, vector<std::string> *> multiFeatures;
    std::map<std::string, pair<string, vector<string> > > featureSets;
    
   protected:
    FeatureOptimizer *fo;
    FilterEngine *fe;
    FSM *fsm;
    vector<Feature *> features;
    vector<std::string> featDescriptions;
    std::map<std::string, Feature *> namedFeatures;
    GlobalVector *params;     //!< Parameter vector
    int usedPCArrIndexes;     //!< number of entries defined for precomputation
    std::string featureHeaders;

    //! the object factory class used to generated on-the-fly Feature objects
    ObjectFactory6
      < Feature,
      int, int, vector<std::string> &, int &,
      FilterEngine *, FeatureEngine *,
      std::string
      > featureFactory;

   public:
    FeatureEngine(
                  FSM &fsm,
                  FilterEngine &fe,
                  std::ifstream &fd
                  ) {

      params = NULL;
      fo = NULL;
      this->fe = &fe;
      this->fsm = &fsm;
      usedPCArrIndexes = 0;
      registerDefaultFeatures();
      readHeaderDefinitions(fd);

    }

    inline FSM *getFSM() {
      return this->fsm;
    }

    /*
     * Feature-related functions 
     */
    inline void freeze(int b, int e) {
      if(e < 0) e = features.size();
      for(unsigned int i = b; i < e; i++)
        features[i]->freeze();
    }

    inline void unfreeze(int b, int e) {
      if(e < 0) e = features.size();
      for(unsigned int i = b; i < e; i++)
        features[i]->unfreeze();
    }

    inline Feature * findFeature(std::string featName) { 
      std::map<std::string, Feature *>::iterator it;
      it = namedFeatures.find(featName);
      
      if(it == namedFeatures.end()) 
        return NULL;
      
      return it->second;
    }
    
    inline Feature * getFeature(std::string featName) { 
      Feature *f = findFeature(featName);
      if(!f) {
        assert(0);
        throw EXCEPTION( PARSE_ERROR, std::string("undefined Feature ") + featName);
      }
      
      return f;
    }
    
    inline void freeze(std::string & name) {
      Feature *f = getFeature(name);
      f->freeze();
    }

    inline void unfreeze(std::string & name) {
      Feature *f = getFeature(name);
      f->unfreeze();
    }

    inline void saveFeatures(std::ofstream &fd) {
      fd << featureHeaders;
    }

    inline void setFeature(std::string &id, std::string &description,
			   Feature *f) {
      assert(f != NULL);
      features.push_back(f);
      featDescriptions.push_back(description);
      namedFeatures[id] = f;
    }

    inline GlobalVector * currentParams() {  return params; }
    inline vector<Feature *> & getFeatures() {  return this->features; }
    inline vector<std::string> & getFeatDescriptions() {
      return this->featDescriptions; 
    }
    inline Feature * getFeature(int index) {  return this->features[index]; }

    inline double dotParamV4Edges(EdgeInst *sig, int frame) {
      return this->fo->dotParamV(sig, TR, frame);
    }

    inline void updParamV4Edges(EdgeInst *sig, int frame, double updVal) {
      if(sig->getParseType() == INVALID_EDGE)
        return;
      this->fo->updParamV(sig, TR, frame, updVal);
    }

    inline double dotParamV4Nodes(NodeInst *b, int frame) {
      return this->fo->dotParamV(b, ST, frame);
    }

    inline void updParamV4Nodes(NodeInst *b, int frame, double updVal) {
      if(b->getParseType() == INVALID_NODE)
        return;
      this->fo->updParamV(b, ST, frame, updVal);
    } 

    inline double dotParamV4Words(WordInst *w, int frame) {
      return this->fo->dotParamV(w, WO, frame);
    }

    inline void updParamV4Words(WordInst *w, int frame, double updVal) {
      if(w->getParseType() == INVALID_WORD)
        return;
      this->fo->updParamV(w, WO, frame, updVal);
    } 

    /* 
     * Param-related functions 
     */
    inline GlobalVector * getParamVector() {
      return this->params;
    }

    inline void setParamVector(GlobalVector *params) {
      this->params = params;
      for(unsigned int i = 0; i < this->features.size(); i++) { 
        Feature *f = this->features[i]; 
        f->setParams(params->values());
      }
    }

    inline void updFeatureCacheParams() {
      for(unsigned int i = 0; i < this->features.size(); i++) { 
        Feature *f = this->features[i];

	if(f->paramInd() < 0)
	  continue;

	f->updCacheParams((const FeatureVector ***)currentParams()->values());
      }
    }

    inline vector<std::string> *multiFeatureEntries(string multiFeatureName) {
      std::map<string, vector<std::string> *>::iterator it;
      it = multiFeatures.find(multiFeatureName);
      
      if(it == multiFeatures.end())
        return NULL;
      
      return it->second;
    } 

    inline pair<string, vector<string> > * findFeatureSet(string featureSetName) {
      std::map<std::string, pair<string, vector<string> > >::iterator it;
      it = featureSets.find(featureSetName);
      
      if(it == featureSets.end())
        return NULL;
      
      return &it->second;
    }

    void registerDefaultFeatures();
    void preComputeFeatures(int beg, int end, TStrand strand);
    void deletePreComputedFeatures();
    void tieTagsToFeature(Feature *f, std::string tagGroups);
    bool readFeatureHeader(int &, std::string &);
    void readHeaderDefinitions(std::ifstream &fd);
    set<Feature *> syncEndFeatures();
    set<Feature *> syncBegFeatures();
    ~FeatureEngine();
    
  };
  
}

#endif
