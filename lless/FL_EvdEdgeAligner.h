/****************************************************************************
* FL_EvdEdgeAligner.h - part of the lless namespace, a general purpose
*                  linear semi-markov structure prediction library
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

#ifndef _EVDEDGE_ALIGNER_H_
#define _EVDEDGE_ALIGNER_H_

#include "Filter.h"
#include <sstream>
#include "FilterEngine.h"
#include "EvdEdges.h"
#include "TagUtils.h"
#include "InpFile.h"

namespace lless {
  class ResourceEngine;

  /** 
   * FL_EvdEdgeAligner
   * is a subclass of TypedFilter<int> that annotates segment over the 
   * input sequence
   **************************************************************************/

  class FL_EvdEdgeAligner : public TypedFilter<EvdEdges> {
   private:
    int *paddedAccWeights [NUM_STRANDS];
    int *accWeights [NUM_STRANDS];
    bool normalize;
    string upSymbols;
    string dwSymbols;
     
   public:
    //! Default constructor
     FL_EvdEdgeAligner(int fInd,             //!< A unique identifier.
		       std::string name,   //!< A unique name.
		       bool normalize = false,
		       string upSymbols = "D",
		       string dwSymbols = "A"
		       ) :
    TypedFilter<EvdEdges>(fInd, name, 
			  1, FT_EDGE) {
      
      for(int strand = 0; strand < NUM_STRANDS; strand++) {
	paddedAccWeights[strand] = NULL;
	accWeights[strand] = NULL;
      }

      this->normalize = normalize;
      this->upSymbols = upSymbols;
      this->dwSymbols = dwSymbols;

    }
    
    /**
     * Constructor from a Header string definition
     */
    FL_EvdEdgeAligner(int fInd,          //!< A unique identifier.
	   /*! The Header string definition, loaded as a vector of strings.
	    * The Header has the following form:\n\n
	    * Filter name Bin filter numKnots [knot_1 .. ]\n\n
	    * numKnots is the number of knots used to define the bins.
	    * knot_i defines the i^th bin's boundaries i.e.
	    * if defined, bin i would range from knot_i+1 to
	    * knot_{i+1}
	    */
	   vector<std::string> &params, 
	   int & offset,       //!< The index for vector params.
	   ResourceEngine *re, //!< A pointer to the ResourceEngine object.
	   FilterEngine *fe    //!< A pointer to the FilterEngine object.
	   ) : 
    TypedFilter<EvdEdges>(fInd, params, 
			   offset, 1, FT_EDGE) {

      for(int strand = 0; strand < NUM_STRANDS; strand++) {
	paddedAccWeights[strand] = NULL;
	accWeights[strand] = NULL;
      }

      normalize = Utils::stringToBoolean(params[offset++]);
      upSymbols = "D";
      dwSymbols = "A";

      if(offset == params.size()) 
	return;
	
      upSymbols = params[offset++];
      dwSymbols = params[offset++];

    }
    
    void freeValArrays(TStrand strand) {
      if(paddedAccWeights[strand])
	delete [] paddedAccWeights[strand];

      paddedAccWeights[strand] = NULL;
      accWeights[strand] = NULL;

    }

    bool isValidEdge(int upS, int dwS) {
      for(int i = 0; i < upSymbols.size(); i++)
	if(upS == upSymbols[i] && dwS == dwSymbols[i])
	  return true;
      
      return false;
    }

    inline void computeVals(char *seq, EvdEdges *vals, int len) {
      vector<int> *edges = (vector<int> *)seq;
      int maxWeight = 0, i;
      int *accW = accWeights[this->defaultStrand];
      
      if(accW)  assert(0);
      else {
	paddedAccWeights[this->defaultStrand] = new int [len + FILTER_PADDING];
	for(i = 0; i < len + FILTER_PADDING; i++)
	  paddedAccWeights[this->defaultStrand][i] = int();
	
	accWeights[this->defaultStrand] = paddedAccWeights[this->defaultStrand] + FILTER_PADDING/2;
	
	accW = accWeights[this->defaultStrand];
      }
      
      for(i = 0; i < len; i++) {

	for(int j = 1; j < edges[i].size(); j += 3) {
	  int back_p = edges[i][j+1];
	  
	  if(back_p < 0)
	    continue;
	  if(!isValidEdge(edges[back_p - 1][0], edges[i][0]))
	    continue;

	  //	  cerr << "edge from " << i + 1 << "(" << edges[i][0] << ") to " << back_p << "(" << edges[back_p - 1][0] << ")" << endl;
	  vals[i + 1].addBackEdge(edges[i][j], back_p, edges[i][j+2]);
	  vals[back_p].addFwdEdge(edges[i][j], i + 1, edges[i][j+2]);
	  accW[i + 1] = edges[i][j+2]/2 > accW[i+1] ?
	    edges[i][j+2]/2 : accW[i+1];
	  accW[back_p] = edges[i][j+2]/2 > accW[back_p] ? 
	    edges[i][j+2]/2 : accW[back_p];
	    
	  if(edges[i][j+2] > maxWeight)
	    maxWeight = edges[i][j+2];

	}
      }

      for(i = 1; i <= len; i++)
	accW[i] += accW[i - 1];

      //      cerr << "normalize factor = " << maxWeight << " strand=" << this->defaultStrand << endl;
      if(normalize && maxWeight) {
	for(i = 1; i <= len; i++) {
	  vals[i].normalize(maxWeight);
	  //	  cerr << "accW[" << i << "]=" << accW[i] << endl;
	}
      }
    }

    /**
     * A member function that returns an alignment code for an edge
     * with a frame, fiveEnd and threeEnd boundaries given as 
     * arguments.
     * The code is 0 for none, 1 for (partial) containment, 2 for fiveP or 
     * threeP and 3 for full exact alignment
     */
    UCHAR _alignmentType(int frame, pair<int,int> coords, TStrand strand) {
      UCHAR code = 0;
      //      if(coords.first == 1390 && coords.second == 2248 && strand == STRAND_FWD)
      //	cerr << "hola\n";
      EvdEdges & fwdEdges = this->value(coords.first, strand);
      if(fwdEdges.bestFwdEdge(frame)) {
	code = 2;
	code += (fwdEdges.edgeEndingAt(frame, coords.second) != NULL);
      }
      //      cerr << "type_info=" << this->getName() << " " << frame << "," << coords.first << "," << coords.second << " " << (int)code << " " << strand << endl;
      if(code == 3)
	return code;

      EvdEdges & backEdges = this->value(coords.second, strand);
      if(backEdges.bestBackEdge(frame)) {
	code = 2;
	//	code += backEdges.edgeStartingAt(frame, fivePend);
	// if edge does not align forward, it should not align backwards
	// either
      }
      if(code == 3 || code == 2)
	return code;

      int *accW = accWeights[strand];
      //      cerr << strand << " accW[" << coords.first - 1 << "]=" << accW[coords.first - 1] << " accW[" << coords.second << "]=" << accW[coords.second] << endl;
      return (accW[coords.second - 1] - accW[coords.first] != 0);
    }

    int _alignmentWeight(int frame, pair<int,int> coords, TStrand strand) {
      //      cerr << "weight_info=" << this->getName() << " " << frame << "," << coords.first << "," << coords.second << "," << strand << endl;
	
      EvdEdges &fwdEdges = this->value(coords.first, strand);
      EvdEdge *edgeF = fwdEdges.edgeEndingAt(frame, coords.second);
      
      if(edgeF)  return edgeF->weight;

      EvdEdges & backEdges = this->value(coords.second, strand);
      edgeF = fwdEdges.bestFwdEdge(frame);
      EvdEdge *edgeB = backEdges.bestBackEdge(frame);

      if(edgeF && edgeB)  
	return Utils::max(edgeF->weight, edgeB->weight);
      if(edgeB) return edgeB->weight;
      if(edgeF) return edgeF->weight;

      int *accW = accWeights[strand];
      return (accW[coords.second - 1] - accW[coords.first]);
    }

    UCHAR alignmentType(Node *node, Tag *tag, int frame) {
      TStrand strand = tag->getStrand();
      UCHAR code = this->_alignmentType(frame, TagUtils::tagCoords(node, tag), strand);
      //      cerr << " " << (int)code << endl;
      return code;
    }

    int alignmentWeight(Node *node, Tag *tag, int frame) {
      TStrand strand = tag->getStrand();
      return this->_alignmentWeight(frame, TagUtils::tagCoords(node, tag), strand);
    }

    void inactivateEdge(int frame, pair<int,int> coords, TStrand strand) {
      this->value(coords.first, strand).inactivateFwdEdge(frame, coords.second);
      this->value(coords.second, strand).inactivateBackEdge(frame, coords.first);
    }

    void inactivateTag(Node *node, Tag *tag, int frame) {
      TStrand strand = tag->getStrand();
      this->inactivateEdge(frame, TagUtils::tagCoords(node, tag), strand);
    }

    ~FL_EvdEdgeAligner() {
    }
  };
}

#endif
