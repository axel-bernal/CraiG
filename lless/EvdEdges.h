/****************************************************************************
* EvdEdges.h - part of the lless namespace, a general purpose
*              linear semi-markov structure prediction library
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

#ifndef _EVD_EDGES_H_
#define _EVD_EDGES_H_

#include "Utils.h"

namespace lless {

  /**
   * The class represents instances of evdEdges in the lattice whose 3prime 
   * end is at a given position in the input sequence.
   ***************************************************************************/
  
  struct EvdEdge {
    int boundary;
    int weight;

    EvdEdge(int _boundary, int _weight) {
      boundary = _boundary;
      weight = _weight;
    }

    ~EvdEdge() {}
  };

  class EvdEdges {
    //   protected:
    /* [0] -> forward edges, [1] -> backward edges **/
   public:
    vector<EvdEdge> edges[2];
    list<int> active[2][4];
    list<int> inactive[2][4];

    /**
     * Default constructor
     */
    EvdEdges() { }

    /**
     * Copy constructor
     */
    EvdEdges(const EvdEdges &fe) {
      *this = fe;
    }
    
    EvdEdge* bestFwdEdge(int phase) {
      if(active[0][phase].size()) {
	int pos = active[0][phase].front();
	return &edges[0][pos];
      }
      return NULL;
    }
    
    EvdEdge* bestBackEdge(int phase) {
      if(active[1][phase].size()) {
	int pos = active[1][phase].front();
	return &edges[1][pos];
      }
      return NULL;
    }

    void addEvdEdge(UCHAR sense, int phase, 
		     int boundary, int weight) {

      active[sense][phase].push_back(edges[sense].size());
      if(phase == 3) 
	for(int i = 0; i < 3; i++)
	  active[sense][i].push_back(edges[sense].size());
      
      edges[sense].push_back(EvdEdge(boundary, weight));
    }

    void addFwdEdge(int phase, int end, int weight) {
      addEvdEdge(0, phase, end, weight);
    }

    void addBackEdge(int phase, int start, int weight) {
      addEvdEdge(1, phase, start, weight);
    }

    /**
     * normalize weights
     */
    void normalize(int norm_factor) {
      for(UCHAR sense = 0; sense < 2; sense++) 
	for(int j = 0; j < edges[sense].size(); j++) {
	  //	  cerr << "\t" << sense << "," << edges[sense][j].boundary << ",";
	  //	  cerr << edges[sense][j].weight << ",";
	  edges[sense][j].weight *= 100;
	  edges[sense][j].weight /= norm_factor;
	  //	  cerr << edges[sense][j].weight;
	}
      //      cerr << endl;
    }

    // apply log scale
    void logScale() {
      for(UCHAR sense = 0; sense < 2; sense++) 
	for(int j = 0; j < edges[sense].size(); j++) 
	  edges[sense][j].weight = edges[sense][j].weight > 1 ?
	    log(edges[sense][j].weight) : 
	    0.001;
    }

    // apply log scale
    void raise2pow(double power) {
      for(UCHAR sense = 0; sense < 2; sense++) 
	for(int j = 0; j < edges[sense].size(); j++) 
	  edges[sense][j].weight = pow(edges[sense][j].weight, power);
    }
    
    
    EvdEdge *edgeEndingAt(int phase, int end) {
      return retrieveEvdEdge(0, phase, end);
    }

    EvdEdge *edgeStartingAt(int phase, int start) {
      return retrieveEvdEdge(1, phase, start);
    }

    EvdEdge *retrieveEvdEdge(UCHAR sense, int phase, int boundary) {
      list<int>::iterator it = active[sense][phase].begin();
      for( ; it != active[sense][phase].end(); it++) {
	if(edges[sense][(*it)].boundary == boundary)
	  return &edges[sense][(*it)];
      }
      return NULL;
    }
    
    vector<EvdEdge> &allFwdEdges() {
      return edges[0];
    }
    
    vector<EvdEdge> &allBackEdges() {
      return edges[1];
    }

    /**
     * Activate evdEdge
     */
    void activateEvdEdge(UCHAR sense, int phase, int boundary) {
      // moving from inactive to active
      
    }

    inline void activateFwdEdge(int phase, int end) {
      activateEvdEdge(0, phase, end);
    }

    inline void activateBackEdge(int phase, int start) {
      activateEvdEdge(1, phase, start);
    }
    
    /**
     * Inactivate evdEdge
     */
    void inactivateEvdEdge(UCHAR sense, int phase, int boundary) {
      // moving from active to inactive
    }

    inline void inactivateFwdEdge(int phase, int end) {
      inactivateEvdEdge(0, phase, end);
    }
    
    inline void inactivateBackEdge(int phase, int start) {
      inactivateEvdEdge(1, phase, start);
    }

    ~EvdEdges() { }
  };
}

#endif
