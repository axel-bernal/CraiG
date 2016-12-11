#include "FL_EvdEdgeAligner.h"
#include "FL_ScoreFile.h"
#include "RSqUtils.h"
#include <numeric>

/****************************************************************************
 * RNASeqUtils.cpp - part of the craig namespace, a genomics library
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

/**
 * A namespace created for all classes which are part of craig, a genomics 
 * library.
 */

extern bool verbose;

namespace craig {
  
  ostream& operator<< (ostream &out, RSqInterval &ival) {
    out << " interval (" << ival.start << "," << ival.end << ") " 
	<< ival.strand << " " << " mean=" << ival.cov_mean << " stdev="
	<< ival.cov_stdev << " quartiles= " << ival.cov_quartiles[0] << " " 
	<<  ival.cov_quartiles[1] << " " << ival.cov_quartiles[2];	 
    return out;
  }

  void RSqInterval::computeCoverageStats(vector<float>::iterator first,
					 vector<float>::iterator last) {
    if(last - first == 0)
      return;

    std::nth_element(first, first + (last - first)/4, last);
    cov_quartiles.push_back(*(first + (last - first)/4));
    std::nth_element(first, first + (last - first)/2, last);
    cov_quartiles.push_back(*(first + (last - first)/2));
    std::nth_element(first, first + 3*(last - first)/4, last);
    cov_quartiles.push_back(*(first + 3*(last - first)/4));

    cov_mean = accumulate(first, last, 0.0f)/(last - first);

    if(last - first <= 1) 
      return;

    vector<float> zero_mean(last - first);
    for(vector<float>::iterator it = first; it != last; it++)
      zero_mean[it - first] = *it;
    transform( zero_mean.begin(), zero_mean.end(), zero_mean.begin(),bind2nd( minus<float>(), cov_mean ) );
    
    cov_stdev = inner_product( zero_mean.begin(),zero_mean.end(), zero_mean.begin(), 0.0f );
    cov_stdev = sqrt( cov_stdev / (last - first - 1 ) );
  }

  bool RSqUtils::interval_match(list<RSqInterval> &I,
				int my_start, int my_end) {
    
    list<RSqInterval>::iterator it = I.size() ? I.begin() : I.end();
    while(it != I.end() && it->start < my_start)
      it++;
    
    if(it->start == my_start && it->end == my_end)
      return true;
    
    return false;
  }


  // if missing = false, checking routine is a bit more relaxed. We don't 
  // want to punish the annotations for 'wrong' introns which just happen
  // not to be in the RNA-Seq naive predicion.. we don't trust the later
  // that much
  bool RSqUtils::validateIntron(int lbegin, int lend, int ibegin, int iend,
				int rbegin, int rend, float junction_support,
				bool missing, vector<float> & cov) {
    
    bool validated = true;
    float intron_cov = (cov[iend] -  cov[ibegin - 1])/(iend - ibegin + 1);
    
    float lexon_cov = (cov[lend] - cov[lbegin - 1])/(lend - lbegin + 1);
    float rexon_cov = (cov[rend] - cov[rbegin - 1])/(rend - rbegin + 1);
    float exon_cov = 0.5*(lexon_cov + rexon_cov);
    float intron_cov_gap = exon_cov - intron_cov - junction_support;
    if(intron_cov_gap == 0) intron_cov_gap = 1;

    if(exon_cov > 500) {
      if(!missing && junction_support < 0.05*exon_cov) {
	if(junction_support/intron_cov_gap < 0.25)
	  validated = false;
      }
      if(missing && junction_support > 0.075*exon_cov) {
	if(junction_support/intron_cov_gap > 0.5)
	  validated = false;
      }
    }
    else if(exon_cov > 50) {
      if(!missing && junction_support < 0.1*exon_cov) {
	if(junction_support/intron_cov_gap < 0.5)
	  validated = false;
      }
      if(missing && exon_cov > 50 && junction_support > 0.15*exon_cov) {
	if(junction_support/intron_cov_gap > 0.75)
	  validated = false;
      }
    }
    else if(exon_cov > 10) {
      if(!missing && junction_support < 0.2*exon_cov) {
	if(junction_support/intron_cov_gap < 0.75)
	  validated = false;
      }
      if(missing && exon_cov > 10 && junction_support > 0.30*exon_cov) {
	if(junction_support/intron_cov_gap > 0.9)
	  validated = false;
      }
    }
    else {
      if(missing && junction_support < 0.2*exon_cov && 
	 junction_support/intron_cov_gap < 1)
	validated = false;
      if(!missing && junction_support > 0.6*exon_cov && 
	 junction_support/intron_cov_gap > 1)
	validated = false;
    }

    cerr << (missing ? "missing": "wrong") << " intron?\n";
    cerr << "intron w (" << ibegin << ", " << iend << ")=" 
	 << junction_support << " " << intron_cov_gap
	 << "\ncov intron " << intron_cov
	 << "\n cov lexon (" << lbegin << ", " << lend << ")=" << lexon_cov 
	 << "\n cov rexon (" << rbegin << ", " << rend << ")=" << rexon_cov
	 << " validated = " << validated << "\n";	      
    
    return validated;
  }
  
  void RSqUtils::insert_interval(list<RSqInterval> &I, 
				 int my_start, int my_end, 
				 float weight, TStrand strand) {
  
    list<RSqInterval>::iterator it = I.size() ? I.begin() : I.end();
    while(it != I.end() && it->start < my_start)
      it++;
    
    RSqInterval rv(my_start, my_end, weight, strand);
    I.insert(it, (const RSqInterval &)rv); 
  }

  bool RSqUtils::hash_interval_match(vector<pair<int, int> > &I, 
				     int my_start, int my_end) {
    
    for(int k = my_start; k < I.size(); k++) {
      if(!I[k].first)
	return false;
      if(I[k].first = my_start && I[k].second == my_end)
	return true;
    }

    return false;
  }	    
  
  bool RSqUtils::insert_hash_interval(pair<int, int> ival,
				      vector<pair<int, int> > &I) {
    
    for(int k = ival.first; k < I.size(); k++) {
      if(!I[k].first) {
	I[k] = ival;
	return true;
      }
    }

    return false;
  }

  void RSqUtils::computeSeqCoverage(Sequence &c,
				    TypedFilter<double> *coverage,
				    vector<float> *seq_cov,
				    bool stranded) {
        
    int loStrand = 0, hiStrand = 1;

    if(stranded)
      hiStrand = NUM_STRANDS;

    for(int s = loStrand; s < hiStrand; s++) {
      TStrand strand = (TStrand)s;
      seq_cov[s] = vector<float>(c.length() + 1);      

      for(int i = 1; i <= c.length(); i++) {
	float cov = coverage->value(i, strand);
	int pos = (strand == STRAND_COMP ? c.length() - i + 1 : i);

	//	cout << strand << " " << pos << " " << cov << endl;
	seq_cov[s][pos] = cov;
      }
    }
  }


  void RSqUtils::computeJunctionCoverage(Sequence &c,
					 FL_EvdEdgeAligner *junctions,
					 vector<float> *junction_cov,
					 bool stranded) {    
    
    for(int s = 0; s < NUM_STRANDS; s++) {
      TStrand strand = (TStrand)s;
      int leftB = 1;
      int rightB;

      if(stranded || !s)
	junction_cov[s] = vector<float>(c.length() + 1);

      for(leftB = 1; leftB <= c.length(); leftB++) {
	EvdEdges &edges = junctions->value(leftB, strand);
	vector<EvdEdge> &fwdEdges = edges.allFwdEdges();
	
	if(!fwdEdges.size())  continue;
	
	int start;
	int end;
	if(strand == STRAND_COMP)
	  end = c.length() - leftB + 1;
	else
	  start = leftB;

	for(int j = 0; j < fwdEdges.size(); j++) {
	  int k;
	  rightB = fwdEdges[j].boundary;
	  float weight = fwdEdges[j].weight;
	  
	  if(strand == STRAND_COMP)
	    start = c.length() - rightB + 1;
	  else
	    end = rightB;
	  
	  //	  cerr << "edge " << strand << " " << start << " " << end << " " << weight << endl;
	  if(stranded)
	    for(k = start; k <= end; k++)
	      junction_cov[s][k] += weight;
	  else
	    for(k = start; k <= end; k++)
	      junction_cov[0][k] += weight;
	}
      }
    }
  }

  void RSqUtils::getRSqIntrons(Sequence &c, FL_EvdEdgeAligner *junctions,
			       list<RSqInterval> &introns, 
			       TStrand origStrand) {

    int loStrand = 0, hiStrand = NUM_STRANDS;
    
    if(origStrand != BOTH_STRANDS) {
      loStrand = origStrand;
      hiStrand = loStrand + 1;
    }
    
    for(int s = loStrand; s < hiStrand; s++) {
      TStrand strand = (TStrand)s;
      int leftB = 1;

      for( ; leftB <= c.length(); leftB++) {
	EvdEdges &edges = junctions->value(leftB, strand);
	vector<EvdEdge> &fwdEdges = edges.allFwdEdges();
	
	if(!fwdEdges.size())  continue;
	
	int start;
	int end;
	if(strand == STRAND_COMP)
	  start = c.length() - leftB + 1;
	else
	  start = leftB;
	
	for(int j = 0; j < fwdEdges.size(); j++) {
	  int k;
	  int rightB = fwdEdges[j].boundary;
	  float weight = fwdEdges[j].weight;
	  if(strand == STRAND_FWD)
	    end = rightB;
	  else
	    end = c.length() - rightB + 1;
	  
	  //	    bool annotated_junction = RSqUtils::hash_interval_match(I, start, end);
	  //	    cout << annotated_junction << " read junction " << start << " " << end << " " << weight << " " << strand << endl;
	  if(strand == STRAND_FWD)
	    RSqUtils::insert_interval(introns, start, end, weight, strand);
	  else // reverting coordinates for sorting COMP strand reads
	    RSqUtils::insert_interval(introns, end, start, weight, strand);
	}
      }
    }
  }

  /* if keep_bestolap is true, only the overlapping introns with the best
     weight will be kept */
  
  void RSqUtils::getRSqIntrons(Sequence &c,
			       FL_EvdEdgeAligner *junctions,
			       bool keep_bestolap,
			       vector<RSqInterval> &introns,
			       TStrand origStrand) {
    
    int loStrand = 0, hiStrand = NUM_STRANDS;
    
    if(origStrand != BOTH_STRANDS) {
      loStrand = origStrand;
      hiStrand = loStrand + 1;
    }
    
    for(int s = loStrand; s < hiStrand; s++) {
      TStrand strand = (TStrand)s;
      int leftB = 1;
      int rightB;
      
      for(leftB = 1; leftB <= c.length(); leftB++) {
	EvdEdges &edges = junctions->value(leftB, strand);
	vector<EvdEdge> &fwdEdges = edges.allFwdEdges();
	
	if(!fwdEdges.size())  continue;
	
	int start;
	int end;
	if(strand == STRAND_COMP)
	  end = c.length() - leftB + 1;
	else
	  start = leftB;
	for(int j = 0; j < fwdEdges.size(); j++) {
	  int k;
	  rightB = fwdEdges[j].boundary;
	  float weight = fwdEdges[j].weight;
	  if(strand == STRAND_FWD)
	    end = rightB;
	  else
	    start = c.length() - rightB + 1;
	  
	  RSqInterval intron(start, end, strand);
	  intron.setWeight(weight);
	  if(keep_bestolap && introns.size()) {
	    int min, max;
	    if(strand == STRAND_COMP) {
	      min = Utils::min(introns.back().start, introns.back().end);
	      max = Utils::max(start, end);
	    }
	    else {
	      max = Utils::max(introns.back().start, introns.back().end);
	      min = Utils::min(start, end);
	    }

	    if(max >= min - 3) { // correcting as these are ss positions, no
	                         // intron lengths
	      if(weight > introns.back().weight)
		introns.back() = intron;
	    }
	    else  introns.push_back(intron);
	  }
	  else  introns.push_back(intron);
	}
      }
    }
  }  
  
  bool RSqUtils::computeSplSCoverage(TypedFilter<EdgeInst> **signals,
				     vector<float> &seq_cov, 
				     int window,
				     int begin, int end,
				     vector<int> &splpos,
				     vector<float> &splcov,
				     float zero_cutoff,
				     float zero_value,
				     TTranscriptRegion region,
				     TStrand strand) {
    
    int i;
    for(i = 0; i <= window; i++) {
      splpos.push_back(0);
      splcov.push_back(0);
    }

    TStrand compStrand = strand;
    if(strand == BOTH_STRANDS) {
      strand = STRAND_FWD;
      compStrand = STRAND_COMP;
    }

    float cov = DOUBLE_INFINITY;
    for(i = begin; i <= end; i++) {
      int sigPos = -1;
      TEdgeId2 sigId = NO_EDGE_INST;;
      EdgeInst *sig = NULL;
      cov = (cov > seq_cov[i] ? seq_cov[i] : cov);

      if(i == begin || i == end) {
	sigPos = i;
      }
      else {
	if(region != TR_DOWNSTREAM) {
	  sigId = (strand == STRAND_FWD ? ACCEPTOR : DONOR);
	  sig = &signals[sigId]->value(i, strand);

	  if(!sig && compStrand != strand) {
	    sigId = (compStrand == STRAND_FWD ? ACCEPTOR : DONOR);  
	    sig = &signals[sigId]->value(i, compStrand);
	  }

	  if(sig && sig->getType() == sigId)
	    sigPos = i;
	}
	if(region != TR_UPSTREAM) {
	  sigId = (strand == STRAND_FWD ? DONOR : ACCEPTOR);
	  sig = &signals[sigId]->value(i, strand);

	  if(!sig && compStrand != strand) {
	    sigId = (compStrand == STRAND_FWD ? DONOR : ACCEPTOR);
	    sig = &signals[sigId]->value(i, compStrand);
	  }

	  if(sig && sig->getType() == sigId)
	    sigPos = i;
	}
      }

      if(sigPos > 0) {
	cov = (cov <= zero_cutoff ? zero_value : cov);
	//	cout << "score " << splcov.size() << " " << sigPos << " " << sigId << " " << cov << endl;
	splpos.push_back(sigPos);
	if(region != TR_DOWNSTREAM) {
	  if(splcov.size())
	    splcov.back() = cov;
	}
	splcov.push_back(cov);
	cov = DOUBLE_INFINITY;
      }
    }

    if(splcov.size() > window + 1) {
      for(i = 0; i <= window; i++) {
	splcov[i] = splcov[window + 1];
	splpos.push_back(splpos.back());
	splcov.push_back(splcov.back());
      } 
      return true;
    }

    return false;
  }

  float RSqUtils::findMean(vector<float> &X, int start, int end) {  
    /*    vector<float>::iterator first = X.begin() + start;
    vector<float>::iterator last = X.begin() + end;
    return accumulate(first, last, 0.0f)/X.size(); */
    float meanX = 0;
    for(int i = start; i < end; i++)
      meanX += X[i];
    
    int len = end - start;
    if(len) return meanX/len;
    return 0;
  }

  float RSqUtils::findSdiff(vector<float> &X, int start, int end) {
    // finding mean
    float meanX = findMean(X, start, end);

    // find cumulative score S
    vector<float> S(1, 0);
    //    cout.precision(4);
    //    cout << "S = ";
    for(int i = start; i < end; i++) {
      S.push_back(S[i - start] + X[i] - meanX);
      //      cout << " " << S.back();
    }

    // find Smax and Smin
    pair<float, float> Smin_max = Utils::findMinMax(S, 0, S.size());
    //    cout << " min,max = " << Smin_max.first << "," << Smin_max.second << endl;    
    return Smin_max.second - Smin_max.first;
  }
  
  float RSqUtils::findMSE(vector<float> &X, int start, int end) {
    // finding mean
    float Xmean = findMean(X, start, end);
    float MSE = 0;
    for(int i = start; i < end; i++)
      MSE += pow(X[i] - Xmean, 2.0);
    
    return MSE;
  }


  float RSqUtils::meanShiftConfidence(vector<float> &X, 
				      int start, int end,
				      int num_bootstraps,
				      int blockLen) {
    
    // check whether blockLen is too big
    vector<pair<int, int> > bindexes;
    int i = 0;

    while(1) {
      if(i >= end - start)  break;
      int bSize = (i + blockLen < end - start) ? blockLen : end - start - i;
      bindexes.push_back(pair<int, int>(i, bSize));
      i += bSize;
    }
    
    // see if there is a change
    int num_changed = 0;

    vector<float> Y(end - start), OY(end - start);
    for(i = start; i < end; i++)
      Y[i - start] = OY[i - start] = X[i];

    // finding Sdiff
    float Sdiff = findSdiff(OY, 0, OY.size());
    
    //    cerr << "bootstrapping for " << start << " " << end  << endl;

    for(i = 0; i < num_bootstraps; i++) {
      vector<pair<int, int> >::iterator first = bindexes.begin(),
	last = bindexes.end();
      
      random_shuffle(first, last);    // randomly generated an ordering of S 
      int blocks = 0;
      for(int l = 0, j = 0; j < bindexes.size(); j++) {
	//	cerr << "bindex " << j << " " << bindexes[j].first << " " << bindexes[j].first + bindexes[j].second - 1 << endl;
	blocks += bindexes[j].second;
	for(int k = 0; k < bindexes[j].second; k++) 
	  Y[l++] = OY[bindexes[j].first + k];
      }
      //      cerr << "blocks sum " << blocks << endl;
      float mySdiff = findSdiff(Y, 0, Y.size());
      //      cerr << Sdiff << " " << mySdiff << " " << endl;
      if(mySdiff < Sdiff)
	num_changed++;
      
    }
    //    cerr << "bootstrapping results " << num_changed << endl;
    return num_changed*100.0/num_bootstraps;    
  }

  int RSqUtils::getIntervalRSqChangePoint(int start, int end, 
					  vector<float> &X,
					  DENSE_HASH<int, int> *allowed_cps) {

    int m = start + 1;
    int min_m = -1;
    float min_MSE = DOUBLE_INFINITY;
    
    for( ; m < end - 1; m++) {
      if(allowed_cps) {
	DENSE_HASH<int, int>::iterator it = allowed_cps->find(m);
	if(it == allowed_cps->end())
	  continue;
      }

      float MSE = findMSE(X, start, m);
      MSE += findMSE(X, m, end);

      if(MSE < min_MSE) {
	min_MSE = MSE;
	min_m = m;
      }
    }
    return min_m;
  }

  void RSqUtils::computeRSqChangePoints(vector<float> &X,
					vector<RSqChangePoint> &cpoints,
					int num_bootstraps,
					int start, int end,
					int blockLen,
					double min_conf,
					vector<int> *positions,
					DENSE_HASH<int, int> *allowed_cps) {
					    
    if(positions && X.size() != positions->size()) {
      cerr << "Input vectors X and positions must have same length";
      return;
    }

    float conf = meanShiftConfidence(X, start, end, num_bootstraps, blockLen);
    //    cerr << "conf " << conf << " for " << start << " " << end << endl;
    if(conf > min_conf) {
      int min_m = getIntervalRSqChangePoint(start, end, X, allowed_cps);
					    
      if(start >= min_m || min_m + 1 >= end)
	return;

      //      cout << level << "change_point " << min_m << "(" << x_i << ") "  << conf << endl;      
      float left_mean = findMean(X, start, min_m);
      float right_mean = findMean(X, min_m, end);

      computeRSqChangePoints(X, cpoints, num_bootstraps, start, min_m, 
			     blockLen, min_conf, positions, allowed_cps);

      computeRSqChangePoints(X, cpoints, num_bootstraps, min_m + 1, end,
			     blockLen, min_conf, positions, allowed_cps);
      
      vector<RSqChangePoint>::iterator it = cpoints.size() ? 
	cpoints.begin() : cpoints.end();

      int x_i = positions ? (*positions)[min_m] : min_m;      
      while(it != cpoints.end() && it->pos < x_i)
	it++;
      
      cpoints.insert(it, RSqChangePoint(x_i, conf, left_mean, right_mean));

    }
  }


  void RSqUtils::refineRSqChangePoints(vector<float> &X,
				       vector<RSqChangePoint> &cpoints,
				       int num_bootstraps,
				       int start, int end,
				       int blockLen,
				       double min_conf,
				       vector<int> *positions,
				       DENSE_HASH<int, int> *allowed_cps) {
    
    if(positions && X.size() != positions->size()) {
      cerr << "Input vectors X and positions must have same length";
      return;
    }
    //    cout << "before refinement\n";
    //    displayRSqChangePoints(cpoints, start, end, positions);
    int A = start, i;
    for(i = 0; i <  cpoints.size(); i++) {
      int B = (i  + 1 < cpoints.size() ? cpoints[i + 1].pos : end);
      RSqChangePoint &cp = cpoints[i];
      // updating X[cp.pos]
      cp.conf = meanShiftConfidence(X, A, B, num_bootstraps, blockLen);
      cp.setMeans(findMean(X, A, cp.pos), findMean(X, cp.pos, B));    
      double min_mean = Utils::min(cp.left_mean, cp.right_mean);
      double variation = fabs(cp.right_mean - cp.left_mean);

      if(verbose)
	cerr << "signal " << (positions ? (*positions)[A] : A) << " " << (positions ? (*positions)[cp.pos] : cp.pos) << " " << (positions ? (*positions)[B] : B) << " " << " " << cp.left_mean << " " << cp.right_mean << " " << variation << " " << cp.conf;
      
      if(cp.conf <= min_conf) {
	// recalculate B
	int C = (i + 2 < cpoints.size() ? cpoints[i + 2].pos : end);	
	if(verbose)
	  cerr << " out " << cp.pos << " and " << B;
	if(B < C) {
	  int min_m = getIntervalRSqChangePoint(A, C, X, allowed_cps);

	  if(A < min_m && min_m + 1 < C) {
	    cpoints[i + 1].pos = min_m;
	    cp.right_mean = findMean(X, cp.pos, cpoints[i + 1].pos);
	    if(verbose)
	      cerr << " for " << min_m;
	  }
	}
	
	vector<RSqChangePoint>::iterator it = cpoints.begin() + i;
	it = cpoints.erase(it);
	i--;	
      }

      if(verbose)
	cerr << endl;
      A = i < 0 ? start : cpoints[i].pos;
    }
    //    cout << "after refinement\n";
    //    displayRSqChangePoints(cpoints, start, end, positions);
    
    for(i = 0; i <  cpoints.size(); i++) {
      RSqChangePoint &cp = cpoints[i];     
      if(positions)
	cp.pos = (*positions)[cp.pos];
    }
  }


  bool RSqUtils::isRealOnset(RSqChangePoint *cp, 
			     vector<RSqChangePoint> &changes,
			     int start, int end,
			     float ovmean_l, float ovmean_r,
			     int level) {
    int i = 0;
    double min_dws = DOUBLE_INFINITY, max_ups = 0, right_cov = -1;
    bool upstream = true, right_increasing = true;
    
    for(; i < changes.size(); i++) {
      RSqChangePoint &tcp = changes[i];
      if(tcp.pos < start || tcp.pos > end)
	continue;
      
      if(upstream)
	max_ups = (tcp.left_mean > max_ups ? tcp.left_mean : max_ups);
      else
	min_dws = (tcp.right_mean < min_dws ? tcp.right_mean : min_dws);
      
      if(right_cov >= 0) {
	cout << "onset previous " << right_cov << " new " << tcp.right_mean << endl;
	if(right_increasing && 0.8*right_cov < tcp.right_mean) 
	  right_cov = tcp.right_mean;
	else
	  right_increasing = false;
      }
      if(tcp.pos == cp->pos) {
	upstream = false;
	right_cov = cp->right_mean;
      }
    }
    
    cout << "onset max_UPS " << max_ups << " right " << cp->right_mean << " min_DWS " << min_dws << " left " << right_cov << "(" << cp->right_mean << ")";
    
    if(max_ups >= right_cov || !min_dws || min_dws <= cp->left_mean) {
      cout << " false" << endl;
      return false;
    }
    cout << " true" << endl;
    return true;
  }

  bool RSqUtils::isRealOffset(RSqChangePoint *cp, 
			      vector<RSqChangePoint> &changes,
			      int start, int end,
			      float ovmean_l, float ovmean_r,
			      int level, bool gene_dws) {
    int i = changes.size() - 1;
    double max_dws = 0, min_ups = DOUBLE_INFINITY, left_cov = -1;
    bool upstream = false, left_increasing = true;

    for( ; i >= 0; i--) {
      RSqChangePoint &tcp = changes[i];
      if(tcp.pos < start || tcp.pos > end)
	continue;

      if(upstream)
	min_ups = (tcp.left_mean < min_ups ? tcp.left_mean : min_ups);
      else
	max_dws = (tcp.right_mean > max_dws ? tcp.right_mean : max_dws);
      
      if(left_cov >= 0) {
	//	cout << "offset previous " << left_cov << " new " << tcp.left_mean << endl;
	if(left_increasing && 0.8*left_cov < tcp.left_mean)
	  left_cov = tcp.left_mean;
	else
	  left_increasing = false;
      }
      if(tcp.pos == cp->pos) {
	upstream = true;
	left_cov = tcp.left_mean;
      }
    }

    cout << "offset min_UPS " << min_ups << " right " << cp->right_mean << " max_DWS " << max_dws << " left " << left_cov << "(" << cp->left_mean << ")"; 
          
    if(!min_ups || min_ups <= cp->right_mean || 
       (!gene_dws && max_dws >= left_cov)) {
      cout << " false\n";
      return false;
    }
    cout << " true\n";
    return true;
  }


  void RSqUtils::changePoints2Histogram(vector<RSqChangePoint> &cps,
					int start, int end, int special_pos,
					char special_char, int offset) {
					
    int max = 0;
    int i;
    cout << "| ";
    for(i = 0; i < cps.size(); i++) {
      int pos = cps[i].pos + offset;
      if(pos < start || pos > end)
	continue;

      if(cps[i].right_mean > max)
	max = (int)(cps[i].right_mean + 1);
      cout << pos << " ";
    }
    cout << endl;
    
    if(!max) return;
    
    for(int j = 2*max; j > 0; j--) {
      cout << "| ";
      for(i = 0; i < cps.size(); i++) {
	int pos = cps[i].pos + offset;
	if(pos < start || pos > end)
	  continue;
	int p = pos, w = 0;
	for(; p > 0; w++)
	  p = p/10;
	for(int k = 0; k < w; k++) {
	  if(j == 2*max && pos == special_pos)  cout << special_char;
	  else if(2*cps[i].right_mean >= j)  cout << "*";
	  else  cout << " ";
	}
	cout << " ";
      }
      cout << endl;
    }
  }

  void RSqUtils::displayRSqChangePoints(vector<RSqChangePoint> &cpoints,
					int start, int end,
					vector<int> *positions) {
    
    int A = start, i;
    for(i = 0; i <  cpoints.size(); i++) {
      RSqChangePoint &cp = cpoints[i];
      int C = (i < cpoints.size() - 1 ? cpoints[i + 1].pos : end - 1);
      int B = cpoints[i].pos;
      if(positions) {
	A = (*positions)[A];
	B = (*positions)[B];
	C = (*positions)[C];
      }
      
      cout << "cp [" << A << ", " << B << ", " 
      	   << C << "] = " << cp.left_mean << " "
	   << cp.right_mean << " " << cp.conf  << endl;
      A = i < 0 ? start : cp.pos;
    }
  }


  void RSqUtils::selectOnsetPoints(vector<RSqChangePoint> &changes, 
				   int start, int end,
				   vector<float> &X, 
				   double left_ovmean, double right_ovmean,
				   double min_conf,
				   RSqChangePoint &onset) {

    
    int i;
    float max_score = 0;
    vector<RSqChangePoint>::iterator it;
    
    for(i = 0; i < changes.size(); i++) {
      RSqChangePoint &cp = changes[i];
      if(cp.pos < start || cp.pos > end)
	continue;
      if(cp.conf < min_conf || cp.left_mean > cp.right_mean)
	continue;

      //      float ldiff = cp.left_mean - left_ovmean, 
      //	rdiff = right_ovmean - cp.right_mean;
      //      float mean_diff = pow(ldiff > 0 ? ldiff : 0, 2.0) +
      //	pow(rdiff > 0 ? rdiff : 0, 2.0);
      
      float ups_mean = (start < 0 ? findMean(X, 0, cp.pos) : 
			findMean(X, start, cp.pos));
      float ds_mean =  (end > X.size() ? findMean(X, cp.pos, X.size()) :
			findMean(X, cp.pos, end));
      float score = ds_mean - ups_mean;
      if(verbose)
	cerr << "signal onset " << cp.pos << " " << cp.conf << " " 
	     << cp.left_mean << " " << cp.right_mean << " " << " " << ups_mean << " " << ds_mean << " " << score << endl;

      //      it = onsets.size() ? onsets.begin() : onsets.end();
      //      while(it != onsets.end() && it->score > score)
      //       	it++;
	
      //      it = onsets.insert(it, cp);
      //      it->score = score;

      if(max_score < score) { // && cp.pos > 1) {
      	max_score = score;
	onset = cp;
      }
    }

    onset.score = max_score;    
    if(max_score && verbose)
      cerr << "best signal onset " << onset.pos << " " 
	   << onset.conf << " " << onset.left_mean << " " 
	   << onset.right_mean << " " << max_score << endl;
    
  }

  void RSqUtils::selectOffsetPoints(vector<RSqChangePoint> &changes,
				    int start, int end,
				    vector<float> &X,
				    double left_ovmean, double right_ovmean,
				    double min_conf,
				    RSqChangePoint &offset) {
    
    int i;
    float max_score = 0;
    vector<RSqChangePoint>::iterator it;

    for(i = 0; i < changes.size(); i++) {
      RSqChangePoint &cp = changes[i];
      if(cp.pos < start || cp.pos > end)
	continue;
      if(cp.conf < min_conf || cp.left_mean < cp.right_mean)
	continue;

      //      float ldiff = left_ovmean - cp.left_mean,
      //	rdiff = cp.right_mean - right_ovmean;
      //      float mean_diff = pow(ldiff > 0 ? ldiff : 0, 2.0) +
      //	pow(rdiff > 0 ? rdiff : 0, 2.0);

      float ups_mean = (start < 0 ? findMean(X, 0, cp.pos) : 
			findMean(X, start, cp.pos));
      float ds_mean =  (end > X.size() ? findMean(X, cp.pos, X.size()) :
			findMean(X, cp.pos, end));
      float score = ups_mean - ds_mean;

      if(verbose)
	cerr << "signal offset " << cp.pos << " " << cp.conf << " " 
	     << cp.left_mean << " " << cp.right_mean << " " << " " << ups_mean << " " << ds_mean << " " << score << endl;
      
      //      it = offsets.size() ? offsets.begin() : offsets.end();
      
      //      while(it != offsets.end() && it->score > score)
      //      	it++;
      
      //      it = offsets.insert(it, cp);
      //      it->score = score;

      if(max_score < score) { // && cp.pos < X.size()) {
	max_score = score;
	offset = cp;
      }
    }

    offset.score = max_score;
    if(max_score && verbose)
      cerr << "best signal offset " << offset.pos << " " 
	   << offset.conf << " " << offset.left_mean << " " 
	   << offset.right_mean << " " << max_score << endl;

  }

  void RSqUtils::changePoints2Filter(int offset,
				     string &symbols, 
				     vector<RSqChangePoint> &changes,
				     vector<float> &X, 			
				     double min_conf,	
				     int clen,
				     vector<char> &sig_filter,
				     vector<vector<int> > &sigscore_filter,
				     TStrand strand) {
    
    int i;
    int start = 0, end = X.size();
    for(i = 0; i < changes.size(); i++) {
      RSqChangePoint &cp = changes[i];
      if(cp.pos > clen)
	assert(0);
      
      if(cp.conf < min_conf)
	continue;
      int pos = cp.pos;
      int pos_l = (i == 0 ? start : changes[i - 1].pos);
      int pos_c = pos;
      int pos_r = (i == changes.size() - 1 ? end : changes[i + 1].pos);
      float mean_l = (pos_l >= 1 ? findMean(X, pos_l, pos_c) : 0.1);
      float mean_r = (pos_r <= clen ? findMean(X, pos_c, pos_r) : 0.1);
      
      int s = (int)strand;

      char symbol = symbols[s*2];
      if(cp.left_mean > cp.right_mean) {
	symbol = symbols[s*2 + 1];
	pos = (cp.pos <= 1 ? cp.pos : cp.pos - 1);
	pos_c = pos;
      }
      
      if(strand == STRAND_COMP) {
	pos = clen - pos + 1;
	int pos_tmp = pos_l;
	pos_l = clen - pos_r + 2;
	pos_c = clen - pos_c + 2;
	pos_r = clen - pos_tmp + 2;
      }
      
      sig_filter[pos + offset] = symbol;
      sigscore_filter[pos + offset] = vector<int>(3);
      sigscore_filter[pos + offset][0] = pos_l + offset;
      sigscore_filter[pos + offset][1] = pos_c + offset;
      sigscore_filter[pos + offset][2] = pos_r + offset;

    }
  }


  // Even for stranded RNA-Seq the complementary strand will be in InOneStrand
  // format. This means the offsets for the complementary strand will need to
  // be complemented before assignment to the filter
  void RSqUtils::onsets2Filter(vector<RSqChangePoint> &changes, char symbol,
			       int start, int end,
			       vector<char> &sig_filter,
			       vector<vector<int> > &sigscore_filter,
			       double left_ovmean, double right_ovmean,
			       double min_conf,
			       int clen,
			       vector<float> &X, TStrand strand) {
    
    int i;
    for(i = 0; i < changes.size(); i++) {
      RSqChangePoint &cp = changes[i];
      if(cp.pos > clen)
	assert(0);
      
      if(cp.conf < min_conf)
	continue;
	 
      int pos = cp.pos;
      int pos_l = (i == 0 ? start : changes[i - 1].pos);
      int pos_c = pos;
      int pos_r = (i == changes.size() - 1 ? end : changes[i + 1].pos);
      float mean_l = (pos_l >= 1 ? findMean(X, pos_l, pos_c) : 0.1);
      float mean_r = (pos_r <= clen ? findMean(X, pos_c, pos_r) : 0.1);

      if(cp.left_mean > cp.right_mean)
	continue;
      
      if(strand == STRAND_COMP) {
	pos = clen - pos + 1;
	int pos_tmp = pos_l;
	pos_l = clen - pos_r + 2;
	pos_c = clen - pos_c + 2;
	pos_r = clen - pos_tmp + 2;
      }
      
      if(verbose)
	cerr << pos << " " << pos_l << " " << pos_c << " " << pos_r;

      sig_filter[pos] = symbol;
      sigscore_filter[pos] = vector<int>(3);
      sigscore_filter[pos][0] = pos_l;
      sigscore_filter[pos][1] = pos_c;
      sigscore_filter[pos][2] = pos_r;

      if(verbose) 
	if((strand == STRAND_FWD && (cp.left_mean != mean_l || cp.right_mean != mean_r)) ||
	   (strand == STRAND_COMP && (cp.left_mean != mean_r || cp.right_mean != mean_l)))
	  cerr << " bad onset " << cp.left_mean << " " << mean_l
	       << " and " << cp.right_mean << " " << mean_r << " " << strand << endl;
	else cerr << "\n";      
      
    }

  }

  // Even for stranded RNA-Seq the complementary strand will be in InOneStrand
  // format. This means the offsets for the complementary strand will need to
  // be complemented before assignment to the filter
  void RSqUtils::offsets2Filter(vector<RSqChangePoint> &changes, char symbol,
				int start, int end,
				vector<char> &sig_filter,
				vector<vector<int> > &sigscore_filter,
				double left_ovmean, double right_ovmean,
				double min_conf,
				int clen, 
				vector<float> &X, TStrand strand) {
    

    int i;
    for(i = 0; i < changes.size(); i++) {
      RSqChangePoint &cp = changes[i];
      if(cp.pos > clen)
	assert(0);
      
      if(cp.conf < min_conf)
	continue;

      int pos = (cp.pos <= 1 ? cp.pos : cp.pos - 1);
      int pos_l = (i == 0 ? start : changes[i - 1].pos);
      int pos_c = pos;
      int pos_r = (i == changes.size() - 1 ? end : changes[i + 1].pos);
      float mean_l = (pos_l >= 1 ? findMean(X, pos_l, pos_c) : 0.1);
      float mean_r = (pos_r <= clen ? findMean(X, pos_c, pos_r) : 0.1);
      
      //      if(cp.left_mean < cp.right_mean)
      if(mean_l < mean_r)
	continue;

      if(strand == STRAND_COMP) {
 	pos = clen - pos + 1;
	int pos_tmp = pos_l;
	pos_l = clen - pos_r + 2;
	pos_c = clen - pos_c + 2;
	pos_r = clen - pos_tmp + 2;
      }

      if(verbose)
	cerr << pos << " " << pos_l << " " << pos_c << " " << pos_r;

      sig_filter[pos] = symbol;
      sigscore_filter[pos] = vector<int>(3);
      sigscore_filter[pos][0] = pos_l;
      sigscore_filter[pos][1] = pos_c;
      sigscore_filter[pos][2] = pos_r;

      if(verbose)
	if((strand == STRAND_FWD && (cp.left_mean != mean_l || cp.right_mean != mean_r)) ||
	   (strand == STRAND_COMP && (cp.left_mean != mean_r || cp.right_mean != mean_l))) 
	  cerr << " bad offset " << cp.left_mean << " " << mean_l
	       << " and " << cp.right_mean << " " << mean_r << " " << strand << endl;
	else cerr << "\n";
    }
  }
  
  void RSqUtils::filter2XFastaFile(vector<char> &filter,
				   ofstream & fd) {

    char last_char = filter[1];
    int times = 1, pos;
    
    for(int i = 2; i < filter.size(); i++) {
      if(filter[i] != last_char) {
	pos =  i - times;
	if(times == 1)
	  fd << last_char;
	else
	  fd << last_char << " x " << times << " " << pos;
	fd << "\n";
	last_char = filter[i];
	times = 1;
      }
      else
	times++;
    }
    
    pos =  filter.size() - times;
    fd << last_char << " x " << times << " " <<  pos << "\n";
    
  }

  void RSqUtils::filter2MultiScoreFile(vector< vector<int> > &filter,
				       int num_cols,
				       ofstream & fd) {
        
    for(int i  = 1; i < filter.size(); i++)
      if(filter[i].size()) {
	fd << i;
	for(int j = 0; j < filter[i].size(); j++)
	  fd << "\t" << filter[i][j];
	fd << endl;
      }
  }
  
  void RSqUtils::storeXcriptSignals(Sequence &c,
				    vector<RSqChangePoint> & pchanges,
				    double intron_mean, double exon_mean,
				    double min_conf,
				    std::ofstream &sigd,
				    std::ofstream &sigscored,
				    vector<float> &X,
				    TStrand strand) {
        
    string symbols("XZZX");
    vector<char> filter;
    vector<vector<int> > score_filter;
    int s = (int)strand;
    filter = vector<char>(c.length() + 1, '-');
    score_filter = vector<vector<int> >(c.length() + 1);
    RSqUtils::onsets2Filter(pchanges, symbols[s*2],
			    0, X.size(),
			    filter, score_filter,
			    intron_mean, exon_mean, min_conf,
			    c.length(), X, strand);
    RSqUtils::offsets2Filter(pchanges, symbols[s*2 + 1],
			     0, X.size(),
			     filter, score_filter,
			     exon_mean, intron_mean, min_conf,
			     c.length(), X, strand);
    
    filter2XFastaFile(filter, sigd);   
    filter2MultiScoreFile(score_filter, 3, sigscored);
  }
  
  /**
   * Converted from MATLAB script at http://billauer.co.il/peakdet.html    
     Currently returns two lists of tuples, but maybe arrays would be better
     function [maxtab, mintab]=peakdet(v, delta, x)
     %PEAKDET Detect peaks in a vector
     %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
     %        maxima and minima ("peaks") in the vector V.
     %        MAXTAB and MINTAB consists of two columns. Column 1
     %        contains indices in V, and column 2 the found values.
     %      
     %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
     %        in MAXTAB and MINTAB are replaced with the corresponding
     %        X-values.
     %
     %        A point is considered a maximum peak if it has the maximal
     %        value, and was preceded (to the left) by a value lower by
     %        DELTA.
     
     % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
     % This function is released to the public domain; Any use is allowed.
  */    
  void RSqUtils::find_peaks(vector<float> & X,
			    vector<pair<int, float> > &min,
			    vector<pair<int, float> > &max, 
			    float mn_delta, float mx_delta, 
			    vector<int> *positions) {
    
    if(positions && positions->size() != positions->size()) {
      cerr << "Input vectors X and positions must have same length";
      return;
    }
    
    if(mx_delta <= 0 || mn_delta <= 0)  {
      cerr << "Input argument delta must be positive";
      return;
    }
    
    float mn = DOUBLE_INFINITY, mx = -DOUBLE_INFINITY;
    int mnpos = -1, mxpos = -1;
    bool lookformax = true;
    
    for(int i = 0; i < X.size(); i++) {
      int x_i = positions ? (*positions)[i] : i;
      float this_x = X[i];
      if(this_x > mx) {
	mx = this_x;
	mxpos = x_i;
      }
      if(this_x < mn) {
	mn = this_x;
	mnpos = x_i;
      }
      
      if(lookformax) {
	if(this_x < mx - mx_delta) {
	  max.push_back(pair<int, float>(mxpos, mx));
	  mn = this_x;
	  mnpos = x_i;
	  lookformax = false;
	}
      }
      else {
	if(this_x > mn + mn_delta) {
	  min.push_back(pair<int, float>(mnpos, mn));
	  mx = this_x;
	  mxpos = x_i;
	  lookformax = true;
	}
      }
    }
  }
  
  void RSqUtils::computeRamps(vector<float> &cov,
			      int window, 
			      float score_cutoff,
			      vector<pair<int, float> > &onsets,
			      vector<pair<int, float> > &offsets,
			      vector<int> *positions) {
    
    float leftAvg, rightAvg, thisAvg;
    vector<float> ronsets, roffsets;
    vector<int> ponsets, poffsets;
    vector<pair<int, float> > min;

    int i, j, leftB = 0, rightB = 2*window;
    int size = cov.size() - window;

    // onsets
    for(i = window + 1; i < size; i++) {
      leftAvg = 0;
      for(j = leftB; j >= 0; j--) {
	thisAvg = (cov[i - 1] - cov[j])/(window + leftB - j);
	if(thisAvg > leftAvg)
	  leftAvg = thisAvg;
      }
      int x_i = positions ? (*positions)[i - 1] : i - 1;
      cout << "onset " << x_i << " " << window + leftB - j << " " << leftAvg;
      rightAvg = DOUBLE_INFINITY;
      for(j = rightB; j < size; j++) {
	thisAvg = (cov[j] - cov[i - 1])/(window + j - rightB);
	if(thisAvg < rightAvg)
	  rightAvg = thisAvg;
      }
      ronsets.push_back(leftAvg ? rightAvg/leftAvg : 0);
      cout << " " << window + j - rightB << " " << rightAvg << " " << ronsets.back() << endl;
      ponsets.push_back(positions ? (*positions)[i] : i);
      leftB++; rightB++;
    }

    // offsets
    leftB = 0, rightB = 2*window;
    size = cov.size() - window - 1;
    for(i = window; i < size; i++) {
      leftAvg = DOUBLE_INFINITY;
      for(j = leftB; j >= 0; j--) {
        thisAvg = (cov[i] - cov[j])/(window + leftB - j);
        if(thisAvg < leftAvg)
          leftAvg = thisAvg;
      }
      int x_i = positions ? (*positions)[i] : i;
      cout << "offset " << x_i << " " << window + leftB - j << " " << leftAvg;
      rightAvg = 0;
      for(j = rightB; j < size; j++) {
	thisAvg = (cov[j] - cov[i])/(window + j - rightB);
	if(thisAvg > rightAvg)
	  rightAvg = thisAvg;
      }
      roffsets.push_back(rightAvg ? leftAvg/rightAvg : 0);
      cout << " " << window + j - rightB << " " << rightAvg << " " << roffsets.back() << endl;
      poffsets.push_back(positions ? (*positions)[i] : i);
      leftB++; rightB++;
    }
    
    find_peaks(ronsets, min, onsets, 0.1, score_cutoff, &ponsets);      
    min.clear();
    find_peaks(roffsets, min, offsets, 0.1, score_cutoff, &poffsets);
  }
   
}
