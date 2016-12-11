/****************************************************************************
* FT_HiddenSeq.h - part of the lless namespace, a general purpose
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

#ifndef F_HIDDEN_SEQ_FEAT_H
#define F_HIDDEN_SEQ_FEAT_H
  
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "FT_PeriodicCountingSegment.h"
#include "NGram.h"
#include "Sequence.h"
#include "Feature.h"
#include "FT_Conjunction.h"

namespace lless {  

  /**
   * The FT_HiddenSeq class  is used to comprehensively account n-grams which
   * are shared by two non-contiguos states.
   * The problem with this implementation is that there is need to 
   * access directly the parameter vectors of the features that count n-grams,
   * because of performance reason (this feature is heavily called throughout 
   * the decoding phase.
   **************************************************************************/
  class FT_HiddenSeq : public TypedFeature<UCHAR> {
  
   protected:
    int maxCodingOrder;
    int period;
    char *hiddenSeq;
    Tag *lastTag;
    int lastEnd;
    int lastFrame;
    XORConsensus *_forbseq;
    vector<FT_CountingSegment<UCHAR> *> features;
  
   public:
    FT_HiddenSeq(
		 int fInd, 
		 int paramInd, 
		 int parsingFrames, 
		 char *name, 
		 FilterEngine *fe, 
		 int maxOrder, 
		 int period, 
		 string signals,          /*!< A string representing the signal
					    that needs to be found, as it appears
					    in the sequence. */
		 Sigma *sigma           //!< An alphabet for signals.
		 )
      : TypedFeature<UCHAR>(fInd, paramInd, parsingFrames, 
			    name, fe, 0, FT_UCHAR) {
      
      maxCodingOrder = maxOrder;
      this->period = period;
      hiddenSeq = new char [2*maxOrder + 4];
      hiddenSeq[0] = '\0';
      lastFrame = -1;
      lastTag = NULL;

      boost::RegEx rExDelim("\\s+");
      std::vector<std::string> vsignals;
      rExDelim.Split(vsignals, signals, 0, 100);
      assert(vsignals.size());
      this->_forbseq = new XORConsensus((char *)vsignals[0].c_str(), sigma);
      for(unsigned int j = 1; j < vsignals.size(); )
        _forbseq->matchPattern((char *)vsignals[j++].c_str());

    }
    
    FT_HiddenSeq(
		 int fInd, 
		 int paramInd, 
		 vector<std::string> & fargs, 
		 int & offset, 
		 FilterEngine *fe, 
		 FeatureEngine *fte
		 )
      : TypedFeature<UCHAR>(
			  fInd, paramInd,
			  fargs, offset, 
			  fe, 0, FT_UCHAR) {
  
      if(!sscanf(fargs[offset++].c_str(), "%d", &maxCodingOrder))
        assert(0);
      if(!sscanf(fargs[offset++].c_str(), "%d", &period))
        assert(0);
  
      int numSignals, numFeaturesToAdd;

      if(!sscanf(fargs[offset++].c_str(), "%d", &numSignals))
        assert(0);
      ResourceEngine *re = fe->getResourceEngine();
      Sigma *sigma  = (Sigma *)re->getResource(fargs[offset + numSignals]);
      this->_forbseq = new XORConsensus((char *)fargs[offset++].c_str(), sigma);
      for(int i = 1; i < numSignals; i++)
        _forbseq->matchPattern((char *)fargs[offset++].c_str());
      offset++; // for sigma

      if(!sscanf(fargs[offset++].c_str(), "%d", &numFeaturesToAdd))
        assert(0);
  
      for(int i = 0; i < numFeaturesToAdd; i++)
        addFeature((FT_CountingSegment<UCHAR> *)fte->getFeature(fargs[offset++]));
  
      hiddenSeq = new char [2*maxCodingOrder + 4];
      hiddenSeq[0] = '\0';
      lastFrame = -1;
      lastTag = NULL;
    }
  
  
    inline void addFeature(FT_CountingSegment<UCHAR> * f) {
      features.push_back(f);
    }
  
    /**
     * Computes the hidden sequence between two exons.
     * ge->getPos() is the beginning of the intron.
     * ge->getPos() + ge->getLen() is the beginning of the next Exon.
     * So the hidden sequence is going to be located joining parts of the end
     * of the last exon and the beginning of the next exon
     */
    inline void computeHiddenSeq(Tag *ge, int frame) {
      computeHiddenSeq(ge, ge->getPos(), ge->getPos() + ge->getLen() - 1, frame);
    }

    inline void computeHiddenSeq(Tag *ge, int beg, int end, int frame) {
      //      cerr << endl << ge << " " << end << " " << lastTag << " " << lastEnd;
      if(ge == lastTag && lastEnd == end && frame == lastFrame)
        return;

      int i;
  
      for(i = 0; i < 2*maxCodingOrder + 3; i++)
        hiddenSeq[i] = 'N';
  
      char *middle = hiddenSeq + maxCodingOrder; // the codon in the middle
      char *seq = this->fe->getSeq(ge->getStrand());
      int seqLen = strlen(seq);
  
      if(beg - frame >= 1 && end - frame + 3 <= seqLen) {
        if(beg - frame - maxCodingOrder >= 1)
          strncpy(hiddenSeq, seq + beg - frame - maxCodingOrder - 1, maxCodingOrder);
        
        strncpy(middle, seq + beg - frame -1, frame);
        strncpy(middle + frame, seq + end, 3 - frame);
        if(end - frame + 3 + maxCodingOrder <= seqLen) 
          strncpy(middle + 3, seq + end + 3 - frame, maxCodingOrder);
  
        hiddenSeq[2*maxCodingOrder + 3] = '\0';
      }
      else
        hiddenSeq[0] = '\0';

      lastTag = ge;
      lastEnd = end;
      lastFrame = frame;
    }
    
    inline bool hiddenSeqHasStop() {
  
      if(!hiddenSeq || !strlen(hiddenSeq))
        return false;
      //      cerr << hiddenSeq << " " << maxCodingOrder << endl;
      return _forbseq->matches(hiddenSeq + maxCodingOrder, _forbseq->size());
    }
  

    inline double dotParamV(Tag *ge, int frame) {
      if(this->isOff())
        return 0;

      if(!frame)
        return 0;

      computeHiddenSeq(ge, frame);
      
      if(!strlen(hiddenSeq) || hiddenSeqHasStop())
        return 0;
  
      FL_Gram<UCHAR> *gram;
      Sigma *sigma;
      double  res = 0;
      int histLen = maxCodingOrder - 1;
  
      for(unsigned int ind = 0; ind < features.size(); ind++) {

        FT_CountingSegment<UCHAR> *f = features[ind];
  
        gram = (FL_Gram<UCHAR> *)f->getFilter();
        sigma = gram->alphabet();
  
        for(int i = 3 + histLen - gram->order() + 1 + f->phase(); i < histLen + frame; i += period) {
          int kmer = sigma->randIndS(hiddenSeq, i, i + gram->order() - 1);
          res += f->posDotParamV((UCHAR &)kmer, ge->getPos(),
				 ge->getStrand(),
				 this->cparams[f->paramInd()]);
	}
      }
      //    cerr << "hidden " << hiddenSeq << "@" << ge->getPos() << " = " << res << endl;
      return res;
    }
    
    inline void updParamV(double updVal, Tag *ge, int frame) {
      if(this->frozen())
        return;
  
      if(!frame)
        return;
  
      computeHiddenSeq(ge, frame);
  
      if(!strlen(hiddenSeq)|| hiddenSeqHasStop())
        return;
  
      FL_Gram<UCHAR> *gram;
      Sigma *sigma;
      int histLen =  maxCodingOrder - 1;
  
      for(unsigned int ind = 0; ind < features.size(); ind++) {
	FT_CountingSegment<UCHAR> *f = (FT_CountingSegment<UCHAR> *)features[ind];	
	gram = (FL_Gram<UCHAR> *)f->getFilter();
        sigma = gram->alphabet();
	
        for(int i = 3 + histLen - gram->order() + 1 + f->phase(); i < histLen + frame; i += period) {
          int kmer = sigma->randIndS(hiddenSeq, i, i + gram->order() - 1);
	  f->posUpdParamV(updVal, (UCHAR &)kmer, ge->getPos(),
			  ge->getStrand(),
			  this->params[f->paramInd()]);
        }
      }
    }
  
    ~FT_HiddenSeq() {
  
      if(hiddenSeq)
        delete [] hiddenSeq;
      hiddenSeq = NULL;
      delete _forbseq;
    }
  };

}
  
#endif
