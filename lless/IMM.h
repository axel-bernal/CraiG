/****************************************************************************
* IMM.h - part of the lless namespace, a general purpose
*         linear semi-markov structure prediction library
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

#ifndef _MMSEQ_H
#define _MMSEQ_H
  
#include <streambuf>
#include <string.h>
#include "NGram.h"
#include "Utils.h"
#include "Sigma.h"
#include <math.h>
#include <assert.h>


namespace lless {  

  /**
   * MM is a base class for implementation of a inhomogeneous
   * markov chain model.
   */
  class MM {
   protected:
      int buildComplemen;
      Sigma *alphabet;
      int numNGrams;
      NGram **ngrams;

   public:
      MM() { ngrams = NULL;}

      MM(Sigma *alphabet, 
         int n = 2, 
         int np = 3, 
         int buildC = 0);    

      MM(double  ***freqMat, 
         Sigma *alphabet = NULL, 
         int n = 2, 
         int np = 3, 
         int buildC = 0);

      inline int compModel() {
        return buildComplemen - 1;
      }
  
      inline int numPos() {
        return numNGrams;
      }
  
      inline NGram * ngram(int num) {
        return ngrams[num];
      }

      ~MM();
          
   protected:
      double getambCp(char *, int, int);
      double getavgCp(char *, int, int);
      void readTable(::ifstream *);
  
  };
  
  /**
   * FMM is a subclass of MM which implements a normal inhomogeneous
   * markov chain of order n and period np.
   */
  class FMM : public MM {
   public:
    FMM(char *fdes = NULL, 
        Sigma *alphabet = NULL, 
        int n=2,
        int np=3, 
        int buildComp = 0);

    inline void buildCp(int order) {
      ngrams[0]->buildCp();
      ngrams[order]->buildCp();
    }

    ~FMM();

  };
  

  /**
   * The IMM class implements an Interpolated Markov Model. Probabilities are
   * computed as in Salzberg's glimmer paper.
   **************************************************************************/ 

  class IMM : public MM {
   private:
    double ****lambda; 
    unsigned int thresholdFreq;
    double ***interpolatedProbs;
    ChiSquare *chisqr;
   public:
    IMM(char *file, 
        Sigma *alphabet, 
        int n, 
        int np=3, 
        unsigned int maxf=400, 
        int buildComp = 0);

    IMM(::ifstream &, Sigma *alphabet);

    IMM(double ***freqMat, 
        Sigma *alphabet, 
        int n=2, 
        int np=3, 
        unsigned int maxf=400, 
        int buildComp = 0);
    
    void computeInterpolatedProbs(int numNGram);
    void emitSequence(char *seq, int frame, int length);
    void probeAllNGramProbs();
    void buildCp(void);
    double imm(unsigned int, int, char *, unsigned int, int);
    double seqScore(char *, int = 1, int initPos = 1, int len = 0);
    void seqScore(double *, char *, int = 1, int initPos = 1, int len = 0);
    void storeModel(ofstream &);
    void deleteLambdas();
    void deleteInterpolatedProbs();
    ~IMM();
  };
  
}
  
#endif
