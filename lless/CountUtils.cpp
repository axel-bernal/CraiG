#include "CountUtils.h"
#include "Utils.h"

/****************************************************************************
* CountUtils.cpp - part of the lless namespace, a general purpose
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



namespace lless {
  
  /** 
   * This count grams of orders minOrder to maxOrder, but it does not use the
   * NGram class because of the fractional counts introduced by ambiguous 
   * characters which give more precise countings.
   */

  void CountUtils::countNGrams(double ****gramCounts, 
                               Sigma *sigma,
                               int numPos, 
                               int minOrder, 
                               int maxOrder,
                               char *seq, 
                               UCHAR *contextVals,
                               TStrand strand) {
    
    assert(contextVals && seq);
    double ***im;
    int seqLen = strlen(seq);

    for(int i = 0; i <= maxOrder-minOrder; i++) {
      int lenFL_Gram = i + minOrder;
      int beg, frame, numFL_Grams;
      int end;
      char *gram =  new char[lenFL_Gram + 1];
      int numPlets = (int)Utils::pow2(sigma->alphabetSize(), lenFL_Gram);

      if(lenFL_Gram > seqLen)
        return;

      for(beg = 0, end = lenFL_Gram - 1, frame = ((lenFL_Gram - 1) % numPos) + 1;end < seqLen; beg++,end++,frame = ((frame % numPos) + 1)) {
        /*
         * choose the gramCount matrix corresponding to current context
         */
        im = gramCounts[contextVals[beg + 1]];
        im[i][numPlets][0]++;
        im[i][numPlets][frame]++;
        strncpy(gram, seq + beg, lenFL_Gram);
        gram[lenFL_Gram] = '\0';
        numFL_Grams = sigma->lenIndS(gram, lenFL_Gram);
        /*
         * Compute count fractions if there are ambiguous characters in the 
         * sequence
         */
        {
          int *bb = new int[numFL_Grams];
          sigma->indS(bb, gram, lenFL_Gram);

          for(int j = 0; j < numFL_Grams; j++) {
            im[i][bb[j]][0] += 1.00/numFL_Grams;
            im[i][bb[j]][frame] += 1.00/numFL_Grams;
          }

          delete [] bb;
        }
      }

      delete [] gram;
    }                                            
  }

  /**
   * allocates and initializes a matrix for storing relative word counts
   */
  int CountUtils::initkcodonMatrix(double ****im, int numPos, int minPlets, int maxPlets, int bC, Sigma *sigma) {
    (*im) = new double ** [maxPlets-minPlets+1];

    if(!*im) return -1;

    for(int i=0;i<=maxPlets-minPlets;i++) {
      int numPlets = (int)Utils::Utils::pow2(sigma->alphabetSize(),i+minPlets);
      (*im)[i] = new double * [numPlets+1];
      if(!(*im)[i]) return -1;
      for(int j=0;j < numPlets + 1;j++) {
        (*im)[i][j] = new double[numPos*(bC+1) + 1];
        if(!(*im)[i][j]) return -1;
        for(int k=0;k < numPos*(bC+1) + 1;k++) {
        (*im)[i][j][k] = 0;
        }
      }
    }
    return 1;
  }                                  

  /*
   * Frees relative word count matrix
   */
  void CountUtils::freekcodonMatrix(double ****im, int minPlets, int maxPlets, Sigma *sigma) {
    for(int i=0;i<=maxPlets-minPlets;i++) {
      int numPlets = (int)Utils::pow2(sigma->alphabetSize(),i+minPlets);
      for(int j=0;j < numPlets + 1;j++)
        delete [] (*im)[i][j];
      delete [] (*im)[i];
    }
    delete [] *im;
  }             

  //! Computes word counts for gram of order lenkCodon
  
  void CountUtils::fillbaseTable(double **imi, int numPos, int lenkCodon, char *seq, int initPos, Sigma *sigma) {
    int beg, frame, numkCodons;
    unsigned int end , len = strlen(seq);
    char *kCodon =  new char[lenkCodon + 1], *s = seq;
    int numPlets = (int)Utils::pow2(sigma->alphabetSize(), lenkCodon);

    if((unsigned int)lenkCodon > len)
      return;

    for(beg = 0, end = lenkCodon - 1, frame = ((initPos + lenkCodon - 2) % numPos) + initPos;end < len;beg++,end++,frame = ((frame % numPos)+initPos)) {
      imi[numPlets][0]++;
      imi[numPlets][frame]++;
      strncpy(kCodon,s + beg, lenkCodon);
      kCodon[lenkCodon] = '\0';
      numkCodons = sigma->lenIndS(kCodon, lenkCodon);
      {
        int *bb = new int[numkCodons];
        sigma->indS(bb, kCodon, lenkCodon);
        for(int j=0;j<numkCodons;j++) {
        imi[bb[j]][0] += 1.00/numkCodons;
        imi[bb[j]][frame] += 1.00/numkCodons;
        }
        delete [] bb;
      }
    }
    delete [] kCodon;
  }                                          

  /**
   * Computes word counts for gram orders between minPlets and maxPlets
   * @see fillbaseTable
   */
  void CountUtils::fillnbaseTable(double ***im, int numPos, int minPlets, int maxPlets, char *seq, int initPos, Sigma *sigma) {
    for(int i=0; i <= maxPlets-minPlets;i++) {
      fillbaseTable(im[i], numPos, i + minPlets, seq, initPos, sigma);
    }
  }            

  //!Computes the position asymmetry index

  double CountUtils::PAS(char *seq, int numPos, int initPos, Sigma *sigma) {
    double pas = 0;
    double **fm = new double *[sigma->alphabetSize()+1];
    int i;
    assert(fm);
    for(i = 0 ; i <= sigma->alphabetSize(); i++) {
      fm[i] = new double [numPos + 1];
      assert(fm[i]);
      for(int j=0;j <= numPos;j++)
        fm[i][j] = 0;
    }
    fillbaseTable(fm, numPos, 1, seq, initPos, sigma);
    for(i=0;i < sigma->alphabetSize();i++) {
      double asym = 0, avg = fm[i][0] / 3.0;
      for(int j=1;j<=numPos;j++) {
        asym += pow(fm[i][j] - avg,2.0);
      }
      pas += asym;
    }
    for(i = 0; i <= sigma->alphabetSize(); i++)
      delete [] fm[i];
    delete [] fm;
    return pas;
  }                    

  double CountUtils:: AMI(char *seq, int numPos, int initPos, Sigma *sigma) {
    int size = sigma->alphabetSize();
    double **fm = new double *[size + 1];
    double **pij[2];
    double *pi = new double [size], Iin, Iout;
    int i;
    assert(pi);
    pij[0]= new double *[size];
    pij[1]= new double *[size];
    assert(pij[0] && pij[1]);
    for(i = 0 ; i <= size; i++) {
      fm[i] = new double [numPos + 1];
      assert(fm[i]);
      for(int j=0;j <= numPos;j++)
        fm[i][j] = 0;
    }
    for(i = 0 ; i < sigma->alphabetSize(); i++) {
      pij[0][i]= new double [size];
      pij[1][i]= new double [size];
      assert(pij[0][i] && pij[1][i]);
      for(int j=0 ; j < size; j++) {
        pij[0][i][j] = 0;
        pij[1][i][j] = 0;
      }
    }
    fillbaseTable(fm, numPos, 1, seq, initPos, sigma);

    for(i=0 ; i < size; i++) {
      for(int j=0 ; j < size; j++) {
        pij[0][i][j] = ((fm[i][1]/fm[size][1])*(fm[j][1]/fm[size][1])+(fm[i][2]/fm[size][2])*(fm[j][2]/fm[size][2])+(fm[i][3]/fm[size][3])*(fm[j][3]/fm[size][3]))/3.0;
        pij[1][i][j] = ((fm[i][1]/fm[size][1])*(fm[j][2]/fm[size][2])+(fm[i][2]/fm[size][2])*(fm[j][3]/fm[size][3])+(fm[i][3]/fm[size][3])*(fm[j][1]/fm[size][1]))/3.0;
      }
    }
    for(i= 0 ; i < size; i++) {
      pi[i] = 0;
      for(int j=1 ; j <= numPos ; j++) {
        pi[i] +=fm[i][j]/fm[size][j];
      }
      pi[i] /= 3.0;
    }
    Iin = MI(pij[0], pi, sigma);
    Iout = MI(pij[1], pi, sigma);

    for(i=0 ; i <= size; i++)
      delete [] fm[i];      
    for(i=0 ; i < size; i++) {
      delete [] pij[0][i];
      delete[] pij[1][i];
    }
    delete [] fm;
    delete [] pij[0];
    delete [] pij[1];
    delete [] pi;
    return (1.0/3.0)*Iin + (2.0/3.0)*Iout;
  }   

  //! Computes the mutual information     
  double CountUtils::MI(double *pij[], double pi[], Sigma *sigma) {
    double res = 0;
    for(int i = 0; i < sigma->alphabetSize(); i++) {
      for(int j = 0; j < sigma->alphabetSize(); j++) {
        res += pij[i][j]*(log(pij[i][j]/(pi[i]*pi[j]))/log(2.0));
      }
    }
    return res;
  }            

  //! Computes the positional Information

  double CountUtils::PI(double posFreq[], double genomFreq[], Sigma *sigma) {
    double res = 0;
    for(int i = 0 ; i < sigma->alphabetSize(); i++) {
      res += (log(posFreq[i]/genomFreq[i])/log(2.0));
    }
    return res;
  }            

}
