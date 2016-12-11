#include <iostream>
#include <streambuf>
#include <fstream>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "Sigma.h"
#include "Utils.h"

/****************************************************************************
* Utils.cpp - part of the lless namespace, a general purpose
*             linear semi-markov structure prediction library
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

  static std::map<std::string, bool> booleanVals;
  static std::map<std::string, TStrand> strands;
  static std::map<TStrand, std::string> strandNames;
  static std::map<std::string, TValType> fTypeVals;

  Utils _ut;

  Utils::Utils() {
    // boolean values
    booleanVals["true"] = true;
    booleanVals["false"] = false;
    // strands
    strands["STRAND_FWD"] = STRAND_FWD;
    strands["STRAND_COMP"] = STRAND_COMP;
    strands["NO_STRAND"] = NO_STRAND;
    strands["BOTH_STRANDS"] = BOTH_STRANDS;

    strandNames[STRAND_FWD] = "STRAND_FWD";
    strandNames[STRAND_COMP] = "STRAND_COMP";
    strandNames[BOTH_STRANDS] = "BOTH_STRANDS";

    fTypeVals["int"] = FT_INTEGER;
    fTypeVals["UCHAR"] = FT_UCHAR;
    fTypeVals["USHORT"] = FT_USHORT;
    fTypeVals["ULONG"] = FT_ULONG;
    fTypeVals["UCHAR#"] = FT_HASH;
    fTypeVals["double"] = FT_DOUBLE;
    fTypeVals["double*"] = FT_DBLVECTOR;
    fTypeVals["int*"] = FT_INTVECTOR;
    fTypeVals["EdgeInst"] = FT_SIGNAL;
    fTypeVals["ChangePoint"] = FT_CHANGEPT;
    fTypeVals["EdgeAlignment"] = FT_EDGE;
  }

  bool Utils::stringToBoolean(const std::string &s) {
    std::map<std::string, bool>::iterator it = booleanVals.find(s);

    if(it == booleanVals.end()) {
      assert(0);
      throw EXCEPTION(PARSE_ERROR, std::string("undefined Boolean Val ")+s);
    }
    return it->second;
  }

  std::string Utils::booleanToString(const bool v) {
    std::string s("false");
    if(v)
      s = std::string("true");
    return s;
  }

  TStrand Utils::stringToTStrand(const std::string &s) {
    std::map<std::string, TStrand>::iterator it = strands.find(s);
    
    if(it == strands.end()) {
      assert(0);
      throw EXCEPTION(PARSE_ERROR, std::string("undefined Strand: ")+s);
    }
    return it->second;
  }

  std::string Utils::tstrandToString(const TStrand strand) {
    std::map<TStrand, std::string>::iterator it = strandNames.find(strand);

    if(it == strandNames.end()) {
      assert(0);
      throw EXCEPTION(PARSE_ERROR, std::string("undefined Strand name "));
    }
    return it->second;
  }

  TValType Utils::stringToTValType(const std::string &s) {
    std::map<std::string, TValType>::iterator it = fTypeVals.find(s);

    if(it == fTypeVals.end()) {
      assert(0);
      throw EXCEPTION(PARSE_ERROR, std::string("undefined Type: ")+s);
    }
    return it->second;
  }
  
  /**
   * A ClassName like ClassName<Type1,.., TypeN> would return
   * TypeN as TValType
   */
  TValType Utils::extractValTypeFromClassName(string &clsId) {
    string::size_type index =  clsId.find('<', 0);
    assert(index != std::string::npos);
    string::size_type index2 =  clsId.rfind(',');
    if(index2 != std::string::npos) 
      index = index2;
    std::string t = clsId.substr(index + 1, clsId.length() - index - 2);
    return Utils::stringToTValType(t);
  }

  int Utils::permutations(int p, int n, vector< vector<int> > &perms) {
    if(p == 1) {
      for(int i = 0; i < n; i++) {
	vector<int> v;
	v.push_back(i);
	perms.push_back(v);
      }
      return n;
    }
    
    if(p == n) {
      vector<int> v;
      for(int i = 0; i < n; i++)
	v.push_back(i);
      perms.push_back(v);
      return 1;
    }
    
    int p1 = permutations(p - 1, n - 1, perms);
    for(int i = 0; i < perms.size(); i++)
      if(perms[i].size() == p - 1)
	perms[i].push_back(n - 1);
    
    int p2 = permutations(p, n - 1, perms); 
    
    return p1 + p2;
    
  }

  void Utils::strtoUpper(char *s) {
    int i=0;
    int len = strlen(s);
    while(i < len) {
      s[i] = toupper(s[i]);
      i++;            
    }
  }


  void Utils::removeBlanks(char *s) {
    int i =0, j = 0;
    int len = strlen(s);
    while(i < len) {
      if(isspace(s[i])||iscntrl(s[i])) {
        i++; continue;
      }
      s[j++] = s[i++];
    }
    s[j] = '\0';
  }


  void Utils::printTableFormat(double ***im, 
                               int numPos, 
                               int minPlets, 
                               int maxPlets, 
                               int bC, 
                               Sigma *sigma, 
                               ostream *cuStream) {

    char * plet = new char[maxPlets + 1];

    for(int i=0; i <= maxPlets-minPlets; i++) {
      int numPlets = (int)pow2(sigma->alphabetSize(),i+minPlets);
      *cuStream << "0";

      for(int k=0;k < numPos*(bC+1) + 1;k++)
        *cuStream << "\t" << setprecision(10) << im[i][numPlets][k];
      *cuStream  << "\n";

      for(int j=0;j < numPlets;j++) {
        sigma->contS(plet, j, i + minPlets);
        *cuStream << plet;

        for(int k=0;k < numPos*(bC+1) + 1;k++)
          *cuStream << "\t" << setprecision(10) << im[i][j][k];
        *cuStream << "\n";
      }

      *cuStream << "//" << "\n";

    }

    delete [] plet;

  }                                       

  void Utils::printArpaFormat(double ***im, 
                              int numPos,
                              int minPlets,
                              int maxPlets, 
                              int bC, 
                              Sigma *sigma,
                              char *prefixFile) {

    char *plet = new char[maxPlets + 1];
    char *outFile = new char[strlen(prefixFile) + 10];

    for(int k=0;k < numPos*(bC+1); k++) {
      sprintf(outFile, "%s%d", prefixFile, k);
      ::ofstream cuStream(outFile);
      cuStream << "\\data\\" << endl;

      for(int i=0; i <= maxPlets-minPlets; i++) {
        if((minPlets + i) % 3)
          continue;

        int numPlets = (int)pow2(sigma->alphabetSize(),i+minPlets);
        cuStream << "ngram " << (minPlets + i)/numPos << "=" << numPlets << endl;
      }

      for(int i=0; i <= maxPlets-minPlets; i++) {
        if((minPlets + i) % 3)
          continue;

        int numPlets = (int)pow2(sigma->alphabetSize(),i+minPlets);
        cuStream << "\\" << (minPlets + i)/numPos << "-grams:" << endl;

        for(int j=0;j< numPlets;j++) {
          if(im[i][j][k] > 0 && im[i][numPlets][k] > 0)
            cuStream << setprecision(10) << log(im[i][j][k]/im[i][numPlets][k])/log(10.0);
          else
            cuStream << "-99";
          sigma->contS(plet, j, i + minPlets);	

          for(unsigned int offs = 0; offs < strlen(plet); offs++) {
            if(!(offs % 3))
              cuStream << " ";

            cuStream << plet[offs];

          }

          if(i < maxPlets - minPlets)
            cuStream << " 0"; 

          cuStream  << "\n";

        }
      }
      cuStream << "\\end\\" << endl;
      cuStream.close();
    }                                       
    delete [] plet;
    delete [] outFile;
  }


  int Utils::posChar(const char *seq,char c) {
    unsigned int pos = 0;

    while(pos < strlen(seq) && seq[pos] != c) {
      pos++;
    }

    if(pos == strlen(seq))
      return -1;

    return pos;
  }

  void Utils::fastRead(char * buff, ifstream **file, unsigned int nbytes) {
    (*file)->read(buff, nbytes);
  }

  std::string & Utils::reverse(std::string & s) {
    for(unsigned int i = 0; i < s.length()/2; i++) {
      char tmp = s[i];
      s[i] = s[s.length() - 1 - i];
      s[s.length() - 1 - i] = tmp;
    }
    return s;
  }


  char * Utils::reverse(char *seq) {
    int i, len = strlen(seq);
    char *_tmp = strdup(seq);
    for(i = 0; i < len; i++) 
      _tmp[i] = seq[len - i - 1];
    return _tmp;
  }
  
  int Utils::substitute(std::string & text, 
                        const std::string & pattern, 
                        std::string & subst, 
                        bool global) {

    int counter = 0;
    string::size_type base = 0, index;

    while(index = text.find(pattern, base),
          index != std::string::npos) {
      // Do replacement 
      text.replace(index, pattern.length(), subst); 
      ++counter; 

      // If not global discontinue loop 
      if(!global) break; 
      base = index + 1;
    } 

    return counter; // number of substitutions made 
  }

  void Utils::revComp(char *dest, char *src, Sigma *sigma) {
    int len = strlen(src);
    for(int i=0;i < len;i++)
      dest[len-(i+1)] = sigma->toComplement(src[i]);
    dest[len] = '\0';
  }

  double Utils::seqscoreIndep(char * seq, int len) {
    double res = 0;
    for(int j = 0; j < len - 3; j++) { // dont want to count the stop codon as part of the seq.
      //  cerr << j << " ";
      res += 0.255;
    }
    return res;
  }                     

  int Utils::isbegSeqCodon(char *seq, const vector<std::string> & seqCodons) {
    if(strlen(seq) < 3) return 0;
    for(unsigned int i = 0; i < seqCodons.size(); i++) {
      if(!strncmp(seqCodons[i].c_str(), seq, 3))
        return i+1;
    }
    return 0;
  }                                                         


  int Utils::findNextCodon(const vector<std::string> & seqCodons, 
                           char *seq, 
                           unsigned int base) {

    assert(!(base%3));
    char *ptr = seq+base;

    if(base >= strlen(seq))
      return -1;

    for(unsigned int i = base ; i < strlen(seq) ; i += 3) {
      assert(ptr);

      if(isbegSeqCodon(ptr, seqCodons)) {
        return i;
      }

      ptr += 3;
    }
    return -1;
  }               

  double Utils::scoreRBSRegion(char * upsReg, 
                               int lenRBS, 
                               double ** weightRBS,
                               double *vectorRBS,
                               vector<unsigned int> &basesFreq, 
                               int *offset, 
                               Sigma *sigma) {

    int windowRBS = strlen(upsReg);
    double *accH = new double[windowRBS], res;
    int j, k;

    for(j=0 ; j < windowRBS - lenRBS +1; j++) {
      accH[j] = 0;

      for(k = 0; k < lenRBS; k++) {
        int baseCode = sigma->indC(upsReg[j+k]);

        if(baseCode < 0) 
          baseCode = sigma->indC(sigma->contC((int)(sigma->alphabetSize()*rand()/(RAND_MAX + 1.0))));

        accH[j] += weightRBS[baseCode][k];
      }     
      accH[j] += vectorRBS[j+k];  
    }

    // getting the maximum value
    for(j=0 ; j < (windowRBS - lenRBS + 1) ; j++)
      if(accH[j] > accH[*offset]) 
        *offset = j;

    res = accH[*offset];
    delete [] accH;

    return res;
  }
      
  void Utils::printProgress(char symbol, int current, int total) {
    static char symbols[] = "|/-\\|/-\\";
    static int indic = 0;
    static int cycle = 8;
    int myCurrent;
    cout.fill(symbol);
    myCurrent = (int)((double)current/total*TICKLINE_WIDTH);
    cout << setw(myCurrent) << symbol;
    cout.fill(' ');
    cout << setw(TICKLINE_WIDTH+1 - myCurrent) << "(" << symbols[indic%cycle] << ")" << setprecision(10) << myCurrent*100.0/60.0 << "%";
    cout.seekp(TICKLINE_WIDTH+4, ios::end);
  }

  void Utils::getNewLimits(int annotSeqLen, TStrand strand, int *beg, int *end) {
    if(*end > *beg) {
      if(strand == STRAND_COMP) {
        int e = *end;
        *end = annotSeqLen - *beg + 1;
        *beg = annotSeqLen - e + 1;
      }
    }
    else {
      if(strand == STRAND_FWD) {
        int e = *end;
        *end = *beg;
        *beg = e;
      }
      else {
        *beg = annotSeqLen - *beg + 1;
        *end = annotSeqLen - *end + 1;
      }
    }
  }

  int Utils::hammingDist(char *seq1, char*seq2, int *posi) {
    int l1 = strlen(seq1), l2 = strlen(seq2);
    if(l1 != l2) {
      assert(0);
      throw new exception();
    }
    int res = 0, i;
    for(i = 0; i < l1; i++) 
      posi[i] = -1;
    for(i = 0; i < l1; i++) {
      if(seq1[i] != seq2[i]) {
        posi[res] = i;
        res++;
      }
    }
    return res;
  }

  /**
   * We need to be able to sum probabilities that are represented as
   * costs (which are -log(probabilities)).  Naively, we would just
   * convert them into probabilities, sum them, and then convert them
   * back into costs.  This would be:
   * 
   * double sumNegLogProb (double a, double b) {
   * return -Math.log (Math.exp(-a) + Math.exp(-b));
   * }
   *
   * This is how this function was originally implemented, but it
   * fails when a or b is too negative.  The machine would have the
   * resolution to represent the final cost, but not the resolution to
   * represent the intermediate exponentiated negative costs, and we
   * would get -infinity as our answer.
   *  
   * What we want is a method for getting the sum by exponentiating a
   * number that is not too large.  We can do this with the following.
   * Starting with the equation above, then:
   *
   * sumNegProb = log (exp(a) + exp(b))
   * exp(sumNegProb) = exp(a) + exp(b)
   * exp(sumNegProb)/exp(a) = 1 + exp(b)/exp(a)
   * exp(sumNegProb-a) = 1 + exp(b-a)
   * sumNegProb-a = log (1 + exp(b-a))
   * sumNegProb = a + log (1 + exp(b-a)).
   *
   * We want to make sure that "b-a" is negative or a small positive
   * number.  We can assure this by noticing that we could have
   * equivalently derived
   * 
   * sumNegProb = b + log (1 + exp(a-b)),
   *
   * and we can simply select among the two alternative equations the
   * one that would have the smallest (or most negative) exponent.
   *
   */
  
  double Utils::sumLogProb (double a, double b) {
    if (isinf(a) == -1 && isinf(b) == -1)
      return -exp((double)DBL_MAX_EXP + 1);
    else if (a > b)
      return a + log (1 + exp(b-a));
    else
      return b + log (1 + exp(a-b));
  }

  int Utils::logScaleIt(int r, double base) {
    if(r >= 1) {
      if(base == 2) {
	int j = 0;
	while(r = r >> 1)
	  j++;
	/*
	while(1) {
	  if(r) j++;
	  else break;
	  r = r >> 1;
	  }*/

	return j;
      }
      r = (int)(log((double)r)/log(base)) + 1;      
    }
    return (int)r;
  }

  int Utils::pickValueAtRandom(int numVals, double *accProbs) {
    if(!accProbs)
      return rand() % numVals;
    double randFreq = accProbs[numVals - 1]*rand()/(RAND_MAX+1.0);
    for(int i = 0; i < numVals; i++)
      if(accProbs[i] >= randFreq) 
        return i;
    assert(0);
    return -1;
  }

  void Utils::permuteArray(vector<int> & arr) {
    srand((unsigned int)time(NULL));
    for(unsigned i = 0; i < arr.size(); i++) {
      int ind = pickValueAtRandom(arr.size());
      int tmp = arr[ind];
      arr[ind] = arr[i];
      arr[i] = tmp;
    }
  }

  
  /**
   *  ChiSquare constructor. It hardcodes the chi-square table.
   */

  ChiSquare::ChiSquare() {
    /*
     * Hardcoded values
     */
    chisqrTable[0][0] = 0.005;   chisqrTable[0][1] = 0.01;   chisqrTable[0][2] = 0.025;   chisqrTable[0][3] = 0.05;   chisqrTable[0][4] = 0.1;   chisqrTable[0][5] = 0.25;   chisqrTable[0][6] = 0.5;   chisqrTable[0][7] = 0.75;   chisqrTable[0][8] = 0.9;   chisqrTable[0][9] = 0.95;   chisqrTable[0][10] = 0.975;   chisqrTable[0][11] = 0.99;   chisqrTable[0][12] = 0.995; 
    chisqrTable[1][0] = 0.00003;   chisqrTable[1][1] = 0.00016;   chisqrTable[1][2] = 0.00098;   chisqrTable[1][3] = 0.0039;   chisqrTable[1][4] = 0.0158;   chisqrTable[1][5] = 0.102;   chisqrTable[1][6] = 0.455;   chisqrTable[1][7] = 1.32;   chisqrTable[1][8] = 2.71;   chisqrTable[1][9] = 3.84;   chisqrTable[1][10] = 5.02;   chisqrTable[1][11] = 6.63;   chisqrTable[1][12] = 7.88; 
    chisqrTable[2][0] = 0.01;   chisqrTable[2][1] = 0.0201;   chisqrTable[2][2] = 0.0506;   chisqrTable[2][3] = 0.103;   chisqrTable[2][4] = 0.211;   chisqrTable[2][5] = 0.575;   chisqrTable[2][6] = 1.39;   chisqrTable[2][7] = 2.77;   chisqrTable[2][8] = 4.61;   chisqrTable[2][9] = 5.99;   chisqrTable[2][10] = 7.38;   chisqrTable[2][11] = 9.21;   chisqrTable[2][12] = 10.6; 
    chisqrTable[3][0] = 0.0717;   chisqrTable[3][1] = 0.115;   chisqrTable[3][2] = 0.216;   chisqrTable[3][3] = 0.352;   chisqrTable[3][4] = 0.584;   chisqrTable[3][5] = 1.21;   chisqrTable[3][6] = 2.37;   chisqrTable[3][7] = 4.11;   chisqrTable[3][8] = 6.25;   chisqrTable[3][9] = 7.81;   chisqrTable[3][10] = 9.35;   chisqrTable[3][11] = 11.3;   chisqrTable[3][12] = 12.8; 
    chisqrTable[4][0] = 0.207;   chisqrTable[4][1] = 0.297;   chisqrTable[4][2] = 0.484;   chisqrTable[4][3] = 0.711;   chisqrTable[4][4] = 1.06;   chisqrTable[4][5] = 1.92;   chisqrTable[4][6] = 3.36;   chisqrTable[4][7] = 5.39;   chisqrTable[4][8] = 7.78;   chisqrTable[4][9] = 9.49;   chisqrTable[4][10] = 11.1;   chisqrTable[4][11] = 13.3;   chisqrTable[4][12] = 14.9; 
    chisqrTable[5][0] = 0.412;   chisqrTable[5][1] = 0.554;   chisqrTable[5][2] = 0.831;   chisqrTable[5][3] = 1.15;   chisqrTable[5][4] = 1.61;   chisqrTable[5][5] = 2.67;   chisqrTable[5][6] = 4.35;   chisqrTable[5][7] = 6.63;   chisqrTable[5][8] = 9.24;   chisqrTable[5][9] = 11.1;   chisqrTable[5][10] = 12.8;   chisqrTable[5][11] = 15.1;   chisqrTable[5][12] = 16.7; 
    chisqrTable[6][0] = 0.676;   chisqrTable[6][1] = 0.872;   chisqrTable[6][2] = 1.24;   chisqrTable[6][3] = 1.64;   chisqrTable[6][4] = 2.2;   chisqrTable[6][5] = 3.45;   chisqrTable[6][6] = 5.35;   chisqrTable[6][7] = 7.84;   chisqrTable[6][8] = 10.6;   chisqrTable[6][9] = 12.6;   chisqrTable[6][10] = 14.4;   chisqrTable[6][11] = 16.8;   chisqrTable[6][12] = 18.5; 
    chisqrTable[7][0] = 0.989;   chisqrTable[7][1] = 1.24;   chisqrTable[7][2] = 1.69;   chisqrTable[7][3] = 2.17;   chisqrTable[7][4] = 2.83;   chisqrTable[7][5] = 4.25;   chisqrTable[7][6] = 6.35;   chisqrTable[7][7] = 9.04;   chisqrTable[7][8] = 12;   chisqrTable[7][9] = 14.1;   chisqrTable[7][10] = 16;   chisqrTable[7][11] = 18.5;   chisqrTable[7][12] = 20.3; 
    chisqrTable[8][0] = 1.34;   chisqrTable[8][1] = 1.65;   chisqrTable[8][2] = 2.18;   chisqrTable[8][3] = 2.73;   chisqrTable[8][4] = 3.49;   chisqrTable[8][5] = 5.07;   chisqrTable[8][6] = 7.34;   chisqrTable[8][7] = 10.2;   chisqrTable[8][8] = 13.4;   chisqrTable[8][9] = 15.5;   chisqrTable[8][10] = 17.5;   chisqrTable[8][11] = 20.1;   chisqrTable[8][12] = 22; 
    chisqrTable[9][0] = 1.73;   chisqrTable[9][1] = 2.09;   chisqrTable[9][2] = 2.7;   chisqrTable[9][3] = 3.33;   chisqrTable[9][4] = 4.17;   chisqrTable[9][5] = 5.9;   chisqrTable[9][6] = 8.34;   chisqrTable[9][7] = 11.4;   chisqrTable[9][8] = 14.7;   chisqrTable[9][9] = 16.9;   chisqrTable[9][10] = 19;   chisqrTable[9][11] = 21.7;   chisqrTable[9][12] = 23.6; 
    chisqrTable[10][0] = 2.16;   chisqrTable[10][1] = 2.56;   chisqrTable[10][2] = 3.25;   chisqrTable[10][3] = 3.94;   chisqrTable[10][4] = 4.87;   chisqrTable[10][5] = 6.74;   chisqrTable[10][6] = 9.34;   chisqrTable[10][7] = 12.5;   chisqrTable[10][8] = 16;   chisqrTable[10][9] = 18.3;   chisqrTable[10][10] = 20.5;   chisqrTable[10][11] = 23.2;   chisqrTable[10][12] = 25.2; 
    chisqrTable[11][0] = 2.6;   chisqrTable[11][1] = 3.05;   chisqrTable[11][2] = 3.82;   chisqrTable[11][3] = 4.57;   chisqrTable[11][4] = 5.58;   chisqrTable[11][5] = 7.58;   chisqrTable[11][6] = 10.3;   chisqrTable[11][7] = 13.7;   chisqrTable[11][8] = 17.3;   chisqrTable[11][9] = 19.7;   chisqrTable[11][10] = 21.9;   chisqrTable[11][11] = 24.7;   chisqrTable[11][12] = 26.8; 
    chisqrTable[12][0] = 3.07;   chisqrTable[12][1] = 3.57;   chisqrTable[12][2] = 4.4;   chisqrTable[12][3] = 5.23;   chisqrTable[12][4] = 6.3;   chisqrTable[12][5] = 8.44;   chisqrTable[12][6] = 11.3;   chisqrTable[12][7] = 14.8;   chisqrTable[12][8] = 18.5;   chisqrTable[12][9] = 21;   chisqrTable[12][10] = 23.3;   chisqrTable[12][11] = 26.2;   chisqrTable[12][12] = 28.3; 
    chisqrTable[13][0] = 3.57;   chisqrTable[13][1] = 4.11;   chisqrTable[13][2] = 5.01;   chisqrTable[13][3] = 5.89;   chisqrTable[13][4] = 7.04;   chisqrTable[13][5] = 9.3;   chisqrTable[13][6] = 12.3;   chisqrTable[13][7] = 16;   chisqrTable[13][8] = 19.8;   chisqrTable[13][9] = 22.4;   chisqrTable[13][10] = 24.7;   chisqrTable[13][11] = 27.7;   chisqrTable[13][12] = 29.8; 
    chisqrTable[14][0] = 4.07;   chisqrTable[14][1] = 4.66;   chisqrTable[14][2] = 5.63;   chisqrTable[14][3] = 6.57;   chisqrTable[14][4] = 7.79;   chisqrTable[14][5] = 10.2;   chisqrTable[14][6] = 13.3;   chisqrTable[14][7] = 17.1;   chisqrTable[14][8] = 21.1;   chisqrTable[14][9] = 23.7;   chisqrTable[14][10] = 26.1;   chisqrTable[14][11] = 29.1;   chisqrTable[14][12] = 31.3; 
    chisqrTable[15][0] = 4.6;   chisqrTable[15][1] = 5.23;   chisqrTable[15][2] = 6.26;   chisqrTable[15][3] = 7.26;   chisqrTable[15][4] = 8.55;   chisqrTable[15][5] = 11;   chisqrTable[15][6] = 14.3;   chisqrTable[15][7] = 18.2;   chisqrTable[15][8] = 22.3;   chisqrTable[15][9] = 25;   chisqrTable[15][10] = 27.5;   chisqrTable[15][11] = 30.6;   chisqrTable[15][12] = 32.8; 
    chisqrTable[16][0] = 5.14;   chisqrTable[16][1] = 5.81;   chisqrTable[16][2] = 6.91;   chisqrTable[16][3] = 7.96;   chisqrTable[16][4] = 9.31;   chisqrTable[16][5] = 11.9;   chisqrTable[16][6] = 15.3;   chisqrTable[16][7] = 19.4;   chisqrTable[16][8] = 23.5;   chisqrTable[16][9] = 26.3;   chisqrTable[16][10] = 28.8;   chisqrTable[16][11] = 32;   chisqrTable[16][12] = 34.3; 
    chisqrTable[17][0] = 5.7;   chisqrTable[17][1] = 6.41;   chisqrTable[17][2] = 7.56;   chisqrTable[17][3] = 8.67;   chisqrTable[17][4] = 10.1;   chisqrTable[17][5] = 12.8;   chisqrTable[17][6] = 16.3;   chisqrTable[17][7] = 20.5;   chisqrTable[17][8] = 24.8;   chisqrTable[17][9] = 27.6;   chisqrTable[17][10] = 30.2;   chisqrTable[17][11] = 33.4;   chisqrTable[17][12] = 35.7; 
    chisqrTable[18][0] = 6.26;   chisqrTable[18][1] = 7.01;   chisqrTable[18][2] = 8.23;   chisqrTable[18][3] = 9.39;   chisqrTable[18][4] = 10.9;   chisqrTable[18][5] = 13.7;   chisqrTable[18][6] = 17.3;   chisqrTable[18][7] = 21.6;   chisqrTable[18][8] = 26;   chisqrTable[18][9] = 28.9;   chisqrTable[18][10] = 31.5;   chisqrTable[18][11] = 34.8;   chisqrTable[18][12] = 37.2; 
    chisqrTable[19][0] = 6.84;   chisqrTable[19][1] = 7.63;   chisqrTable[19][2] = 8.91;   chisqrTable[19][3] = 10.1;   chisqrTable[19][4] = 11.7;   chisqrTable[19][5] = 14.6;   chisqrTable[19][6] = 18.3;   chisqrTable[19][7] = 22.7;   chisqrTable[19][8] = 27.2;   chisqrTable[19][9] = 30.1;   chisqrTable[19][10] = 32.9;   chisqrTable[19][11] = 36.2;   chisqrTable[19][12] = 38.6; 
    chisqrTable[20][0] = 7.43;   chisqrTable[20][1] = 8.26;   chisqrTable[20][2] = 9.59;   chisqrTable[20][3] = 10.9;   chisqrTable[20][4] = 12.4;   chisqrTable[20][5] = 15.5;   chisqrTable[20][6] = 19.3;   chisqrTable[20][7] = 23.8;   chisqrTable[20][8] = 28.4;   chisqrTable[20][9] = 31.4;   chisqrTable[20][10] = 34.2;   chisqrTable[20][11] = 37.6;   chisqrTable[20][12] = 40; 
    chisqrTable[21][0] = 8.03;   chisqrTable[21][1] = 8.9;   chisqrTable[21][2] = 10.3;   chisqrTable[21][3] = 11.6;   chisqrTable[21][4] = 13.2;   chisqrTable[21][5] = 16.3;   chisqrTable[21][6] = 20.3;   chisqrTable[21][7] = 24.9;   chisqrTable[21][8] = 29.6;   chisqrTable[21][9] = 32.7;   chisqrTable[21][10] = 35.5;   chisqrTable[21][11] = 38.9;   chisqrTable[21][12] = 41.4; 
    chisqrTable[22][0] = 8.64;   chisqrTable[22][1] = 9.54;   chisqrTable[22][2] = 11;   chisqrTable[22][3] = 12.3;   chisqrTable[22][4] = 14;   chisqrTable[22][5] = 17.2;   chisqrTable[22][6] = 21.3;   chisqrTable[22][7] = 26;   chisqrTable[22][8] = 30.8;   chisqrTable[22][9] = 33.9;   chisqrTable[22][10] = 36.8;   chisqrTable[22][11] = 40.3;   chisqrTable[22][12] = 42.8; 
    chisqrTable[23][0] = 9.26;   chisqrTable[23][1] = 10.2;   chisqrTable[23][2] = 11.7;   chisqrTable[23][3] = 13.1;   chisqrTable[23][4] = 14.8;   chisqrTable[23][5] = 18.1;   chisqrTable[23][6] = 22.3;   chisqrTable[23][7] = 27.1;   chisqrTable[23][8] = 32;   chisqrTable[23][9] = 35.2;   chisqrTable[23][10] = 38.1;   chisqrTable[23][11] = 41.6;   chisqrTable[23][12] = 44.2; 
    chisqrTable[24][0] = 9.89;   chisqrTable[24][1] = 10.9;   chisqrTable[24][2] = 12.4;   chisqrTable[24][3] = 13.8;   chisqrTable[24][4] = 15.7;   chisqrTable[24][5] = 19;   chisqrTable[24][6] = 23.3;   chisqrTable[24][7] = 28.2;   chisqrTable[24][8] = 33.2;   chisqrTable[24][9] = 36.4;   chisqrTable[24][10] = 39.4;   chisqrTable[24][11] = 43;   chisqrTable[24][12] = 45.6; 
    chisqrTable[25][0] = 10.5;   chisqrTable[25][1] = 11.5;   chisqrTable[25][2] = 13.1;   chisqrTable[25][3] = 14.6;   chisqrTable[25][4] = 16.5;   chisqrTable[25][5] = 19.9;   chisqrTable[25][6] = 24.3;   chisqrTable[25][7] = 29.3;   chisqrTable[25][8] = 34.4;   chisqrTable[25][9] = 37.7;   chisqrTable[25][10] = 40.6;   chisqrTable[25][11] = 44.3;   chisqrTable[25][12] = 46.9; 
    chisqrTable[26][0] = 11.2;   chisqrTable[26][1] = 12.2;   chisqrTable[26][2] = 13.8;   chisqrTable[26][3] = 15.4;   chisqrTable[26][4] = 17.3;   chisqrTable[26][5] = 20.8;   chisqrTable[26][6] = 25.3;   chisqrTable[26][7] = 30.4;   chisqrTable[26][8] = 35.6;   chisqrTable[26][9] = 38.9;   chisqrTable[26][10] = 41.9;   chisqrTable[26][11] = 45.6;   chisqrTable[26][12] = 48.3; 
    chisqrTable[27][0] = 11.8;   chisqrTable[27][1] = 12.9;   chisqrTable[27][2] = 14.6;   chisqrTable[27][3] = 16.2;   chisqrTable[27][4] = 18.1;   chisqrTable[27][5] = 21.7;   chisqrTable[27][6] = 26.3;   chisqrTable[27][7] = 31.5;   chisqrTable[27][8] = 36.7;   chisqrTable[27][9] = 40.1;   chisqrTable[27][10] = 43.2;   chisqrTable[27][11] = 47;   chisqrTable[27][12] = 49.6; 
    chisqrTable[28][0] = 12.5;   chisqrTable[28][1] = 13.6;   chisqrTable[28][2] = 15.3;   chisqrTable[28][3] = 16.9;   chisqrTable[28][4] = 18.9;   chisqrTable[28][5] = 22.7;   chisqrTable[28][6] = 27.3;   chisqrTable[28][7] = 32.6;   chisqrTable[28][8] = 37.9;   chisqrTable[28][9] = 41.3;   chisqrTable[28][10] = 44.5;   chisqrTable[28][11] = 48.3;   chisqrTable[28][12] = 51; 
    chisqrTable[29][0] = 13.1;   chisqrTable[29][1] = 14.3;   chisqrTable[29][2] = 16;   chisqrTable[29][3] = 17.7;   chisqrTable[29][4] = 19.8;   chisqrTable[29][5] = 23.6;   chisqrTable[29][6] = 28.3;   chisqrTable[29][7] = 33.7;   chisqrTable[29][8] = 39.1;   chisqrTable[29][9] = 42.6;   chisqrTable[29][10] = 45.7;   chisqrTable[29][11] = 49.6;   chisqrTable[29][12] = 52.3; 
    chisqrTable[30][0] = 13.8;   chisqrTable[30][1] = 15;   chisqrTable[30][2] = 16.8;   chisqrTable[30][3] = 18.5;   chisqrTable[30][4] = 20.6;   chisqrTable[30][5] = 24.5;   chisqrTable[30][6] = 29.3;   chisqrTable[30][7] = 34.8;   chisqrTable[30][8] = 40.3;   chisqrTable[30][9] = 43.8;   chisqrTable[30][10] = 47;   chisqrTable[30][11] = 50.9;   chisqrTable[30][12] = 53.7; 
    chisqrTable[31][0] = 14.5;   chisqrTable[31][1] = 15.7;   chisqrTable[31][2] = 17.5;   chisqrTable[31][3] = 19.3;   chisqrTable[31][4] = 21.4;   chisqrTable[31][5] = 25.4;   chisqrTable[31][6] = 30.3;   chisqrTable[31][7] = 35.9;   chisqrTable[31][8] = 41.4;   chisqrTable[31][9] = 45;   chisqrTable[31][10] = 48.2;   chisqrTable[31][11] = 52.2;   chisqrTable[31][12] = 55; 
    chisqrTable[32][0] = 15.1;   chisqrTable[32][1] = 16.4;   chisqrTable[32][2] = 18.3;   chisqrTable[32][3] = 20.1;   chisqrTable[32][4] = 22.3;   chisqrTable[32][5] = 26.3;   chisqrTable[32][6] = 31.3;   chisqrTable[32][7] = 37;   chisqrTable[32][8] = 42.6;   chisqrTable[32][9] = 46.2;   chisqrTable[32][10] = 49.5;   chisqrTable[32][11] = 53.5;   chisqrTable[32][12] = 56.3; 
    chisqrTable[33][0] = 15.8;   chisqrTable[33][1] = 17.1;   chisqrTable[33][2] = 19;   chisqrTable[33][3] = 20.9;   chisqrTable[33][4] = 23.1;   chisqrTable[33][5] = 27.2;   chisqrTable[33][6] = 32.3;   chisqrTable[33][7] = 38.1;   chisqrTable[33][8] = 43.7;   chisqrTable[33][9] = 47.4;   chisqrTable[33][10] = 50.7;   chisqrTable[33][11] = 54.8;   chisqrTable[33][12] = 57.6; 
    chisqrTable[34][0] = 16.5;   chisqrTable[34][1] = 17.8;   chisqrTable[34][2] = 19.8;   chisqrTable[34][3] = 21.7;   chisqrTable[34][4] = 24;   chisqrTable[34][5] = 28.1;   chisqrTable[34][6] = 33.3;   chisqrTable[34][7] = 39.1;   chisqrTable[34][8] = 44.9;   chisqrTable[34][9] = 48.6;   chisqrTable[34][10] = 52;   chisqrTable[34][11] = 56.1;   chisqrTable[34][12] = 59; 
    chisqrTable[35][0] = 17.2;   chisqrTable[35][1] = 18.5;   chisqrTable[35][2] = 20.6;   chisqrTable[35][3] = 22.5;   chisqrTable[35][4] = 24.8;   chisqrTable[35][5] = 29.1;   chisqrTable[35][6] = 34.3;   chisqrTable[35][7] = 40.2;   chisqrTable[35][8] = 46.1;   chisqrTable[35][9] = 49.8;   chisqrTable[35][10] = 53.2;   chisqrTable[35][11] = 57.3;   chisqrTable[35][12] = 60.3; 
    chisqrTable[36][0] = 17.9;   chisqrTable[36][1] = 19.2;   chisqrTable[36][2] = 21.3;   chisqrTable[36][3] = 23.3;   chisqrTable[36][4] = 25.6;   chisqrTable[36][5] = 30;   chisqrTable[36][6] = 35.3;   chisqrTable[36][7] = 41.3;   chisqrTable[36][8] = 47.2;   chisqrTable[36][9] = 51;   chisqrTable[36][10] = 54.4;   chisqrTable[36][11] = 58.6;   chisqrTable[36][12] = 61.6; 
    chisqrTable[37][0] = 18.6;   chisqrTable[37][1] = 20;   chisqrTable[37][2] = 22.1;   chisqrTable[37][3] = 24.1;   chisqrTable[37][4] = 26.5;   chisqrTable[37][5] = 30.9;   chisqrTable[37][6] = 36.3;   chisqrTable[37][7] = 42.4;   chisqrTable[37][8] = 48.4;   chisqrTable[37][9] = 52.2;   chisqrTable[37][10] = 55.7;   chisqrTable[37][11] = 59.9;   chisqrTable[37][12] = 62.9; 
    chisqrTable[38][0] = 19.3;   chisqrTable[38][1] = 20.7;   chisqrTable[38][2] = 22.9;   chisqrTable[38][3] = 24.9;   chisqrTable[38][4] = 27.3;   chisqrTable[38][5] = 31.8;   chisqrTable[38][6] = 37.3;   chisqrTable[38][7] = 43.5;   chisqrTable[38][8] = 49.5;   chisqrTable[38][9] = 53.4;   chisqrTable[38][10] = 56.9;   chisqrTable[38][11] = 61.2;   chisqrTable[38][12] = 64.2; 
    chisqrTable[39][0] = 20;   chisqrTable[39][1] = 21.4;   chisqrTable[39][2] = 23.7;   chisqrTable[39][3] = 25.7;   chisqrTable[39][4] = 28.2;   chisqrTable[39][5] = 32.7;   chisqrTable[39][6] = 38.3;   chisqrTable[39][7] = 44.5;   chisqrTable[39][8] = 50.7;   chisqrTable[39][9] = 54.6;   chisqrTable[39][10] = 58.1;   chisqrTable[39][11] = 62.4;   chisqrTable[39][12] = 65.5; 
    chisqrTable[40][0] = 20.7;   chisqrTable[40][1] = 22.2;   chisqrTable[40][2] = 24.4;   chisqrTable[40][3] = 26.5;   chisqrTable[40][4] = 29.1;   chisqrTable[40][5] = 33.7;   chisqrTable[40][6] = 39.3;   chisqrTable[40][7] = 45.6;   chisqrTable[40][8] = 51.8;   chisqrTable[40][9] = 55.8;   chisqrTable[40][10] = 59.3;   chisqrTable[40][11] = 63.7;   chisqrTable[40][12] = 66.8; 
    chisqrTable[41][0] = 21.4;   chisqrTable[41][1] = 22.9;   chisqrTable[41][2] = 25.2;   chisqrTable[41][3] = 27.3;   chisqrTable[41][4] = 29.9;   chisqrTable[41][5] = 34.6;   chisqrTable[41][6] = 40.3;   chisqrTable[41][7] = 46.7;   chisqrTable[41][8] = 52.9;   chisqrTable[41][9] = 56.9;   chisqrTable[41][10] = 60.6;   chisqrTable[41][11] = 65;   chisqrTable[41][12] = 68.1; 
    chisqrTable[42][0] = 22.1;   chisqrTable[42][1] = 23.7;   chisqrTable[42][2] = 26;   chisqrTable[42][3] = 28.1;   chisqrTable[42][4] = 30.8;   chisqrTable[42][5] = 35.5;   chisqrTable[42][6] = 41.3;   chisqrTable[42][7] = 47.8;   chisqrTable[42][8] = 54.1;   chisqrTable[42][9] = 58.1;   chisqrTable[42][10] = 61.8;   chisqrTable[42][11] = 66.2;   chisqrTable[42][12] = 69.3; 
    chisqrTable[43][0] = 22.9;   chisqrTable[43][1] = 24.4;   chisqrTable[43][2] = 26.8;   chisqrTable[43][3] = 29;   chisqrTable[43][4] = 31.6;   chisqrTable[43][5] = 36.4;   chisqrTable[43][6] = 42.3;   chisqrTable[43][7] = 48.8;   chisqrTable[43][8] = 55.2;   chisqrTable[43][9] = 59.3;   chisqrTable[43][10] = 63;   chisqrTable[43][11] = 67.5;   chisqrTable[43][12] = 70.6; 
    chisqrTable[44][0] = 23.6;   chisqrTable[44][1] = 25.1;   chisqrTable[44][2] = 27.6;   chisqrTable[44][3] = 29.8;   chisqrTable[44][4] = 32.5;   chisqrTable[44][5] = 37.4;   chisqrTable[44][6] = 43.3;   chisqrTable[44][7] = 49.9;   chisqrTable[44][8] = 56.4;   chisqrTable[44][9] = 60.5;   chisqrTable[44][10] = 64.2;   chisqrTable[44][11] = 68.7;   chisqrTable[44][12] = 71.9; 
    chisqrTable[45][0] = 24.3;   chisqrTable[45][1] = 25.9;   chisqrTable[45][2] = 28.4;   chisqrTable[45][3] = 30.6;   chisqrTable[45][4] = 33.4;   chisqrTable[45][5] = 38.3;   chisqrTable[45][6] = 44.3;   chisqrTable[45][7] = 51;   chisqrTable[45][8] = 57.5;   chisqrTable[45][9] = 61.7;   chisqrTable[45][10] = 65.4;   chisqrTable[45][11] = 70;   chisqrTable[45][12] = 73.2; 
  }

  /**
   * @return the chi-square value for table ct[rows, cols]
   */
  double ChiSquare::chisqrValue(double ** ct, int rows, int cols) {
    double **aux = new double * [rows+1], res = 0;
    double *totalc = new double[cols], *totalr = new double[rows], total = 0;
    int i,j;
    // calculating totals
    aux[0] = new double[cols];
    for(i=0;i<rows;i++) {
      totalr[i] = 0;
        aux[i+1] = new double[cols];
    }  
    for(j=0;j<cols;j++) {
      totalc[j] = 0;
      for(i=0;i<rows;i++) {
        totalr[i] += ct[i][j];
        totalc[j] += ct[i][j];
        total += ct[i][j];
      }
    }
    for(j=0;j<cols;j++) {
      double factor, diff;
      factor = totalc[j]/total;
      aux[rows][j] = 0;
      for(i=0;i<rows;i++) {
        diff = ct[i][j] - totalr[i]*factor;
        res += pow(diff,2)/(totalr[i]*factor);
        //      cout << j << " " << i << " " << pow(diff,2) << " "<<totalr[i]*factor << "\n";
        aux[rows][j] += diff;
      }
    }
    for(j=0;j < cols;j++) {
      //    cout << aux[rows][j] << "\t";
      assert((int)aux[rows][j] == 0);
    }
    for(i = 0; i <= rows;i++) 
        delete [] aux[i];
    delete [] totalc;
    delete [] totalr;
    delete [] aux;
    return res;
  }

  /**
   * @return the chi-square confidence value for value chisqrValue
   * and df degrees of freedom.
   */
  double ChiSquare::chisqrProbDst(double chisqrValue, int df) {
    int found = 0;
    // find the extreme values;
    assert(df > 0 && df < MAX_DFS);

    if(chisqrValue < chisqrTable[df][0])
      return 0;
    if(chisqrValue > chisqrTable[df][NUM_CHITABLE_COLS - 1])
      return 1;
    for(int i=1;i < NUM_CHITABLE_COLS;i++) {
      if(chisqrValue == chisqrTable[df][i]) {
        return chisqrTable[0][i];
      }
      else if(chisqrValue < chisqrTable[df][i]) {
        found = i;
        break;
      }
    }
    if(found) // found the interval: interpolate in [found-1], [found]
      return chisqrTable[0][found-1] + (chisqrTable[0][found] - chisqrTable[0][found-1])*(chisqrValue - chisqrTable[df][found-1])/(chisqrTable[df][found] - chisqrTable[df][found-1]);
    else
      throw EXCEPTION(NOT_ANNOTATED, "An error occurred while retrieving chi square values");
    return 1;                                                   
  }

}
