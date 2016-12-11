/****************************************************************************
* StatUtils.h - part of the lless namespace, a general purpose
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

#ifndef _STAT_UTILS_H_
#define _STAT_UTILS_H_
#include "Utils.h"

namespace lless {
  class StatUtils {
   public:
    StatUtils() {
    }
    
    static double findMean(vector<double> &X, int start, int end) {
      double meanX = 0;
      for(int i = start; i < end; i++)
	meanX += X[i];
      
      int len = end - start;
      if(len) return meanX/len;
      return 0;
    }

    static double findCUSUMdiff(vector<double> &X, int start, int end) {
      float meanX = findMean(X, start, end);
      
      // find cumulative score S
      vector<double> S(1, 0);
      //    cout.precision(4);
      //    cout << "S = ";
      for(int i = start; i < end; i++) {
	S.push_back(S[i - start] + X[i] - meanX);
	//      cout << " " << S.back();
      }
      
      // find Smax and Smin
      pair<double, double> Smin_max = Utils::findMinMax(S, 0, S.size());
      //    cout << " min,max = " << Smin_max.first << "," << Smin_max.second << endl;    
      return Smin_max.second - Smin_max.first;
      
    }

    static double findMSE(vector<double> &X, int start, int end) {
      // finding mean
      float Xmean = findMean(X, start, end);
      float MSE = 0;
      for(int i = start; i < end; i++)
	MSE += pow(X[i] - Xmean, 2.0);
      
      return MSE;
    }

    ~StatUtils() {
    }    
  };
}

#endif
