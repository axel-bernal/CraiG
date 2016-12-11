/****************************************************************************
* CountUtils.h - part of the lless namespace, a general purpose
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

#ifndef _COUNT_UTILS_H_
#define _COUNT_UTILS_H_

#include "Filter.h"
#include "Sigma.h"


namespace lless {

  /**
   * CountUtils
   *
   * The CountUtils class contains suborutines to count n-grams which then
   * are using by counting features or filters.
   ***************************************************************************/

  class CountUtils {
   public:
    static void countNGrams(double ****gramCounts, Sigma *sigma, int numPos,
                            int minOrder, int maxOrder, char *seq, 
                            UCHAR *contextVals, TStrand);

    static int initkcodonMatrix(double ****im, int numPos, int minPlets,
                                int maxPlets, int, Sigma *);
    static void freekcodonMatrix(double ****im, int minPlets, int maxPlets, 
                                 Sigma *);
    static void fillbaseTable(double **imi, int numPos, int lenkCodon, 
                              char *seq, int, Sigma *);
    static void fillnbaseTable(double ***im, int numPos, int minLenFL_Gram, 
                               int maxLenFL_Gram, char *seq, int, Sigma *);

    static double PAS(char *seq, int numPos, int, Sigma *);
    static double AMI(char *seq, int numPos, int, Sigma *);
    static double MI(double *pij[], double pi[], Sigma *);
    static double PI(double posFreq[], double genomFreq[], Sigma *); 

    CountUtils();
    ~CountUtils() { }
  };

}

#endif
