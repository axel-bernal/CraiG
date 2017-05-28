/****************************************************************************
 * ChangePtUtils.h - part of the craig namespace, a genomics library
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

#ifndef _CHANGEPT_UTILS_H_
#define _CHANGEPT_UTILS_H_
#include "FeatureVector.h"
#include "StatUtils.h"
#include <numeric>

namespace lless {

    struct ChangePoint {
        int pos;
        double conf;
        double left_mean;
        double right_mean;
        double score;
        ChangePoint(int pos = 0, double conf = -DOUBLE_INFINITY,
                    double left_mean = 0, double right_mean = 0) {
            this->pos = pos;
            this->conf = conf;
            this->left_mean = left_mean;
            this->right_mean = right_mean;
            this->score = 0;
        }

        void set_means(double lmean, double rmean) {
            this->left_mean = lmean;
            this->right_mean = rmean;
        }

        ~ChangePoint() {}
    };

    class ChangePtUtils {

    public:
        ChangePtUtils() {

        }

        static void onsets2Filter(vector<ChangePoint> &changes, char symbol,
                                  char *fvals,
                                  double left_ovmean, double right_ovmean,
                                  int clen, TStrand strand);

        static void offsets2Filter(vector<ChangePoint> &changes, char symbol,
                                   char *fvals,
                                   double left_ovmean, double right_ovmean,
                                   int clen, TStrand strand);


        static double meanShiftConfidence(vector<double> &X,
                                          int start, int end,
                                          int num_bootstraps,
                                          int block_len);
        static int findChangePoint(int start, int end,
                                   vector<double> &X);

        static void computeChangePoints(vector<double> &X,
                                        vector<ChangePoint> &cpoints,
                                        int num_bootstraps,
                                        int start, int end,
                                        int block_len = 1,
                                        double min_conf = 80,
                                        vector<int> *positions = NULL);

        static void refineChangePoints(vector<double> &X,
                                       vector<ChangePoint> &cpoints,
                                       int num_bootstraps,
                                       int start, int end,
                                       int block_len = 1,
                                       double min_conf = 80,
                                       vector<int> *positions = NULL);
    };
}

#endif
