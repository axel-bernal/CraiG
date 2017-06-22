#include "ChangePtUtils.h"

/****************************************************************************
 * ChangePtUtils.cpp - part of the craig namespace, a genomics library
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

    void ChangePtUtils::onsets2Filter(vector<ChangePoint> &changes, char symbol,
                                      char *fvals,
                                      double left_ovmean, double right_ovmean,
                                      int clen, TStrand strand) {

        int i;
        for(i = 0; i < changes.size(); i++) {
            ChangePoint &cp = changes[i];
            if(cp.pos > clen)
                assert(0);

            if(cp.left_mean > cp.right_mean)
                continue;

            int pos = (strand == STRAND_FWD ? cp.pos : clen - cp.pos + 1);
            fvals[pos] = symbol;
        }
    }

    void ChangePtUtils::offsets2Filter(vector<ChangePoint> &changes, char symbol,
                                       char *fvals,
                                       double left_ovmean, double right_ovmean,
                                       int clen, TStrand strand) {


        int i;
        for(i = 0; i < changes.size(); i++) {
            ChangePoint &cp = changes[i];
            if(cp.pos > clen)
                assert(0);

            if(cp.left_mean < cp.right_mean)
                continue;

            int pos = (strand == STRAND_FWD ? cp.pos : clen - cp.pos + 1);
            fvals[pos] = symbol;
        }
    }

    double ChangePtUtils::meanShiftConfidence(vector<double> &X,
                                              int start, int end,
                                              int num_bootstraps,
                                              int block_len) {

        // check whether block_len is too big
        vector<pair<int, int> > bindexes;
        int i = 0;

        while(1) {
            if(i >= end - start)  break;
            int bSize = (i + block_len < end - start) ? block_len : end - start - i;
            bindexes.push_back(pair<int, int>(i, bSize));
            i += bSize;
        }

        // see if there is a change
        int num_changed = 0;

        vector<double> Y(end - start), OY(end - start);
        for(i = start; i < end; i++)
            Y[i - start] = OY[i - start] = X[i];

        // finding Sdiff
        double Sdiff = StatUtils::findCUSUMdiff(OY, 0, OY.size());

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
            double mySdiff = StatUtils::findCUSUMdiff(Y, 0, Y.size());
            //      cerr << Sdiff << " " << mySdiff << " " << endl;
            if(mySdiff < Sdiff)
                num_changed++;

        }
        //    cerr << "bootstrapping results " << num_changed << endl;
        return num_changed*100.0/num_bootstraps;
    }

    int ChangePtUtils::findChangePoint(int start, int end,
                                       vector<double> &X) {

        int m = start + 1;
        int min_m = -1;
        double min_MSE = DOUBLE_INFINITY;

        for( ; m < end - 1; m++) {
            double MSE = StatUtils::findMSE(X, start, m);
            MSE += StatUtils::findMSE(X, m, end);

            if(MSE < min_MSE) {
                min_MSE = MSE;
                min_m = m;
            }
        }
        return min_m;
    }

    void ChangePtUtils::computeChangePoints(vector<double> &X,
                                            vector<ChangePoint> &cpoints,
                                            int num_bootstraps,
                                            int start, int end,
                                            int block_len,
                                            double min_conf,
                                            vector<int> *positions) {

        if(positions && X.size() != positions->size())
            throw EXCEPTION(BAD_USAGE, "Input vectors X and positions must have same length");

        double conf = meanShiftConfidence(X, start, end, num_bootstraps, block_len);
        //    cerr << "conf " << conf << " for " << start << " " << end << endl;
        if(conf >= min_conf) {
            int min_m = ChangePtUtils::findChangePoint(start, end, X);

            if(start >= min_m || min_m + 1 >= end)
                return;

            //      cout << level << "change_point " << min_m << "(" << x_i << ") "  << conf << endl;
            double left_mean = StatUtils::findMean(X, start, min_m);
            double right_mean = StatUtils::findMean(X, min_m, end);

            ChangePtUtils::computeChangePoints(X, cpoints, num_bootstraps,
                                               start, min_m,
                                               block_len, min_conf, positions);

            ChangePtUtils::computeChangePoints(X, cpoints, num_bootstraps,
                                               min_m + 1, end,
                                               block_len, min_conf, positions);

            vector<ChangePoint>::iterator it = cpoints.size() ?
                cpoints.begin() : cpoints.end();

            int x_i = positions ? (*positions)[min_m] : min_m;
            while(it != cpoints.end() && it->pos < x_i)
                it++;

            cpoints.insert(it, ChangePoint(x_i, conf, left_mean, right_mean));

        }
    }


    void ChangePtUtils::refineChangePoints(vector<double> &X,
                                           vector<ChangePoint> &cpoints,
                                           int num_bootstraps,
                                           int start, int end,
                                           int block_len,
                                           double min_conf,
                                           vector<int> *positions) {

        if(positions && X.size() != positions->size())
            throw EXCEPTION(BAD_USAGE, "Input vectors X and positions must have same length");

        //    cout << "before refinement\n";
        //    displayChangePoints(cpoints, start, end, positions);
        int A = start, i;
        for(i = 0; i <  cpoints.size(); i++) {
            int B = (i  + 1 < cpoints.size() ? cpoints[i + 1].pos : end);
            ChangePoint &cp = cpoints[i];
            // updating X[cp.pos]
            cp.conf = meanShiftConfidence(X, A, B, num_bootstraps, block_len);
            cp.set_means(StatUtils::findMean(X, A, cp.pos),
                         StatUtils::findMean(X, cp.pos, B));
            double min_mean = Utils::min(cp.left_mean, cp.right_mean);
            double variation = fabs(cp.right_mean - cp.left_mean);

            //      cerr << (positions ? (*positions)[A] : A) << " " << (positions ? (*positions)[cp.pos] : cp.pos) << " " << (positions ? (*positions)[B] : B) << " " << " " << cp.left_mean << " " << cp.right_mean << " " << variation << " " << cp.conf;

            if(cp.conf < min_conf) {
                // recalculate B
                int C = (i + 2 < cpoints.size() ? cpoints[i + 2].pos : end);
                //	cerr << " out " << cp.pos << " and " << B;
                if(B < C) {
                    int min_m = ChangePtUtils::findChangePoint(A, C, X);

                    if(A < min_m && min_m + 1 < C) {
                        cpoints[i + 1].pos = min_m;
                        //	    cerr << " for " << min_m;
                    }
                }

                vector<ChangePoint>::iterator it = cpoints.begin() + i;
                it = cpoints.erase(it);
                i--;
            }
            //      cerr << endl;
            A = i < 0 ? start : cpoints[i].pos;
        }
        //    cout << "after refinement\n";
        //    displayChangePoints(cpoints, start, end, positions);

        for(i = 0; i <  cpoints.size(); i++) {
            ChangePoint &cp = cpoints[i];
            if(positions)
                cp.pos = (*positions)[cp.pos];
        }
    }
}
