/****************************************************************************
 * RSqUtils.h - part of the craig namespace, a genomics library
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

#ifndef _RNASEQ_UTILS_H_
#define _RNASEQ_UTILS_H_
#include "FeatureVector.h"
#include "FL_EvdEdgeAligner.h"
#include <numeric>

namespace craig {

    typedef enum {TR_UPSTREAM, TR_DOWNSTREAM, TR_ANYWHERE} TTranscriptRegion;
    /**
     * The RNASeqUtils class implements routines for converting Tag to Gene
     * objects and visceversa. It also contains utility functions for
     * Gene, Transcripts and Exon objects
     ***************************************************************************/

    struct RSqChangePoint {
        int pos;
        float conf;
        float left_mean;
        float right_mean;
        float score;
        RSqChangePoint(int pos = 0, float conf = -DOUBLE_INFINITY,
                       float left_mean = 0, float right_mean = 0) {
            this->pos = pos;
            this->conf = conf;
            this->left_mean = left_mean;
            this->right_mean = right_mean;
            this->score = 0;
        }

        void setMeans(double lmean, double rmean) {
            this->left_mean = lmean;
            this->right_mean = rmean;
        }

        ~RSqChangePoint() {}
    };

    class RSqInterval {

    public:
        int start;
        int end;
        TStrand strand;
        float weight;
        double cov_mean;
        vector<double> cov_quartiles;
        double cov_stdev;
        RSqInterval(int _start = 0, int _end = 0, TStrand _strand = STRAND_FWD) {
            start = _start;
            end = _end;
            weight = -1;
            strand = _strand;
        }
        RSqInterval(int _start, int _end, double _weight, TStrand _strand) {
            start = _start;
            end = _end;
            weight = _weight;
            strand = _strand;
        }

        friend ostream& operator<< (ostream &out, RSqInterval &);
        void computeCoverageStats(vector<float>::iterator first,
                                  vector<float>::iterator last);
        inline void setWeight(float weight) {
            this->weight = weight;
        }
        bool operator< (const RSqInterval& re);
    };

    class RSqUtils {

    public:
        RSqUtils() {

        }

        static bool validateIntron(int lbegin, int lend, int ibegin, int iend,
                                   int rbegin, int rend, float junction_support,
                                   bool missing, vector<float> & cov);

        static bool interval_match(list<RSqInterval> &I,
                                   int my_start, int my_end);

        static void insert_interval(list<RSqInterval> &I,
                                    int my_start, int my_end,
                                    float weight, TStrand strand);

        // For the hashing functions, I is a hash structure. Collisions
        // are handled by using the next available spot in the array

        static bool hash_interval_match(vector<pair<int, int> > &I,
                                        int my_start, int my_end);

        static bool insert_hash_interval(pair<int, int> ival,
                                         vector<pair<int, int> > &I);

        static void computeSeqCoverage(Sequence &c,
                                       TypedFilter<double> *coverage,
                                       vector<float> *seq_cov,
                                       bool stranded);

        static void computeJunctionCoverage(Sequence &c,
                                            FL_EvdEdgeAligner *junctions,
                                            vector<float> *junction_cov,
                                            bool stranded);

        static void getRSqIntrons(Sequence &c,
                                  FL_EvdEdgeAligner *junctions,
                                  list<RSqInterval> &introns,
                                  TStrand origStrand = BOTH_STRANDS);

        static void getRSqIntrons(Sequence &c,
                                  FL_EvdEdgeAligner *junctions,
                                  bool keep_bestolap,
                                  vector<RSqInterval> &introns,
                                  TStrand origStrand = BOTH_STRANDS);

        static bool computeSplSCoverage(TypedFilter<EdgeInst> **signals,
                                        vector<float> &seq_cov,
                                        int window,
                                        int begin, int end,
                                        vector<int> &splpos,
                                        vector<float> &splcov,
                                        float zero_cutoff,
                                        float zero_value,
                                        TTranscriptRegion,
                                        TStrand strand);

        static float findMean(vector<float> &X, int start, int end);
        static float findSdiff(vector<float> &X, int start, int end);
        static float findMSE(vector<float> &X, int start, int end);
        static float meanShiftConfidence(vector<float> &X, int start, int end,
                                         int num_bootstraps, int blockLen = 1);
        static int getIntervalRSqChangePoint(int start, int end,
                                             vector<float> &X,
                                             DENSE_HASH<int, int> * = NULL);

        static void computeRSqChangePoints(vector<float> &X,
                                           vector<RSqChangePoint> &cpoints,
                                           int num_bootstraps,
                                           int start, int end,
                                           int blockLen = 1,
                                           double min_conf = 80,
                                           vector<int> *positions = NULL,
                                           DENSE_HASH<int,int> * = NULL);

        static void refineRSqChangePoints(vector<float> &X,
                                          vector<RSqChangePoint> &cpoints,
                                          int num_bootstraps,
                                          int start, int end,
                                          int blockLen = 1,
                                          double min_conf = 80,
                                          vector<int> *positions = NULL,
                                          DENSE_HASH<int, int> * = NULL);

        static void displayRSqChangePoints(vector<RSqChangePoint> &cpoints,
                                           int start, int end,
                                           vector<int> *positions);

        static bool isRealOnset(RSqChangePoint *cp,
                                vector<RSqChangePoint> &changes,
                                int start, int end,
                                float ovmean_l, float ovmean_r,
                                int level);

        static bool isRealOffset(RSqChangePoint *cp,
                                 vector<RSqChangePoint> &changes,
                                 int start, int end,
                                 float ovmean_l, float ovmean_r,
                                 int level, bool = false);

        static void changePoints2Histogram(vector<RSqChangePoint> &cps,
                                           int start, int end, int pos,
                                           char = '+', int = 0);

        static void selectOnsetPoints(vector<RSqChangePoint> &changes,
                                      int start, int end,
                                      vector<float> &X,
                                      double lovmean, double rovmean,
                                      double min_conf,
                                      RSqChangePoint &onset);

        static void selectOffsetPoints(vector<RSqChangePoint> &changes,
                                       int start, int end,
                                       vector<float> &X,
                                       double lovmean, double rovmean,
                                       double min_conf,
                                       RSqChangePoint &offset);

        static void changePoints2Filter(int offset,
                                        string &symbols,
                                        vector<RSqChangePoint> &changes,
                                        vector<float> &X,
                                        double min_conf,
                                        int clen,
                                        vector<char> &sig_filter,
                                        vector<vector<int> > &sigscore_filter,
                                        TStrand strand);

        static void onsets2Filter(vector<RSqChangePoint> &changes, char symbol,
                                  int start, int end,
                                  vector<char> &filter,
                                  vector<vector<int> > &sigscore_filter,
                                  double lovmean, double rovmean,
                                  double min_conf,
                                  int clen,
                                  vector<float> &X, TStrand strand);

        static void offsets2Filter(vector<RSqChangePoint> &changes, char symbol,
                                   int start, int end,
                                   vector<char> &filter,
                                   vector<vector<int> > &sigscore_filter,
                                   double lovmean, double rovmean,
                                   double min_conf,
                                   int clen,
                                   vector<float> &X, TStrand strand);

        static void filter2XFastaFile(vector<char> &filter,
                                      ofstream & fd);
        static void filter2MultiScoreFile(vector<vector<int> > &filter,
                                          int num_cols,
                                          ofstream & fd);

        static void storeXcriptSignals(Sequence &c,
                                       vector<RSqChangePoint> & pchanges,
                                       double intron_mean, double exon_mean,
                                       double min_conf,
                                       std::ofstream &sigd,
                                       std::ofstream &sigscored,
                                       vector<float> &, TStrand);

        static void find_peaks(vector<float> &v,
                               vector<pair<int, float> > &min,
                               vector<pair<int, float> > &max,
                               float mn_delta, float mx_delta,
                               vector<int> *positions = NULL);

        static void computeRamps(vector<float> &cov,
                                 int window,
                                 float score_cutoff,
                                 vector<pair<int, float> > &onsets,
                                 vector<pair<int, float> > &offsets,
                                 vector<int> *positions = NULL);
    };
}

#endif
