/****************************************************************************
 * Utils.h - part of the lless namespace, a general purpose
 *           linear semi-markov structure prediction library
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

#ifndef _UTILS_H__
#define _UTILS_H__

#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include <cfloat>
#include <map>
#include <streambuf>
#include "GenLibExcept.h"
#include <vector>
#include "limits.h"

/****************************************************************************
 * Utils
 *
 * The Utils class provides some needed routines and definitions that are used
 * through the whole lless library
 ****************************************************************************/


#define UINT unsigned int
#define UCHAR unsigned char
#define USHORT short unsigned int
#define ULONG unsigned long

#define MAX_LINE 200000
#define NAMEFILE_LEN 200
#define TICKLINE_WIDTH 10
#define MAX_DFS 46
#define NUM_CHITABLE_COLS 13
#define DOUBLE_INFINITY exp((double)DBL_MAX_EXP - 1)
#define PRINT_TICKS(symbol, counter, total) if(!((counter)%((total + TICKLINE_WIDTH)/TICKLINE_WIDTH))) cerr << symbol;
#define NUM_STRANDS 2
#define BOTH_STRANDS NO_STRAND

using namespace std;

namespace lless {

    //! an enumerate for types
    typedef enum {
        FT_INTEGER,
        FT_DOUBLE,
        FT_UCHAR,
        FT_ULONG,
        FT_USHORT,
        FT_SIGNAL,
        FT_CHANGEPT,
        FT_HASH,
        FT_DBLVECTOR,
        FT_INTVECTOR,
        FT_SEGMENT,
        FT_EDGE
    } TValType;

    //! an enumerate for strands
    typedef enum {
        STRAND_FWD,
        STRAND_COMP,
        NO_STRAND
    } TStrand;

    class Sigma;

    /**
     * A Set of utility functions
     */
    class Utils {
    public:
        Utils();
        ~Utils() { }

        static bool stringToBoolean(const std::string &s);
        static std::string booleanToString(const bool v);
        static std::string tstrandToString(const TStrand strand);
        static TStrand stringToTStrand(const std::string &s);
        static TValType stringToTValType(const std::string &s);
        static TValType extractValTypeFromClassName(string &clsId);

        static int permutations(int p, int n, vector< vector<int> > &perms);

        static void strtoUpper(char *s);
        static void removeBlanks(char *s);

        static void printTableFormat(double ***im, int numPos, int minPlets, int maxPlets, int, Sigma *, ostream * cuStream = &cout);
        static void printArpaFormat(double ***im, int numPos, int minPlets, int maxPlets, int bC, Sigma *, char *prefixFile);

        static int posChar(const char *seq,char c);
        static std::string &reverse(std::string &s);
        static void revComp(char *dest, char *src, Sigma *sigma);
        static char * reverse(char *srcdest);

        /**
         * @see  int substitute(std::string &, const char *, const char *, bool)
         */
        static int substitute(std::string &, const std::string &, std::string &, bool);
        /**
         * A member function to split s in as many parts as there are commas
         * @param v the splitted version of s
         * @param s the input string with format  "part_1seppart_2,.."
         * @param sep field separator
         */
        template <class TClass>
            static void stringSplit(vector<TClass> & v, std::string & s, std::string sep) {
            string::size_type base = 0, index;
            TClass cls;

            while(index = s.find(sep, base),
                  index != std::string::npos) {

                std::istringstream clstr(s.substr(base, index - base));
                clstr >> cls;
                v.push_back(cls);
                base = index + 1;
            }

            std::istringstream clstr(s.substr(base, index - base));
            clstr >> cls;
            v.push_back(cls);

        }

        static void fastRead(char *buf, ifstream **file, unsigned int);

        static double seqscoreIndep(char * seq, int len);
        static int isbegSeqCodon(char *seq, const vector<std::string> & seqCodons);
        static int findNextCodon(const vector<string> &seqCodons, char *seq, unsigned int base);

        static double scoreRBSRegion(char * upsReg, int lenRBS, double ** weightRBS, double *vectorRBS,vector<unsigned int> & basesFreq, int *offset, Sigma *);
        static void printProgress(char symbol, int current, int total);

        static inline void getNewLimits(int, TStrand, int *, int *);

        static int hammingDist(char *seq1, char *seq2, int *);
        static double sumLogProb (double a, double b);
        static int logScaleIt(int r, double base);

        static int pickValueAtRandom(int numVals, double *accProbs = NULL);
        static void permuteArray(vector<int> & arr);


        /*
         * inline functions
         */

        template<class TClass>
            static pair<TClass, TClass> findMinMax(vector<TClass> &v,
                                                   int start = 0, int end = -1) {
            pair<TClass, TClass> min_max;
            int index;
            int n = end - start; //n: the number of elements to be sorted, assuming n>0
            if (n % 2 != 0) { // if n is odd
                min_max.first = v[start];
                min_max.second = v[start];
                index = start + 1;
            }
            else {// n is even
                if(v[start] < v[start + 1]) {
                    min_max.first = v[start];
                    min_max.second = v[start + 1];
                }
                else {
                    min_max.first = v[start + 1];
                    min_max.second = v[start];
                }
                index = start + 2;
            }

            TClass big, small;
            for(int i = index; i < n - 1; i = i + 2) {
                if(v[i] < v[i + 1]) { //one comparison
                    small = v[i];
                    big = v[i + 1];
                }
                else{
                    small = v[i + 1];
                    big = v[i];
                }
                if(min_max.first > small) { //one comparison
                    min_max.first = small;
                }
                if(min_max.second < big) { //one comparison
                    min_max.second = big;
                }
            }
            return min_max;
        }

        static inline  void removeBlanks(std::string & s) {
            removeBlanks((char *)s.c_str());
        }

        /**
         * substitutes string subs for string pattern in text.
         * @param global substitute all occurrences if true, only the first one
         * if false.
         */
        static inline int substitute(std::string & text, const char *pattern, const char *subst, bool global) {
            std::string _pattern = std::string(pattern);
            std::string _subst = std::string(subst);
            return substitute(text, _pattern, _subst, global);
        }

        /**
         * This pow function only works with bases that are multiples of 2
         */
        static inline int pow2(int b, int p) {
            if(p == 0)
                return 1;
            return b << ((p - 1)*(b/2));
        }

        static inline double min(double a, double b) {
            return a > b ? b : a;
        }

        static inline double max(double a, double b) {
            return a > b ? a : b;
        }

        static inline void revComp(char *srcdest, Sigma *sigma) {
            char *_tmp = strdup(srcdest);
            revComp(srcdest, _tmp, sigma);
            delete _tmp;
        }

        static inline void swap(char &a, char &b) {
            char tmp = b;
            b = a;
            a = tmp;
        }

        static inline void fixWithinBounds(int &beg, int &end, int length) {
            if(beg < 1) beg = 1;
            if(end > length) end = length;
        }

        /**
         * A function to smooth(average) scores, using a sliding window
         * of size smoothenWSize. smootheWSize and smoothenWSize/2 must be
         * multiple of period
         * @param filterVals the array containing scores of type TClass
         * @param len the length to consider. It could be shorter than the
         * @param period the period
         * @param smoothWindow the windowSize to smooth out scores. Must
         *        be multiple of period
         * actual length of the input sequence seq.
         */
        template <class TClass>
            static void smoothScores(vector<TClass> &vals, int period,
                                     int smoothWindow, int len) {

            TClass *fvals = &vals.front();
            smoothScores(fvals, period, smoothWindow, len);
        }

        template <class TClass>
            static void smoothScores(TClass *filterVals, int period,
                                     int smoothWindow, int len) {
            int i = 0, j;
            int wSize = smoothWindow;

            if(wSize < 0 || len < wSize)
                wSize = len;

            if(wSize < 2)
                return;

            vector<double> sc(period), scores(len + 1);

            // 5 prime
            for(i = 1; i <= wSize; i += period)
                for(j = 0; j < period; j++)
                    if(i + j <= len)
                        sc[j] += filterVals[i+j];
            for(i = 1; i <= wSize/2; i += period)
                for(j = 0; j < period; j++)
                    if(i+j <= len)
                        scores[i+j] = sc[j]*period/wSize;

            // middle
            for(; i <= len - wSize/2; i += period) {
                for(j = 0; j < period; j++) {
                    if(i+j <= len) {
                        sc[j] += filterVals[i+j+wSize/2] - filterVals[i+j-wSize/2];
                        scores[i+j] = sc[j]*period/wSize;
                    }
                }
            }

            // 3 prime
            for(;i <= len; i += period)
                for(j = 0; j < period; j++)
                    if(i+j <= len)
                        scores[i+j] = sc[j]*period/wSize;

            // updating score arrays
            for(i = 1; i <= len; i++)
                filterVals[i] = (TClass)scores[i];

        }

        template <class TClass>
            static void sequentialKMeans(vector<TClass> &vals, double lambda,
                                         vector<TClass> &means, vector<int> t,
                                         int len) {

            TClass *fvals = &vals.front();
            sequentialKMeans(fvals, lambda, means, t, len);
        }

        template <class TClass>
            static void sequentialKMeans(TClass *filterVals, double lambda,
                                         vector<TClass> &means, vector<int> t,
                                         int len) {

            for(int i = 1; i <= len; i++) {
                //	  cout << last << " " << i << " " << coverage[i] << endl;
                int closest_mean = 0;
                double min_error = DOUBLE_INFINITY;

                for(int j = 0; j < means.size(); j++) {
                    double error = fabs(means[j] - filterVals[i]);
                    if(error < min_error) {
                        min_error = error;
                        closest_mean = j;
                        t[closest_mean]++;
                    }
                }

                means[closest_mean] += (1/lambda*t[closest_mean])*(filterVals[i] - means[closest_mean]);

            }
        }
    };

    /**
     * ChiSquare is a class that implements a chi-square test.
     * The table with the distribution probabilities has been
     * hardcoded within the class.
     */

    class ChiSquare {
    private:
        double chisqrTable[MAX_DFS][NUM_CHITABLE_COLS];
        int chstblnumProb;
        int chstblDfs;
    public:
        ChiSquare();
        ~ChiSquare() {}
        double chisqrValue(double ** ct, int rows, int cols);
        double chisqrProbDst(double chisqrValue, int df);
    };

    /**
     * Comparator class for integers. Integers will appear in desceding
     * order if a list is ordered using this function as a comparator.
     */
    struct intDesCmp {
        bool operator()( const int i1, const int i2 ) const {
            return i1 > i2;
        }
    };

    /**
     * Comparator class for integers. Integers will appear in ascending
     * order if a list is ordered using this function as a comparator.
     */
    struct intCmp {
        bool operator() (const int a, const int b) const {
            return a < b;
        }
    };

    /**
     * Pair of fields of class T. Type T's operator= must be implemented.
     */
    template <class T> class Pair {
    public:
        T f; //!< first component
        T s; //!< second component
        Pair() {
        }
        Pair(T _f, T _s) {
            updVals(_f, _s);
        }
        Pair(const Pair<T> & pi) {
            f = ((Pair<T> &)pi).f;
            s = ((Pair<T> &)pi).s;
        }
        inline Pair<T> & updVals(T _f, T _s) {
            f = _f;
            s = _s;
            return *this;
        }
        inline bool operator==(const Pair<T> & pi) {
            return (this->f == ((Pair<T> &)pi).f && this->s == ((Pair<T> &)pi).s);
        }
        ~Pair<T>() {}
    };

    class DisjSetsFast : public vector<int> {
    public:
        DisjSetsFast(int sz) {
            this->resize(sz, -1);
        }

        DisjSetsFast() {
        }

        void reset(int sz = -1) {
            sz = (sz >= 0) ? sz : this->size();
            this->resize(sz, -1);

        }

        void dsUnion(int root1, int root2 ) {
            vector<int> &v = *this;
            if(v[root2] < v[root1] )  // root2 is deeper
                v[root1] = root2;        // Make root2 new root
            else {
                if(v[root1] == v[root2])
                    v[root1]--;          // Update height if same
                v[root2] = root1;        // Make root1 new root
            }
        }

        int dsFind(int x) {
            vector<int> &v = *this;
            if(v[x] < 0 )
                return x;
            else
                return v[x] = dsFind(v[x]);
        }
    };

}

#endif
