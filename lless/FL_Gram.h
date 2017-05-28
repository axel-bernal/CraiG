/****************************************************************************
 * FL_Gram.h - part of the lless namespace, a general purpose
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

#ifndef _FT_GRAM_H
#define _FT_GRAM_H

#include <streambuf>
#include <string>
#include <string.h>
#include "Utils.h"
#include "Sigma.h"
#include "Filter.h"
#include <math.h>
#include <assert.h>

#define MAX_ORDER 8
#define MIN_FREQ 1e-50
#define MIN_PROB exp((double)-DBL_MAX_EXP)
#define MAX_GRAM_ORDER 9

namespace lless {

    class ResourceEngine;
    class FilterEngine;

    template <class TClass1, class TClass2>
        class FL_BaseGram : public TypedFilter<TClass2> {
    protected:
        unsigned int _order;
        int _numHistories;
        unsigned int _alphabetSize;

    public:
    FL_BaseGram(int fInd,
                std::string & name,
                TValType type,
                int order,
                int alphabetSize,
                int period = 0) :
        TypedFilter<TClass2>(fInd,
                             name,
                             (int)(pow((double)alphabetSize, order)),
                             type) {

            _order = order;
            _alphabetSize = alphabetSize;
            _numHistories = (pow((double)alphabetSize, _order - 1));
            this->_period = period;

        }

    FL_BaseGram(int fInd,
                vector<std::string> & params,
                int & offset) :
        TypedFilter<TClass2>(fInd,
                             params,
                             offset,
                             0) {

            if(!sscanf(params[offset++].c_str(), "%u", &_order))
                assert(0);

        }

        FL_BaseGram(FL_BaseGram &g);

        virtual int wmap(TClass1 *seq, int beg, int end) = 0;
        virtual int wmap(TClass1 *seq, int wlen) = 0;

        inline int alphabetSize() {
            return _alphabetSize;
        }

        /**
         * @return the total number of words of length order - 1
         */
        inline int numHistories() {
            return _numHistories;
        }

        inline unsigned int order() {
            return _order;
        }

        /**
         * A member function for computing filter values.
         *
         * The value at each position of the input sequence is the gram found in the
         * sequence at the position.
         *
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */

        virtual void computeVals(char *inputSeq, TClass2 *vals, int len) {
            TClass1 *seq = (TClass1 *)inputSeq;
            unsigned int i = 1;
            assert(order() > 0);

            for( ; i + order() - 1 <= (unsigned int)len; i++)
                vals[i] = (TClass2)wmap(seq, i - 1, i + order() - 2);

        }

        virtual ~FL_BaseGram() {  }

    };

    /**
     * The FL_Gram class computes a distribution of order-gram words over the
     * input sequence.
     * This class  derives from TypedFilter<int> and as such it sets the filter\
     * value at any position p equal to the gram that appears at p in the
     * sequence.
     */

    template <class TClass>
        class FL_Gram : public FL_BaseGram<char, TClass> {
    protected:
        Sigma *_alphabet;
        unsigned int totFreq;
        double *transMat;
        vector<unsigned int> _freqs;

    public:
    FL_Gram(int fInd,                 //!< A unique identifier.
            std::string & name,       //!< A unique name.
            Sigma *alphabet,          //!< The FL_Gram's alphabet.
            TValType type,
            int order,                //!< The FL_Gram's order.
            int period = 0
        ) :
        FL_BaseGram<char, TClass>(fInd,
                                  name,
                                  type,
                                  order,
                                  alphabet->alphabetSize(),
                                  period) {

            _alphabet = alphabet;
            initGram();
        }


    FL_Gram(int fInd,           //!< A unique identifier.
            /*! The Header string definition, loaded as a vector of strings.
             * The Header has the following form:\n\n
             * Filter name FL_Gram arrInd sigma order compactHomogBlocks\n\n
             * The description of the fields could be found in
             * the other constructor(s)
             */
            vector<std::string> & params,
            int & offset,       //!< The index for vector params.
            ResourceEngine *re, //!< A pointer to the ResourceEngine object.
            FilterEngine *fe    //!< A pointer to the FilterEngine object.
        ) :
        FL_BaseGram<char, TClass>(fInd, params, offset) {

            _alphabet = (Sigma *)re->getResource(params[offset++]);
            this->_period = 0;

            if(offset < params.size())
                if(!sscanf(params[offset++].c_str(), "%u", &this->_period))
                    assert(0);

            this->_alphabetSize = _alphabet->alphabetSize();
            this->maxFilterVals = (int)(pow((double)this->alphabetSize(), this->order()));
            this->_numHistories = (int)(pow((double)this->alphabetSize(), this->order() - 1));
            initGram();

        }

        /**
         * Copy constructor
         * @param g a reference to a FL_Gram object
         */
    FL_Gram(FL_Gram<TClass> &g) :
        FL_BaseGram<char, TClass>(g.ind(),
                                  g.getName(),
                                  g.valType(),
                                  g.order(),
                                  g.alphabetSize()) {

            *this = g;
            transMat = NULL;
            buildCp();
        }

        inline int wmap(char *seq, int beg, int end) {
            int ind = alphabet()->randIndS(seq, beg, end);
            assert(ind < this->maxNumFilterValues());
            return ind;
        }

        inline int wmap(char *seq, int wlen) {
            int ind = alphabet()->indS(seq, wlen);
            assert(ind < this->maxNumFilterValues());
            return ind;
        }

        inline int wmap(int suffix, int word) {
            int ind = _alphabet->addSuffix(word, suffix);
            assert(ind < this->maxNumFilterValues());
            return ind;
        }

        void addWord(char *word, unsigned int freq = 0) {
            int index = this->wmap(word, this->order());
            _freqs[index] = freq;
        }

        inline unsigned int getWordFreq(int ind) {
            return _freqs[ind];
        }

        inline void setWordFreq(TClass ind, unsigned int value) {
            _freqs[ind] = value;
        }

        inline vector<unsigned int> & wordCounts() {
            return _freqs;
        }

        inline unsigned int getWordFreq() {
            return totFreq;
        }

        inline void setWordFreq(unsigned int value) {
            totFreq = value;
        }

        inline Sigma *alphabet() {
            return _alphabet;
        }

        /**
         * @return the relative frequency of symbol i  when it appears immediately
         * after word j
         * @param i a symbol of length 1
         * @param j a word of length order - 1
         */
        inline double getCp(int i, int j) {
            return transMat[this->wmap(i, j)];
        }

        /**
         * A member function that sets the relative frequency of symbol
         * i when it appears immediately after word j.
         * @param i a symbol of length 1
         * @param j a word of length order - 1
         * @param val the relative frequency
         */
        inline void setCp(int i, int j, double val) {
            transMat[this->wmap(i, j)] = val;
        }

        /**
         * @return true if table with relative word frequencies has been built,
         * false otherwise.
         */
        inline bool cpDone() {
            if(transMat)
                return true;
            else return false;
        }

        /**
         * A member function that initializes *this
         */
        void initGram() {
            transMat = NULL;
            _freqs.resize(this->maxNumFilterValues());

            for(int i = 0; i < this->maxNumFilterValues(); i++)
                _freqs.push_back(0);
        }

        /**
         * A member function that builds the table with relative frequencies of the
         * words which form part of *this.
         * The resulting table has dimensions [alphaSize, numHistores].
         */
        void buildCp() {
            int i, j, k, f;

            if(transMat)
                return;

            transMat = new double [this->maxNumFilterValues()]; // r+1 state
            assert(transMat);

            for(i = 0; i < this->maxNumFilterValues(); i++)
                transMat[i] = 0;

            for(f = 0; f < this->maxNumFilterValues(); f++) {
                int ind = f;
                unsigned int freq = _freqs[ind];
                transMat[ind] = freq;
            }

            for(j = 0; j < this->_numHistories; j++) {
                double freqStr;
                freqStr = 0;

                for(k = 0; k < this->alphabetSize(); k++)
                    freqStr += getCp(k, j);

                for(k = 0; k < this->alphabetSize(); k++) {
                    if(freqStr)
                        setCp(k, j, getCp(k, j) / freqStr);
                    else
                        setCp(k, j, 0);

                }
            }
        }

        /**
         * A member function that stores the relative frequencies of all words
         * that form part of *this
         * @param fd the output file stream
         */
        void storeParams(::ofstream **fs) {
            ::ofstream *file = *fs;
            file->write((char *)&totFreq, sizeof(unsigned int));
            file->write((char *)&this->maxFilterVals, sizeof(int));
            //  cerr << totFreq << " " << distNextWord << " " << maxNumFilterValues() << endl;

            for(int i = 0; i < this->maxNumFilterValues(); i++) {
                file->write((char *)&_freqs[i], sizeof(unsigned int));
                file->write((char *)&transMat[i], sizeof(double));
            }
        }



        /**
         * A member operator which adds the relative frequencies of the words
         * that are common to *this and g
         * @param g a FL_Gram object which needs to have been built using the same
         * alphabet and the same order as *this.
         */
        FL_Gram & operator+=(FL_Gram<TClass> &g) {
            assert(*_alphabet == (*g.alphabet()) && this->order() == g.order());
            /*
             * Merging _freqs...
             */

            for(unsigned int i = 0; i < g.maxNumFilterValues(); i++)
                _freqs[i] += g.getWordFreq(i);

            if(transMat) {
                delete transMat;
                transMat = NULL;
                buildCp();
            }

            return *this;

        }

        ~FL_Gram() {
            if(transMat)
                delete [] transMat;
        }

    };


    template <class TClass1, class TClass2>
        class FL_FastGram : public FL_BaseGram<TClass1, TClass2> {
    protected:

    public:
    FL_FastGram(int fInd,
                std::string & name,
                int order,
                int alphabetSize,
                int period = 0) :
        FL_BaseGram<TClass1, TClass2>(fInd,
                                      name,
                                      order,
                                      alphabetSize,
                                      period) {

        }

    FL_FastGram(int fInd,
                vector<std::string> & params,
                int & offset,
                ResourceEngine *re,
                FilterEngine *fe) :
        FL_BaseGram<TClass1, TClass2>(fInd, params, offset) {

            if(!sscanf(params[offset++].c_str(), "%d", &this->_alphabetSize))
                assert(0);

            if(offset < params.size())
                if(!sscanf(params[offset++].c_str(), "%u", &this->_period))
                    assert(0);

            this->maxFilterVals = (int)(pow((double)this->alphabetSize(), this->order()));
            this->_numHistories = (int)(pow((double)this->alphabetSize(), this->order() - 1));

        }

        void computeVals(char *inputSeq, SPARSE_HASH<UCHAR, int> *vals, int len) {
            TClass1 *seq = (TClass1 *)inputSeq;
            int i = 1;
            assert(this->order() > 0);
            typename SPARSE_HASH<UCHAR, int>::iterator it;

            for( ; i + (int)this->order() - 1 <= len; i++) {
                UCHAR val = (UCHAR)wmap(seq, i - 1, i + this->order() - 2);
                it = vals[i].find(val);
                if(it ==  vals[i].end())
                    vals[i][val] = 1;
                else
                    it->second++;

                //	cerr << "\t" << i << "->" << (int)val << "=" << vals[i][val] << endl;

            }
        }

        inline int concat(int word1, TClass1 *word2, int len2) {
            assert(word1 >= 0);
            int res = word1;

            for(int i = 0; i < len2; i++)
                res = res*this->alphabetSize() + word2[i];

            assert(res < this->maxNumFilterValues());
            return res;
        }

        inline int wmap(TClass1 *seq, int beg, int end) {
            return concat(0, seq + beg, end - beg + 1);
        }

        inline int wmap(TClass1 *seq, int wlen) {
            return concat(0, seq, wlen);
        }

        ~FL_FastGram() {}
    };


    template <class TClass> struct HBlock {
        int min;
        int max;
        TClass val;
        int length;
        HBlock() {
            min = max = length = 0;
            val = TClass();
        }
        HBlock(int min, int max, TClass val, int length) {
            this->min = min;
            this->max = max;
            this->val = val;
            this->length = length;
        }
    };

    /*! A Gram that counts only grams formed
     * from different characters and organizes them into blocks
     */
    template <class TClass1, class TClass2>
        class FL_SparseGram : public FL_FastGram<TClass1, TClass2> {
    protected:
        int _maxPeriodicity;
        vector<HBlock<TClass2> > hblocks[NUM_STRANDS];
        typename vector<HBlock<TClass2> >::iterator curr_iter[NUM_STRANDS];

    public:
    FL_SparseGram(int fInd,
                  std::string & name,
                  int order,
                  int alphabetSize,
                  int maxPeriodicity = 0)
        :    FL_FastGram<TClass1, TClass2>(fInd,
                                           name,
                                           order,
                                           alphabetSize,
                                           0) {
            _maxPeriodicity = maxPeriodicity;

        }

    FL_SparseGram(int fInd,
                  vector<std::string> & params,
                  int & offset,
                  ResourceEngine *re,
                  FilterEngine *fe)
        : FL_FastGram<TClass1, TClass2>(fInd, params, offset,
                                        re, fe) {

            _maxPeriodicity = this->_period;
            this->_period = 0; // doesn't make sense to accumulate a sparse filter
        }

        inline int maxPeriodicity() {
            return _maxPeriodicity;
        }

        /**
         * A member function for computing filter values.
         *
         * The value at each position of the input sequence is the gram found in the
         * sequence at the position.
         *
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */

        inline void computeVals(char *inputSeq, TClass2 *vals, int len) {
            constructHBlocks(inputSeq, len);
        }

        inline typename vector<HBlock<TClass2> >::iterator findHBlock(int pos,
                                                                      TStrand strand) {

            while(curr_iter[strand] != hblocks[strand].begin() &&
                  pos < curr_iter[strand]->max)
                --curr_iter[strand];

            while(curr_iter[strand] != hblocks[strand].end() &&
                  pos > curr_iter[strand]->max)
                ++curr_iter[strand];

            //      cerr << "found block for " << pos << " (" << (ret.first)->first.min << "," << (ret.first)->first.max << ")= " << (int)(ret.first)->second.val << " " << (ret.first)->second.length << endl;
            return  curr_iter[strand];

        }

        inline typename vector<HBlock<TClass2> >::iterator end(TStrand strand) {
            return hblocks[strand].end();
        }

        /**
         * Construct homogenous blocks of grams.
         */
        inline void constructHBlocks(char *inputSeq, int len) {
            TClass1 *seq = (TClass1 *)inputSeq;
            unsigned int i = 1;
            assert(this->order() > 0);

            vector<int> pb(this->maxNumFilterValues(), -1);

            for( ; i + this->order() - 1 <= (unsigned int)len; i++)  {
                TClass2 fval = (TClass2)this->wmap(seq, i - 1, i + this->order() - 2);

                if(pb[fval] >= 0) {
                    HBlock<TClass2> &b = hblocks[this->defaultStrand][pb[fval]];

                    if(i - 1 > b.max + this->maxPeriodicity()) {
                        //    cerr << "ending " << (int)fval << "-block (" << b.min << "," << b.max << ")= " << b.length << endl;
                        pb[fval] = -1;
                    }
                    else {
                        b.max = i;
                        b.length++;
                    }
                }

                if(pb[fval] < 0) {
                    pb[fval] = hblocks[this->defaultStrand].size();
                    hblocks[this->defaultStrand].push_back(HBlock<TClass2>(i, i, fval, 1));

                    //  cerr << "inserting " << (int)fval << "-block #" << hblocks[this->defaultStrand].size() - 1 << " at pos " << b[fval]->min << endl;
                }
            }

            //      for(int i = 0; i < hblocks[this->defaultStrand].size(); i++) {
            //HBlock<TClass2> & b = hblocks[this->defaultStrand][i];
            //cerr << (int)b.val << "-block #" << i << " (" << b.min << "," << b.max << ") = " << b.length << endl;
            //      }
            curr_iter[this->defaultStrand] = hblocks[this->defaultStrand].begin();

        }

        void freeValArrays(TStrand strand) {
            hblocks[strand].clear();
        }

        ~FL_SparseGram() {}
    };


    /*************************************************************************
     * FL_ConsGram class
     *************************************************************************/

    template <class TClass>
        class FL_ConsGram : public FL_BaseGram<char, TClass> {
    protected:
        Sigma *_alphabet;
        int _numWords;
        int _numGapWords;
        int _numUAWords;
        char _uaSymbol;
        Sigma *_gapSigma;

    public:
    FL_ConsGram(int fInd,                 //!< A unique identifier.
                std::string & name,       //!< A unique name.
                Sigma *alphabet,          //!< The FL_ConsGram's alphabet.
                Sigma *gapSigma,                //!< gap alphabet
                char uaSymbol,               //!< unaligned symbol
                int order,                //!< The FL_ConsGram's order.
                int period = 0
        ) :
        FL_BaseGram<char, TClass>(fInd,
                                  name,
                                  order,
                                  alphabet->alphabetSize(),
                                  false,
                                  period) {

            _alphabet = alphabet;
            _gapSigma = gapSigma;
            _uaSymbol = uaSymbol;
            initConsGram();
        }


    FL_ConsGram(int fInd,           //!< A unique identifier.
                /*! The Header string definition, loaded as a vector
                 * of strings.
                 * The Header has the following form:\n\n
                 * Filter name FL_ConsGram arrInd sigma order
                 * compactHomogBlocks\n\n
                 * The description of the fields could be found in
                 * the other constructor(s)
                 */
                vector<std::string> & params,
                int & offset,       //!< The index for vector params.
                ResourceEngine *re, //!< A pointer to the ResourceEngine object.
                FilterEngine *fe    //!< A pointer to the FilterEngine object.
        ) :
        FL_BaseGram<char, TClass>(fInd, params, offset) {

            _alphabet = (Sigma *)re->getResource(params[offset++]);
            _gapSigma = (Sigma *)re->getResource(params[offset++]);
            if(!sscanf(params[offset++].c_str(), "%c", &_uaSymbol))
                assert(0);

            this->_period = 0;
            if(offset < params.size())
                if(!sscanf(params[offset++].c_str(), "%u", &this->_period))
                    assert(0);

            this->_alphabetSize = _alphabet->alphabetSize();
            this->maxFilterVals = (int)(pow((double)this->alphabetSize(), this->order()));
            this->_numHistories = (int)(pow((double)this->alphabetSize(), this->order() - 1));

            initConsGram();

        }

        inline int wmap(char *seq, int beg, int end) {
            return wmap(seq + beg, end - beg + 1);
        }

        inline int wmap(char *seq, int wlen) {
            int ind = alphabet()->strictIndS(seq, wlen);
            if(ind >= 0)
                return ind;

            ind = _gapSigma->strictIndS(seq, wlen);
            if(ind >= 0)
                return _numWords + ind;
            ind = _numWords + _numGapWords;
            int uaInd = 0;
            if(seq[0] != _uaSymbol)
                uaInd++;
            if(wlen > 1 && seq[wlen - 1] != _uaSymbol)
                uaInd += 2;

            return ind + uaInd;
        }

        inline Sigma *alphabet() {
            return _alphabet;
        }


        /**
         * A member function that initializes *this
         */
        void initConsGram() {

            int gapSigmaSize = _gapSigma->alphabetSize();

            _numWords = this->maxFilterVals;
            _numGapWords = (int)(pow((double)gapSigmaSize, this->_order));
            _numUAWords = 2 + 2*(this->_order > 1);

            this->_alphabetSize += gapSigmaSize + 2;
            this->maxFilterVals += _numGapWords + _numUAWords;
            this->_numHistories += (int)(pow((double)gapSigmaSize, this->_order - 1));
            this->_numHistories += 1 + 1*(this->_order > 1);
        }

        ~FL_ConsGram() {}

    };

}

#endif
