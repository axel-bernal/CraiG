/****************************************************************************
 * FL_CountingScore.h - part of the lless namespace, a general purpose
 *                      linear semi-markov structure prediction library
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

#ifndef _FILTER_COUNTINGSCORE_H_
#define _FILTER_COUNTINGSCORE_H_

#include "Utils.h"
#include "FL_Score.h"
#include "ContextIMM.h"
#include "FilterEngine.h"
#include "ResourceEngine.h"
#include "IMM.h"
#include "Motif.h"

namespace lless {

    /**
     *
     * FL_CountingScore is a subclass of FL_Score. The scoring function here
     * counts a single pattern when found in the input sequence.
     *
     */
    class FL_CountingScore : public FL_Score<int> {

    protected:
        int wIndex;            // symbol to look <-- sigma->ind(pattern, order)
        unsigned int order; /* the ind refers to patterns of length order
                               in the original seq. */
        Sigma *sigma;

    public:
        //! Default contructor
    FL_CountingScore(int fInd,            //!< A unique identifier.
                     std::string & name,  //!< A unique name.
                     Sigma *sigma,        //!< The alphabet used for pattern
                     char *pattern        //!< The word to count
        )
        :  FL_Score<int>(fInd,
                         name,
                         1,
                         0, FT_INTEGER
            ) {

            this->sigma = sigma;
            this->order = strlen(pattern);
            this->wIndex = sigma->indS(pattern, this->order);
        }

        /**
         * Constructor from a Header string definition
         */
    FL_CountingScore(int fInd,        //!< A unique identifier.
                     /*! The Header string definition, loaded as a vector
                      * of strings.
                      * The Header has the following form:\n\n
                      * Filter  name CountingScore arrInd period
                      * smoothenWSize sigma "pattern"|ind\n\n
                      * The last parameter could be either a pattern enclosed
                      * in quotes or a symbol of order 1.
                      * The description of the other fields could be found in
                      * the other constructor(s)
                      */
                     vector<std::string> & params,
                     int & offset,       //!< The index for vector params.
                     ResourceEngine *re, /*!< A pointer to the
                                           ResourceEngine object. */
                     FilterEngine *fe    /*!< A pointer to the FilterEngine
                                           object. */
        )
        : FL_Score<int>(fInd,
                        params,
                        offset,
                        fe,
                        FT_INTEGER) {

            this->sigma  = (Sigma *)re->getResource(params[offset++]);
            assert(sigma);

            if(params.size() == (unsigned)offset)
                assert(0);

            if(!sscanf(params[offset].c_str(), "%d", &wIndex)) {
                assert(params[offset][0] == '\"'
                       && params[offset][params[offset].length() - 1] == '\"');
                this->order = params[offset].length() - 2;
                this->wIndex = sigma->indS((char *)params[offset].c_str() + 1, this->order);
            }
            else
                order = 1;

            offset++;
        }

        /**
         * The implemented scoring function. At each position p the score is 1
         * if sigma->ind(seq[p..p+order]) matches ind or zero otherwise.
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        void scoreInputSeq(char *seq, int *filterVals, int len) {
            //      cerr << this->getName() << " " << wIndex << "\n";
            for(int i = 0; i < len; i++) {
                filterVals[i + 1] = (sigma->indS(seq + i, order)  == wIndex) ? 1 : 0;
                //        if(wIndex == 7 && filterVals[i + 1])
                //          cerr << i << " " << seq[i] << " " << filterVals[i + 1] << endl;
            }
        }

        ~FL_CountingScore() {

        }

    };


    /**
     *
     * FL_MultiCountingScore is a subclass of FL_Score. The scoring function here
     * counts occurrences of multiple pattern when found in the input sequence.
     *
     */
    class FL_MultiCountingScore : public FL_Score<int> {

    protected:
        int *wIndexes;            // symbols to look <-- sigma->ind(pattern, order)
        unsigned int order;       /* the ind refers to patterns of length order
                                     in the original seq. */
        Sigma *sigma;

    public:
        //! Default contructor
    FL_MultiCountingScore(int fInd,            //!< A unique identifier.
                          std::string & name,  //!< A unique name.
                          Sigma *sigma,        //!< The alphabet used for pattern
                          vector<string> patterns        //!< The word to count
        )
        :  FL_Score<int>(fInd,
                         name,
                         1,
                         0, FT_INTEGER
            ) {

            this->sigma = sigma;
            assert(patterns.size() > 0);
            this->order = patterns[0].size();
            initWordIndexes();
            for(int i = 0; i < patterns.size(); i++)
                this->wIndexes[sigma->indS((char *)patterns[i].c_str(), this->order)] = 1;
        }

        /**
         * Constructor from a Header string definition
         */
    FL_MultiCountingScore(int fInd,        //!< A unique identifier.
                          /*! The Header string definition, loaded as a vector
                           * of strings.
                           * The Header has the following form:\n\n
                           * Filter  name CountingScore arrInd period
                           * smoothenWSize sigma "pattern"|ind\n\n
                           * The last parameter could be either a pattern enclosed
                           * in quotes or a symbol of order 1.
                           * The description of the other fields could be found in
                           * the other constructor(s)
                           */
                          vector<std::string> & params,
                          int & offset,       //!< The index for vector params.
                          ResourceEngine *re, /*!< A pointer to the
                                                ResourceEngine object. */
                          FilterEngine *fe    /*!< A pointer to the FilterEngine
                                                object. */
        )
        : FL_Score<int>(fInd,
                        params,
                        offset,
                        fe,
                        FT_INTEGER) {

            this->sigma  = (Sigma *)re->getResource(params[offset++]);
            assert(sigma);

            int count;
            if(!sscanf(params[offset++].c_str(), "%d", &count))
                assert(0);

            if(params.size() == (unsigned)offset)
                assert(0);
            this->wIndexes = NULL;

            for(int i = 0; i < count; i++) {
                int val;
                if(!sscanf(params[offset].c_str(), "%d", &val)) {
                    assert(params[offset][0] == '\"'
                           && params[offset][params[offset].length() - 1] == '\"');
                    this->order = params[offset].length() - 2;
                    val = sigma->indS((char *)params[offset].c_str() + 1, this->order);
                }
                else
                    order = 1;

                if(!wIndexes)
                    initWordIndexes();

                wIndexes[val] = 1;
                offset++;
            }
        }

        void initWordIndexes() {
            int numSymbols = sigma->alphabetSize()^this->order;
            wIndexes = new int [numSymbols];

            for(int i = 0; i < numSymbols; i++)
                this->wIndexes[i] = 0;
        }

        /**
         * The implemented scoring function. At each position p the score is
         * wIndexes[sigma->ind(seq[p..p+order])]
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        void scoreInputSeq(char *seq, int *filterVals, int len) {
            //      cerr << this->getName() << " " << ind << "\n";
            for(int i = 0; i < len; i++) {
                filterVals[i + 1] = wIndexes[sigma->indS(seq + i, order)];
                //        cerr << seq[i] << " " << filterVals[i + 1] << endl;
            }
        }

        ~FL_MultiCountingScore() {
            delete [] wIndexes;
            wIndexes = NULL;
        }

    };


    /**
     *
     * FL_ContainScore is a subclass of FL_Score. The scoring function here
     * counts all words of size wordLen which contain pattern or its
     * associated symbol ind as defined in alphabet sigma, both passed
     * as parameter.
     *
     */
    class FL_ContainScore : public FL_CountingScore {

    private:
        unsigned int wordSize; //!< the context length used for searching pattern

    public:
        //! Default contructor
    FL_ContainScore(int fInd,            //!< A unique identifier.
                    std::string & name,  //!< A unique name.
                    Sigma *sigma,        //!< The alphabet used for pattern
                    char *pattern,        //!< The word to count
                    int wordSize          //!< The size of the words to look
        )
        :  FL_CountingScore(fInd,
                            name,
                            sigma,
                            pattern
            ) {

            this->wordSize = wordSize;
        }

        /**
         * Constructor from a Header string definition
         */
    FL_ContainScore(int fInd,        //!< A unique identifier.
                    /*! The Header string definition, loaded as a vector
                     * of strings.
                     * The Header has the following form:\n\n
                     * Filter  name ContainScore arrInd period
                     * smoothenWSize sigma "pattern"|ind wordSize\n\n
                     * The second to last parameter could be either a
                     * pattern enclosed in quotes or a symbol of order 1.
                     * The description of the other fields could be found in
                     * the other constructor(s)
                     */
                    vector<std::string> & params,
                    int & offset,       //!< The index for vector params.
                    ResourceEngine *re, /*!< A pointer to the
                                          ResourceEngine object. */
                    FilterEngine *fe    /*!< A pointer to the FilterEngine
                                          object. */
        )
        : FL_CountingScore(fInd,
                           params,
                           offset,
                           re,
                           fe
            ) {

            if(!sscanf(params[offset++].c_str(), "%d", &this->wordSize))
                assert(0);

        }

        /**
         * The implemented scoring function. At each position p the score is 1
         * if seq[p..p+order] contains symbol
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        void scoreInputSeq(char *seq, int *filterVals, int len) {
            //      cerr << this->getName() << " " << ind << "\n";
            int accum = 0, i = 0;
            int limit = Utils::min(wordSize, len);
            for( ; i < limit; i++)
                accum += (sigma->indS(seq + i, order)  == wIndex);
            for(i = 0; i < len; i++) {
                filterVals[i + 1] = (accum > 0);
                accum -= (sigma->indS(seq + i, order)  == wIndex);
                if(wordSize + i + order <= strlen(seq))
                    accum += (sigma->indS(seq + wordSize + i, order)  == wIndex);
                assert(accum >= 0);
                //        cerr << seq[i] << " " << filterVals[i + 1] << endl;
            }

        }

        ~FL_ContainScore() {

        }

    };

}

#endif
