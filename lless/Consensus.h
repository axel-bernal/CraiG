/****************************************************************************
* Consensus.h - part of the lless namespace, a general purpose
*               linear semi-markov structure prediction library
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

#ifndef _PATTERN_CONSENSUS_H
#define _PATTERN_CONSENSUS_H
#include "Utils.h"
#include "Sigma.h"
#include "PrefixTree.h"
#include "NGram.h"
#include <map>


namespace lless {

  /**
   * The Consensus class represents a consensus sequence. The class will match 
   * any pattern that falls in the OR space defined by the matched(added)
   * patterns.
   *
   ***************************************************************************/

  class Consensus {

   private:
    bool ***cache;
    bool useCache;
    vector<bool *> mat; //!< the consensus sequence
    Sigma *sigma; 
    int order;
    int maxNumValues;

   public:
    /** 
     * Constructor
     * @param pattern the first consensus pattern that is to be matched(added)
     * to the consensus sequence. It can contain ambiguous characters.
     * @param useCache caches results from function matches(int, int, int)
     * if true. Use it only if the gram's alphabet is small enough that 
     * generating all possible words of the same length as the consensus
     * sequence is not prohibitively expensive.
     */
    Consensus(char *pattern, Sigma *sigma,
	      int order, int maxNumValues,
	      bool useCache = true) {

      unsigned int length = strlen(pattern);
      assert(length >= (unsigned)order);
      this->useCache = useCache;
      this->sigma = sigma;
      this->order = order;
      this->maxNumValues = maxNumValues;
      initializeCache(length);
      matchPattern(pattern);

    }

    /**
     * A member function that creates and initializes a cache table to store
     * all possible words which have the same length as the consensus 
     * sequence. 
     */
    inline void initializeCache(unsigned int length) {
      unsigned int i,j;
      cache = new bool ** [length];

      for(i = 0; i < length; i++) {
        mat.push_back(new bool [sigma->alphabetSize()]);
        for(j = 0; (int)j < sigma->alphabetSize(); j++)
          mat[i][j] = false;

        cache[i] = new bool *[length];
        for(j = 0; j < length; j++)
          cache[i][j] = NULL;

        if(!useCache)
          continue;

        cache[i][order - 1] = new bool [maxNumValues];
        for(int k = 0; k < maxNumValues; k++)
          cache[i][order - 1][k] = false;
      }
    }

    /**
     * A member function that matchs(adds) pattern to the consensus sequence; 
     * patterm must have the same length as any other previously
     * added patterns. If usecache is true, the function updates the cache
     * with all the words that pattern matches. The resulting consensus 
     * sequence is the OR of all matched(added) patterns.
     */ 
    void matchPattern(char *pattern) {
      assert(strlen(pattern) == (unsigned)size());
      int i;

      int *ba = new int[sigma->alphabetSize()];

      for(i = 0; i < sigma->alphabetSize(); i++)
        ba[i] = 0;

      for(i = 0; i < size(); i++) {
        int numAmbChars = sigma->indC(ba, pattern[i]), j;

        for(j = 0; j < numAmbChars; j++)
          mat[i][ba[j]] = true;

      }

      if(useCache) {
        for(i = 0; i < size(); i++) {      
          for(int k = 0; k < maxNumValues; k++)
            cache[i][order - 1][k] = _matches(k, i, order);
        }
      }
      delete [] ba;
    }

    inline int size() {
      return mat.size();
    }

    /**
     * @return true if word matches the consensus sequence defined by all 
     * patterns previously added through function matchPattern(char *),
     * false otherwise.
     */
    inline bool matches(int word, int pos, int wordLen) {
      assert(wordLen != 0);

      if(cache[pos][wordLen - 1])
        return cache[pos][wordLen - 1][word];

      return _matches(word, pos, wordLen);
    }

   protected:
    inline bool _matches(int word, int pos, int wordLen) {
      assert(pos >= 0);
      if(pos + wordLen > (int)mat.size())
        return false;
      for(int i = pos + wordLen - 1 ; i >=  pos; i--) {
        if(!mat[i][sigma->wordTail(word)])
          return false;
        word = sigma->wordHead(word);
      }
      assert(word == 0);
      return true;
    }                                                                            

   public:
    ~Consensus() {
      for(unsigned i = 0; i < mat.size(); i++) {
        if(mat[i])
          delete [] mat[i];
        mat[i] = NULL;

        for(unsigned j = 0; j < mat.size(); j++) {

          if(cache[i][j])
            delete [] cache[i][j];

          cache[i][j] = NULL;
        }
        if(cache[i])
          delete [] cache[i];
      }
      if(cache)
        delete [] cache;
      cache = NULL;
    }
  };


  /**
   * The XORConsensus class will match any pattern that coincides with one 
   * and only one matched pattern. It uses a PrefixTree representation for
   * storing the matched patterns.
   **************************************************************************/

  class XORConsensus {
   protected:
    PrefixTree *root;
    Sigma *sigma;
    int _size;
   public:
    /**
     * Constructor
     * @param pattern a pattern to add to the consensus. It may contain
     * ambiguous characters.
     * @param sigma an alphabet compatible with pattern.
     */
    XORConsensus(char *pattern, Sigma *sigma) {
      this->sigma = sigma;
      root = new PrefixTree(sigma->alphabetSize());
      _size = strlen(pattern);
      matchPattern(pattern);
    }

    inline int size() {  return _size; }


    /**
     * A member function that matchs(adds) pattern to the consensus sequence; 
     * patterm may have a different length from  other previously
     * added patterns. The resulting consensus sequence is the XOR of 
     * all matched(added) patterns.
     */ 
    inline void matchPattern(char *pattern) {
      _matchPattern(root, pattern);
    }

    /**
     * A function that matches(adds) pattern to the consensus by recursively
     * going from left to right through the pattern and adding a child to the
     * PreffixTree object at each recursion level.
     */
    void _matchPattern(PrefixTree *node, char *pattern) {
      int j;

      if(!strlen(pattern))
        return;

      int *ba = new int[sigma->alphabetSize()];
      int numAmbChars = sigma->indC(ba, pattern[0]);

      for(j = 0; j < numAmbChars; j++) 
        _matchPattern(node->addChild(ba[j]), pattern + 1);
      delete [] ba;
    }

    /**
     * @return 1 if word matches any of the patterns previously added
     * through function matchPattern(char *), 0 otherwise.
     */
    int matches(char *word, int wordLen) {
      if(wordLen == 0)
        return 0;
      PrefixTree *ptr = root;
      int indC;
      for(int i = 0; i < wordLen; i++) {
        //      cerr << "->" << word[i];
        indC = sigma->indC(word[i]);
        if(ptr == NULL || indC < 0)
          return 0;
        ptr = (*ptr)[indC];
      }
      return ptr == NULL ? 0 : 1;
    }

    ~XORConsensus() {
      delete root;
    }

  };

}
#endif
