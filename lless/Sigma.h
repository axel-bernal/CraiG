/****************************************************************************
 * Sigma.h - part of the lless namespace, a general purpose
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

#ifndef _SIGMA_H_
#define _SIGMA_H_
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include "ResourceEngine.h"


#define AA_ALPHABET_SIZE 21
#define DNA_ALPHABET_SIZE 4
#define DNA_CHARS "ACGTMRWSYKBDHVNXacgtmrwsykbdhvnx"
#define COMP_DNA_CHARS "TGCAKYWSRMVHDBNXtgcakywsrmvhdbnx"
#define AA_CHARS "*ACDEFGHIKLMNPQRSTVWY*acdefghiklmnpqrstvwy"
#define COMP_AA_CHARS "*ACDEFGHIKLMNPQRSTVWY*acdefghiklmnpqrstvwy"

using namespace lless;

namespace lless {

    class ResourceEngine;

    /**
     * Sigma is a generic alphabet implementation. The user can define
     * sets of characters in which all members share the same symbol. One
     * important shortcomming is that the class only works efficiently on
     * alphabets with not too many symbols (<= 8)
     *
     ***************************************************************************/

    class Sigma : public Resource {
    protected:
        int ind[128];            //!< The index array, any entry != 0 is a symbol
        int _alphabetSize;        //!< The alphabet size
        std::string *charsPerSymbol;       //!< The characters allowed per symbol
        std::string symbols;     //!< the alphabet of valid characters
        std::string compSymbols; //!< the complementary valid alphabet characters

    public:
    Sigma() : Resource() {
            charsPerSymbol = NULL;
        }

    Sigma(int alphabetSize,       //!< The alphabet size.
          //! An ordered set of characters, regarded as valid by the alphabet.
          const char *symbols,
          /*! the ordered complementary set of characters regarded as valid
            by the alphabet. */
          const char *compSymbols,
	  const char *chrs4Symbol = NULL
	  /*  Defines the set of characters that each alphabet
           *  symbol encodes for, as they appear in the input sequence.
	   */
        ) : Resource() {

            initializeSigma(alphabetSize, string(symbols), string(compSymbols));

            if(chrs4Symbol) {
                string s(chrs4Symbol);
                boost::RegEx rExDelim("\\s+");
                std::vector<std::string> vs;
                rExDelim.Split(vs, s, 0, 100);
                assert(vs.size() == alphabetSize);

                for(int i = 0; i < alphabetSize; i++) {
                    for(unsigned int j = 0; j < vs[i].length(); j++)
                        addInd(vs[i][j], i);
                }
            }
        }

        /**
         * Constructor from a Header string definition which is not used.
         */
    Sigma(
        /*! The Header string definition, loaded as a vector of strings.
          it should be empty.
        */
        std::vector<std::string> &params,
        int & offset,       //!< The index for vector params.
        int alphabetSize,   //!< The alphabet size.
        //! An ordered set of characters, regarded as valid by the alphabet.
        const char *symbols,
        /*! the ordered complementary set of characters regarded as valid
          by the alphabet. */
        const char *compSymbols
        ) : Resource(params, offset, false) {

            initializeSigma(alphabetSize, string(symbols), string(compSymbols));

        }

        /**
         * Constructor from a Header string definition.
         */
    Sigma(
        /*! The Header string definition, loaded as a vector of strings.
         *  The Header has the following form:\n
         *  Resource name Sigma alphabetSize symbols  compSymbols
         *  [chars4Symbol_1 ..] \n
         *  The field(s) chars4Symbol_i, define the set of characters that
         *  symbol i encodes for, as they appear in the input sequence.
         */
        std::vector<std::string> &params,
        int & offset,       //!< The index for vector params.
        ResourceEngine *re  //!< A pointer to the ResourceEngine object.
        ) : Resource(params, offset, false) {

            if(!sscanf(params[offset++].c_str(), "%d", &_alphabetSize))
                assert(0);

            initializeSigma(_alphabetSize, params[offset], params[offset+1]);
            offset += 2;
            for(int i = 0; i < alphabetSize(); i++) {
                for(unsigned int j = 0; j < params[offset].length(); j++)
                    addInd(params[offset][j], i);
                offset++;
            }
        }

        inline std::string & getSymbols() {
            return this->symbols;
        }

        /**
         * A member function to initialize
         * @param alphabetSize the alphabet size
         * @param symbols the string containing the legal alphabet symbols
         * @param compSymbols the string containing the legal complementary
         * symbols
         */
        inline void initializeSigma(int alphabetSize,
                                    string symbols,
                                    string compSymbols) {

            int i;
            _alphabetSize = alphabetSize;
            for(i = 0; i < 128; i++)
                ind[i] = -1;
            charsPerSymbol = new std::string [alphabetSize + 1];
            this->symbols = symbols;
            this->compSymbols = compSymbols;
        }

        /**
         * A member function which adds to the alphabet validated character c
         * as symbol id
         * @param c the character to be added
         * @param id the symbol which encodes character c
         */
        void addInd(char c, int id) {
            if(c == '.') {
                for(unsigned int i = 0; i < symbols.length(); i++)
                    ind[(int)symbols[i]] = id;
            }
            else
                ind[(int)c] = id;
            charsPerSymbol[id].push_back(c);
        }

        /**
         * A member function which saves the Header definition.
         * @see Resource::saveHeader(std::ofstream &)
         */
        inline void saveHeader(std::ofstream &fd) {
            Resource::saveHeader(fd);

            fd << " " << alphabetSize() << " " << symbols << " " << compSymbols;
            for(int i = 0; i < alphabetSize(); i++)
                fd << " " << charsPerSymbol[i];
            fd << endl;
        }

        /**
         * @return the complement of c
         * @param c the character which we need the complement from
         */
        inline char toComplement(char c) {
            string::size_type pos = symbols.find(c);
            if(pos == string::npos) {
                assert(0);
                throw EXCEPTION(STRANGE_CHAR, string(1, c));
            }
            return compSymbols[pos];
        }

        bool operator==(Sigma &s) {
            bool sameSig = true;
            for(int i = 0; i < alphabetSize(); i++)
                sameSig = sameSig && (charsPerSymbol[i].compare(s.getCharsPerSymbol(i)) == 0);

            return sameSig && (alphabetSize() == s.alphabetSize());
        }

        inline int alphabetSize() {
            return _alphabetSize;
        }

        /**
         * @param word a word defined in *this alphabet of unknown length
         * @return the last symbol of word
         */
        inline int wordTail(int word) {
            return word % alphabetSize();
        }

        /**
         * @param word a word defined in *this alphabet of unknown length
         * @return the first symbol(s) of word
         */
        inline int wordHead(int word) {
            return word / alphabetSize();
        }

        /**
         * @param word a word defined in *this alphabet of unknown length
         * @param suffix the symbol to add at the end of the word
         * @return a new word which encodes the concatenation of word and suffix
         */
        inline int addSuffix(int word, int suffix) {
            return alphabetSize()*word + suffix;
        }

        /**
         * @param word a word defined in *this alphabet of unknown length
         * @param order the length of the word.
         * @return a new word which results from deleting prefix from word;
         * prefix is the alphabet symbol appearing at the beginning of the
         * word given as parameter.
         */
        inline int removePrefix(int word, int order) {
            if(!order) return 0;
            int quanta = (int)pow((double)alphabetSize(), (double)order - 1);
            return word - quanta*(word/quanta);
        }

        /**
         * @param word a word defined in *this alphabet of unknown length.
         * @param order the length of the word.
         * @param prefix is the word of unknown length appearing at the beginning
         * of word.
         * @return a new word which results from deleting prefix from word.
         */
        inline int removePrefix(int word, int prefix, int order) {
            if(!order) return 0;
            int quanta = (int)pow((double)alphabetSize(), (double)order - 1);
            return word - quanta*prefix;
        }

        /**
         * @return the set of characters which are encoded by symbol i.
         * @param i the symbol we want to inquiry about.
         */
        inline std::string & getCharsPerSymbol(int i) {
            return charsPerSymbol[i];
        }

        /**
         * @return the alphabet symbol encoding character ch.
         * @param ch the input character.
         */
        inline int indC(char ch) {
            return ind[(int)ch];
        }

        /**
         * @return the word encoding string str.
         * @param str the input string.
         * @param len the length of the string str.
         */
        inline int randIndS(char *str, unsigned int len) {
            assert(len <= strlen(str));
            return concat(0, str, len);
        }

        /**
         * @return the word that encodes str[beg..end].
         * @param str the input string.
         * @param beg the starting point in str
         * @param end the end point in str.
         */
        inline int randIndS(char *str, int beg, int end) {
            assert(beg >= 0);
            return concat(0, str + beg, (unsigned int)(end - beg + 1));
        }

        /**
         * A member function that concatenates words word1 and word2 of length
         * len2.
         * @param word1 the word that goes at the beginning.
         * @param word2 the word that goes at the end.
         * @param len2 the length of word2.
         * @return the word  that encodes the concatenation
         * of contS(word1) + word2.
         */
        inline int concat(int word1, char *word2, unsigned int len2) {
            assert(len2 <= strlen(word2) && word1 >= 0);
            unsigned int i;
            int res = word1;

            for(i = 0; i < len2; i++)
                if(word2 + i) {
                    int baseCode =  indC(word2[i]);
                    if(baseCode < 0)
                        baseCode = indC(charsPerSymbol[(int)(alphabetSize()*rand()/(RAND_MAX + 1.0))][0]);
                    //      cout << str <<" " << i << "\n";
                    res = res*alphabetSize() + baseCode;
                }
            //    words[s] = res;
            //    str[end + 1] = saved;
            return res;
        }


        /**
         * @return the word encoding string str.
         * @param str the input string.
         * @param len the length of the string str.
         */
        inline int strictIndS(char *str, unsigned int len) {
            assert(len <= strlen(str));
            return strictConcat(0, str, len);
        }

        /**
         * @return the word that encodes str[beg..end].
         * @param str the input string.
         * @param beg the starting point in str
         * @param end the end point in str.
         */
        inline int strictIndS(char *str, int beg, int end) {
            assert(beg >= 0);
            return strictConcat(0, str + beg, (unsigned int)(end - beg + 1));
        }

        /**
         * A member function that concatenates words word1 and word2 of length
         * len2.
         * @param word1 the word that goes at the beginning.
         * @param word2 the word that goes at the end.
         * @param len2 the length of word2.
         * @return the word  that encodes the concatenation
         * of contS(word1) + word2.
         */
        inline int strictConcat(int word1, char *word2, unsigned int len2) {
            assert(len2 <= strlen(word2) && word1 >= 0);
            unsigned int i;
            int res = word1;

            for(i = 0; i < len2; i++)
                if(word2 + i) {
                    int baseCode =  indC(word2[i]);
                    if(baseCode < 0)
                        return -1;
                    //      cout << str <<" " << i << "\n";
                    res = res*alphabetSize() + baseCode;
                }
            //    words[s] = res;
            //    str[end + 1] = saved;
            return res;
        }

        /**
         * @return the word the encodes str
         * @param str the input string
         * @param len the input string's length
         */
        inline int indS(char *str, unsigned int len) {
            unsigned int i, res=0;
            if(!len) return 0;
            assert(len <= strlen(str));
            for(i = 0; i < len; i++)
                if(str + i) {
                    int baseCode = indC(str[i]);
                    if(baseCode < 0)
                        baseCode = indC(charsPerSymbol[(int)(alphabetSize()*rand()/(RAND_MAX + 1.0))][0]);

                    res = res*alphabetSize() + baseCode;
                }
            return res;
        }

        /**
         * @param ind the symbol we want more information about.
         * @return the character which is encoded by symbol ind. If
         * symbol is encoded by more than one character, the first character in the
         * set is chosen (since it is an ordered set)
         */
        char contC(int ind) {
            assert(ind >= 0 && ind < alphabetSize());
            return charsPerSymbol[ind][0];
        }

        /**
         * @return the string of characters encoded by the word ind of length
         * len.
         * @param s the output string which is encoded by input word ind.
         * @param ind the input word.
         * @param len the length of ind.
         */
        char * contS(char *s, int ind, int len) {
            s[len] = '\0';
            if(!len) return NULL;
            for(int i=0;i < len;i++) {
                s[len - (i + 1)] = contC(wordTail(ind));
                ind = wordHead(ind);
            }
            return s;
        }

        /**
         * @return true if c is a valid character, false otherwise
         */
        inline bool checkSymbol(char c) {
            string::size_type pos = symbols.find(c);

            if(pos == string::npos)
                return false;
            return true;
        }

        /**
         * @return true if s is a valid string, false otherwise
         */
        virtual bool checkSeq(char *s) {
            int i =0;
            int len = strlen(s);
            while(i < len)
                if(!checkSymbol(s[i++]))
                    return false;
            return true;
        }

        /**
         * @see checkSeq(char *)
         */
        virtual bool checkSeq(std::string & s) {
            return checkSeq((char *)s.c_str());
        }

        /**
         * A member function that deals with an ambiguos character which have more
         * than one symbol in the alphabet that could encode for it.
         * @param ba the ordered set of all possible symbols which encode
         * ambiguous character ch.
         * @param ch the input ambiguous character.
         * @return the number of symbols which encode the ambiguous
         * character ch.
         */
        virtual inline int indC(int *ba, char ch) {
            ba[0] = indC(ch);
            return 1;
        }

        /**
         * A member function that deals with strings which have at least one
         * ambiguous character.
         * @param sa the ordered set of all possible words which could encode str
         * @param str the input ambiguous string
         * @param len the length of str
         * @return the number of ambiguous strings
         */
        virtual int indS(int *sa, char *str, unsigned int len) {
            sa[0] = indS(str, len);
            return 1;
        }

        /**
         * @return the number of ambiguous strings
         * @see indS(int *, char *, unsigned int)
         */
        virtual int lenIndS(char *str, unsigned int len) {
            return 1;
        }

        inline bool operator==(const Sigma &sgm) {
            Sigma &sigma = (Sigma &)sgm;
            return alphabetSize() == sigma.alphabetSize() &&
            !symbols.compare(sigma.getSymbols());
        }

        virtual ~Sigma() {
            assert(charsPerSymbol);
            delete [] charsPerSymbol;
            charsPerSymbol = NULL;
        }
    };

    /**
     * DNASigma is a hardcoded Sigma object with all the properties of an
     * alphabet to recognize DNA sequences, in particular it deals with
     * IUPAC ambiguous characters.
     */
    class DNASigma : public Sigma {
    public:
    DNASigma() : Sigma(DNA_ALPHABET_SIZE, DNA_CHARS, COMP_DNA_CHARS) {
            initialize();
        }

    DNASigma(std::vector<std::string> &params, int & offset, ResourceEngine *re) : Sigma(params, offset, DNA_ALPHABET_SIZE, DNA_CHARS, COMP_DNA_CHARS) {
            initialize();
        }

        /**
         * A member function that initializes an alphabet for DNA sequences
         */
        void initialize() {
            addInd('A', 0); addInd('a', 0);
            addInd('C', 1); addInd('c', 1);
            addInd('G', 2); addInd('g', 2);
            addInd('T', 3); addInd('t', 3);
        }


        /**
         * A member function that deals with IUPAC ambiguous characters.
         * @param ba the ordered set of all possible symbols which could encode
         * character ch.
         * @param ch the input ambiguous character.
         * @return the number of symbols which encode the ambiguous
         * character ch.
         * @see int Sigma::indC(int *, char)
         */
        inline int indC(int *ba, char ch) {
            char c = toupper(ch);
            switch(c){
            case 'A': case 'a' : ba[0] = 0; return 1;
            case 'C': case 'c' : ba[0] = 1; return 1;
            case 'G': case 'g' : ba[0] = 2; return 1;
            case 'T': case 't' : ba[0] = 3; return 1;
            case 'R': case 'r' : ba[0] = 2; ba[1] = 0; return 2;
            case 'Y': case 'y' : ba[0] = 3; ba[1] = 1; return 2;
            case 'K': case 'k' : ba[0] = 2; ba[1] = 3; return 2;
            case 'M': case 'm' : ba[0] = 0; ba[1] = 1; return 2;
            case 'S': case 's' : ba[0] = 2; ba[1] = 1; return 2;
            case 'W': case 'w' : ba[0] = 0; ba[1] = 3; return 2;
            case 'B': case 'b' : ba[0] = 2; ba[1] = 3; ba[2] = 1; return 3;
            case 'D': case 'd' : ba[0] = 2; ba[1] = 0; ba[2] = 3; return 3;
            case 'H': case 'h' : ba[0] = 0; ba[1] = 1; ba[2] = 3; return 3;
            case 'V': case 'v' : ba[0] = 2; ba[1] = 1; ba[2] = 0; return 3;
            case 'N': case 'n' : ba[0] = 0; ba[1] = 2; ba[2] = 1; ba[3] = 3; return 4;
            case 'X': case 'x' : ba[0] = 0; ba[1] = 2; ba[2] = 1; ba[3] = 3; return 4;
            default: return -1;
            }
            return -1;
        }

        /**
         * @see Resource::saveHeader(std::ofstream &)
         */
        inline void saveHeader(std::ofstream &fd) {
            Resource::saveHeader(fd);
            fd << endl;
        }

        inline int indC(char ch) {
            return ind[(int)ch];
        }

        /**
         * A member function that deals with dna strings containing IUPAC
         * ambiguous characters.
         * @see lenIndS(char *str, unsigned int len)
         */
        int lenIndS(char *str, unsigned int len) {
            int *ba = new int[alphabetSize() + 1];
            int salen = 1;
            for(unsigned int i=0; i<len; i++) {
                int numChars  = indC(ba,str[i]);
                assert(numChars > 0);
                salen *= numChars;
            }
            delete [] ba;
            return salen;
        }

        /**
         * A member function that deals with dna strings containing IUPAC
         * ambiguous characters.
         * @param sa the ordered set of all possible words which could encode
         * str
         * @param str the input ambiguous string
         * @param len the length of str
         * @return the number of all possible words which can encode str.
         * @see int Sigma::indS(int *, char *, unsigned int)
         */
        int indS(int *sa, char *str, unsigned int len) {
            unsigned int i;
            int j,salen = 1, *ba = new int[alphabetSize() + 1];
            assert(len && len <= strlen(str));
            salen = lenIndS(str, len);
            for(j = 0;j < salen;j++) sa[j] = 0;
            for(i = 0; i < len; i++) {
                if(str + i) {
                    int numChars = indC(ba,str[i]);
                    assert(numChars > 0);
                    for(int k=0;k<salen;k++) {
                        sa[k] = sa[k] * alphabetSize() + ba[k % numChars];
                        //      cout << str <<" " << i << "\n";
                    }
                }
            }
            delete [] ba;
            return salen;
        }

        /**
         * @param str the input string
         * @param len the input string's length
         * @return the word which encodes str[0..len]
         * @see int Sigms::indS(char *, unsigned int)
         */
        inline int indS(char *str, unsigned int len) {
            unsigned int i, res=0;
            if(!len) return 0;
            assert(len <= strlen(str));
            for(i = 0; i < len; i++)
                if(str + i) {
                    int baseCode = indC(str[i]);
                    if(baseCode < 0)
                        baseCode = indC(charsPerSymbol[(int)(alphabetSize() * rand()/(RAND_MAX + 1.0))][0]);

                    res = res * alphabetSize() + baseCode;
                }
            return res;
        }

        ~DNASigma() {}

    };


    /**
     * AASigma is a hardcoded Sigma object with all the properties of an
     * alphabet to recognize amino acid sequences(proteins).
     */
    class AASigma : public Sigma {
    protected:
        int dna2aa[64];
        DNASigma s;

    public:
    AASigma() : Sigma(AA_ALPHABET_SIZE, AA_CHARS, COMP_AA_CHARS) {
            initialize();
        }

    AASigma(std::vector<std::string> &params,
            int & offset,
            ResourceEngine *re) :
        Sigma(params,
              offset,
              AA_ALPHABET_SIZE,
              AA_CHARS,
              COMP_AA_CHARS) {

            initialize();

        }

        /**
         * A member function that initializes an alphabet for protein sequences
         */
        void initialize() {
            int i;
            for(i = 0; i < 64; i++)
                dna2aa[i] = 0;

            addInd('*', 0);
            dna2aa[s.indS((char *)"TAA", 3)] = 0; dna2aa[s.indS((char *)"TAG", 3)] = 0; dna2aa[s.indS((char *)"TGA", 3)] = 0;
            addInd('A', 1);
            dna2aa[s.indS((char *)"GCA", 3)] = 1; dna2aa[s.indS((char *)"GCC", 3)] = 1; dna2aa[s.indS((char *)"GCG", 3)] = 1; dna2aa[s.indS((char *)"GCT", 3)] = 1;
            addInd('C', 2);
            dna2aa[s.indS((char *)"TGC", 3)] = 2; dna2aa[s.indS((char *)"TGT", 3)] = 2;
            addInd('D', 3);
            dna2aa[s.indS((char *)"GAC", 3)] = 3; dna2aa[s.indS((char *)"GAT", 3)] = 3;
            addInd('E', 4);
            dna2aa[s.indS((char *)"GAA", 3)] = 4; dna2aa[s.indS((char *)"GAG", 3)] = 4;
            addInd('F', 5);
            dna2aa[s.indS((char *)"TTC", 3)] = 5; dna2aa[s.indS((char *)"TTT", 3)] = 5;
            addInd('G', 6);
            dna2aa[s.indS((char *)"GGA", 3)] = 6; dna2aa[s.indS((char *)"GGC", 3)] = 6; dna2aa[s.indS((char *)"GGG", 3)] = 6; dna2aa[s.indS((char *)"GGT", 3)] = 6;
            addInd('H', 7);
            dna2aa[s.indS((char *)"CAC", 3)] = 7; dna2aa[s.indS((char *)"CAT", 3)] = 7;
            addInd('I', 8);
            dna2aa[s.indS((char *)"ATA", 3)] = 8; dna2aa[s.indS((char *)"ATC", 3)] = 8; dna2aa[s.indS((char *)"ATT", 3)] = 8;
            addInd('K', 9);
            dna2aa[s.indS((char *)"AAA", 3)] = 9; dna2aa[s.indS((char *)"AAG", 3)] = 9;
            addInd('L', 10);
            dna2aa[s.indS((char *)"CTA", 3)] = 10; dna2aa[s.indS((char *)"CTC", 3)] = 10; dna2aa[s.indS((char *)"CTG", 3)] = 10; dna2aa[s.indS((char *)"CTT", 3)] = 10; dna2aa[s.indS((char *)"TTA", 3)] = 10;
            dna2aa[s.indS((char *)"TTG", 3)] = 10;
            addInd('M', 11);
            dna2aa[s.indS((char *)"ATG", 3)] = 11;
            addInd('N', 12);
            dna2aa[s.indS((char *)"AAC", 3)] = 12; dna2aa[s.indS((char *)"AAT", 3)] = 12;
            addInd('P', 13);
            dna2aa[s.indS((char *)"CCA", 3)] = 13; dna2aa[s.indS((char *)"CCC", 3)] = 13; dna2aa[s.indS((char *)"CCG", 3)] = 13; dna2aa[s.indS((char *)"CCT", 3)] = 13;
            addInd('Q', 14);
            dna2aa[s.indS((char *)"CAA", 3)] = 14; dna2aa[s.indS((char *)"CAG", 3)] = 14;
            addInd('R', 15);
            dna2aa[s.indS((char *)"AGA", 3)] = 15; dna2aa[s.indS((char *)"AGG", 3)] = 15; dna2aa[s.indS((char *)"CGA", 3)] = 15; dna2aa[s.indS((char *)"CGC", 3)] = 15; dna2aa[s.indS((char *)"CGG", 3)] = 15; dna2aa[s.indS((char *)"CGT", 3)] = 15;
            addInd('S', 16); dna2aa[s.indS((char *)"AGC", 3)] = 16; dna2aa[s.indS((char *)"AGT", 3)] = 16; dna2aa[s.indS((char *)"TCA", 3)] = 16; dna2aa[s.indS((char *)"TCC", 3)] = 16; dna2aa[s.indS((char *)"TCG", 3)] = 16; dna2aa[s.indS((char *)"TCT", 3)] = 16;
            addInd('T', 17);
            dna2aa[s.indS((char *)"ACA", 3)] = 17; dna2aa[s.indS((char *)"ACC", 3)] = 17; dna2aa[s.indS((char *)"ACG", 3)] = 17; dna2aa[s.indS((char *)"ACT", 3)] = 17;
            addInd('V', 18);
            dna2aa[s.indS((char *)"GTA", 3)] = 18; dna2aa[s.indS((char *)"GTC", 3)] = 18; dna2aa[s.indS((char *)"GTG", 3)] = 18; dna2aa[s.indS((char *)"GTT", 3)] = 18;
            addInd('W', 19);
            dna2aa[s.indS((char *)"TGG", 3)] = 19;
            addInd('Y', 20);
            dna2aa[s.indS((char *)"TAC", 3)] = 20; dna2aa[s.indS((char *)"TAT", 3)] = 20;
        }

        /**
         * @return the amino acid which the codon represents
         */
        inline int AAIndex(int codon) {
            assert(codon >= 0 && codon < 64);
            return dna2aa[codon];
        }

        /**
         * @see Resource::saveHeader(std::ofstream &)
         */
        inline void saveHeader(std::ofstream &fd) {
            Resource::saveHeader(fd);
            fd << endl;
        }

        ~AASigma() {
        }
    };

}

#endif
