#include "IMM.h"
#include <string.h>
#include <math.h>
#include <time.h>

/****************************************************************************
 * IMM.cpp - part of the lless namespace, a general purpose
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


/**
 * Main Constructor
 * @param alphabet the alphabet for the input sequences which were/are used
 * to train this IMM.
 * @param n is maximum order of the markov chain
 * @param np period, i.e. number of subchains that will model any
 * sequence.
 * @param buildC deprecated
 */
MM::MM(Sigma *alphabet, int n, int np, int buildC) {

    this->alphabet = alphabet;
    buildComplemen = buildC + 1;
    numNGrams = np*buildComplemen;
    ngrams = new NGram* [numNGrams];
    assert(ngrams);

    for(int i = 0; i < numNGrams; i++) {
        ngrams[i] = new NGram(alphabet, n);
        assert(ngrams[i]);
    }
}

/**
 * Constructor that read wrd frequencies from array freqMat
 * @param n order
 * @param np period, i.e. np subchains will be used to model a sequene
 * one for each position within the period.
 */

MM::MM(double  ***freqMat, Sigma *alphabet, int n, int np, int buildC) {
    buildComplemen = buildC + 1;
    numNGrams = np*buildComplemen;
    this->alphabet = alphabet;
    ngrams = new NGram* [numNGrams];
    assert(ngrams);

    /*
     * position 0 in the frequency vector for each gram/word has the total count
     * for the row.
     */
    for(int i = 0; i < numNGrams; i++) {

        ngrams[i] = new NGram(alphabet, n);
        assert(ngrams[i]);

        for(int j = 0; j < ngrams[i]->maxOrder(); j++) {
            FL_Gram<int> *gram = ngrams[i]->gram(j);
            gram->setWordFreq((unsigned int)freqMat[j][gram->maxNumFilterValues()][i + 1]);

            for(int k = 0; k < gram->maxNumFilterValues(); k++)
                gram->setWordFreq(k, (unsigned int)freqMat[j][k][i + 1]);
        }
    }
}

/**
 * Reads word frequencies. Deprecated
 */
void MM::readTable(::ifstream *file) {
    int i = 0, j, dummy;
    char name[20];
    if(!file || !(*file))
        return;
    double freq;
    for(j = 0; j < ngrams[i]->maxOrder(); i++) {
        if(((*file) >> dummy), !dummy) {
            for(i = 0;i <= numNGrams;i++) {
                //      cerr << numPositions << "\n";
                (*file) >> freq;
                ngrams[i]->gram(j)->setWordFreq((unsigned int)freq);
            }
        }

        while(((*file) >> name), (name && strcmp(name,"//"))) {
            int indFreq = alphabet->indS(name, strlen(name));

            for(i = 0; i < numNGrams; i++) {
                assert(strlen(name) == ngrams[i]->gram(j)->order());
                (*file) >> freq;
                ngrams[i]->gram(j)->setWordFreq(indFreq, (unsigned int)freq);
            }
        }
        cerr << ".";
    }
}

/**
 * @return probability of seq, using order-degree interpolated markov model
 * and subchain numNGram in the first position of seq. If any ambiguous
 * character is found in seq, the model averages the probabilities of all
 * possible sequences that are generated under differente hypothesis, each
 * of them correpsponding to a differente character value. Only works well with
 * alphabets with small number of symbols. Too expensive for other alphabets.
 */
double MM::getavgCp(char *seq, int order, int numNGram) {
    double res = 0;
    DNASigma *dnaAlphabet = (DNASigma *)alphabet;
    int *ba = new int[dnaAlphabet->alphabetSize() + 1];
    int numChars =  dnaAlphabet->indC(ba,seq[order]);
    int numStrings =  dnaAlphabet->lenIndS(seq,order);
    assert(numChars > 0);
    int *bb = new int[numStrings];
    dnaAlphabet->indS(bb,seq,order);
    //    cerr << ptr <<  " " << bb[0] << " " << order;

    for(int i=0;i<numChars;i++) {
        for(int j=0;j<numStrings;j++)
            res += ngrams[numNGram]->gram(order)->getCp(ba[i],bb[j]);
    }

    delete [] bb;
    delete [] ba;

    return res/(numChars*numStrings);
}

/**
 * @return probability of seq, using order-degree interpolated markov model
 * and subchain numNGram in the first position of seq. If ambiguous
 * characters are found in the sequence, the model generates a random
 * character instead using information about the relative frequencies of
 * symbols present in the sequence.
 */
double MM::getambCp(char *seq, int order, int numNGram) {
    double res = 0;
    int baseCode, stringCode;
    assert(ngrams[numNGram]->gram(0));
    baseCode = alphabet->indC(seq[order]);

    if(baseCode < 0)
        baseCode = alphabet->indC(alphabet->contC((int)(alphabet->alphabetSize()*rand()/(RAND_MAX + 1.0))));

    stringCode = alphabet->randIndS(seq, order);
    res = ngrams[numNGram]->gram(order)->getCp(baseCode,stringCode);
    return res;
}

MM::~MM() {
    if(!ngrams)
        return;
    for(int i = 0; i < numNGrams; i++) {
        if(ngrams[i])
            delete ngrams[i];
    }
    delete [] ngrams;
}


/*
 * Main Constructor
 * @param fdes file name of the table containing relative frequencies of words
 * @param alphabet the alphabet for the input sequences which were/are used
 * to train this IMM.
 * @param n is maximum order of the markov chain
 * @param np number of subchains
 * @param buildComp deprecated
 */
FMM::FMM(char *fdes,
         Sigma *alphabet,
         int n,
         int np,
         int buildComp) : MM(alphabet,
                             n,
                             np,
                             buildComp) {

    if(!fdes)
        throw string("must provide file containing table in FMM\n");

    ::ifstream streamDes(fdes);
    assert(streamDes);
    readTable(&streamDes);
}


/*
 * Main Constructor
 * @param fdes file name of the table containing relative frequencies of words
 * @param alphabet the alphabet for the input sequences which were/are used
 * to train this IMM.
 * @param n is maximum order of the markov chain
 * @param np number of subchains
 * @param maxf threshold to indicate the minimum number of counts of each word.
 * If the word w, has frequency that is too low to trust then apply chi-square
 * test between distributions of words of same length as w and words with
 * shorter length (given by word prefixes) If the test is positive
 * use the prefix word instead as it should have higher frequency, if the test
 * fails the model cannot estimate reliable the probability of w and it's
 * assigned a value of zero.
 * @param buildComp deprecated
 */

IMM::IMM(char *fdes,
         Sigma *alphabet,
         int n,
         int np,
         unsigned int maxf,
         int buildComp) : MM(alphabet,
                             n,
                             np,
                             buildComp) {

    if(!fdes)
        throw string("must provide file containing table in IMM\n");

    ::ifstream streamDes(fdes);
    assert(streamDes);
    readTable(&streamDes);
    thresholdFreq = maxf;
    assert(n -1 <= MAX_ORDER);
    lambda = new double *** [numNGrams];
    assert(lambda);
    interpolatedProbs = NULL;

    for(int i = 0; i < numNGrams; i++) {
        lambda[i] = new double ** [ngrams[i]->maxOrder()];
        assert(lambda[i]);

        for(int j = 0; j < ngrams[i]->maxOrder(); j++) {
            FL_Gram<int> *gram = ngrams[i]->gram(j);
            assert(gram);
            int limit = gram->maxNumFilterValues();
            lambda[i][j] = new double * [gram->alphabetSize()];
            assert(lambda[i][j]);

            for(int k = 0; k < gram->alphabetSize(); k++) {
                lambda[i][j][k] = new double [limit];
                assert(lambda[i][j][k]);

                for(int l = 0; l < limit; l++)
                    lambda[i][j][k][l] = 0;
            }
        }
    }
    // Load chi-square table
    chisqr = new ChiSquare();
    assert(chisqr);
}

/*
 * Loading constructor
 * @param params the input stream containing the IMM parameters
 * @param the alphabet
 */
IMM::IMM(::ifstream & params, Sigma *alphabet) {
    int i,j;
    this->alphabet = alphabet;

    std::string dummy;
    params >> buildComplemen;
    params >> dummy >> numNGrams;
    ngrams = new NGram* [numNGrams];
    assert(ngrams);
    int alphabetSize = alphabet->alphabetSize();
    interpolatedProbs = new double ** [numNGrams];
    lambda = NULL;

    /*
     * Initializing arrays
     */
    for(i = 0; i < numNGrams; i++) {
        int maxOrder;
        params >> dummy;
        params >> maxOrder;
        ngrams[i] = new NGram(maxOrder);
        assert(ngrams[i]);
        interpolatedProbs[i] = new double * [maxOrder];
        assert(interpolatedProbs[i]);

        for(j = 0; j < maxOrder; j++) {
            params >> dummy;
            int numWords = (int)(Utils::pow2(alphabetSize, j+1));
            interpolatedProbs[i][j] = new double [numWords];
            assert(interpolatedProbs[i][j]);

            for(int k = 0; k < numWords; k++)
                params >> interpolatedProbs[i][j][k];

        }
    }

    params >> dummy;
    assert(dummy.compare("//") == 0);

    /*
      params.read((char *)&buildComplemen,sizeof(int));
      params.read((char *)&numNGrams,sizeof(int));
      ngrams = new NGram* [numNGrams];
      assert(ngrams);
      int alphabetSize = alphabet->alphabetSize();
      interpolatedProbs = new double ** [numNGrams];
      lambda = NULL;
      // initializing arrays
      for(i = 0; i < numNGrams; i++) {
      int maxOrder;
      params.read((char *)&maxOrder, sizeof(int));
      ngrams[i] = new NGram(maxOrder);
      assert(ngrams[i]);
      interpolatedProbs[i] = new double * [maxOrder];
      assert(interpolatedProbs[i]);
      for(j = 0; j < maxOrder; j++) {
      int numWords = (int)(Utils::pow2(alphabetSize, j+1));
      interpolatedProbs[i][j] = new double [numWords];
      assert(interpolatedProbs[i][j]);
      for(int k = 0; k < numWords; k++)
      params.read((char *)&interpolatedProbs[i][j][k], sizeof(double));
      }
      }
    */
    chisqr = NULL;
}

/*
 * Constructor from an array of word frequencies freqMat
 */
IMM::IMM(double ***freqMat,
         Sigma *alphabet,
         int n,
         int np,
         unsigned int maxf,
         int buildComp) : MM(freqMat,
                             alphabet,
                             n,
                             np,
                             buildComp) {

    thresholdFreq = maxf;
    assert(n - 1 <= MAX_ORDER);
    lambda = new double *** [numNGrams];
    assert(lambda);
    interpolatedProbs = NULL;
    int alphabetSize = alphabet->alphabetSize();

    for(int i = 0; i < numNGrams; i++) {
        lambda[i] = new double ** [ngrams[i]->maxOrder()];
        assert(lambda[i]);

        for(int j = 0;j < ngrams[i]->maxOrder(); j++) {
            FL_Gram<int> *gram = ngrams[i]->gram(j);
            assert(gram);
            lambda[i][j] = new double *[alphabetSize];
            assert(lambda[i][j]);

            for(int k = 0; k < alphabetSize; k++) {
                lambda[i][j][k] = new double [gram->numHistories()];
                assert(lambda[i][j][k]);

                for(int l = 0; l < gram->numHistories(); l++)
                    lambda[i][j][k][l] = 0;
            }
        }
    }

    // Load chi-square table
    chisqr = new ChiSquare();
    assert(chisqr);
}

/**
 * A member function for reporting probabilities of each gram
 */

void IMM::probeAllNGramProbs() {
    assert(lambda);
    int alphabetSize = 4;

    for(int numNGram = 0; numNGram < numNGrams; numNGram++) {
        char *s = new char[ngrams[numNGram]->maxOrder() + 1];

        for(int order = 0; order < ngrams[numNGram]->maxOrder(); order++) {
            FL_Gram<int> *gram = ngrams[numNGram]->gram(order);
            Sigma *sigma = gram->alphabet();
            int limit = (int)pow((double)alphabetSize, order + 1);

            for(int i = 0; i < limit; i++) {
                sigma->contS(s, i, order + 1);
                double res = imm(order, numNGram, s, order + 1, 1);
                cerr << numNGram << " " << order << " " << s  <<  " " << res << " " << lambda[numNGram][order][i % 4][i / 4] << " " << gram->getCp(i%4, i/4) << " (" << i/4 << " " << i%4 << ")" <<endl;
            }
        }

        delete [] s;
    }
}

/**
 * Emits(generates) a sequence seq of length len, using subchain frame to
 * generate the first character of seq.
 */
void IMM::emitSequence(char *seq, int frame, int len)  {
    assert(interpolatedProbs);
    int word = 0, i, suffix;
    char *ptr;
    int limit, numNGram;
    assert(frame >= 1 && frame <= numNGrams);
    int alphabetSize = alphabet->alphabetSize();
    double *suffixProbs = new double[alphabetSize];

    for(int j = 0; j < len; j++)  {
        numNGram = (j + frame - 1) % numNGrams + 3*(frame > 3);
        limit = Utils::min(j + 1, ngrams[numNGram]->maxOrder()) - 1;
        ptr = seq + j - limit - 1;

        if(j - limit > 0)
            word = alphabet->removePrefix(word, limit + 1);

        suffixProbs[0] = interpolatedProbs[numNGram][limit][alphabet->addSuffix(word, 0)];

        for(i = 1; i < alphabetSize; i++)
            suffixProbs[i] = suffixProbs[i - 1] + interpolatedProbs[numNGram][limit][alphabet->addSuffix(word, i)];

        suffix = Utils::pickValueAtRandom(alphabetSize, suffixProbs);
        assert(suffix >= 0 && suffix < alphabetSize);
        ptr[limit + 1] = alphabet->contC(suffix);
        word = alphabet->addSuffix(word, suffix);
    }

    seq[len] = '\0';
    delete [] suffixProbs;
}

/**
 * Compute interpolated probabilities (the lambda values)
 */
void IMM :: buildCp(void) {
    assert(lambda);

    for(int i = 0; i < numNGrams; i++)
        ngrams[i]->buildCp();

    // Calculating lambda's
    interpolatedProbs = new double **[numNGrams];
    int counter = 0;

    for(int k = 0; k < numNGrams; k++) {
        interpolatedProbs[k] = new double * [ngrams[k]->maxOrder()];

        for(int i = 0; i < ngrams[k]->maxOrder(); i++) {
            PRINT_TICKS(".", counter++, ngrams[k]->maxOrder()*numNGrams)
                FL_Gram<int> *gram = ngrams[k]->gram(i);
            Sigma *sigma = gram->alphabet();
            int limit = gram->maxNumFilterValues();
            interpolatedProbs[k][i] = new double [limit];

            for(int j = 0; j < limit; j++) {
                unsigned int freq;
                int l, m, b;
                l = j % ngrams[k]->alphabetSize(); // r+1 state
                m = j / ngrams[k]->alphabetSize(); // r states

                if(i - 1 < 0)
                    freq = gram->getWordFreq(sigma->removePrefix(j, i + 1));
                else
                    freq = (ngrams[k]->gram(i - 1))->getWordFreq(sigma->removePrefix(j, i + 1));
                double prob;
                prob = gram->getCp(l,m);

                if(freq >= thresholdFreq)
                    lambda[k][i][l][m] = 1.0;
                else { // oops, frequency is too low to trust in it apply chi-square
                    char *context = new char[gram->order() + 1];
                    double ** contingencyTable = new double * [ngrams[k]->alphabetSize()];
                    double q;
                    assert(contingencyTable);
                    for(b=0;b < ngrams[k]->alphabetSize();b++) { // get the frequencies
                        // of the following base
                        contingencyTable[b] = new double [2];
                        assert(contingencyTable[b]);
                        for(int e=0;e < 2;e++) {
                            if(e)
                                q = gram->getWordFreq(m*ngrams[k]->alphabetSize() + b);
                            else {
                                assert(i-2 >= 0);
                                sigma->contS(context, sigma->removePrefix(m, i)*ngrams[k]->alphabetSize() + b, gram->order() - 1);
                                q = imm(i-1, k, context, gram->order() - 1, 1)*ngrams[k]->gram(i - 2)->getWordFreq(m / ngrams[k]->alphabetSize());
                                /* adjust IMM for countings instead of probabilities,
                                 * otherwise chi-square test does not apply.. revise
                                 * this in Salzberg paper.
                                 */
                            }
                            contingencyTable[b][e] = q?q:MIN_FREQ;
                        }
                    }

                    q = chisqr->chisqrValue(contingencyTable,
                                            ngrams[k]->alphabetSize(), 2);

                    if((q=chisqr->chisqrProbDst(q, (ngrams[k]->alphabetSize() - 1)*(2 - 1))) >= 0.5) {
                        lambda[k][i][l][m] = (q/thresholdFreq)*freq;
                    }
                    else
                        lambda[k][i][l][m] = 0;

                    for(b=0;b < ngrams[k]->alphabetSize();b++)
                        delete [] contingencyTable[b];

                    delete [] contingencyTable;
                    delete [] context;
                }
                assert(lambda[k][i][l][m] <= 1.0 && lambda >= 0);
            }
        }
    }
    for(int k = 0; k < numNGrams; k++)
        computeInterpolatedProbs(k);
}


/*
//  This is the new buildCp, which may not be correct. Takes the prefix out (instead of the suffix)
void IMM :: buildCp(void) {
cerr << "Building codon table CP's ";
for(int i = 0; i < numNGrams; i++) {
ngrams[i]->buildCp();
cerr << ".";
}
cerr << "Done!\n";
// Calculating lambda's
cerr << "Calculating Lambda parameters";
for(int k = 0; k < numNGrams; k++) {
for(int i = 0; i < ngrams[k]->maxOrder(); i++) {
cerr << ".";
FL_Gram *gram = ngrams[k]->gram(i);
Sigma *sigma = gram->alphabet();
int limit = gram->maxNumFilterValues();
for(int j = 0; j < limit; j++) {
unsigned int freq;
int l, m, b;
l = sigma->wordTail(j); // r+1 state
m = sigma->wordHead(j); // r states
if(i - 1 < 0)
freq = gram->getWordFreq(m);
else
freq = (ngrams[k]->gram(i - 1))->getWordFreq(m);
double prob;
prob = gram->getCp(l,m);
if(freq >= thresholdFreq)
lambda[k][i][l][m] = 1.0;
else { // oops, frequency is too low to trust in it apply chi-square
char *context = new char[gram->order() + 1];
double ** contingencyTable = new double * [ngrams[k]->alphabetSize()];
double q;
assert(contingencyTable);
assert(i - 2 >= 0);
int shorter = sigma->removePrefix(m, gram->order() - 1);
for(b = 0; b < ngrams[k]->alphabetSize(); b++) { // get the frequencies of the following base
contingencyTable[b] = new double [2];
assert(contingencyTable[b]);
for(int e=0;e < 2;e++) {
if(e)
q = gram->getWordFreq(sigma->addSuffix(m, b));
else {
gram->alphabet()->contS(context, sigma->addSuffix(shorter, b), gram->order() - 1);
q = imm(i - 1, k, context, gram->order() - 1, 1)*ngrams[k]->gram(i - 2)->getWordFreq(shorter);
// adjust IMM for countings instead of probabilities, otherwise
// chi-square test does not apply.. revise this in Salzberg pap
//   er.

}
contingencyTable[b][e] = q?q:MIN_FREQ;
}
}
q = chisqr->chisqrValue(contingencyTable, ngrams[k]->alphabetSize(), 2);
if((q=chisqr->chisqrProbDst(q, (ngrams[k]->alphabetSize() - 1)*(2 - 1))) >= 0.5)
lambda[k][i][l][m] = (q/thresholdFreq)*freq;
else
lambda[k][i][l][m] = 0;
for(b = 0; b < ngrams[k]->alphabetSize(); b++)
delete [] contingencyTable[b];
delete [] contingencyTable;
delete [] context;
}
assert(lambda[k][i][l][m] <= 1.0 && lambda >= 0);
}
}
}
cerr << "Done!\n";
}
*/

/**
 * Gets rid of storage space and make a much more efficient scoring pass over
 * a sequence, by using the lambda parameters directly to compute the final
 * probability for each gram that appears in the sequence.
 */
void IMM::computeInterpolatedProbs(int numNGram) {
    assert(lambda);
    int baseCode, stringCode;
    int alphabetSize = ngrams[numNGram]->alphabetSize();
    double l;
    FL_Gram<int> *gram = ngrams[numNGram]->gram(0);
    /*
     * Initializing the first simple letters first
     */
    Sigma *sigma;

    for(int order = 0; order < ngrams[numNGram]->maxOrder(); order++) {
        gram = ngrams[numNGram]->gram(order);
        sigma = gram->alphabet();
        int limit = (int)pow((double)alphabetSize, order + 1);

        for(int i = 0; i < limit; i++) {
            baseCode = sigma->wordTail(i);
            stringCode = sigma->wordHead(i);
            l = lambda[numNGram][order][baseCode][stringCode];
            interpolatedProbs[numNGram][order][i] = l*gram->getCp(baseCode, stringCode);

            if(order)
                interpolatedProbs[numNGram][order][i] += (1 - l)*
                    interpolatedProbs[numNGram][order - 1][sigma->removePrefix(i, order + 1)];

            if(interpolatedProbs[numNGram][order][i] == 0)
                interpolatedProbs[numNGram][order][i] = MIN_PROB;

            /*      char s[order + 1];
                    sigma->contS(s, i, order + 1);
                    cerr << numNGram << " " << order << " " << s << " " << interpolatedProbs[numNGram][order][i] << " " << l << " " << gram->getCp(baseCode, stringCode) << " (" << stringCode << " " << baseCode <<")" <<endl; */
        }
    }
}


/**
 * This is a slower version for computing the probability of seq. This
 * function uses lambdas and table frequencies for each NGram object and
 * it is deprecated because of IMM::computeInterpolatedProbs(int).
 */
double IMM::imm(unsigned int order,
                int numNGram,
                char *seq,
                unsigned int x,
                int amb)  {

    assert(lambda);
    double res = 0, l;
    char * ptr;
    FL_Gram<int> *gram;
    Sigma *sigma;
    assert(amb);
    double factor = 1;
    int baseCode, stringCode;
    int limit = Utils::min(x - 1, order);

    for(int ord = limit; ord >= 0; ord--) {
        gram = ngrams[numNGram]->gram(ord);
        sigma = gram->alphabet();
        ptr = seq + x - (ord + 1);
        baseCode = sigma->indC(ptr[ord]);

        if(baseCode < 0) {
            int randIndex  = (int)(sigma->alphabetSize()*rand()/(RAND_MAX + 1.0));
            baseCode = sigma->indC(sigma->contC(randIndex));
        }

        stringCode = sigma->randIndS(ptr,ord);
        l = lambda[numNGram][ord][baseCode][stringCode];
        res += factor*l*gram->getCp(baseCode,stringCode);

        if(l == 1) return res;
        factor *= (1-l);
    }

    if(res < 0)  {
        cerr << "Error computing imm value\n";
        assert(0);
    }

    if(!res)
        res = 0.25;

    return res;
}

/**
 * Scores sequence seq using subchain frame in the first position of seq.
 * @return the probability of sequence seq as a whole.
 */

double IMM::seqScore(char *seq,
                     int amb,
                     int frame,
                     int len)  {

    assert(amb && interpolatedProbs);
    int word = 0, suffix;
    char * ptr;
    int limit, numNGram;
    double res = 0;
    assert(frame >= 1 && frame <= numNGrams);

    for(int j = 0; j < len; j++)  {
        numNGram = (j + frame - 1) % numNGrams + 3*(frame > 3);
        limit = Utils::min(j + 1, ngrams[numNGram]->maxOrder()) - 1;
        ptr = seq + j - limit - 1;

        if(j - limit > 0)
            word = alphabet->removePrefix(word, limit + 1);

        suffix = alphabet->indC(ptr[limit + 1]);

        if(suffix < 0) {
            int randInd = (int)(alphabet->alphabetSize()*rand()/(RAND_MAX + 1.0));
            suffix = alphabet->indC(alphabet->contC(randInd));
        }
        word = alphabet->addSuffix(word, suffix);
        res += log(interpolatedProbs[numNGram][limit][word]);
    }
    return res;
}

/**
 * Scores sequence seq using subchain frame in the first position of seq.
 * Stores the results for each position in array res
 */
void IMM::seqScore(double *res,
                   char *seq,
                   int amb,
                   int frame,
                   int len)  {

    assert(amb && interpolatedProbs);
    int word = 0, suffix;
    char * ptr;
    int limit, numNGram;
    assert(frame >= 1 && frame <= numNGrams);

    for(int j = 0; j < len; j++)  {
        numNGram = (j + frame - 1) % numNGrams + 3*(frame > 3);
        limit = Utils::min(j + 1, ngrams[numNGram]->maxOrder()) - 1;
        ptr = seq + j - limit - 1;

        if(j - limit > 0)
            word = alphabet->removePrefix(word, limit + 1);

        suffix = alphabet->indC(ptr[limit + 1]);

        if(suffix < 0) {
            int randInd = (int)(alphabet->alphabetSize()*rand()/(RAND_MAX + 1.0));
            suffix = alphabet->indC(alphabet->contC(randInd));
        }

        word = alphabet->addSuffix(word, suffix);

        // comment out when modifying probabilities in training
        if(!interpolatedProbs[numNGram][limit][word])
            res[j] = -10;
        else
            res[j] = log(interpolatedProbs[numNGram][limit][word]);
    }
}

/**
 * Stores the gram probabilities
 */
void IMM::storeModel(::ofstream & params)  {
    assert(interpolatedProbs);

    params << buildComplemen << endl;
    params << "num-periods= " << numNGrams << endl;
    //  cerr << buildComplemen << " " << numNGrams << "\n";

    for(int i = 0; i < numNGrams; i++)  {
        int maxOrder = ngrams[i]->maxOrder();
        params << "period=" << i << endl;
        params << maxOrder << endl;

        for(int j = 0; j < maxOrder; j++)  {
            params << "gram-order=" << j << endl;
            int numWords = (int)(Utils::pow2(alphabet->alphabetSize(), j+1));

            for(int k = 0; k < numWords; k++) {
                params << interpolatedProbs[i][j][k];
                if(k < numWords - 1)
                    params << " ";
            }
            params << endl;
        }
    }
    params << "//" << endl;
}

void IMM::deleteLambdas() {
    if(!lambda)
        return;

    for(int h = 0; h < numNGrams; h++)  {
        if(!lambda[h])
            continue;

        for(int i = 0; i < ngrams[h]->maxOrder(); i++)  {
            if(!lambda[h][i])
                continue;

            for(int j = 0; j < alphabet->alphabetSize(); j++) {
                if(!lambda[h][i][j])
                    continue;

                delete [] lambda[h][i][j];
                lambda[h][i][j] = NULL;
            }

            delete [] (lambda[h][i]);
            lambda[h][i] = NULL;
            ngrams[h]->removeFL_Gram(i+1);
        }

        delete [] lambda[h];
        lambda[h] = NULL;
    }

    delete [] lambda;
    lambda = NULL;
}

void IMM::deleteInterpolatedProbs() {
    if(!interpolatedProbs)
        return;

    for(int h = 0; h < numNGrams; h++)  {
        if(!interpolatedProbs[h])
            continue;

        for(int i = 0; i < ngrams[h]->maxOrder(); i++) {
            if(!interpolatedProbs[h][i])
                continue;
            delete [] interpolatedProbs[h][i];
            interpolatedProbs[h][i] = NULL;
        }

        delete [] interpolatedProbs[h];
        interpolatedProbs[h] = NULL;
    }

    delete [] interpolatedProbs;
    interpolatedProbs = NULL;
}

IMM::~IMM() {
    deleteLambdas();
    deleteInterpolatedProbs();
    if(chisqr)
        delete chisqr;
}
