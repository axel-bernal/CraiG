#include "ContextIMM.h"
#include "CountUtils.h"
#include "SequenceUtils.h"
#include "FL_Context.h"
#include "InpFile.h"

/****************************************************************************
 * ContextIMM.cpp - part of the lless namespace, a general purpose
 *                  linear semi-markov structure prediction library
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


#include "Utils.h"

namespace lless {
    /**
     * @see Resource::saveContents(::ofstream &)
     */
    void ContextIMM::saveContents(::ofstream &fd) {
        fd << "CIMM Parameters " << numContexts << endl;
        for(int i = 0; i < numContexts; i++) {
            fd << "Context= " << i << endl;
            imm[i]->storeModel(fd);
        }
    }

    /**
     * A member function that loads the cimm model from an ifstream object
     * @param fd the input ifstream
     */
    void ContextIMM::loadModel(::ifstream & fd) {
        int i, context;
        std::string dummy;

        fd >> dummy;
        fd >> dummy;
        fd >> numContexts;

        if(numDesiredContexts < 0)
            numDesiredContexts = numContexts;

        imm = new IMM * [numDesiredContexts];

        for(i = 0; i < numDesiredContexts; i++)
            imm[i] = NULL;

        //  fd.read((char *)&numContexts, sizeof(unsigned int));
        assert(numContexts > 0);

        for(i = 0; i < numContexts; i++) {
            fd >> dummy >> context;
            //    fd.read((char *)&context, sizeof(int));
            imm[i] = new IMM(fd, sigma);
        }

        // if there wasn't enough contexts, fill the rest with last context
        for(; i < numDesiredContexts; i++) {
            imm[i] = imm[numContexts - 1];
        }

    }

    /**
     * A member function for training the IMM models on a set of fasta input
     * sequences read from an ifstream
     * @param aseqStream the input ifstream
     */
    void ContextIMM::preProcess(std::ifstream & aseqStream) {
        int i;
        unsigned int thresholdFreq = 400;
        std::vector<int> ctxLevRanges;
        Utils::stringSplit<int>(ctxLevRanges, ctxString, ",");

        std::string cName("ctx class");
        FL_Context contexts(0, cName, &ctxLevRanges, sldWindow);

        numContexts = ctxLevRanges.size();
        numDesiredContexts = numContexts;

        cerr << "Loading sequences\n";
        InpFile fasta("default", aseqStream, FASTA, sigma, true);
        list<Sequence *> & annotSeqs = (list<Sequence *> &)fasta.sequences();
        // Computing the imm model

        numContexts = numDesiredContexts;
        double **** nbaseMatrix = new double ***[numContexts];
        imm = new IMM *[numContexts];

        for(i = 0; i < numContexts; i++) {
            nbaseMatrix[i] = NULL;
            imm[i] = NULL;
        }

        for(i = 0; i < numContexts; i++)
            if(CountUtils::initkcodonMatrix(&nbaseMatrix[i],
                                            numPos,
                                            minOrd, maxOrd,
                                            0, sigma) < 0)

                throw EXCEPTION(OUT_OF_MEMORY, "Cannot initialize codon Matrix");

        cerr << "Building frequency tables";
        list<Sequence *>::iterator cit = annotSeqs.begin();
        for( ; cit != annotSeqs.end(); cit++) {
            Sequence &c = *(*cit);
            UCHAR *values = new UCHAR [c.length() + 100];
            contexts.computeVals(c.getSeq(STRAND_FWD), values + 50, c.length());

            CountUtils::countNGrams(nbaseMatrix,
                                    sigma,
                                    numPos, minOrd, maxOrd,
                                    c.getSeq(STRAND_FWD),
                                    values + 50, STRAND_FWD);
            delete [] values;
        }

        //computing distances between CTX models and augment counts accordingly
        ChiSquare *chiSqr = new ChiSquare();
        int e;
        double *cTable[2];
        int alphabetSize = sigma->alphabetSize();
        int cTableRows = alphabetSize*numPos;
        int ctxClass4i, ctxClass4j;

        cTable[0] =  new double[cTableRows];
        cTable[1] =  new double[cTableRows];

        for(i = 0; i < numContexts; i++) {
            double *weights = new double[numContexts];
            double totalWeight = 0;
            ctxClass4i = ctxLevRanges[i];
            // fill annotSeqency table for i
            double min = -1;
            //      cerr << ctxClass4i << endl;
            for(e = 0; e < cTableRows; e++) {
                cTable[0][e] = nbaseMatrix[i][0][e % alphabetSize][e / alphabetSize + 1];
                //        cerr << e << " " << cTable[0][e] << endl;
                if(min < 0 || min > cTable[0][e])
                    min = cTable[0][e];
            }
            assert(min);
            for(e = 0; e < cTableRows; e++)
                cTable[0][e] = cTable[0][e]*2/min;

            int j;
            for(j = 0; j < numContexts; j++) {
                ctxClass4j = ctxLevRanges[j];
                // fill annotSeqency table for j

                for(e = 0; e < cTableRows; e++)
                    cTable[1][e] = nbaseMatrix[j][0][e % alphabetSize][e / alphabetSize + 1]*2/min;

                weights[j] = exp(-chiSqr->chisqrValue(cTable, 2, cTableRows));
                //        cerr << "dists " << ctxClass4i << " " << ctxClass4j << " " << weights[j] << " " << chiSqr->chisqrValue(cTable, 2, cTableRows) << endl;
                //        weights[j] = chiSqr->chisqrProbDst(q, cTableRows - 1);
                totalWeight += weights[i];
            }
            for(j = 0; j < numContexts; j++) {
                ctxClass4j = ctxLevRanges[j];
                weights[j] = numContexts*weights[j]/totalWeight;

                for(int k = 0; k < maxOrd - minOrd; k++) {

                    int numFL_Grams = Utils::pow2(alphabetSize, k + minOrd);
                    for(int l = 0; l <= numFL_Grams; l++)
                        for(int n = 0; n <= numPos; n++)
                            ; //nbaseMatrix[i][k][l][n] += weights[j]*nbaseMatrix[j][k][l][n];
                }
            }
            delete [] weights;
        }

        delete [] cTable[0];
        delete [] cTable[1];

        //computing interpolation, using new augmented counts
        cerr << "Computing interpolation";

        for(i = 0; i < numContexts; i++) {
            ctxClass4i = ctxLevRanges[i];

            for(int j=0; j < alphabetSize; j++) // fixing threshold in case genome is too small
                for(int k=1; k <= numPos; k++)
                    if(thresholdFreq > nbaseMatrix[i][0][j][k]) {
                        thresholdFreq = (unsigned int)nbaseMatrix[i][0][j][k];
                    }

            imm[i] = new IMM(nbaseMatrix[i],
                             sigma,
                             maxOrd,
                             numPos,
                             thresholdFreq,
                             0);

            imm[i]->buildCp(); // building lambda values from the codon tables

            CountUtils::freekcodonMatrix(&nbaseMatrix[i],
                                         minOrd, maxOrd,
                                         sigma);

        }
        cerr << "Done!\n";

    }

    ContextIMM::~ContextIMM() {
        for(int i = 0; i < numContexts; i++)
            if(imm[i]) {
                delete imm[i];
                imm[i] = NULL;
            }
        if(imm)
            delete [] imm;
        imm = NULL;
    }

}
