#include <time.h>
#include <stdio.h>
#include "Organism.h"
#include "Gene.h"
#include <sys/stat.h>
#include <iostream>
#include <string>
#include <stdexcept>
#include "CountUtils.h"

/****************************************************************************
 * Organism.cpp - part of the craig namespace, a genomics library
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


using namespace std;

namespace craig {

    Organism::Organism(Sigma &sigma) {
        this->sigma = &sigma;
        setDefaults();
    }

    void Organism::setDefaults() {
        countgenomicBases.reserve(4);
        for(int i = 0; i < sigma->alphabetSize(); i++) {
            weightRBS[i] = NULL;
            countgenomicBases[i] = 0;
        }
        vectorRBS = NULL;
    }

    void Organism::getcodonTables(double ***im,
                                  int numPos,
                                  int minPlets,
                                  int maxPlets,
                                  int bC,
                                  TGene cond,
                                  TSetType geneSet) {

        list<Gene> *listO = &genes[geneSet];
        int counter = 0,initPos;

        for(list<Gene>::iterator it=listO->begin();it != listO->end();it++) {
            if(!(it->state() & cond))
                continue;

            for(unsigned int i = 0; i < it->transcripts().size(); i++) {
                Transcript &t = it->transcripts()[i];
                char * tmpSeq = (char *)t.codingSeq();
                initPos = 1;

                if(bC && t.getStrand() == STRAND_COMP) // complementary frame
                    initPos = 4;

                PRINT_TICKS(".", counter++, listO->size())
                    CountUtils::fillnbaseTable(im, numPos, minPlets, maxPlets, tmpSeq, initPos, sigma);
            }
        }
        cerr << endl << "# of Genes used for training (after filtering): " << counter << endl;
    }

    void Organism::genesInOneStrand(TSetType geneSet) {
        for(list<Gene>::iterator it = genes[geneSet].begin(); it != genes[geneSet].end(); it++)
            it->toOneStrand();
    }


    void Organism::attachGenes2Sequences( TSetType geneSet) {
        list<Gene> *listO = &genes[geneSet];
        listO->sort();
        list<Gene>::iterator it = listO->begin();
        list<Sequence *>::iterator cit = annotSeqs.begin();

        for(; cit != annotSeqs.end(); cit++) {
            Sequence &c = *(*cit);

            while(it != listO->end() && !c.id().compare(it->getSequenceId())) {
                it->setSequence(&c);
                it++;
            }
        }

        if(it != listO->end() || cit != annotSeqs.end())
            cerr << "BAD:" << it->getSequenceId() << " " << (*cit)->id() << "\n";

        assert(it == listO->end() && cit == annotSeqs.end());
    }


    void Organism::filterGenesByLength(int minGeneLen,  TSetType geneSet) {
        list<Gene> *listO = &genes[geneSet];
        list<Gene>::iterator it = listO->begin();

        while(it != listO->end()) {
            //      cerr << geneLen << " "<< it->getId() << " " << it->getBegin() << " " << it->getEnd() << endl << seq << endl;
            for(vector<Transcript>::iterator tit = it->transcripts().begin(); tit != it->transcripts().end();) {
                Transcript &t = (Transcript &)(*tit);

                if(t.transcriptLen() + 1 < (unsigned int)minGeneLen) {
                    tit = it->transcripts().erase(tit);
                    continue;
                }
                tit++;
            }
            if(it->transcripts().size() == 0) {
                it = listO->erase(it);
                continue;
            }
            it++;
        }
    }

    void Organism::getannotSeqGaps(int mingapLen, TSetType geneSet, TSetType gapSet) {
        list<Gene> *listO = &genes[geneSet];
        list<Gene> *listT = &genes[gapSet];
        assert(!listT->size());
        listO->sort();
        list<Gene>::iterator it = listO->begin();
        list<Sequence *>::iterator cit = annotSeqs.begin();

        for( ; cit != annotSeqs.end(); cit++) {
            Sequence &c = *(*cit);
            int i = 0;
            int initGap = 1, endGap;
            ostringstream geneId;

            while(it != listO->end() && !c.id().compare(it->getSequenceId())) {
                endGap = it->begin() - 1;
                geneId.str("");
                geneId <<  c.id() << i;

                if(endGap - initGap + 1 > mingapLen) {
                    std::string geneId2 = geneId.str();
                    listT->push_back(Gene(geneId2, &c, true, STRAND_FWD));
                    Gene &gene = listT->back();
                    Transcript tmp(geneId2, &gene, false, false, 0);
                    Transcript &t = gene.addTranscript(tmp);
                    t.addExon(initGap, endGap, 0);
                }
                initGap = it->end() + 1;
                it++;
                i++;
            }
            // process up to the end of the annotSeq
            endGap = c.length();
            geneId.str("");
            geneId <<  c.id() << i;

            if(endGap - initGap + 1 > mingapLen) {
                std::string geneId2 = geneId.str();
                listT->push_back(Gene(geneId2, &c, true, STRAND_FWD));
                Gene &gene = listT->back();
                Transcript tmp(geneId2, &gene, false, false, 0);
                Transcript &t = gene.addTranscript(tmp);
                t.addExon(initGap, endGap, 0);
            }
        }

        if(it != listO->end()) {
            cerr << "BAD " << it->getSequenceId() << "_" << it->begin() << "_" << it->end() << "\n";
        }

        if(cit != annotSeqs.end()) {
            cerr << (*cit)->id() << " \n";
        }
        assert(it == listO->end() && cit == annotSeqs.end());
    }

    /*
      void Organism::processSIMS(::ifstream *simStream,  TSetType geneSet) {// list must be ordered
      list<Gene> *listO = &genes[geneSet];
      int counter = 0;
      assert(*simStream);
      char line[MAX_LINE], oldGene[100] = "", gene1[100], gene2[100],d;
      int len1, len2, perc, regini1, regini2, regend1, regend2;
      double score;
      ::ifstream *fd = simStream;
      listO->sort();
      fd->getline(line, MAX_LINE);
      //  cout << line << "\n";

      while((!fd->eof())) {
      sscanf(line, "%[^,] %[,] %d %[,] %[^,] %[,] %d %[,] %d %[,] %lf %[,] %d %[,] %d %[,] %d %[,] %d", gene1, &d, &len1, &d, gene2, &d, &len2, &d, &perc, &d, &score, &d, &regini1, &d, &regend1, &d, &regini2, &d, &regend2);
      assert(regend1 && regend2);
      strcpy(oldGene, gene1);
      int i = 0;

      do {
      Sim s(gene1, gene2, len1, len2, perc, score, regini1, regend1, regini2, regend2);

      listSims.push_back(s);
      fd->getline(line, MAX_LINE);
      //	cout << line << "\n";
      i++;

      if(fd->eof())
      break;

      sscanf(line, "%[^,] %[,] %d %[,] %[^,] %[,] %d %[,] %d %[,] %lf %[,] %d %[,] %d %[,] %d %[,] %d", gene1, &d, &len1, &d, gene2, &d, &len2, &d, &perc, &d, &score, &d, &regini1, &d, &regend1, &d, &regini2, &d, &regend2);
      } while(!strcmp(gene1,oldGene));

      //	cout << "pushed\n";
      Gene * gene = Gene::findGene(oldGene, genes[geneSet]);
      assert(gene);

      for(int sims = 1; sims <= i; sims++) {
      assert(listSims.size() - sims >= 0);
      Sim & s = listSims[listSims.size() - sims];
      gene->addSimFact(&s);
      //      cout << oldGene << " " << gene->getBegin() << "_" << gene->getEnd() << s.getPercentage() << " " << s.getScore() << " " << s.getRegini2() << " " << s.getRegend2() << "\n";
      }

      if(fd->eof())
      break;

      PRINT_TICKS(".", counter++, listO->size())
      }
      cerr << "\n" << listSims.size() << " Similarities processed.....";
      }
    */
    void Organism::buildRBSVector(int numGenes, int *offsetRBS) {
        int *freqDistance = new int[windowRBS], max = 0, i;

        for(i = 0 ; i < windowRBS; i++ )
            freqDistance[i] = 0;

        for(i=0;i < numGenes; i++) {

            if(offsetRBS[i] < 0)
                continue;

            freqDistance[offsetRBS[i]]++;
        }

        for(i = 0 ; i < windowRBS ; i++)
            if(freqDistance[i] > freqDistance[max])
                max = i;

        for(i=0;i < windowRBS; i++)
            vectorRBS[i] = log((freqDistance[i] + 0.5)/(freqDistance[max] + 0.5))/log(2.0);

        delete [] freqDistance;
    }


    void Organism::buildRBSweightMatrix(char **upsRegs, int *offsetRBS, int numGenes, double **weightFreq)  {
        char * sdBox = new char[lenRBS + 1];
        int *maxBase = new int[lenRBS], i;

        for(i=0 ; i < numGenes ; i++) {
            if(offsetRBS[i] < 0)
                continue;

            strncpy(sdBox, upsRegs[i] + offsetRBS[i],lenRBS);
            sdBox[lenRBS] = '\0';
            CountUtils::fillbaseTable(weightFreq,lenRBS,1,sdBox, 1, sigma);
        }

        for(i = 0 ; i < lenRBS ; i++) {
            maxBase[i] = 1;

            for(int j = 1; j <= sigma->alphabetSize(); j++)
                if(weightFreq[j][i] > weightFreq[maxBase[i]][i])
                    maxBase[i] = j;
        }

        for(i=0;i < sigma->alphabetSize();i++) {
            for(int j = 0 ; j < lenRBS; j++) {
                weightRBS[i][j] = log((weightFreq[i][j] + 0.5)/(weightFreq[maxBase[i]][j] + 0.5))/log(2.0);
            }
        }

        delete [] sdBox;
        delete  [] maxBase;
    }

    void Organism::refineRBSweightMatrix(char **upsRegs, int *offsetRBS, int numGenes, double **weightFreq, int percentage)  {
        double *scoresRBS = new double[numGenes];
        int firstStageDone = 0;
        int numCurrentgenes;
        percentage = percentage;

        for(;;) {
            int i;
            numCurrentgenes = 0;
            int offsetsModified = 0;

            for(i=0;i < numGenes;i++) {
                if(offsetRBS[i] < 0)
                    continue;

                int maxSD = offsetRBS[i];
                // Using the weight matrix to calculate the offsets
                scoresRBS[i] = Utils::scoreRBSRegion(upsRegs[i], lenRBS, weightRBS, vectorRBS, countgenomicBases, &offsetRBS[i], sigma);
                //      cout << upsRegs[i] + offsetRBS[i] << "\t" << offsetRBS[i] << "\t" << scoresRBS[i] << "\n";

                if(maxSD != offsetRBS[i])
                    offsetsModified = 1;

                numCurrentgenes++;
            }

            if(!offsetsModified) {
                if(firstStageDone)
                    break;
                else
                    firstStageDone = 1;
            }
            //    cout << "Iteration : number of genes = " << numCurrentgenes << "\n";
            for(i=0;i < numGenes - 1;i++) {
                for(int j = i + 1; j < numGenes; j++) {
                    if((offsetRBS[i] >= 0) && (offsetRBS[j] >= 0) && (scoresRBS[i] < scoresRBS[j])) {
                        char *tempc = upsRegs[j];
                        upsRegs[j] = upsRegs[i];
                        upsRegs[i] = tempc;
                        int tempi = offsetRBS[j];
                        offsetRBS[j] = offsetRBS[i];
                        offsetRBS[i] = tempi;
                        double tempd = scoresRBS[j];
                        scoresRBS[j] = scoresRBS[i];
                        scoresRBS[i] = tempd;
                    }
                }
            }
            int j = 0;
            for(i = 0; i < numGenes; i++) {
                if(offsetRBS[i] < 0)
                    continue;
                if(j > (numCurrentgenes*percentage)/100.0)
                    offsetRBS[i] = - 1; // do not take them for the next iteration
                j++;
            }

            if(firstStageDone)
                buildRBSVector(numGenes, offsetRBS);
            buildRBSweightMatrix(upsRegs, offsetRBS, numGenes, weightFreq);
        }
        delete [] scoresRBS;
    }


    void Organism::computegenomFreqs() {
        int total = 0;
        for(list<Sequence *>::iterator cit = annotSeqs.begin(); cit != annotSeqs.end(); cit++) {
            Sequence &c = *(*cit);
            char *annotSeqSeq = c.getSeq(STRAND_FWD);
            for(int i=0;i < c.length();i++) {
                int chr = sigma->indC(annotSeqSeq[i]);

                if(chr >=0) {
                    countgenomicBases[chr]++;
                    total++;
                }
            }
        }
    }


    void Organism::computeRBSprofile(int lenSD, TGene cond, int windowLen, int percentage,  TSetType geneSet) {
        list<Gene> *listO = &genes[geneSet];

        // Initialization
        int i,j;
        int numGenes = listO->size();
        lenRBS = lenSD;
        windowRBS = windowLen;
        double **posFreq = new double*[sigma->alphabetSize() + 1];
        double **weightFreq = new double*[sigma->alphabetSize() + 1];
        double *H = new double[windowRBS], *accH = new double[windowRBS];
        int *offsetRBS = new int[numGenes];
        char **upsRegs = new char*[numGenes + 1];
        int maxHPos = 0;
        assert(H && posFreq && weightFreq && upsRegs);
        for(i=0;i < sigma->alphabetSize();i++) {
            weightRBS[i] = new double[lenRBS];
            assert(weightRBS[i]);
        }
        vectorRBS = new double[windowRBS];
        assert(vectorRBS);
        for(i=0 ; i< windowRBS ; i++)
            vectorRBS[i] = 0;

        for(i=0;i<=sigma->alphabetSize();i++) {
            weightFreq[i] = new double[lenRBS + 1];
            posFreq[i] = new double [windowRBS + 1];
            assert(posFreq[i] && weightFreq[i]);

            for(j=0;j <= windowRBS;j++)
                posFreq[i][j] = 0;

            for(int j =0 ; j <= lenRBS ; j++)
                weightFreq[i][j] = 0;
        }

        for(i=0;i < windowRBS;i++) {
            accH[i] = 0;
            H[i] = 0;
        }
        // Computing the frequency table
        i = 0;
        char *upsReg = new char[windowRBS + 1];

        for(list<Gene>::iterator it = listO->begin(); it != listO->end();it++, i++) {

            if(!(it->state() == cond)) {
                upsRegs[i] = NULL;
                continue;
            }

            (it->getSequence())->getLocation(upsReg, windowRBS, -abs(it->begin() - it->end()) - 1, it->begin(), it->end());

            if(strlen(upsReg) < (unsigned int)windowRBS) {
                upsRegs[i] = NULL;
                continue;
            }

            upsRegs[i] = strdup(upsReg);
            CountUtils::fillbaseTable(posFreq,windowRBS,1,upsRegs[i], 1, sigma);
        }

        for(i=0;i<= sigma->alphabetSize();i++) {
            for(j= 0;j <= windowRBS;j++)
                cout << posFreq[i][j]  << " ";
            cout << "\n";
        }

        // Using PI(positional information) to get the initial alignment
        {
            double *posfreqBase = new double[sigma->alphabetSize()];
            double *totalfreqBase = new double[sigma->alphabetSize()];
            for(j= 1; j <=sigma->alphabetSize(); j++)
                totalfreqBase[j - 1] = posFreq[j][0];

            for(i=1 ; i <= windowRBS; i++) {
                for(j= 1; j <=sigma->alphabetSize(); j++)
                    posfreqBase[j - 1] = posFreq[j][i];

                H[i-1] = CountUtils::PI(posfreqBase, totalfreqBase, sigma);

                for(j = ((i - lenRBS) <= 0 ? 0 : (i - lenRBS)); j < i;j++)
                    accH[j] += H[i-1];
            }
        }
        cout << "\nH";
        for(i =0 ; i < windowRBS; i++) {
            cout << H[i] << " ";
        }
        cout << "\nACCH";
        for(i =0 ; i < windowRBS; i++) {
            cout << accH[i] << " ";
        }
        cout << "\n";

        // getting the maximum value

        for(i=0 ; i < (windowRBS - lenRBS + 1) ; i++) {
            if(accH[i] < accH[maxHPos])
                maxHPos = i;
        }
        for(i = 0; i < numGenes;i++)
            if(!upsRegs[i])
                offsetRBS[i] = -1;
            else
                offsetRBS[i] = maxHPos;

        cout << "Initial Offset = " << maxHPos << "\n";
        // Calculating the initial weight matrix

        buildRBSweightMatrix(upsRegs, offsetRBS, numGenes, weightFreq);
        refineRBSweightMatrix(upsRegs, offsetRBS, numGenes, weightFreq, percentage);

        // freeing memory

        for(i=0;i<=sigma->alphabetSize();i++) {
            delete [] posFreq[i];
            delete [] weightFreq[i];
        }
        for(i = 0; i < numGenes; i++) {
            if(upsRegs[i]) delete [] upsRegs[i];
        }

        delete [] posFreq;
        delete [] weightFreq;
        delete [] upsRegs;
        delete [] upsReg;
        delete [] offsetRBS;
        delete [] H;
        delete [] accH;
    }


    void Organism::storegenomFreqs(::ofstream *statStream) {
        ::ofstream *fd = statStream;
        for(int i=0;i<sigma->alphabetSize();i++) {
            fd->write((char *)&countgenomicBases[i], sizeof(double));
        }
    }

    void Organism::storeRBSprofile(::ofstream *statStream) {
        ::ofstream *fd = statStream;
        int i;
        fd->write((char *)&lenRBS, sizeof(int));
        fd->write((char *)&windowRBS, sizeof(int));

        for(i=0;i<windowRBS;i++)
            fd->write((char *)&vectorRBS[i], sizeof(double));

        for(i=0;i<sigma->alphabetSize();i++) {
            for(int j=0;j<lenRBS;j++)
                fd->write((char *)&weightRBS[i][j], sizeof(double));
        }
    }

    void Organism::storeStats(::ofstream * statStream) {
        storegenomFreqs(statStream);
        storeRBSprofile(statStream);
    }

    void Organism::printStats() {
        int i,j;

        cerr << "\nOffset weight vector for " << windowRBS << " bases upstream the start codon\n";
        for(i=0;i<windowRBS;i++)
            cerr << i-windowRBS << " ";
        cerr << "\n";
        for(i=0;i<windowRBS;i++)
            cerr << vectorRBS[i] << " ";

        cerr << "\nPosition Weight matrix using " << lenRBS << " as the expected SD signal length\n";
        for(j=0;j<lenRBS;j++)
            cerr << " " <<j-lenRBS ;
        cerr << "\n";
        for(i=0;i<sigma->alphabetSize();i++) {
            cerr << sigma->contC(i);
            for(j=0;j<lenRBS;j++)
                cerr << " " << weightRBS[i][j];
            cerr << "\n";
        }

        //  cerr << "Percentages for each start\n";
        //  for(i=0; i < (int)startCodons.size(); i++)
        //    cerr << startCodons[i] << " " << countstartCodons[i] << "\n";

        cerr << "GC content in frequency\n";
        for(i=0;i<sigma->alphabetSize();i++)
            cerr << sigma->contC(i) << " " << countgenomicBases[i] << "\n";
    }

    void Organism::retrievegenomFreqs(::ifstream *statStream) {
        //  cout << "FREQ\n";
        ::ifstream *fd = statStream;
        for(int i=0;i<sigma->alphabetSize();i++) {
            fd->read((char *)&countgenomicBases[i], sizeof(double));
            //    cout << countgenomicBases[i] << "\t";
        }
        //  cout << "\n";
    }


    void Organism::retrieveRBSprofile(::ifstream *statStream) {
        int i;
        ::ifstream *fd = statStream;
        fd->read((char *)&lenRBS, sizeof(int));
        fd->read((char *)&windowRBS, sizeof(int));
        //  cout << windowRBS<< "\n";
        vectorRBS = new double[windowRBS];
        assert(vectorRBS);

        for(i=0;i<windowRBS;i++)
            fd->read((char *)&vectorRBS[i], sizeof(double));

        for(i=0;i<sigma->alphabetSize();i++) {
            weightRBS[i] = new double[lenRBS];
            assert(weightRBS[i]);

            for(int j=0;j<lenRBS;j++)
                fd->read((char *)&weightRBS[i][j], sizeof(double));
        }
    }

    void Organism::retrieveStats(::ifstream *statStream) {
        //  retrievestartFRQ(statStream);
        retrievegenomFreqs(statStream);
        retrieveRBSprofile(statStream);
    }


    void Organism::changeStateGenes(TGene state,  TSetType geneSet) {
        list<Gene> *listO = &genes[geneSet];

        for(list<Gene>::iterator it = listO->begin(); it != listO->end() ; it++)
            it->setState(state);
    }


    Organism::~Organism() {
        unsigned int i;

        if(vectorRBS)
            delete [] vectorRBS;

        for(i=0; (int)i < sigma->alphabetSize();i++)
            if(weightRBS[i]) delete [] weightRBS[i];

        annotSeqs.clear();
        listOverlaps.clear();

        for(int i = 0; i < NUM_SET_TYPES; i++)
            genes[i].clear();
    }

}
