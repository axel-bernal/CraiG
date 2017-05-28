/****************************************************************************
 * SequenceUtils.h - part of the lless namespace, a general purpose
 *                   linear semi-markov structure prediction library
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

#ifndef _ANNOT_SEQ_UTILS_H_
#define _ANNOT_SEQ_UTILS_H_

#include "Sequence.h"
#include <string>
#include <fstream>
#include <boost/regex.hpp>

using namespace lless;

namespace lless {
    class Sigma;

    /**
     * SequenceUtils is a class that contains utility subroutines for Sequence
     * type objects.
     * The only formats that are allowed for an AnnotSe object so far are FASTA
     * and XFASTA.The latter one is a compact representation of the
     * former which reduces loading times tremendously in those cases in which
     * the input sequence contains many repeated characters/symbols.
     ***************************************************************************/

    class SequenceUtils {

    public:

        SequenceUtils();

        static TSetType stringToTSetType(const std::string &s);

        static Sequence* loadFastaSequence(std::ifstream & annotSeqStream,
                                           Sigma *sigma,
                                           std::string &line);
        static Sequence* loadTagSequence(std::ifstream & annotSeqStream, FSM &fsm,
                                         std::string &line);
        static Sequence* loadExtFastaSequence(std::ifstream & annotSeqStream,
                                              Sigma *sigma,
                                              std::string &line);
        static Sequence* loadMultizSequence(std::ifstream & annotSeqStream,
                                            Sigma *sigma,
                                            std::string &line);
        static EdgeAnnotSeq* loadEdgeAnnotSequence(std::ifstream & annotSeqStream,
                                                   Sigma *sigma,
                                                   std::string &line);

        /**
         * A static member function that load scores associated with
         * a sequence. The scores are specified at each position of the
         * sequence starting from 1.
         * @param fd the file containing the scores.
         * @param line the string containing the fasta header which
         * needs to be parsed to obtain the sequence id
         */
        template <class TClass>
            static ScoreSeq<TClass> * loadExtScoreSequence(std::ifstream &fd,
                                                           std::string &line) {

            boost::RegEx rExScoreId("^>(\\S+)\\s(\\d+)\\s*$");
            boost::RegEx rExCompSeqMark("^\\s*\\/\\/\\s*$");

            TStrand strand = STRAND_FWD;
            ScoreSeq<TClass> *c = NULL;
            int seqLen = 0;

            while(1) {

                if(rExScoreId.Match(line) || fd.eof())  {
                    if(c) {
                        assert(seqLen);

                        if(strand == STRAND_FWD)
                            c->setNatComplementSeq();

                        return c;
                    }

                    if(fd.eof())
                        break;

                    std::string cId = rExScoreId[1];
                    if(!sscanf((char *)rExScoreId[2].c_str(), "%d", &seqLen))
                        assert(0);

                    c = new ScoreSeq<TClass>(cId, seqLen);
                    strand = STRAND_FWD;
                }
                else if(rExCompSeqMark.Match(line))
                    strand = STRAND_COMP;
                else { // it must match a line="pos score" or "pos1..pos2 score"
                    int start, end;
                    double score;
                    if(sscanf(line.c_str(), "%d..%d %lf", &start, &end, &score) != 3) {
                        if(sscanf(line.c_str(), "%d %lf", &start, &score) != 2) {
                            assert(0);
                            throw EXCEPTION(STRANGE_CHAR, line);
                        }
                        end = start;
                    }
                    for(int i = start; i <= end; i++)
                        (*c)(strand, i) = (TClass)score;
                }
                std::getline(fd, line);
            }

            return NULL;
        }

        template <class TClass>
            static MultiScoreSeq<TClass> *loadExtMultiScoreSequence(std::ifstream &fd,
                                                                    std::string &line) {

            boost::RegEx rExMultiScoreId("^>(\\S+)\\s+(\\d+)\\s+(\\d+)\\s*$");
            boost::RegEx rExCompSeqMark("^\\s*\\/\\/\\s*$");

            TStrand strand = STRAND_FWD;
            MultiScoreSeq<TClass> *c = NULL;
            vector<ScoreSeq<TClass> *> subcs;

            int seqLen = 0, numColumns = 1, i;

            while(1) {

                if(rExMultiScoreId.Match(line) || fd.eof())  {

                    if(c) {
                        assert(subcs.size() == numColumns && seqLen);

                        for(i = 0; i < numColumns; i++)
                            c->setSeq(i, subcs[i]);

                        if(strand == STRAND_FWD)
                            c->setNatComplementSeq();

                        return c;
                    }

                    if(fd.eof())
                        break;

                    std::string cId = rExMultiScoreId[1];
                    if(!sscanf((char *)rExMultiScoreId[2].c_str(), "%d", &seqLen))
                        assert(0);
                    if(!sscanf((char *)rExMultiScoreId[3].c_str(), "%d", &numColumns))
                        assert(0);

                    c = new MultiScoreSeq<TClass>(cId, numColumns);
                    subcs = vector<ScoreSeq<TClass> *>(numColumns);

                    for(i = 0; i < numColumns; i++)
                        subcs[i] = new ScoreSeq<TClass>(cId, seqLen);

                    strand = STRAND_FWD;
                }
                else if(rExCompSeqMark.Match(line))
                    strand = STRAND_COMP;
                else { // it must match a line="pos score"
                    int pos;
                    std::istringstream lineStream(line);
                    lineStream >> pos;

                    for(i = 0; i < numColumns; i++) {
                        double score;

                        if(lineStream.eof()) {
                            assert(0);
                            throw EXCEPTION(STRANGE_CHAR, line);
                        }

                        lineStream >> score;
                        (*subcs[i])(strand, pos) = (TClass)score;

                    }
                }
                std::getline(fd, line);
            }

            return NULL;
        }


        /**
         * Not supported yet.
         * \todo Implement this function
         */
        template <class TClass>
            static ScoreSeq<TClass>* loadScoreSequence(std::ifstream &fd,
                                                       std::string &line) {
            return NULL;
        }

        static void loadSeqIndexes(std::ifstream &fd,
                                   map<std::string, streampos> & seqIndexes);
        ~SequenceUtils() { }
    };
}

#endif
