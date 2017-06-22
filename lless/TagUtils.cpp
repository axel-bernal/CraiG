#include "TagUtils.h"
#include "Utils.h"
#include "FSM.h"
#include "FeatureEngine.h"
#include "boost/regex.hpp"
#include "NodeInst.h"
#include "EdgeInst.h"

/****************************************************************************
 * TagUtils.cpp - part of the lless namespace, a general purpose
 *                linear semi-markov structure prediction library
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
namespace lless {

    static std::map<std::string, TTagClass> tagClasses;

    TagUtils _tmpTagUtilsObj;

    TagUtils::TagUtils() {
        tagClasses["NODE_INST"] = NODE_INST;
        tagClasses["EDGE_INST"] = EDGE_INST;
        tagClasses["WORD_INST"] = WORD_INST;
    }

    int TagUtils::tagType(Tag *tag) {
        switch(tag->getGEClass()) {
        case EDGE_INST: return ((EdgeInst*)tag)->getType();
        case NODE_INST: return ((NodeInst*)tag)->getType();
        case WORD_INST: return 0;
        }
        assert(0);
    }

    TTagClass TagUtils::stringToTTagClass(const std::string &s) {
        std::map<std::string, TTagClass>::iterator it;
        it = tagClasses.find(s);
        if(it == tagClasses.end())
            throw EXCEPTION( PARSE_ERROR, std::string("undefined TTagClass ")+s);
        return it->second;
    }

    void TagUtils::loadSeqTags(std::string tagFile,
                               list<Sequence *> & listC,
                               FSM &fsm,
                               TSetType set, bool add2List) {

        ::ifstream tagStream(tagFile.c_str());

        if(!tagStream) {
            assert(0);
            throw EXCEPTION( BIOFEATS_UNAVAILABLE, string(tagFile));
        }

        loadTags((istream &)tagStream, listC, fsm, set, add2List);

    }


    void TagUtils::sloadSeqTags(std::string &tagString,
                                list<Sequence *> & listC,
                                FSM &fsm,
                                TSetType set, bool add2List) {

        std::istringstream tagStream(tagString);
        loadTags((istream &)tagStream, listC, fsm, set, add2List);

    }

    void TagUtils::loadTags(std::istream &tagStream,
                            list<Sequence *> & listC,
                            FSM &fsm,
                            TSetType set, bool add2List) {

        std::string line;
        Sequence *annotSeqObj = NULL;
        int pos;
        int len;

        boost::RegEx rExSequence("^>(\\S+)\\s*(\\d*)$");
        boost::RegEx rExLabel("^Label\\s+Set\\s+(\\d.*)$");
        boost::RegEx rExNodeInst("^NODE_INST\\s+(\\S+)\\s+(\\d+)\\s+(\\d.*)$");
        boost::RegEx rExEdgeInst("^EDGE_INST\\s+(\\S+)\\s+(\\d+)\\s*$");

        boost::RegEx rExDelim("\\s+");

        std::string seqId;
        int rank = -1;
        double score = 0;
        int initPhase = 0, endPhase = 0;
        int tagStats[2] = {0,0};
        Tag *tag;
        SeqTags *currTags = NULL;

        while(std::getline(tagStream, line), !tagStream.eof()) {
            tag = NULL;
            int phaseBreak = -1;

            if(rExNodeInst.Match(line)) {
                std::string nodeName(rExNodeInst[1]);
                Node *node = fsm.node(nodeName);

                if(!sscanf(rExNodeInst[2].c_str(), "%d", &pos))
                    assert(0);
                int numScans = sscanf(rExNodeInst[3].c_str(), "%d %d",
                                      &len, &phaseBreak);

                assert(numScans);
                if(numScans < 2)  phaseBreak = -1;

                tag = new NodeInst(node->id3(), node->id(), pos, node->strand(), len);
                tagStats[0]++;
            }
            else if(rExEdgeInst.Match(line)) {
                std::string edgeName(rExEdgeInst[1]);
                Edge *edge = fsm.edge(edgeName);

                if(!sscanf(rExEdgeInst[2].c_str(), "%d", &pos))
                    assert(0);

                tag = new EdgeInst(edge->id2(), edge->id(), pos, edge->strand());
                tagStats[1]++;
            }
            else if(rExLabel.Match(line)) {

                int numScans = sscanf(rExLabel[1].c_str(), "%d %lf %d %d",
                                      &rank, &score, &initPhase, &endPhase);

                assert(numScans);
                currTags = &annotSeqObj->addTags(rank, set);
                if(numScans == 4)
                    currTags->initialize(score, initPhase, endPhase, rank);

                continue;
            }
            else if(rExSequence.Match(line))  {
                seqId = rExSequence[1];
                int length = 0;
                if(rExSequence[2].length())
                    sscanf(rExSequence[2].c_str(), "%d", &length);

                annotSeqObj = (Sequence *)BasicSeq::findSequence(seqId, (list<BasicSeq *> &)listC);

                if(!annotSeqObj && add2List) {
                    annotSeqObj = new Sequence(seqId);
                    annotSeqObj->setLength(length);
                    listC.push_back(annotSeqObj);
                }
                else length = annotSeqObj->length();

                if(!annotSeqObj || !length || length != annotSeqObj->length())
                    throw EXCEPTION(TAGFILE_FMTERROR, seqId);

                continue;
            }
            else if(line.length() > 0) {
                assert(0);
                throw EXCEPTION(PARSE_ERROR, line);
            }
            else continue;

            if(!annotSeqObj) {
                assert(tag && annotSeqObj && rank >= 0 && currTags);
                throw EXCEPTION(CONTIG_UNAVAILABLE, seqId);
            }

            currTags->tpush_back(tag);
            if(phaseBreak > 0)
                currTags->insertPhaseBreak(currTags->back(), phaseBreak);

        }

        if(!add2List) {
            cerr << " Tags loaded for Set " << set << ".....\n";
            cerr  << "\t" << tagStats[0] << " Node Instances\n";
            cerr  << "\t" << tagStats[1] << " Edge Instances\n";
        }
    }


    void TagUtils::saveSeqTags(std::ofstream &fd,
                               Sequence & c,
                               FSM &fsm,
                               TSetType set) {

        fd << ">" << c.id() << " " << c.length() << endl;
        saveTags((ostream &)fd, c.getTags(set), fsm);
    }

    void TagUtils::ssaveSeqTags(std::string &tagString,
                                Sequence & c,
                                FSM &fsm,
                                TSetType set) {

        std::ostringstream tagStream(tagString);

        tagStream << ">" << c.id() << " " << c.length() << endl;
        saveTags((ostream &)tagStream, c.getTags(set), fsm);
        tagString = tagStream.str();
    }

    void TagUtils::saveTags(std::ostream &fd,
                            vector<SeqTags> & tags,
                            FSM &fsm) {

        for(unsigned int i = 0; i < tags.size(); i++) {
            fd << "Label Set " << tags[i].rank() << " " << tags[i].score() << " " <<
                tags[i].initPhase() << " " << tags[i].endPhase() << endl;

            for(SeqTags::iterator git = tags[i].begin(); git != tags[i].end(); git++) {

                if((*git)->getGEClass() == NODE_INST) {
                    NodeInst *ni = (NodeInst *)(*git);
                    Node *node = fsm.node((TParseNode)ni->getParseType());
                    fd << "NODE_INST " << node->name() << " " << ni->getPos() << " " << ni->getLen();
                    int phaseBreak = tags[i].phaseBreak(*git);

                    if(phaseBreak >= 0)
                        fd << " " << phaseBreak;

                }
                else {
                    EdgeInst *ei = (EdgeInst *)(*git);
                    Edge *edge = fsm.edge((TParseEdge)ei->getParseType());
                    fd << "EDGE_INST " << edge->name() << " " << ei->getPos();
                }
                fd << endl;
            }
        }
    }

}
