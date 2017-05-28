/****************************************************************************
 * TagUtils.h - part of the lless namespace, a general purpose
 *              linear semi-markov structure prediction library
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

#ifndef _TAG_UTILS_H_
#define _TAG_UTILS_H_
#include <string>
#include "Sequence.h"
#include <fstream>


namespace lless {
    class FSM;
    class FeatureEngine;
    class Tag;
    /**
     * The TagUtils class contains utilities for the Tag class
     ***************************************************************************/

    class TagUtils {
    public:
        TagUtils();
        static TTagClass stringToTTagClass(const std::string &s);
        static int tagType(Tag *tag);
        static void loadSeqTags(std::string tagFile,
                                list<Sequence *> & listC,
                                FSM &fsm,
                                TSetType set = DEFAULT_SET,
                                bool add2List = false);


        static void sloadSeqTags(std::string &tagString,
                                 list<Sequence *> & listC,
                                 FSM &fsm,
                                 TSetType set = DEFAULT_SET,
                                 bool add2List = false);

        static void loadTags(std::istream &tagStream,
                             list<Sequence *> &,
                             FSM &fsm,
                             TSetType = DEFAULT_SET,
                             bool add2List = false);

        static void saveSeqTags(std::ofstream &fd,
                                Sequence & c,
                                FSM &fsm,
                                TSetType set = DEFAULT_SET);

        static void ssaveSeqTags(std::string &tagString,
                                 Sequence & c,
                                 FSM &fsm,
                                 TSetType set = DEFAULT_SET);

        static void saveTags(std::ostream &fd,
                             vector<SeqTags> & tags,
                             FSM &fsm);

        static pair<int,int> tagCoords(Node *node, Tag *tag) {
            TStrand strand = tag->getStrand();
            if(strand == STRAND_COMP)
                node = node->complementNode();

            int fivePend = tag->getPos() + node->prevEdgePos();
            int threePend = tag->getPos() + tag->getLen() + node->nextEdgePos();
            //      cerr << "tagcoords " << node->name() << " " << fivePend << " " << threePend << endl;
            return pair<int, int>(fivePend, threePend);

        }

        ~TagUtils() { }
    };

}

#endif
