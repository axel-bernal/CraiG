/****************************************************************************
 * WordInst.h - part of the lless namespace, a general purpose
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

#ifndef _WORD_INST_H_
#define _WORD_INST_H_

#include "Utils.h"
#include "Tag.h"
#include "FSM.h"

namespace lless {

    /**
     * The WordInst class represents instances of lattice nodes in the sequence.
     * it derives from the Tag class the parseType and position among other
     * things. It adds up to this, a length and an internal type.
     ***************************************************************************/
    class WordInst : public Tag {
    protected:
        SeqTags *stObj;
        SeqTags::iterator bwInst, ewInst;

    public:
        /**
         * Main constructor
         * @param pType the parsing type of the word, one of TParseWord enumerate.
         * @param pos the starting position of the word in the sequence.
         * @param strand the strand of the word.
         * @param len the length of the word.
         */
    WordInst(Word &word,
             SeqTags &stObj,
             SeqTags::iterator &it,
             FSM *fsm
        ) : Tag((TParseWord)word.id()) {

            this->stObj = &stObj;
            bwInst = (word.prevNode() == SYNC_STATE) ?
                it : stObj.end();
            ewInst = stObj.end();
            this->len = 0;
            SeqTags::iterator last_it = it;

            for( ; it != stObj.end(); it++) {
                if((*it)->getGEClass() == NODE_INST) {
                    NodeInst &b = (NodeInst &)(**it);
                    TNodeId2 bId2 = fsm->node(b.getParseType())->id2();

                    if(word.prevNode() == bId2)
                        if(bwInst == stObj.end())
                            bwInst = ++it; // jump to edge

                    if(word.nextNode() == bId2) {
                        ewInst = last_it;
                        break;
                    }

                    if(bwInst != stObj.end())
                        this->len += b.getLen();

                }
                else
                    last_it = it;
            }

            if(word.nextNode() == SYNC_STATE)
                ewInst = last_it;

            if(bwInst == stObj.end() || ewInst == stObj.end()) {
                this->pType = INVALID_WORD;
                return;
            }

            this->pos = (*bwInst)->getPos();
            this->strand = (*bwInst)->getStrand();

            this->geClass = WORD_INST;
        }

        /**
         * Copy constructor
         */
        WordInst(const WordInst & st) {
            *this = st;
        }

        inline SeqTags * parse() {
            return stObj;
        }

        ~WordInst() {
        }
    };

}

#endif
