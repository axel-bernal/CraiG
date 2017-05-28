/****************************************************************************
 * TagPrinter.h - part of the lless namespace, a general purpose
 *          linear semi-markov structure prediction library
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

#ifndef _TAG_PRINTER_H_
#define _TAG_PRINTER_H_
#include "FSM.h"
#include "Sequence.h"


namespace lless {

    // forward declarations
    class FeatureEngine;

    /**
     * The TagPrinter class contains the abstract methods to display tags. These
     * methods are to be overriden by derived classes.
     ***************************************************************************/

    class TagPrinter {
    public:
        TagPrinter() {;}

        virtual void displayTags(std::string & format,
                                 std::ofstream & os,
                                 Sequence & c,
                                 SeqTags & seqTags,
                                 FSM &fsm,
                                 const char *idPrefix,
                                 int offset = 0,
                                 int upperLimit = INT_MAX) = 0;

        virtual void displayTags(std::string & format,
                                 std::ofstream & os,
                                 Sequence & c,
                                 FSM &fsm,
                                 const char *idPrefix,
                                 TSetType set,
                                 int offset = 0,
                                 int upperLimit = INT_MAX) = 0;


        virtual ~TagPrinter() {;}

    };

}

#endif
