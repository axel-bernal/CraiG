/****************************************************************************
 * GenLibExcept.h - part of the lless namespace, a general purpose
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

#ifndef _GEN_EXCEPTION_H_
#define _GEN_EXCEPTION_H_

#include <string>
#include <sstream>
#include <string.h>
#include <stdexcept>
#include <signal.h>
#include <execinfo.h>

using namespace std;

#define EXCEPTION(name, msg) GenLibExcept((name), (msg), __FILE__, __LINE__)

namespace lless {

    /**
     * The enumerate that registers the annotated types of error found in craig
     * and less
     */
    typedef enum {BAD_LOCATION_FORMAT,
                  INVALID_CONTIG_COORDS,
                  CONTIG_UNAVAILABLE,
                  BIOFEATS_UNAVAILABLE,
                  STRANGE_CHAR,
                  BAD_CONTIG_SCORE_HEADER,
                  INCONSISTENT_CONTIG,
                  STATES_NOT_GENERATED,
                  INVARIANT,
                  GCCLASS_UNDEFINED,
                  LATTICE_DECODING,
                  CONTEXT_TOO_LONG,
                  FORB_VIRTUAL_FUNCTION,
                  PARSE_ERROR,
                  INCOMPAT_GENE,
                  INCOMPAT_TRANSCRIPT,
                  MISSING_ARGUMENT,
                  OUT_OF_MEMORY,
                  FILE_UNAVAILABLE,
                  TAGFILE_FMTERROR,
                  BAD_RESOURCE_HEADER,
                  BAD_FILTER_HEADER,
                  BAD_FEATURE_HEADER,
                  INSUFFICIENT_ARGS,
                  BAD_USAGE,
                  NOT_SUPPORTED,
                  NOT_ANNOTATED
    } TErrorName;

    static int NUM_ERRORS = 27;
    static string errorMsg[] = {
        "bad format in BioFeature object's location ",
        "invalid sequence coordinates ",
        "sequence not available ",
        "BioFeature objects not available ",
        "unrecognized symbol in sequence ",
        "bad format for Sequence object ",
        "sequence information is bogus ",
        "must generate states at",
        "invariant condition violated ",
        "gc-class is not defined in training set ",
        "lattice decoding general error ",
        "attempting to use a markov context longer than allowed ",
        "attempting to use a virtual function that is forbidden to execute ",
        "parse general error ",
        "coordinates of gene are in wrong format ",
        "coordinates of transcript are in wrong format ",
        "required argument has not bee provided to subroutine ",
        "out of memory ",
        "cannot find/open/delete file ",
        "error reading format in tag file ",
        "malformed resource header ",
        "malformed filter header ",
        "malformed feature header ",
        "insufficient arguments ",
        "bad usage ",
        "not supported ",
        ""
    };


    /**
     * The GenLibExcept class, an exception class for lless and craig
     ***************************************************************************/


    class GenLibExcept : public exception {
    protected:
        int errorIndex;
        string excepMsg;
    public:

        GenLibExcept(TErrorName errorName,
                     std::string addMsg,
                     const char *file = "NO FILE SPECIFIED",
                     const long line = 0) {

            errorIndex  = errorName;

            if(errorIndex < 0 || errorIndex >= NUM_ERRORS)
                errorIndex = NUM_ERRORS - 1;

            std::ostringstream buffer;
            buffer << file << ":" << line << " " << errorMsg[errorIndex] << ". \"" << addMsg << "\"\n";
            show_stackframe(buffer);
            excepMsg = buffer.str();
        }

        void show_stackframe(ostream &ost) {
            void *trace[16];
            char **messages = (char **)NULL;
            int i, trace_size = 0;

            trace_size = backtrace(trace, 16);
            messages = backtrace_symbols(trace, trace_size);
            ost << "[bt] Execution path:\n";
            for (i=0; i<trace_size; ++i)
                ost << "[bt] " << messages[i] << endl;
        }

        inline int error() {
            return errorIndex;
        }

        GenLibExcept() {
            errorIndex = -1;
        }

        ~GenLibExcept() throw() {
        }

        /**
         * @ return the rror description
         */
        virtual const char* what() const throw() {
            return excepMsg.c_str();
        }
    };

}

#endif
