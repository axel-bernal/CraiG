/****************************************************************************
 * Sim.h - part of the craig namespace, a genomics library
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

#ifndef _SIMS_H
#define _SIMS_H

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <vector>


namespace craig {

    /**
     * The Sim class represents an exon similarity.
     ***************************************************************************/

    class Sim {

    private:
        char *geneId1, *geneId2; // the gene id's
        int lenId1, lenId2;
        int percentage;
        double score;
        int regIni1, regEnd1;
        int regIni2, regEnd2;

    public:
        Sim(char *g1, char *g2,
            int l1, int l2,
            int perc, double sc,
            int rIni1, int rEnd1,
            int rIni2, int rEnd2);

        Sim(const Sim & s);

        Sim& operator=(Sim & s);

        inline char *getGeneid1() {
            return geneId1;
        }

        inline char *getGeneid2() {
            return geneId2;
        }

        inline int getLen1() {
            return lenId1;
        }

        inline int getLen2() {
            return lenId2;
        }

        inline int getPercentage() {
            return percentage;
        }

        inline double getScore() {
            return score;
        }
        inline int getRegini1() {
            return regIni1;
        }
        inline int getRegend1() {
            return regEnd1;
        }

        inline int getRegini2() {
            return regIni2;
        }

        inline int getRegend2() {
            return regEnd2;
        }

        inline bool operator< ( Sim& s) {
            int b = strcmp(geneId1, s.getGeneid1());
            if(b)
                return b < 0;
            b = strcmp(geneId2, s.getGeneid2());
            if(b)
                return b < 0;
            return score < s.getScore();
        }

        inline bool operator== (Sim& s) {
            return (!strcmp(geneId1,s.getGeneid1()) &&
                    !strcmp(geneId2, s.getGeneid2()) &&
                    score == s.getScore());;
        }

        ~Sim();

    };
}

#endif
