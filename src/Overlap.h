/****************************************************************************
 * Overlap.h - part of the craig namespace, a genomics library
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

#ifndef _OVERLAP_H_
#define _OVERLAP_H_
#include "IMM.h"

namespace craig {

    /**
     * The Overlap class represents the physical overlap between two exons.
     ***************************************************************************/

    class Overlap {
    private:
        char *olapregLOrf, *olapregROrf;
        int iniposLOrf, iniposROrf;
        int frameLOrf, frameROrf;
        double score[7];

    public:
        Overlap(char *orRO, char *orLO,
                int posRO, int posLO,
                int frROrf, int frLOrf);

        Overlap();
        Overlap(Overlap & ov);

        void scoreolapRegion(IMM *imm, int amb, int indep);

        Overlap& operator= (Overlap& ov);

        inline char *getolapregLO() {
            return olapregLOrf;
        }
        inline char *getolapregRO() {
            return olapregROrf;
        }

        inline int getposLO() {
            return iniposLOrf;
        }
        inline int getposRO() {
            return iniposROrf;
        }
        inline int getframeLO() {
            return frameLOrf;
        }
        inline int getframeRO() {
            return frameROrf;
        }

        inline double getolapScore(int frame) {
            int fr = abs(frame) + 3*(frame < 0);
            assert(fr >= 0 && fr <= 6);
            return score[fr];
        }

        inline void setolapScore(double val, int frame) {
            int fr = abs(frame) + 3*(frame < 0);
            assert(fr >= 0 && fr <= 6);
            score[fr] = val;
        }

        inline bool operator< (Overlap& ov) {
            int b = strcmp(olapregROrf, ov.getolapregRO());
            return b <= 0;
        }

        inline bool operator== (Overlap& ov) {
            return !strcmp(olapregROrf, ov.getolapregRO());
        }

        ~Overlap();

    };
}

#endif
