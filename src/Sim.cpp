#include "Utils.h"
#include "Sim.h"

/****************************************************************************
* Sim.cpp - part of the craig namespace, a genomics library
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

namespace craig {

  Sim::Sim(char *g1, char *g2,
           int l1, int l2, 
           int perc, double sc, 
           int rIni1, int rEnd1, 
           int rIni2, int rEnd2) {

    geneId1 = strdup(g1);
    geneId2 = strdup(g2);
    assert(geneId1 && geneId2);
    lenId1 = l1;
    lenId2 = l2;
    percentage = perc;  
    score = sc;
    regIni1 = rIni1;
    regEnd1 = rEnd1;
    regIni2 = rIni2;
    regEnd2 = rEnd2;
  }

  Sim::Sim(const Sim & sim) {
    Sim & s = (Sim &)sim;
    char *g1 = s.getGeneid1(), *g2 = s.getGeneid2();
    assert(g1 && g2);
    geneId1 = strdup(g1);
    geneId2 = strdup(g2);
    lenId1 = s.getLen1();
    lenId2 = s.getLen2();
    percentage = s.getPercentage();  
    score = s.getScore();
    regIni1 = s.getRegini1();
    regEnd1 = s.getRegend1();
    regIni2 = s.getRegini2();
    regEnd2 = s.getRegend2();
  }                  

  Sim& Sim::operator=(Sim& s) {
    char *g1 = s.getGeneid1(), *g2 = s.getGeneid2();
    assert(g1 && g2);
    geneId1 = strdup(g1);
    geneId2 = strdup(g2);
    lenId1 = s.getLen1();
    lenId2 = s.getLen2();
    percentage = s.getPercentage();  
    score = s.getScore();
    regIni1 = s.getRegini1();
    regEnd1 = s.getRegend1();
    regIni2 = s.getRegini2();
    regEnd2 = s.getRegend2();
    return *this;
  }                 

  Sim::~Sim() {
    if(geneId1)
      delete [] geneId1;

    if(geneId2)
      delete [] geneId2;
  }                 
}

