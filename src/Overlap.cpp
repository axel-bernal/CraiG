#include "Overlap.h"

/****************************************************************************
* Overlap.cpp - part of the craig namespace, a genomics library
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

  Overlap::Overlap() {
    olapregROrf = NULL;
    olapregLOrf = NULL;
  }

  Overlap::Overlap(char *orRO, char *orLO,
                   int posRO, int posLO, 
                   int frROrf, int frLOrf) {

    olapregROrf = strdup(orRO);
    olapregLOrf = strdup(orLO);
    assert(olapregROrf && olapregLOrf);
    iniposLOrf = posLO;
    iniposROrf = posRO;
    frameLOrf = frLOrf;
    frameROrf = frROrf;

    for(int i=0;i< 7;i++)
      score[i] = -1; 
  }

  Overlap::Overlap(Overlap & ov) {
    *this = ov;
  }                  

  Overlap& Overlap::operator= (Overlap& ov) {
    if(ov.getolapregRO() && ov.getolapregRO()) {
      olapregROrf = strdup(ov.getolapregRO());
      olapregLOrf = strdup(ov.getolapregLO());
      assert(olapregROrf && olapregLOrf);
    }
    else {
      olapregROrf = NULL;
      olapregLOrf = NULL;
    }
    iniposLOrf = ov.getposLO();
    iniposROrf = ov.getposRO();
    frameLOrf = ov.getframeLO();
    frameROrf = ov.getframeRO();

    for(int i=0;i < 7;i++)
      score[i] = ov.getolapScore(i);

    return *this;
  }                 

  void Overlap::scoreolapRegion(IMM *imm, int amb, int indep) {
    assert(olapregLOrf && olapregROrf);
    int frLO = abs(frameLOrf) + 3*(frameLOrf < 0) , 
      frRO = abs(frameROrf) + 3*(frameROrf < 0), frInd = 0;

    assert(frLO >= 0 && frLO <= 6 && frRO >= 0 && frRO <= 6);
    if(score[frLO] <= 0)
      score[frLO] =  imm->seqScore(olapregLOrf, amb, 
                                   iniposLOrf, 
                                   strlen(olapregLOrf));
    if(score[frRO] <= 0)
      score[frRO] =  imm->seqScore(olapregROrf, amb, 
                                   iniposROrf, 
                                   strlen(olapregROrf));
       
    if(indep && score[frInd] <= 0)
      score[frInd] = Utils::seqscoreIndep(olapregROrf,  strlen(olapregROrf));

  }


  Overlap::~Overlap() {

    if(olapregLOrf)
      delete olapregLOrf;

    if(olapregROrf) 
      delete olapregROrf;
  }
}
