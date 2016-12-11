/****************************************************************************
* LengthSegment.h - part of the lless namespace, a general purpose
*             linear semi-markov structure prediction library
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

#ifndef _FT_LENGTH_SEGMENT_H_
#define _FT_LENGTH_SEGMENT_H_

#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <set>
#include "Filter.h"
#include "FL_Signal.h"
#include "Utils.h"
#include "FilterEngine.h"
#include <sstream>


namespace lless {

  /**
   * FT_LengthSegment is a subclass of TypedFeature<int> and its value is the
   * length of a segment (or node).
   */
  class FT_LengthSegment : public TypedFeature<int> {

   public:
    /**
     * Default constructor
     * @param fInd A feature index.
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * @param name A unique feature identifier.
     * @param fe A pointer to a FilterEngine object.
     * @param maxNumFeatVals The number of possible feature values. 
     * If greater than zero then the Feature object's possible values are
     * countable.
     */
    FT_LengthSegment(
                  int fInd,
                  int paramInd,
                  int parsingFrames,
                  char *name,
                  FilterEngine *fe,
                  int maxNumFeatVals = 0
                  ) 
      : TypedFeature<int>(fInd, paramInd, parsingFrames, 
                          name, fe, 
                          maxNumFeatVals, FT_INTEGER) {
      ;
    }

    /**
     * Constructor from a Header string definition. The feature value type
     * is provided as parameter.
     * @param fInd A feature index.
     * @param paramInd feature's parameter vector index. If negative, the 
     * feature does not have any ties to Tag objects.
     * object to which this feature is tied to.
     * @param fargs The Header string definition loaded as a vector of strings.
     * The Header has the following form:\n\n Feature name LengthSegment 
     * numParsingFrames\n\n
     * The description of the fields could be found in the other
     * constructor(s)
     * @param offset The index for vector fargs.
     * @param fe A pointer to a FilterEngine object.
     * @param fte A pointer to a FeatureEngine object.
     */
    FT_LengthSegment(
                  int fInd, 
                  int paramInd, 
                  vector<std::string> & fargs, 
                  int & offset, 
                  FilterEngine *fe, 
                  FeatureEngine *fte
                  ) 
      : TypedFeature<int>(fInd, paramInd,
                          fargs, offset, 
                          fe, 1, FT_INTEGER) {

    }

    /**
     * @param ge Tag object to which this feature is tied.
     * @param frame Tag object's current phase(or frame).
     * @return the length of the Tag object given as parameter.
     */
    inline double featValue(Tag *ge, int frame) {
      return ge->getLen();
    }
  };

}

#endif
