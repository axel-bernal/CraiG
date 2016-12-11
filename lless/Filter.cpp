#include "Filter.h"

/****************************************************************************
* Filter.cpp - part of the lless namespace, a general purpose
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


namespace lless {
  //! Default constructor
  Filter::Filter(int fInd,           //!< A unique filter identifier
                 std::string & name, //!< A unique name.
                 int maxFilterVals   /*!< The number of possible filter values.
                                       It is zero if filter values are not 
                                       countable. */
                 ) {
    this->fInd = fInd;
    this->arrInd = -1;
    this->name = name;
    this->maxFilterVals = maxFilterVals; 
    this->defaultStrand = STRAND_FWD;
    this->sparse = false;
    this->logScaled = false;
    this->exponent = 1;
    this->_period = 0;
  }
  
  /**
   * Constructor from a Header string definition
   */
  Filter::Filter(int fInd,         //!< A unique filter identifier
                 /*! The Header string definition, loaded as a vector 
                   of strings. */
                 vector<std::string> & params, 
                 int & offset,     //!< The index for vector params.
                 int maxFilterVals /*!< The number of possible filter values. 
                                     It is zero if filter values are not 
                                     countable. */
                 ) {
    this->fInd = fInd;
    this->arrInd = -1;
    this->name = params[offset++];
    this->maxFilterVals = maxFilterVals; 
    this->defaultStrand = STRAND_FWD;
    this->sparse = false;
    this->logScaled = false;
    this->exponent = 1;
    this->_period = 0;
  }
  
}
