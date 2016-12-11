#include <time.h>
#include <stdio.h>
#include "NGram.h"
#include "ResourceEngine.h"
#include "FilterEngine.h"

/****************************************************************************
* NGram.cpp - part of the lless namespace, a general purpose
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

using namespace std;
#include <string>
#include <stdexcept>


/**
 * Constructor
 * @param alphabet the alphabet used by the input sequence
 * @param maxOrder the maximum order of any contained FL_Gram object
 */
NGram::NGram(Sigma *alphabet, int maxOrder) {
  _initialize();
  _maxOrder = maxOrder;

  for(int i = 0; i < _maxOrder; i++) {
    std::string name = "FL_Gram";
    grams[i] = new FL_Gram<int>(0, name, alphabet, FT_INTEGER, i + 1);
  }

}

//! Copy constructor
NGram::NGram(NGram & ng) {
  _initialize();
  _maxOrder = ng.maxOrder();

  for(int i = 0; i < _maxOrder; i++) {
    if(ng.gram(i)) {
      grams[i] = new FL_Gram<int>(*ng.gram(i));
    }
  }

}

/**
 * Frequency table constructor
 * @param fs input file stream containing frequency table
 * @param alphabet the alphabet of the input sequences
 */
NGram::NGram(::ifstream *fs, Sigma *alphabet) {
  _initialize();
  std::ifstream &file = *fs;
  string::size_type index;
  int order = 0;
  _maxOrder = 0;
  std::string str;

  do {
    std::getline(file, str);
    index = str.find("-",0);

    if(str[0] == '\\' && index != string::npos)
      break;

  } while(!file.eof());
  
  while(!file.eof() && str[0] == '\\' && index != string::npos) {
    std::getline(file, str);

    if(file.eof()) break;
    index = str.find("-",0);
    order = 0;
    int lastOrder = 0;

    while(!file.eof() && (str[0] != '\\' || index == string::npos)) {
      char word[100] = "";
      int j = 0;

      for(unsigned int i = 0; i < str.length(); i++) {
        if(alphabet->indC(str[i]) >= 0)
          word[j++] = str[i];
      }
      word[j] = '\0';
      order = strlen(word);

      if(order) { // we found a valid word
        assert(!lastOrder || (lastOrder == order));
        lastOrder = order;

        if(gram(order - 1) == NULL) {
          std::string name = "FL_Gram";
          grams[order - 1] = new FL_Gram<int>(0, name, alphabet, FT_INTEGER, order, false);
          _maxOrder = order > _maxOrder ? order : this->_maxOrder;
        }
        grams[order - 1]->addWord(word);
      }
      std::getline(file, str);
      
      if(file.eof()) break;
      index = str.find("-", 0);
    }
  }
}

/**
 * A member function that stores the contents of *this in a file and calls
 * each contained FL_Gram object's own storeParams function.
 * @see FL_Gram::storeParams(::ofstream **)
 */
void NGram::storeParams(::ofstream **fd) {
  ::ofstream *file = *fd;
  file->write((char *)&_maxOrder, sizeof(int));
//  cerr << "order " <<  _maxOrder << " " << endl;
  for(int i = 0; i < _maxOrder; i++) {
    int wordSize = 0;

    if(grams[i])
      wordSize = grams[i]->order();

//    cerr << "wordsize " << i << " " << wordSize << endl;
    file->write((char *)&wordSize, sizeof(unsigned int));

    if(grams[i])
      grams[i]->storeParams(fd);
  }
}

/**
 * @see FL_Gram::operator+=(FL_Gram &)
 */
NGram & NGram::operator+=(NGram & ng) {
  _maxOrder = ng.maxOrder();
  for(int i = 0; i < ng.maxOrder(); i++) {

    if(ng.gram(i)) {
      if(!gram(i)) {
        assert(i < ng.maxOrder());
        grams[i] = new FL_Gram<int>(*(ng.gram(i)));
      }
      else {
        assert(grams[i]->maxNumFilterValues() == ng.gram(i)->maxNumFilterValues());
        (*grams[i]) += *(ng.gram(i));
      }
    }
  }
  return *this;
}

/**
 * @see FL_Gram::computeVals(char *, int *, int) 
 */
void NGram::generateSeqKmers(char *seq, int **seqFL_Grams, int len) {
  for(int i = 0; i < _maxOrder; i++) {
    if(!grams[i])
      continue;
    grams[i]->computeVals(seq, seqFL_Grams[i], len);
  }
}

FL_Gram<int> * NGram::gram(int ord) {
  assert(ord >= 0);
  return this->grams[ord];
}
  

NGram::~NGram() {
  for(int i = 0; i < _maxOrder; i++)
    if(grams[i])
      delete grams[i];
}

