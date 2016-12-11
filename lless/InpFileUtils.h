#ifndef _INPFILE_UTILS_H_
#define _INPFILE_UTILS_H_

#include "InpFile.h"

/****************************************************************************
* InpFile.cpp - part of the lless namespace, a general purpose
*               linear semi-markov structure prediction library
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

  static std::map<std::string, TFileFormat> fileTypes;
  static std::map<TFileFormat, std::string> fileNames;

  class InpFileUtils {
   public:
    InpFileUtils() {
      // boolean values
      fileTypes["FASTA"] = FASTA;
      fileTypes["SEQTAGS"] = SEQTAGS;
      fileTypes["LDSCORE"] = LDSCORE;
      fileTypes["XFASTA"] = XFASTA;
      fileTypes["MULTIZ"] = MULTIZ;
      fileTypes["EDGEANNOT_FASTA"] = EDGEANNOT_FASTA;
      fileTypes["EXTENDED_LDSCORE"] = EXTENDED_LDSCORE;
      fileTypes["EXTENDED_MULTI_INTSCORE"] = EXTENDED_MULTI_INTSCORE;
      fileTypes["EXTENDED_MULTI_LDSCORE"] = EXTENDED_MULTI_LDSCORE;

      fileNames[FASTA] = "FASTA";
      fileNames[SEQTAGS] = "SEQTAGS";
      fileNames[LDSCORE] = "LDSCORE";
      fileNames[XFASTA] = "XFASTA";
      fileNames[MULTIZ] = "MULTIZ";
      fileNames[EDGEANNOT_FASTA] = "EDGEANNOT_FASTA";
      fileNames[EXTENDED_LDSCORE] = "EXTENDED_LDSCORE";
      fileNames[EXTENDED_MULTI_INTSCORE] = "EXTENDED_MULTI_INTSCORE";
      fileNames[EXTENDED_MULTI_LDSCORE] = "EXTENDED_MULTI_LDSCORE";
    }
    
    static TFileFormat stringToTFileFormat(const std::string &s) {
      std::map<std::string, TFileFormat>::iterator it;
      it = fileTypes.find(s);
      if(it == fileTypes.end()) {
        assert(0);
        throw EXCEPTION( PARSE_ERROR, std::string("undefined TFileFormat: ")+s);
      }
      return it->second;
    }
    
    static std::string tInpFileToString(const TFileFormat fType) {
      std::map<TFileFormat, std::string>::iterator it;
      it = fileNames.find(fType);
      if(it == fileNames.end()) {
        assert(0);
        throw EXCEPTION( PARSE_ERROR, std::string("undefined TFileFormat name "));
      }
      return it->second;
    }
  };
  
  InpFileUtils _tmpfileobj;

}

#endif
