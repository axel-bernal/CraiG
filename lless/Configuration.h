/****************************************************************************
* Configuration.h - part of the lless namespace, a general purpose
*                   linear semi-markov structure prediction library
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

#ifndef _CONFIGURATION_H_
#define _CONFIGURATION_H_
  
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Utils.h"
#include <vector>
#include "Resource.h"


namespace lless {
  
  /**
   * The Configuration class is a class for storing filenames of training
   * and validation data and other details.
   * The following is an example of a configuration file:\n
   * <TABLE>
   * <TR><TD>  ValidationSequences </TD><TD> fasta_for_validation</TD>
   * <TR><TD>  ValidationTags      </TD><TD> tags_for_validation</TD>
   * <TR><TD>  Sequences           </TD><TD> fasta_for_training</TD>
   * <TR><TD>  Tags                </TD><TD> tags_for_training</TD>
   * <TR><TD>  Path                </TD><TD> where_this_file_is_located</TD>
   * <TR><TD>  Name                </TD><TD> model_name</TD>
   * <TR><TD>  PrefixFiles         </TD><TD> where input files should be 
   *                                         accessed
   * </TABLE>\n
   * It basically contains all the information needed by the trainer to 
   * start the learning process. Sequences and ValidationSequences are
   * the file names of the input sequence files for training and validation
   * respectively. ValidationTags and Tags are Tag annotations of the 
   * input sequence.
   *
   ***************************************************************************/
  
  class Configuration {
   protected:
    // general
    std::string _name;
    std::string _domain;
    std::string _trainingSequences;
    std::string _validationSequences;
    std::string _trainingEntries;
    std::string _validationEntries;
    std::string _prefixFiles;
   public:
    Configuration(ifstream &fd)  {
      std::string propName;
      _prefixFiles = "";
      
      while(!fd.eof()) {
        fd >> propName;
        if(!propName.compare("Name")) 
          fd >> _name;
        else if(!propName.compare("Domain"))
          fd >> _domain;
        else if(!propName.compare("Tags"))
          fd >> this->_trainingEntries;
        else if(!propName.compare("ValidationTags"))
          fd >> this->_validationEntries;
        else if(!propName.compare("Sequences"))
          fd >> this->_trainingSequences;
        else if(!propName.compare("ValidationSequences"))
          fd >> this->_validationSequences;
	else if(!propName.compare("PrefixFiles"))
	  fd >> this->_prefixFiles;
      }
    }
  
    inline std::string & name() {  return _name; }
    inline std::string & domain() {  return _domain; }
    inline std::string & trainingSequences() {  return this->_trainingSequences; }
    inline std::string & validationSequences() {  return this->_validationSequences; }
    inline std::string & trainingSet() {  return this->_trainingEntries; }
    inline std::string & validationSet() {  return this->_validationEntries; }
    inline std::string & prefixFiles() {  return this->_prefixFiles; }
  };

}
  
#endif
  
  
  
  
  
