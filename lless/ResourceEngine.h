/****************************************************************************
* ResourceEngine.h - part of the lless namespace, a general purpose
*                    linear semi-markov structure prediction library
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

#ifndef _RESOURCE_ENGINE_H_
#define _RESOURCE_ENGINE_H_

#include <ctype.h>
#include <assert.h>
#include <iostream>
#include <cfloat>
#include "ObjectFactory.h"
#include <list>  
#include <vector>
#include <boost/regex.hpp>
#include "Resource.h"


namespace lless {

  /**
   * The class FilterEngine read Resource Header/Content definitions and 
   * creates Resource objects dynamically.
   * In lless, resource definitions can be retrieved from an input  resource
   * file in order to create Resource objects dynamically. The ResourceEngine 
   * class takes care of all this by using an instance of a customized 
   * ObjectFactory which receives these definitions and creates objects of type
   * Resource.\n
   * Within the resource file, a Resource object is defined in at most two 
   * parts: the Header and possibly the Contents. Each line in a resource file
   * that starts with the keyword Resource is a Header. The format is the 
   * following:\n\n
   * Resource  resource-name   ResourceClass   space-separated-parameters \n\n
   * ResourceClass is either a Resource or any Class derived from Resource, 
   * whereas resource-name is a unique identifier.\n
   * Subsequent lines in the file possibly specify the Resource's Contents 
   * unless the Resource class member _contentsInline is false.\n
   * The reserved keyword ~ is an acronym for the model path, which is
   * equivalent to the absolute path to the directory containing the model(s)
   * concatenated with the model name. This acronym is mostly used to access 
   * resources whose Contents are too large or cannot be provided inline.
   *
   ***************************************************************************/

  class ResourceEngine {

   protected:
    vector<Resource *> resources;
    vector<string> resourceHeaders;
    std::map<std::string, Resource *> namedResources;
    ObjectFactory3< Resource, vector<std::string> &, int &, ResourceEngine *, std::string > resourceFactory;

   public:
    ResourceEngine() {
      registerDefaultResources();
    }

    ResourceEngine(std::ifstream & fd, 
		   map<std::string, std::string> *rsubs = NULL) {
      registerDefaultResources();
      retrieveResources(fd, rsubs);
    }

    /*    ResourceEngine(std::ifstream &fd, const char *modelPath)  {
      registerDefaultResources();
      readHeaderDefinitions(fd, modelPath);
      }*/

    void saveResources(std::ofstream &fd) {
      boost::RegEx rExPrefixFiles("PREFIX_EVIDENCE");
      for(unsigned int i = 0; i < resources.size(); i++) {
	// preserve the orig header if it contains PREFIX_EVIDENCE
	if(rExPrefixFiles.Search(resourceHeaders[i]))
	  fd << resourceHeaders[i] << endl;
	else 
	  resources[i]->saveHeader(fd);

	resources[i]->saveContents(fd);

      }
      fd << "//" << endl;
    }

    inline void retrieveResources(std::ifstream &fd,
				  map<string, string> *rsubs = NULL) {
      readHeaderDefinitions(fd, rsubs);
    }

    inline void setResource(Resource *r) {
      assert(r != NULL);
      resources.push_back(r);
      namedResources[r->getName()] = r;
    }

    inline Resource * getResource(std::string resourceName) { 
      std::map<std::string, Resource *>::iterator it;
      it = namedResources.find(resourceName);

      if(it == namedResources.end()) {
        assert(0);
        throw EXCEPTION(PARSE_ERROR, std::string("undefined Resource ") + resourceName);
      }
      return it->second;    
    }  

    void registerDefaultResources();
    void readHeaderDefinitions(std::ifstream &reStream,
			       map<std::string, std::string> *rsubs = NULL);
    
    void deleteResourceSeqContents() {
      for(unsigned int i = 0; i < resources.size(); i++)
	resources[i]->releaseSeqContents();
    }

    ~ResourceEngine() {  
      for(unsigned int i = 0; i < resources.size(); i++)
        delete resources[i];
    }

  };

}

#endif
