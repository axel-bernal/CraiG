/****************************************************************************
* ContextIMM.h - part of the lless namespace, a general purpose
*                linear semi-markov structure prediction library
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

#ifndef _CONTEXT_IMM_H_
#define _CONTEXT_IMM_H_

#include <fstream>
#include <vector>
#include "ResourceEngine.h"
#include "Utils.h"
#include "Sigma.h"
#include "IMM.h"



namespace lless {

  /**
   * ContextIMM is a subtype of Resource. The class contains 
   * different IMM models, each of them corresponding to a different context. 
   * The method for training the model are also provided and it's called when
   * the Resource Contents need preprocessing - in which case the Contents 
   * are a list of input sequences in FASTA format.
   *
   ***************************************************************************/

  class ContextIMM : public Resource {
   private:
    int numContexts;        //!< Number of contexts
    int numDesiredContexts; /*!< This number could be larger than the actual 
                              number of contexts. */
    Sigma *sigma;
    IMM **imm;
    // only for training
    int numPos;
    int minOrd, maxOrd;
    int sldWindow;
    std::string ctxString;

   public:
    //! Default constructor for testing the cimm
    ContextIMM(::ifstream &fd,   //!< ifstream containing the cimm definition.
               Sigma *sigma,     /*!< alphabet used by the input training
                                   sequences. */
               
               int contexts = -1 //!< number of desired contexts.
               ) :
      Resource(false) {
      this->sigma = sigma;
      this->numDesiredContexts = contexts;
      retrieveContents(fd);
    }

    //! Default constructor for training the cimm
    ContextIMM(const char *inputSeqs,  /*!< the filename containing the input
                                         sequences for training. */
               std::string &ctxString, /*!< the string which represents the
                                         the context levels used by this cimm.
                                       */
               int numPos,             /*!< the number of subchains for each
                                         imm submodel. */
               
               int minOrd,             //!< minimum order for each imm submodel
               int maxOrd,             //!< maximum order for each imm submodel
               int sldWindow,          /*!< the sliding window to compute 
                                         context levels */
               
               Sigma *sigma) :

      Resource(true) {

      this->numPos = numPos;
      this->minOrd = minOrd;
      this->maxOrd = maxOrd;
      this->ctxString = ctxString;
      this->sldWindow = sldWindow;
      this->sigma = sigma;
      std::ifstream fd(inputSeqs, ios::in);
      retrieveContents(fd);
    }


    /**
     * Constructor from a Header string definition.
     */   
    ContextIMM(
               /**
                * The Header string definition, loaded as a vector of strings.
                * The Header has the following form:\n\n
                * Resource id ContextIMM needsTrain fileName contexts sigma\n\n
                * The field needsTrain is equivalent to needsPreProcessing in
                * Resource; fileName is the name of file containing either the
                * input sequences for training (if needsTrain is set to true)
                * or the cimm definition (if needsTrain is set to false);
                * If fileName equals "null", then the cimm Contents are
                * to be read right after the Header.
                
                * contexts and sigma are defined as in the other
                * constructors.
                */
                std::vector<std::string> &params,
                int & offset,       //!<The index for vector params.
                ResourceEngine *re  //!<A pointer to the ResourceEngine object.
                ) 
      : Resource(params, offset) {

      _needsPreProcessing = Utils::stringToBoolean(params[offset++]);
      std::string & file = params[offset++];

      if(!sscanf(params[offset++].c_str(), "%d", &numDesiredContexts))
        assert(0);

      sigma = (Sigma *)re->getResource(params[offset++]);

      if(file.compare("null") == 0)
        _contentsInline = true; // model to be read inline @ ResourceEngine scope
      else {
        std::ifstream fd(file.c_str());
        assert(fd);
        retrieveContents(fd);
      }
    }
    
    /**
     * @see Resource::retrieveContents(std::ifstream &)
     */
    void retrieveContents(::ifstream & fd) {

      if(needsPreProcessing()) { 
        preProcess(fd);
        return;
      }

      loadModel(fd);
    }

    inline IMM **getIMMModels() {
      return imm;
    }
    
    /**
     * @see Resource::saveHeader(std::ofstream &)
     */
    inline void saveHeader(::ofstream &fd) {
      Resource::saveHeader(fd);
      fd << " false null " << numDesiredContexts << " " << sigma->getName() << endl;
    }

    /**
     * A member function which generates a sequence of characters using the 
     * the cimm emission probabilities
     * @param context the input context used.
     * @param seq the generated sequence.
     * @param frame the initial frame which should be used to generate the 
     * first character of sequence seq. This parameter is different from 1
     * when the imm submodel associated with context is inhomogeneous.
     * @param length the length of the generated sequence.
     */
    inline void emitSequence(int context, char *seq, int frame, int length) {
      imm[context]->emitSequence(seq, frame , length);
    }

    void saveContents(::ofstream &fd);
    void loadModel(::ifstream & fd);  
    void preProcess(std::ifstream & fd);

    ~ContextIMM();

  };

}
#endif
