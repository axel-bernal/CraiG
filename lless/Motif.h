/****************************************************************************
* Motif.h - part of the lless namespace, a general purpose
*           linear semi-markov structure prediction library
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

#ifndef _DB_SOFT_MOTIF_
#define _DB_SOFT_MOTIF_

#include "ResourceEngine.h"
#include "Sigma.h"
#include <sstream>


namespace lless {

  /**
   * The Motif class is a subtype of Resource. It is an abstract base class
   * which could be used to represent a motif, a signal or a pattern through
   * a motif profile. Derived classes may opt to read the motif profile
   * directly from a file or read sequence data first and then compute the
   * motif profile afterwards. The first scenario would be described as a
   * Resource whose Contents do not need preprocessing and the second scenario
   * is best described as a Resource which need preprocessing. In this case, 
   * the preprocessing is needed to compute the motif profile.
   *
   ***************************************************************************/

  class Motif : public Resource {
   protected:
    Sigma *sigma;

   public:
    /**
     * Minimal constructor
     * @param sigma the alphabet for the motif
     * @param needsPreProcess true if motif profile needs to be computed,
     * false otherwise.
     */
    Motif(
          Sigma *sigma,
          bool needsPreProcess) 
      : Resource(needsPreProcess) {

      this->sigma = sigma;

    }

    /**
     * Constructor from a Header string definition.
     * @param params The Header string definition, loaded as a vector of 
     * strings.
     * @param offset The index for vector params
     * @param re A pointer to the ResourceEngine object.
     */
    Motif(
          std::vector<std::string> &params,
          int & offset, 
          ResourceEngine *re) 
      : Resource(params, offset) {
      
      sigma = (Sigma *)re->getResource(params[offset++]);
    }

    virtual inline int motifLength() = 0;
    virtual double score(double threshold, char *word) = 0;
    
    virtual ~Motif() {
    }

  };

  /**
   * TransfacMotif is a subclass of Motif in which the motif profile is
   * implemented as Transfac Motif definition that is read from an input file.
   * Since the Transfac definition is already a motif profile, there is no 
   * need to preprocess it.
   * 
   */
  class TransfacMotif : public Motif {

   protected:
    vector<double *> mat;
    std::string transfacMatStr;

   public:
    /**
     * Minimal constructor
     * @param fd an input file stream containing the Transfac motif 
     * definition. 
     * @param sigma a motif alphabet
     */
    TransfacMotif(
             std::ifstream & fd,
             Sigma *sigma
             ) : Motif(sigma, false) {
      
      retrieveContents(fd);
    }

    /**
     * Constructor from a Header string definition.
     * @param params The Header string definition, loaded as a vector of 
     * strings. The Header has the following form:\n\n
     * Resource name  TransfacMotif sigma transfacFile\n\n
     * If transfacFile equals "null", then the motif Contents are to be read
     * right after the Header.
     * The description of the fields could be found in 
     * the other constructor(s)     
     * @param offset The index for vector params
     * @param re A pointer to the ResourceEngine object.
     */
    TransfacMotif(
             std::vector<std::string> &params,
             int & offset, 
             ResourceEngine *re) 
      : Motif(params, offset, re) {
      
      _needsPreProcessing = true;

      std::string file = params[offset++];
      
      if(file.compare("null") == 0)
        _contentsInline = true; // model to be read inline @ 
                                // ResourceEngine scope
      else {
        std::ifstream fd(file.c_str());
        assert(fd);
        retrieveContents(fd);
      }
    }

    inline int motifLength() {
      return (int)mat.size();
    }    

    /**
     * @return the motif score if it is greater than threshold, 0 otherwise. 
     * The score is computed from the probabilities read from the Transfac
     * definition file.
     */
    double score(double threshold, char *word) {
      double result = 0;
      if(motifLength() > (int)strlen(word)) 
        return result;
      
      for(int i = 0; i < motifLength(); i++) 
        result += mat[i][sigma->indC(word[i])];

      if(result < threshold)
        return 0;
      return result;
    }
    
    /**
     * @see Resource::retrieveContents(std::ifstream &)
     */
    void retrieveContents(::ifstream & fd) {

      std::string line;
      int i = 0, j = 0, k;
      int *header = new int [sigma->alphabetSize()];
      vector<double> totals;
      char val[100];

      /*
       * Initializing header
       */
      for(k = 0; k < sigma->alphabetSize(); k++)
        header[k] = 0;

      /*
       * Get the first line of the matrix
       */
      getline(fd, line);
      transfacMatStr += line + std::string("\n");
      
      while(!fd.eof() && line.substr(0, 2).compare("PO") != 0) {
        getline(fd, line);
        transfacMatStr += line + std::string("\n");
      }
      /*
       * Processing header
       */
      char *ptr;
      ptr = (char *)line.c_str() + 2;
      while(sscanf(ptr, "%s", val) == 1 && j < sigma->alphabetSize()) {
        header[j++] = sigma->indC(val[0]);
        ptr = strchr(ptr, val[0]) + 1;
      }
      
      /*
       *  Processing matrix values
       */
      getline(fd, line);
      transfacMatStr += line + std::string("\n");

      while(!fd.eof() && line.substr(0, 2).compare("XX") != 0) {
        j = 0;
        ptr = (char *)line.c_str();
        double *row = new double [sigma->alphabetSize()];
        mat.push_back(row);
        totals.push_back(0);

        while(sscanf(ptr, "%s", val) == 1) {
          if(j == 0) {
            assert(atoi(val) == i + 1);
          }
          else if(j == sigma->alphabetSize() + 1)
            ;
          else {
            mat[i][header[j - 1]] = atof(val);
            totals[i] += mat[i][header[j - 1]];
          }
          j++;
          ptr = strstr(ptr, val) + strlen(val);
        }
        i++;
        getline(fd, line);
        transfacMatStr += line + std::string("\n");
      }
      while(!fd.eof() && line.substr(0, 2).compare("//") != 0) {
        getline(fd, line);
        transfacMatStr += line + std::string("\n");
      }
      delete [] header;
      
      for(i = 0; (unsigned)i < totals.size(); i++)
        for(j = 0; j < sigma->alphabetSize(); j++) 
          mat[i][j] = mat[i][j]/totals[i];
    }

    inline void saveHeader(::ofstream &fd) {
      Resource::saveHeader(fd);
      fd << " " << sigma->getName() << " null" << endl;
    }

    inline void saveContents(std::ofstream & fd) {
      fd << transfacMatStr;
    }
    
    ~TransfacMotif() {
      for(unsigned int i = 0; i < mat.size(); i++) {
        delete mat[i];
        mat[i] = NULL;
      }
      mat.clear();
    }
    
  };

  /**
   * The class DBSoftMatchMotif builds a motif definition based on soft matches
   * against a database composed of multiple sequences, each of which is a
   * motif occurrence.
   * A soft match is a motif match with at most maxHammingDist mismatches
   * in any motif position. This Motif class cannot work for maxHammingDist 
   * greater than 1, as the memory grows combinatorialy with the number of 
   * allowed mismatches.
   */
  class DBSoftMatchMotif : public Motif {

   protected:
    int **database; //!< the database containing the different motif strings
    int column;
    int len;
    int maxHammingDist;
    int numWords[2];

   public:
    /**
     * Default constructor. 
     * @param fd the input file stream containing the database of motif
     * occurrences
     * @param column the column in file fd at which the motif starts
     * @param len the length of the motif
     * @param sigma the motif's alphabet
     */
    DBSoftMatchMotif(
                     std::ifstream & fd,
                     int column, 
                     int len, 
                     Sigma *sigma
                     ) : Motif(sigma, true) {
      
      this->maxHammingDist = 1; 
      this->column = column;
      this->len = len;
      initializeDB();
      retrieveContents(fd);
    }
    

    /**
     * Constructor from a Header string definition.
     * @param params The Header string definition, loaded as a vector of 
     * strings. The Header has the following form:\n\n
     * Resource name DBSoftMatchMotif sigma needsPreProcessing contentsFile
     * column len \n\n
     * If contentsFile equals "null", then the motif Contents are to be read
     * right after the Header.
     * The description of the fields could be found in 
     * the other constructor(s)
     * @param offset The index for vector params
     * @param re A pointer to the ResourceEngine object.
     */
    DBSoftMatchMotif(
                     std::vector<std::string> &params,
                     int & offset, 
                     ResourceEngine *re) 
      : Motif(params, offset, re) {

      _needsPreProcessing = Utils::stringToBoolean(params[offset++]);
      std::string file = params[offset++];

      if(!sscanf(params[offset++].c_str(), "%d", &column))
        assert(0);

      if(!sscanf(params[offset++].c_str(), "%d", &len))
        assert(0);

      maxHammingDist = 1;

      initializeDB();

      if(file.compare("null") == 0)
        _contentsInline = true; // model to be read inline @ 
                                // ResourceEngine scope
      else {
        std::ifstream fd(file.c_str());
        assert(fd);
        retrieveContents(fd);
      }
    }

    /**
     * A member function that allocates memory and initializes the table 
     * that will contain all motif occurrences found in the input file
     */
    void initializeDB() {
      int i, j;
      int alphSize = sigma->alphabetSize();

      database = new int *[ maxHammingDist + 1];

      for(i = 0; i < maxHammingDist + 1; i++) {
        database[i] = new int [Utils::pow2(alphSize, len)];

        for(j = 0; j < Utils::pow2(alphSize, len); j++)
	  database[i][j] = 0;
      }
    }

    inline int maxHammingDistance() {
      return maxHammingDist;
    }

    inline int motifLength() {
      return len;
    }    

    /**
     * @return the motif score if it is greater than threshold, 0 otherwise. 
     * The score is computed as the relative frequency of the word's soft 
     * matches against the database of motif occurrences.
     */
    inline double score(double threshold, char *word) {
      if(len > (int)strlen(word))
        return 0;
      
      int ind = sigma->randIndS(word, len);
      double result = ((double)softFreq(0, ind))/numWords[0];

      if(result < threshold)
        return 0;

      return result;
    }

    /**
     * @return the absolute frequency of soft matches of word against the
     * database of motif occurrences.
     */
    inline int softFreq(int hammingDistance, int word) {
      assert(hammingDistance <= maxHammingDist);
      return database[hammingDistance][word];
    }

    inline Sigma * alphabet() {
      return sigma;
    }

    /**
     * @see Resource::retrieveContents(std::ifstream &)
     */    
    void retrieveContents(::ifstream & fd) {

      assert(fd);
      if(needsPreProcessing()) { 
        preProcess(fd);
        return;
      }

      loadDatabase(fd);
    }

    /**
     * A member function that loads the database of motif occurrences into
     * memory from input file stream fd.
     */
    void loadDatabase(std::ifstream & fd) {
      std::string dummy;
      fd >> dummy >> len;
      int i, j, nonZero;

      for(i = 0; i <= maxHammingDist; i++) {
        fd >> dummy >> dummy >> numWords[i];
        nonZero = numWords[i];

        while(nonZero) {
          fd >> j;
	  fd >> database[i][j];
          nonZero--;
        }
      }    
      fd >> dummy;
      assert(dummy.compare("//") == 0);
    }

    /**
     * A member function which preprocess the database of motif occurrences
     * which is in input file stream fd.
     * The resulting table contains all possible words of the same length
     * as the motif and their respective number of soft matches in the 
     * database
     */
    void preProcess(std::ifstream & fd) {
      int i, j, baseCode;
      std::string line;
      int alphSize = sigma->alphabetSize();
      numWords[0] = 0;
      numWords[1] = 0;

      while(!fd.eof()) {
        std::getline(fd, line);

        if(!strlen(line.c_str())) 
          continue;
        
        /*
         * hamming distance  = 0
         */
        char *sigSeq = (char *)line.c_str() + column;
        baseCode = sigma->randIndS(sigSeq, len);
        database[0][baseCode]++;
        numWords[0]++;

        /*
         * hamming distance = 1
         */
        int *s = new int[len + 1];
        assert(column + len <= (int)strlen(line.c_str()));

        for(i = 1; i <= len; i++)
          s[i]  =  sigma->randIndS(sigSeq, i - 1, i -1);

        int prefix;

        for(i = 0; i < len; i++) {
          for(j = 0; j < alphSize; j++)
            /*
             * build baseCode for sigSeq[0..i-1],j,sigSeq[i+1..len]
             */
            if(j != s[i+1]) { 
              prefix = sigma->randIndS(sigSeq, i);
              baseCode = sigma->concat(prefix*alphSize+j, sigSeq+i+1, len-i-1);
              if(!database[1][baseCode])
                numWords[1]++;
	      
              database[1][baseCode]++;
              //            sigma->contS(seq, baseCode, len);
              //            cerr << "\t" << seq << endl;
            }
        }

        delete [] s;
      }

      fd.close();    
    }

    /**
     * @see Resource::saveHeader(::ofstream &fd)
     */
    inline void saveHeader(::ofstream &fd) {
      Resource::saveHeader(fd);
      fd << " " << sigma->getName() << " false null " << column << " " << len << " " << endl;
    }

    /**
     * @see Resource::saveContents(::ofstream &fd)
     */
    inline void saveContents(std::ofstream & fd) {
      fd << "len= " << len << endl;
      int alphSize = sigma->alphabetSize();

      for(int i = 0; i <= maxHammingDist; i++) {
        int nonZero = 0;
	int numW = (i == 0 ? 1 : len + 1);

        for(int j = 0; j < Utils::pow2(alphSize, len); j++)
          if(database[i][j])
            nonZero++;

        fd << "HD= " << i << " " << nonZero << endl;

        for(int j = 0; j < Utils::pow2(alphSize, len); j++) {
          if(database[i][j]) 
	    fd << j << " " << database[i][j] << endl;
	}
      }
      fd << "//" << endl;
    }

    ~DBSoftMatchMotif() {
      int alphSize = sigma->alphabetSize();
      for(int i = 0; i <= maxHammingDist; i++) {
        if(this->database[i]) 
	  delete [] this->database[i];
	this->database[i] = NULL;
      }
      database = NULL;    
    }

  };

}

#endif
