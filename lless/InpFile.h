/****************************************************************************
* InpFile.h - part of the lless namespace, a general purpose
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

#ifndef _INP_FILE_H_
#define _INP_FILE_H_
#include <map>
#include "Sequence.h"

namespace lless {

  typedef enum {FASTA, SEQTAGS, XFASTA, MULTIZ, EDGEANNOT_FASTA, EXTENDED_INTSCORE, LDSCORE, EXTENDED_LDSCORE, EXTENDED_MULTI_INTSCORE, EXTENDED_MULTI_LDSCORE} TFileFormat;
 
  class IOSource : public Resource {
   protected:
    TFileFormat _type;
    std::string filename;
    bool mapIds;
    bool loadAllSeqs;
    Sigma *sigma;  
    std::map<std::string, BasicSeq *> seqs;
    list<BasicSeq *> _seqs;
    std::vector<std::string> _seqNames;    
    int _numCollapsedVals;
    FSM *fsm;
    bool _lockContents;
   public:
    
    /**
     * Constructor 
     * @param name A unique resource name
     * @param type File type
     * for more details
     * @param mapIds An auxiliary hash table with the sequence ids as
     * keys is built if this parameter is true. This implies that 
     * unique identifiers must be supplied for the sequences
     * @param loadAllSeqs load the file completely into memory. Only
     * recommended if file is small enough
     */
    IOSource(std::string name, 
             TFileFormat type,
	     Sigma *sigma = NULL,
             bool mapIds = false,
             bool loadAllSeqs = true,
	     FSM *fsm = NULL) : Resource (name) {
      
      initialize("null", type, mapIds, sigma, loadAllSeqs, fsm);
      
    }

    /**
     * Default constructor
     * @param name A unique resource name
     * @param filename The name of the file containing the sequences in the
     * format required
     * @param type File type
     * @param sigma An alphabet for the contained sequences in the input file.
     * Its value must be initialized to NULL if the sequence is made of 
     * floating-point numbers
     * @param mapIds An auxiliary hash table with the sequence ids as
     * keys is built if this parameter is true
     * @param loadAllSeqs load the file completely into memory. Only
     * recommended if file is small enough
     */
    IOSource(std::string name, 
             std::string filename,
             TFileFormat type, 
             Sigma *sigma = NULL,
             bool mapIds = false,
             bool loadAllSeqs = true,
	     FSM *fsm = NULL) : Resource (name) {

      initialize(filename, type, mapIds, sigma, loadAllSeqs, fsm);
      
    }
    
    /**
     * Default constructor
     * @param name A unique resource name
     * @param filestream The file stream containing the sequences in the
     * format required
     * @param type File type
     * @param sigma An alphabet for the contained sequences in the input file.
     * Its value must be initialized to NULL if the sequence is made of 
     * floating-point numbers
     * @param mapIds An auxiliary hash table with the sequence ids as
     * keys is built if this parameter is true
     * @param loadAllSeqs load the file completely into memory. Only 
     * recommended if file is small enoughif file is small enough
     */
    IOSource(std::string name, 
             std::ifstream &filestream,
             TFileFormat type,
             Sigma *sigma = NULL,
             bool mapIds = false,
             bool loadAllSeqs = true,
	     FSM *fsm = NULL) : Resource (name) {
      
      initialize("null", type, mapIds, sigma, loadAllSeqs, fsm);

    }

    /**
     * Auxiliary Header Constructor 
     * @param name A unique resource name
     * @param type File type
     * for more details
     * @param mapIds An auxiliary hash table with the sequence ids as
     * keys is built if this parameter is true
     * @param loadAllSeqs load the file completely into memory. Only
     * recommended if file is small enough
     */
    IOSource(std::vector<std::string> &params,
             int & offset,
             TFileFormat type,
             bool mapIds,
             bool loadAllSeqs,
	     FSM *fsm = NULL) : Resource (params, offset) {
      
      initialize("null", type, mapIds, NULL, loadAllSeqs, fsm);

    }

    /**
     * Constructor for the Collapser/Merger classes
     */
    IOSource(std::vector<std::string> &params,
             int & offset,
	     ResourceEngine *re,
             TFileFormat type,
             bool mapIds) : Resource (params, offset) {
      
      std::string sigmaId = params[offset++];
      Sigma *sigma = NULL;

      if(sigmaId.compare("null"))
	sigma  = (Sigma *)re->getResource(sigmaId);

      bool loadAllSeqs = Utils::stringToBoolean(params[offset++]);
      initialize("null", type, mapIds, sigma, loadAllSeqs);

    }

    IOSource(std::vector<std::string> &params,
             int & offset,
             ResourceEngine *re);
    
    inline void initialize(std::string filename, TFileFormat type, 
                           bool mapIds, Sigma *sigma,
                           bool loadAllSeqs,
			   FSM *fsm = NULL) {
      
      this->filename = filename;
      this->_type = type;
      this->mapIds = mapIds;
      this->sigma = sigma;
      this->loadAllSeqs = loadAllSeqs;
      this->_numCollapsedVals = 1;
      this->fsm = fsm;
      this->_lockContents = false;
    }
        
    // sequences are being edited and they won't be released
    // but reading new sequences is still allowed
    void lockContents() {
      this->_lockContents = true;
    }

    void unlockContents() {
      this->_lockContents = false;
    }
    
    int numCollapsedVals() {
      return _numCollapsedVals;
    }

    virtual std::string & collapsedfNames(int permNo = 0) {
      return filename;
    }

    inline Sigma *alphabet() {
      return sigma;
    }

    virtual inline int alphabetSize() {
      if(sigma)
	return sigma->alphabetSize();
      return 1;
    }

    virtual inline void print(BasicSeq *) {
      ;
    }

    inline TFileFormat type() {
      return _type;
    }

    inline bool isFasta() {
      return (_type == FASTA || _type == XFASTA ||
	      _type == MULTIZ || _type == EDGEANNOT_FASTA);
    }

    inline bool isMultiScore() {
      return _type == EXTENDED_MULTI_INTSCORE || _type == EXTENDED_MULTI_LDSCORE;
    }

    inline std::vector<std::string> & seqNames() {
      if(!mapIds)
        throw EXCEPTION(BAD_USAGE, filename + "'s ids need to be mapped first\n");
      
      return _seqNames;
    }

    virtual inline std::string & fileName() {
      return filename;
    }
    
    /**
     * @return true if fasta file has been loaded in memory, false otherwise
     */
    inline bool loaded() {
      return (seqs.size() != 0 || !loadAllSeqs);
    }
    
    /**
     * @return the list of BasicSeq objects 
     */
    inline list<BasicSeq *> & sequences() {
      return _seqs;
    }
    
    BasicSeq* loadSequence(std::ifstream &annotSeqStream, 
			   std::string &line);
    
    BasicSeq* findSeq(std::string id);
    BasicSeq* cacheSeq(BasicSeq *c);   
    
    // virtual functions

    virtual void releaseSeqContents();
    virtual void saveHeader(std::ofstream &fd);
    virtual BasicSeq *loadFileSequence(std::string id) = 0;
    virtual void loadFileSequences() {  }
    virtual void mapInputSeqs() = 0;
    
    virtual ~IOSource() {
      list<BasicSeq *>::iterator it = _seqs.begin(); 
      for( ; it != _seqs.end(); it++)
        delete (*it);

      _seqs.clear();
      seqs.clear();
    }
    
  };

  /**
   * The InpFile class is a subtype of Resource, which allows the user to 
   * upload a file and keep it in memory. There are three types of files 
   * that can be uploaded:  FASTA, XFASTA and EXTENDED_LDSCORE. 
   * If the FASTA type is  specified, then the class expects every entry 
   * in the file to have the following format:
   * >header
   * sequence_of_characters
   * ..
   * When XFASTA is specified the format file should be:
   * >header
   * symbol1 x times_symbol1_appears_in_sequence
   * symbol2 x times_symbol2_appears_in_sequence
   * ..
   *
   * It is in fact a very basic compression and it would prove to be very 
   * useful when uploading fasta sequences which contain mostly repeated
   * characters.
   * The third type, EXTENDED_LDSCORE, is specified for files with the
   * following format:
   * >header
   * pos double_val
   * ..
   * The EXTENDED_INTSCORE type is exactly the same as EXTENDED_LDSCORE,
   * except each entry point to an integer argument.
   * For the score-type formats, supplying an alphabet is not necessary, i.e.
   * the member sigma in this case is assumed to be NULL.
   *
   * All these files can be in a per-sequence id basis within a directory,
   * in which case, the variable is_dir will be set to true
   *
   ***************************************************************************/

  class InpFile : public IOSource {
   private:
    std::ifstream _fs;
   protected:
    std::map<std::string, streampos> seqIndexes;
    std::ifstream *fs;
    
   public:

    /**
     * Constructor 
     * @param name A unique resource name
     * @param type File type
     * for more details
     * @param mapIds An auxiliary hash table with the sequence ids as
     * keys is built if this parameter is true. This implies that 
     * unique identifiers must be supplied for the sequences
     * @param loadAllSeqs load the file completely into memory. Only
     * recommended if file is small enough
     */
    InpFile(std::string name, 
            TFileFormat type, 
            bool mapIds = false,
            bool loadAllSeqs = true,
	    FSM *fsm = NULL) : 
    IOSource(name, type, NULL, mapIds, loadAllSeqs, fsm) {

    }
    
    /**
     * Default constructor
     * @param name A unique resource name
     * @param filename The name of the file containing the sequences in the
     * format required
     * @param type File type
     * @param sigma An alphabet for the contained sequences in the input file.
     * Its value must be initialized to NULL if the sequence is made of 
     * floating-point numbers
     * @param mapIds An auxiliary hash table with the sequence ids as
     * keys is built if this parameter is true
     * @param loadAllSeqs load the file completely into memory. Only
     * recommended if file is small enough
     */
    InpFile(std::string name, 
            std::string filename,
            TFileFormat type, 
            Sigma *sigma = NULL,
            bool mapIds = false,
            bool loadAllSeqs = true,
	    FSM *fsm = NULL) : 
    IOSource(name, filename, type, sigma, 
	     mapIds, loadAllSeqs, fsm) {
      
      _fs.open(filename.c_str(), ios::in);
      this->fs = &_fs;
      
      if(!fs->is_open())
        throw EXCEPTION(FILE_UNAVAILABLE, filename);
      
      loadFileSequences();
      mapInputSeqs();
      
    }
    
    /**
     * Default constructor
     * @param name A unique resource name
     * @param filestream The file stream containing the sequences in the
     * format required
     * @param type File type
     * @param sigma An alphabet for the contained sequences in the input file.
     * Its value must be initialized to NULL if the sequence is made of 
     * floating-point numbers
     * @param mapIds An auxiliary hash table with the sequence ids as
     * keys is built if this parameter is true
     * @param loadAllSeqs load the file completely into memory. Only 
     * recommended if file is small enoughif file is small enough
     */
    InpFile(std::string name, 
            std::ifstream &filestream,
            TFileFormat type,
            Sigma *sigma = NULL,
            bool mapIds = false,
            bool loadAllSeqs = true,
	    FSM *fsm = NULL) : 
    IOSource(name, filestream, type, 
	     sigma, mapIds, loadAllSeqs, fsm) {
      
      this->fs = &filestream;
      loadFileSequences();
      mapInputSeqs();

    }

    InpFile(std::vector<std::string> &params,
            int & offset,
            ResourceEngine *re) : 
    IOSource(params, offset, re) {
      
      this->_fs.open(filename.c_str(), ios::in);
      this->fs = &_fs;
      
      if(!fs->is_open())
        throw EXCEPTION(FILE_UNAVAILABLE, filename);
      
      loadFileSequences();
      mapInputSeqs();

    }

    BasicSeq* loadFileSequence(std::string id);
    void loadFileSequences();
    void mapInputSeqs();
    
    ~InpFile() {
      
    }
    
  };

  class InpDir : public IOSource {
   public:

    /**
     * Constructor 
     * @param name A unique resource name
     * @param type File type
     * for more details
     * @param mapIds An auxiliary hash table with the sequence ids as
     * keys is built if this parameter is true. This implies that 
     * unique identifiers must be supplied for the sequences
     * @param loadAllSeqs load the file completely into memory. Only
     * recommended if file is small enough
     */
    InpDir(std::string name, 
           TFileFormat type)
      : IOSource(name, type, NULL, true, false) {
      
    }

    /**
     * Default constructor
     * @param name A unique resource name
     * @param filename The name of the file containing the sequences in the
     * format required
     * @param type File type
     * @param sigma An alphabet for the contained sequences in the input file.
     * Its value must be initialized to NULL if the sequence is made of 
     * floating-point numbers
     * @param mapIds An auxiliary hash table with the sequence ids as
     * keys is built if this parameter is true
     * @param loadAllSeqs load the file completely into memory. Only
     * recommended if file is small enough
     */
    InpDir(std::string name, 
           std::string filename,
           TFileFormat type, 
           Sigma *sigma = NULL) :
      IOSource(name, filename, type, sigma, 
               true, false) {

      mapInputSeqs();      

    }

    /**
     * Auxiliary Header Constructor 
     * @param name A unique resource name
     * @param type File type
     * for more details
     * @param mapIds An auxiliary hash table with the sequence ids as
     * keys is built if this parameter is true
     * @param loadAllSeqs load the file completely into memory. Only
     * recommended if file is small enough
     */
    InpDir(std::vector<std::string> &params,
           int & offset,
           TFileFormat type) :
      IOSource(params, offset, type,
               true, false) {
      
    }
    
    InpDir(std::vector<std::string> &params,
           int & offset,
           ResourceEngine *re) : 
      IOSource(params, offset, re) {
      
      mapInputSeqs();
      
    }
    
    BasicSeq* loadFileSequence(std::string id);
    void mapInputSeqs();
    
    ~InpDir() {
      
    }
    
  };

  /**
   * The InpFileMerger class is a subtype of InpFile, which allows the user
   * to merge multiple InpFile type objects which contain sequences that are
   * defined through the same alphabet, i.e.the merge process will not work 
   * on files containing score sequences, only on objects of type Sequence. 
   * Also, all merging files must contain sequence based on the same alphabet. 
   *
   ***************************************************************************/

  class InpFileMerger : public IOSource {
   private:
    vector<InpFile *> files2Merge;
    char fillChar;
    bool overrideSigma;

   public:
    /**
     * Default constructor
     * @param name A unique resource name
     * @param fillChar character which identifies regions in sequences with
     * no useful information. These characters are the only ones that could be
     * yanked in a merge.
     */
    InpFileMerger(std::string name, 
		  bool loadAllSeqs,
                  char fillChar, 
		  Sigma *sigma = NULL)
      : IOSource(name, FASTA, sigma, true, loadAllSeqs) {
      
      this->fillChar = fillChar;
      this->overrideSigma = (this->sigma ? true : false);
      
      loadFileSequences();
      mapInputSeqs();

    }

    InpFileMerger(std::vector<std::string> &params,
		  int & offset,
		  ResourceEngine *re);

    void releaseSeqContents();
    void saveHeader(std::ofstream &fd);    
    BasicSeq* loadFileSequence(std::string id);
    void loadFileSequences();
    void mapInputSeqs();

    ~InpFileMerger() {}
  };


  /**
   * The InpFileCollapser class is a subtype of InpFile, which allows the 
   * user to collapse multiple InpFile type objects which contain objects 
   * of type Sequence
   *
   ***************************************************************************/
  template<class TClass1, class TClass2>
  class InpFileCollapser : public IOSource {
   protected:
    vector<InpFile *> files2Collapse;
    int _alphabetSize;
    int _gramLength;
    vector<int> *_indexes;
    vector<std::string> _collapsedfNames;
    bool overrideSigma;

   public:
    /**
     * Default constructor
     * @param name A unique resource name
     */
    InpFileCollapser(std::string name, 
		     int gramLength, 
		     bool loadAllSeqs,
		     Sigma *sigma = NULL) :
    IOSource(name, EXTENDED_MULTI_INTSCORE,
	     sigma, true, loadAllSeqs) {
      
      _gramLength = gramLength;
      _indexes = NULL;
      _alphabetSize = 1;
      this->overrideSigma = (this->sigma ? true : false);

    }
    
    /**
     * Constructor from a Header string definition.
     * @param params The Header string definition, loaded as a vector of 
     * strings.
     * @param offset The index for vector params
     * @param re A pointer to the ResourceEngine object.
     */    
    InpFileCollapser(std::vector<std::string> &params,
                     int & offset,
                     ResourceEngine *re) :
    IOSource(params, offset, re,
	     EXTENDED_MULTI_INTSCORE, true) {
      
      int numInpFiles;
      _indexes = NULL;

      if(!sscanf(params[offset++].c_str(), "%d", &this->_gramLength))
        assert(0);
      
      if(!sscanf(params[offset++].c_str(), "%d", &numInpFiles))
        assert(0);
      
      for(int i = 0; i < numInpFiles; i++)
        files2Collapse.push_back((InpFile *)re->getResource(params[offset++]));

      assert(_gramLength <= files2Collapse.size());

      computeNumCollapsedVals();
      
      _alphabetSize = 1;
      this->overrideSigma = (this->sigma ? true : false);
      
    }
	  
    inline std::string & collapsedfNames(int permutation = 0) {
      return _collapsedfNames[permutation];
    }

    inline int alphabetSize() {
      return _alphabetSize;
    }

    void computeAlphabetSize(int numSymbols) {
      for(int i = 0; i < _gramLength; i++)
        _alphabetSize = _alphabetSize*numSymbols;
    } 

    void releaseSeqContents() {
      // freeing contained files
      for(int i = 0; i < files2Collapse.size(); i++)
	files2Collapse[i]->releaseSeqContents();
      
      IOSource::releaseSeqContents();
    }

    /**
     * @see Resource::saveHeader(::ofstream &)
     */    
    virtual void saveHeader(std::ofstream &fd) {
      Resource::saveHeader(fd);
      
      fd << " " << (this->overrideSigma ? this->sigma->getName() : "null")
	 << (this->loadAllSeqs ? " true " : " false ")
	 << _gramLength << " " << files2Collapse.size();

      for(int i = 0; i < files2Collapse.size(); i++) 
        fd  << " " << files2Collapse[i]->getName();
    }

    /**
     * A member function that collapses all id sequences from files into
     * *this. At each position within this' sequence a value of
     * 1 (or 0) is shifted depending on weather any character contained 
     * within chars4BitOn appears in file's sequence. This process is 
     * carried for each file sequence.
     */
    BasicSeq *loadFileSequence(std::string id) {
      TStrand strand;
      int fileNo, i, t;
      TClass1 *target = NULL;
      
      for(fileNo = 0; fileNo < files2Collapse.size(); fileNo++) {
	Sequence *c = (Sequence *)files2Collapse[fileNo]->findSeq(id);
	vector<int> &fidxs = _indexes[fileNo];

	if(!files2Collapse[fileNo]->alphabet()) {
	  assert(0);
	  throw EXCEPTION( BAD_USAGE, files2Collapse[fileNo]->fileName() + " needs alphabet");
	}
	
	if(!overrideSigma)
	  this->sigma = files2Collapse[fileNo]->alphabet();
	
	if(!c) { 
	  if(target)  delete target;
	  return NULL;
	}
	
	if(!target) {
	  target = new TClass1(c->id(), _numCollapsedVals);
	  vector<TClass2 *> thisSeq(_numCollapsedVals);
	  
	  for(i = 0; i < _numCollapsedVals; i++)
	    thisSeq[i] = new TClass2(c->id(), c->length());
	  
	  for(t = 0; t < NUM_STRANDS; t++) {
	    strand = (TStrand)t;
	    
	    char *fileSeq = c->getSeq(strand);
	    
	    for(i = 0; i < c->length(); i++) {
	      int cval = seqValue(c, fileSeq, i, strand);

	      //	      if(!target->id().compare("61270.2204") && i == 36)
	      //		cerr << "initial " << t << "@" << cval << endl;

	      if(cval)
		for(int j = 0; j < fidxs.size(); j++)
		  (*thisSeq[fidxs[j]])(strand, i + 1) = cval;

	      /*	      if(!target->id().compare("61270.2204") && i == 36) {
		for(int j = 0; j < fidxs.size(); j++)
		  cerr << t << "@" << (int)(*thisSeq[fidxs[j]])(strand, i + 1) << " ";
		cerr << endl;
	      }
	      */
	    }
	  }
	  
	  for(i = 0; i < _numCollapsedVals; i++)
	    target->setSeq(i, thisSeq[i]);
	  
	  continue;
	}
	
	const TClass1 & ctarget = (const TClass1 &)(*target);
	
	for(t = 0; t < NUM_STRANDS; t++) {
	  strand = (TStrand)t;
	  
	  char *fileSeq = c->getSeq(strand);
	  
	  assert(target->length(strand) == c->length(strand));
	  
	  for(i = 0; i < c->length(); i++) {
	    int cval = seqValue(c, fileSeq, i, strand);
	    
	    for(int j = 0; j < fidxs.size(); j++) {
	      int nval = ctarget(strand, fidxs[j], i + 1);
	      //	      if(!target->id().compare("61270.2204") && i == 36)
	      //		cerr << "upd to " << fidxs[j] << " " << t << "@" << cval << "+" << nval << "=";	      
	      nval = collapsedValue(nval, cval);
	      //	      if(!target->id().compare("61270.2204") && i == 36)
	      //		cerr << nval << endl;
	      if(nval)
		(*target)(strand, fidxs[j], i + 1) = nval; 

	    }
	    /*
	    if(!target->id().compare("61270.2204") && i == 36) {
	      for(int j = 0; j < this->numCollapsedVals(); j++)
		cerr << t << "@" << (int)((const TClass1 &)*target)(strand, j, i + 1) << " ";
	      cerr << endl;
	    }
	    */
	  }
	}
      }
      
      return target;
      
    }  
    
    /**
     * A member function that collapses all sequences of all files into *this'.
     * as shifted bits.
     */
    void loadFileSequences() {
      
      assert(files2Collapse.size());

      if(!overrideSigma)
	this->sigma = files2Collapse[0]->alphabet();

      if(!loadAllSeqs)
	return;
      
      vector<std::string> &ids = files2Collapse[0]->seqNames();
      
      if(!ids.size()) {
	assert(0);
	throw new EXCEPTION(BAD_USAGE, files2Collapse[0]->fileName() + " must have mapped ids");
      }
      
      for(int i = 0; i < ids.size(); i++) {
	BasicSeq *c = loadFileSequence(ids[i]);
	if(!c) {
	  assert(0);
	  throw EXCEPTION(CONTIG_UNAVAILABLE, ids[i] + " not found for collapsing");
	}
	sequences().push_back(c);
      }
      
      cerr  << sequences().size() << " BasicSeq* objects collapsed\n";
    }  
    
  
    /**
     * A member function that creates a hashmap, with the sequence identifiers
     * as they appear in the file in the header information, and the sequence 
     * themselves as values. the mapping is created only if the class member 
     * mapIds is true.
     */
    void mapInputSeqs() {
      if(!mapIds)
	return;

      _seqNames = files2Collapse[0]->seqNames();

      if(loadAllSeqs) {
	list<BasicSeq *>::iterator it = _seqs.begin();
	std::map<std::string, BasicSeq*>::iterator it2;
	
	for( ; it != _seqs.end(); it++) {
	  BasicSeq *c = *it;
	  it2 = seqs.find(c->id());
	  
	  if(it2 != seqs.end())
	    throw EXCEPTION( NOT_SUPPORTED, string("duplicated identifiers for file ") + filename);
	  
	  seqs[c->id()] = c;
	}
	
	return;
	
      }
      
    }

    
    /* 
     * This function computes the number of permutations of gramLength in
     * files.size() without repetition and without order. This number
     * will be the size of the vector of integers at each position
     */
    
    void computeNumCollapsedVals() {

      vector< vector<int> > perms;
      _numCollapsedVals = Utils::permutations(_gramLength, files2Collapse.size(), perms);
      _indexes = new vector<int> [files2Collapse.size()];
      
      assert(_numCollapsedVals == perms.size());
      
      for(int i = 0; i < perms.size(); i++) {
	assert(_gramLength == perms[i].size());

	std::string collapsed_fns = "";
	for(int j = 0; j < perms[i].size(); j++)  {
	  _indexes[perms[i][j]].push_back(i);
	  collapsed_fns += files2Collapse[perms[i][j]]->fileName() + " ";
	}
	//	cerr << " names " << _collapsedfNames.size() << " " << collapsed_fns << endl;
	_collapsedfNames.push_back(collapsed_fns);
      }
    }    
    
    void print(BasicSeq *c) {
      
      TClass1 &target = (TClass1 &)*c;

      //      if(target.id().compare("61270.2204"))
      //	return;
	
      cerr << ">" << target.id() << endl;
      
      for(int fileNo = 0; fileNo < files2Collapse.size(); fileNo++) {
	cerr << fileNo << " " << _indexes[fileNo].size();
	for(int k = 0; k < _indexes[fileNo].size(); k++)
	  cerr << " " << _indexes[fileNo][k];
	cerr << endl;
      }
      
      for(int t = 0; t < NUM_STRANDS; t++) {
	for(int i = 0; i < target.length(); i++) {
	  cerr << t << "|" << i;
	  for(int j = 0; j < target.numSequences(); j++) 
	    cerr << " " << (int)((const TClass1 &)target)((TStrand)t, j, i);;
	  cerr << endl;
	}
      }
    }
    


    virtual int collapsedValue(int word, int suffix) = 0;
    virtual int seqValue(BasicSeq *c, char *fileSeq,
                         int pos, TStrand strand) = 0;
    
    virtual ~InpFileCollapser() { 
      if(_indexes)
        delete [] _indexes;
    }
    
  };

  /**
   * The InpFileCollapser class is a subtype of InpFile, which allows the 
   * user to collapse multiple InpFile type objects which contain objects 
   * of type Fasta Sequence that have an alphabet definition.
   * Also, all collapsing files must contain sequence based on the same alphabet.
   *
   ***************************************************************************/
  template<class TClass1, class TClass2>
    class InpFastaFileCollapser : public InpFileCollapser<TClass1, TClass2> {
   public:

    /**
     * Default constructor
     * @param name A unique resource name
     * @param chars4BitOn if any of these character appears in file sequence, 
     * the On value (1) is bit-shifted into the character of this' sequence 
     * appearing at the same position
     */
    InpFastaFileCollapser(std::string name, 
                          int gramLength, bool loadAllSeqs, 
			  Sigma *sigma = NULL)
      : InpFileCollapser<TClass1, TClass2>(name, gramLength, 
					   loadAllSeqs, sigma) {
      
    }
    
    /**
     * Constructor from a Header string definition.
     * @param params The Header string definition, loaded as a vector of 
     * strings. The Header has the following form:\n\n Resource name  
     * BitFileCollapser numInpFiles InpFileobj1 [..]\n\n
     * The description of the fields could be found in the other
     * constructor(s)
     * @param offset The index for vector params
     * @param re A pointer to the ResourceEngine object.
     */    
    InpFastaFileCollapser(std::vector<std::string> &params,
                          int & offset,
                          ResourceEngine *re) :
    InpFileCollapser<TClass1, TClass2>(params, offset, re) {

      this->loadFileSequences();
      this->mapInputSeqs();
      this->computeAlphabetSize(this->alphabet()->alphabetSize());
    }

    inline int collapsedValue(int word, int suffix) {
      return this->sigma->addSuffix(word, suffix);
    }
    
    inline int seqValue(BasicSeq *c, char *fileSeq,
                        int pos, TStrand strand) {
      
      return this->sigma->indC(fileSeq[pos]);
    }


    /**
     * @see Resource::saveHeader(::ofstream &)
     */
    void saveHeader(std::ofstream &fd) {      
      InpFileCollapser<TClass1, TClass2>::saveHeader(fd);
      fd << endl;
    }
   
    virtual ~InpFastaFileCollapser() { }
    
  };


  template<class TClass1, class TClass2>
    class InpPhaseFastaFileCollapser : public InpFileCollapser<TClass1, TClass2> {
   public:
    /**
     * Default constructor
     * @param name A unique resource name
     */
    InpPhaseFastaFileCollapser(std::string name,
                               int gramLength, bool loadAllSeqs,
			       Sigma *sigma = NULL)
      : InpFileCollapser<TClass1, TClass2>(name, gramLength, loadAllSeqs) {

    }    

    InpPhaseFastaFileCollapser(std::vector<std::string> &params,
			       int & offset,
			       ResourceEngine *re) 
      : InpFileCollapser<TClass1, TClass2>(params, offset, re) {

      this->loadFileSequences();
      this->mapInputSeqs();
      this->computeAlphabetSize(3*this->alphabet()->alphabetSize());

    }    

    inline int collapsedValue(int word, int suffix) {
      return this->sigma->alphabetSize()*3*word + suffix;
    }    

    inline int seqValue(BasicSeq *c, char *fileSeq,
			int pos, TStrand strand) {
      
      int phase = ((EdgeAnnotSeq *)c)->activePhase(pos, strand);
      int value = this->sigma->indC(fileSeq[pos]);
      return value*3 + phase;

    }

    /**
     * @see Resource::saveHeader(::ofstream &)
     */
    void saveHeader(std::ofstream &fd) {      
      InpFileCollapser<TClass1, TClass2>::saveHeader(fd);
      fd << endl;
    }

    ~InpPhaseFastaFileCollapser() { }
  };

  /**
   * The InpBitFileCollapser class is a subtype of InpFile, which allows the 
   * user to bit-collapse multiple InpFile type objects which contain objects 
   * of type Sequence.
   * Also, all collapsing files must contain sequence based on the same alphabet.
   *
   ***************************************************************************/
  template<class TClass1, class TClass2>
    class InpBitFileCollapser : public InpFileCollapser<TClass1, TClass2> {
   protected:
    std::string chars4BitOn;

   public:
    /**
     * Default constructor
     * @param name A unique resource name
     * @param chars4BitOn if any of these character appears in file sequence, 
     * the On value (1) is bit-shifted into the character of this' sequence 
     * appearing at the same position
     */
    InpBitFileCollapser(std::string name,
                        int gramLength,
			bool loadAllSeqs,
                        std::string chars4BitOn,
			Sigma *sigma = NULL)
      : InpFileCollapser<TClass1, TClass2>(name, gramLength,
					   loadAllSeqs, sigma) {
      
      this->chars4BitOn = chars4BitOn;

    }    

    InpBitFileCollapser(std::vector<std::string> &params,
                        int & offset,
                        ResourceEngine *re)
      : InpFileCollapser<TClass1, TClass2>(params, offset, re) {

      chars4BitOn = params[offset++];
      this->loadFileSequences();
      this->mapInputSeqs();
      this->computeAlphabetSize(2);
    
    }

    inline int collapsedValue(int word, int suffix) {
      return 2*word + suffix;
    }

    inline int seqValue(BasicSeq *c, char *fileSeq,
				int pos, TStrand strand) {
      
      return (chars4BitOn.find(fileSeq[pos]) != string::npos) ? 1 : 0;
    }

    void saveHeader(std::ofstream &fd) {      
      InpFileCollapser<TClass1, TClass2>::saveHeader(fd);
      fd << " " << chars4BitOn << endl; 
    }

    virtual ~InpBitFileCollapser() { }
  };

  template<class TClass1, class TClass2>  
    class InpPhaseBitFileCollapser : public InpFileCollapser<TClass1, TClass2> {
   protected:
    std::string chars4BitOn;

   public:
    /**
     * Default constructor
     * @param name A unique resource name
     * @param chars4BitOn if any of these character appears in file sequence, 
     * the On value (1) is bit-shifted into the character of this' sequence 
     * appearing at the same position
     */
    InpPhaseBitFileCollapser(std::string name,
                             int gramLength,
			     bool loadAllSeqs,
                             std::string chars4BitOn,
			     Sigma *sigma = NULL)
      : InpFileCollapser<TClass1, TClass2>(name, gramLength,
					   loadAllSeqs, sigma) {

      this->chars4BitOn = chars4BitOn;
      
    }
    

    InpPhaseBitFileCollapser(std::vector<std::string> &params,
                             int & offset,
                             ResourceEngine *re) 
      : InpFileCollapser<TClass1, TClass2>(params, offset, re) {
      
      chars4BitOn = params[offset++];
      this->loadFileSequences();
      this->mapInputSeqs();
      this->computeAlphabetSize(6);
      
    }
 
    inline int collapsedValue(int word, int suffix) {
      return 6*word + suffix;
    }

    inline int seqValue(BasicSeq *c, char *fileSeq,
			int pos, TStrand strand) {

      int phase = ((EdgeAnnotSeq *)c)->activePhase(pos, strand);
      int value = (chars4BitOn.find(fileSeq[pos]) != string::npos);
      return 3*value + phase;
    }

    void saveHeader(std::ofstream &fd) {      
      InpFileCollapser<TClass1, TClass2>::saveHeader(fd);
      fd << " " << chars4BitOn << endl; 
    }    

    ~InpPhaseBitFileCollapser() { }
  };

}

#endif
