#include "SequenceUtils.h"
#include <sys/stat.h>
#include <iostream>
#include "Sigma.h"
#include <string>
#include <ctype.h>

/****************************************************************************
* SequenceUtils.cpp - part of the lless namespace, a general purpose
*                     linear semi-markov structure prediction library
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

  static boost::RegEx rExFastaId("^>(\\S+)\\s*");
  static boost::RegEx rExSeqTagId("^>(\\S+)\\s+(\\d+)$");
  static boost::RegEx rExMultizId("^>(\\S+)\\s+(\\d+)\\s+(\\S)\\s*$");
  static boost::RegEx rExEdgeAnnotId("^>(\\S+)\\s+(\\d+)\\s+(\\d+)\\s*$");
  static boost::RegEx rExCompSeqMark("^\\s*\\/\\/\\s*$");
  static boost::RegEx rExLabel("^Label\\s+Set\\s+(\\d.*)$");
  static boost::RegEx rExNodeInst("^NODE_INST\\s+(\\S+)\\s+(\\d+)\\s+(\\d.*)$");
  static boost::RegEx rExEdgeInst("^EDGE_INST\\s+(\\S+)\\s+(\\d+)\\s*$");

  static std::map<std::string, TSetType> setTypes;

  SequenceUtils::SequenceUtils _tmpSeqUtilsObj;
  
  SequenceUtils::SequenceUtils() {
    setTypes["DEFAULT_SET"] = DEFAULT_SET;
    setTypes["PRED_SET"] = PRED_SET;
    setTypes["TRAIN_SET"] = TRAIN_SET;
    setTypes["NO_SET"] = NO_SET;
  }

  TSetType SequenceUtils::stringToTSetType(const std::string &s) {
    std::map<std::string, TSetType>::iterator it;
    it = setTypes.find(s);
    if(it == setTypes.end())
      throw EXCEPTION( PARSE_ERROR, std::string("undefined TSetType ")+s);
    return it->second;
  }
  
  Sequence * SequenceUtils::loadFastaSequence(std::ifstream & fd, 
                                              Sigma *sigma,
                                              std::string &line) {
    
    std::string annotSeq[NUM_STRANDS];
    TStrand strand = STRAND_FWD;
    Sequence *c = NULL;

    while(1) {    
      
      if(rExFastaId.Search(line) || fd.eof())  {
        
        if(c) {
          for(int t = 0; t <= (int)strand; t++) {
            c->setSeq((char *)annotSeq[t].c_str(), (TStrand)t);
            annotSeq[t].clear();                       
          }

          /*
           * Sequence in forward strand should always be present
           * if backward strand isn't present, compute it as a
           * natural complement of the forward strand
           */
          if(strand == STRAND_FWD) 
            c->setNatComplementSeq();

          return c;
        }

        if(fd.eof())
          break;

        std::string cId = rExFastaId[1];
        c = new Sequence(cId, NULL, sigma);
        strand = STRAND_FWD;
      }
      else if(rExCompSeqMark.Match(line)) 
        strand = STRAND_COMP;
      else { // it's a sequence
        Utils::removeBlanks(line);
        
        if(!sigma->checkSeq(line)) {
          assert(0);
          throw EXCEPTION(STRANGE_CHAR, line);
        }
        
        annotSeq[strand] += line;

      }
      
      std::getline(fd, line);
    }
    
    return NULL;

  }         

  Sequence* SequenceUtils::loadTagSequence(std::ifstream & fd, 
					   FSM &fsm,
					   std::string &line) {
    
    int seqLen;
    int pos;
    Sequence *c = NULL;
    std::string seqId;
    int rank = -1;
    double score = 0;
    int initPhase = 0, endPhase = 0;
    SeqTags *currTags = NULL;
    
    while(1) {    
      int phaseBreak = -1;
      Tag *tag = NULL;
    
      if(rExSeqTagId.Match(line) || fd.eof())  {        
        if(c) {
	  if(!currTags)
	    throw EXCEPTION(CONTIG_UNAVAILABLE, seqId);
          return c;
	}	

        if(fd.eof())
          break;
	
        std::string cId = rExSeqTagId[1];
        if(!sscanf((char *)rExSeqTagId[2].c_str(), "%d", &seqLen)) 
          assert(0);
	
	c = new Sequence(cId);
	c->setLength(seqLen);
      }
      else if(rExLabel.Match(line)) {
	int numScans = sscanf(rExLabel[1].c_str(), "%d %lf %d %d", 
			      &rank, &score, &initPhase, &endPhase);
	
	assert(numScans);
	currTags = &c->addTags(rank);
	if(numScans == 4)
	  currTags->initialize(score, initPhase, endPhase, rank);
      }
      else {
	if(rExNodeInst.Match(line)) {
	  std::string nodeName(rExNodeInst[1]);
	  Node *node = fsm.node(nodeName);
	  int len;
	  
	  if(!sscanf(rExNodeInst[2].c_str(), "%d", &pos))
	    assert(0);  
	  int numScans = sscanf(rExNodeInst[3].c_str(), "%d %d", 
				&len, &phaseBreak);
	  
	  assert(numScans);         
	  if(numScans < 2)  phaseBreak = -1;
	  
	  tag = new NodeInst(node->id3(), node->id(), pos, node->strand(), len);
	}
	else if(rExEdgeInst.Match(line)) {
	  std::string edgeName(rExEdgeInst[1]);
	  Edge *edge = fsm.edge(edgeName);
	  
	  if(!sscanf(rExEdgeInst[2].c_str(), "%d", &pos))
	    assert(0);
	  
	  tag = new EdgeInst(edge->id2(), edge->id(), pos, edge->strand());
	}
	else if(line.length() > 0) {
	  assert(0);
	  throw EXCEPTION(PARSE_ERROR, line);
	}
	else {
	  std::getline(fd, line);
	  continue;
	}

	currTags->tpush_back(tag);

	if(phaseBreak > 0)
	  currTags->insertPhaseBreak(currTags->back(), phaseBreak);       
	
      }
      std::getline(fd, line);
    }

    return NULL;
  }

  // for handling .state evidence files
  Sequence * SequenceUtils::loadExtFastaSequence(std::ifstream & fd,
                                                 Sigma *sigma,
                                                 std::string &line) {

    std::string annotSeq[NUM_STRANDS];
    TStrand strand = STRAND_FWD;
    char symbol;
    int len;
    Sequence *c = NULL;

    while(1) {    
      
      if(rExFastaId.Search(line) || fd.eof())  {
        
        if(c) {
          for(int t = 0; t <= (int)strand; t++) {
            c->setSeq((char *)annotSeq[t].c_str(), (TStrand)t);
	    //	    cerr << c->id() << " " << t << endl << annotSeq[t] << endl;
            annotSeq[t].clear();           
          }

          /*
           * Sequence in forward strand should always be present
           * if backward strand isn't present, compute it as a
           * natural complement of the forward strand
           */
          if(strand == STRAND_FWD) 
            c->setNatComplementSeq();

          return c;
        }

        if(fd.eof())
          break;

        std::string cId = rExFastaId[1];
        c = new Sequence(cId, NULL, sigma);
        strand = STRAND_FWD;
      }
      else if(rExCompSeqMark.Match(line)) 
        strand = STRAND_COMP;
      else { 
        // it must match a line="symbol x length abs_pos"
        // or it can be a symbol line
	if(sscanf(line.c_str(), "%c x %d", &symbol, &len) == 2) {
          if(!sigma->checkSymbol(symbol)) {
            assert(0);
            throw EXCEPTION(STRANGE_CHAR, line);
          }

          for(int i = 0; i < len; i++) 
            annotSeq[strand].push_back(symbol);

        }
        else {
	  //          Utils::removeBlanks(line);
	           
          if(!sigma->checkSeq(line)) {
            assert(0);
            throw EXCEPTION(STRANGE_CHAR, line);
          }
          
          annotSeq[strand] += line;          
        }
      }
      
      std::getline(fd, line);
    }
    
    return NULL;
    
  }

  // for handling .cons Conservation files
  Sequence * SequenceUtils::loadMultizSequence(std::ifstream &fd, 
                                                  Sigma *sigma,
                                                  std::string &line) {
    
    char *annotSeq[NUM_STRANDS];
    TStrand strand = STRAND_FWD;
    Sequence *c = NULL;
    int seqLen = 0;
    char fillChar = '.';

    for(int t = 0; t < NUM_STRANDS; t++) 
      annotSeq[t] = NULL;

    while(1) {    
      
      if(rExMultizId.Match(line) || fd.eof())  {
        
        if(c) {
          assert(seqLen);
          
          for(int t = 0; t <= (int)strand; t++) {
            if(!annotSeq[t]) {
              annotSeq[t] = new char[seqLen + 1];
              annotSeq[t][seqLen] = '\0';
              for(int i = 0; i < seqLen; i++)
                annotSeq[t][i] = fillChar;
            }
            
            c->setSeq(annotSeq[t], (TStrand)t);

            delete [] annotSeq[t];
            annotSeq[t] = NULL;
            
          }

          if(strand == STRAND_FWD) 
            c->setNatComplementSeq();

          return c;
        }
        
        if(fd.eof())
          break;
        
        if(!sscanf((char *)rExMultizId[2].c_str(), "%d", &seqLen)) 
          assert(0);
        
        if(!sscanf((char *)rExMultizId[3].c_str(), "%c", &fillChar)) 
          assert(0);
        
        std::string cId = rExMultizId[1];
        c = new Sequence(cId, NULL, sigma);
        strand = STRAND_FWD;
      }
      else if(rExCompSeqMark.Match(line)) 
        strand = STRAND_COMP;
      else { // it must match a line="pos score"

        int sp = 0;
        while(isblank(line[sp]) && sp < line.length()) sp++;

        if(sp < line.length()) {
          int pos;
          
          if(!annotSeq[strand]) {
            annotSeq[strand] = new char[seqLen + 1];
            annotSeq[strand][seqLen] = '\0';          
            for(int i = 0; i < seqLen; i++)
              annotSeq[strand][i] = fillChar;
          }
          
          if(sscanf(line.c_str(), "%d", &pos) != 1) {
            assert(0);
            throw EXCEPTION(STRANGE_CHAR, line);
          }
          
          sp = 0;
          
          while(isblank(line[sp]) || isdigit(line[sp])) sp++;
          assert(pos + line.length() - sp <= seqLen);
          
          while(sp < line.length()) 
            annotSeq[strand][pos++] = line[sp++];
          
        }
      }
      std::getline(fd, line);
    }
    
    return NULL;
  }
  
  // for handling .signal files
  EdgeAnnotSeq * SequenceUtils::loadEdgeAnnotSequence(std::ifstream & fd,
                                                      Sigma *sigma,
                                                      std::string &line) {

    
    std::string tmpSeq[NUM_STRANDS];
    vector<int> *annotSeq[NUM_STRANDS];
    TStrand strand = STRAND_FWD;
    EdgeAnnotSeq *c = NULL;
    boost::RegEx rExDelim("\\s+");
    int seqLen = 0;
    int numEdgePhases = 0;

    for(int t = 0; t < NUM_STRANDS; t++) 
      annotSeq[t] = NULL;
    
    while(1) {    
      //      cerr << line << " " << tmpSeq[strand].length() << endl;
      if(rExEdgeAnnotId.Match(line) || fd.eof())  {
        
        if(c) {
          for(int t = 0; t <= (int)strand; t++) {
            TStrand mstrand = (TStrand)t;
            c->setSeq(annotSeq[t], tmpSeq[t].length(), mstrand);
            assert(c->length(mstrand) == seqLen);
            annotSeq[t] = NULL;
	    tmpSeq[t].clear();
          }

          /*
           * Sequence in forward strand should always be present
           * if backward strand isn't present, compute it as a
           * natural complement of the forward strand
           */
          if(strand == STRAND_FWD) 
            c->setNatComplementSeq();

          return c;
        }

        if(fd.eof())
          break;

        if(!sscanf((char *)rExEdgeAnnotId[2].c_str(), "%d", &seqLen))
          assert(0);

        if(!sscanf((char *)rExEdgeAnnotId[3].c_str(), "%d", &numEdgePhases))
          assert(0);

        std::string cId = rExEdgeAnnotId[1];
        c = new EdgeAnnotSeq(cId, numEdgePhases, NULL, sigma);
        strand = STRAND_FWD;
      }
      else if(rExCompSeqMark.Match(line)) 
        strand = STRAND_COMP;
      else { // it must match a line="symbol x length back_pos phase"
	char symbol;
	int len, pos = tmpSeq[strand].length() + 1, weight = 0, offset = 3;
	
        if(sscanf(line.c_str(), "%c x %d", &symbol, &len) == 2) {
	  if(!sigma->checkSymbol(symbol)) {
            assert(0);
            throw EXCEPTION(STRANGE_CHAR, line);
          }
        }
        else {
	  symbol = line[0];
          if(!sigma->checkSymbol(symbol)) {
            assert(0);
            throw EXCEPTION(STRANGE_CHAR, line);
          }          
	  offset = len = 1;
        }
	
	vector<std::string> detail;
	rExDelim.Split(detail, line, 0, 99);
	int i;

	//        cerr << c->id() << " " << len << endl;
	
	if(detail.size() < 3) 
	  throw EXCEPTION(BAD_USAGE, "format in edge file is wrong\n");

	if(!annotSeq[strand])
	  annotSeq[strand] = new vector<int> [seqLen + 1];
	
        for(int j = 0; j < len; j++) 
	  if(!annotSeq[strand][pos + j].size())
	    annotSeq[strand][pos + j].push_back((int)symbol);
	
	for(i = offset; i < detail.size(); i += 3) {
	  int phase, backPos;
	  
	  if(!sscanf(detail[i].c_str(), "%d", &backPos))
	    assert(0);
	  if(!sscanf(detail[i + 1].c_str(), "%u", &phase))
	    assert(0);
	  if(!sscanf(detail[i + 2].c_str(), "%d", &weight))
	    assert(0);
        	  
	  for(int j = 0; j < len; j++)  {
	    //	    cerr << "edges[" << tmpSeq[strand].length() + 1 + j << "] = " << pos + j << " " << phase << " " << backPos << " " << weight << " " << strand << endl;
	    annotSeq[strand][pos + j].push_back(phase);
	    annotSeq[strand][pos + j].push_back(backPos);
	    annotSeq[strand][pos + j].push_back(weight);
	  }
	}

        for(int i = 0; i < len; i++) 
          tmpSeq[strand].push_back(symbol);

      }
      
      std::getline(fd, line);

    }
    
    return NULL;
    
  }


  /**
   * A static member function to prepare a file for random access. offsets 
   * for each sequence id will be loaded in seqIndexes
   * @param fd the input file to be readied.
   * @param seqIndexes the hash table containing sequence ids and offsets
   */

  void SequenceUtils::loadSeqIndexes(std::ifstream &fd, 
                                     map<std::string, streampos> & seqIndexes) {

    char last_c, c;
    long offset = 0;
    std::string id;
    while(1) {
      fd.ignore(INT_MAX, '>');
      offset += fd.gcount();
      fd >> id;

      if(fd.eof())
        return;

      seqIndexes[id] = offset - 1;
      offset += id.length();
    }

  }

}
