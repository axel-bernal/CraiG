#include "FSM.h"
#include "TypeDefs.h"

/****************************************************************************
* FSM.cpp - part of the lless namespace, a general purpose
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


#include <string>

namespace lless {

  /**
   * Minimal constructor. 
   * @param strand the strand in which the nodes and edges of this FSM object
   * are defined.
   */
  FSM::FSM(TStrand strand) : Resource() {
    unsigned int i;
    this->numContexts = 0;
    this->parseStrand = strand;

    _numParseWords = 0;
    _numParseNodes = 0;
    _numParseEdges = 0;

    for(i = 0; i < MAX_NUMPARSE_EDGES + 1; i++) 
      _edges[i] = NULL;

    for(i = 0; i < MAX_NUMPARSE_NODES + 2; i++)
      _nodes[i] = NULL;

    for(i = 0; i < MAX_NUMPARSE_WORDS + 1; i++)
      _words[i] = NULL;

  }

  /**
   * Constructor from a FSM topology stored in a file
   * @param topology the FSM nodes, edges and node length definitions
   * @param numContexts the number of context levels present in the sequence.
   * Usually variable length nodes have context-dependent lengths.
   * @param strand the strand in which the nodes and edges of this FSM object
   * are defined.
   */
  FSM::FSM(::ifstream &topology, 
           int numContexts, 
           TStrand strand) : Resource() {

    this->numContexts = numContexts;
    this->parseStrand = strand;
    retrieveContents(topology);
  }

  /**
   * Constructor from a Header string definition.
   * @param params The Header string definition, loaded as a vector of strings.
   * The Header has the following form:\n\n
   * Resource name FSM topologyFile numContexts strand \n\n
   * The description of the fields could be found in 
   * the other constructor(s)
   * @param offset The index for vector params
   * @param re A pointer to the ResourceEngine object.
   *
   **/
  FSM::FSM(std::vector<std::string> &params, 
           int &offset, 
           ResourceEngine *re) : Resource(params, offset, false) {

    std::string file = params[offset++];

    if(!sscanf(params[offset++].c_str(), "%d", &numContexts))
      assert(0);

    parseStrand = Utils::stringToTStrand(params[offset++]);

    if(file.compare("null") == 0)
      // model to be read inline @ ResourceEngine scope
      _contentsInline = true; 
    else {
      std::ifstream fd(file.c_str());
      assert(fd);
      retrieveContents(fd);
    }  
  }


  /**
   * @see Resource::retrieveContents(std::ifstream &)
   */
  void FSM::retrieveContents(std::ifstream &topology) {
    
    assert(topology);
    std::string line;
    Node *thisNode;

    unsigned int i;

    /*
     * Initialization
     */
    _numParseWords = 0;
    _numParseNodes = 0;
    _numParseEdges = 0;

    for(i = 0; i < MAX_NUMPARSE_EDGES; i++) 
      _edges[i] = NULL;

    for(i = 0; i < MAX_NUMPARSE_NODES; i++) 
      _nodes[i] = NULL;

    for(i = 0; i < MAX_NUMPARSE_WORDS + 1; i++)
      _words[i] = NULL;

    /*
     * Creating the begin and end-of-sequence  synchronization nodes
     */
    _nodes[SYNC_END_STATE] = new VarLenNode(SYNC_END_STATE, 
                                            std::string("$"), 
                                            SYNC_STATE, 
                                            NO_NODE_INST, 
                                            1, 3, 3,
                                            STRAND_FWD);

    namedNodes["$"] = _nodes[SYNC_END_STATE];

    /*
     * Creating different expressions which are to be recognized in the input
     */
    boost::RegEx rExComment("^\\s*\\#.*");
    boost::RegEx rExEndofFile("^\\s*\\/\\/\\s*$");
    boost::RegEx rExEdge("Edge\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s*(\\S*)\\s*$");
    boost::RegEx rExWord("Word\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s*(\\S*)\\s*$");
    boost::RegEx rExVarLenNode("VLNode\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s*(\\S*)\\s*$");
    boost::RegEx rExVarPhaseNode("VPNode\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s+(\\S+)\\s*(\\S*)\\s*$");
    boost::RegEx rExFixedLenNode("FLNode\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\-?\\d+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s*(\\S*)\\s*$");
    boost::RegEx rExTransition("Transition\\s+from\\s+(\\S+)\\s+to\\s+\\[\\s*(.+)\\s*\\]\\s*");
    boost::RegEx rExLengths("Lengths\\s+(\\S+)\\s+\\[\\s*(.+)\\s*\\]\\s*");

    while(std::getline(topology, line), !topology.eof()) {

      fsmContents += line + std::string("\n");

      if(rExComment.Match(line)) // line is a comment
        continue; 

      if(rExEndofFile.Match(line))
        break;

      if(rExEdge.Search(line)) // an edge
        createEdge(rExEdge);
      else if(rExWord.Search(line)) // a word
        createWord(rExWord);
      else if(rExVarLenNode.Match(line)) // a variable-length node
        createNode(rExVarLenNode, VARLEN_NODE, numContexts);
      else if(rExVarPhaseNode.Match(line)) // a variable-length node
        createNode(rExVarPhaseNode, VARPHASE_NODE, numContexts);
      else if(rExFixedLenNode.Match(line)) // a fixed-length node
        createNode(rExFixedLenNode, FIXEDLEN_NODE, numContexts);
      else if(rExTransition.Match(line)) // a transition
        createTransition(rExTransition);
      else if(rExLengths.Match(line)) {
        std::string nodeName(rExLengths[1]);
        thisNode = node(nodeName);
        boost::RegEx rExDelim("\\s+");
        ::vector<std::string> lengthSplits;
        ::string lengthString(rExLengths[2]);
        rExDelim.Split(lengthSplits, lengthString, 0, numContexts);
        int context;
        int length1, length2;
        boost::RegEx rExLenEntry("\\s*(\\d+)\\s*\\:\\s*(\\d+)\\s*\\->\\s*(\\-?\\d+)\\s*");

        for(i = 0; i < lengthSplits.size(); i++) {
          bool matched;
          matched = rExLenEntry.Match(lengthSplits[i]);
          assert(matched);
          sscanf(rExLenEntry[1].c_str(),"%d", &context);
          sscanf(rExLenEntry[2].c_str(),"%d", &length1);
          sscanf(rExLenEntry[3].c_str(),"%d", &length2);
          thisNode->setMinLength(length1, context);

          if(length2 < 0) 
            length2 = INT_MAX;

          thisNode->setMaxLength(length2, context);
        }
      }
      else if(line.length() != 0)
        cerr << "Warning! unknown parsing line : " << line << endl;
    }

    assert(_nodes[SYNC_BEG_STATE]->name().compare("^") == 0);

    /*
     * Construct the model with reverse flow
     */
    for(int i = 0; i < numParseNodes(); i++) {
      Node *curr = node((TParseNode)i);

      for(unsigned int j = 0; j < curr->nextNodes().size(); j++) {
        Node *next = curr->nextNodes()[j];
        Edge *e = edge(curr->id(), next->id());
        next->addPrevNode(-e->nextNodePos(next->id()), e, curr);
      }

    }
  }

  /** 
   * A member function that creates an Edge object out of a string that
   * matches the regular expression for Edge object definitions.
   * @param r the input regular expression
   */
  void FSM::createEdge(boost::RegEx & r) {
    assert(_numParseEdges < MAX_NUMPARSE_EDGES);
    TParseEdge idt = (TParseEdge)_numParseEdges;
    TStrand strand = Utils::stringToTStrand(r[5]);
    _edges[idt] = new Edge(idt, 
                           r[1],
                           TypeDefs::stringToTEdgeId2(r[2]),
                           Utils::stringToBoolean(r[3]), 
                           Utils::stringToBoolean(r[4]), 
                           strand);
    if(r[6].length() > 0) {
      std::string compEdge(r[6]);
      _edges[idt]->setComplementEdge(edge(compEdge));
      edge(compEdge)->setComplementEdge(_edges[idt]);
    }

    namedEdges[r[1]] = _edges[idt];
    _numParseEdges++;
  }

  /** 
   * A member function that creates an Word object out of a string that
   * matches the regular expression for Word object definitions.
   * @param r the input regular expression
   */
  void FSM::createWord(boost::RegEx & r) {
    assert(_numParseWords < MAX_NUMPARSE_WORDS);
    TParseWord ids;
    int numPhases, maxInPhases, maxOutPhases;

    ids = (TParseWord)_numParseWords;
    
    if(!sscanf(r[4].c_str(),"%d", &numPhases))
      assert(0);
    if(!sscanf(r[5].c_str(),"%d", &maxInPhases))
      assert(0);
    if(!sscanf(r[6].c_str(),"%d", &maxOutPhases))
      assert(0);
    
    TStrand strand = Utils::stringToTStrand(r[7]);
    _words[ids] = new Word(ids, 
                           r[1], 
                           TypeDefs::stringToTNodeId2(r[2]), 
                           TypeDefs::stringToTNodeId2(r[3]),                            
                           numPhases, 
                           maxInPhases, 
                           maxOutPhases,
                           strand);

    std::string compWord = std::string(r[8]);

    if(compWord.length() > 0) {
      _words[ids]->setComplementWord(word(compWord));
      word(compWord)->setComplementWord(_words[ids]);
    }
    
    _numParseWords++;
    namedWords[r[1]] = _words[ids];
  }

  /** 
   * A member function that creates an Node object out of a string that
   * matches the regular expression for Node object definitions.
   * @param r the input regular expression
   * @param t the enumerate type of the Node object to be created.
   * @param numContexts the number of context levels present in the sequence.
   */
  void FSM::createNode(boost::RegEx & r, TNodeType t, int numContexts) {
    assert(_numParseNodes < MAX_NUMPARSE_NODES);
    TParseNode ids;
    int numPhases, maxInPhases, maxOutPhases, period;

    ids = (TParseNode)_numParseNodes;

    if(!sscanf(r[4].c_str(),"%d", &numPhases))
      assert(0);
    if(!sscanf(r[5].c_str(),"%d", &maxInPhases))
      assert(0);
    if(!sscanf(r[6].c_str(),"%d", &maxOutPhases))
      assert(0);

    TStrand strand = Utils::stringToTStrand(r[7]);
    std::string compNode;

    switch(t) {
    case FIXEDLEN_NODE:      
      _nodes[ids] = new FixedLenNode(ids, 
                                     r[1],
                                     TypeDefs::stringToTNodeId2(r[2]), 
                                     TypeDefs::stringToTNodeId3(r[3]), 
                                     numPhases, 
                                     maxInPhases, 
                                     maxInPhases,
                                     strand);
      compNode = std::string(r[8]);
      break;

    case VARLEN_NODE: 
      _nodes[ids] = new VarLenNode(ids, 
                                   r[1], 
                                   TypeDefs::stringToTNodeId2(r[2]), 
                                   TypeDefs::stringToTNodeId3(r[3]), 
                                   numPhases, 
                                   maxInPhases, 
                                   maxOutPhases,
                                   strand);
      compNode = std::string(r[8]);
      break;

    case VARPHASE_NODE: {
      TParseEdge disruptor = INVALID_EDGE;
      if(r[8].compare("INVALID_EDGE"))
	disruptor = edge(r[8])->id();
      
      if(!sscanf(r[6].c_str(),"%d", &maxOutPhases))
        assert(0);

      _nodes[ids] = new VarPhaseNode(ids, 
                                     r[1],
                                     TypeDefs::stringToTNodeId2(r[2]), 
                                     TypeDefs::stringToTNodeId3(r[3]), 
                                     numPhases, 
                                     maxInPhases, 
                                     maxOutPhases,
                                     strand,
                                     disruptor);
      compNode = std::string(r[9]);
      break;
    }      
    default: break;
    }

    if(compNode.length() > 0) {
      _nodes[ids]->setComplementNode(node(compNode));
      node(compNode)->setComplementNode(_nodes[ids]);
    }

    _numParseNodes++;
    namedNodes[r[1]] = _nodes[ids];
  }

  /** 
   * A member function that establishes as transition between two nodes in the
   * FSM object. The transition definition is read out of a string that
   * matches the regular expression for transition definitions.
   * @param r the input regular expression
   */
  void FSM::createTransition(boost::RegEx & r) {
    unsigned int i;
    Edge *nextEdge;
    Node *nextNode;
    int nextEdgePos;

    boost::RegEx rExDelim("\\s+");
    std::string nodeName(r[1]);
    Node *thisNode = node(nodeName); 
    std::string nextNodeString(r[2]);

    std::vector<std::string> nextNodeSplits;
    rExDelim.Split(nextNodeSplits, nextNodeString, 0, MAX_NUMPARSE_NODES);

    boost::RegEx rExNextNode("(\\-?\\d+)\\:(\\S+)\\:(\\S+)");

    for(i = 0; i < nextNodeSplits.size(); i++) {
      bool matched;
      matched = rExNextNode.Match(nextNodeSplits[i]);
      std::string nextNodeName(rExNextNode[3]), edgeName(rExNextNode[2]);
      assert(matched);
      sscanf(rExNextNode[1].c_str(),"%d", &nextEdgePos);
      nextEdge = edge(edgeName);
      nextNode = node(nextNodeName);
      nextEdge->setNextNodePos(nextNode->id(), -nextEdgePos);

      //        cerr << ids << " " << nextEdgePos << " " << nextEdge << " " << nextNode << endl;
      if(thisNode->name().compare("^") == 0)
        nextNode->setSyncBegFlag(true);  

      if(nextNode->name().compare("$") == 0) {
        thisNode->addNextEdge(nextEdgePos, nextEdge, nextNode);
        thisNode->setSyncEndFlag(true);
      }
      else
        thisNode->addNextNode(nextEdgePos, nextEdge, nextNode);
    }

    thisNode->setNumContexts(numContexts);
  }

  void FSM::restrictSyncBegStates(TParseNode nsyncNode) {
    vector<Node *> nextNodes = node(nsyncNode)->nextNodes();
    for(int j = 0; j < nextNodes.size(); )
      if(!nextNodes[j]->isSyncBeg())
	nextNodes.erase(nextNodes.begin() + j, nextNodes.begin() + j + 1);
      else
	j++;
    
    node(SYNC_BEG_STATE)->setNextNodes(nextNodes);

  }

  /**
   * A member function that removes the Edge object between prevNode and
   * nextNode.
   */
  void FSM::removeEdge(TParseNode prevNode, TParseNode nextNode) {
    vector<Node *> & nextNodes = node(prevNode)->nextNodes();

    for(unsigned int j = 0; j < nextNodes.size();) {

      if(nextNode == nextNodes[j]->id())
        nextNodes.erase(nextNodes.begin() + j, nextNodes.begin() + j + 1);
      else
        j++;  
    }
  }

  /**
   * A member function that removes any node with TNodeId3 type id3 
   */

  void FSM::removeId2Node(TNodeId2 id2) {
    int i;
    unsigned int j;  

    for(i = 0; i < numParseNodes(); i++) {
      Node *nod = node((TParseNode)i);
      if(!nod) 
	continue;

      /*
       * Remove nodes even if they belong to another strand
       */
      vector<Node *> & nextNodes = nod->nextNodes();
      for(j = 0; j < nextNodes.size();)
        if(nextNodes[j]->id2() == id2)
          nextNodes.erase(nextNodes.begin() + j, nextNodes.begin() + j + 1);
        else
          j++;

      vector<Node *> & prevNodes = nod->prevNodes();
      for(j = 0; j < prevNodes.size();)
        if(prevNodes[j]->id2() == id2)
          prevNodes.erase(prevNodes.begin() + j, prevNodes.begin() + j + 1);
        else
          j++;
      
    }
  }
  
  /**
   * @return the Edge object whose name is edgeName, if not found it throws
   * an EXCEPTION
   */
  Edge *FSM::edge(std::string & edgeName) {
    Edge *thisEdge = findEdge(edgeName);

    if(thisEdge == NULL) {
      assert(0);
      throw EXCEPTION(PARSE_ERROR, 
                          std::string("undefined Edge ") + edgeName);
    }
    return thisEdge;
  }

  /**
   * @return the Edge object whose name is edgeName, if not found it returns
   * NULL
   */
  Edge *FSM::findEdge(std::string & edgeName) {
    std::map<std::string, Edge *>::iterator it;
    it = namedEdges.find(edgeName);

    if(it == namedEdges.end())
      return NULL;

    return it->second;    
  }

  /**
   * @return the Node object whose name is nodeName, if not found it throws
   * an EXCEPTION
   */
  Node *FSM::node(std::string & nodeName) {
    Node *thisNode = findNode(nodeName);

    if(!thisNode) {
      assert(0);
      throw EXCEPTION( PARSE_ERROR, 
                             std::string("undefined Node ") + nodeName);
    }
    return thisNode;
  }

  /**
   * @return the Node object whose name is nodeName, if not found it returns
   * NULL
   */
  Node *FSM::findNode(std::string & nodeName) {
    std::map<std::string, Node *>::iterator it;
    it = namedNodes.find(nodeName);

    if(it == namedNodes.end())
      return NULL;

    return it->second;    
  }

  /**
   * @return the Word object whose name is wordName, if not found it throws
   * an EXCEPTION
   */
  Word *FSM::word(std::string & wordName) {
    Word *thisWord = findWord(wordName);

    if(!thisWord) {
      assert(0);
      throw EXCEPTION( PARSE_ERROR, 
                           std::string("undefined Word ") + wordName);
    }
    return thisWord;
  }

  /**
   * @return the Word object whose name is wordName, if not found it returns
   * NULL
   */
  Word *FSM::findWord(std::string & wordName) {
    std::map<std::string, Word *>::iterator it;
    it = namedWords.find(wordName);

    if(it == namedWords.end())
      return NULL;

    return it->second;    
  }

  /**
   * A member function that sets the parse strand of the FSM object to be
   * strand. If any entry in the list of nodes and edges which form part of
   * the finite state automata are in a strand that is incompatible(they
   * belong to another strand) with this function's parameter then they 
   * are removed immediately.
   * @param strand the new parse strand.
   */
  void FSM::setParseStrand(TStrand strand) {
    int i;
    unsigned int j;  

    if(strand == parseStrand)
      return;

    if(strand == BOTH_STRANDS)
      throw EXCEPTION(NOT_ANNOTATED, "Cannot restore FSM's state to BOTH_STRANDS");

    parseStrand = strand;

    for(i = 0; i < numParseNodes(); i++) {
      Node *nod = node((TParseNode)i);
      if(!nod)
        continue;

      vector<Node *> & nextNodes = nod->nextNodes();
      for(j = 0; j < nextNodes.size();)
        /*
         * Remove nodes if they belong to another strand
         */
        if(nextNodes[j]->strand() != strand)
          nextNodes.erase(nextNodes.begin() + j, nextNodes.begin() + j + 1);
        else
          j++;

      vector<Node *> & prevNodes = nod->prevNodes();
      for(j = 0; j < prevNodes.size();)

        if(prevNodes[j]->strand() != strand)
          prevNodes.erase(prevNodes.begin() + j, prevNodes.begin() + j + 1);
        else
          j++;
      
    }
    
    node(SYNC_BEG_STATE)->setStrand(strand);
  }
  
  FSM::~FSM() {
    int i;

    for(i = 0; i < numParseNodes(); i++) {
      
      if(_nodes[i])
        delete _nodes[i];
      
      _nodes[i] = NULL;
    }

    delete _nodes[SYNC_END_STATE];

    for(i = 0; i < numParseEdges(); i++) {

      if(_edges[i])
        delete _edges[i];

      _edges[i] = NULL;
    }

  }
  
}
