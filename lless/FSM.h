/****************************************************************************
 * FSM.h - part of the lless namespace, a general purpose
 *         linear semi-markov structure prediction library
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

#ifndef _FSM_H_
#define _FSM_H_

#include <ctype.h>
#include <string.h>
#include <fstream>
#include <stdio.h>
#include <map>
#include <vector>
#include "Utils.h"
#include <TypeDefs.h>
#include <boost/regex.hpp>
#include "ResourceEngine.h"


namespace lless {

    /**
     * An enumerate for describing the different node types.
     */
    typedef enum {
        FIXEDLEN_NODE,
        VARLEN_NODE,
        VARPHASE_NODE,
        VARPFIXEDLEN_NODE
    } TNodeType;

    /**
     * An enumerate for describing the two ways in which transitions flow within
     * the finite state automata
     */
    typedef enum {
        FSM_NEXT,
        FSM_PREV
    } TFlow;

    /**
     * The Edge class represents an edge or a transition in a FSM object.
     */
    class Edge {
    protected:
        TParseEdge _id;
        std::string _name;
        TEdgeId2 _id2;
        /**
         * Relative position of the beginning of the next state
         * or end of the previous from the current signal's position
         */
        int _nextNodePos[MAX_NUMPARSE_NODES+1];
        int _nextNodeP;
        bool _hasSignal;
        bool _toSelf;              //!< Edge goes to same node it departed from
        TStrand _strand;
        Edge *_complementEdge;      //!< This Edge in the complementary strand

    public:
        Edge(TParseEdge id,
             std::string name,
             TEdgeId2 id2,
             bool hasSignal,
             bool toSelf,
             TStrand strand) {

            int i;
            _id = id;
            _name = name;
            _id2 = id2;
            _toSelf = toSelf;
            _hasSignal = hasSignal;
            _strand = strand;

            for(i = 0; i <= MAX_NUMPARSE_NODES; i++)
                _nextNodePos[i] = INT_MIN;

            _nextNodeP = INT_MIN;
            _complementEdge = NULL;
        }

        inline TParseEdge id() {  return _id; }
        inline TEdgeId2 id2() {  return _id2; }
        inline std::string & name() {  return _name; }
        inline bool hasSignal() {  return _hasSignal; }
        inline bool toSelf() {  return _toSelf; }
        inline Edge *complementEdge() {  return _complementEdge; }
        inline TStrand strand() { return _strand; }

        inline int nextNodePos(int node) {
            assert(node >= 0 && node <= MAX_NUMPARSE_NODES && _nextNodePos[node] != INT_MIN);
            return _nextNodePos[node];
        }
        inline int nextNodePos() {
            assert(_nextNodeP != INT_MIN);
            return _nextNodeP;
        }
        inline void setNextNodePos(int node, int d) {
            assert(node >=0 && node <= MAX_NUMPARSE_NODES);
            _nextNodePos[node] = d;
            _nextNodeP = d;
        }

        inline void setComplementEdge(Edge *edge) {
            assert(edge != NULL);
            _complementEdge = edge;
        }

        ~Edge() {;}

    };


    /**
     * The Node base class represents an node or a state in a FSM object.
     */
    class Node {
    protected:
        TParseNode _id;
        std::string _name;
        TNodeId2 _id2;
        TNodeId3 _id3;
        int _numPhases;       //!< Number of phases of Node
        int _maxInPhases;     //!< Max number of phases at the beginning of Node
        int _maxOutPhases;    //!< Max number of phases at the end of Node
        TStrand _strand;
        bool _syncBeg, _syncEnd;  //!< true if sequence can start(end) in Node
        Pair<int> *_length;                   //!< context-length pairs
        TNodeType _nodeType;
        int _defaultAdjacentNode[2];     //!< the default
        int _edgePos[2][MAX_NUMPARSE_NODES+1]; /*!< edge positions between Node
                                                 and any adjacent node */
        Edge* _edge[2][MAX_NUMPARSE_NODES+1];

        vector<Node *> _adjacentNodes[2];
        Node *_complementNode;

    public:
        Node(TParseNode id,
             std::string name,
             TNodeId2 id2,
             TNodeId3 id3,
             int numPhases,
             int maxInPhases,
             int maxOutPhases,
             TStrand strand) {

            this->_id = id;
            this->_name = name;
            this->_id2 = id2;
            this->_id3 = id3;
            this->_numPhases = numPhases;
            this->_maxInPhases = maxInPhases;
            this->_maxOutPhases = maxOutPhases;
            this->_strand = strand;
            this->_syncEnd = false;
            this->_syncBeg = false;
            _length = NULL;

            for(int j = 0; j < 2; j++ ) {
                for(int i = 0; i <= MAX_NUMPARSE_NODES; i++) {
                    _edgePos[j][i] = INT_MIN;
                    _edge[j][i] = NULL;
                }
                _defaultAdjacentNode[j] = INT_MIN;
            }
            _complementNode = NULL;
        }

        inline TParseNode id() {  return _id; }
        inline TNodeId2 id2() {  return _id2; }
        inline TNodeId3 id3() {  return _id3; }
        inline std::string & name() { return _name; }
        inline TNodeType nodeType() {  return _nodeType; }
        inline void setSyncBegFlag(bool syncBeg) {  _syncBeg = syncBeg; }
        inline void setSyncEndFlag(bool syncEnd) {  _syncEnd = syncEnd; }
        inline TStrand strand() {  return _strand; }
        inline void setStrand(TStrand strand) {  _strand = strand; }
        inline int numPhases() { return _numPhases; }
        inline int maxInPhases() { return _maxInPhases; }
        inline int maxOutPhases() { return _maxOutPhases; }
        inline bool isSyncBeg() {  return _syncBeg; }
        inline bool isSyncEnd() {  return _syncEnd; }
        inline Node *complementNode() { return _complementNode; }
        inline ::vector<Node *> & nextNodes() {  return _adjacentNodes[FSM_NEXT]; }
        inline int nextPhase(int phase, int len) {
            return (numPhases() > 1) ? (phase + len) % numPhases() : phase;
        }

        virtual inline TParseEdge phaseDisruptor() {
            return INVALID_EDGE;
        }
        virtual inline int disruptingPhase(int beg, int phase, TStrand strand) {
            return -1;
        }

        inline void setComplementNode(Node *compNode) {
            assert(compNode != NULL);
            _complementNode = compNode;
        }

        inline void setNextNodes(vector<Node *> & nextNodes) {
            _adjacentNodes[FSM_NEXT] = nextNodes;
        }

        inline void addNextNode(int edgePos, Edge *edge, Node *node) {
            assert(node->id() >= 0 && node->id() < MAX_NUMPARSE_NODES);
            _adjacentNodes[FSM_NEXT].push_back(node);
            addNextEdge(edgePos, edge, node);
        }

        inline void addNextEdge(int edgePos, Edge *edge, Node *node) {
            assert(node->id() <= MAX_NUMPARSE_NODES && edge->id() < MAX_NUMPARSE_EDGES && _edge[FSM_NEXT][node->id()] == NULL);
            _edgePos[FSM_NEXT][node->id()] = edgePos;
            _edge[FSM_NEXT][node->id()] = edge;
            if(edge->hasSignal() && edge->strand() == strand())
                _defaultAdjacentNode[FSM_NEXT] = node->id();
        }

        inline int nextEdgePos(int node) {
            assert(node >= 0 && node <= MAX_NUMPARSE_NODES && _edgePos[FSM_NEXT][node] != INT_MIN);
            return _edgePos[FSM_NEXT][node];
        }

        inline int nextEdgePos() {
            assert(_defaultAdjacentNode[FSM_NEXT] != INT_MIN);
            return _edgePos[FSM_NEXT][_defaultAdjacentNode[FSM_NEXT]];
        }

        inline Edge *nextEdge(int node) {
            assert(_edge[FSM_NEXT][node] != NULL);
            return _edge[FSM_NEXT][node];
        }

        inline Edge *findNextEdge(int node) {
            return _edge[FSM_NEXT][node];
        }

        inline Edge *nextEdge() {
            assert(_defaultAdjacentNode[FSM_NEXT] != INT_MIN);
            return _edge[FSM_NEXT][_defaultAdjacentNode[FSM_NEXT]];
        }

        /**
         * @return The reverse flow connections, useful to go backwards at
         * decoding time
         */
        inline ::vector<Node *> & prevNodes() {  return _adjacentNodes[FSM_PREV]; }

        inline void addPrevNode(int edgePos, Edge *edge, Node *node) {
            assert(node->id() >= 0 && node->id() < MAX_NUMPARSE_NODES);
            _adjacentNodes[FSM_PREV].push_back(node);
            addPrevEdge(edgePos, edge, node);
        }

        inline void addPrevEdge(int edgePos, Edge *edge, Node *node) {
            assert(node->id() <= MAX_NUMPARSE_NODES && edge->id() < MAX_NUMPARSE_EDGES && _edge[FSM_PREV][node->id()] == NULL);
            _edgePos[FSM_PREV][node->id()] = edgePos;
            _edge[FSM_PREV][node->id()] = edge;
            if(edge->hasSignal() && edge->strand() == strand())
                _defaultAdjacentNode[FSM_PREV] = node->id();

        }

        inline int prevEdgePos(int node) {
            assert(node >= 0 && node <= MAX_NUMPARSE_NODES && _edgePos[FSM_PREV][node] != INT_MIN);
            return _edgePos[FSM_PREV][node];
        }

        inline int prevEdgePos() {
            assert(_defaultAdjacentNode[FSM_PREV] != INT_MIN);
            return _edgePos[FSM_PREV][_defaultAdjacentNode[FSM_PREV]];
        }

        inline Edge *prevEdge(int node) {
            assert(_edge[FSM_PREV][node] != NULL);
            return _edge[FSM_PREV][node];
        }

        inline Edge *findPrevEdge(int node) {
            return _edge[FSM_PREV][node];
        }

        inline Edge *prevEdge() {
            assert(_defaultAdjacentNode[FSM_PREV] != INT_MIN);
            return _edge[FSM_PREV][_defaultAdjacentNode[FSM_PREV]];
        }

        inline void setNumContexts(int numContexts = 1) {
            assert(!_length);
            _length = new Pair<int> [numContexts];
            for(int i = 0; i < numContexts; i++) {
                _length[i].f = INT_MIN;
                _length[i].s = INT_MIN;
            }
        }

        inline int minLength(int context = 0) {
            assert(_length[context].f != INT_MIN);
            return _length[context].f;
        }
        inline int maxLength(int context = 0) {
            assert(_length[context].s != INT_MIN);
            return _length[context].s;
        }
        inline void setMinLength(int length = 1, int context = 0) {
            assert(_length[context].f == INT_MIN);
            _length[context].f = length;
        }
        inline void setMaxLength(int length = 1, int context = 0) {
            assert(_length[context].s == INT_MIN);
            _length[context].s = length;
        }

        virtual ~Node() {
            if(_length)
                delete [] _length;
        }

    };

    /**
     * FixedLenNode is a subclass of Node that represents a node or a state
     * in a FSM object which a fixed length
     */
    class FixedLenNode : public Node {
    public:
    FixedLenNode(TParseNode id,
                 std::string name,
                 TNodeId2 id2,
                 TNodeId3 id3,
                 int numPhases,
                 int maxInPhases,
                 int maxOutPhases,
                 TStrand strand
        ) : Node(id, name,
                 id2, id3,
                 numPhases,
                 maxInPhases, maxOutPhases,
                 strand) {

            _nodeType = FIXEDLEN_NODE;

        }

        ~FixedLenNode() { }

    };

    /**
     * VarLenNode is a subclass of Node that represents a node or a state
     * in a FSM object which has variable length.
     */
    class VarLenNode : public Node {
    public:
    VarLenNode(TParseNode id,
               std::string name,
               TNodeId2 id2,
               TNodeId3 id3,
               int numPhases,
               int maxInPhases,
               int maxOutPhases,
               TStrand strand
        ) : Node(id, name,
                 id2, id3,
                 numPhases, maxInPhases,
                 maxOutPhases,
                 strand) {

            _nodeType = VARLEN_NODE;

        }
        ~VarLenNode() { }
    };

    /**
     * VarPhaseNode is a subclass of Node that represents a node or a state
     * in a FSM object which has variable length and variable phase. These are
     * the only class of Nodes which can be disrupted by Signal phase
     * disruptors.
     */
    class VarPhaseNode : public Node {
        TParseEdge _phaseDisruptor;
    public:
    VarPhaseNode(TParseNode id,
                 std::string name,
                 TNodeId2 id2,
                 TNodeId3 id3,
                 int numPhases,
                 int maxInPhases,
                 int maxOutPhases,
                 TStrand strand,
                 TParseEdge disruptor
        ) : Node(id, name,
                 id2, id3,
                 numPhases,
                 maxInPhases,
                 maxOutPhases,
                 strand) {

            _nodeType = VARPHASE_NODE;
            _phaseDisruptor = disruptor;
        }

        /**
         * @param beg the beginning of this node
         * @param phase the phase at the beginning of this node
         * @param strand the strand in which the node is defined
         * @return the phase which disrupts the continuity of this node. If
         * a signal which is a phase disruptor for this node is found at any time
         * within the node and has the same phase as the returning value of this
         * function then the node's end cannot be beyond the disrupting signal's
         * position.
         */
        inline int disruptingPhase(int beg, int phase, TStrand strand) {
            int p = (beg  + numPhases() - phase) % numPhases();
            if(strand == STRAND_COMP)
                p = (p + 2) % numPhases();
            return p;
        }

        /**
         * @return the id of the signal which disrupts this node
         */
        inline TParseEdge phaseDisruptor() {
            return _phaseDisruptor;
        }

        ~VarPhaseNode() {   }
    };

    class Word {
    protected:
        TParseWord _id;
        std::string _name;
        int _numPhases;       //!< Number of phases of Word
        int _maxInPhases;     //!< Max number of phases at the beginning of Word
        int _maxOutPhases;    //!< Max number of phases at the end of Word
        TStrand _strand;
        TNodeId2 _prevNode, _nextNode;
        Word *_complementWord;

    public:
        Word(TParseWord id,
             std::string name,
             TNodeId2 prevNode,
             TNodeId2 nextNode,
             int numPhases,
             int maxInPhases,
             int maxOutPhases,
             TStrand strand) {

            this->_id = id;
            this->_name = name;
            _prevNode = prevNode;
            _nextNode = nextNode;
            this->_numPhases = numPhases;
            this->_maxInPhases = maxInPhases;
            this->_maxOutPhases = maxOutPhases;
            this->_strand = strand;
        }

        inline TParseWord id() {  return _id; }
        inline std::string & name() { return _name; }
        inline TStrand strand() {  return _strand; }
        inline void setStrand(TStrand strand) {  _strand = strand; }
        inline int numPhases() { return _numPhases; }
        inline int maxInPhases() { return _maxInPhases; }
        inline int maxOutPhases() { return _maxOutPhases; }
        inline Word *complementWord() { return _complementWord; }
        inline TNodeId2 nextNode() {  return _nextNode; }

        inline void setComplementWord(Word *compWord) {
            assert(compWord != NULL);
            _complementWord = compWord;
        }
        inline TNodeId2 prevNode() {  return _prevNode; }

        ~Word() { }
    };

    /**
     * The FSM class is a subtype of resource which models a finite state
     * automata with variable-length and/or variable-phase nodes.
     * It allows for edges to optionally be associated with signals at
     * positions where the patterns that represent the signals are found in the
     * sequence. The edges and nodes defined by this class are used within the
     * decoder to retrieve important information for nodes (lengths, phases,
     * next and previous nodes, next and previous signals, and so on) and for
     * edges (signal information, self-jumping, and so on)
     ***************************************************************************/

    class FSM : Resource {
    protected:
        int numContexts;
        std::string fsmContents;
        Edge *_edges[MAX_NUMPARSE_EDGES];  //!< list of edges
        Node *_nodes[MAX_NUMPARSE_NODES];  //!< list of nodes
        Word *_words[MAX_NUMPARSE_WORDS];  //!< list of words
        int _numParseNodes, _numParseEdges, _numParseWords;
        std::map<std::string, Edge *> namedEdges;
        std::map<std::string, Node *> namedNodes;
        std::map<std::string, Word *> namedWords;
        TStrand parseStrand;

    public:
        FSM(TStrand = BOTH_STRANDS);
        FSM(::ifstream &  fd, int numContexts, TStrand =  BOTH_STRANDS);
        FSM(std::vector<std::string> &params, int &offset, ResourceEngine *re);

        inline int  numParseWords() {
            return _numParseWords;
        }

        inline int  numParseNodes() {
            return _numParseNodes;
        }

        inline int numParseEdges() {
            return _numParseEdges;
        }

        inline TStrand parsingStrand() {
            return parseStrand;
        }

        inline void saveHeader(std::ofstream &fd) {
            Resource::saveHeader(fd);
            fd << " null " << numContexts << " " << Utils::tstrandToString(parseStrand) << endl;
        }

        inline  void saveContents(std::ofstream & fd) {
            fd << fsmContents;
        }

        inline Word **words() {  return _words; }
        inline Node **nodes() {  return _nodes; }
        inline Edge **edges() {  return _edges; }

        inline void removeSyncEndNode(TParseNode st) {
            _nodes[st]->setSyncEndFlag(false);
        }

        inline Word *word(TParseWord w) {
            assert(w >= 0 && w <= numParseWords() + 1);
            return _words[w];
        }

        inline Node *node(int nd) {
            return node((TParseNode)nd);
        }

        inline Node *node(TParseNode nd) {
            assert(nd >= 0 && nd <= numParseNodes() + 1);
            return _nodes[nd];
        }

        inline Word *word(TParseWord w, TStrand strand) {
            if(_words[w]->strand() == strand)
                return _words[w];
            return _words[w]->complementWord();
        }

        inline Node *node(TParseNode nd, TStrand strand) {
            if(_nodes[nd]->strand() == strand)
                return _nodes[nd];
            return _nodes[nd]->complementNode();
        }

        inline Word *word(const char * wordName) {
            std::string s(wordName);
            return word(s);
        }

        inline Node *node(const char * nodeName) {
            std::string s(nodeName);
            return node(s);
        }

        inline Edge *edge(TParseEdge ed) {
            assert(ed >= 0 && ed <= numParseEdges());
            return _edges[ed];
        }

        inline Edge *edge(TParseNode prev, TParseNode next) {
            return node(prev)->nextEdge(next);
        }

        inline Edge *findEdge(TParseNode prev, TParseNode next) {
            return node(prev)->findNextEdge(next);
        }

        inline Edge *edge(TParseEdge e, TStrand strand) {
            if(_edges[e]->strand() == strand)
                return _edges[e];
            return _edges[e]->complementEdge();
        }

        inline Edge *edge(const char *edgeName) {
            std::string s(edgeName);
            return edge(s);
        }

        inline vector<Node *> & retrieveSyncBegStates() {
            return node(SYNC_BEG_STATE)->nextNodes();
        }

        inline void restoreSyncBegStates(vector<Node *> &nextNodes) {
            node(SYNC_BEG_STATE)->setNextNodes(nextNodes);
        }

        void retrieveContents(std::ifstream &);
        void createEdge(boost::RegEx & r);
        void createNode(boost::RegEx & r, TNodeType t, int numContexts);
        void createWord(boost::RegEx & r);
        void createTransition(boost::RegEx & r);
        void restrictSyncBegStates(TParseNode nsyncNode);

        void removeEdge(TParseNode prevNode, TParseNode nextNode);
        void removeId2Node(TNodeId2 node);

        inline Edge *edge(const std::string & edgeName) {
            return edge((std::string &)edgeName);
        }
        Edge *edge(std::string & edgeName);
        Edge *findEdge(std::string & edgeName);

        inline Word *word(const std::string & wordName) {
            return word((std::string &)wordName);
        }

        inline Node *node(const std::string & nodeName) {
            return node((std::string &)nodeName);
        }

        Word *word(std::string & wordName);
        Node *node(std::string & nodeName);

        Word *findWord(std::string & wordName);
        Node *findNode(std::string & nodeName);

        void setParseStrand(TStrand strand);

        virtual void storeInfo(std::ofstream & fd) {;}
        virtual void retrieveInfo(std::ifstream & fd) {;}

        virtual ~FSM();

    };

}

#endif
