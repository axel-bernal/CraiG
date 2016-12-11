#include "TypeDefs.h"
#include "GenLibExcept.h"
#include <map>

/****************************************************************************
* TypeDefs.cpp - part of the craig namespace, a genomics library
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


static std::map<std::string, TLossType> lossFunctions;
static std::map<std::string, TTrainMethod> trainMethods;
static std::map<std::string, TParseNode> parseNodes;
static std::map<std::string, TNodeId2>   id2Nodes;
static std::map<std::string, TNodeId3>   id3Nodes;
static std::map<std::string, TEdgeId2> id2Edges;
static std::map<std::string, TMultiUpd> multiUpdMethods;
static std::map<std::string, TOracleUpd> oracleUpdMethods;
static std::map<std::string, TAvgMethod> avgMethods;
static std::map<std::string, TCombMethod> combMethods;
static std::map<TEdgeId2, std::string> id2EdgeNames;
using namespace lless;

TypeDefs TypeDefs::initializer;

TypeDefs::TypeDefs() {

  //averaging methods
  avgMethods["none"] = AVG_NONE;
  avgMethods["all"] = AVG_ALL;
  avgMethods["last"] = AVG_LAST;

  // combination methods
  combMethods["euclid"] = COMB_EUCLID; 
  combMethods["k-l"] = COMB_KL;
  
  //multilabel update methods
  multiUpdMethods["ALL"] = ML_ALL;
  multiUpdMethods["MIN"] = ML_MIN;
  multiUpdMethods["EXP"] = ML_EXP;
  multiUpdMethods["SEP"] = ML_SEP;
  multiUpdMethods["LONGEST"] = ML_LONGEST;
  multiUpdMethods["MINLONGER"] = ML_MINLONGER;
  
  // orcale updates
  oracleUpdMethods["TOP"] = OC_TOP;
  oracleUpdMethods["ALL"] = OC_ALL;
  oracleUpdMethods["BETTER"] = OC_BETTER;
  oracleUpdMethods["NONE"] = OC_NONE;

  //loss functions

  lossFunctions["NONE"] = LF_NONE;
  lossFunctions["EDGE"] = LF_EDGE;
  lossFunctions["SOFT_EDGE"] = LF_SOFT_EDGE;
  lossFunctions["HAMMING"] = LF_HAMMING;
  lossFunctions["SEGMENT"] = LF_SEGMENT;
  lossFunctions["CORR_COEF"] = LF_CORR_COEF;
  lossFunctions["ZERO_ONE"] = LF_ZERO_ONE;
  lossFunctions["F_SCORE"] = LF_FSCORE;
  
  // id3Nodes
  id3Nodes["INTRON"] = INTRON;
  id3Nodes["INTERGENIC"] = INTERGENIC;
  id3Nodes["EXON"] = EXON;
  id3Nodes["UTR"] = UTR;
  id3Nodes["NO_NODE_INST"] = NO_NODE_INST;
  
  // id2Nodes
  id2Nodes["ANY_INTRON"] = ANY_INTRON;
  id2Nodes["ANY_INTERGENIC"] = ANY_INTERGENIC;
  id2Nodes["ANY_UTR_INTRON"] = ANY_UTR_INTRON;
  id2Nodes["ANY_5UTR"] = ANY_5UTR;
  id2Nodes["ANY_3UTR"] = ANY_3UTR;
  id2Nodes["ANY_UTR"] = ANY_UTR;
  id2Nodes["INIT_EXON"] = INIT_EXON;
  id2Nodes["INTERNAL_EXON"] = INTERNAL_EXON;
  id2Nodes["LAST_EXON"] = LAST_EXON;
  id2Nodes["SINGLE_EXON"] = SINGLE_EXON;
  id2Nodes["SYNC_STATE"] = SYNC_STATE;
  
  id2Edges["START"] = START;
  id2Edges["STOP"] = STOP;
  id2Edges["DONOR"] = DONOR;
  id2Edges["ACCEPTOR"] = ACCEPTOR;
  id2Edges["TSS"] = TSS;
  id2Edges["PAS"] = PAS;
  id2Edges["NO_EDGE_INST"] = NO_EDGE_INST;

  id2EdgeNames[START] = "START";
  id2EdgeNames[STOP] = "STOP";
  id2EdgeNames[DONOR] = "DONOR";
  id2EdgeNames[ACCEPTOR] = "ACCEPTOR";
  id2EdgeNames[TSS] = "TSS";
  id2EdgeNames[PAS] = "PAS";
  id2EdgeNames[NO_EDGE_INST] = "NO_EDGE_INST";
  
  trainMethods["MIRA"] = MIRA;
  trainMethods["PEGASOS"] = PEGASOS;
  trainMethods["CWL"] = CWL;
  trainMethods["PERCEPTRON"] = PERCEPTRON;
  trainMethods["ARROW"] = ARROW;

}

TAvgMethod TypeDefs::stringToTAvgMethod(const std::string &s) {
  std::map<std::string, TAvgMethod>::iterator it;
  it = avgMethods.find(s);
  if(it == avgMethods.end())
    throw EXCEPTION( PARSE_ERROR, std::string("undefined TAvgMethod ")+s);
  return it->second;
}

TCombMethod TypeDefs::stringToTCombMethod(const std::string &s) {
  std::map<std::string, TCombMethod>::iterator it;
  it = combMethods.find(s);
  if(it == combMethods.end())
    throw EXCEPTION( PARSE_ERROR, std::string("undefined TCombMethod ")+s);
  return it->second;
}

TMultiUpd TypeDefs::stringToTMultiUpd(const std::string &s) {
  std::map<std::string, TMultiUpd>::iterator it;
  it = multiUpdMethods.find(s);
  if(it == multiUpdMethods.end())
    throw EXCEPTION( PARSE_ERROR, std::string("undefined TMultiUpd ")+s);
  return it->second;
}

TOracleUpd TypeDefs::stringToTOracleUpd(const std::string &s) {
  std::map<std::string, TOracleUpd>::iterator it;
  it = oracleUpdMethods.find(s);
  if(it == oracleUpdMethods.end())
    throw EXCEPTION( PARSE_ERROR, std::string("undefined TOracleUpd ")+s);
  return it->second;
}

TLossType TypeDefs::stringToTLossType(const std::string &s) {
  std::map<std::string, TLossType>::iterator it;
    it = lossFunctions.find(s);
    if(it == lossFunctions.end())
      throw EXCEPTION( PARSE_ERROR, std::string("undefined TLossType ")+s);
    return it->second;
}

TTrainMethod TypeDefs::stringToTTrainMethod(const std::string &s) {
  std::map<std::string, TTrainMethod>::iterator it;
    it = trainMethods.find(s);
    if(it == trainMethods.end())
      throw EXCEPTION( PARSE_ERROR, std::string("undefined TTrainMethod ")+s);
    return it->second;
}

TNodeId2 TypeDefs::stringToTNodeId2(const std::string &s) {
  std::map<std::string, TNodeId2>::iterator it;
    it = id2Nodes.find(s);
    if(it == id2Nodes.end())
      throw EXCEPTION( PARSE_ERROR, std::string("undefined TNodeId2 ")+s);
    return it->second;
}

TNodeId3 TypeDefs::stringToTNodeId3(const std::string &s) {
  std::map<std::string, TNodeId3>::iterator it;
  it = id3Nodes.find(s);
  if(it == id3Nodes.end())
    throw EXCEPTION( PARSE_ERROR, std::string("undefined TNodeId3 ")+s);
  return it->second;
}

TEdgeId2 TypeDefs::stringToTEdgeId2(const std::string &s) {
  std::map<std::string, TEdgeId2>::iterator it;
  it = id2Edges.find(s);
  if(it == id2Edges.end())
    throw EXCEPTION( PARSE_ERROR, std::string("undefined EdgeId2 ")+s);
  return it->second;
}

std::string & TypeDefs::tEdgeId2ToString(const TEdgeId2 id) {
  std::map<TEdgeId2, std::string>::iterator it;
  it = id2EdgeNames.find(id);
  if(it == id2EdgeNames.end())
    throw EXCEPTION( PARSE_ERROR, std::string("undefined EdgeId2 "));
  return it->second;
}
