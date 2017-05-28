/****************************************************************************
 * TypeDefs.h - part of the craig namespace, a genomics library
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


#ifndef _TYPE_DEFS_H_
#define _TYPE_DEFS_H_

#include <string>
#include "Utils.h"

/*
 * Constants for CRAIG
 */


#define MAX_NUM_EXONS 60
#define MAX_NUMPARSE_WORDS 10
#define MAX_NUMPARSE_EDGES 63
#define BITS4P_MAXNODELEN 20
#define BITS4P_ID2EDGES 4
#define BITS4P_EDGES 6
#define BITS4P_WORDS 1
#define MAX_NUMPARSE_NODES 63
#define BITS4P_NODES 6
#define BITS4P_ID3NODES 3
#define BITS4P_NODEPHASE 2
#define BITS4P_STRAND 1
#define BITS4P_TAG_CLASS 2
#define BITS4P_NUMEXONS 7
#define MAX_NUMPARSE_ELEMS ((MAX_NUMPARSE_NODES > MAX_NUMPARSE_EDGES) ? MAX_NUMPARSE_NODES : MAX_NUMPARSE_EDGES)
#define BITS4P_TAG (BITS4P_NODES > BITS4P_EDGES ? BITS4P_NODES : BITS4P_EDGES)

//! enumerate for different types of parsing nodes
typedef enum {
    SYNC_BEG_STATE=0,
    SYNC_END_STATE=MAX_NUMPARSE_NODES-1,
    INVALID_NODE=MAX_NUMPARSE_NODES
} TParseNode;

//! enumerate for different types of parsing edges
typedef enum {
    INVALID_EDGE=MAX_NUMPARSE_EDGES
} TParseEdge;

//! enumerate for different types of parsing words, for global features
typedef enum {
    INVALID_WORD=MAX_NUMPARSE_WORDS
} TParseWord;


#define NUM_BIOTR_STATES 10
//! enumerate for different types of secondary id for nodes.
typedef enum {
    ANY_INTRON,
    ANY_INTERGENIC,
    ANY_UTR_INTRON,
    ANY_5UTR,
    ANY_3UTR,
    ANY_UTR,
    INIT_EXON,
    INTERNAL_EXON,
    LAST_EXON,
    SINGLE_EXON,
    SYNC_STATE
} TNodeId2;

#define NUM_NODE_INSTS 4
//! enumerate for different types of tertiary id for nodes
typedef enum {
    INTRON,
    INTERGENIC,
    EXON,
    UTR,
    NO_NODE_INST
} TNodeId3;

// Biosignal-related type definitions
//! enumerate for different types of secondary id for edges
#define NUM_EDGE_INSTS 6
typedef enum {
    START,
    STOP,
    DONOR,
    ACCEPTOR,
    TSS,
    PAS,
    NO_EDGE_INST
} TEdgeId2;

//! defined type for enumerating the available training methods
typedef enum {
    PERCEPTRON,
    MIRA,
    PEGASOS,
    CWL,
    ARROW
} TTrainMethod;

/**
 * Different Types of loss functions. This has to be updated
 * if newer types are to be added. HAMMING definition must be
 * presentin all of them
 */

typedef enum {
    LF_NONE,
    LF_EDGE,
    LF_SOFT_EDGE,
    LF_HAMMING,
    LF_SEGMENT,
    LF_CORR_COEF,
    LF_ZERO_ONE,
    LF_FSCORE
} TLossType;

typedef enum {
    ML_ALL,
    ML_MIN,
    ML_EXP,
    ML_SEP,
    ML_LONGEST,
    ML_MINLONGER
} TMultiUpd;

typedef enum {
    OC_TOP,
    OC_ALL,
    OC_BETTER,
    OC_NONE
} TOracleUpd;

typedef enum {
    AVG_NONE,
    AVG_ALL,
    AVG_LAST
} TAvgMethod;

typedef enum {
    COMB_KL,
    COMB_EUCLID,
} TCombMethod;

/**
 * The TypeDefs class contains all the type definitions and constant
 * declarations needed by CRAIG.
 ****************************************************************************/

class TypeDefs  {
private:
    static TypeDefs initializer;
public:
    TypeDefs();

    static TAvgMethod stringToTAvgMethod(const std::string &s);
    static TCombMethod stringToTCombMethod(const std::string &s);
    static TMultiUpd stringToTMultiUpd(const std::string &s);
    static TOracleUpd stringToTOracleUpd(const std::string &s);
    static TLossType stringToTLossType(const std::string &s);
    static TTrainMethod stringToTTrainMethod(const std::string &s);
    static TParseNode stringToTParseNode(const std::string &s);
    static TNodeId2 stringToTNodeId2(const std::string &s);
    static TNodeId3 stringToTNodeId3(const std::string &s);
    static TEdgeId2 stringToTEdgeId2(const std::string &s);
    static std::string &tEdgeId2ToString(const TEdgeId2 id);
    ~TypeDefs() {}
};

typedef int TFilter;
typedef int TFeature;


#endif
