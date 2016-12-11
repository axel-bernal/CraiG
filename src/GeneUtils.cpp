#include <iostream>
#include <streambuf>
#include <fstream>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "GeneUtils.h"
#include "Gene.h"
#include "FeatureEngine.h"
#include "NodeInst.h"
#include "FSM.h"
#include <boost/regex.hpp>

/****************************************************************************
 * GeneUtils.cpp - part of the craig namespace, a genomics library
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

/**
 * A namespace created for all classes which are part of craig, a genomics 
 * library.
 */
namespace craig {

  static boost::RegEx rExLocation1("([^\\[]+)\\[(\\S*)\\s*$");
  static boost::RegEx rExLocation2("([^\\]]*)\\](\\S+)\\s*$");
  static boost::RegEx rExLocation3("(\\d+)?<(\\S*)\\s*$");
  static boost::RegEx rExLocation4("^(\\S*)>");
  static boost::RegEx rExSep("\\s*\\;\\s*");
  static boost::RegEx rExLoc("^(\\S+)_(\\d+)_(\\d+)\\s*$");

  void GeneUtils::loadIGenicRegions(map<string, list<BioFeature> > &igenics, 
				    string & igenics_fn) {

    map<string, list<BioFeature> >::iterator igit;
    ifstream igStream(igenics_fn.c_str());
    if(igStream) {
      string line;
      while(std::getline(igStream, line), !igStream.eof()) {
	string cid, id, location;
	vector<Pair<int> > ve;
	std::istringstream istr(line);
	istr >> id >> location;
	id.erase(0); // removing the >
	GeneUtils::readExonLocations(location, cid, ve);
	
	if((igit = igenics.find(cid)) == igenics.end())
	  igenics[cid] = list<BioFeature>();
	
	igenics[cid].push_back(BioFeature(id, NULL, ve[0].f, ve[0].s, 
					  true, STRAND_FWD));
      }
    }
  }
  
  void GeneUtils::getGenicRegions(list<BioFeature> &igenics, 
				  list<BioFeature> &genics,
				  int percent_igenic) {

    igenics.sort();
    list<BioFeature>::iterator it = igenics.begin();
    list<BioFeature>::iterator igenicL = it++;
    int cend = -1;
    for( ; it != igenics.end(); it++) {
      list<BioFeature>::iterator igenicR = it;

      if(igenicR->getStrand() != igenicL->getStrand())
	throw EXCEPTION(BAD_USAGE, " intergenics are in one strand");

      int cbeg = cend < 0 ? percent_igenic*(igenicL->begin() + igenicL->end())/100 + 1 :
	cend + 1;
      cend = percent_igenic*(igenicR->begin() + igenicR->end())/100 - 1;
      BioFeature nbf = *igenicR;
      nbf.setBegin(cbeg); nbf.setEnd(cend);
      genics.push_back(nbf);
      igenicL = igenicR;
    }
  }
  

  void GeneUtils::genes2transcripts(Sequence &c, TSetType typeSet,
				    list<Transcript *> &transcripts) {

    list<Gene *> &listGenes = (list<Gene *> &)c.getAnnotBioFeats(typeSet);
    list<Gene *>::iterator git = listGenes.begin();
    for( ; git != listGenes.end(); git++) {
      Gene &gene = *(*git);
      for(int j = 0; j < gene.transcripts().size(); j++) {
	Transcript &t = gene.transcripts()[j];
	transcripts.push_back(&t);
      }
    }
    
    transcripts.sort(greater<Transcript *>());
  }

  /**
   * Loads gene headers and sequences in locs format. This is the only format
   * which craig is able to read for the time being.
   */
  void GeneUtils::loadGenes(list<Gene> & listGenes,
                           const char *geneFile, 
                           bool buildstartInfo, 
                           bool readfastaSequence, 
                           list<Sequence *> &annotSeqs, 
                           Sigma *sigma, 
                           TSetType typeSet) {

    ::ifstream geneStream(geneFile);

    if(!geneStream) {
      assert(0);
      throw EXCEPTION(BIOFEATS_UNAVAILABLE, string(geneFile));
    }
    
    std::string annotSeq = "", geneId, transcriptId;
    std::string locationHeader, id;
    std::string exonSeq = "", line;
    std::ifstream &fd = geneStream;
    TStrand strand;
    Sequence *annotSeqObj;
    
    Gene *geneObj;
    bool locationIsPresent = false;
    int tstats[3] = {0,0,0};
    int altTranscripts = 0;

    /*
     * vectors where we will load gene exons
     */
    vector<Pair<int> > ve5;
    vector<Pair<int> > ve3;
    vector<Pair<int> > ve; 

    string::size_type begTransId;
    
    fd >> id;

    do {
      if(fd.eof())
        break;

      bool hasStart = true, hasStop = true, hasPeptide = false;
      int phase = 0;
      strand = STRAND_FWD;

      if(id[0] != '>') {
        assert(0);
        throw EXCEPTION(BAD_LOCATION_FORMAT, id + " must start with >");
      }
      
      if(!id.compare(1, 4, "PEP-")) {
        hasPeptide = true;
      }

      fd >> locationHeader;
      std::string remnant;
      std::getline(fd, remnant);
      locationHeader += remnant;
      locationIsPresent = (locationHeader.length() != 0);
      /*
       * Check if it's an alternative splice
       */
      begTransId = id.rfind("|");
      
      if(begTransId != std::string::npos) {
        geneId = id.substr(1, begTransId - 1);
        transcriptId = id.substr(begTransId + 1);
      }
      else {
        geneId = id.substr(1);
        transcriptId = geneId;
      }

      if(readfastaSequence) { 
        /*
         * Reading the fasta sequence that follows the fasta header
         */
        exonSeq.clear();
        fd >> line;

        while(!fd.eof() && line.c_str()[0] != '>') {
          exonSeq = exonSeq + line;
          std::getline(fd, line);
        }

        id = line;

        if(exonSeq.length() == 0)
          continue;
      }
      else 
        fd >> id;

      if(locationIsPresent)
        GeneUtils::parseSequenceLocation(locationHeader, annotSeq, ve5, ve, ve3, hasStart, hasStop, phase);
      else { 
        /*
         * Location not present. Take default
         */
        if(!readfastaSequence) { 
          /*
           * if sequence is not present there is a format error
           */
          assert(0);
          throw EXCEPTION(CONTIG_UNAVAILABLE, std::string(id));
        }      

        Pair<int> pi(1, exonSeq.length());
        ve.push_back(pi);
        annotSeq = geneId;
      }

      int gene_start = ve5.size() ? ve5[0].f : 
	ve.size() ? ve[0].f : 
	ve3.size() ? ve3[0].f : 0;
      int gene_end = ve3.size() ? ve3[ve3.size() - 1].s :
	ve.size() ? ve[ve.size() - 1].s : 
	ve5.size() ? ve5[ve5.size() - 1].s : 0;

      assert(gene_start && gene_end);

      if(gene_end < gene_start)
	strand = STRAND_COMP;    
      
      if(readfastaSequence) { 
        /*
         * We create a new annotSeq for each sequence
         */
        sigma->checkSeq((char *)exonSeq.c_str());
        annotSeqObj = new Sequence(annotSeq, (char *)exonSeq.c_str(), false);
        annotSeqs.push_back(annotSeqObj);
      }
      else  // find the annotSeq entry
        annotSeqObj = (Sequence *)BasicSeq::findSequence(annotSeq, (list<BasicSeq *> &)annotSeqs);

      if(!annotSeqObj) {
        assert(0);
        throw EXCEPTION(CONTIG_UNAVAILABLE, string(annotSeq));
      }

      /*
       * Creating the geneObj
       */
      geneObj = Gene::findGene(geneId, listGenes);
      
      if(!geneObj) {
        // genes are loaded in OneStrand format by default
        listGenes.push_front(Gene(geneId, annotSeqObj, true, strand));
        geneObj = &(listGenes.front());
        annotSeqObj->addBioFeature(&(*geneObj), typeSet);
      }
      
      if(geneObj->getSequence()->id().compare(annotSeqObj->id())) {
        assert(0);
        throw EXCEPTION(INCOMPAT_TRANSCRIPT, geneObj->getSequence()->id());
      }

      /*
       * Creating a new transcript
       */
      Transcript tmp(transcriptId, geneObj, hasStart, hasStop, hasPeptide, phase);
      Transcript &t = geneObj->addTranscript(tmp);

      /*
       * Add Exons
       */
      unsigned int i, currentLen = 0;
      for(i = 0; i < ve5.size(); i++)  {
        if(strand == STRAND_COMP)
          ve5[i].updVals(ve5[i].s, ve5[i].f);
        t.add5Exon(ve5[i].f, ve5[i].s);
      }

      for(unsigned int i=0; i < ve.size(); i++)  {
        int exonPhase = currentLen % 3;

        if(strand == STRAND_COMP)
          ve[i].updVals(ve[i].s, ve[i].f);

        //        cerr << geneId << " " << transcriptId << " " << strand << " " << annotSeq << " " << i << " " << ve[i].f << " " << ve[i].s << endl;
        currentLen += ve[i].s - ve[i].f + 1;
        t.addExon(ve[i].f, ve[i].s, exonPhase);
      }

      for(i = 0; i < ve3.size(); i++)  {
        if(strand == STRAND_COMP)
          ve3[i].updVals(ve3[i].s, ve3[i].f);

        t.add3Exon(ve3[i].f, ve3[i].s);
      }
      
      if(!t.hasStart())
        tstats[1]++;
      else if(!t.hasStop())
        tstats[2]++;
      else
        tstats[0]++;

      if(geneObj->transcripts().size() == 2)
        altTranscripts += 2;
      else if(geneObj->transcripts().size() > 2)
        altTranscripts++;
      ve.clear();
      ve5.clear();
      ve3.clear();
    } while(1);

    cerr  << listGenes.size() << " Genes loaded for Set " << typeSet << ".....\n";
    cerr  << "\t" << tstats[0] << " Complete transcripts\n";
    cerr  << "\t" << tstats[1] << " 5' truncated transcripts\n";
    cerr  << "\t" << tstats[2] << " 3' truncated transcripts\n";
    cerr  << "\t" << altTranscripts << " alternative transcripts\n";
  }

  
  void GeneUtils::extractFrameOrfs(vector<pair<int,int> > &orfs,
				   Sequence &c, 
				   TypedFilter<EdgeInst> **signals,
				   int frame, TStrand strand,
				   int mingeneLen,
				   bool reportTruncated,
				   bool useAllStarts)  {
    
    assert(frame > 0 );
    vector<int> stops;
    int i, lenSeq = c.length(strand);
    char *tmp = NULL, *seq = c.getSeq(strand);
    EdgeInst *sig = NULL;

    if(reportTruncated)
      stops.push_back(frame - 3);

    for(i = frame; i <= lenSeq; i += 3) {
      sig = &signals[STOP]->value(i, strand);
      if(sig && sig->getType() == STOP) 
	stops.push_back(i); 
    }

    if(c.isCircular()) {
      tmp = new char[2*lenSeq + 2];
      strcpy(tmp, seq);
      strcat(tmp,seq);
      seq = tmp;
    } 
    
    while(i <= lenSeq) {
      sig = &signals[STOP]->value((i - 1) % c.length(strand) + 1, strand);
      if(sig && sig->getType() == STOP) {
	break;
      }
      i += 3;
    }
    if(i > c.length(strand) + 3 || reportTruncated)
      stops.push_back(i);

    for(i = 0; i < (int)stops.size() - 1; i++) {
      vector<int> starts;

      for(int j = stops[i] + 3; j < stops[i + 1]; j += 3) {
	sig = &signals[START]->value((j - 1) % c.length(strand) + 1, strand);
	if(sig && sig->getType() == START) {
	  starts.push_back(j);
	  if(!useAllStarts)
	    break;
	}
      }
      
      //      cerr << c.id() << " " << frame << " " << strand << " " << stops[i] << " " << stops[i+1] << " " << c.length(strand) << "\n";
      //      if(stops[i+1] >= c.length(strand) && strand == STRAND_FWD)
      //	cerr << "bad stop\n";

      if(!starts.size())
        continue;

      for(int j = 0; j < starts.size(); j++) {
	int start = starts[j];
	int geneLen = (strand == STRAND_FWD) ? stops[i+1] - start :
	  start - stops[i];
	
	if(geneLen <= mingeneLen) 
	  continue;

	int end = (strand == STRAND_FWD) ? stops[i+1] - 1 :
	  end = stops[i] + 1;

	if(start > c.length(strand))
	  start -= c.length(strand);
	
	if(end > c.length(strand))
	  end = (end <= c.length(strand) + 3) ? c.length(strand) :
	    end - c.length(strand);
	if(end <= 0) end = 1;
	
	orfs.push_back(pair<int,int>(start, end));
      }
    }

    if(tmp)
      delete [] tmp;

  }


  void GeneUtils::findNextGene(Sequence &c, 
			       SeqTags &st, SeqTags::iterator it, 
			       int &beg, int &end, TStrand &strand) {
    
    // find next gene's strand and boundaries
    beg = -1;
    end = -1;
    strand = NO_STRAND;
    bool gene_was_observed = false;
    for(; it != st.end(); it++) {
      if((*it)->getGEClass() == NODE_INST) {
	NodeInst *b = (NodeInst *)(*it);
 	
	if(b->getType() == INTERGENIC) {
	  if(!gene_was_observed)
	    continue;

	  if(strand == STRAND_FWD)
	    end = b->getPos();
	  else
	    beg = c.length() - b->getPos() + 1;
	  break;
	}
	else if(b->getType() == EXON || b->getType() == UTR || 
		b->getType() == INTRON) {
	  if(b->getStrand() == STRAND_FWD && beg < 0) {
	    beg = b->getPos();
	    strand = STRAND_FWD;
	  }
	  if(b->getStrand() == STRAND_COMP && end < 0) {
	    end = c.length() - b->getPos() + 1;
	    strand = STRAND_COMP;
	  }
	  gene_was_observed = true;
	}
      }
    }
    if(gene_was_observed) {
      end = (end < 0 ? c.length() : end);
      beg = (beg < 0 ? 1 : beg);
    }
  }    

  void GeneUtils::updateTranscript(Transcript *t,
				   bool hasStart, bool hasStop,
				   bool hasPeptide, int phase) {
    
    t->setStart(hasStart);
    t->setStop(hasStop);
    t->setPeptide(hasPeptide);
    
    if(t->getStrand() == STRAND_COMP)
      t->setPhase((3 - phase) % 3);
  }

  /**
   * Computes listGenes, the set of genes associated with the list of Tags
   * listTags. Genes contain at most one transcript and are in OneStrand fmt
   */
  void GeneUtils::tags2Genes(Sequence &c,
                             SeqTags & listTags,
                             const char *idPrefix,
                             FSM & fsm,
                             list<Gene> & listGenes,
                             int offset,
			     bool reportNonCodingGenes) {
    
    NodeInst *b;
    EdgeInst *sig = NULL;
    TEdgeId2 sigType = NO_EDGE_INST;;

    if(!listTags.size())
      return;

    Transcript transcript;
    int phase = listTags.initPhase();    
    Transcript *t = NULL;
    bool hasStart = false, hasStop = false, hasPeptide = false;
    SeqTags::iterator it;
    static int gcounter = 1, tcounter = 1;
    std::ostringstream geneId2, transcriptId2;
    geneId2 << idPrefix << "gen_" << gcounter++;
    transcriptId2 << "tr_" << tcounter++;
    
    std::string geneId = geneId2.str(), transcriptId = transcriptId2.str();
    geneId2.str(""); transcriptId2.str("");  

#ifndef NDEBUG
    cerr << transcriptId << endl;
    cerr << "init phase " << phase << endl;
    int length = 0;
#endif

    for(it = listTags.begin(); it != listTags.end();it++) {

      if((*it)->getGEClass() == NODE_INST) {
        b = (NodeInst *)(*it);
        Node *node = fsm.node((TParseNode)b->getParseType());
        int phBreak = listTags.phaseBreak(*it);

        if(phBreak > 0) // phaseBreak only could happen in between genes
          phase = phBreak;
#ifndef NDEBUG
	cerr << node->name() << b->getPos() << " " << b->getLen() << " " << phBreak << endl;
#endif
        if(b->getType() == INTERGENIC) {	  
	  if(t && 
	     (t->exons().size() || 
	      (reportNonCodingGenes && 
	       (t->TPexons().size() ||t->FPexons().size())))) {
	    
	    /*
	     * All exons from a transcript have been read. Complete gene info
	     * and reinitialize Transcript variable t
	     */
	    updateTranscript(t, hasStart, hasStop, hasPeptide, phase);
	    
	    geneId2 << idPrefix << "gen_" << gcounter++;
	    transcriptId2 << "tr_" << tcounter++;
            
	    geneId = geneId2.str();
	    transcriptId = transcriptId2.str();
	    geneId2.str(""); transcriptId2.str("");
          }
          t = NULL;
          hasStart = false;
          hasStop = false;
          hasPeptide = false;
#ifndef NDEBUG
          length = 0;
#endif
        }
        else if(b->getType() == EXON || b->getType() == UTR) {

          int exonBeg = b->getPos(), exonEnd = b->getPos() + b->getLen() - 1;
          
          if(!t) { // create a new transcript
            listGenes.push_back(Gene(geneId, &c, true, b->getStrand()));
            Gene & gene = listGenes.back();
            transcript = Transcript(transcriptId, &gene, false, false, false, 0);
            t = &gene.addTranscript(transcript);
            
            if(t->getStrand() == STRAND_FWD)
              t->setPhase(phase);

          }

	  if(b->getType() == EXON) {
	    t->addExon(exonBeg, exonEnd, 0);
	  }
	  else {
	    if(node->id2() == ANY_5UTR || node->id2() == ANY_UTR)  
	      t->add5Exon(exonBeg, exonEnd, 0);
	    else if(node->id2() == ANY_3UTR) 
	      t->add3Exon(exonBeg, exonEnd, 0);
	    else throw EXCEPTION(BAD_USAGE, "something really wrong\n");
	  }
#ifndef NDEBUG
	  length += b->getLen();
	  cerr << "transcript length = " << length << endl;
#endif
	}

        phase = node->nextPhase(phase, b->getLen());

#ifndef NDEBUG
        cerr << "phase = " << phase << endl;
#endif

      }
      else if((*it)->getGEClass() == EDGE_INST) {
        sig = (EdgeInst *)(*it);
        sigType = (TEdgeId2)sig->getType();
	//	cout << fsm.edge((TParseEdge)sig->getParseType())->name() << sig->getPos() << endl;
        if(sigType == START) {
          hasStart = true;
          std::string & name = fsm.edge((TParseEdge)sig->getParseType())->name();
          if(!name.compare("PEPSTART_F") || !name.compare("PEPSTART_B"))
            hasPeptide = true;
        }
        else if(sigType == STOP)
          hasStop = true;
      }
    }

    // updating last transcript if there was no downstream intergenic region
    if(t && 
       (t->exons().size() || 
	(reportNonCodingGenes && 
	 (t->TPexons().size() ||t->FPexons().size()))))
      updateTranscript(t, hasStart, hasStop, hasPeptide, phase);
    
#ifndef NDEBUG
    cerr << "final phase = " << phase << " " << listTags.endPhase() << endl;    
#endif

    assert(phase == listTags.endPhase());
  }

  void GeneUtils::updateAltTranscript(list<Gene> &listGenes, Gene &gene, 
                                      Transcript *t, bool hasStart,
                                      bool hasStop, bool hasPeptide,
                                      int phase) {
                                      

    updateTranscript(t, hasStart, hasStop, hasPeptide, phase);
    
    bool isAltSpliced = false;
    
    list<Gene>::iterator git = listGenes.begin();
    for( ; git != listGenes.end(); git++) {
      if(git->isAlternativeTranscript(*t)) {
        git->addTranscript(*t);
        isAltSpliced = true;
        break;
      }
    }
    
    if(!isAltSpliced) {
      listGenes.push_back(gene);
      (listGenes.back()).addTranscript(*t);
    }

  }

  /**
   * Computes tags, the list of tags associated with the list of Genes, 
   * listGenes. Genes 
   * can contain multiple transcripts. Overlapping transcripts are assumed to 
   * belong to the same gene as alternative transcripts.
   * Genes are in OneStrand format.
   */
  void GeneUtils::tags2Genes(Sequence &c, 
                             vector<SeqTags> & tags, 
                             const char *idPrefix,
                             FSM & fsm,
                             list<Gene> & listGenes,
                             int offset,
			     bool reportNonCodingGenes) {
    
    NodeInst *b;
    EdgeInst *sig = NULL;
    TEdgeId2 sigType = NO_EDGE_INST;

    if(!tags.size())
      return;

    SeqTags::iterator it, next_it;
    Gene gene;
    Transcript transcript;

    static int gcounter = 1, tcounter = 1;
    std::ostringstream geneId2, transcriptId2;
    geneId2 << idPrefix << "gen_" << gcounter++;
    transcriptId2 << "tr_" << tcounter++;
    
    std::string geneId = geneId2.str(), transcriptId = transcriptId2.str();
    geneId2.str(""); transcriptId2.str("");  

    for(unsigned int i = 0; i < tags.size(); i++) {

      Transcript *t = NULL;
      int phase = tags[i].initPhase();
      bool hasStart = false, hasStop = false, hasPeptide = false;
      SeqTags & seqTags = tags[i];
#ifndef NDEBUG
      cerr << transcriptId2 << endl;
      cerr << "init phase " << phase << endl;
      int length = 0;
#endif
      for(it = seqTags.begin(); it != seqTags.end(); it++) {

        if((*it)->getGEClass() == NODE_INST) {
          b = (NodeInst *)(*it);
          Node *node = fsm.node((TParseNode)b->getParseType());
          //          cerr << "phase @ " << b->getPos() << " " << phase << endl;
          int phBreak = seqTags.phaseBreak(*it);
          
          if(phBreak > 0)
            phase = phBreak;
	  
          if(b->getType() == INTERGENIC) {
	    if(t && 
	       (t->exons().size() || 
		(reportNonCodingGenes && 
		 (t->TPexons().size() || t->FPexons().size())))) {
	      
	      updateAltTranscript(listGenes, gene, t,
				  hasStart, hasStop, 
				  hasPeptide, phase);
	      
	      geneId2 << idPrefix << "gen_" << gcounter++;
	      transcriptId2 << "tr_" << tcounter++;
	      
	      geneId = geneId2.str();
	      transcriptId = transcriptId2.str();
	      geneId2.str(""); transcriptId2.str("");
            }

            t = NULL;
            hasStart = false;
            hasStop = false;
            hasPeptide = false;

          }
          else if(b->getType() == EXON || b->getType() == UTR) {
#ifndef NDEBUG
            cerr  << fsm.node((TParseNode)b->getParseType())->name() << " " << b->getPos() << " " << b->getLen() << endl;
#endif
            int exonBeg = b->getPos(), exonEnd = b->getPos() + b->getLen() - 1;
            
            if(!t) {
              gene = Gene(geneId, &c, true, b->getStrand());
              transcript = Transcript(transcriptId, &gene, 
                                      false, false, false, 0);
              t = &transcript;

              if(t->getStrand() == STRAND_FWD)
                t->setPhase(phase);

            }
            
	    if(b->getType() == EXON)
	      t->addExon(exonBeg, exonEnd, 0);
	    else {
	      if(node->id2() == ANY_5UTR || node->id2() == ANY_UTR) 
		t->add5Exon(exonBeg, exonEnd, 0);
	      else if(node->id2() == ANY_3UTR)  
		t->add3Exon(exonBeg, exonEnd, 0);
	      else throw EXCEPTION(BAD_USAGE, "something really wrong\n");
	    }

#ifndef NDEBUG
            length += b->getLen();
            cerr << "coding length = " << length << endl;
#endif
          }
          
          phase = node->nextPhase(phase, b->getLen());
          
#ifndef NDEBUG
          cerr << "phase = " << phase << endl;
#endif
          
        }
        else if((*it)->getGEClass() == EDGE_INST) {
#ifndef NDEBUG
          cerr  << fsm.edge((TParseEdge)(*it)->getParseType())->name() << " " << (*it)->getPos()  << endl;
#endif
          sig = (EdgeInst *)(*it);
          sigType = (TEdgeId2)sig->getType();
          
          if(sigType == START) {
            hasStart = true;
            std::string & name = fsm.edge((TParseEdge)sig->getParseType())->name();
            if(!name.compare("PEPSTART_F") || !name.compare("PEPSTART_B"))
              hasPeptide = true;
          }
          else if(sigType == STOP)
            hasStop = true;
        }
      }

      if(t && 
	 (t->exons().size() || 
	  (reportNonCodingGenes && 
	   (t->TPexons().size() || t->FPexons().size()))))
        updateAltTranscript(listGenes, gene, t, hasStart, 
                            hasStop, hasPeptide, phase);
      
      assert(phase == tags[i].endPhase());
#ifndef NDEBUG
      cerr << "final phase = " << phase << " " << tags[i].endPhase() << endl;    
#endif
    }
  }
  
  EdgeInst *GeneUtils::_getSignal(Edge *fsmEdge, int sigPos, 
                                  TypedFilter<EdgeInst> **signals,
                                  string msg) {

    EdgeInst *sig = NULL;
    
    if(fsmEdge->id2() != NO_EDGE_INST)
      sig = &signals[fsmEdge->id2()]->value(sigPos, fsmEdge->strand());
    
    if(!sig || sig->getType() == NO_EDGE_INST) {
      sig = new EdgeInst(fsmEdge->id2(), fsmEdge->id(), 
                         sigPos, fsmEdge->strand());
      if(msg.length())
        cerr << "Added Signal " <<  fsmEdge->strand() << " " << sigPos << " " << sig->getType() << " " << fsmEdge->name() << msg << endl;
    }
    else 
      sig = new EdgeInst(*sig);
    
    sig->setParseType(fsmEdge->id());

#ifndef NDEBUG
    cerr << "signal down " << fsmEdge->name() << " " << sig->getPos() << " " << sig->getStrand() << endl;
#endif

    return sig;
  }

  /**
   * Transforms a generic exon structure into a SeqTags object
   */
  void GeneUtils::exonStructure2Tags(vector<Exon> & v,
                                     FSM *fsm, 
                                     TypedFilter<UCHAR> *contexts,
                                     TypedFilter<EdgeInst> **signals,
                                     Sequence &c,
                                     Edge **fsmEdge,
                                     Node **fsmNode,
                                     std::string &intronNode,
                                     std::string &lintronNode,
                                     const char *singleExonNode,
                                     const char *initExonNode,
                                     const char *innerExonNode,
                                     const char *lastExonNode,
                                     SeqTags &listTags,
                                     int lExonEnd,
                                     int nExonBeg,
                                     TStrand strand
                                     ) {
    
    EdgeInst *sig = NULL;
    NodeInst *b = NULL;
    Node *tmpNode = NULL;
    int sigPos = -1;
    int to = c.length(); 
    int gcClass = contexts->value(1, strand);
    Node *longIntronNode = fsm->findNode(lintronNode);

    for(unsigned int i = 0; i < v.size(); i++) {
      if(i > 0) { // inserting utr intron
        if(v[i].begin() > v[i - 1].end() + 1) {
          tmpNode = fsm->node(intronNode);
          if(longIntronNode && v[i].begin() - v[i - 1].end() + 1 > (*fsmNode)->maxLength(gcClass)) 
            tmpNode = longIntronNode;
          
          *fsmEdge = fsm->edge((*fsmNode)->id(), tmpNode->id());
          sigPos = v[i - 1].end() + 1 - (*fsmEdge)->nextNodePos();        
          sig = _getSignal(*fsmEdge, sigPos, signals,
                           string(" utrdonor+|utracceptor- for ") + c.id());
          listTags.tpush_back((Tag *)sig); 

          *fsmNode = tmpNode;                      
          b = new NodeInst((*fsmNode)->id3(), (*fsmNode)->id(), 
                           v[i - 1].end() + 1, strand, 
                           v[i].begin() - v[i - 1].end() - 1);
          
          listTags.tpush_back((Tag *)b); 

#ifndef NDEBUG
          cerr << "intron-st " << (*fsmNode)->name() << " " << b->getPos() << " " << b->getLen() << " " << b->getStrand() << endl;
#endif
        }
      }

      tmpNode = fsm->node(innerExonNode);
      if(v.size() == 1)
        tmpNode = fsm->node(singleExonNode);
      else {
        if(i == v.size() - 1)
          tmpNode = fsm->node(lastExonNode);
        else if(i == 0)
          tmpNode = fsm->node(initExonNode);
      }
      
      if(i == 0)
        if(lExonEnd > 0 && v[0].begin() - lExonEnd > 1)
          if(v.size() == 1)
            tmpNode = fsm->node(lastExonNode);
          else
            tmpNode = fsm->node(innerExonNode);

      if(i == v.size() - 1)
        if(nExonBeg > 0 && nExonBeg - v[v.size() - 1].end() > 1)
          if(v.size() == 1)
            tmpNode = fsm->node(initExonNode);
          else
            tmpNode = fsm->node(innerExonNode);         

      Edge *fsmEdge2 = fsm->edge((*fsmNode)->id(), tmpNode->id());
      if(fsmEdge2->id2() != START)
	*fsmEdge = fsmEdge2;
      sigPos = v[i].begin() - (*fsmEdge)->nextNodePos();
      
      if(v[i].end() >= v[i].begin()) { // insert utr
        *fsmNode = tmpNode;
        
        sig = _getSignal(*fsmEdge, sigPos, signals, 
                         string(" utracceptor+|utrdonor- for ") + c.id()); 
	
        listTags.tpush_back((Tag *)sig);
               
        b = new NodeInst((*fsmNode)->id3(), (*fsmNode)->id(), 
                         v[i].begin(), strand, 
                         v[i].end() - v[i].begin() + 1);

#ifndef NDEBUG
        cerr << "exon-st " << (*fsmNode)->name() << " " << b->getPos() << " " << b->getLen() << " " << b->getStrand() << endl;
#endif
        listTags.tpush_back((Tag *)b);

      }
    }
  }

  /**
   * Routine to transform an sequence c which has Gene objects in 
   * set typeSet into a vector of SeqTags objects. Each element in the vector
   * corresponds to an alternative transcript which belongs to some gene.
   * @todo This routine is very complex.. it needs to be replaced sometime later:
   * the correct approach should be to use a tree structure and a queue to "match"
   * the states parsed in the sequence with types TNodeId2 and TEdgeId2 with
   * the TParseTypes that the FSM knows.
   * As usual, input genes must be in OneStrand format.
   * @todo insert phaseBreaks... not needed now, since training data needs to be 
   * clean of phase breaks
   */
  void GeneUtils::annotSeq2Tags(vector<SeqTags> & listTags, 
				FSM *fsm, 
				TypedFilter<UCHAR> *contexts,
				TypedFilter<EdgeInst> **signals,
				Sequence &c, 
				TSetType typeSet) {
    
    SeqTags currSeqTags;
    list<Gene *> &listGenes = (list<Gene *> &)c.getAnnotBioFeats(typeSet);
    list<Gene *>::iterator it;
    SeqTags::iterator git;
    list<Transcript *>::iterator tit;

    list<Transcript *> transcripts;
    vector<int> transcriptEnds;
    vector<Transcript *> fivePts;
    vector<Transcript *> threePts;
    
    unsigned int i, j;

    // get all transcripts in one list
    for(it = listGenes.begin(); it != listGenes.end(); it++) {
      Gene &gene = *(*it);
      for(j = 0; j < gene.transcripts().size(); j++) {
        Transcript &t = gene.transcripts()[j];
        transcripts.push_back(&t);
      }
    }

    transcripts.sort(greater<Transcript *>());

    NodeInst *b = NULL;
    Node *tmpNode = NULL;
    Edge *fsmEdge = NULL, *fsmLastEdge = NULL;
    TStrand strand = NO_STRAND;
    int sigPos = -1, lastSigPos = -1, intronLen = 0;
    SeqTags *chosen = NULL;
    vector<Exon> *utrs = NULL;
    int chosenInd = -1;

    for(tit = transcripts.begin(); tit != transcripts.end(); tit++) {
      Transcript &t = *(*tit);
      strand = t.getStrand();

      std::string utr5PIntronNode = (strand == STRAND_FWD) ? 
        "5UTR_INTRON_F" : "5UTR_INTRON_B";
      std::string utr3PIntronNode = (strand == STRAND_FWD) ?
        "3UTR_INTRON_F" : "3UTR_INTRON_B";
      std::string utr5PLIntronNode = (strand == STRAND_FWD) ? 
        "5UTR_LONG_INTRON_F" : "5UTR_LONG_INTRON_B";
      std::string utr3PLIntronNode = (strand == STRAND_FWD) ? 
        "3UTR_LONG_INTRON_F" : "3UTR_LONG_INTRON_B";
      std::string igenicNode = (strand == STRAND_FWD) ? 
        "INTERGENIC_F" : "INTERGENIC_B";
      std::string lintronNode = (strand == STRAND_FWD) ? 
        "LONG_INTRON_F" : "LONG_INTRON_B"; 
      std::string intronNode = (strand == STRAND_FWD) ? 
        "INTRON_F" : "INTRON_B"; 
      std::string interExonNode = (strand == STRAND_FWD) ? 
        "INTERNAL_EXON_F" : "INTERNAL_EXON_B"; 
      std::string linterExonNode = (strand == STRAND_FWD) ?
        "LIN_INTERNAL_EXON_F" : "LIN_INTERNAL_EXON_B"; 
      std::string sinterExonNode = (strand == STRAND_FWD) ? 
        "SIN_INTERNAL_EXON_F" : "SIN_INTERNAL_EXON_B"; 
      std::string initExonNode = (strand == STRAND_FWD) ? 
        "INIT_EXON_F" : "INIT_EXON_B"; 
      std::string singleExonNode = (strand == STRAND_FWD) ? 
        "SINGLE_EXON_F" : "SINGLE_EXON_B"; 
      std::string lastExonNode = (strand == STRAND_FWD) ? 
        "LAST_EXON_F" : "LAST_EXON_B"; 

      Node *fsmNode = fsm->node(SYNC_BEG_STATE), *fsmLastNode = NULL;
      EdgeInst *sigUp = NULL, *sigDw = NULL;
      Node *nextIntron = NULL;

      vector<Exon> & ve5 = t.FPexons();
      vector<Exon> & ve = t.exons();
      vector<Exon> & ve3 = t.TPexons();

      int begInterg = (chosenInd >= 0) ? transcriptEnds[chosenInd] + 1 : 1;
      int endInterg = t.begin() - 1;
      int gcClass = contexts->value(begInterg, strand);
      Node *longIntronNode = fsm->findNode(lintronNode);
      Node *fivePlintronNode = fsm->findNode(utr5PLIntronNode);
      Node *threePlintronNode = fsm->findNode(utr3PLIntronNode);

      if(begInterg > endInterg + 1) {
        begInterg = 1;
        chosenInd = -1;
        chosen = NULL;
        for(unsigned int i = 0; i < listTags.size(); i++) {

          if(transcriptEnds[i] < endInterg) {
            begInterg = transcriptEnds[i] + 1;
            assert(c.length() > transcriptEnds[i]);
            listTags[i].pop_back();
            listTags[i].pop_back();
            chosenInd = i;
            chosen = &listTags[i];
            break;
          }

        }
      }
      else {
        if(chosenInd >= 0) { // poping the last intergenic region
          listTags[chosenInd].pop_back();
          listTags[chosenInd].pop_back();
        }
      }
      assert(endInterg + 1 >= begInterg);
      // dealing with upstream intergenic regions 
      if(endInterg >= begInterg) {
        tmpNode = fsm->node(igenicNode);
        fsmLastNode = fsmNode;
        fsmNode = tmpNode;

        if(begInterg == 1) { // first gene

          if((!t.hasStart() && strand == STRAND_FWD) || 
             (!t.hasStop() && strand == STRAND_COMP))  { // check if it's actually an intron
            fsmLastEdge = fsm->edge((strand == STRAND_FWD) ? "ACCEPTOR_F" : "DONOR_B");
            lastSigPos = endInterg + 1 - fsmLastEdge->nextNodePos();
            sigDw = &signals[fsmLastEdge->id2()]->value(lastSigPos, strand);

            if(sigDw->getType() == fsmLastEdge->id2()) {
              tmpNode = fsm->node(intronNode);

              if(longIntronNode != NULL &&
                 endInterg - begInterg > longIntronNode->minLength(gcClass))
                tmpNode = longIntronNode;

              fsmNode = tmpNode;
              
              /*
               * Checking next intron's type
               */
              nextIntron = NULL;

              if(ve.size() > 1) {
                if(longIntronNode != NULL &&
                   ve[1].begin() - ve[0].end() > longIntronNode->minLength(gcClass))
                  nextIntron = longIntronNode;
                else
                  nextIntron = tmpNode;
              }
            }
            else
              sigDw = NULL; // not an splice site, but no start either
          }
          
          fsmEdge = fsm->edge(fsmLastNode->id(), fsmNode->id());
          sigPos = begInterg + fsmLastNode->nextEdgePos(fsmNode->id());
          sigUp = _getSignal(fsmEdge, sigPos, signals, string(""));
          currSeqTags.tpush_back((Tag *)sigUp);
        }
        
        b = new NodeInst(fsmNode->id3(), fsmNode->id(), 
                         begInterg, strand, 
                         endInterg - begInterg + 1);
        
        currSeqTags.tpush_back((Tag *)b);
#ifndef NDEBUG
        cerr << "ig-state " << fsmNode->name() << " " << b->getPos() << " " << b->getLen() << " " <<  b->getStrand() << endl;
#endif
      }
      
      // dealing with utr exons
      utrs = &ve5;
      if(strand == STRAND_COMP)
        utrs = &ve3;

      if(utrs->size() > 0) {
        exonStructure2Tags(*utrs, fsm, contexts,
                           signals, c,
                           &fsmEdge,
                           &fsmNode,
                           (strand == STRAND_FWD) ?
                           utr5PIntronNode : utr3PIntronNode,
                           (strand == STRAND_FWD) ? 
                           utr5PLIntronNode : utr3PLIntronNode,
                           (strand == STRAND_FWD) ?
                           "5SINGLE_UTR_F" : "3SINGLE_UTR_B",
                           (strand == STRAND_FWD) ? 
                           "5INIT_UTR_F" : "3LAST_UTR_B",
                           (strand == STRAND_FWD) ? 
                           "5INTERNAL_UTR_F" : "3INTERNAL_UTR_B",
                           (strand == STRAND_FWD) ? 
                           "5LAST_UTR_F" : "3INIT_UTR_B",
                           currSeqTags,
                           -1,
                           (ve.size() > 0 ? ve[0].begin() : -1),
                           strand
                           );
	
        fsmLastEdge = fsmEdge;
        fsmLastNode = fsmNode;

        if(strand == STRAND_FWD && ve.size() > 0) {
          intronLen = ve[0].begin() - (*utrs)[utrs->size() - 1].end() - 1;
          if(intronLen > 0) {
            // there is a UTR intron
            
            tmpNode = fsm->node(utr5PIntronNode);
            if(fivePlintronNode && intronLen > tmpNode->maxLength(gcClass))
              tmpNode = fivePlintronNode;
            
            fsmEdge = fsm->edge(fsmNode->id(), tmpNode->id());
            sigPos = (*utrs)[utrs->size() - 1].end() + 1 - fsmEdge->nextNodePos();
              
            
            fsmNode = tmpNode;
            sigUp = _getSignal(fsmEdge, sigPos, signals, 
                               string(" utrdonor+ for ") + c.id());
            currSeqTags.tpush_back((Tag *)sigUp);

            b = new NodeInst(fsmNode->id3(), fsmNode->id(), 
                             (*utrs)[utrs->size() - 1].end() + 1, strand, 
                             intronLen);
            
            currSeqTags.tpush_back((Tag *)b);          
#ifndef NDEBUG
            cerr << "intron-utr " << fsmNode->name() << " " << b->getPos() << " " << b->getLen() << " " <<  b->getStrand() << endl;
#endif

          }
        }       
      }

      utrs = &ve3;
      if(strand == STRAND_COMP)
        utrs = &ve5;

      int coding_len = 0;
      for(i = 0; i < ve.size(); i++) 
	coding_len += ve[i].end() - ve[i].begin() + 1;
    
      for(i = 0; i < ve.size(); i++) {
        fsmLastNode = fsmNode;

        if(i > 0) { // insert intron
          intronLen = ve[i].begin() - ve[i - 1].end() - 1;
          if(intronLen > 0) {
            tmpNode = fsm->node(intronNode);
            
            if(longIntronNode != NULL &&
               intronLen >= longIntronNode->minLength(gcClass))
              tmpNode = longIntronNode;
            
            fsmNode = tmpNode;
            
            // next intron's type
            nextIntron = NULL;

            if(i < ve.size() - 1) {
              if(longIntronNode != NULL && 
                 ve[i+1].begin() - ve[i].end() > longIntronNode->minLength(gcClass))
                nextIntron = longIntronNode;
              else
                nextIntron = tmpNode;
            }
            
            b = new NodeInst(fsmNode->id3(), fsmNode->id(), 
                             ve[i - 1].end() + 1, strand, 
                             ve[i].begin() - ve[i - 1].end() - 1);
            
            // inserting intron
            currSeqTags.tpush_back((Tag *)b); 
            
#ifndef NDEBUG
            cerr << "intron-st " << fsmNode->name() << " " << b->getPos() << " " << b->getLen() << " " << b->getStrand() << endl;
#endif
          }
          fsmLastNode = fsmNode;
        }
        
        tmpNode = fsm->node(interExonNode);

        // normal flow, no truncated genes
        if(fsmNode->id3() == INTERGENIC || fsmNode->id3() == UTR) {
          // 5' edge
          fsmLastEdge = fsm->edge((strand == STRAND_COMP) ? 
                                  "STOP_B" : 
                                  (t.hasPeptide() ? "PEPSTART_F" : "START_F"));
          lastSigPos = ve[i].begin() - fsmLastEdge->nextNodePos();
          sigUp = &signals[fsmLastEdge->id2()]->value(lastSigPos, strand);

          if(sigUp->getType() != fsmLastEdge->id2())
            sigUp = NULL;

          if(ve.size() == 1 &&
             ((t.hasStop() && strand == STRAND_FWD) ||
              (t.hasStart() && strand == STRAND_COMP)))
            
              tmpNode = fsm->node(singleExonNode);
          else {
            if(strand == STRAND_COMP)
              tmpNode = fsm->node(lastExonNode);
            else
              tmpNode = fsm->node(initExonNode);
          }
        }
        //5' intron-truncated gene
        else if(fsmNode->id3() == INTRON) {
          fsmLastEdge = fsm->edge((strand == STRAND_FWD) ? 
                                  "ACCEPTOR_F" : "DONOR_B");
          lastSigPos = ve[i].begin() - fsmLastEdge->nextNodePos();
          sigUp = &signals[fsmLastEdge->id2()]->value(lastSigPos, strand);

          if(sigUp->getType() != fsmLastEdge->id2())
            sigUp = NULL;

          if(i == ve.size() - 1) { 
            if(strand == STRAND_COMP && t.hasStart())
              tmpNode = fsm->node(initExonNode);
            else if(strand == STRAND_FWD && t.hasStop())
              tmpNode = fsm->node(lastExonNode);
          }
        }

        // 5' exon-truncated gene
        else if(fsmNode->id() == SYNC_BEG_STATE) { 
          assert(ve[i].begin() == 1);
          if(ve.size() == 1 &&
             ((t.hasStop() && strand == STRAND_FWD) ||
              (t.hasStart() && strand == STRAND_COMP)))
            
            tmpNode = fsm->node(singleExonNode);
          else { // multiple-exon gene
            if(strand == STRAND_FWD)
              tmpNode = fsm->node(initExonNode);
            else
              tmpNode = fsm->node(lastExonNode);
          }
          
          fsmLastEdge = fsm->edge(SYNC_BEG_STATE, tmpNode->id());
          lastSigPos = 0;
          sigUp = NULL;      
        }

        // looking for the downstream signal, a splice site or translation
        if(i == ve.size() - 1 &&
           ((t.hasStop() && strand == STRAND_FWD) ||
            (t.hasStart() && strand == STRAND_COMP))) {
	  
          fsmEdge = fsm->edge((strand == STRAND_COMP) ?
                              (t.hasPeptide() ? "PEPSTART_B" : "START_B") :
                              "STOP_F");
	  //	  cerr << c.id() << " " << ve[i].end() + 1 - fsmEdge->nextNodePos() << " " << fsmEdge->name() << " " << strand << " " <<  t.hasPeptide() << endl;
	}
        else
          fsmEdge = fsm->edge((strand == STRAND_COMP) ?
                              "ACCEPTOR_B" : "DONOR_F");
	
        sigPos = ve[i].end() + 1 - fsmEdge->nextNodePos();
        sigDw = &signals[fsmEdge->id2()]->value(sigPos, strand);
	
        if(sigDw->getType() != fsmEdge->id2()) {
          sigDw = NULL;
          if(i == ve.size() - 1) {
            if(ve.size() == 1)  // single-exon gene
              tmpNode = fsm->node(singleExonNode);
            else { // multiple-exon gene
              if(strand == STRAND_FWD)
                tmpNode = fsm->node(lastExonNode);
              else
                tmpNode = fsm->node(initExonNode);
            }
            
            if(ve[i].end() == c.length()) { // maybe 3' exon-truncated gene
              fsmEdge = fsm->edge(tmpNode->id(), SYNC_END_STATE);
              sigPos = ve[i].end() + 1;
            }
          }
        }

        if(tmpNode->id2() == INTERNAL_EXON && nextIntron != NULL) {
          Node *tmpNode2 = tmpNode;

          if(fsmLastNode->id() == nextIntron->id()) {
            if(fsmLastNode->name().compare("LONG_INTRON_F") == 0 || fsmLastNode->name().compare("LONG_INTRON_B") == 0) {
              tmpNode = fsm->findNode(linterExonNode);
              if(!tmpNode)
                tmpNode = tmpNode2;
            }
            else {
              tmpNode = fsm->findNode(sinterExonNode);
              if(!tmpNode)
                tmpNode = tmpNode2;
            }          

            Edge *tmpEdge = fsm->findEdge(fsmLastNode->id(), tmpNode->id());
            if(tmpEdge)
              fsmLastEdge = tmpEdge;
            else
              tmpNode = tmpNode2;

            tmpEdge = fsm->findEdge(tmpNode->id(), nextIntron->id());
            if(tmpEdge)
              fsmEdge = tmpEdge;
            else
              tmpNode = tmpNode2;

          }
        }

        if(ve[i].end() >= ve[i].begin()) { // insert exon
          fsmLastNode = fsmNode;

          fsmNode = tmpNode;
          if(!sigUp) {
            sigUp = new EdgeInst(fsmLastEdge->id2(), fsmLastEdge->id(), 
                                 lastSigPos, strand);

            cerr << "Added Signal " << strand << " " << i << " " << ve[i].begin() << " " << sigUp->getType() << " " << fsmLastEdge->name() << " start+|acceptor+|donor-|stop- for " << c.id() << "\n";
          }
          else
            sigUp = new EdgeInst(*sigUp);
          
          sigUp->setParseType(fsmLastEdge->id());

#ifndef NDEBUG
          cerr << "signal up " << fsmLastEdge->name() << " " << sigUp->getPos() << " " << sigUp->getStrand() << endl;
#endif

          b = new NodeInst(fsmNode->id3(), fsmNode->id(), 
                           ve[i].begin(), strand, 
                           ve[i].end() - ve[i].begin() + 1);

          if(!sigDw) {
            sigDw = new EdgeInst(fsmEdge->id2(), fsmEdge->id(), sigPos, strand);
            cerr << "Added Signal " << strand << " " << i << " " << ve[i].end() << " " << sigDw->getType() << " " << fsmEdge->name() << " stop+|donor+|start-|acceptor- for " << c.id() << "\n";
          }
          else
            sigDw = new EdgeInst(*sigDw);

          sigDw->setParseType(fsmEdge->id());
          currSeqTags.tpush_back((Tag *)sigUp); // the upstream signal
          currSeqTags.tpush_back((Tag *)b);   // the state
	  
	  if(i == 0) {
	    if(strand == STRAND_FWD && t.phase())
	      currSeqTags.insertPhaseBreak(b, t.phase());
	    if(strand == STRAND_COMP) {
	      int comp_phase = (t.phase() + coding_len) % 3;
	      if(comp_phase) 
		currSeqTags.insertPhaseBreak(b, (3 - comp_phase) % 3);
	    }
	  }

#ifndef NDEBUG
          cerr << "exon-state " << fsmNode->name() << " " << b->getPos() << " " << b->getLen() << " " << b->getStrand() << endl;
          cerr << "signal dw " << fsmEdge->name() << " " << sigDw->getPos() << " " << sigDw->getStrand() << endl;
#endif

          if(!utrs->size() || (utrs->size() > 0 && i < ve.size() - 1))
            currSeqTags.tpush_back((Tag *)sigDw); // the downstream signal
	  else 
	    delete sigDw;
          
        }
        else
          assert(0);
      }

      // dealing with utr exons
      if(utrs->size() > 0 && ve.size() > 0 && strand == STRAND_COMP) {
        intronLen = (*utrs)[0].begin() - ve[ve.size() - 1].end() - 1;
        if(intronLen > 0) {
	  throw EXCEPTION(BAD_USAGE, "can't support this");
          tmpNode = fsm->node(utr5PIntronNode);
          if(fivePlintronNode && intronLen > tmpNode->maxLength(gcClass)) 
            tmpNode = fivePlintronNode;
          
          fsmEdge = fsm->edge(fsmNode->id(), tmpNode->id());
          sigPos = ve[ve.size() - 1].end() + 1 - fsmEdge->nextNodePos();
          
          fsmNode = tmpNode;
	  sigUp = _getSignal(fsmEdge, sigPos, signals, 
			     string(" utracceptor- for ") + c.id());
	  currSeqTags.tpush_back((Tag *)sigUp);

          b = new NodeInst(fsmNode->id3(), fsmNode->id(), 
                           ve[ve.size() - 1].end() + 1, strand,
                           intronLen);
          
          currSeqTags.tpush_back((Tag *)b);
#ifndef NDEBUG
          cerr << "intron-utr " << fsmNode->name() << " " << b->getPos() << " " << b->getLen() << " " <<  b->getStrand() << endl;
#endif

        }
      }
      if(utrs->size()) {
        exonStructure2Tags(*utrs, fsm, contexts,
                           signals, c,
                           &fsmEdge,
                           &fsmNode,
                           (strand == STRAND_FWD) ? 
                           utr3PIntronNode : utr5PIntronNode,
                           (strand == STRAND_FWD) ? 
                           utr3PLIntronNode : utr5PLIntronNode,
                           (strand == STRAND_FWD) ?
                           "3SINGLE_UTR_F" : "5SINGLE_UTR_B",
                           (strand == STRAND_FWD) ? 
                           "3INIT_UTR_F" : "5LAST_UTR_B",
                           (strand == STRAND_FWD) ? 
                           "3INTERNAL_UTR_F" : "5INTERNAL_UTR_B",
                           (strand == STRAND_FWD) ? 
                           "3LAST_UTR_F" : "5INIT_UTR_B",
                           currSeqTags,
                           (ve.size() > 1 ? ve[ve.size() - 1].end() : -1),
                           -1,
                           strand
                           );
        
        fsmLastEdge = fsmEdge;
        fsmLastNode = fsmNode;

        fsmEdge = fsm->edge(fsmNode->id(),
                            fsm->node((strand == STRAND_FWD) ? 
                                      "INTERGENIC_F" : "INTERGENIC_B")->id());
        if((*utrs)[utrs->size() - 1].end() == c.length())  // maybe 3' exon-truncated gene
          fsmEdge = fsm->edge(fsmNode->id(), SYNC_END_STATE);
        
        sigPos = (*utrs)[utrs->size() - 1].end() + 1 - fsmEdge->nextNodePos();
        sigDw = _getSignal(fsmEdge, sigPos, signals,
                           string(" tes+|tss- for ") + c.id());

        currSeqTags.tpush_back((Tag *)sigDw);
        
      }
      // dealing with downstream intergenic regions
      if(c.length() > t.end()) {
        if(!fsmEdge->name().compare("DONOR_F") ||
           !fsmEdge->name().compare("ACCEPTOR_B")) {
          
          tmpNode = fsm->node(intronNode);
          if(longIntronNode && c.length() - t.end() + 1 > tmpNode->maxLength(gcClass))
            tmpNode = longIntronNode;
        }
        else
          tmpNode = fsm->node(igenicNode);          

        fsmLastNode = fsmNode;
        fsmNode = tmpNode;

        b = new NodeInst(fsmNode->id3(), fsmNode->id(), 
                         t.end() + 1, strand, 
                         c.length() - t.end());

        currSeqTags.tpush_back((Tag *)b);  
#ifndef NDEBUG
        cerr << "ig-state " << fsmNode->name() << " " << b->getPos() << " " << b->getLen() << " " << b->getStrand() << endl;
#endif
        fsmLastNode = fsmNode;
        fsmEdge = fsm->edge(fsmLastNode->id(), SYNC_END_STATE);
        sigPos = c.length() + 1 + fsmLastNode->nextEdgePos(SYNC_END_STATE);
        sigDw = new EdgeInst(fsmEdge->id2(), fsmEdge->id(), sigPos, strand);
        currSeqTags.tpush_back((Tag *)sigDw);
#ifndef NDEBUG
        cerr << "signal " << fsmEdge->name() << " " << sigDw->getPos() << " " << sigDw->getStrand() << endl;
#endif
      }

      /*
       * Check whether to merge current list in an exisiting Label Set
       * or create a new Label Set.
       */
      if(chosen) {
        assert(chosenInd >= 0);
        transcriptEnds[chosenInd] = t.end();
        threePts[chosenInd] = *tit;
      }
      else {
        SeqTags tmpSeqTags; 
        listTags.push_back(tmpSeqTags);
        chosen = &(listTags.back());
        chosen->setRank(listTags.size() - 1);
        transcriptEnds.push_back(t.end());
        fivePts.push_back(*tit);
        threePts.push_back(*tit);
        chosenInd = transcriptEnds.size() - 1;
      }

      for(git = currSeqTags.begin(); git != currSeqTags.end(); ) {
        chosen->push_back(*git);
        if((*git)->getGEClass() == NODE_INST) {
          int phaseBreak = currSeqTags.phaseBreak(*git);
          if(phaseBreak > 0) // phaseBreak only could happen in between genes
            chosen->insertPhaseBreak(chosen->back(), phaseBreak);
	}
        git = currSeqTags.erase(git);
      }
    }
    
    // setting phases
    assert(fivePts.size() == threePts.size());
    int currLen;

    for(i = 0; i < transcriptEnds.size(); i++) {
      if(fivePts[i]) {
	int coding_len = 0;
	vector<Exon> &fve = fivePts[i]->exons();
	for(int j = 0; j < fve.size(); j++) 
	  coding_len += fve[j].end() - fve[j].begin() + 1;
	
	if(fivePts[i]->getStrand() == STRAND_FWD) {
	  int fwd_phase = fivePts[i]->phase();
	  listTags[i].setInitPhase(fivePts[i]->phase());
	//	  listTags[i].setEndPhase((currLen + fivePts[i]->phase()) % 3);
	}
	else {
	  int comp_phase = (fivePts[i]->phase() + coding_len) % 3;
	  listTags[i].setInitPhase((3 - comp_phase) % 3);
	}
      }
      if(threePts[i]) {
	int coding_len = 0;
	vector<Exon> &tve = threePts[i]->exons();
	for(int j = 0; j < tve.size(); j++) 
	  coding_len += tve[j].end() - tve[j].begin() + 1;

	if(threePts[i]->getStrand() == STRAND_FWD) {
	  int fwd_phase = (threePts[i]->phase() + coding_len) % 3;
	  listTags[i].setEndPhase(fwd_phase);
	//	  listTags[i].setEndPhase((currLen + fivePts[i]->phase()) % 3);
	}
	else {
	  int comp_phase = threePts[i]->phase();
	  listTags[i].setEndPhase((3 - comp_phase) % 3);
	}

	  //        listTags[i].setEndPhase((3 - threePts[i]->phase()) % 3);
	  //        int threeLen = currLen + threePts[i]->phase();
	  //        listTags[i].setInitPhase((3 - threeLen % 3) % 3);
	
      }      
    }
    
    /* 
     * If annotSeq does not contain any transcripts, then insert an intergenic
     * region
     */
    
    if(!listTags.size()) {
      Node *fsmNode = fsm->node("INTERGENIC_F");
      fsmLastEdge = fsm->edge(SYNC_BEG_STATE, fsmNode->id());
      fsmEdge = fsm->edge(fsmNode->id(), SYNC_END_STATE);

      EdgeInst *sigUp = new EdgeInst(fsmEdge->id2(), fsmLastEdge->id(), 0, STRAND_FWD);
      b = new NodeInst(fsmNode->id3(), fsmNode->id(), 
                       1, STRAND_FWD, 
                       c.length());
      EdgeInst *sigDw = new EdgeInst(fsmEdge->id2(), fsmEdge->id(), c.length() + 1, STRAND_COMP);

      listTags.push_back(currSeqTags);
      chosen = &(listTags.back());
      chosen->setRank(listTags.size() - 1);
      chosen->tpush_back((Tag *)sigUp); // the upstream signal
      chosen->tpush_back((Tag *)b);   // the state
      chosen->tpush_back((Tag *)sigDw); // the downstream signal
      
      
    }
  }

  /**
   * Function that parses sequence header which specifies the exon locations.
   * This only applies to sequences in locs format.
   */
  void GeneUtils::parseSequenceLocation(std::string &location, 
                                        std::string &annotSeqId, 
                                        vector<Pair<int> > & ve5, 
                                        vector<Pair<int> > & ve, 
                                        vector<Pair<int> > & ve3, 
                                        bool & hasStart, 
                                        bool & hasStop,
                                        int & phase) {

    hasStart = true;
    hasStop = true;
    phase = 0;
    annotSeqId = "";
    /*
     * Process 5utr regions
     */    
    if(rExLocation1.Match(location)) {
      location = rExLocation1[2];
      readExonLocations(rExLocation1[1], annotSeqId, ve5);
    }

    /*
     * Process 3utr regions
     */    
    if(rExLocation2.Match(location)) {
      location = rExLocation2[1];
      readExonLocations(rExLocation2[2], annotSeqId, ve3);
    }

    /*
     * Check for truncated genes
     */
    if(rExLocation3.Match(location)) {
        if(rExLocation3[1].length()) 
          sscanf(rExLocation3[1].c_str(),"%d", &phase);

        location = rExLocation3[2];
        hasStart = false;
    }

    /*
     * Process coding regions
     */
    if(rExLocation4.Match(location)) {
      location = rExLocation4[1];
      hasStop = false;
    }
       
    readExonLocations(location, annotSeqId, ve);
  }

  /**
   * Function called from parseSequenceLocation in order to extract the exon
   * boundaries which are part of the sequence header.
   * @see parseSequenceLocation
   */
  void GeneUtils::readExonLocations(std::string location, std::string & annotSeqId, vector<Pair<int> > & v) {

    int beg, end;
    std::vector<std::string> locs;
    rExSep.Split(locs, location, 0, 100);

    for(unsigned int i = 0; i < locs.size(); i++) {
      bool matched = rExLoc.Match(locs[i]);

      if(!matched) {
        assert(0);
        throw EXCEPTION(BAD_LOCATION_FORMAT, location);
      }
      
      if(annotSeqId.length()) {
        if(annotSeqId.compare(rExLoc[1])) {
          assert(0);
          throw EXCEPTION(BAD_LOCATION_FORMAT, "Location spawns multiple contig sequences");
        }
      }
      else annotSeqId = rExLoc[1];

      sscanf(rExLoc[2].c_str(),"%d", &beg);
      sscanf(rExLoc[3].c_str(),"%d", &end);
      
      Pair<int> pi(beg, end);
      v.push_back(pi);
      
    }
  }


  void GeneUtils::computeRepeatedElems(list<RepeatedElem> & reps, 
                                      Sequence &c,
                                      int minLen, TStrand strand) {

    int i, cLen = c.length();
    char *seq = c.getSeq(strand);

    if(!minLen || cLen < minLen)
      minLen = cLen;

    int currLen = seq[0] == 'N' ? 1 : 0;

    for( i = 1; i < cLen; i++) {
      if(seq[i] == 'N')
        currLen++;

      if(seq[i-1] != 'N' && (i == cLen - 1 || seq[i + 1] != 'N')) {
        if(currLen >= minLen) {
          ostringstream reId;
          reId << c.id() << i;
          std::string reId2 = reId.str();

          if(strand == STRAND_FWD)
            reps.push_back(RepeatedElem(reId2, i - currLen + 2, i + 1, &c, true, strand));
          else
            reps.push_back(RepeatedElem(reId2, cLen - i + currLen - 1, cLen - i, &c, true, strand));

        }
        currLen = 0;
      }
    }
  }

}
  
  /*
  void FilterEngine::computeCodingDifferential(FeatureEngine *params, 
                                                   TSetType geneSet, 
                                                   TStrand strand) {

    if(!getSequence()->hasAnnotBioFeats(geneSet))
      return;

    list<Gene *> & listGenes = (list<Gene *> &)getSequence()->getAnnotBioFeats(geneSet);
    listGenes.sort();
    double sc[3] = {0,0,0};
    vector<Exon> *exons;

    for(list<Gene *>::iterator it = listGenes.begin(); it != listGenes.end(); it++) {
      Gene & expGene = *((*it));
      TStrand thisStrand = expGene.getStrand();

      if(thisStrand != strand)
        continue;

      for(unsigned int t = 0; t < expGene.transcripts().size(); t++) {
        int currTranscriptLen = 0;
        exons = &expGene.transcript(t).exons();

        for(unsigned int j = 0; j < exons->size(); j++) {
          double scPhase = 0, scNonPhase = 0;
          TParseNode pType = INTERNAL_EXON_F;

          if(j == 0) 
            pType = INIT_EXON_F;

          if(j == exons->size() - 1)
            pType = LAST_EXON_F;

          if(exons->size() == 1)
            pType = SINGLE_EXON_F;

          pType = BioStMap::switchParsNodeStrand(pType, strand);
          NodeInst b(EXON, pType, (*exons)[j].begin(), strand, (*exons)[j].end() - (*exons)[j].begin() + 1);

          for(int phase = 0; phase < 3; phase++) {
            sc[phase] += params->stFeatSum(&b, phase);

            if(phase != currTranscriptLen % 3)
              scNonPhase += sc[phase];

          }
          scPhase += sc[currTranscriptLen % 3];
          currTranscriptLen += abs((*exons)[j].end() - (*exons)[j].begin()) + 1;

          if(scPhase != 0 && scNonPhase != 0)
            cerr << getSequence()->id() << "\t" << gcClass->value(2, STRAND_FWD) << "\t" << scPhase*2.0 << "\t" << scNonPhase << "\n";
        }
      }
    }
  }
  */

  /**
   * Computes the difference in score between true positive and true negative 
   * signals. Genes must be in OneStrand format. The results are formatted so
   * that a histogram can be easily made from them.
   */
  /*
  void GeneUtils::computeSignalDifferential(FSM *fsm, FilterEngine *fe, 
                                           TypedFilter<int> *contexts,
                                           TypedFilter<EdgeInst> **signals,
                                           FeatureEngine *params, 
                                           TSetType geneSet, 
                                           TStrand strand) {

    std::map<int, int> geneSignals;
    Sequence *c = fe->getSequence();

    if(!c->hasAnnotBioFeats(geneSet))
      return;

    SeqTags expTags;
    int i = 0;
    EdgeInst *sig;
    int phase = 0;
    vector<Exon> *exons;
    list<Gene *> & listGenes = (list<Gene *> &)c->getAnnotBioFeats(geneSet);

    list<Gene *>::iterator git = listGenes.begin();
    for( ; git != listGenes.end(); git++) {
      Gene & expGene = *((*git));
      TStrand thisStrand = expGene.getStrand();

      if(thisStrand != strand)
        continue;

      expGene.toTwoStrand();

      for(unsigned int t = 0; t < expGene.transcripts().size(); t++) {
        int currTranscriptLen = 0;
        exons = &expGene.transcript(t).exons();

        for(unsigned int j = 0; j < exons->size(); j++) {
          int pos3P = -2, pos5P = 1;

          if(j == 0)  
            pos3P = 0;

          if(j == exons->size() - 1)
            pos5P = 1;

          geneSignals[(*exons)[j].begin() + pos3P] = currTranscriptLen % 3;
          currTranscriptLen += abs((*exons)[j].end() - (*exons)[j].begin()) + 1;
          geneSignals[(*exons)[j].end() + pos5P] = currTranscriptLen % 3;
        }
      }

      expGene.toOneStrand();

    }

    std::map<int, int>::iterator it;
    bool truePos;
    int sigPos;

    for(i = 1; i < c->length(); i++) {
      sigPos = i;

      if(strand == STRAND_COMP) 
        sigPos = c->length() - i + 1;

      sig = &signals->value(sigPos, strand);

      if(sig->getType() == NO_EDGE_INST) 
        continue;
      
      it = geneSignals.find(sigPos);
      phase = 0;
      truePos = false;

      if(it != geneSignals.end()) {
        phase = it->second;
        truePos = true;
      }

      Edge *e = fsm->edge("START_F");

      switch(sig->getType()) {
      case STOP: e = fsm->edge("STOP_F"); break;
      case DONOR: e = fsm->edge("DONOR_F"); break;
      case ACCEPTOR: e = fsm->edge("ACCEPTOR_F"); break;
      default: break;
      }

      if(strand != STRAND_FWD)
        e = e->complementEdge();
      sig->setParseType(e->id());

      cout << truePos << "\t" << c->id() << "\t" << sig->getType() << "\t" << sigPos << "\t" << contexts->value(2, STRAND_FWD) << "\t" << params->dotParamV4Edges(sig, phase) << endl;
    }
  }
  */


  /* 
   * In development.
   * GFF goes like this : <seqname> <source> <feature> <start> <end> <score>
   * <strand> <frame> [attributes] [comments], [attributes] usually have gene 
   * "gene_id" as first attribute
   */
  /*
  void GeneUtils::loadGTF(list<Gene> & listGenes, 
                           const char *geneFile, 
                           bool buildstartInfo, 
                           bool readfastaSequence, 
                           list<Sequence> &annotSeqs, 
                           Sigma *sigma, 
                           TSetType typeSet) {

    std::ifstream geneStream(geneFile);
    std::ostringstream ostr; 

    int $counter = 0;
    if(!geneStream) {
      assert(0);
      throw EXCEPTION( BIOFEATS_UNAVAILABLE, string(geneFile));
    }

    boost::RegEx rExDelim("\\s+");
    boost::RegEx rExComment("^\\s*\\#.*");
    boost::RegEx r1("gene",  boost::regex::perl|boost::regex::icase);
    boost::RegEx r2("gene_id",  boost::regex::perl|boost::regex::icase);
    boost::RegEx r3("\"?([^\"]+)\"?\;?",  boost::regex::perl|boost::regex::icase);
    boost::RegEx r4("transcript_id",  boost::regex::perl|boost::regex::icase);
    boost::RegEx r5("trans.+\=([^\";]+)", boost::regex::perl|boost::regex::icase);
    boost::RegEx r6("^(\S+)(\.gb\.fa)?", boost::regex::perl|boost::regex::icase);

    std::string line;
    ::vector<std::string> gtfArgs;
    bool matched;

    std::map<std::string, vector<Exon> >::iterator it;
    map<std::string, std::vector<Exon> > exons;
    std::string geneId, transcriptId;
    vector<Exon> & theseExons;  

    while(std::getline(geneStream, line), !geneStream.eof()) {
      if(rExComment.Match(line)) // a comment
        continue; 

      rExDelim.Split(gtfArgs, line, 0, 100);
      if(r1.Match(gtfArgs[2])) {
        counter++;
      }
      if(gtfArgs[8].length() > 0) { 
        if(r2.Match(gtfArgs[8])) { // twinscan
          matched = r3.Match(gtfArgs[9]);
          assert(matched);
          geneId = r3[1];
          transcriptId = r3[1];
        }
        else if(r4.Match(gtfArgs[8])) { // augustus and encode
          transcriptId = gtfArgs[9];
          if(gtfArgs[10].length() > 0 && r2.Match(gtfArgs[10])) 
            geneId = gtfArgs[11];
          else
            geneId = transcriptId;

          matched = r3.Match(transcriptId);
          assert(matched);
          transcriptId = r3[1];
          matched = r3.Match(geneId);
          assert(matched);
          ostr << counter;
          geneId = r3[1] + ostr.str();
        }
        else if(r5.Match(gtfArgs[8]))  {  // genezilla
          geneId = r5[1];
          transcriptId = r5[1];
        }
        else {
          geneId = gtfArgs[8];
          transcriptId = gtfArgs[8];
        }
      }
      else {
        geneId = gtfArgs[10];
        transcriptId = gtfArgs[10];
      }
      matched = r6.Match(gtfArgs[0]);
      assert(matched);
      annotSeqId = r6[1];
      if(oldGeneId.compare("") != 0 && oldGeneId.compare(geneId) == 0) {
        // insert alternative transcripts
        std::string tId;
        if(onlyOneTranscript) {
          int maxExons;
          std::string chosentId;

          for(it = exons.begin(); it != exons.end(); it++) {
            vector<Exon> & theseExons = it->second;
            if(theseExons.size() > maxExons) {
              chosentId =it->first;
              maxExons = it->second.size();
            }
          }
          for(it = exons.begin(); it != exons.end(); it++) {
            if(chosentId.length() > 0 && tId.compare(chosentId) != 0)
              it->second.clear();
          }
        }
        for(it = exons.begin(); it != exons.end(); it++) {
          std::string transId = tId;

          if(oldGeneId.compare("transId") == 0) 
            transId = oldGeneId + std::string("|") + tId;


          if(
        }
      }
    }
  }

  */
