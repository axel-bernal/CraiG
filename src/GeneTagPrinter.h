/****************************************************************************
 * GeneTagPrinter.h - part of the craig namespace, a genomics library
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
 ****************************************************************************/


#ifndef _GENE_TAGPRINTER_H_
#define _GENE_TAGPRINTER_H_
#include "TagPrinter.h"
#include "Sequence.h"
#include "Gene.h"
#include "GeneUtils.h"
#include "FSM.h"
#include "TagUtils.h"
#include "ArgParseUtils.h"
#include "IntervalSet.h"


namespace craig {

  /**
   * This class overrides the methods of class TagPrinter in order to display
   * genes as lists of exons and using known display/output formats, such
   * as GTF.
   ***************************************************************************/

  class GeneTagPrinter : public TagPrinter {
    bool reportNonCoding;
   public:
    GeneTagPrinter() {
      reportNonCoding = false;
    }
    
    inline void reportNonCodingGenes() {
      reportNonCoding = true;
    }

    void displayHeader(std::string &fmt, ::ofstream &fd);
    void displayTranscript(std::string &fmt,
			   std::ofstream &fd,
			   std::string &gid,
			   Transcript &t,
			   int offset = 0,
			   int upperLimit = INT_MAX,
			   IntervalSet<double> *tag_scores = NULL);

    void prepareExons4Display(Transcript &t,
			      int offset,
			      int upperLimit,
			      vector<Exon> & exons,
			      vector<pair<int, int> > &v,
			      vector<double> &v_scores,
			      IntervalSet<double> *tag_scores = NULL);
      
      
    void displayCoordsAsGTF(std::ofstream &fd,
			    Transcript &t, 
			    int initPhase,
			    std::string & gtf_detail,
			    vector<pair<int, int> > & v,
			    vector<double> &v_scores);
    
  
    void displayCoordsAsLocs(std::ofstream &fd,
			     Transcript &t, 
			     vector<pair<int, int> > & v);
          
    inline void displayGenes(std::string  &fmt, std::ofstream &fd, 
			     list<Gene> &genes,
			     int offset = 0,
			     int upperLimit = INT_MAX,
			     IntervalSet<double> *tag_scores = NULL) {

      for(list<Gene>::iterator it = genes.begin(); it != genes.end(); it++) {
	Gene & gene = (*it);
	displayGene(fmt, fd, gene, offset, upperLimit, tag_scores);
      }

    }

    inline void displayGenes(std::string  &fmt, std::ofstream &fd, 
			     list<Gene *> &genes,
			     int offset = 0,
			     int upperLimit = INT_MAX,
			     IntervalSet<double> *tag_scores = NULL) {

      for(list<Gene *>::iterator it = genes.begin(); it != genes.end(); it++) {
	Gene & gene = *(*it);
	displayGene(fmt, fd, gene, offset, upperLimit, tag_scores);
      }

    }

    inline void displayGene(std::string  &fmt, std::ofstream &fd, 
			    Gene &gene,
			    int offset = 0,
			    int upperLimit = INT_MAX,
			    IntervalSet<double> *tag_scores = NULL) {
      
      for(unsigned i = 0; i < gene.transcripts().size(); i++) {
        Transcript & t = gene.transcripts()[i];
	displayTranscript(fmt, fd, gene.getId(), t, 
			  offset, upperLimit, tag_scores);
      }

    }

    inline void displayTags(std::string & fmt,
                            std::ofstream & os,
                            Sequence & c,
                            SeqTags & seqTags, 
                            FSM &fsm,
                            const char *idPrefix,
                            int offset = 0,
                            int upperLimit = INT_MAX) {
      
      list<Gene> genes;
      GeneUtils::tags2Genes(c, seqTags, idPrefix, fsm, 
			    genes, offset, reportNonCoding);
      displayGenes(fmt, os, genes, offset, upperLimit);
      os.flush();
    }

    inline void displayTags(std::string & fmt,
                            std::ofstream & os,
                            Sequence & c,
                            FSM &fsm,
                            const char *idPrefix,
                            TSetType set,
                            int offset = 0,
                            int upperLimit = INT_MAX) {
      
      if(!fmt.compare("tags")) {
        TagUtils::saveSeqTags(os, c, fsm, set);
        return;
      }

      vector<SeqTags> & seqTags = c.getTags(set);
      
      list<Gene> genes;
      GeneUtils::tags2Genes(c, seqTags, idPrefix, fsm, 
			    genes, offset, reportNonCoding);
      displayGenes(fmt, os, genes, offset, upperLimit);
      os.flush();
      
    }
    
    ~GeneTagPrinter() {;}
    
  };
}

#endif
