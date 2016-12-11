#include <iostream>
#include <streambuf>
#include <fstream>
#include <string.h>
#include "GeneUtils.h"
#include "Gene.h"
#include "FSM.h"
#include "GeneTagPrinter.h"


/****************************************************************************
* GeneTagPrinter.cpp - part of the craig namespace, a genomics library
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


namespace craig {
  
  void GeneTagPrinter::displayHeader(std::string &fmt, ::ofstream &fd) {
    if(fmt.compare("gtf") == 0) {
      fd << "#gtf-format\n";
      fd << "#Output generated with CRAIG v " << CRAIG_VERSION << ".\n";
      fd << "#Copyright (C) 2003-2007 Axel E. Bernal (abernal@seas.upenn.edu)\n\n";
    }
  }

  void GeneTagPrinter::displayTranscript(std::string &fmt,
					 std::ofstream &fd,
					 std::string &gid,
					 Transcript &t,
					 int offset,
					 int upperLimit,
					 IntervalSet<double> *tag_scores) {
    
    bool fmt_is_locs = false;
    // if(tag_scores) {
    //   IntervalSet<double>::iterator it = tag_scores->begin();
    //   for( ; it != tag_scores->end(); it++) {
    // 	cout << (int)it->first.first << " " << (int)it->first.second << " " << it->second << endl;
    //   }
    // }
    if(!fmt.compare("locs"))
      fmt_is_locs = true;
    string gtf_detail = "";
    if(fmt_is_locs) {
      fd << ">";
      if(t.hasPeptide())
        fd << "PEP-";
      fd << gid << "|" << t.getId() << "\t";
    }
    else
      gtf_detail = "\t0\tgene_id \"" + gid + "\"; transcript_id \"" + gid + 
	"|" + t.getId() + "\";\n";
    
    vector<pair<int, int> > v5, v, v3;
    IntervalSet<double>::iterator vit;
    vector<double> v5scores, vscores, v3scores;

    prepareExons4Display(t, offset, upperLimit, t.FPexons(), v5,
			 v5scores, tag_scores);
    prepareExons4Display(t, offset, upperLimit, t.exons(), v,
			 vscores, tag_scores);
    prepareExons4Display(t, offset, upperLimit, t.TPexons(), v3, 
			 v3scores, tag_scores);

    if(v5.size()) {
      // check if 5UTR finished in an intron
      if(v.size() && v5scores.size()) {
	pair<int, int> intron_coords(t.FPexons().back().end() + 1, 
				     t.exons()[0].begin() - 1);

	if(intron_coords.second - intron_coords.first >= 0) {
	  vit = tag_scores->find(intron_coords);
	  if(vit == tag_scores->end())
	    throw EXCEPTION(INCONSISTENT_CONTIG, "cannot find score for 5UTR intron");
	  v5scores.back() += vit->second; // adding the intron score
	}
      }

      if(fmt_is_locs) {
	displayCoordsAsLocs(fd, t, v5);
	fd << "[";
      }
      else 
	displayCoordsAsGTF(fd, t, -1, gtf_detail, v5, v5scores);
    }
    
    if(v.size()) {
      if(fmt_is_locs) {
	if(!t.hasStart()) {
	  if(t.phase() != 0)
	    fd << t.phase();
	  fd << "<";
	}
	displayCoordsAsLocs(fd, t, v);
	if(!t.hasStop())
	  fd << ">";
      }
      else {
	int startBeg = v[0].first, startEnd = v[0].first + 2;
	
	if(t.getStrand() == STRAND_COMP) {
	  startBeg = v[0].first - 2; 
	  startEnd = v[0].first;
	}
	
	if(t.hasPeptide())
	  fd << t.getSequenceId() << "\tCRAIG\tsignal_peptide\t.\t" << 
	    startBeg << "\t" << startEnd << "\t.\t" << gtf_detail;

	if(t.hasStart())
	  fd << t.getSequenceId() << "\tCRAIG\tstart_codon\t0\t" << 
	    startBeg << "\t" << startEnd << "\t.\t" << gtf_detail;
	
	displayCoordsAsGTF(fd, t, t.phase(), gtf_detail, v, vscores);

	int stopBeg = v[v.size() - 1].second + 1;
	int stopEnd = v[v.size() - 1].second + 3;
	
	if(t.getStrand() == STRAND_COMP) {
	  stopBeg = stopEnd - 3; 
	  stopEnd = v[v.size() - 1].second - 1;
	}
	
	fd << t.getSequenceId() << "\tCRAIG\tstop_codon\t0\t" << 
	  stopBeg << "\t" << stopEnd << "\t.\t" << gtf_detail;
      }
    }

    if(v3.size())  {
      // check if 3UTR started with an intron
      if(v.size() && v3scores.size()) {
	pair<int, int> intron_coords(t.exons().back().end() + 1,
				     t.TPexons()[0].begin() - 1);
	
	if(intron_coords.second - intron_coords.first >= 0) {
	  vit = tag_scores->find(intron_coords);
	  if(vit == tag_scores->end())
	    throw EXCEPTION(INCONSISTENT_CONTIG, "cannot find score for 3UTR intron");	  
	  v3scores[0] += vit->second; // adding the intron score
	
	}
      }

      if(fmt_is_locs) {
	fd << "]";
	displayCoordsAsLocs(fd, t, v3);
      }
      else
	displayCoordsAsGTF(fd, t, -1, gtf_detail, v3, v3scores);
    }

    if(fmt_is_locs) {
      int i = 0;    
      if(v5scores.size() || vscores.size() || v3scores.size())
	fd << "\t";
      for( ; i < v5scores.size(); i++)
	fd << (i > 0 ? ";" : "") << v5scores[i];
      fd << (v5scores.size() ? "[" : "");

      for(i = 0; i < vscores.size(); i++)
	fd << (i > 0 ? ";" : "") << vscores[i];      

      fd << (v3scores.size() ? "]" : "");
      for(i = 0; i < v3scores.size(); i++)
	fd << (i > 0 ? ";" : "") << v3scores[i];

      fd << endl;
    }
  }
  
  void GeneTagPrinter::prepareExons4Display(Transcript &t,
					    int offset,
					    int upperLimit,
					    vector<Exon> & exons,
					    vector<pair<int, int> > &v,
					    vector<double> &vscores,
					    IntervalSet<double> *tag_scores) {

    
    int beg, end;
    unsigned int numExons = exons.size();    
    for(unsigned int j = 0; j < numExons; j++) {
      beg = exons[j].begin();
      end = exons[j].end();
      
      if(tag_scores) {
	pair<int, int> exon_coords(beg, end);
	IntervalSet<double>::iterator vit = tag_scores->find(exon_coords);
	if(vit == tag_scores->end())
	  throw EXCEPTION(INCONSISTENT_CONTIG, "cannot find score for exon");

	vscores.push_back(vit->second);

	if(j > 0 && numExons > 1) {
	  pair<int, int> intron_coords(exons[j-1].end() + 1, beg - 1);
	  vit = tag_scores->find(intron_coords);
	  if(vit == tag_scores->end())
	    throw EXCEPTION(INCONSISTENT_CONTIG, "cannot find score for intron");	  
	  vscores.back() += vit->second; // adding the intron score
	}
      }

      if(t.getStrand() == STRAND_COMP) {
	if(t.isInOneStrand()) {
	  beg = exons[numExons - j - 1].end();
	  end = exons[numExons - j - 1].begin();
	}
	else {
	  beg = (t.getSequence())->length() - beg + 1;
	  end = (t.getSequence())->length() - end + 1;
	}
      }
      
      beg += offset;
      end += offset;
      
      if(beg > upperLimit || end > upperLimit)
	break;
      
      v.push_back(pair<int, int>(beg, end));
    }

    assert(!vscores.size() || vscores.size() == v.size());

  }

  
  void GeneTagPrinter::displayCoordsAsGTF(std::ofstream &fd,
					  Transcript &t, 
					  int initPhase,
					  std::string & gtf_detail,
					  vector<pair<int, int> > & v,
					  vector<double> & vscores) {
    
    unsigned int tLen = initPhase;
    string strand = (t.getStrand() == STRAND_COMP ? "-" : "+");

    for(int i = 0; i < v.size(); i++) {
      fd << t.getSequenceId() << "\tCRAIG\t" << 
	(initPhase >= 0 ? "CDS" : "UTR") << "\t" << 
	((t.getStrand() == STRAND_COMP) ? v[i].second : v[i].first) << "\t" << 
	((t.getStrand() == STRAND_COMP) ? v[i].first : v[i].second) << 
	"\t";
      if(vscores.size())  fd << vscores[i];  else fd << ".";
      fd << "\t" << strand << "\t";

      if(initPhase >= 0)  fd << (3 - (tLen % 3)) % 3;  else fd << ".";
      fd << gtf_detail;	
      
      if(initPhase >= 0)
	tLen = v[i].second - v[i].first + 1;
    }
  }
  
  void GeneTagPrinter::displayCoordsAsLocs(std::ofstream &fd,
					   Transcript &t, 
					   vector<pair<int, int> > & v) {

    for(int i = 0; i < v.size(); i++) {    
      if(i > 0)
	fd << ";";
      fd << t.getSequenceId() << "_" << v[i].first << "_" << v[i].second;

    }
  }
}
