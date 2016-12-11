#include "ResourceEngine.h"
#include "Motif.h"
#include <boost/regex.hpp>
#include "ContextIMM.h"
#include "FSM.h"
#include "ParamModel.h"
#include "InpFile.h"

/****************************************************************************
* ResourceEngine.cpp - part of the lless namespace, a general purpose
*                      linear semi-markov structure prediction library
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

  /**
   * Registers all Resource objects that can be created dynamically in the 
   * Resource Header definition file.
   */
  void ResourceEngine::registerDefaultResources() {
    if(!resourceFactory.Register(std::string("Sigma"), Type2Type<Sigma>()))
      assert(0);
    if(!resourceFactory.Register(std::string("DNASigma"), Type2Type<DNASigma>()))
      assert(0);
    if(!resourceFactory.Register(std::string("AASigma"), Type2Type<AASigma>()))
      assert(0);
    if(!resourceFactory.Register(std::string("ContextIMM"), Type2Type<ContextIMM>()))
        assert(0);
    if(!resourceFactory.Register(std::string("DBSoftMatchMotif"), Type2Type<DBSoftMatchMotif>()))
      assert(0);
    if(!resourceFactory.Register(std::string("TransfacMotif"), Type2Type<TransfacMotif>()))
      assert(0);
    if(!resourceFactory.Register(std::string("FSM"), Type2Type<FSM>()))
      assert(0);
    if(!resourceFactory.Register(std::string("File"), Type2Type<InpFile>()))
      assert(0);
    if(!resourceFactory.Register(std::string("Dir"), Type2Type<InpDir>()))
      assert(0);
    if(!resourceFactory.Register(std::string("FileMerger"), Type2Type<InpFileMerger>()))
      assert(0);
    if(!resourceFactory.Register(std::string("BitFileCollapser<int>"), Type2Type<InpBitFileCollapser<MultiScoreSeq<int>, ScoreSeq<int> > >()))
      assert(0);
    if(!resourceFactory.Register(std::string("PhaseBitFileCollapser<int>"), Type2Type<InpPhaseBitFileCollapser<MultiScoreSeq<int>, ScoreSeq<int> > >()))
      assert(0);
    if(!resourceFactory.Register(std::string("FastaFileCollapser<int>"), Type2Type<InpFastaFileCollapser<MultiScoreSeq<int>, ScoreSeq<int> > >()))
      assert(0);
    if(!resourceFactory.Register(std::string("FastaFileCollapser<UCHAR>"), Type2Type<InpFastaFileCollapser<MultiScoreSeq<UCHAR>, ScoreSeq<UCHAR> > >()))
      assert(0);
    if(!resourceFactory.Register(std::string("FastaFileCollapser<USHORT>"), Type2Type<InpFastaFileCollapser<MultiScoreSeq<USHORT>, ScoreSeq<USHORT> > >()))
      assert(0);
    if(!resourceFactory.Register(std::string("SparseFastaFileCollapser<UCHAR>"), Type2Type<InpFastaFileCollapser<MultiSparseSeq<UCHAR>, SparseSeq<UCHAR> > >()))
      assert(0);
    if(!resourceFactory.Register(std::string("SparseFastaFileCollapser<USHORT>"), Type2Type<InpFastaFileCollapser<MultiSparseSeq<USHORT>, SparseSeq<USHORT> > >()))
      assert(0);
    if(!resourceFactory.Register(std::string("SparseFastaFileCollapser<int>"), Type2Type<InpFastaFileCollapser<MultiSparseSeq<int>, SparseSeq<int> > >()))
      assert(0);
    if(!resourceFactory.Register(std::string("PhaseFastaFileCollapser<int>"), Type2Type<InpPhaseFastaFileCollapser<MultiScoreSeq<int>, ScoreSeq<int> > >()))
      assert(0);
    if(!resourceFactory.Register(std::string("ParamModel"), Type2Type<ParamModel>()))
      assert(0);

  }

  /**
   * Reads Resource Header/Content  definitions from input file stream 
   * reStream
   * @param replace_res string to be replaced during header reading
   * @param with_res string to replace replace_with  with
   * During training time, with_res is a string representing the the
   * absolute path of the directory containing the model concatenated
   * with the name of the model, whereas with_res is the "~" string
   */
  void ResourceEngine::readHeaderDefinitions(std::ifstream &reStream,
					     map<string, string> *rsubs) {

    assert(reStream);
    Resource *r;
    boost::RegEx rExDelim("\\s+");
    boost::RegEx rExComment("^\\s*\\#.*");
    boost::RegEx rExEndofFile("^\\s*\\/\\/\\s*$");
    boost::RegEx rExResource("^\\s*Resource\\s+(\\S.+\\S)\\s*$");
    std::vector<std::string> reArgs;
    std::string line;
    int offset;

    while(std::getline(reStream, line), !reStream.eof()) {
      r = NULL;
      offset = 0;

      if(rExComment.Match(line)) // a comment
        continue; 

      if(line.length() == 0)
        continue;

      if(rExEndofFile.Match(line))
        break;

      if(rExResource.Match(line)) {
        ::string args(rExResource[1]);
        reArgs.clear();

	if(rsubs) {
	  map<string, string>::iterator it = rsubs->begin();
	  for( ; it != rsubs->end(); it++) 
	    Utils::substitute(args, it->first, it->second, true);
	}

	rExDelim.Split(reArgs, args, 0, 100);
        r = resourceFactory.Create(reArgs[1], reArgs, offset, this);

        if(!r) {
          cerr << "Warning! unknown resource name " << reArgs[1] << endl;
          continue;
        }
        if(r->contentsAreInline())
          r->retrieveContents(reStream);
      }
      else {
          cerr << "Resource file: Warning! unknown parsing line : " << line << endl;
          continue;
      }

      if(r == NULL || (unsigned)offset != reArgs.size()) {
        assert(0);
        throw EXCEPTION(PARSE_ERROR, line);
      }

      resourceHeaders.push_back(line);
      setResource(r);
    }
  }

}
