#include "FilterEngine.h"
#include <boost/regex.hpp>
#include "FL_Context.h"
#include "FL_Signal.h"
#include "FL_Score.h"
#include "FL_CountingScore.h"
#include "FL_DiffCountingScore.h"
#include "FL_ImmScore.h"
#include "FL_MotifScore.h"
#include "FL_SpacerMotifScore.h"
#include "FL_FastaxFilter.h"
#include "FL_ScoreFile.h"
#include "FL_Reranker.h"
#include "FL_EvdEdgeAligner.h"
#include "FL_Coverage.h"
#include "FL_OrderStatistics.h"
#include "FL_MovingQuantile.h"
#include "FL_ChangePoint.h"
#include "FL_ChangePtSignal.h"
#include "FL_ChangePtScore.h"

/****************************************************************************
 *   FilterEngine.cpp - part of the lless namespace, a general purpose
 *                    linear semi-markov structure prediction library
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

    FilterEngine::FilterEngine(ResourceEngine &re, std::ifstream & fd,
                               std::string **sig_info, int num_sigfilters) {

        assert(fd);
        int i;
        vector<std::string> flines;
        std::string fline;
        vector<std::string> ffields;
        while(std::getline(fd, fline), !fd.eof()) {
            istringstream ifline(fline);
            if(!fline.length())  continue;

            std::string first_tok, second_tok;
            ifline >> first_tok >> second_tok;
            if(!first_tok.compare("//"))  {
                flines.push_back(fline);
                break;
            }
            if(first_tok[0] == '#')  continue;

            bool skip = false;
            for(i = 0; i < num_sigfilters; i++) {
                //	cerr << sig_info[0][i] << " " << second_tok+"-input\n";
                if(sig_info[2][i].length() && !second_tok.compare(sig_info[0][i]+"-input"))
                    skip = true;
            }

            //      cerr << fline << " " << skip << endl;

            if(skip)
                continue;

            for(i = 0; i < num_sigfilters; i++) {
                if(!sig_info[2][i].length() || sig_info[0][i].compare(second_tok))
                    continue;

                IOSource *sig_res = (IOSource *)re.getResource(sig_info[2][i]);
                assert(sig_res);
                std::string &sigma_name = sig_res->alphabet()->getName();
                std::string new_line = std::string("LazyFilter\t")+sig_info[0][i]+"\tEdgeInst\t"+TypeDefs::tEdgeId2ToString((TEdgeId2)i)+" "+sigma_name+" 1 "+sig_info[1][i];
                flines.push_back(new_line);
                fline = std::string("Filter\t")+sig_info[0][i]+"-input\tFastaxFilter<char,EdgeInst> "+sig_info[2][i]+" "+sig_info[0][i];
            }

            flines.push_back(fline);
        }

        initialize(re, flines);
        //    for(int k = 0; k < flines.size(); k++)
        //      cerr << flines[k] << endl;

    }


    /**
     * Registers all Filter objects that can be created dynamically in the Filter
     * Header definition file
     */
    void FilterEngine::registerDefaultFilters() {
        if(!filterFactory.Register(std::string("Content<double>"), Type2Type<FL_Content<double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FilterContent<double>"), Type2Type<FL_FilterContent<double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("Content<UCHAR>"), Type2Type<FL_Content<UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("Bin<int>"), Type2Type<FL_Bin<int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("Bin<double>"), Type2Type<FL_Bin<double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("Context"), Type2Type<FL_Context>()))
            assert(0);
        if(!filterFactory.Register(std::string("ImmScore"), Type2Type<FL_ImmScore>()))
            assert(0);
        if(!filterFactory.Register(std::string("Gram<int>"), Type2Type<FL_Gram<int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("Gram<UCHAR>"), Type2Type<FL_Gram<UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("Gram<USHORT>"), Type2Type<FL_Gram<USHORT> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ConsGram<int>"), Type2Type<FL_ConsGram<int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ConsGram<UCHAR>"), Type2Type<FL_ConsGram<UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastGram<int,int>"), Type2Type<FL_FastGram<int, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastGram<int,USHORT>"), Type2Type<FL_FastGram<int, USHORT> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastGram<UCHAR,UCHAR>"), Type2Type<FL_FastGram<UCHAR, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastGram<UCHAR,int>"), Type2Type<FL_FastGram<UCHAR, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastGram<UCHAR,USHORT>"), Type2Type<FL_FastGram<UCHAR, USHORT> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastGram<USHORT,USHORT>"), Type2Type<FL_FastGram<USHORT, USHORT> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastGram<USHORT,int>"), Type2Type<FL_FastGram<USHORT, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastGram<USHORT,ULONG>"), Type2Type<FL_FastGram<USHORT, ULONG> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseGram<UCHAR,UCHAR>"), Type2Type<FL_SparseGram<UCHAR, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseGram<USHORT,USHORT>"), Type2Type<FL_SparseGram<USHORT, USHORT> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastGram<UCHAR,UCHAR#>"), Type2Type<FL_FastGram<UCHAR, SPARSE_HASH<UCHAR, int> > >()))
            assert(0);
        if(!filterFactory.Register(std::string("EdgeInst"), Type2Type< FL_Signal<EdgeInst> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FilterEdgeInst"), Type2Type< FL_FilterSignal<EdgeInst> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseFastaxFilter<int*,int>"), Type2Type< FL_FastaxFilter<MultiSparseSeq<int>, SparseSeq<int>, int, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseFastaxFilter<UCHAR*,int>"), Type2Type< FL_FastaxFilter<MultiSparseSeq<UCHAR>, SparseSeq<UCHAR>, int, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseFastaxFilter<UCHAR*,UCHAR>"), Type2Type< FL_FastaxFilter<MultiSparseSeq<UCHAR>, SparseSeq<UCHAR>, UCHAR, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseFastaxFilter<UCHAR*,USHORT>"), Type2Type< FL_FastaxFilter<MultiSparseSeq<UCHAR>, SparseSeq<UCHAR>, USHORT, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseFastaxFilter<USHORT*,USHORT>"), Type2Type< FL_FastaxFilter<MultiSparseSeq<USHORT>, SparseSeq<USHORT>, USHORT, USHORT> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseFastaxFilter<USHORT*,int>"), Type2Type< FL_FastaxFilter<MultiSparseSeq<USHORT>, SparseSeq<USHORT>, int, USHORT> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseFastaxFilter<USHORT*,ULONG>"), Type2Type< FL_FastaxFilter<MultiSparseSeq<USHORT>, SparseSeq<USHORT>, ULONG, USHORT> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseFastaxFilter<int*,USHORT>"), Type2Type< FL_FastaxFilter<MultiSparseSeq<int>, SparseSeq<int>, USHORT, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseFastaxFilter<UCHAR*,UCHAR#>"), Type2Type< FL_FastaxFilter<MultiSparseSeq<UCHAR>, SparseSeq<UCHAR>, SPARSE_HASH<UCHAR, int>, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseFastaxFilter<int*,UCHAR>"), Type2Type< FL_FastaxFilter<MultiSparseSeq<int>, SparseSeq<int>, UCHAR, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<char,UCHAR>"), Type2Type< FL_FastaxFilter<Sequence, Sequence, UCHAR, char> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<int,int>"), Type2Type< FL_FastaxFilter<ScoreSeq<int>, ScoreSeq<int>, int, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<char,int>"), Type2Type< FL_FastaxFilter<ScoreSeq<int>, ScoreSeq<int>, int, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<char,EdgeInst>"), Type2Type<FL_FastaxFilter<Sequence, Sequence, EdgeInst, char> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<int*,int>"), Type2Type< FL_FastaxFilter<MultiScoreSeq<int>, ScoreSeq<int>, int, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<UCHAR*,int>"), Type2Type< FL_FastaxFilter<MultiScoreSeq<UCHAR>, ScoreSeq<UCHAR>, int, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<UCHAR*,UCHAR>"), Type2Type< FL_FastaxFilter<MultiScoreSeq<UCHAR>, ScoreSeq<UCHAR>, UCHAR, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<UCHAR*,USHORT>"), Type2Type< FL_FastaxFilter<MultiScoreSeq<UCHAR>, ScoreSeq<UCHAR>, USHORT, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<USHORT*,USHORT>"), Type2Type< FL_FastaxFilter<MultiScoreSeq<USHORT>, ScoreSeq<USHORT>, USHORT, USHORT> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<int*,USHORT>"), Type2Type< FL_FastaxFilter<MultiScoreSeq<int>, ScoreSeq<int>, USHORT, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<int*,UCHAR>"), Type2Type< FL_FastaxFilter<MultiScoreSeq<int>, ScoreSeq<int>, UCHAR, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<UCHAR*,UCHAR#>"), Type2Type< FL_FastaxFilter<MultiScoreSeq<UCHAR>, ScoreSeq<UCHAR>, SPARSE_HASH<UCHAR, int>, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ScoreFile<int*,int>"), Type2Type< FL_ScoreFile<MultiScoreSeq<int>, ScoreSeq<int>, int, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ScoreFile<UCHAR*,UCHAR>"), Type2Type< FL_ScoreFile<MultiScoreSeq<UCHAR>, ScoreSeq<UCHAR>, UCHAR, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ScoreFile<USHORT*,USHORT>"), Type2Type< FL_ScoreFile<MultiScoreSeq<USHORT>, ScoreSeq<USHORT>, USHORT, USHORT> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ScoreFile<UCHAR*,int>"), Type2Type< FL_ScoreFile<MultiScoreSeq<UCHAR>, ScoreSeq<UCHAR>, int, UCHAR> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ScoreFile<double*,double>"), Type2Type< FL_ScoreFile<MultiScoreSeq<double>, ScoreSeq<double>, double, double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ScoreFile<int,int>"), Type2Type< FL_ScoreFile<ScoreSeq<int>, ScoreSeq<int>, int, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("SparseScoreFile<int*,int>"), Type2Type< FL_ScoreFile<MultiSparseSeq<int>, SparseSeq<int>, int, int> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ScoreFile<double,double>"), Type2Type< FL_ScoreFile<ScoreSeq<double>, ScoreSeq<double>, double, double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("Coverage<double>"), Type2Type< FL_Coverage<double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("NormalizedScoreFile<double,double>"), Type2Type< FL_NormalizedScoreFile<ScoreSeq<double>, ScoreSeq<double>, double, double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("CountingScore"), Type2Type< FL_CountingScore >()))
            assert(0);
        if(!filterFactory.Register(std::string("MultiCountingScore"), Type2Type< FL_MultiCountingScore >()))
            assert(0);
        if(!filterFactory.Register(std::string("ContainScore"), Type2Type< FL_ContainScore >()))
            assert(0);
        if(!filterFactory.Register(std::string("DiffCountingScore"), Type2Type< FL_DiffCountingScore >()))
            assert(0);
        if(!filterFactory.Register(std::string("RelPeptideOffset"), Type2Type< FL_RelPeptideOffset >()))
            assert(0);
        if(!filterFactory.Register(std::string("RelPeptideScore"), Type2Type< FL_RelPeptideScore >()))
            assert(0);
        if(!filterFactory.Register(std::string("MotifScore"), Type2Type< FL_MotifScore >()))
            assert(0);
        if(!filterFactory.Register(std::string("SpacerMotifScore"), Type2Type< FL_SpacerMotifScore >()))
            assert(0);
        if(!filterFactory.Register(std::string("RerankOffset<double>"), Type2Type< FL_RerankOffset<double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("RerankScore<double>"), Type2Type< FL_RerankScore<double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("EdgeAligner"), Type2Type< FL_EvdEdgeAligner >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<EdgeAlignment,EdgeAlignment>"), Type2Type<FL_FastaxFilter<EdgeAnnotSeq, EdgeAnnotSeq, EvdEdges, vector<int> > >()))
            assert(0);
        if(!filterFactory.Register(std::string("OrderStatistics"), Type2Type< FL_OrderStatistics >()))
            assert(0);
        if(!filterFactory.Register(std::string("MovingQuantile<double>"), Type2Type< FL_MovingQuantile<double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<double,double*>"), Type2Type< FL_FastaxFilter<ScoreSeq<double>, ScoreSeq<double>, vector<double>, double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("FastaxFilter<double,int*>"), Type2Type< FL_FastaxFilter<ScoreSeq<double>, ScoreSeq<double>, vector<int>, double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ChangePoint<double>"), Type2Type< FL_ChangePoint<double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ChangePtSignal<double>"), Type2Type< FL_ChangePtSignal<double> >()))
            assert(0);
        if(!filterFactory.Register(std::string("ChangePtScore<double>"), Type2Type< FL_ChangePtScore<double> >()))
            assert(0);

    }

    /**
     * Reads a Filter Header definition from single line
     */
    bool FilterEngine::readFilterHeader(int & id, std::string & line) {

        boost::RegEx rExLogFilter("^\\s*log\\s*\\(\\s*(\\S+)\\s*\\)(.+)$");
        boost::RegEx rExPowFilter("^\\s*pow\\s*\\(\\s*(\\S+)\\,([0-9,\\.,\\-]+)\\)(.+)$");
        boost::RegEx rExAccFilter("^\\s*acc\\s*\\(\\s*(\\S+)\\,([0-9,\\.,\\-]+)\\)(.+)$");
        boost::RegEx rExMultiFilter("^Multi(\\S+)\\((\\S+)\\)(\\s+\\S.+)\\s*$");
        boost::RegEx rExDelim("\\s+");
        boost::RegEx rExFilter("^Filter\\s+(\\S.+\\S)\\s*$");
        boost::RegEx rExFilterName("^\\S+\\s+(\\S+)\\s+.+$");
        boost::RegEx rExSparseFilter("^\\s*Sparse(\\S.+\\S)\\s*$");
        boost::RegEx rExLazyFilter("^Lazy(\\S.+\\S)\\s*$");

        Filter *f = NULL;
        vector<std::string> headers;
        string header = line;

        bool acc_match = rExAccFilter.Match(header);
        if(acc_match)
            header = rExAccFilter[1] + rExAccFilter[3];

        bool pow_match = rExPowFilter.Match(header);
        if(pow_match)
            header = rExPowFilter[1] + rExPowFilter[3];

        bool log_match = rExLogFilter.Match(header);
        if(log_match)
            header = rExLogFilter[1] + rExLogFilter[2];

        bool multi_match = rExMultiFilter.Match(header);

        if(multi_match) {
            IOSource *ios = (IOSource *)re->getResource(rExMultiFilter[2]);

            if(!ios) {
                assert(0);
                throw EXCEPTION( PARSE_ERROR, std::string("unknown Resource ") + rExMultiFilter[2]);
            }

            header = rExMultiFilter[1] + rExMultiFilter[3];

            if(!rExFilterName.Match(header))
                throw EXCEPTION(PARSE_ERROR, header);

            vector<std::string> collapsed_fnames;

            for(int j = 0; j < ios->numCollapsedVals(); j++) {
                std::ostringstream ioString;
                ioString << j;
                std::string myHeader = header;
                Utils::substitute(myHeader, "$$", ioString.str().c_str(), true);
                headers.push_back(myHeader);
                collapsed_fnames.push_back(ios->collapsedfNames(j));
            }

            multiFilters[rExFilterName[1]] = collapsed_fnames;

        }
        else
            headers.push_back(header);

        int i = 0;
        for( ; i < headers.size(); i++) {
            bool lazy_match = rExLazyFilter.Match(headers[i]);
            bool sparse_match = rExSparseFilter.Match(headers[i]);

            if(lazy_match)
                headers[i] = rExLazyFilter[1];
            else if(sparse_match)
                headers[i] = rExSparseFilter[1];

            bool match = rExFilter.Match(headers[i]);
            if(!match)
                break;

            int offset = 0;
            ::vector<std::string> feHeaders;
            string filtDef = rExFilter[1];
            rExDelim.Split(feHeaders, filtDef, 0, 100);
            f = filterFactory.Create(feHeaders[1], id, feHeaders, offset, re, this);

            if(!f) {
                assert(0);
                throw EXCEPTION( PARSE_ERROR, std::string("unknown Filter ") + headers[i]);
            }

            if((unsigned)offset != feHeaders.size()) {
                throw EXCEPTION(PARSE_ERROR, headers[i]);
                assert(0);
            }

            double power = 1;
            if(pow_match)
                if(!sscanf(rExPowFilter[2].c_str(), "%lf", &power))
                    assert(0);

            int period = 0;
            if(acc_match)
                if(!sscanf(rExAccFilter[2].c_str(), "%d", &period))
                    assert(0);

            setFilter(feHeaders[0], f, id, power, period,
                      lazy_match, sparse_match, log_match);
            id++;
        }

        if(i && i != headers.size()) {
            assert(0);
            throw EXCEPTION( PARSE_ERROR, std::string("undefined Filter ") + headers[i]);
        }

        return (i == 0 ? false : true);

    }

    /**
     * Reads Filter Header definitions from input file stream
     * feStream
     */
    void FilterEngine::readHeaderDefinitions(vector<std::string> & lines) {

        boost::RegEx rExComment("^\\s*\\#.*");
        boost::RegEx rExDelim("\\s+");
        boost::RegEx rExEndofFile("^\\s*\\/\\/\\s*$");
        boost::RegEx rExFilterSet("^\\s*\\@define\\s+FilterSet\\s+(\\S+)\\s+(\\S+)\\s*$");
        boost::RegEx rExEndFilterSet("^\\s*\\@end\\s+FilterSet\\s*$");

        int id = 0;
        numArrayEntries = 0;

        /*
         * Load filters
         */
        for(int k = 0; k < lines.size(); k++) {
            std::string line = lines[k];
            filterHeaders += line + std::string("\n");

            if(line.length() == 0 || rExComment.Match(line))
                continue;

            if(rExEndofFile.Match(line))
                break;

            if(rExFilterSet.Match(line)) {
                pair<string, vector<string> > filterSet;
                filterSet.first = rExFilterSet[2];
                filterSets[rExFilterSet[1]] = filterSet;

                while(k++, k < lines.size()) {
                    line = lines[k];
                    filterHeaders += line + std::string("\n");

                    if(line.length() == 0 || rExComment.Match(line))
                        continue;

                    if(rExEndofFile.Match(line))
                        throw EXCEPTION(PARSE_ERROR, line);

                    if(rExEndFilterSet.Match(line)) {
                        break;
                    }
                    filterSets[rExFilterSet[1]].second.push_back(line);
                }

                continue;
            }

            bool isHeader = readFilterHeader(id, line);

            if(!isHeader) {
                ::vector<std::string> filterSetVect;

                rExDelim.Split(filterSetVect, line, 0, 100);
                map<string, pair<string, vector<string> > >
                    ::iterator it = filterSets.find(filterSetVect[0]);

                if(it != filterSets.end()) {
                    pair<string, vector<string> > &fs = it->second;
                    for(int i = 0; i < fs.second.size(); i++) {
                        std::string myline = fs.second[i];
                        Utils::substitute(myline, fs.first, filterSetVect[1], true);

                        if(!readFilterHeader(id, myline)) {
                            assert(0);
                            throw EXCEPTION(PARSE_ERROR, string("after ...") + filterHeaders.substr((int)Utils::max(0, filterHeaders.size()) - 200, 200));
                        }
                    }
                }
                else {
                    assert(0);
                    throw EXCEPTION(PARSE_ERROR, string("after ...") + filterHeaders.substr((int)Utils::max(0, filterHeaders.size()) - 200, 200));
                }
            }
        }
    }

    void FilterEngine::setFilter(std::string &name, Filter *f,
                                 int id, double power, int period,
                                 bool lazy_match,
                                 bool sparse_match,
                                 bool log_match) {

        assert(f != NULL);
        filters.push_back(f);
        namedFilters[name] = f;
        filters[id]->setPaddedValArrays(this->paddedFilterValArrays);

        if(!lazy_match && !sparse_match) {
            filters[id]->setArrayIndex(numArrayEntries);
            numArrayEntries++;
        }

        if(sparse_match)
            f->makeSparse();

        if(log_match)
            f->makeLogScaled();

        if(power != 1)
            f->setPower(power);

        if(period)
            f->setAccPeriod(period);

    }


    /*
     * Routine to compute all defined filters.
     */
    void FilterEngine::computeSeqFilters(TStrand strand) {
        assert(C);
        unsigned int i;

        for(i = 0; i < this->filters.size(); i++) {
            Filter *f = this->filters[i];

            f->allocValArrays(getSequence(), strand);
            f->computeVals(getSequence(), strand);
        }
    }


    void FilterEngine::deleteSeqFilters() {
        unsigned int i;
        if(!this->filters.size())
            return;

        for(int strand = 0; strand < NUM_STRANDS; strand++) {

            if(!this->paddedFilterValArrays[strand])
                continue;

            for(i = 0; i < this->filters.size(); i++) {
                Filter *f = this->filters[i];
                f->freeValArrays((TStrand)strand);
            }
        }
    }


    void FilterEngine::releaseSequence() {
        C = NULL;
        this->seqLen = -1;
    }

    FilterEngine::~FilterEngine() {

        if(filterHeaders.length())
            this->deleteSeqFilters();

        releaseSequence();

        for(int strand = 0; strand < NUM_STRANDS; strand++) {
            if(paddedFilterValArrays[strand])
                delete [] paddedFilterValArrays[strand];
        }

        resetFilterValArrays();

        if(filterHeaders.length())
            for(unsigned int i = 0; i < this->filters.size(); i++) {
                delete this->filters[i];
                this->filters[i] = NULL;
            }
    }
}
