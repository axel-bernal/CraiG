#include "InpFile.h"
#include "InpFileUtils.h"
#include "SequenceUtils.h"
#include "dirent.h"

/****************************************************************************
 * InpFile.cpp - part of the lless namespace, a general purpose
 *               linear semi-markov structure prediction library
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
     * Constructor from a Header string definition.
     * @param params The Header string definition, loaded as a vector of
     * strings. The Header has the following form:\n\n Resource name File
     * type filename [sigma] mapIds\n\n
     * The description of the fields could be found in the other
     * constructor(s)
     * @param offset The index for vector params
     * @param re A pointer to the ResourceEngine object.
     */
    IOSource::IOSource(std::vector<std::string> &params,
                       int & offset,
                       ResourceEngine *re)
        : Resource(params, offset) {

        this->_type = InpFileUtils::stringToTFileFormat(params[offset++]);
        this->filename = params[offset++];
        sigma = NULL;

        if(isFasta()) {
            this->sigma  = (Sigma *)re->getResource(params[offset++]);
            assert(sigma);
        }

        this->mapIds = Utils::stringToBoolean(params[offset++]);
        this->loadAllSeqs = Utils::stringToBoolean(params[offset++]);
        this->_numCollapsedVals = 1;
        this->fsm = NULL;
        this->_lockContents = false;

        if(offset < params.size()) {
            if(_type != SEQTAGS) {
                if(!sscanf(params[offset++].c_str(), "%d", &this->_numCollapsedVals))
                    assert(0);
            }
            else
                this->fsm  = (FSM *)re->getResource(params[offset++]);
        }
        if(_numCollapsedVals > 1 && !isMultiScore())
            throw EXCEPTION( NOT_SUPPORTED, string("Cannot collapse more than one value into ") + this->filename);

        if(_type == SEQTAGS && !fsm)
            throw EXCEPTION(BAD_RESOURCE_HEADER, this->getName());

    }

    /**
     * \todo add sparse sequence loading routines
     */
    BasicSeq* IOSource::loadSequence(std::ifstream & annotSeqStream,
                                     std::string &line) {

        BasicSeq *c = NULL;

        switch(_type) {
        case FASTA:
            c = SequenceUtils::loadFastaSequence(annotSeqStream, sigma, line);
            break;
        case XFASTA:
            c = SequenceUtils::loadExtFastaSequence(annotSeqStream, sigma, line);
            break;
        case SEQTAGS:
            c = SequenceUtils::loadTagSequence(annotSeqStream, *fsm, line);
            break;
        case MULTIZ:
            c = SequenceUtils::loadMultizSequence(annotSeqStream, sigma, line);
            break;
        case EXTENDED_INTSCORE:
            c = SequenceUtils::loadExtScoreSequence<int>(annotSeqStream, line);
            break;
        case EDGEANNOT_FASTA:
            c = SequenceUtils::loadEdgeAnnotSequence(annotSeqStream, sigma, line);
            break;
        case EXTENDED_LDSCORE:
            c = SequenceUtils::loadExtScoreSequence<double>(annotSeqStream, line);
            break;
        case EXTENDED_MULTI_INTSCORE:
            c = SequenceUtils::loadExtMultiScoreSequence<int>(annotSeqStream, line);
            break;
        case EXTENDED_MULTI_LDSCORE:
            c = SequenceUtils::loadExtMultiScoreSequence<double>(annotSeqStream, line);
            break;
        default:
            throw EXCEPTION( NOT_SUPPORTED, string("file type cannot be loaded for ") + filename);
        }

        return c;
    }

    /**
     * @see Resource::saveHeader(::ofstream &)
     */
    void IOSource::saveHeader(std::ofstream &fd) {
        Resource::saveHeader(fd);

        fd << " " << InpFileUtils::tInpFileToString(_type) << " " << filename;

        if(sigma)
            fd << " " << sigma->getName();

        fd  << " " << Utils::booleanToString(mapIds) << " "
            << Utils::booleanToString(loadAllSeqs) << " "
            << this->_numCollapsedVals << endl;

    }

    /*
     * InpFile routines
     */

    /**
     * @param id the identification of the BasicSeq object we want.
     * @return a pointer to a BasicSeq object with identification id.
     */
    BasicSeq* IOSource::findSeq(std::string id) {

        std::map<std::string, BasicSeq*>::iterator it;
        it = seqs.find(id);

        if(it != seqs.end())
            return it->second;

        if(loadAllSeqs)
            return NULL;

        /*
         * Keep only one sequence in memory, succesive queries about
         * sequences are usually on the same sequence id so we cache it
         */
        BasicSeq *c = loadFileSequence(id);
        if(!c)  return c;

        cacheSeq(c);

        return c;

    }

    BasicSeq * IOSource::cacheSeq(BasicSeq *c) {
        assert(c);
        _seqs.push_back(c);

        if(mapIds)
            seqs[c->id()] = _seqs.back();

        return c;
    }

    void IOSource::releaseSeqContents() {

        // if all sequences are loaded or being edited, no need to uncache.
        if(loadAllSeqs || this->_lockContents)
            return;

        list<BasicSeq *>::iterator it = _seqs.begin();
        for( ; it != _seqs.end(); it++) {
            if(*it) {
                //	cerr << "unloading " << (*it)->id() << " for " << this->getName() << endl;
                delete (*it);
            }
            *it = NULL;
        }

        _seqs.clear();
        seqs.clear();

    }

    /**
     * A member function that loads a sequence contained in the fasta file
     * into memory. The sequence is created and inserted into a list and
     * mapped to a hash table if mapIds is true.
     * @param id the identifier of the sequence
     */

    BasicSeq *InpFile::loadFileSequence(std::string id) {
        std::map<std::string, streampos>::iterator it;
        it = seqIndexes.find(id);

        if(it == seqIndexes.end())
            return NULL;

        std::string line;
        fs->clear();   // unset fail flags
        fs->seekg(it->second);
        std::getline(*fs, line);
        //    cerr << "loading " << id << " for " << this->getName() << endl;
        return loadSequence(*fs, line);

    }


    /**
     * A member function that loads the sequences contained in the fasta file
     * into memory. For every read sequence, a pointer to a BasicSeq object is
     * created and inserted into a list. It won't do anything unless
     * loadAllSeqs has been set to true.
     */
    void InpFile::loadFileSequences() {

        if(!loadAllSeqs)
            return;

        std::string line;
        std::getline(*fs, line);

        while(1) {
            BasicSeq *c = loadSequence(*fs, line);
            if(!c)  break;
            sequences().push_back(c);
        }

        cerr  << sequences().size() << " BasicSeq* objects loaded from file "
              << filename << "\n";

    }



    /**
     * A member function that creates a hashmap, with the sequence identifiers
     * as they appear in the file in the header information, and the sequence
     * themselves as values. the mapping is created only if the class member
     * mapIds is true.
     */
    void InpFile::mapInputSeqs() {
        if(!mapIds)
            return;

        if(loadAllSeqs) {
            list<BasicSeq *>::iterator it = _seqs.begin();
            std::map<std::string, BasicSeq*>::iterator it2;

            for( ; it != _seqs.end(); it++) {
                BasicSeq *c = *it;
                it2 = seqs.find(c->id());

                if(it2 != seqs.end())
                    throw EXCEPTION( NOT_SUPPORTED, string("duplicated identifiers for file ") + filename);

                seqs[c->id()] = c;
                _seqNames.push_back(c->id());
            }

            return;

        }

        SequenceUtils::loadSeqIndexes(*fs, seqIndexes);

        map<string, streampos>::iterator it = seqIndexes.begin();
        for( ;it != seqIndexes.end(); it++)
            _seqNames.push_back(it->first);

    }


    /*
     * InpDir routines
     */

    /**
     * A member function that loads a sequence contained in the fasta file
     * into memory. The sequence is created and inserted into a list and
     * mapped to a hash table if mapIds is true.
     * @param id the identifier of the sequence
     */

    BasicSeq *InpDir::loadFileSequence(std::string id) {
        std::string seqFile(filename + "/");
        seqFile += id;
        std::ifstream fs(seqFile.c_str());

        if(!fs)
            return NULL;

        std::string line;
        std::getline(fs, line);

        return loadSequence(fs, line);

    }

    /**
     * A member function that creates a hashmap, with the sequence identifiers
     * as they appear in the file in the header information, and the sequence
     * themselves as values. the mapping is created only if the class member
     * mapIds is true.
     */
    void InpDir::mapInputSeqs() {
        struct dirent **dirList;
        int numFiles = scandir(filename.c_str(), &dirList, 0, 0);

        if (numFiles < 0)
            throw EXCEPTION(FILE_UNAVAILABLE, filename + "'s scandir failed\n");


        while(numFiles--) {
            std::string entry(dirList[numFiles]->d_name);
            free(dirList[numFiles]);

            if(!entry.compare(".") || !entry.compare(".."))
                continue;

            _seqNames.push_back(entry);

        }

        free(dirList);

    }

    /*
     * FileMerger routines
     */

    /**
     * Constructor from a Header string definition.
     * @param params The Header string definition, loaded as a vector of
     * strings. The Header has the following form:\n\n Resource name
     * FileMerger fillChar numInpFiles InpFileobj1 [..]\n\n
     * The description of the fields could be found in the other
     * constructor(s)
     * @param offset The index for vector params
     * @param re A pointer to the ResourceEngine object.
     */
    InpFileMerger::InpFileMerger(std::vector<std::string> &params,
                                 int & offset,
                                 ResourceEngine *re) :
        IOSource(params, offset, re, FASTA, true) {

        if(!sscanf(params[offset++].c_str(), "%c", &this->fillChar))
            assert(0);

        int numInpFiles;

        if(!sscanf(params[offset++].c_str(), "%d", &numInpFiles))
            assert(0);

        for(int i = 0; i < numInpFiles; i++) {
            InpFile *file = (InpFile *)re->getResource(params[offset++]);
            files2Merge.push_back(file);
        }

        this->overrideSigma = (this->sigma ? true : false);

        loadFileSequences();
        mapInputSeqs();

    }


    void InpFileMerger::releaseSeqContents() {
        // freeing contained files
        for(int i = 0; i < files2Merge.size(); i++)
            files2Merge[i]->releaseSeqContents();

        IOSource::releaseSeqContents();
    }

    /**
     * @see Resource::saveHeader(::ofstream &)
     */
    void InpFileMerger::saveHeader(std::ofstream &fd) {
        Resource::saveHeader(fd);

        fd << " " << (this->overrideSigma ? this->sigma->getName() : "null")
           << (this->loadAllSeqs ? " true " : " false ") << fillChar
           << " " << files2Merge.size();

        for(int i = 0; i < files2Merge.size(); i++)
            fd  << " " << files2Merge[i]->getName();

        fd << endl;

    }


    /**
     * A member function that merges the sequences with name id into
     * *this. For each sequence, and position within the sequence the character
     * from this' sequence is yanked if it is a fillChar and file's
     * character at that position is not. Disagreements between characters
     * different from fillChar will throw an exception
     */
    BasicSeq *InpFileMerger::loadFileSequence(std::string id) {
        TStrand strand;
        int fileNo, i, t;
        Sequence *target = NULL;

        for(fileNo = 0; fileNo < files2Merge.size(); fileNo++) {
            Sequence *c = (Sequence *)files2Merge[fileNo]->findSeq(id);

            if(!c) {
                if(target)  delete target;
                return NULL;
            }

            if(!target) {
                target = new Sequence(*c);
                continue;
            }

            //      cerr << files[fileNo]->fileName() << endl;

            for(t = 0; t < NUM_STRANDS; t++) {
                strand = (TStrand)t;
                char *thisSeq = target->getSeq(strand);
                char *fileSeq = c->getSeq(strand);

                assert(target->length(strand) == c->length(strand));

                for(i = 0; i < c->length(); i++) {
                    if(fileSeq[i] != fillChar) {
                        if(thisSeq[i] != fillChar && thisSeq[i] != fileSeq[i]) {
                            //              assert(0);
                            ostringstream istring; istring << t << ":" << i;
                            cerr << files2Merge[fileNo]->fileName() << " contains discordance " <<
                                fileSeq[i] << " instead of " << thisSeq[i] << " for " << c->id()
                                 << " " << t << " " << i + 1<< endl;
//              throw EXCEPTION(STRANGE_CHAR, string("InpFileMerger::merging ") + file->fileName() + string(":") + c->id() + string(" @pos ") + istring.str());
                        }
                        else {
                            // yank thisSeq[i]
                            //	      cerr << c->id() << " " << t << " " << i + 1 << " " << fileSeq[i] << endl;
                            thisSeq[i] = fileSeq[i];
                        }
                    }
                }
            }
        }

        return target;

    }

    /**
     * A member function that merges sequences contained in files into *this'.
     */
    void InpFileMerger::loadFileSequences() {

        assert(files2Merge.size());
        int i;

        if(!overrideSigma) {
            for(i = 0; i < files2Merge.size(); i++) {
                if(!this->alphabet())
                    this->sigma = files2Merge[i]->alphabet();

                if(!this->alphabet() || this->alphabet() != files2Merge[i]->alphabet()) {
                    assert(0);
                    throw EXCEPTION( BAD_USAGE, string("InpFileMerger::collapse(alphabet)"));
                }
            }
        }

        if(!loadAllSeqs)
            return;

        vector<std::string> &ids = files2Merge[0]->seqNames();

        if(!ids.size()) {
            assert(0);
            throw new EXCEPTION(BAD_USAGE, files2Merge[0]->fileName() + " must have mapped ids");
        }

        for(i = 0; i < ids.size(); i++) {
            BasicSeq *c = loadFileSequence(ids[i]);
            if(!c) {
                assert(0);
                throw EXCEPTION(CONTIG_UNAVAILABLE, ids[i] + " not found for merging");
            }
            sequences().push_back(c);
        }

        cerr  << sequences().size() << " BasicSeq* objects merged\n";
    }


    /**
     * A member function that creates a hashmap, with the sequence identifiers
     * as they appear in the file in the header information, and the sequence
     * themselves as values. the mapping is created only if the class member
     * mapIds is true.
     */
    void InpFileMerger::mapInputSeqs() {
        if(!mapIds)
            return;

        _seqNames = files2Merge[0]->seqNames();

        if(loadAllSeqs) {
            list<BasicSeq *>::iterator it = _seqs.begin();
            std::map<std::string, BasicSeq*>::iterator it2;
            for( ; it != _seqs.end(); it++) {
                BasicSeq *c = *it;
                it2 = seqs.find(c->id());
                if(it2 != seqs.end())
                    throw EXCEPTION( NOT_SUPPORTED, string("duplicated identifiers for file ") + filename);

                seqs[c->id()] = c;
            }
            return;
        }
    }
}
