/****************************************************************************
 * ParamModel.h - part of the lless namespace, a general purpose
 *                   linear semi-markov structure prediction library
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

#ifndef _PARAM_MODEL_H_
#define _PARAM_MODEL_H_

#include "Utils.h"
#include "FeatureEngine.h"
#include <stdio.h>

namespace lless {

    class ParamModel : public Resource {
    private:
        ::ifstream ipStream;
        std::ostringstream tmpFile;

    protected:
        ResourceEngine *re;
        FilterEngine *fe;
        GlobalVector *params;
        FeatureEngine *fte;
        ::ifstream *paramStream;
        bool readFeaturesOnly, saveFeaturesOnly, preComputed[NUM_STRANDS];
        std::string modelPath;
        std::string prefixFiles;
    public:
    ParamModel(
        std::string & paramFile,
        std::string & modelPath,
        std::string & prefixFiles,
        bool readFeaturesOnly = false,
        bool saveFeaturesOnly = false
        ) : Resource(false) {

            this->modelPath = modelPath;
            this->prefixFiles = prefixFiles;
            this->readFeaturesOnly = readFeaturesOnly;
            this->saveFeaturesOnly = saveFeaturesOnly;
            initialize(paramFile);

        }


    ParamModel(
        std::vector<std::string> &params,
        int & offset,       //!<The index for vector params.
        ResourceEngine *re  //!<A pointer to the ResourceEngine object.
        )
        : Resource(params, offset) {

            std::string paramFile = params[offset++];
            this->modelPath = params[offset++];
            this->prefixFiles = params[offset++];
            this->readFeaturesOnly = Utils::stringToBoolean(params[offset++]);
            this->saveFeaturesOnly = Utils::stringToBoolean(params[offset++]);
            initialize(paramFile);

        }


        void initialize(std::string & paramFile) {

            preComputed[STRAND_FWD] = preComputed[STRAND_COMP] = false;
            re = NULL;
            fe = NULL;
            tmpFile.str("");

            if(paramFile.compare("null") == 0)
                _contentsInline = true; // model to be read inline @ ResourceEngine scope
            else {
                ipStream.open(paramFile.c_str(), ios::in);

                if(!ipStream.is_open())
                    throw EXCEPTION(FILE_UNAVAILABLE, paramFile);

                paramStream = &ipStream;

                if(!readFeaturesOnly) {
                    map<string, string> rsubs;
                    rsubs["PREFIX_EVIDENCE"] = prefixFiles;
                    rsubs["~"] = modelPath;
                    this->re = new ResourceEngine(ipStream, &rsubs);
                    this->fe = new FilterEngine(*re, ipStream);
                }
            }
        }


        inline void readParams(FSM *fsm, FilterEngine *fe = NULL) {
            if(!paramStream->is_open()) // parameters were already read
                return;

            FilterEngine *currfe = (fe ? fe : this->fe);

            if(!currfe)
                throw EXCEPTION(BAD_USAGE, getName() + " needs a FilterEngine* object");

            this->fte = new FeatureEngine(*fsm, *currfe, *paramStream);
            this->params = new GlobalVector(fte->getFeatures(), paramStream);
            this->fte->setParamVector(params);

            paramStream->close();
            if(tmpFile.str().length()) {
                if(remove(tmpFile.str().c_str()))
                    throw EXCEPTION(FILE_UNAVAILABLE, tmpFile.str());
                tmpFile.str("");
            }
        }

        /**
         * @see Resource::saveHeader(std::ofstream &)
         */
        inline void saveHeader(::ofstream &fd) {
            Resource::saveHeader(fd);
            fd << " null ~ " << Utils::booleanToString(saveFeaturesOnly) << " "
               << Utils::booleanToString(saveFeaturesOnly) << " " << endl;
        }

        /**
         * @see Resource::saveContents(::ofstream &)
         */
        void saveContents(::ofstream &fd) {
            if(!saveFeaturesOnly) {
                if(!re || !fe)
                    throw EXCEPTION(BAD_USAGE, getName() + " cannot save NULL content");

                re->saveResources(fd);
                fe->saveFilters(fd);

            }

            fte->saveFeatures(fd);
            params->store(fd);
            fd << "//" << endl;
        }

        /**
         * @see Resource::retrieveContents(std::ifstream &)
         */
        inline void retrieveContents(std::ifstream &fd) {

            tmpFile << modelPath << this->getName() << "." << getpid();
            ::ofstream opStream(tmpFile.str().c_str(), ios::out);

            if(!opStream.is_open())
                throw EXCEPTION(FILE_UNAVAILABLE, tmpFile.str());

            boost::RegEx rExEndofFile("^\\s*\\/\\/\\s*$");
            std::string line;
            bool endParamFile = false;
            while(std::getline(fd, line), !fd.eof()) {

                if(rExEndofFile.Match(line))
                    if(endParamFile)  break;
                    else endParamFile = true;

                opStream << line << endl;

            }

            opStream.close();

            ipStream.open(tmpFile.str().c_str(), ios::in);
            paramStream = &ipStream;

        }

        inline void preComputeFeatures(int beg, int end, TStrand strand) {

            if(preComputed[strand])
                return;

            fte->preComputeFeatures(beg, end, strand);
            preComputed[strand] = true;

        }

        inline void deletePreComputedFeatures() {

            if(!preComputed[STRAND_FWD] &&
               !preComputed[STRAND_COMP])
                return;

            fte->deletePreComputedFeatures();
            preComputed[STRAND_FWD] = preComputed[STRAND_COMP] = false;

        }

        inline FeatureEngine * featureEngine() {
            return fte;
        }

        ~ParamModel() {
            if(re)  delete re;
            if(fe)  delete fe;

            delete params;
            delete fte;
        }
    };
}

#endif
