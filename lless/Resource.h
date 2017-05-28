/****************************************************************************
 * Resource.h - part of the lless namespace, a general purpose
 *              linear semi-markov structure prediction library
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

#ifndef _RESOURCE_H_
#define _RESOURCE_H_

#include "Utils.h"

namespace lless {

    /**
     * The Resource class is used for managing input data in lless.
     * Almost anything which should be read/gathered from either the filesystem
     * or the user should be represented by a Resource. lless contains a few
     * classes of resources, which were implemented as needs arised during
     * development, but this list is by no means comprehensive.
     *
     * Resources provide an elegant representation of any type of input to the
     * program and could be used as black boxes by filters, features or in any
     * other environment in a transparent way.
     *
     * Every Resource contains two fields which control its behaviour: the field
     * _needsPreProcessing controls whether the Resource Contents to be
     * preProcessed before they get stored or not. When a Resource object
     * is stored back in the filesystem, Contents are always preProcessed. The
     * field _contentsInline is explained in detail in class ResourceEngine
     *
     * Examples of Resource-derived objects are FMS objects, Fasta objects,
     * ContextIMM objects and so on.
     *
     *
     ***************************************************************************/

    class Resource {

    protected:
        std::string name;         //!< Resource Name
        std::string uniqueClsId;  //!< unique identifier for Object Factory
        bool _contentsInline;     //!< whether Contents are inline in the file
        bool _needsPreProcessing; //!< whether Contents need preprocessing or not

    public:
        //! Default constructor. The default name is "default-resource".
        Resource(
            //! checks whether Contents need preprocessing or not
            bool needsPreProcessing = false,
            /*! checks whether Resource Contents are right after the Header
              definition, in the resource definition file or not. */
            bool contentsInline = false
            ) {

            this->name = std::string("default-resource");
            this->uniqueClsId = std::string("Resource");
            this->_needsPreProcessing = needsPreProcessing;
            this->_contentsInline = contentsInline;
        }

        //! Contructor with Resource name.
        Resource(std::string &name,               //!< A unique resource name
                 //! checks whether Contents need preprocessing or not
                 bool needsPreProcessing = false,
                 /*! checks whether Resource Contents are right after the Header
                   definition, in the resource definition file or not. */
                 bool contentsInline = false
            ) {

            this->name = name;
            this->uniqueClsId = std::string("Resource");
            this->_needsPreProcessing = needsPreProcessing;
            this->_contentsInline = contentsInline;
        }

        /**
         * Constructor from a Header string definition.
         */
        Resource(
            //! The Header string definition, loaded as a vector of strings.
            vector<std::string> & params,
            int & offset,       //!< The index for vector params.
            /*! checks whether Resource Contents are right after the Header
              definition, in the resource definition file or not. */
            bool contentsInline = false
            ) {

            this->name = params[offset++];
            this->uniqueClsId = params[offset++];
            this->_contentsInline = contentsInline;
            this->_needsPreProcessing = false;

        }

        /**
         * @return true if Resource Contents and inline, false otherwise.
         */
        bool contentsAreInline() {
            return _contentsInline;
        }

        /**
         * @return true if Resource Contents need preprocessing, false otherwise.
         */
        bool needsPreProcessing() {
            return  _needsPreProcessing;
        }

        /**
         * @return Resource unique name
         */
        inline std::string & getName() {
            return name;
        }

        /**
         * A member function for saving the Header definition in some output file
         * Retrieving the Header is made through the constructor
         * @see Resource(std::vector<::string> &, int, bool)
         * because information about the type of the derived Resource class
         * must be known beforehand.
         *
         * @param fd the output stream
         */
        virtual inline void saveHeader(::ofstream & fd) {
            // print header
            fd << "Resource\t" << name << "\t" << uniqueClsId;
        }

        /**
         * A member function for saving the Resource Contents in some output file.
         * The Contents can only be stored once they have been preprocessed.
         * @param fd the output stream
         */
        virtual inline void saveContents(::ofstream & fd) {

        }

        /**
         * A member function for retrieving the Resource Contents from some input
         * file.
         * The Contents may need preprocessing
         * @param fd the input stream
         */
        virtual inline void retrieveContents(::ifstream & fd) {

        }

        /**
         * A member function that does 'something' with the original Contents
         * before they become usable.
         */
        virtual inline void preProcess(::ifstream &fd) {

        }

        virtual void releaseSeqContents() {

        }

        virtual ~Resource() {}
    };

}

#endif
