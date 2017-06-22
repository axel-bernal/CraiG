/****************************************************************************
 * FL_FastaxFilter.h - part of the lless namespace, a general purpose
 *                     linear semi-markov structure prediction library
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

#ifndef _FL_FASTAxFILTER_H_
#define _FL_FASTAxFILTER_H_

#include "Filter.h"
#include "InpFile.h"
#include <boost/regex.hpp>

/**
 * A namespace created for classes belonging to lless, a general purpose linear
 * semi-markov structure prediction library.
 */
namespace lless {

    /**
     * FL_InpFile is
     * a subclass of Filter which wraps an input file so that any filter could
     * be computed on the file's content.
     *
     * In all cases, Filter objects can only compute their values on the default
     * input sequence. The FL_InpFile class makes it possible for any
     * filter to compute its values on arbitrary input sequences contained in
     * an external file, as long as these sequences are in a recognized format.
     * Hence, this class can behave as a wrapper, so that any filter could be
     * computed over alternative input sequences.
     *
     ***************************************************************************/
    template <class TClass1, class TClass2>
        class FL_InpFile : public TypedFilter<TClass2> {
    protected:
        IOSource *file;
        int _seqNo;
    public:
        //! Default constructor
    FL_InpFile(
        int fInd,          //!< A unique identifier.
        std::string &name, //!< A unique name.
        IOSource *file,      //!< An IOSource* object
        int maxFilterVals,  //!< number of filter values
        int seqNo,          //!< sequence number inside file
        TValType type = FT_DOUBLE //!< Filter value type.
        ) : TypedFilter<TClass2>(fInd,
                                 name,
                                 maxFilterVals,
                                 type) {

            this->file = file;
            this->_seqNo = seqNo;
            assert(seqNo < file->numCollapsedVals());
        }

        /**
         * Constructor from Header string definition
         */
    FL_InpFile(
        int fInd,       //!< A unique identifier.
        /*! The Header string definition, loaded as a vector of
         * strings.
         * The Header has the following form:\n\n
         * Filter name InputFile<type> arrInd file\n\n
         * The meaning of these fields are as in the default
         * constructor.
         *
         * @see FL_InpFile(int,int,std::string &,InpFile*, Filter*)
         */
        vector<std::string> & params,
        int & offset,       //!< The index for vector params.
        int maxFilterVals,  //!< Max number of filter values
        ResourceEngine *re /*!< A pointer to the ResourceEngine
                             object. */
        ) : TypedFilter<TClass2>(fInd,
                                 params,
                                 offset,
                                 maxFilterVals) {

            this->file = (IOSource *)re->getResource(params[offset++]);
            this->_seqNo = 0;
            if(this->file->type() == EXTENDED_MULTI_INTSCORE ||
               this->file->type() == EXTENDED_MULTI_LDSCORE)
                if(!sscanf(params[offset++].c_str(), "%d", &this->_seqNo))
                    assert(0);

            assert(file && this->_seqNo < file->numCollapsedVals());

        }

        /**
         * Constructor from Header string definition
         */
    FL_InpFile(
        int fInd,       //!< A unique identifier.
        /*! The Header string definition, loaded as a vector of
         * strings.
         * The Header has the following form:\n\n
         * Filter name InputFile<type> arrInd file\n\n
         * The meaning of these fields are as in the default
         * constructor.
         *
         * @see FL_InpFile(int,int,std::string &,InpFile*, Filter*)
         */
        vector<std::string> & params,
        int & offset,       //!< The index for vector params.
        int maxFilterVals,  //!< Max number of filter values
        TValType type,      //!< Filter value type
        ResourceEngine *re /*!< A pointer to the ResourceEngine
                             object. */
        ) : TypedFilter<TClass2>(fInd,
                                 params,
                                 offset,
                                 maxFilterVals,
                                 type) {

            this->file = (IOSource *)re->getResource(params[offset++]);
            this->_seqNo = 0;
            if(this->file->type() == EXTENDED_MULTI_INTSCORE ||
               this->file->type() == EXTENDED_MULTI_LDSCORE)
                if(!sscanf(params[offset++].c_str(), "%d", &this->_seqNo))
                    assert(0);

            assert(file && this->_seqNo < file->numCollapsedVals());

        }


        BasicSeq *findSeqInFile(std::string &id, bool &isSubSeq) {

            TClass1 *fileSeq = (TClass1 *)file->findSeq(id);
            isSubSeq = false;

            if(fileSeq) {
                //	file->print((BasicSeq *)fileSeq);
                return (*fileSeq)[_seqNo];
            }

            /*
             * Check if it is a subSequence
             */
            boost::RegEx rExscont("^(\\S+)_(\\d+)_(\\d+)$");
            isSubSeq = rExscont.Match(id);

            if(!isSubSeq)  // unknown annotSeq object has been requested
                return NULL;

            /*
             * It's a subSequence object, subSequence objects have an id which is
             * SequenceId_beg_end, where beg and end mark the coordinates in the
             * real SequenceId object.
             */
            std::string superCId(rExscont[1]);
            int beg, end;

            if(!sscanf(rExscont[2].c_str(), "%d", &beg))
                assert(0);
            if(!sscanf(rExscont[3].c_str(), "%d", &end))
                assert(0);

            TClass1 *superC = (TClass1 *)file->findSeq(superCId);

            if(!superC) {
                isSubSeq = false;
                return NULL;
            }

            //      file->print((BasicSeq *)superC);
            return (*superC)[_seqNo]->getSubSequence(beg, end);

        }

        void freeValArrays(TStrand strand) {
            TypedFilter<TClass2>::freeValArrays(strand);
            //      file->uncacheSeqs(); // uncachingSeq occurs within the IOSource::cacheSeq routine now
        }

        /**
         * A member function for computing filter values.
         *
         * The value at each position of c is whatever value the private member
         * filter computes at said position.
         *
         * @param c a pointer to a Sequence object which has been read from the
         * fasta file.
         * @param strand the strand to use for accessing the c's string sequence.
         */
        virtual void computeVals(Sequence *c, TStrand strand) = 0;

        virtual ~FL_InpFile() {

        }
    };


    /**
     * FL_FastaxFilter is a subclass of FL_InpFile which wraps a fasta file
     * so that any filter could be computed on the fasta's contained sequences.
     *
     * Filter objects can only compute their values on the default input
     * sequence. The FL_FastaxFilter class makes it possible for any filter to
     * compute its values on arbitrary input sequences, as long as these
     * sequences are in fasta format.
     * Hence, this class behaves as a wrapper, so that any filter could be
     * computed over alternative input sequences stored in fasta format.
     *
     ***************************************************************************/

    template <class TClass1, class TClass2, class TClass3, class TClass4>
        class FL_FastaxFilter : public FL_InpFile<TClass1, TClass3> {
    protected:
        TypedFilter<TClass3> *filter;

    public:
        //! Default constructor
    FL_FastaxFilter(
        int fInd,          //!< A unique identifier.
        std::string name, //!< A unique name.
        IOSource *file,      //!< A Fasta file
        /*! A filter which is to be computed on the input
         * sequences contained in fasta
         */
        TypedFilter<TClass3> *filter,
        int seqNo          //!< sequence number inside file
        ) :
        FL_InpFile<TClass1, TClass3>(fInd,
                                     name,
                                     file,
                                     filter->maxNumFilterValues(),
                                     seqNo,
                                     filter->valType()
            ) {

            assert(this->file->type() != LDSCORE || this->file->type() != EXTENDED_LDSCORE);

            this->filter = filter;
            this->_period = filter->period();

        }

        /**
         * Constructor from Header string definition
         */
    FL_FastaxFilter(
        int fInd,       //!< A unique identifier.
        /*! The Header string definition, loaded as a vector of
         * strings.
         * The Header has the following form:\n\n
         * Filter name FastaxFilter<type> arrInd file filter\n\n
         * The meaning of these fields are as in the default
         * constructor
         * @see FL_FastaxFilter()
         */
        vector<std::string> & params,
        int & offset,       //!< The index for vector params.
        ResourceEngine *re, /*!< A pointer to the ResourceEngine
                              object. */
        FilterEngine *fe    /*!< A pointer to the FilterEngine
                              object. */
        )
        : FL_InpFile<TClass1, TClass3>(fInd,
                                       params,
                                       offset,
                                       0,
                                       re) {

            assert(this->file->type() != LDSCORE || this->file->type() != EXTENDED_LDSCORE);
            this->filter = (TypedFilter<TClass3> *)fe->getFilter(params[offset++]);
            assert(filter);

            this->maxFilterVals = filter->maxNumFilterValues();
            assert(this->type == filter->valType());

            this->_period = filter->period();
        }

        /**
         * A member function to set the default strand for the filter values.
         * @param strand the default strand
         */
        inline void setDefaultStrand(TStrand strand) {
            this->defaultStrand = strand;
            filter->setDefaultStrand(strand);
        }

        /**
         * A member function for computing filter values.
         *
         * The value at each position of c is whatever value the private member
         * filter computes at said position.
         *
         * @param c a pointer to a Sequence object which has been read from the
         * fasta file.
         * @param strand the strand to use for accessing the c's string sequence.
         * @see FL_InpFile::computeVals(Sequence *, TStrand)
         */
        inline void computeVals(Sequence *c, TStrand strand) {
            if(this->arrayInd() < 0 && !this->isSparse())
                return;

            if(!this->file->loaded())
                return;

            this->setDefaultStrand(strand);

            int first_seqNo = this->_seqNo, last_seqNo = this->_seqNo + 1,
                orig_seqNo = this->_seqNo, fileLength = 0;

            if(this->_seqNo <  0) {
                first_seqNo = 0;
                last_seqNo = this->file->numCollapsedVals();
            }

            for(this->_seqNo = first_seqNo; this->_seqNo < last_seqNo; this->_seqNo++) {
                bool isSubSeq = false;

                TClass2 *fileSeq = (TClass2 *)this->findSeqInFile(c->id(), isSubSeq);

                if(!fileSeq) {
                    assert(0);
                    throw EXCEPTION(BAD_USAGE, c->id() + string(" not found in file ") + this->getName());
                }

                int zero = (fileSeq->isZeroBased() ? 0 : 1);
                TClass4 *clone = (TClass4 *)fileSeq->cloneSeq(strand);
                fileLength = fileSeq->length(strand);
                //	cerr << this->getName() << "_" << this->_seqNo << " " << c->id() << " " << strand << endl;
                this->filter->computeVals((char *)((TClass4 *)(clone + zero)),
                                          this->filterVals[strand],
                                          fileLength);

                fileSeq->deleteClone(clone);

                if(isSubSeq)
                    delete fileSeq;

            }

            this->_seqNo = orig_seqNo;

            //      assert(fileLength);

            /*
             * Set the contained filter values variable to point to this
             * object's filter values.
             */
            this->filter->setValues(this->filterVals[strand], strand);

            if(this->isSparse())
                return;

            if(this->isLogScaled())
                this->filter->logScaleVals(this->filterVals[strand], fileLength);

            if(this->power() != 1)
                this->filter->raiseVals2Pow(this->filterVals[strand], fileLength,
                                            this->power());

            if(this->period())
                this->filter->accumVals(this->filterVals[strand],
                                        this->accFilterVals[strand],
                                        fileLength);

            this->filter->setAccValues(this->accFilterVals[strand], strand);

        }

        ~FL_FastaxFilter() {

        }
    };
}

#endif
