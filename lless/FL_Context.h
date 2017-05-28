/****************************************************************************
 * FL_Context.h - part of the lless namespace, a general purpose
 *                  linear semi-markov structure prediction library
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

#ifndef _FILTER_CONTEXT_H_
#define _FILTER_CONTEXT_H_

#include "Filter.h"
#include <sstream>
#include "FilterEngine.h"

namespace lless {
    class ResourceEngine;

    /**
     * FL_Bin
     * is a subclass of TypedFilter<int> which computes a bin for the contained
     * filter at each positon along the input sequence.
     *
     **************************************************************************/

    template <class TClass> class FL_Bin : public TypedFilter<UCHAR> {
    protected:
        vector<TClass> knots;
        TypedFilter<TClass> *filter;
    public:
        //! Default constructor
    FL_Bin(int fInd,                //!< A unique identifier.
	   std::string & name,      //!< A unique name.
	   TypedFilter<TClass> *filter, //!< Contained filter.
	   vector<TClass> &knots             //!< bin control points.
        ) :
        TypedFilter<UCHAR>(fInd, name,
                           knots.size() - 1, FT_UCHAR) {

            this->filter = filter;
            this->knots = knots;

        }

        /**
         * Constructor from a Header string definition
         */
    FL_Bin(int fInd,          //!< A unique identifier.
	   /*! The Header string definition, loaded as a vector of strings.
	    * The Header has the following form:\n\n
	    * Filter name Bin filter numKnots [knot_1 .. ]\n\n
	    * numKnots is the number of knots used to define the bins.
	    * knot_i defines the i^th bin's boundaries i.e.
	    * if defined, bin i would range from knot_i+1 to
	    * knot_{i+1}
	    */
	   vector<std::string> &params,
	   int & offset,       //!< The index for vector params.
	   ResourceEngine *re, //!< A pointer to the ResourceEngine object.
	   FilterEngine *fe    //!< A pointer to the FilterEngine object.
        ) :
        TypedFilter<UCHAR>(fInd, params,
                           offset, 1, FT_UCHAR) {

            double knot;
            int num_knots;

            filter = (TypedFilter<TClass> *)fe->getFilter(params[offset++]);

            if(!sscanf(params[offset++].c_str(), "%d", &num_knots))
                assert(0);

            for(int i = 0; i < num_knots; i++) {
                if(!sscanf(params[offset++].c_str(), "%lf", &knot))
                    assert(0);

                knots.push_back((TClass)knot);
            }

            this->maxFilterVals = num_knots - 1;

        }

        /**
         * A member function for finding the bin number.
         * @param fVal the filter value.
         * @return The computed bin number.
         */
        UCHAR val2bin(TClass fVal) {
            int bin = 0;

            for(int j = 1; j <= this-> maxFilterVals; j++) {
                if(fVal > knots[j - 1] && fVal <= knots[j]) {
                    bin = j - 1;
                    break;
                }
            }

            return (UCHAR)bin;

        }

        inline vector<TClass> & ctrlPoints() {
            return knots;
        }

        /**
         * A member function for computing filter values.
         *
         * The value at each position of the input sequence is the context level
         * obtained through computation of the relative abundance of ctxSymbols
         * in a sliding window centered at said position.
         *
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        inline void computeVals(char *seq, UCHAR *filterVals, int len) {

            for(int i = -5; i <= len + 5; i++)
                filterVals[i] = val2bin(filter->value(i, defaultStrand));

        }

        ~FL_Bin() {
        }
    };


    /**
     * FL_Content
     * is a subclass of Filter which computes context levels at each sequence
     * position.
     *
     * The FL_Content class computes the context level at
     * each position in the sequence. A context is defined by a group of
     * characters which may appear in the input sequence and whose relative
     * abundance expressed in terms of levels, may define some key properties
     * of the sequence structure.
     *
     **************************************************************************/

    template<class TClass>
        class FL_Content: public TypedFilter<TClass> {
    private:
        int sldWindSize;
        bool ctxSymbols[128];

    public:
        //! Default constructor
    FL_Content(int fInd,                //!< A unique identifier.
               std::string & name,      //!< A unique name.
               int sldWindSize          /*!< A sliding window's size.
                                          If equal to zero, the window size
                                          is set automatically to the length
                                          of the current input sequence.
                                        */
        ) :
        TypedFilter<TClass>(fInd, name, 1) {

            this->sldWindSize = sldWindSize;

            for(int i = 0; i < 128; i++)
                ctxSymbols[i] = false;
        }

        /**
         * Constructor from a Header string definition
         */
    FL_Content(int fInd,           //!< A unique identifier.
               /*! The Header string definition, loaded as a vector of strings.
                * The Header has the following form:\n\n
                * Filter name Context arrInd ctxSymbols sldWindSize
                *
                */
               vector<std::string> &params,
               int & offset,       //!< The index for vector params.
               ResourceEngine *re, //!< A pointer to the ResourceEngine object.
               FilterEngine *fe    //!< A pointer to the FilterEngine object.
        ) :
        TypedFilter<TClass>(fInd, params, offset, 1) {

            for(int i = 0; i < 128; i++)
                ctxSymbols[i] = false;

            addCtxSymbols(params[offset++]);

            if(!sscanf(params[offset++].c_str(), "%d", &sldWindSize))
                assert(0);

        }

        inline void addCtxSymbols(std::string &symbols) {
            for(unsigned int i = 0; i < symbols.length(); i++)
                ctxSymbols[(int)symbols[i]] = true;
        }

        inline bool isCtxSymbol(char symbol) {
            return ctxSymbols[(int)symbol];
        }

        inline void setSlidingWindow(int window) {
            sldWindSize = window;
        }


        /**
         * A member function for computing filter values.
         *
         * The value at each position of the input sequence is the context level
         * obtained through computation of the relative abundance of ctxSymbols
         * in a sliding window centered at said position.
         *
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        inline void computeVals(char *seq, TClass *filterVals, int len) {
            assert(len);
            int i = 0;
            int ctxc = 0;
            int wSize = sldWindSize;

            if(!wSize || len < wSize)
                wSize = len;

            // 5 prime
            for(i = 0; i < wSize; i++)
                if(isCtxSymbol(seq[i]))
                    ctxc++;

            for(i = -5; i < wSize/2; i++)
                filterVals[i] = (UCHAR)(ctxc*100/wSize);

            // middle
            for( ; i < len - wSize/2; i++) {
                if(isCtxSymbol(seq[i-wSize/2]))
                    ctxc--;

                if(isCtxSymbol(seq[i+wSize/2]))
                    ctxc++;

                filterVals[i] = (TClass)(ctxc*100/wSize);
            }

            // 3 prime
            for( ; i <= len + 5; i++)
                filterVals[i] = (TClass)(ctxc*100/wSize);
        }

        ~FL_Content() {
        }
    };

    template<class TClass>
        class FL_FilterContent: public TypedFilter<TClass> {
    private:
        int sldWindSize;
        TypedFilter<TClass> *filter;

    public:
        //! Default constructor
    FL_FilterContent(int fInd,                //!< A unique identifier.
		     std::string & name,      //!< A unique name.
		     int sldWindSize,          /*!< A sliding window's size.
                                                 If equal to zero, the window size
                                                 is set automatically to the length
                                                 of the current input sequence.
					       */
                     TypedFilter<TClass> *filter
        ) :
        TypedFilter<TClass>(fInd, name, 1) {
            this->filter = filter;
            this->sldWindSize = sldWindSize;

        }

        /**
         * Constructor from a Header string definition
         */
    FL_FilterContent(int fInd,           //!< A unique identifier.
		     /*! The Header string definition, loaded as a vector of strings.
		      * The Header has the following form:\n\n
		      * Filter name Context arrInd ctxSymbols sldWindSize
		      *
		      */
		     vector<std::string> &params,
		     int & offset,       //!< The index for vector params.
		     ResourceEngine *re, //!< A pointer to the ResourceEngine object.
		     FilterEngine *fe    //!< A pointer to the FilterEngine object.
        ) :
        TypedFilter<TClass>(fInd, params, offset, 1) {

            if(!sscanf(params[offset++].c_str(), "%d", &sldWindSize))
                assert(0);
            this->filter = (TypedFilter<TClass> *)fe->getFilter(params[offset++]);
        }


        inline void setSlidingWindow(int window) {
            sldWindSize = window;
        }

        inline void setDefaultStrand(TStrand strand) {
            this->defaultStrand = strand;
            filter->setDefaultStrand(strand);

        }

        /**
         * A member function for computing filter values.
         *
         * The value at each position of the input sequence is the context level
         * obtained through computation of the relative abundance of ctxSymbols
         * in a sliding window centered at said position.
         *
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        inline void computeVals(char *seq, TClass *filterVals, int len) {
            assert(len);
            int i = 0;
            int wSize = sldWindSize;
            TClass acc_vals = TClass();
            if(!wSize || len < wSize)
                wSize = len;

            // 5 prime
            for(i = 0; i < wSize; i++)
                acc_vals += filter->value(i);

            for(i = -5; i < wSize/2; i++)
                filterVals[i] = acc_vals/wSize;

            // middle
            for( ; i < len - wSize/2; i++) {
                acc_vals -= filter->value(i-wSize/2);
                acc_vals += filter->value(i+wSize/2);
                filterVals[i] = acc_vals/wSize;
            }

            // 3 prime
            for( ; i <= len + 5; i++)
                filterVals[i] = acc_vals/wSize;
        }

        ~FL_FilterContent() {
        }
    };


    class FL_Context: public TypedFilter<UCHAR> {
    private:
        int sldWindSize;
        vector<int> _ctxLevRanges;
        bool ctxSymbols[128];

    public:
        //! Default constructor
    FL_Context(int fInd,                //!< A unique identifier.
               std::string & name,      //!< A unique name.
               vector<int> *ctxLevRanges, //!< Context level ranges.
               int sldWindSize          /*!< A sliding window's size.
                                          If equal to zero, the window size
                                          is set automatically to the length
                                          of the current input sequence.
                                        */
        ) :
        TypedFilter<UCHAR>(fInd, name, ctxLevRanges->size(), FT_UCHAR) {

            this->sldWindSize = sldWindSize;
            this->_ctxLevRanges = *ctxLevRanges;

            for(int i = 0; i < 128; i++)
                ctxSymbols[i] = false;
        }

        /**
         * Constructor from a Header string definition
         */
    FL_Context(int fInd,           //!< A unique identifier.
               /*! The Header string definition, loaded as a vector of strings.
                * The Header has the following form:\n\n
                * Filter name Context arrInd ctxSymbols sldWindSize
                * numCtxLevels [ctxLevel_1 .. ]\n\n
                * numCtxLevels is the number of levels that appear in this
                * context; ctxLevel_i defines the upper range for level i, i.e.
                * if defined, level i+1 would range from ctxLevel_i+1 to
                * ctxLevel_{i+1}
                *
                */
               vector<std::string> &params,
               int & offset,       //!< The index for vector params.
               ResourceEngine *re, //!< A pointer to the ResourceEngine object.
               FilterEngine *fe    //!< A pointer to the FilterEngine object.
        ) :
        TypedFilter<UCHAR>(fInd, params, offset, 1, FT_UCHAR) {

            int ctxClass;

            for(int i = 0; i < 128; i++)
                ctxSymbols[i] = false;

            addCtxSymbols(params[offset++]);

            if(!sscanf(params[offset++].c_str(), "%d", &sldWindSize))
                assert(0);

            if(!sscanf(params[offset++].c_str(), "%d", &maxFilterVals))
                assert(0);

            for(int i = 0; i < maxFilterVals; i++) {

                if(!sscanf(params[offset++].c_str(), "%d", &ctxClass))
                    assert(0);

                _ctxLevRanges.push_back(ctxClass);

            }
        }

        inline void addCtxSymbols(std::string &symbols) {
            for(unsigned int i = 0; i < symbols.length(); i++)
                ctxSymbols[(int)symbols[i]] = true;
        }

        inline bool isCtxSymbol(char symbol) {
            return ctxSymbols[(int)symbol];
        }

        inline void setSlidingWindow(int window) {
            sldWindSize = window;
        }

        /**
         * A member function for computing the context level.
         * @param ctxContent the relative content(abundance).
         * @return The computed context level.
         */
        UCHAR ctxLevel(double ctxContent) {
            int lastCTXClass = -1, ctxClassIndex = 0;

            for(unsigned int j = 0; j < _ctxLevRanges.size(); j++) {
                if(ctxContent >= lastCTXClass && ctxContent <= _ctxLevRanges[j]) {
                    ctxClassIndex = j;
                    break;
                }
                lastCTXClass = _ctxLevRanges[j];
            }

            return (UCHAR)ctxClassIndex;

        }

        inline vector<int> & ctxLevRanges() {
            return _ctxLevRanges;
        }

        /**
         * A member function for computing filter values.
         *
         * The value at each position of the input sequence is the context level
         * obtained through computation of the relative abundance of ctxSymbols
         * in a sliding window centered at said position.
         *
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        inline void computeVals(char *seq, UCHAR *filterVals, int len) {
            assert(len);
            int i = 0;
            double ctxc = 0;
            int wSize = sldWindSize;

            if(!wSize || len < wSize)
                wSize = len;

            // 5 prime
            for(i = 0; i < wSize; i++)
                if(isCtxSymbol(seq[i]))
                    ctxc++;

            int ctxl = ctxLevel(ctxc*100/wSize);

            for(i = -5; i < wSize/2; i++)
                filterVals[i] = ctxl;

            // middle
            for( ; i < len - wSize/2; i++) {
                if(isCtxSymbol(seq[i-wSize/2]))
                    ctxc--;

                if(isCtxSymbol(seq[i+wSize/2]))
                    ctxc++;

                ctxl = ctxLevel(ctxc*100/wSize);
                filterVals[i] = ctxl;
            }

            // 3 prime
            for( ; i <= len + 5; i++)
                filterVals[i] = ctxl;
        }

        ~FL_Context() {
        }
    };

}

#endif
