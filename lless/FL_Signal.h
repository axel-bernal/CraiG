/****************************************************************************
 * FL_Signal.h - part of the lless namespace, a general purpose
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

#ifndef _FILTER_SIGNAL_H_
#define _FILTER_SIGNAL_H_

#include "Filter.h"
#include "Consensus.h"
#include "ResourceEngine.h"
#include "FilterEngine.h"



namespace lless {

    /**
     *
     * FL_Signal class models different types of signals, i.e. short
     * stretches of DNA with non-random or irregular patterns. It does by using
     * an auxiliary class (XORconsensus) which allows the user to store many
     * different sequence patterns in a finite state automata. The FL_Signal
     * class will take the XOR of all of those stored patterns.
     * All FL_Signal-type objects are in OneStrand Format, as they need to be
     * heavily used during decoding (and decoding works more efficiently with
     * Tags in OneStrand format.
     * When wrapped up as a Tag object and used within a dotParamV or updParamV
     * feature functions, coordinates of the object need to change to TwoStrand
     * Format
     ***************************************************************************/

    template <class TClass> class FL_Signal: public TypedFilter<TClass> {
    protected:
        int sigType;
        XORConsensus *_consensus;

    public:
        //!< Default constructor
    FL_Signal(int fInd,              //!< A unique identifier.
              std::string name,    //!< A unique name.
              int sigType,           /*!< A signal type of the objects to
                                       be created as filter values. */
	      string signals,          /*!< A string representing the signal
                                         that needs to be found, as it appears
                                         in the sequence. */
              Sigma *sigma           //!< An alphabet for signals.
        )
        : TypedFilter<TClass>(fInd, name,
                              1, FT_SIGNAL) {

            this->sigType = sigType;
            boost::RegEx rExDelim("\\s+");
            std::vector<std::string> vsignals;
            rExDelim.Split(vsignals, signals, 0, 100);
            assert(vsignals.size());
            this->_consensus = new XORConsensus((char *)vsignals[0].c_str(), sigma);
            for(unsigned int j = 1; j < vsignals.size(); )
                _consensus->matchPattern((char *)vsignals[j++].c_str());

        }

        /**
         * Constructor from a Header string definition
         */
    FL_Signal(int fInd,           //!< A unique identifier.
              /*! The Header string definition, loaded as a vector of strings.
               * The Header has the following form:\n\n
               * Filter name  EdgeInst arrInd sigType sigma numSignals
               * [signal1 ..] \n\n
               * Instantiations of TClass = EdgeInst are the only ones to be
               * used with this constructor; sigType should be string which
               * represents a TEdgeId2 type. This string is internally
               * transformed into an integer.
               * The description of the other fields could be found in
               * the other constructor(s)
               */
              vector<std::string> &params,
              int & offset,       //!< The index for vector params.
              ResourceEngine *re, //!< A pointer to the ResourceEngine object.
              FilterEngine *fe    //!< A pointer to the FilterEngine object.
        )
        : TypedFilter<TClass>(fInd,
                              params,
                              offset,
                              1,
                              FT_SIGNAL) {

            sigType = TypeDefs::stringToTEdgeId2(params[offset++]);
            Sigma *sigma  = (Sigma *)re->getResource(params[offset++]);
            int numSignals;

            if(!sscanf(params[offset++].c_str(), "%d", &numSignals))
                assert(0);
            this->_consensus = new XORConsensus((char *)params[offset++].c_str(), sigma);
            for(int i = 1; i < numSignals; i++)
                _consensus->matchPattern((char *)params[offset++].c_str());
        }

        /**
         * A member function that matches (adds) a new patterns to the set of
         * patterns that currently match a signal occurrence.
         * @param pattern the pattern that identifies a signal occurrence
         */
        inline void matchSignal(char *pattern) {
            _consensus->matchPattern(pattern);
        }

        /**
         * A member function
         * @return the length of the _consensus of all matched signals (they all
         * have to be the same length
         */
        inline int len() {
            return _consensus->size();
        }

        inline int type() {
            return sigType;
        }

        inline XORConsensus *consensus() {
            return _consensus;
        };

        /**
         * A member function
         * @param seq the input sequence to be tested.
         * @return true if seq matches the signal, false otherwise.
         */
        inline virtual bool matches(char *seq) {
            return _consensus->matches(seq, len());
        }

        /**
         * A member function for computing filter values.
         *
         * The value at each position of the input sequence is an object of
         * class TClass which in most cases should be an EdgeInst (subclass of
         * a Tag object), but could eventually be used for potentially new objects
         * such as Regulators.
         *
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        inline void computeVals(char *seq, TClass *filterVals, int length) {
            int index = 1;
            int sigPos;
            int upperLimit = length - len() + 1;
            //      cerr << this->getName() << " " << this->defaultStrand << endl;
            while(index <= upperLimit) {
                sigPos = index;
                if(this->defaultStrand == STRAND_COMP)
                    sigPos = length - index + 1;
                if(matches(seq + index - 1)) {
                    filterVals[sigPos] = TClass(type(), sigPos, this->defaultStrand);
                    //	  cerr << sigPos << " " << type() << endl;
                }
                index++;
            }
        }

        virtual ~FL_Signal() {
            assert(_consensus);
            delete _consensus;
            _consensus = NULL;
        }

    };

    /**
     *
     * FL_Signal class models different types of signals, i.e. short
     * stretches of DNA with non-random or irregular patterns. It does by using
     * an auxiliary class (XORconsensus) which allows the user to store many
     * different sequence patterns in a finite state automata. The FL_Signal
     * class will take the XOR of all of those stored patterns.
     * All FL_Signal-type objects are in OneStrand Format, as they need to be
     * heavily used during decoding (and decoding works more efficiently with
     * Tags in OneStrand format.
     * When wrapped up as a Tag object and used within a dotParamV or updParamV
     * feature functions, coordinates of the object need to change to TwoStrand
     * Format
     ***************************************************************************/

    template <class TClass> class FL_FilterSignal: public FL_Signal<TClass> {
    protected:
        TypedFilter<UCHAR> *filter;
    public:
        //!< Default constructor
    FL_FilterSignal(int fInd,              //!< A unique identifier.
		    std::string & name,    //!< A unique name.
		    int sigType,           /*!< A signal type of the objects to
					     be created as filter values. */
		    char *signal,          /*!< A string representing the signal
					     that needs to be found, as it appears
					     in the sequence. */
		    Sigma *sigma,          //!< An alphabet for signals.
		    TypedFilter<UCHAR> *filter
        )
        : FL_Signal<TClass>(fInd, name, sigType,
                            signal, sigma) {

            this->filter = filter;
        }

        /**
         * Constructor from a Header string definition
         */
    FL_FilterSignal(int fInd,           //!< A unique identifier.
		    /*! The Header string definition, loaded as a vector of strings.
		     * The Header has the following form:\n\n
		     * Filter name  EdgeInst arrInd sigType sigma numSignals
		     * [signal1 ..] \n\n
		     * Instantiations of TClass = EdgeInst are the only ones to be
		     * used with this constructor; sigType should be string which
		     * represents a TEdgeId2 type. This string is internally
		     * transformed into an integer.
		     * The description of the other fields could be found in
		     * the other constructor(s)
		     */
		    vector<std::string> &params,
		    int & offset,       //!< The index for vector params.
		    ResourceEngine *re, //!< A pointer to the ResourceEngine object.
		    FilterEngine *fe    //!< A pointer to the FilterEngine object.
        )
        : FL_Signal<TClass>(fInd,
                            params,
                            offset,
                            re, fe) {

            this->filter = (TypedFilter<UCHAR> *)fe->getFilter(params[offset++]);

        }

        /**
         * A member function that extends the functionality of the original
         * function defined in Filter.h to change the contained filter object's
         * default strand to strand
         * @param strand the new default strand
         */
        inline void setDefaultStrand(TStrand strand) {
            this->defaultStrand = strand;
            filter->setDefaultStrand(strand);
        }

        /**
         * A member function for computing filter values.
         *
         * The value at each position of the input sequence is an object of
         * class TClass which in most cases should be an EdgeInst (subclass of
         * a Tag object), but could eventually be used for potentially new objects
         * such as Regulators.
         *
         * @param seq the input sequence.
         * @param filterVals the array of filter values.
         * @param len the length to consider. It could be shorter than the
         * actual length of the input sequence seq.
         */
        inline void computeVals(char *seq, TClass *filterVals, int length) {
            UCHAR *fseq = this->filter->values(this->defaultStrand) + 1;
            // need to add 1 since filters have a 1-based index
            FL_Signal<TClass>::computeVals((char *)fseq, filterVals, length);
        }

        ~FL_FilterSignal() { }

    };

}

#endif
