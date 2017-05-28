/****************************************************************************
 * Filter.h - part of the lless namespace, a general purpose
 *            linear semi-markov structure prediction library
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

#ifndef _FILTER_H_
#define _FILTER_H_
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Utils.h"
#include "Sequence.h"
#include "EvdEdges.h"
#include "ChangePtUtils.h"
#include <map>

#define FILTER_PADDING 500

namespace lless {

    /**
     *
     * Filter is property of the input sequence. In particular it can be seen
     * as a function that defines a value (of any type) at each position
     * of the input sequence. All filters by default, only need to be computed
     * prior to sequence decoding and freed afterwards; in that way, memory
     * usage is kept in check. Input sequences could be viewed in either forward
     * way (forward strand) or in reverse complemented way (complementary strand)
     * so potentially there are two possible values associated with a single
     * position of a sequence, depending on the strand which is used.
     *
     ***************************************************************************/

    class Filter {

    protected:
        int fInd;
        int arrInd;
        std::string name;
        int maxFilterVals;
        TValType type;
        TStrand defaultStrand;
        bool sparse;
        bool logScaled;
        double exponent;
        int _period;
    public:

        Filter(int fInd,
               std::string & name,
               int maxFilterVals = 1
            );

        Filter(int fInd,
               vector<std::string> & params,
               int & offset,
               int maxFilterVals = 1
            );

        inline int ind() {
            return fInd;
        }

        /**
         * @return the filter values' array index.
         */
        inline int arrayInd() {
            return arrInd;
        }

        /**
         * @return the filter's unique name.
         */
        inline std::string & getName() {
            return name;
        }

        /**
         * @return the filter values type
         */
        inline TValType valType() {
            return type;
        }

        inline void makeSparse() {
            this->sparse = true;
        }

        inline bool isSparse() {
            return this->sparse;
        }

        inline void makeLogScaled() {
            logScaled = true;
        }

        inline int isLogScaled() {
            return logScaled;
        }

        inline void setPower(double power) {
            exponent = power;
        }

        inline double power() {
            return exponent;
        }

        inline void setAccPeriod(int period) {
            this->_period = period;
        }

        /**
         * A member function that return the filter's period. This will be
         * different from zero if the filter can accumulate
         */
        virtual int period() {
            return _period;
        }

        virtual inline bool setArrayIndex(int arrayInd) {
            return this->arrInd = arrayInd;
        }

        /**
         * A member function that sets the filter's maximum number of different
         * filter values.
         * @param maxfilterVals the maximum number of different filter values
         */
        inline void setMaxNumFilterValues(int maxFilterVals) {
            this->maxFilterVals = maxFilterVals;
        }

        /**
         * @return the maximum number of different filter values
         */
        virtual inline int maxNumFilterValues() {
            return maxFilterVals;
        }

        /**
         * A member function to set the filter's default strand
         */
        virtual inline void setDefaultStrand(TStrand strand) {
            this->defaultStrand = strand;
        }

        virtual void allocValArrays(Sequence *c, TStrand strand) = 0;
        virtual void setPaddedValArrays(void ***filterVals) = 0;
        virtual void computeVals(Sequence *c, TStrand strand) = 0;
        virtual void freeValArrays(TStrand strand) = 0;

        virtual ~Filter() { }

    };

    template <class TClass>
        void logScaleFilterVals(TClass *vals, int len) {
        for(int i = 1; i <= len; i++)
            vals[i] = vals[i] > 1 ? log(vals[i]) : (TClass)0.001;
    }


    template <class TClass>
        void raiseVals2Pow(TClass *vals, int len, double power) {
        for(int i = 1; i <= len; i++)
            vals[i] = pow(vals[i], power);
    }

    template <class TClass>
        void accumFilterVals(TClass *vals, TClass **accumVals,
                             int len, int period, int maxVals) {

        assert(period);

        int i = 1;

        if(maxVals > 1) { // filter values can be used as indexes
            TClass fVal = vals[i - period];
            accumVals[i - period][(int)fVal]++;
            for( ; i <= len + period; i++) {
                fVal = vals[i];
                for(int j = 0; j < maxVals; j++)
                    accumVals[j][i] = accumVals[j][i - period] + (fVal == j);
            }
        }
        else {
            for(int i = 0; i <= len + period; i++)
                vals[i] += vals[i - period];
            accumVals[0] = vals;
        }
    }

    /**
     * A TypedFilter is a subclass of Filter that adds parametric type TClass
     * to filter values
     *
     ***************************************************************************/

    template <class TClass> class TypedFilter : public Filter {
    protected:
        TClass **filterVals;
        TClass ***paddedFilterVals;
        TClass ***paddedAccFilterVals;
        TClass ***accFilterVals;
    public:
    TypedFilter(int fInd,            //!< A unique identifier.
                std::string & name,  //!< A unique name.
                int maxFilterVals = 1,   /*!< The number of possible filter values.
                                           It is one if filter values are not
                                           countable. */

                TValType type = FT_INTEGER  //!< The value type
        )
        : Filter(fInd,
                 name,
                 maxFilterVals) {

            this->type = type;
            this->_period = 0;
            this->logScaled = false;
            this->exponent = 1;
            initFilterVals();

        }

        /**
         * Constructor from a Header string definition. Here the Header does not
         * contain type information and therefore the type needs to be passed as
         * an additional parameter
         */

    TypedFilter(int fInd,           //!< A unique filter identifier
                /*! The Header string definition, loaded as a vector
                  of strings. */
                vector<std::string> & params,
                int & offset,       //!< The index for vector params.
                int maxFilterVals,  /*!< The number of possible filter values.
                                      It is one if filter values are not
                                      countable. */
                TValType type       //!< The value type
        )
        : Filter(fInd,
                 params,
                 offset,
                 maxFilterVals) {

            this->type = type;
            offset++; // uniqueClsId
            this->_period = 0;
            this->logScaled = false;
            this->exponent = 1;

            initFilterVals();

        }

        /**
         * Constructor from a Header string definition.
         */
    TypedFilter(int fInd,           //!< A unique filter identifier
                /*! The Header string definition, loaded as a vector
                  of strings. */
                vector<std::string> & params,
                int & offset,       //!< The index for vector params.
                int maxFilterVals = 1 /*!< The number of possible filter values.
                                        It is one if filter values are not
                                        countable. */
        )
        : Filter(fInd,
                 params,
                 offset,
                 maxFilterVals) {

            /*
             * There is need to determine type from fargs[offset]
             */
            std::string & uClsId = params[offset++];
            type = Utils::extractValTypeFromClassName(uClsId);
            this->_period = 0;
            this->logScaled = false;
            this->exponent = 1;
            initFilterVals();

        }

        /**
         * A member function which initializes (to NULL) the array(s) containing
         * the filter values
         */
        void initFilterVals() {
            paddedFilterVals = NULL;

            filterVals = new TClass* [NUM_STRANDS];
            paddedAccFilterVals = new TClass** [NUM_STRANDS];
            accFilterVals = new TClass** [NUM_STRANDS];

            for(int j = 0; j < NUM_STRANDS; j++) {
                filterVals[j] = NULL;
                paddedAccFilterVals[j] = NULL;
                accFilterVals[j] = NULL;
            }
        }

        /**
         * A member function to set the filter values array
         * @param vals the new filter values array.
         */
        inline void setPaddedValArrays(void ***vals) {
            this->paddedFilterVals = (TClass ***)vals;
        }

        inline TClass*** accValues() {
            return accFilterVals;
        }

        inline TClass* values(TStrand strand) {
            return this->filterVals[strand];
        }

        inline void setValues(TClass *vals, TStrand strand) {
            this->filterVals[strand] = vals;
        }

        inline void setAccValues(TClass **accVals, TStrand strand) {
            this->accFilterVals[strand] = accVals;
        }

        inline TClass & value(int pos, TStrand strand) {
            assert(filterVals[strand] != NULL);
            return filterVals[strand][pos];
        }


        inline TClass & accValue(int pos, int fVal, TStrand strand) {
            assert(accFilterVals[strand][fVal] != NULL);
            return accFilterVals[strand][fVal][pos];
        }

        /**
         * @return the filter value at position pos in the default strand
         * @param pos the position in the input sequence.
         */
        inline TClass & value(int pos) {
            assert(filterVals[defaultStrand] != NULL);
            return filterVals[defaultStrand][pos];
        }

        /**
         * A member function to allocate the array for storing filter values that
         * appear in strand in the input sequence.
         * @param c the input Sequence object.
         * @param strand the strand used.
         */
        void allocValArrays(Sequence *c, TStrand strand) {
            assert(this->paddedFilterVals[strand]);

            if(arrayInd() < 0)
                return;

            TClass *array = this->paddedFilterVals[strand][arrayInd()];

            if(!array) {
                array = new TClass [c->length() + FILTER_PADDING];
                assert(array);

                for(int i = 0; i < c->length() + FILTER_PADDING; i++)
                    array[i] = TClass();

                this->paddedFilterVals[strand][arrayInd()] = array;

                if(period()) {
                    this->paddedAccFilterVals[strand] = new TClass *[maxNumFilterValues()];
                    this->accFilterVals[strand] = new TClass *[maxNumFilterValues()];
                    this->paddedAccFilterVals[strand][0] = NULL;
                    this->accFilterVals[strand][0] = NULL;

                    if(maxNumFilterValues() > 1) {
                        for(int j = 0; j < maxNumFilterValues(); j++) {
                            TClass *accArray = new TClass [c->length() + FILTER_PADDING];
                            for(int k = 0; k < c->length() + FILTER_PADDING; k++)
                                accArray[k] = TClass();

                            this->paddedAccFilterVals[strand][j] = accArray;
                            this->accFilterVals[strand][j] = accArray + FILTER_PADDING/2;
                        }
                    }
                }

                this->filterVals[strand] = array + FILTER_PADDING/2;
            }
        }

        /**
         * A member function to deallocate the array(s) that contains the filter
         * values in a given strand of the input sequence.
         * @param strand the strand
         */
        virtual void freeValArrays(TStrand strand) {

            if(this->arrayInd() < 0)
                return;

            if(this->paddedFilterVals[strand]
               && this->paddedFilterVals[strand][arrayInd()]) {

                delete [] this->paddedFilterVals[strand][arrayInd()];
            }

            this->paddedFilterVals[strand][arrayInd()] = NULL;
            this->filterVals[strand] = NULL;

            if(period()) {
                if(maxNumFilterValues() > 1) {
                    if(this->paddedAccFilterVals[strand])
                        for(int j = 0; j < maxNumFilterValues(); j++) {
                            if(this->paddedAccFilterVals[strand][j])
                                delete [] this->paddedAccFilterVals[strand][j];
                            this->paddedAccFilterVals[strand][j] = NULL;
                            this->accFilterVals[strand][j] = NULL;
                        }
                }
                else {
                    if(this->paddedAccFilterVals[strand])
                        this->paddedAccFilterVals[strand][0] = NULL;
                    if(this->accFilterVals[strand])
                        this->accFilterVals[strand][0] = NULL;
                }

                if(this->accFilterVals[strand])
                    delete [] this->accFilterVals[strand];
                this->accFilterVals[strand] = NULL;

                if(this->paddedAccFilterVals[strand])
                    delete [] this->paddedAccFilterVals[strand];
                this->paddedAccFilterVals[strand] = NULL;

            }
        }

        virtual inline void computeVals(char *seq, TClass *fVals, int len) {};

        /**
         * A member function that wraps function computeVals(seq, fVals, len)
         * as a function the receives a Sequence object and the strand which
         * is used to access the object's string sequence.
         * @param c a pointer to the input Sequence object.
         * @param strand the strand to use.
         */
        virtual inline void computeVals(Sequence *c, TStrand strand) {

            if(arrayInd() < 0 && !this->isSparse())
                return;

            //      cerr << c->id() << " " << this->getName() << " " << strand << endl;

            setDefaultStrand(strand);
            computeVals(c->getSeq(strand), filterVals[strand], c->length());

            if(this->isSparse())
                return;
            // post-processing
            if(isLogScaled())
                logScaleVals(filterVals[strand], c->length());

            if(this->power() != 1)
                raiseVals2Pow(this->filterVals[strand], c->length(), this->power());

            if(period())
                accumVals(filterVals[strand], accFilterVals[strand], c->length());

        }


        inline void accumVals(EdgeInst *vals, EdgeInst **accumVals, int len) {
            throw EXCEPTION(BAD_USAGE, "EdgeInst objects cannot accumulate");
        }

        inline void accumVals(ChangePoint *vals, ChangePoint **accumVals, int len) {
            throw EXCEPTION(BAD_USAGE, "ChangePoint objects cannot accumulate");
        }

        inline void accumVals(typename SPARSE_HASH<UCHAR, int> *vals,
                              typename SPARSE_HASH<UCHAR, int> **accumVals,
                              int len) {
            throw EXCEPTION(BAD_USAGE, "Hash objects cannot accumulate");
        }

        inline void accumVals(EvdEdges *vals, EvdEdges **accumVals, int len) {

            throw EXCEPTION(BAD_USAGE, "EvdEdges objects cannot accumulate");
        }

        inline void accumVals(vector<double> *vals,
                              vector<double> **accumVals, int len) {
            throw EXCEPTION(BAD_USAGE, "Vector objects cannot accumulate");
        }

        inline void accumVals(vector<int> *vals,
                              vector<int> **accumVals, int len) {
            throw EXCEPTION(BAD_USAGE, "Vector objects cannot accumulate");
        }

        inline void accumVals(UCHAR *vals, UCHAR **accumVals, int len) {
            accumFilterVals(vals, accumVals, len, period(), maxNumFilterValues());
        }

        inline void accumVals(USHORT *vals, USHORT **accumVals, int len) {
            accumFilterVals(vals, accumVals, len, period(), maxNumFilterValues());
        }

        inline void accumVals(int *vals, int **accumVals, int len) {
            accumFilterVals(vals, accumVals, len, period(), maxNumFilterValues());
        }

        inline void accumVals(ULONG *vals, ULONG **accumVals, int len) {
            accumFilterVals(vals, accumVals, len, period(), maxNumFilterValues());
        }

        inline void accumVals(double *vals, double **accumVals, int len) {
            assert(period());

            if(maxNumFilterValues() > 1)
                throw EXCEPTION(BAD_USAGE, "double values cannot be accumulated");

            for(int i = 0; i <= len + period(); i++)
                vals[i] += vals[i - period()];
            accumVals[0] = vals;

        }


        inline void logScaleVals(EdgeInst *vals, int len) {
            throw EXCEPTION(BAD_USAGE, "EdgeInst objects cannot be log scaled");
        }

        inline void logScaleVals(ChangePoint *vals, int len) {
            throw EXCEPTION(BAD_USAGE, "ChangePoint objects cannot be log scaled");
        }

        inline void logScaleVals(typename SPARSE_HASH<UCHAR, int> *vals,
                                 int len) {
            throw EXCEPTION(BAD_USAGE, "Hash objects cannot be log scaled");
        }

        inline void logScaleVals(EvdEdges *vals,
                                 int len) {
            for(int i = 1; i <= len; i++)
                vals[i].logScale();
        }

        inline void logScaleVals(vector<double> *vals,
                                 int len) {
            throw EXCEPTION(BAD_USAGE, "Vector objects cannot be log scaled");
        }

        inline void logScaleVals(vector<int> *vals,
                                 int len) {
            throw EXCEPTION(BAD_USAGE, "Vector objects cannot be log scaled");
        }

        inline void logScaleVals(UCHAR *vals, int len) {
            logScaleFilterVals(vals, len);
        }

        inline void logScaleVals(USHORT *vals, int len) {
            logScaleFilterVals(vals, len);
        }

        inline void logScaleVals(int *vals, int len) {
            logScaleFilterVals(vals, len);
        }

        inline void logScaleVals(ULONG *vals, int len) {
            logScaleFilterVals(vals, len);
        }

        inline void logScaleVals(double *vals, int len) {
            logScaleFilterVals(vals, len);
        }

        inline void raiseVals2Pow(EdgeInst *vals, int len, double power) {
            throw EXCEPTION(BAD_USAGE, "EdgeInst objects cannot be log scaled");
        }

        inline void raiseVals2Pow(ChangePoint *vals, int len, double power) {
            throw EXCEPTION(BAD_USAGE, "ChangePoint objects cannot be log scaled");
        }

        inline void raiseVals2Pow(typename SPARSE_HASH<UCHAR, int> *vals,
                                  int len, double power) {
            throw EXCEPTION(BAD_USAGE, "Hash objects cannot be log scaled");
        }

        inline void raiseVals2Pow(EvdEdges *vals,
                                  int len, double power) {
            for(int i = 1; i <= len; i++)
                vals[i].raise2pow(power);
        }

        inline void raiseVals2Pow(vector<double> *vals,
                                  int len, double power) {
            throw EXCEPTION(BAD_USAGE, "Vector objects cannot be log scaled");
        }

        inline void raiseVals2Pow(vector<int> *vals,
                                  int len, double power) {
            throw EXCEPTION(BAD_USAGE, "Vector objects cannot be log scaled");
        }

        inline void raiseVals2Pow(UCHAR *vals, int len, double power) {
            raiseVals2Pow(vals, len, power);
        }

        inline void raiseVals2Pow(USHORT *vals, int len, double power) {
            raiseVals2Pow(vals, len, power);
        }

        inline void raiseVals2Pow(int *vals, int len, double power) {
            raiseVals2Pow(vals, len, power);
        }

        inline void raiseVals2Pow(ULONG *vals, int len, double power) {
            raiseVals2Pow(vals, len, power);
        }

        inline void raiseVals2Pow(double *vals, int len, double power) {
            raiseVals2Pow(vals, len, power);
        }



        virtual ~TypedFilter() {

            delete [] filterVals;
            filterVals = NULL;
            delete [] accFilterVals;
            accFilterVals = NULL;
            delete [] paddedAccFilterVals;
            paddedAccFilterVals = NULL;
        }

    };

}

#endif
