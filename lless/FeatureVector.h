/****************************************************************************
* FeatureVector.h - part of the lless namespace, a general purpose
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

#ifndef _FEATURE_VECTOR_H_
#define _FEATURE_VECTOR_H_
#include "VectorUtils.h"

namespace lless {

  /***************************************************************************
   * The Vector abstract base class defines the interface of a vector, which 
   * can be implemented as a sparse vector or dense vector. 
   ***************************************************************************/

  struct FeatVectorInd  {
    FeatVectorInd(int ind, int first, ULONG second) {
      this->ind = ind;
      this->first = first;
      this->second = second;
    }
    int ind;
    int first;
    ULONG second;
    ~FeatVectorInd() {}
  };

  template <class TClass> class Vector {
   protected:
    TClass _defaultValue;

   public:
    
    Vector(TClass dValue = TClass()) {
      this->_defaultValue = dValue;
    }

    virtual void reset(TClass = TClass()) = 0;

    inline TClass defaultValue() {
      return _defaultValue;
    }

    inline TClass defaultValue() const {
      return _defaultValue;
    }

    inline void setDefaultValue(TClass val) {
      _defaultValue = val;
    }

    virtual inline bool isSparse() {
      return false;
    }

    virtual inline bool isSparse() const {
      return false;
    }

    virtual TClass & operator[](ULONG i) = 0;
    virtual const TClass operator[](ULONG i) const = 0;

    virtual void maxValues(list<pair<FeatVectorInd, double> >&, 
			   int, int, int) const = 0;
    virtual double averageValue() const = 0;
    virtual Vector<TClass> & operator=(TClass) = 0; 
    virtual Vector<TClass> & product(const Vector<TClass> &fv) = 0;
    virtual Vector<TClass> & productInverse(const Vector<TClass> &fv) = 0;
    virtual Vector<TClass> & operator+=(const Vector<TClass> &fv) = 0;
    virtual Vector<TClass> & operator-=(const Vector<TClass> &fv) = 0;
    virtual Vector<TClass> & substractInverse(const Vector<TClass> &fv) = 0;
    virtual Vector<TClass> & highFilter(const Vector<TClass> &fv, TClass) = 0;
    virtual Vector<TClass> & operator=(const Vector<TClass> &fv) = 0;
    virtual double operator*(const Vector<TClass> &fv) const = 0;
    virtual Vector<TClass> & average(int, const Vector<TClass> &fv) = 0;

    virtual bool operator==(const Vector<TClass> &fv) = 0;
    virtual Vector<TClass> & operator/(double denom) = 0;
    virtual Vector<TClass> & operator*(double factor) = 0;

    virtual unsigned int size() = 0;
    virtual unsigned int size() const = 0;
    virtual unsigned int effectiveSize() = 0;

    virtual void print(std::ostream &ost, int offset = 0) const = 0;
    virtual void store(std::ofstream &ost, bool isBinary) = 0;
    virtual void retrieve(std::ifstream &ist, bool isBinary) = 0;

    virtual ~Vector() {}
  };
  

  template <class TClass1, class TClass2> 
    class TypedVector : public Vector<TClass2> {
   protected:
    TClass1 vals;
    ULONG _maxSize;

   public:
    TypedVector(ULONG maxSize, TClass2 dValue = TClass2()) 
      : Vector<TClass2>(dValue) { 
      this->_maxSize = maxSize;

    }
    
    TypedVector &operator=(const TypedVector<TClass1, TClass2> &fv) {
      vals = fv.values();
      this->_defaultValue = fv.defaultValue();
      _maxSize = fv.maxSize();
    }

    inline ULONG maxSize() {
      return _maxSize;
    }

    inline ULONG maxSize() const {
      return _maxSize;
    }

    inline TClass1 & values() {
      return vals;
    }

    inline const TClass1 & values() const {
      return vals;
    }
    
    inline typename TClass1::iterator begin() {
      return vals.begin();
    }

    inline typename TClass1::iterator end() {
      return vals.end();
    }

    inline typename TClass1::const_iterator begin() const {
      return vals.begin();
    }

    inline typename TClass1::const_iterator end() const {
      return vals.end();
    }
    

    inline void maxValues(list<pair<FeatVectorInd, double> > &maxV,
			  int ind, int i, int numMaxVals) const {
      
      typename TClass1::const_iterator cit = vals.begin();
      for( ; cit != vals.end(); cit++) {
	ULONG offset = position(cit);
	double v = value(cit);
	pair<FeatVectorInd, double> p(FeatVectorInd(ind, i, offset), v);
	
	if(!maxV.size()) {
	  maxV.push_back(p);
	  continue;
	}

	list<pair<FeatVectorInd, double> >::iterator it = maxV.begin();
	while(it != maxV.end() && p.second < it->second)
	  it++;
	maxV.insert(it, p);
	
	if(maxV.size() > numMaxVals)
	  maxV.pop_back();

      }
    }

    inline double averageValue() const {
      
      if(this->defaultValue())
	throw EXCEPTION(BAD_USAGE, "can't product NonZeroSparse TypedVector");

      double avg = 0;
      typename TClass1::const_iterator cit = vals.begin();
      for( ; cit != vals.end(); cit++)
	avg += value(cit);
      
      if(size()) avg /= size();
      
      return avg;
    }

    inline TypedVector<TClass1, TClass2> & 
      product(const Vector<TClass2> &fv) {

      if(this->defaultValue())
	throw EXCEPTION(BAD_USAGE, "can't product NonZeroSparse TypedVector");

      typename TClass1::iterator it = vals.begin();

      for( ; it != vals.end(); it++)
	value(it) *= fv[position(it)];
	
      return *this;
    }

    inline TypedVector<TClass1, TClass2> &
      productInverse(const Vector<TClass2> &fv) {

      if(this->defaultValue())
	throw EXCEPTION(BAD_USAGE, "can't productInverse NonZeroSparse TypedVector");

      typename TClass1::iterator it = vals.begin();
      
      for( ; it != vals.end(); it++)
	value(it) = value(it)*(1/fv[position(it)]);
	
      return *this;
    }

    
    inline TypedVector<TClass1, TClass2>  & operator+=(const Vector<TClass2> &fv) {
      if(fv.isSparse())
	VectorUtils::
	  add(*this, (const TypedVector<DENSE_HASH<ULONG, TClass2>, TClass2> &)fv);
      else
	VectorUtils::
	  add(*this, (const TypedVector<vector<TClass2>, TClass2> &)fv);
			 

      return *this;

    }

    inline TypedVector<TClass1, TClass2> & operator-=(const Vector<TClass2> &fv) {
      if(fv.isSparse())
	VectorUtils::
	  substract(*this, (const TypedVector<DENSE_HASH<ULONG, TClass2>, TClass2> &)fv);
      else
	VectorUtils::
	  substract(*this, (const TypedVector<vector<TClass2>, TClass2> &)fv);

      return *this;

    }

    inline TypedVector<TClass1, TClass2> & substractInverse(const Vector<TClass2> &fv) {
      if(fv.isSparse())
	VectorUtils::
	  substractInverse(*this, (const TypedVector<DENSE_HASH<ULONG, TClass2>, TClass2> &)fv);
      else
	VectorUtils::
	  substractInverse(*this, (const TypedVector<vector<TClass2>, TClass2> &)fv);
				      

      return *this;

    }

    inline TypedVector<TClass1, TClass2> & operator=(const Vector<TClass2> &fv) {
      if(fv.isSparse())
	VectorUtils::
	  assign(*this, (const TypedVector<DENSE_HASH<ULONG, TClass2>, TClass2> &)fv);
      else
	VectorUtils::
	  assign(*this, (const TypedVector<vector<TClass2>, TClass2> &)fv);

      return *this;

    }

    inline double operator*(const Vector<TClass2> &fv) const {
      if(fv.isSparse())
	return VectorUtils::
	  dot_product(*this, (const TypedVector<DENSE_HASH<ULONG, TClass2>, TClass2> &)fv);
      else
	return VectorUtils::
	  dot_product(*this, (const TypedVector<vector<TClass2>, TClass2> &)fv);
    }
    
    inline TypedVector<TClass1, TClass2> & average(int numIterations,
						   const Vector<TClass2> &fv) {
      if(fv.isSparse())
	VectorUtils::
	  average(numIterations, *this,
		  (const TypedVector<DENSE_HASH<ULONG, TClass2>, TClass2> &)fv);
      else
	VectorUtils::
	  average(numIterations, *this, 
		  (const TypedVector<vector<TClass2>, TClass2> &)fv);

      return *this;

    }

    inline TypedVector<TClass1, TClass2> & highFilter(const Vector<TClass2> &fv, TClass2 v) {      
      typename TClass1::iterator it = vals.begin();
      
      for( ; it != vals.end(); it++) {
	TClass2 tv = fv[position(it)];
	if(tv > v)
	  vals[position(it)] = 0;
      }

      return *this;
    }



    inline bool operator==(const Vector<TClass2> &fv) {      

      if(this->defaultValue())
	throw EXCEPTION(BAD_USAGE, "can't == NonZeroSparse TypedVector");

      typename TClass1::iterator it = vals.begin();

      for( ; it != vals.end(); it++)
	if(value(it) != fv[position(it)])
	  return false;

      return true && (fv.size() == this->size());
    }

    inline TypedVector<TClass1, TClass2> & operator/(double denom) {      

      if(this->defaultValue())
	throw EXCEPTION(BAD_USAGE, "can't == NonZeroSparse TypedVector");

      typename TClass1::iterator it = vals.begin();

      for( ; it != vals.end(); it++)
	value(it) /= denom;

      return *this;
    }

    inline TypedVector<TClass1, TClass2> & operator*(double factor) {

      if(this->defaultValue())
	throw EXCEPTION(BAD_USAGE, "can't == NonZeroSparse TypedVector");

      typename TClass1::iterator it = vals.begin();

      for( ; it != vals.end(); it++)
	value(it) *= factor;

      return *this;
    }

    inline void print(std::ostream &ost, int offset = 0) const {
      typename TClass1::const_iterator cit = vals.begin();

      for( ; cit != vals.end(); cit++) {
	TClass2 val = value(cit); 
	if(val)
	  ost << "\t\tbins[" << offset << "][" << position(cit) << "] = " << 
	    val << "\n";
      }
    }

    inline unsigned int size() {
      return vals.size();
    }

    inline unsigned int size() const {
      return vals.size();
    }

    virtual inline TClass2 & operator[](ULONG i) {
      return vals[i];
    }

    virtual const TClass2 operator[](ULONG i) const = 0;

    virtual ULONG position(typename TClass1::const_iterator &) const  = 0;
    virtual ULONG position(typename TClass1::iterator &) = 0;
    virtual TClass2 value(typename TClass1::const_iterator &) const = 0;
    virtual TClass2 & value(typename TClass1::iterator &) = 0;

    virtual ~TypedVector() {}
  };

  template <class TClass>
    class SparseVector : public TypedVector< DENSE_HASH<ULONG, TClass>, TClass> {
   protected:

   public:
    SparseVector(ULONG maxSize = 0, TClass defaultValue = TClass())
      : TypedVector<DENSE_HASH<ULONG, TClass>, TClass>(maxSize, defaultValue) { 

      this->vals.set_empty_key(ULONG_MAX);

    }

    void reset(TClass val = TClass()) {
      this->vals.clear();
      this->setDefaultValue(val);
    }

    inline bool isSparse() {
      return true;
    }

    inline bool isSparse() const {
      return true;
    }

    inline TClass & operator[](ULONG i) {      
      if(this->defaultValue()) {
	typename DENSE_HASH<ULONG, TClass>::iterator it = this->vals.find(i);
	if(it == this->vals.end()) 
	  this->vals[i] = this->defaultValue();
      }
      return this->vals[i];
    }

    inline const TClass operator[](ULONG i) const {
      typename DENSE_HASH<ULONG, TClass>::const_iterator cit = this->vals.find(i);
      return (cit == this->vals.end() ? this->defaultValue(): cit->second);
    }    

    ULONG position(typename DENSE_HASH<ULONG, TClass>::const_iterator & it) const {
      return it->first;
    }

    ULONG position(typename DENSE_HASH<ULONG, TClass>::iterator & it) {
      return it->first;
    }

    TClass value(typename DENSE_HASH<ULONG, TClass>::const_iterator & it) const {
      return it->second;
    }

    TClass & value(typename DENSE_HASH<ULONG, TClass>::iterator & it) {
      return it->second;
    }
    
    inline SparseVector<TClass> & operator=(TClass val) {
      this->vals.clear();      
      this->_defaultValue = val;
      
      return *this;
    }

    inline unsigned int effectiveSize() {
      return this->vals.size();
    }

    inline void store(std::ofstream &ost, bool isBinary) {
      
      unsigned int nz = effectiveSize();
      typename DENSE_HASH<ULONG, TClass>::const_iterator cit = this->vals.begin();
      TClass val;
      int stored = 0;
      
      ost << this->defaultValue() << " " << nz << endl;

      for( ; cit != this->vals.end(); cit++) {

	val = cit->second;
	ost << cit->first << "\t";
	
	if(isBinary)
	  ost.write((char *)&val, sizeof(val));
	else
	  ost << val;
	
	ost << endl;
	stored++;
      }
      
      if(stored != nz)
	throw EXCEPTION(BAD_USAGE, string("wrong #params in file"));
			
    }
    
    inline void retrieve(std::ifstream &ist, bool isBinary) {
      int nz;
      ULONG key;
      TClass val;
      ist >> this->_defaultValue >> nz;
      this->vals.clear();

      while(nz-- > 0) {
	ist >> key;
	if(isBinary) {
	  ist.get(); //tab
	  ist.read((char *)&val, sizeof(val));
	}
	else
	  ist >> val;
	
	ist.get(); //newline
	this->vals[key] = val;
	//	cerr << "\t" << key << "\t" << val << endl;
      }
    }

    ~SparseVector() {
      this->vals.clear();
    }

  };


  template <class TClass>
    class DenseVector : public TypedVector< vector<TClass>, TClass> {

   public:
    DenseVector(ULONG maxSize, TClass dValue = TClass()) 
      : TypedVector< vector<TClass>, TClass>(maxSize) {
      
      this->vals = vector<TClass>(maxSize, dValue);
    }
    
    void reset(TClass val = TClass()) {
      for(unsigned i = 0; i < this->vals.size(); i++) 
        this->vals[i] = val;
    }

    inline const TClass operator[](ULONG i) const {
      return this->vals[i];
    }    

    ULONG position(typename vector<TClass>::const_iterator & it) const {
      return it - this->values().begin();
    }

    ULONG position(typename vector<TClass>::iterator & it) {
      return it - this->values().begin();
    }

    TClass value(typename vector<TClass>::const_iterator & it) const {
      return (*it);
    }

    TClass & value(typename vector<TClass>::iterator & it) {
      return (*it);
    }

    inline DenseVector<TClass> & operator=(TClass val) {
      reset(val);
      return *this;
    }

    inline unsigned int effectiveSize() {
      unsigned nz = 0;
      for(unsigned i = 0; i < this->vals.size(); i++) 
        if(this->vals[i])
	  nz++;

      return nz;
    }

    void store(std::ofstream &ost, bool isBinary) {

      int stored = 0;
      TClass val;
      bool sparseStore = false;
      unsigned int nz = effectiveSize();

      if(nz < this->vals.size()/2) {
	ost << nz << "\t1" << endl;
	sparseStore = true;
      }
      else
	ost << this->vals.size() << "\t0" << endl;

      for(unsigned i = 0; i < this->vals.size(); i++) {	
	val = this->vals[i];
        if(sparseStore) {
          if(!val)
            continue;

	  ost << i << "\t";
        }

	if(isBinary)
	  ost.write((char *)&val, sizeof(val));
	else
	  ost << val;
	
	ost << endl;
	stored++;

      }

    }
    
    void retrieve(std::ifstream &ist, bool isBinary) {
      
      int stored = 0;
      double val;
      size_t key = -1;
      bool sparseStore;
      this->reset();

      ist >> stored >> sparseStore;
      ist.get(); //newline

      while(stored-- > 0) {
        if(sparseStore)
	  ist >> key;
	else 
	  key++;

	if(isBinary) {
	  if(sparseStore)
	    ist.get(); //tab

	  ist.read((char *)&val, sizeof(val));
	}
	else
	  ist >> val;
	
	ist.get(); //newline
	assert(key < this->vals.size());
	this->vals[key] = val;
	//	cerr << "\t" << key << "\t" << val << endl;
      }
    }

    ~DenseVector() {
      this->vals.clear();
    }

  };

  typedef Vector<double> FeatureVector;
  typedef SparseVector<double> SparseFeatureVector;
  typedef DenseVector<double> DenseFeatureVector;

}

#endif
