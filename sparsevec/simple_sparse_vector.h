//=============================================================================
// File Name: simple_sparse_vector.h
// Written by: Shai Shalev-Shwartz (23.03.06)
// implement a very simple sparse vector using stl map
//=============================================================================
#ifndef _SHAI_SIMPLE_SPARSE_VECTOR_
#define _SHAI_SIMPLE_SPARSE_VECTOR_

//*****************************************************************************
// Included Files
//*****************************************************************************
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <map>
#include <infra.h>



/** Implements a simple sparse vector based on stl map.
    The key is a string and the value is a double.
    @author Shai Shalev-Shwartz (shais@cs.huji.ac.il)
*/
class simple_sparse_vector {
 public:

  /** Default Constructor. Allocates an all zeros sparse vector
  */
  inline simple_sparse_vector() {  }

  /** Constructor. Read a line from a file of pairs of string and value.
      Construct a vector from this line.
      @param is The input stream to read from.
  */
  inline simple_sparse_vector(std::istream& is);
  
  /** Performs:  this += s*other
      @param other The other simple_sparse_vector
      @param s A double scalar
  */
  inline void add(simple_sparse_vector& other, double s = 1.0);

  /** Returns the squared l_2 norm of the vector
      @return the squared l_2 norm of the vector
  */
  inline double snorm();

  /** Convert the vector to a binary vector that just indicates 
      which elements are non-zero
  */
  inline void make_binary();


  std::map<std::string,double> my_vec;

};


//-----------------------------------------------------------------------------
/** Operator * for vector-vector multiplication
    @param u A reference to a simple_sparse_vector
    @param v A reference to a simple_sparse_vector
    @return The product (double)
*/
double operator* (simple_sparse_vector& u, simple_sparse_vector& v);


//-----------------------------------------------------------------------------
/** Operator * for "matrix"-vector multiplication.
    Note: we assume that the vector is much sparser than the rows of the matrix
    @param A A reference to a std vector of simple_sparse_vectors 
    @param v A reference to a simple_sparse_vector
    @return An infra::vector_base of the product
*/
infra::vector_base operator* (std::vector<simple_sparse_vector>& A, 
			      simple_sparse_vector& v);



#endif
