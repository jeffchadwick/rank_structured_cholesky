// VECTOR.h: interface for the VECTOR class.
//
//////////////////////////////////////////////////////////////////////

#ifndef VECTOR_H
#define VECTOR_H

#include <rschol/SETTINGS.h>
#include <map>
#include <iostream>
#include <cstdio>
#include <cmath>

using namespace std;

//////////////////////////////////////////////////////////////////////
// An arbitrary dimension vector class
//////////////////////////////////////////////////////////////////////
class VECTOR {

public:
  VECTOR();
  VECTOR(int size);
  VECTOR(const char* filename);
  VECTOR(int size, const Real *data);
  VECTOR(const VECTOR& v);
  VECTOR(FILE* file);
  ~VECTOR();

  inline Real& operator()(int index) { return _vector[index]; };
  inline Real operator()(int index) const { return _vector[index]; };
  inline Real& operator[](int index) { return _vector[index]; };
 
  int& size() { return _size; };
  int size() const { return _size; };

  // wipe the whole vector
  void clear();

  // write the vector to a binary file
  // everything is always written as a double
  void write(const char* filename);
  void write(FILE* file);

  // read the vector to a binary file
  // everything is always read in as a double, then
  // converted if necessary
  // Returns true if successfully read.
  bool read(const char* filename);

  // resize the vector and wipe to zero
  void resizeAndWipe(int size);

  // overloaded operators
  VECTOR& operator=(VECTOR m);
  VECTOR& operator+=(const VECTOR& v);
  VECTOR& operator-=(const VECTOR& v);
  VECTOR& operator*=(const Real& alpha);

  // (stevenan) The trailing underscore is to make this windows-compatible.
  // http://polingplace.blogspot.com/2007/04/stdmin-and-stdmax-versus-visual-studio.html
  // Thanks, Microsoft.
  Real max_();
  Real min_();

  Real absmax();
  Real absmin();

  // 2 norm
  Real norm2();

	// Infinity norm
	Real normInf();

	Real rms();

  // dot product
  Real operator*(const VECTOR& vector);
	Real operator*(const Real *vector);

  // swap contents with another vector --
  // it's your job to ensure that they are the same size
  void swap(VECTOR& vector);

  // raw data pointer
  Real* data() { return _vector; };
  const Real* data() const { return _vector; }

  // BLAS axpy operation: y += alpha * x, where y is this vector
  void axpy(Real alpha, VECTOR& x);

  // same as axpy above, but vector contents are stomped as well
  void clearingAxpy(Real alpha, VECTOR& x);

  // sum of all the elements
  Real sum();

  // in-place copy, since operator= must allocate a new VECTOR
  void copyInplace(VECTOR& vector);
	void equals(VECTOR& vector) { copyInplace(vector); }

	// In place copy from one vector to a specified location
	// in this vector.
	void copyInplace(VECTOR &vector, int startLoc);

	// Copies a sub vector from the input vector in to the
	// given location in this vector.
	void copySubVectorInplace(const VECTOR &vector,
														int startLoc_in, int startLoc_out,
														int size);

  // mean of all the entries
  Real mean() { return sum() / _size; };

  // max of all the elements
  Real maxValue();

  // min of all the elements
  Real minValue();

  // take the absolute value of all entires
  void fabs();

  static Real norm2( const Real *v, int size );

  // Vector dot product
  static Real dot( int n, const Real *x, const Real *y,
                   int incX = 1, int incY = 1 );

  // Absolute sum of vector components
  static Real absSum( int n, const Real *x, int incX = 1 );

private:
  int _size;
  Real* _vector;
};

//////////////////////////////////////////////////////////////////////
// dump vector to iostream
//////////////////////////////////////////////////////////////////////
inline ostream &operator<<(ostream &out, VECTOR& vector)
{
  out << "[" << endl; 
  for (int x = 0; x < vector.size(); x++)
    out << vector(x) << endl;
  out << "]" << endl;
  return out;
}

// overloaded operators
VECTOR operator-(VECTOR& x, VECTOR& y);
VECTOR operator+(VECTOR& x, VECTOR& y);
VECTOR operator*(VECTOR& x, Real& scalar);
VECTOR operator*(Real& scalar, VECTOR& x);

// x^T . y
Real operator^(VECTOR& x, VECTOR& y);

#endif
