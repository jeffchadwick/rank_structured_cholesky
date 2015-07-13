// VECTOR.h: interface for the VECTOR class.
//
//////////////////////////////////////////////////////////////////////

#include "VECTOR.h"
#ifdef USE_MKL
#include <mkl_types.h>
#include <mkl_cblas.h>
#elif defined USING_ATLAS
extern "C" {
#include <cblas.h>
}
#else
// Use standard cblas
extern "C" {
#include <cblas.h>
}
#endif

#include <iostream>
#include <memory.h>

#include <rschol/linearalgebra/MATRIX.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor for the full vector
//////////////////////////////////////////////////////////////////////
VECTOR::VECTOR() :
  _size(0)
{
  _vector = NULL;
}

VECTOR::VECTOR(int size) :
  _size(size)
{
  _vector = new Real[_size];
  clear();
}

VECTOR::VECTOR(const char* filename)
{
	_vector = NULL;
  read(filename);
}

VECTOR::VECTOR(const VECTOR& v) 
{
  _size = v._size;
  _vector = new Real[_size];
  for (int x = 0; x < _size; x++)
    _vector[x] = v._vector[x];
}

VECTOR::VECTOR(int size, const Real *data)
{
  _size = size;
  _vector = new Real[_size];
  for (int x= 0; x < _size; x++) {
    _vector[x] = data[x];
  }
}

VECTOR::VECTOR(FILE* file)
{
  size_t bytes_read;

  // read dimensions
  bytes_read = fread((void*)&_size, sizeof(int), 1, file);

  // read data
  _vector = new Real[_size];
  if (sizeof(Real) == sizeof(float)) {
    double* vecDouble = new double[_size];
    bytes_read = fread((void*)vecDouble, sizeof(double), _size, file);
    for (int x = 0; x < _size; x++)
      _vector[x] = vecDouble[x];
    delete[] vecDouble;
  } else {
    bytes_read = fread((void*)_vector, sizeof(double), _size, file);
  }
}

VECTOR::~VECTOR()
{
  delete[] _vector;
}

//////////////////////////////////////////////////////////////////////
// dot product with another vector
//////////////////////////////////////////////////////////////////////
Real VECTOR::operator*(const VECTOR& vector)
{
#if BOUNDS_CHECKING_ENABLED
  if (vector._size != _size)
  {
    cout << __FILE__ << " " << __LINE__ << " VECTOR dot sizes do not match!: " << endl;
    return 123456.0f;
  }
#endif
#ifdef SINGLE_PRECISION
	return cblas_sdot(_size, _vector, 1, vector._vector, 1);
#else
	return cblas_ddot(_size, _vector, 1, vector._vector, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// dot product with another vector
//////////////////////////////////////////////////////////////////////
Real VECTOR::operator*(const Real *vector)
{
#ifdef SINGLE_PRECISION
	return cblas_sdot(_size, _vector, 1, vector, 1);
#else
	return cblas_ddot(_size, _vector, 1, vector, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// scale vector by a constant
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator*=(const Real& alpha)
{
#ifdef SINGLE_PRECISION
  cblas_sscal(_size, alpha, _vector, 1);
#else
  cblas_dscal(_size, alpha, _vector, 1);
#endif
  return *this;
}

//////////////////////////////////////////////////////////////////////
// wipe the whole vector
//////////////////////////////////////////////////////////////////////
void VECTOR::clear()
{
  memset(_vector, 0, _size * sizeof(Real));
}

//////////////////////////////////////////////////////////////////////
// resize and wipe the vector
//////////////////////////////////////////////////////////////////////
void VECTOR::resizeAndWipe(int size)
{
  if (_size != size)
    delete[] _vector;
  _size = size;

  _vector = new Real[_size];
  clear();
}

//////////////////////////////////////////////////////////////////////
// write the vector to a file
//////////////////////////////////////////////////////////////////////
void VECTOR::write(const char* filename)
{
  FILE* file;
  file = fopen(filename, "wb");

  if( file == NULL )
  {
	  printf( "** WARNING ** Could not write vector to %s\n", filename );
	  return;
  }

  // write dimensions
  fwrite((void*)&_size, sizeof(int), 1, file);

  // write data
  if (sizeof(Real) == sizeof(float))
  {
    double* vecDouble = new double[_size];
    for (int x = 0; x < _size; x++)
      vecDouble[x] = _vector[x];
    fwrite((void*)vecDouble, sizeof(double), _size, file);
    delete[] vecDouble;
  } 
  else
    fwrite((void*)_vector, sizeof(Real), _size, file);
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// write the vector to a file
//////////////////////////////////////////////////////////////////////
void VECTOR::write(FILE* file)
{
  // write dimensions
  fwrite((void*)&_size, sizeof(int), 1, file);

  // write data
  if (sizeof(Real) == sizeof(float))
  {
    double* vecDouble = new double[_size];
    for (int x = 0; x < _size; x++)
      vecDouble[x] = _vector[x];
    fwrite((void*)vecDouble, sizeof(double), _size, file);
    delete[] vecDouble;
  } 
  else
    fwrite((void*)_vector, sizeof(Real), _size, file);
}

//////////////////////////////////////////////////////////////////////
// read vector from a file
//////////////////////////////////////////////////////////////////////
bool VECTOR::read(const char* filename)
{
  FILE* file;
  size_t bytes_read;
  file = fopen(filename, "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : File " << filename << " not found! " << endl;
    return false;
  }

  // read dimensions
  bytes_read = fread((void*)&_size, sizeof(int), 1, file);

  // read data
  delete[] _vector;
  _vector = new Real[_size];

  if (sizeof(Real) == sizeof(float)) {
    double* vecDouble = new double[_size];
    bytes_read = fread((void*)vecDouble, sizeof(double), _size, file);
    for (int x = 0; x < _size; x++)
      _vector[x] = vecDouble[x];
    delete[] vecDouble;
  } else {
    bytes_read = fread((void*)_vector, sizeof(Real), _size, file);
  }
  fclose(file);

  return true;
}

//////////////////////////////////////////////////////////////////////
// Deep copy equality operator
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator=(VECTOR m)
{
  if (m.size() != _size)
  {
    delete[] _vector;
    _size= m.size();
    _vector= new Real[_size];
  }
	memcpy (_vector, m._vector, _size * sizeof(Real));
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Self-add operator
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator+=(const VECTOR& m)
{
#if BOUNDS_CHECKING_ENABLED
  if (m._size != _size)
    cout << __FILE__ << " " << __LINE__ << " : Vector sizes don't match! " << endl;
#endif
#ifdef SINGLE_PRECISION
	cblas_saxpy (_size, 1.0f, m._vector, 1, _vector, 1);
#else
	cblas_daxpy (_size, 1.0, m._vector, 1, _vector, 1);
#endif
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Self-subtract operator
//////////////////////////////////////////////////////////////////////
VECTOR& VECTOR::operator-=(const VECTOR& m)
{
#if BOUNDS_CHECKING_ENABLED
  if (m._size != _size)
  {
    cout << __FILE__ << " " << __LINE__ << " : Vector sizes don't match! " << endl;
    cout << " this: " << _size << endl;
    cout << " input: " << m._size << endl;
  }
#endif
#ifdef SINGLE_PRECISION
	cblas_saxpy (_size, -1.0f, m._vector, 1, _vector, 1);
#else
	cblas_daxpy (_size, -1.0, m._vector, 1, _vector, 1);
#endif
  return *this;
}

//////////////////////////////////////////////////////////////////////
// compute the 2 norm
//////////////////////////////////////////////////////////////////////
Real VECTOR::norm2()
{
#ifdef SINGLE_PRECISION
	return cblas_snrm2( _size, _vector, 1 );
#else
	return cblas_dnrm2( _size, _vector, 1 );
#endif
}

//////////////////////////////////////////////////////////////////////
// largest element in vector
//////////////////////////////////////////////////////////////////////
Real VECTOR::max_()
{
	Real vmax = _vector[0];

	for ( int i = 1; i < _size; i++ )
	{
		if ( _vector[i] > vmax )
			vmax = _vector[i];
	}

	return vmax;
}

//////////////////////////////////////////////////////////////////////
// smallest element in vector
//////////////////////////////////////////////////////////////////////
Real VECTOR::min_()
{
	Real vmin = _vector[0];

	for ( int i = 1; i < _size; i++ )
	{
		if ( _vector[i] < vmin )
			vmin = _vector[i];
	}

	return vmin;
}

//////////////////////////////////////////////////////////////////////
// largest absolute element in vector
//////////////////////////////////////////////////////////////////////
Real VECTOR::absmax()
{
	Real vmax = std::fabs( _vector[0] );

	for ( int i = 1; i < _size; i++ )
	{
		if ( std::fabs( _vector[i] ) > vmax )
			vmax = std::fabs( _vector[i] );
	}

	return vmax;
}

//////////////////////////////////////////////////////////////////////
// smallest absolute element in vector
//////////////////////////////////////////////////////////////////////
Real VECTOR::absmin()
{
	Real vmin = std::fabs( _vector[0] );

	for ( int i = 1; i < _size; i++ )
	{
		if ( std::fabs( _vector[i] ) < vmin )
			vmin = std::fabs( _vector[i] );
	}

	return vmin;
}

//////////////////////////////////////////////////////////////////////
// swap contents with another vector
//////////////////////////////////////////////////////////////////////
void VECTOR::swap(VECTOR& vec)
{
  Real* temp = _vector;
  _vector = vec._vector;
  vec._vector = temp;
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: y += alpha * x, where y is this vector
//////////////////////////////////////////////////////////////////////
void VECTOR::axpy(Real alpha, VECTOR& x)
{
#if BOUNDS_CHECKING_ENABLED
  if (_size != x._size)
    cout << __FILE__ << " " << __LINE__ << " : Vector size do not match! " << endl;
#endif
#ifdef SINGLE_PRECISION
	cblas_saxpy (_size, alpha, x._vector, 1, _vector, 1);
#else
	cblas_daxpy (_size, alpha, x._vector, 1, _vector, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// same as axpy above, but vector contents are stomped as well
//////////////////////////////////////////////////////////////////////
void VECTOR::clearingAxpy(Real alpha, VECTOR& x)
{
#if BOUNDS_CHECKING_ENABLED
  if (_size != x._size)
    cout << __FILE__ << " " << __LINE__ << " : Vector size do not match! " << endl;
#endif
  memset(_vector, 0, _size * sizeof(Real));
#ifdef SINGLE_PRECISION
	cblas_saxpy (_size, alpha, x._vector, 1, _vector, 1);
#else
	cblas_daxpy (_size, alpha, x._vector, 1, _vector, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// return the sum of the vector
//////////////////////////////////////////////////////////////////////
Real VECTOR::sum()
{
#ifdef SINGLE_PRECISION
  return cblas_sasum(_size, _vector, 1);
#else
  return cblas_dasum(_size, _vector, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// in-place copy, since operator= must allocate a new VECTOR
//////////////////////////////////////////////////////////////////////
void VECTOR::copyInplace(VECTOR& vector)
{
#if BOUNDS_CHECKING_ENABLED
  if (_size != vector._size)
    cout << __FILE__ << " " << __LINE__ << " : Vector sizes do not match! " << endl;
#endif
	memcpy (_vector, vector._vector, _size * sizeof(Real));
}

//////////////////////////////////////////////////////////////////////
// in-place copy, starting from a given index
//////////////////////////////////////////////////////////////////////
void VECTOR::copyInplace(VECTOR& vector, int start)
{
#if BOUNDS_CHECKING_ENABLED
  if ( (_size - start) < vector._size)
    cout << __FILE__ << " " << __LINE__ << " : Vector too small for copy! " << endl;
#endif
	memcpy (_vector + start, vector._vector, vector._size * sizeof(Real));
}

//////////////////////////////////////////////////////////////////////
// Copies a sub vector from the input vector in to the
// given location in this vector.
//////////////////////////////////////////////////////////////////////
void VECTOR::copySubVectorInplace(const VECTOR &vector,
																	int startLoc_in, int startLoc_out,
																	int size)
{
	if ( size <= 0 )
		return;

#if BOUNDS_CHECKING_ENABLED
	if ( startLoc_in + size > vector._size )
	{
    cout << __FILE__ << " " << __LINE__ << " : Index out of bounds in"
			"input vector! " << endl;
	}
	if ( startLoc_out + size > _size )
	{
    cout << __FILE__ << " " << __LINE__ << " : Index out of bounds in"
			"output vector! " << endl;
	}
#endif
	memcpy(_vector + startLoc_out, vector._vector + startLoc_in,
				 size * sizeof(Real));
}

//////////////////////////////////////////////////////////////////////
// add two vectors
//////////////////////////////////////////////////////////////////////
VECTOR operator+(VECTOR& x, VECTOR& y) 
{
  VECTOR z(x.size());

  for (int i = 0; i < x.size(); i++)
    z(i) = x(i) + y(i);
  return z;
}

//////////////////////////////////////////////////////////////////////
// subtract two vectors
//////////////////////////////////////////////////////////////////////
VECTOR operator-(VECTOR& x, VECTOR& y) 
{
  VECTOR z(x.size());

  for (int i = 0; i < x.size(); i++)
    z(i) = x(i) - y(i);
  return z;
}

//////////////////////////////////////////////////////////////////////
// scale a vector
//////////////////////////////////////////////////////////////////////
VECTOR operator*(VECTOR& x, Real& scalar) 
{
  VECTOR z(x.size());

  for (int i = 0; i < x.size(); i++)
    z(i) = x(i) * scalar;
  return z;
}

//////////////////////////////////////////////////////////////////////
// scale a vector
//////////////////////////////////////////////////////////////////////
VECTOR operator*(Real& scalar, VECTOR& x) 
{
  VECTOR z(x.size());

  for (int i = 0; i < x.size(); i++)
    z(i) = x(i) * scalar;
  return z;
}

//////////////////////////////////////////////////////////////////////
// compute the 2 norm
//////////////////////////////////////////////////////////////////////
Real operator^(VECTOR& x, VECTOR& y)
{
  Real total = 0;
  for (int i = 0; i < x.size(); i++)
    total += x.data()[i] * y.data()[i];
  return total;
}

//////////////////////////////////////////////////////////////////////
// max of the vector
//////////////////////////////////////////////////////////////////////
Real VECTOR::maxValue()
{
  if (_size == 0) return 0;
  
  Real maxFound = _vector[0];
  for (int x = 0; x < _size; x++)
    if (_vector[x] > maxFound) maxFound = _vector[x];

  return maxFound;
}

//////////////////////////////////////////////////////////////////////
// min of the vector
//////////////////////////////////////////////////////////////////////
Real VECTOR::minValue()
{
  if (_size == 0) return 0;
  
  Real minFound = _vector[0];
  for (int x = 0; x < _size; x++)
    if (_vector[x] < minFound) minFound = _vector[x];

  return minFound;
}

//////////////////////////////////////////////////////////////////////
// Take the absolute value of all entries
//////////////////////////////////////////////////////////////////////
void VECTOR::fabs()
{
  for (int x = 0; x < _size; x++)
    _vector[x] = std::fabs(_vector[x]);
}

//////////////////////////////////////////////////////////////////////
// compute the root mean squared
//////////////////////////////////////////////////////////////////////
Real VECTOR::rms()
{
  /*
  Real sumSq = 0.0;
  for (int x = 0; x < _size; x++)
    sumSq += _vector[x] * _vector[x];
  return sqrt(sumSq / _size);
  */
  return norm2() / sqrt((Real)_size);
}

//////////////////////////////////////////////////////////////////////
// compute the infinity norm
//////////////////////////////////////////////////////////////////////
Real VECTOR::normInf()
{
  if (_size == 0) return 0.0;
  Real max = std::fabs(_vector[0]);
  for (int x = 1; x < _size; x++)
    if (std::fabs(_vector[x]) > max)
      max = std::fabs(_vector[x]);
  return max;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
Real VECTOR::norm2( const Real *v, int size )
{
#ifdef DO_FLOP_COUNT
  MATRIX::MATRIX_FLOP_COUNT += 2 * size;
#endif

#ifdef SINGLE_PRECISION
	return cblas_snrm2( size, v, 1 );
#else
	return cblas_dnrm2( size, v, 1 );
#endif
}

//////////////////////////////////////////////////////////////////////
// Vector dot product
//////////////////////////////////////////////////////////////////////
Real VECTOR::dot( int n, const Real *x, const Real *y, int incX, int incY )
{
#ifdef DO_FLOP_COUNT
  // Only count the flops we are directly responsible for
  MATRIX::MATRIX_FLOP_COUNT += 2 * n;
#endif

#ifdef SINGLE_PRECISION
  return cblas_sdot( n, x, incX, y, incY );
#else
  return cblas_ddot( n, x, incX, y, incY );
#endif
}

//////////////////////////////////////////////////////////////////////
// Absolute sum of vector components
//////////////////////////////////////////////////////////////////////
Real VECTOR::absSum( int n, const Real *x, int incX )
{
#ifdef SINGLE_PRECISION
  return cblas_sasum( n, x, incX );
#else
  return cblas_dasum( n, x, incX );
#endif
}
