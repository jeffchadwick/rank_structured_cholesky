// MATRIX.h: interface for the MATRIX class.
//
// Requires MKL version 10.3 or above, since we need mkl_lapacke.h,
// which stores C-like LAPACK routines.
//
//////////////////////////////////////////////////////////////////////

#include "MATRIX.h"

#ifdef USE_MKL
#include <mkl_lapack.h>
#include <mkl_lapacke.h>
#include <mkl_types.h>
#include <mkl_cblas.h>
#else
// Just use standard blas libraries
extern "C" {
#include <cblas.h>
#include <lapacke.h>
}
#endif

#include <util/IO.h>
#include <util/STLUtil.h>

#define TEST_ACCESS 1

size_t MATRIX::MATRIX_FLOP_COUNT = 0;

//////////////////////////////////////////////////////////////////////
// Constructor for the full matrix
//////////////////////////////////////////////////////////////////////
MATRIX::MATRIX() :
  _rows(0), _cols(0)
{
  _matrix = NULL;
}

MATRIX::MATRIX(int rows, int cols) :
  _rows(rows), _cols(cols)
{
  _matrix = new Real[_rows * _cols];
  clear();
}

MATRIX::MATRIX(int rows, int cols, const Real* data) :
  _rows(rows), _cols(cols)
{
  _matrix = new Real[_rows * _cols];
	memcpy(_matrix, data, _rows * _cols * sizeof(Real));
}

MATRIX::MATRIX(int rows, int cols, int lda, const Real *data) :
  _rows(rows), _cols(cols)
{
  _matrix = new Real[_rows * _cols];

  lda = max( lda, cols );

  for ( int i = 0; i < rows; i++ )
  for ( int j = 0; j < cols; j++ )
  {
    _matrix[ i * cols + j ] = data[ i * lda + j ];
  }
}

MATRIX::MATRIX(const char* filename)
{
	_matrix = NULL;
  read(filename);
}

MATRIX::MATRIX(const MATRIX& m)
{
  _cols = m._cols;
  _rows = m._rows;

  _matrix = new Real[_rows * _cols];
	memcpy(_matrix, m._matrix, _rows * _cols * sizeof(Real));
}

MATRIX::MATRIX(VECTOR& vec)
{
  _cols = vec.size();
  _rows = vec.size();

  _matrix = new Real[_rows * _cols];
  clear();

  for (int x = 0; x < vec.size(); x++)
    (*this)(x,x) = vec(x);
}

MATRIX::MATRIX(MATRIX3& matrix3)
{
  _cols = 3;
  _rows = 3;

  _matrix = new Real[9];
  clear();

  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++)
      (*this)(x,y) = matrix3(x,y);
}

MATRIX::~MATRIX()
{
  delete[] _matrix;
}

//////////////////////////////////////////////////////////////////////
// wipe the whole matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::clear()
{
  memset(_matrix, 0, _rows * _cols * sizeof(Real));
}

//////////////////////////////////////////////////////////////////////
// resize and wipe the matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::resizeAndWipe(int rows, int cols)
{
  if (_rows != rows || _cols != cols) {
    delete[] _matrix;

    _rows = rows;
    _cols = cols;

    _matrix = new Real[_rows * _cols];
  }
  clear();
}

//////////////////////////////////////////////////////////////////////
// write the matrix to a file
//////////////////////////////////////////////////////////////////////
void MATRIX::write(const char* filename)
{
  FILE* file;
  file = fopen(filename, "wb");

  // write dimensions
  fwrite((void*)&_rows, sizeof(int), 1, file);
  fwrite((void*)&_cols, sizeof(int), 1, file);

  // always write out as a double
  if (sizeof(Real) != sizeof(double))
  {
    double* matrixDouble = new double[_rows * _cols];
    for (int x = 0; x < _rows * _cols; x++)
      matrixDouble[x] = _matrix[x];

    fwrite((void*)matrixDouble, sizeof(double), _rows * _cols, file);
    delete[] matrixDouble;
    fclose(file);
  }
  else
    fwrite((void*)_matrix, sizeof(Real), _rows * _cols, file);
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// read matrix from a file
//////////////////////////////////////////////////////////////////////
void MATRIX::read(const char* filename)
{
  FILE* file;
  size_t bytes_read;
  file = fopen(filename, "rb");
  if (file == NULL)
  {
    cout << __FILE__ << " " << __LINE__ << " : File " << filename << " not found! " << endl;
    return;
  }

  // read dimensions
  bytes_read = fread((void*)&_rows, sizeof(int), 1, file);
  bytes_read = fread((void*)&_cols, sizeof(int), 1, file);

  // read data
  if( _matrix != NULL ) delete[] _matrix;

  _matrix = new Real[_rows * _cols];

  // always read in a double
  if (sizeof(Real) != sizeof(double)) {
    double* matrixDouble = new double[_rows * _cols];
    bytes_read = fread((void*)matrixDouble, sizeof(double), _rows * _cols, file);
    for (int x = 0; x < _rows * _cols; x++)
      _matrix[x] = matrixDouble[x];
    delete[] matrixDouble;
  } else {
    bytes_read = fread((void*)_matrix, sizeof(Real), _rows * _cols, file);
  }
  fclose(file);
}

//////////////////////////////////////////////////////////////////////
// return transpose of current matrix
//////////////////////////////////////////////////////////////////////
MATRIX MATRIX::transpose()
{
  MATRIX toReturn(_cols, _rows);

  for (int y = 0; y < _cols; y++)
    for (int x = 0; x < _rows; x++)
      toReturn(y,x) = (*this)(x,y);

  return toReturn;
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply
//////////////////////////////////////////////////////////////////////
VECTOR operator*(MATRIX& A, VECTOR& x) 
{
  VECTOR y(A.rows());

#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		CblasNoTrans, 
		A.rows(), A.cols(),
		1.0f, A.data(), A.cols(), x.data(), 1, 
		0.0f, y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		CblasNoTrans, 
		A.rows(), A.cols(),
		1.0, A.data(), A.cols(), x.data(), 1, 
		0.0, y.data(), 1);
#endif

  return y;
}

//////////////////////////////////////////////////////////////////////
// Scale matrix
//////////////////////////////////////////////////////////////////////
MATRIX operator*(MATRIX& A, Real alpha) 
{
  MATRIX y(A.rows(), A.cols());

  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < A.cols(); j++)
      y(i,j) = A(i, j) * alpha;

  return y;
}

//////////////////////////////////////////////////////////////////////
// scale matrix by a constant
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator*=(const Real& alpha)
{
#ifdef SINGLE_PRECISION
  cblas_sscal(_cols * _rows, alpha, _matrix, 1);
#else
  cblas_dscal(_cols * _rows, alpha, _matrix, 1);
#endif
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Matrix-Matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator*(MATRIX& A, MATRIX& B) 
{
  MATRIX y(A.rows(), B.cols());

#ifdef SINGLE_PRECISION
	cblas_sgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			A.rows(), 
      B.cols(),
      A.cols(), 
			1.0f,
			A.data(), A.cols(),
			B.data(), B.cols(),
			0.0f,
			y.data(), y.cols());
#else
	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			A.rows(), 
      B.cols(),
      A.cols(), 
			1.0,
			A.data(), A.cols(),
			B.data(), B.cols(),
			0.0,
			y.data(), y.cols());
#endif
  return y;
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply where A is transposed
//////////////////////////////////////////////////////////////////////
VECTOR operator^(MATRIX& A, VECTOR& x) 
{
  VECTOR y(A.cols());

#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		CblasTrans, 
		A.rows(), A.cols(), 
		1.0f, A.data(), A.cols(), x.data(), 1, 
		0.0f, y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		CblasTrans, 
		A.rows(), A.cols(), 
		1.0, A.data(), A.cols(), x.data(), 1, 
		0.0, y.data(), 1);
#endif
  return y;
}

//////////////////////////////////////////////////////////////////////
// Matrix^T -Matrix multiply
//////////////////////////////////////////////////////////////////////
MATRIX operator^(MATRIX& A, MATRIX& B) 
{
  MATRIX y(A.cols(), B.cols());

#ifdef SINGLE_PRECISION
	cblas_sgemm(
			CblasRowMajor,
			CblasTrans,
			CblasNoTrans,
			A.cols(), 
      B.cols(), 
      A.rows(),
			1.0f,
			A.data(), A.cols(),
			B.data(), B.cols(),
			0.0f,
			y.data(), y.cols());
#else
	cblas_dgemm(
			CblasRowMajor,
			CblasTrans,
			CblasNoTrans,
			A.cols(), 
      B.cols(), 
      A.rows(),
			1.0,
			A.data(), A.cols(),
			B.data(), B.cols(),
			0.0,
			y.data(), y.cols());
#endif
  return y;
}

//////////////////////////////////////////////////////////////////////
// Print matrix to stream
//////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream &out, MATRIX& matrix)
{
  out << "[" << endl; 
  for (int row = 0; row < matrix.rows(); row++)
  {
    for (int col = 0; col < matrix.cols(); col++)
      out << matrix(row, col) << " ";
    out << endl;
  }
  out << "]" << endl;
  return out;
}

//////////////////////////////////////////////////////////////////////
// Deep copy equality operator
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator=(const MATRIX m)
{
  if (m._cols != _cols || m._rows != _rows)
  {
    delete[] _matrix;
    _cols = m._cols;
    _rows = m._rows;

    _matrix = new Real[_rows * _cols];
  }
	memcpy(_matrix, m._matrix, _rows * _cols * sizeof(Real));
  return *this;
}

//////////////////////////////////////////////////////////////////////
// self minus
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator-=(const MATRIX& m)
{
  if (m._cols != _cols || m._rows != _rows)
  {
    delete[] _matrix;
    _cols = m._cols;
    _rows = m._rows;

    _matrix = new Real[_rows * _cols];
  }
#ifdef SINGLE_PRECISION
	cblas_saxpy(_cols * _rows, -1.0f, m._matrix, 1, _matrix, 1);
#else
	cblas_daxpy(_cols * _rows, -1.0, m._matrix, 1, _matrix, 1);
#endif
  return *this;
}

//////////////////////////////////////////////////////////////////////
// self plus
//////////////////////////////////////////////////////////////////////
MATRIX& MATRIX::operator+=(const MATRIX& m)
{
  if (m._cols != _cols || m._rows != _rows)
  {
    delete[] _matrix;
    _cols = m._cols;
    _rows = m._rows;

    _matrix = new Real[_rows * _cols];
  }
#ifdef SINGLE_PRECISION
	cblas_saxpy(_cols * _rows, 1.0f, m._matrix, 1, _matrix, 1);
#else
	cblas_daxpy(_cols * _rows, 1.0, m._matrix, 1, _matrix, 1);
#endif
  return *this;
}

//////////////////////////////////////////////////////////////////////
// Return the matrix diagonal
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::diagonal()
{
  int minDim = (_rows > _cols) ? _cols : _rows;
  VECTOR diag(minDim);
  for (int x = 0; x < minDim; x++)
    diag(x) = (*this)(x,x);

  return diag;
}

//////////////////////////////////////////////////////////////////////
// stomp the current matrix with the given matrix starting at "row". 
// It is your responsibility to ensure that you don't fall off the 
// end of this matrix.
//////////////////////////////////////////////////////////////////////
void MATRIX::setSubmatrix(MATRIX& matrix, int row)
{
  int totalSize = matrix.rows() * matrix.cols();
  int index = row * _cols;

  for (int x = 0; x < totalSize; x++, index++)
    _matrix[index] = matrix._matrix[x];
}

//////////////////////////////////////////////////////////////////////
// Stomp the current matrix with the given vector,
// starting at (row,col). So the first elemenet of the vector will
// be copied into (row,col), the next into (row+1, col), etc.
//////////////////////////////////////////////////////////////////////
void MATRIX::setVector( VECTOR& vector, int row, int col )
{
	for( int j = 0; j < vector.size(); j++ ) {
		(*this)(row+j, col) = vector(j);
	}
}

//////////////////////////////////////////////////////////////////////
// This assumes row-major storage
//////////////////////////////////////////////////////////////////////
void MATRIX::copyRowFrom( MATRIX& src, int srcRow, int row )
{
	memcpy( &(*this)(row,0), &src(srcRow,0), sizeof(Real)*_cols );
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: B += alpha * A, where B is this matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::axpy(Real alpha, MATRIX& A)
{
#ifdef SINGLE_PRECISION
	cblas_saxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#else
	cblas_daxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// BLAS axpy operation: B = alpha * A, where B is this matrix, and 
// current contents of B are stomped
//////////////////////////////////////////////////////////////////////
void MATRIX::clearingAxpy(Real alpha, MATRIX& A)
{
  memset(_matrix, 0, _rows * _cols * sizeof(Real));
#ifdef SINGLE_PRECISION
	cblas_saxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#else
	cblas_daxpy(_cols * _rows, alpha, A._matrix, 1, _matrix, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// BLAS gemm operation: C += alpha * A * B where C is this matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::gemm(Real alpha, MATRIX& A, MATRIX& B)
{
#ifdef SINGLE_PRECISION
	cblas_sgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			A.rows(), 
      B.cols(),
      A.cols(), 
			1.0f,
			A.data(), A.cols(),
			B.data(), B.cols(),
			0.0f,
			_matrix, _cols);
#else
	cblas_dgemm(
			CblasRowMajor,
			CblasNoTrans,
			CblasNoTrans,
			A.rows(), 
      B.cols(),
      A.cols(), 
			1.0,
			A.data(), A.cols(),
			B.data(), B.cols(),
			0.0,
			_matrix, _cols);
#endif
}

/*
 * Untested -- don't uncomment until a test case comes up
 * 
//////////////////////////////////////////////////////////////////////
// BLAS gemv operation: y += A * x where A is this matrix
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::gemv(VECTOR& x)
{
  if (x.size() != _rows)
    cout << __FILE__ << " " << __LINE__ << " : Matrix-Vector dimensions do not match!" << endl;

  VECTOR y(x.size());
  for (int j = 0; j < _cols; j++)
    for (int i = 0; i < _rows; i++)
      y(j) += (*this)(i,j) * x(j);

  return y;
}
*/

//////////////////////////////////////////////////////////////////////
// BLAS gemv operation: y += A * x where A is this matrix
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::gemv(VEC3F& x)
{
  VECTOR y(_rows);
#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		1.0f, _matrix, _cols, x, 1, 
		0.0f, y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		1.0, _matrix, _cols, x, 1, 
		0.0, y.data(), 1);
#endif
  return y;
}

//////////////////////////////////////////////////////////////////////
// BLAS gemv operation: y += alpha * A * x where A is this matrix
//////////////////////////////////////////////////////////////////////
VECTOR MATRIX::gemv(Real alpha, VEC3F& x)
{
  VECTOR y(_rows);
#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		alpha, _matrix, _cols, x, 1, 
		0.0f, y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		alpha, _matrix, _cols, x, 1, 
		0.0, y.data(), 1);
#endif
  return y;
}

//////////////////////////////////////////////////////////////////////
// BLAS gemm operation: C = alpha * A * B where C is this matrix and
// current contents of C are stomped
//////////////////////////////////////////////////////////////////////
void MATRIX::clearingGemm(Real alpha, MATRIX& A, MATRIX& B,
													bool transposeA, bool transposeB)
{
#ifdef SINGLE_PRECISION
	cblas_sgemm(
			CblasRowMajor,
			transposeA ? CblasTrans : CblasNoTrans,
			transposeB ? CblasTrans : CblasNoTrans,
			transposeA ? A.cols() : A.rows(), 
      transposeB ? B.rows() : B.cols(),
      transposeA ? A.rows() : A.cols(), 
			alpha,
			A.data(),
			A.cols(),
			//transposeA ? A.rows() : A.cols(),
			B.data(),
			B.cols(),
			//transposeB ? B.rows() : B.cols(),
			0.0f,
			_matrix,
			_cols); /// Should be rows if tranposed B?
			//transposeB ? _rows : _cols); /// Should be rows if tranposed B?
#else
	cblas_dgemm(
			CblasRowMajor,
			transposeA ? CblasTrans : CblasNoTrans,
			transposeB ? CblasTrans : CblasNoTrans,
			transposeA ? A.cols() : A.rows(), 
      transposeB ? B.rows() : B.cols(),
      transposeA ? A.rows() : A.cols(), 
			alpha,
			A.data(),
			A.cols(),
			//transposeA ? A.rows() : A.cols(),
			B.data(),
			B.cols(),
			//transposeB ? B.rows() : B.cols(),
			0.0,
			_matrix,
			_cols); /// Should be rows if tranposed B?
			//transposeB ? _rows : _cols); /// Should be rows if tranposed B?
#endif
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply
//
// Do *NOT* let x == y!
//////////////////////////////////////////////////////////////////////
void MATRIX::multiplyInplace(VECTOR& x, VECTOR& y, Real alpha, bool add) 
{
#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		//1.0f, _matrix, _cols, x.data(), 1, 
		alpha, _matrix, _cols, x.data(), 1, 
		add ? 1.0f : 0.0f,
		y.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		//1.0, _matrix, _cols, x.data(), 1, 
		alpha, _matrix, _cols, x.data(), 1, 
		add ? 1.0 : 0.0,
		y.data(), 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiply
//
// Do *NOT* let x == y!
//////////////////////////////////////////////////////////////////////
void MATRIX::multiplyInplace(Real *x, Real *y, Real alpha, bool add)
{
#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		//1.0f, _matrix, _cols, x, 1, 
		alpha, _matrix, _cols, x, 1, 
		add ? 1.0f : 0.0f,
		y, 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		CblasNoTrans,
    _rows, _cols,
		//1.0, _matrix, _cols, x, 1, 
		alpha, _matrix, _cols, x, 1, 
		add ? 1.0 : 0.0,
		y, 1);
#endif
}

void MATRIX::subMatrixMultiplyInplace( VECTOR& x, VECTOR& prod,
                                       int subRows, int subCols,
                                       bool transpose )
{
#ifdef SINGLE_PRECISION
	cblas_sgemv(
		CblasRowMajor, 
		transpose ? CblasTrans : CblasNoTrans,
		subRows, subCols,
		1.0f, _matrix, _cols, x.data(), 1, 
		0.0f, prod.data(), 1);
#else
	cblas_dgemv(
		CblasRowMajor, 
		transpose ? CblasTrans : CblasNoTrans,
		subRows, subCols,
		1.0, _matrix, _cols, x.data(), 1, 
		0.0, prod.data(), 1);
#endif
}

void MATRIX::uppertriMultiplyInplace( VECTOR& x, VECTOR& prod )
{
#ifdef SINGLE_PRECISION
	cblas_ssymv(
		CblasRowMajor, 
		CblasUpper, _rows,
		1.0, _matrix, _cols, x.data(), 1,
		0.0, prod.data(), 1 );
#else
	cblas_dsymv(
		CblasRowMajor, 
		CblasUpper, _rows,
		1.0, _matrix, _cols, x.data(), 1,
		0.0, prod.data(), 1 );
#endif
}

//////////////////////////////////////////////////////////////////////
// solve the linear system Ax = b, return x in the passed in b
// Here, A is a general invertible system.
//////////////////////////////////////////////////////////////////////
void MATRIX::solve(VECTOR& b)
{
  char uplo = 'U';
  int nrhs = 1;
  int info = 0;
  int R = _rows;
	//int ipiv[_rows];

	lapack_int ipiv[_rows];

  lapack_int result;

  // call the general LU solver from MKL
#ifdef SINGLE_PRECISION
  //sgesv(&R, &nrhs, _matrix, &R, ipiv, b.data(), &R, &info);
  result = LAPACKE_sgesv( LAPACK_ROW_MAJOR,
                          R, nrhs, _matrix, R, ipiv, b.data(), 1 );
#else
  //dgesv(&R, &nrhs, _matrix, &R, ipiv, b.data(), &R, &info);
  result = LAPACKE_dgesv( LAPACK_ROW_MAJOR,
                          R, nrhs, _matrix, R, ipiv, b.data(), 1 );
#endif
}

//////////////////////////////////////////////////////////////////////
// solve the linear system Ax = b, return x in the passed in b
// Here, A is a general invertible system.
//////////////////////////////////////////////////////////////////////
void MATRIX::solveSPD(VECTOR& b)
{
  char uplo = 'L';
  int nrhs = 1;
  int info = 0;
  int R = _rows;

  lapack_int result;

  // call the SPD solver from MKL
#ifdef SINGLE_PRECISION
  //sposv(&uplo, &R, &nrhs, _matrix, &R, b.data(), &R, &info);
  result = LAPACKE_sposv( LAPACK_ROW_MAJOR, uplo,
                          R, nrhs, _matrix, R, b.data(), 1 );
#else
  //dposv(&uplo, &R, &nrhs, _matrix, &R, b.data(), &R, &info);
  result = LAPACKE_dposv( LAPACK_ROW_MAJOR, uplo,
                          R, nrhs, _matrix, R, b.data(), 1 );
#endif
}

//////////////////////////////////////////////////////////////////////
// solve the linear system AX = B, return x in the passed in b
// Here, A is a general invertible system.
//
// Note: This version involves taking tranposes (eg. dynamic
// allocation).
//////////////////////////////////////////////////////////////////////
void MATRIX::solveSPD(MATRIX &B)
{
  char uplo = 'L';
  int nrhs = B.cols();
  int info = 0;
  int R = _rows;

#if 0
  // Create column major ordering of data
  Real *rhsData = new Real[ B.rows() * B.cols() ];

  int idx = 0;
  for ( int j = 0; j < B.cols(); j++ )
  {
    for ( int i = 0; i < B.rows(); i++ )
    {
      rhsData[ idx ] = B( i, j );
      idx++;
    }
  }
#endif

  lapack_int result;

  // call the SPD solver from MKL
#ifdef SINGLE_PRECISION
  //sposv(&uplo, &R, &nrhs, _matrix, &R, rhsData, &R, &info);
  result = LAPACKE_sposv( LAPACK_ROW_MAJOR, uplo,
                          R, nrhs, _matrix, R, B.data(), B.cols() );
#else
  //dposv(&uplo, &R, &nrhs, _matrix, &R, rhsData, &R, &info);
  result = LAPACKE_dposv( LAPACK_ROW_MAJOR, uplo,
                          R, nrhs, _matrix, R, B.data(), B.cols() );
#endif

  if ( result < 0 )
  {
    cout << "ERROR: dposv failed!  Invalid argument" << endl;
    abort();
  }
  else if ( result > 0 )
  {
    cout << "ERROR: dposv failed!  Matrix not positive definite" << endl;
    abort();
  }

#if 0
  if ( info < 0 )
  {
    cout << "ERROR: dposv failed!  Invalid argument" << endl;
    abort();
  }
  else if ( info > 0 )
  {
    cout << "ERROR: dposv failed!  Matrix not positive definite" << endl;
    abort();
  }
#endif

#if 0
  // Copy back in to B
  idx = 0;
  for ( int j = 0; j < B.cols(); j++ )
  {
    for ( int i = 0; i < B.rows(); i++ )
    {
      B( i, j ) = rhsData[ idx ];
      idx++;
    }
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Inverts this matrix in place
//////////////////////////////////////////////////////////////////////
void MATRIX::invert()
{
  cerr << "MATRIX::invert not implemented" << endl;
  abort();

#if 0
  // basic error checking
  if (_rows != _cols)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Matrix must be square to be inverted! " << endl;
    return;
  }

  int R = _rows;
	int ipiv[_rows];
	int info = 0;
	Real work_query;
	int l_work = -1;

#ifdef USE_MKL
#ifdef SINGLE_PRECISION
	sgetrf( &R, &R, _matrix, &R, ipiv, &info );

	// Query for workspace size
	sgetri( &R, _matrix, &R, ipiv, &work_query, &l_work, &info );

	l_work = (int)work_query;
	Real work[l_work];

	sgetri( &R, _matrix, &R, ipiv, work, &l_work, &info );
#else
	dgetrf( &R, &R, _matrix, &R, ipiv, &info );

	// Query for workspace size
	dgetri( &R, _matrix, &R, ipiv, &work_query, &l_work, &info );

	l_work = (int)work_query;
	Real work[l_work];

	dgetri( &R, _matrix, &R, ipiv, work, &l_work, &info );
#endif
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// Returns the inverse of this matrix in A
//////////////////////////////////////////////////////////////////////
void MATRIX::inverse( MATRIX &A )
{
  cerr << "MATRIX::inverse not implemented" << endl;
  abort();

#if 0
  // basic error checking
  if (_rows != _cols)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Matrix must be square to be inverted! " << endl;
    return;
  }

	A = *this;

  int R = _rows;
	int ipiv[_rows];
	int info = 0;
	Real work_query;
	int l_work = -1;

#ifdef USE_MKL
#ifdef SINGLE_PRECISION
	sgetrf( &R, &R, A._matrix, &R, ipiv, &info );

	// Query for workspace size
	sgetri( &R, A._matrix, &R, ipiv, &work_query, &l_work, &info );

	l_work = (int)work_query;
	Real work[l_work];

	sgetri( &R, A._matrix, &R, ipiv, work, &l_work, &info );
#else
	dgetrf( &R, &R, A._matrix, &R, ipiv, &info );

	// Query for workspace size
	dgetri( &R, A._matrix, &R, ipiv, &work_query, &l_work, &info );

	l_work = (int)work_query;
	Real work[l_work];

	dgetri( &R, A._matrix, &R, ipiv, work, &l_work, &info );
#endif
#endif
#endif
}

//////////////////////////////////////////////////////////////////////
// solve for the eigensystem of the matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem(VECTOR& eigenvalues, MATRIX& eigenvectors)
{
  cerr << "MATRIX::eigensystem not implemented" << endl;

#if 0
  // basic error checking
  if (_rows != _cols)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " Matrix must be square to get eigenvalues! " << endl;
    return;
  }

  // resize result space
  eigenvalues.resizeAndWipe(_rows);
  eigenvectors.resizeAndWipe(_rows, _rows);

  int rowsize = _rows;
  int worksize = 5 * _rows;
  Real* work = new Real[worksize];
  Real* valuesReal = new Real[2 * _rows];
  Real* valuesImag = valuesReal + _rows;
  Real* vectors = new Real[_rows * _rows];

  // the actual LAPACK call
  int error;
#ifdef SINGLE_PRECISION
  sgeev("V","N", &rowsize, _matrix, &rowsize, 
        valuesReal, valuesImag, 
        vectors, &rowsize, NULL, &rowsize,
        work, &worksize, &error);
#else
  dgeev("V","N", &rowsize, _matrix, &rowsize, 
        valuesReal, valuesImag, 
        vectors, &rowsize, NULL, &rowsize,
        work, &worksize, &error);
#endif

  if (error != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " eigenvalue solver bombed!" << endl;
  }

  // copy out results
  for (int x = 0; x < _rows; x++)
    eigenvalues(x) = valuesReal[x];

  for (int x = 0; x < _rows; x++)
    for (int y = 0; y < _rows; y++)
      eigenvectors(x,y) = vectors[x + y * _cols];
 
  // cleanup
  delete[] work;
  delete[] valuesReal;
  delete[] vectors;
#endif
}

//////////////////////////////////////////////////////////////////////
// Performs QR factorization of the given system.  Other calls
// will be needed following this to work the factorization,
// since this is just a wrapper for the basic LAPACK call.
//
// No pivoting is used here
//////////////////////////////////////////////////////////////////////
int MATRIX::qr( Real *A, int nRows, int nCols,
                Real *extraData, Real *workspace,
                int workSize )
{
  int                  lwork = workSize;
  int                  info;

#if 0
#ifdef SINGLE_PRECISION
  sgeqrf( &nRows, &nCols, A, &nCols, extraData, workspace, &lwork, &info );
#else
  dgeqrf( &nRows, &nCols, A, &nCols, extraData, workspace, &lwork, &info );
#endif
#endif
#ifdef SINGLE_PRECISION
  info = LAPACKE_sgeqrf( LAPACK_ROW_MAJOR, nRows, nCols, A, nCols, extraData );
#else
  info = LAPACKE_dgeqrf( LAPACK_ROW_MAJOR, nRows, nCols, A, nCols, extraData );
#endif

  return info;
}

//////////////////////////////////////////////////////////////////////
// Performs QR factorization with column pivoting on the given
// system.  Other calls will be needed following this to work with
// the factorization, since this is just a wrapper for the basic
// LAPACK call.
//////////////////////////////////////////////////////////////////////
int MATRIX::qrPivot( Real *A, int *pivots, int nRows, int nCols,
                     Real *extraData, Real *workspace,
                     int workSize )
{
  int                  lwork = workSize;
  int                  info;

  IntArray             newPivots;

  if ( pivots == NULL ) {
    newPivots.resize( nCols );
    pivots = newPivots.data();
  }

#ifdef SINGLE_PRECISION
  info = LAPACKE_sgeqp3( LAPACK_ROW_MAJOR, nRows, nCols, A, nCols,
                         pivots, extraData );
#else
  info = LAPACKE_dgeqp3( LAPACK_ROW_MAJOR, nRows, nCols, A, nCols,
                         pivots, extraData );
#endif
}

//////////////////////////////////////////////////////////////////////
// Extracts the leading columns of a QR factor computed
// via the qr function.  Parameters should be the same.
//////////////////////////////////////////////////////////////////////
int MATRIX::extractQRfactor( Real *A, int nRows, int nCols,
                             Real *extraData, Real *workspace,
                             int workSize )
{
  int                  lwork = workSize;
  int                  info;

#if 0
#ifdef SINGLE_PRECISION
  sorgqr( &nRows, &nCols, &nCols, A, &nCols, extraData,
          workspace, &workSize, &info );
#else
  dorgqr( &nRows, &nCols, &nCols, A, &nCols, extraData,
          workspace, &workSize, &info );
#endif
#endif
#ifdef SINGLE_PRECISION
  info = LAPACKE_sorgqr( LAPACK_ROW_MAJOR, nRows, nCols, nCols,
                         A, nCols, extraData );
#else
  info = LAPACKE_dorgqr( LAPACK_ROW_MAJOR, nRows, nCols, nCols,
                         A, nCols, extraData );
#endif

  return info;
}

//////////////////////////////////////////////////////////////////////
// copy this matrix to MATRIX3 type
//////////////////////////////////////////////////////////////////////
void MATRIX::copiesInto(MATRIX3& matrix3)
{
  if (_rows != 3 || _cols != 3)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << "Trying to copy a MATRIX that is not 3x3 to a MATRIX3!" << endl;
    return;
  }
  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++)
      matrix3(x,y) = (*this)(x,y);
}

Real MATRIX::differenceFrobeniusSq( MATRIX& B )
{
	// Not MKL optimized
	Real frobSq = 0.0;
	for( int i = 0; i < _rows; i++ ) {
		for( int j = 0; j < _cols; j++ ) {
			Real diff = (*this)(i,j) - B(i,j);
			frobSq += diff*diff;
		}
	}
	return frobSq;
}

//////////////////////////////////////////////////////////////////////
// Static routines for handling array-based matrices.
// Use with case - nothing here does any bound checking
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Copy one matrix in to another (A <- B)
//////////////////////////////////////////////////////////////////////
void MATRIX::copy( Real *A, const Real *B, int rows, int cols )
{
  memcpy( A, B, rows * cols * sizeof( Real ) );
}

//////////////////////////////////////////////////////////////////////
// Copy operation in which the leading dimension of the matrix
// changes
//////////////////////////////////////////////////////////////////////
void MATRIX::copy( Real *A, const Real *B, int rows, int cols,
                   int ldaA, int ldaB )
{
  const Real                *inputRow;
  Real                      *outputRow;

  TRACE_ASSERT( ldaA >= cols && ldaB >= cols, "Invalid leading dimension" );

  for ( int row_idx = 0; row_idx < rows; row_idx++ )
  {
    inputRow = B + row_idx * ldaB;
    outputRow = A + row_idx * ldaA;

#ifdef SINGLE_PRECISION
  cblas_scopy( cols, inputRow, 1, outputRow, 1 );
#else
  cblas_dcopy( cols, inputRow, 1, outputRow, 1 );
#endif
  }
}

//////////////////////////////////////////////////////////////////////
// Copy a 3x3 matrix
//////////////////////////////////////////////////////////////////////
void MATRIX::copy( Real *A, const MATRIX3 &B )
{
  for ( int i = 0; i < 3; i++ )
  for ( int j = 0; j < 3; j++ )
  {
    access( A, 3, 3, i, j ) = B( i, j );
  }
}

//////////////////////////////////////////////////////////////////////
// Copies a row (or part of it) from one matrix to another
//////////////////////////////////////////////////////////////////////
void MATRIX::copyRow( const Real *A, Real *B,
                      int rowA, int rowB,
                      int nColsA, int nColsB,
                      int nCopyCols )
{
  if ( nCopyCols < 0 )
  {
    nCopyCols = nColsA;
  }

  memcpy( (void *)( B + rowB * nColsB ), (void *)( A + rowA * nColsA ),
          nCopyCols * sizeof( Real ) );
}

//////////////////////////////////////////////////////////////////////
// Adds one row in B to another in A
//////////////////////////////////////////////////////////////////////
void MATRIX::addRow( const Real *A, Real *B,
                     int rowA, int rowB,
                     int nColsA, int nColsB,
                     int nCopyCols,
                     Real alpha )
{
  if ( nCopyCols < 0 )
  {
    nCopyCols = nColsA;
  }

  MATRIX::axpy( B + rowB * nColsB, // Row to add to
                A + rowA * nColsA, // Row to copy from
                1, nCopyCols, // Dimensions of row to copy
                alpha,
                // Both leading dimensions are 1
                1, 1 );
}

//////////////////////////////////////////////////////////////////////
// Copies a set of rows r0, r1, r2, ... from A to rows 0, 1, 2, ... of B
//////////////////////////////////////////////////////////////////////
void MATRIX::copyRows( const Real *A, Real *B,
                       const IntArray &rowsA, int nCols )
{
  for ( int row_idx = 0; row_idx < rowsA.size(); row_idx++ )
  {
    MATRIX::copyRow( A, B, rowsA[ row_idx ], row_idx, nCols, nCols );
  }
}

//////////////////////////////////////////////////////////////////////
// Copies a subset of the rows r0, r1, r2, specified by the given
// index range
//////////////////////////////////////////////////////////////////////
void MATRIX::copyRows( const Real *A, Real *B,
                       const IntArray &rowsA, int nCols,
                       const IndexRange &rowRange,
                       int offset )
{
  for ( int row_idx = rowRange.first; row_idx <= rowRange.second; row_idx++ )
  {
#ifdef TEST_ACCESS
  TRACE_ASSERT( row_idx - rowRange.first >= 0 );
  TRACE_ASSERT( rowsA[ row_idx ] - offset >= 0 );
#endif
    MATRIX::copyRow( A, B,
                     rowsA[ row_idx ] - offset,
                     row_idx - rowRange.first,
                     nCols, nCols );
  }
}

//////////////////////////////////////////////////////////////////////
// Adds rows r0, r1, r2, ... of A to rows 0, 1, 2, ... of B
//////////////////////////////////////////////////////////////////////
void MATRIX::copyAddRows( const Real *A, Real *B,
                          const IntArray &rowsA, int nCols,
                          const IndexRange &rowRange,
                          int offset, Real alpha )
{
  for ( int row_idx = rowRange.first; row_idx <= rowRange.second; row_idx++ )
  {
#ifdef TEST_ACCESS
  TRACE_ASSERT( row_idx - rowRange.first >= 0 );
  TRACE_ASSERT( rowsA[ row_idx ] - offset >= 0 );
#endif
    MATRIX::addRow( A, B,
                    // Row to copy from
                    rowsA[ row_idx ] - offset,
                    // Row to copy to
                    row_idx - rowRange.first,
                    nCols, nCols,
                    // Columns to copy
                    nCols,
                    alpha );
  }
}

//////////////////////////////////////////////////////////////////////
// Copies rows 0, 1, 2, of A in to r0, r1, r2, of B
//////////////////////////////////////////////////////////////////////
void MATRIX::scatterRows( const Real *A, Real *B,
                          const IntArray &rowsB, int nCols )
{
  for ( int row_idx = 0; row_idx < rowsB.size(); row_idx++ )
  {
    MATRIX::copyRow( A, B, row_idx, rowsB[ row_idx ], nCols, nCols );
  }
}

//////////////////////////////////////////////////////////////////////
// Scatters rows 0, 1, ... of A to a subset of r0, r1, r2 of B specified
// with the given index range.
//////////////////////////////////////////////////////////////////////
void MATRIX::scatterRows( const Real *A, Real *B,
                          const IntArray &rowsB, int nCols,
                          const IndexRange &rowRange,
                          // Optional offset for row indices
                          int offset )
{
  for ( int row_idx = rowRange.first; row_idx <= rowRange.second; row_idx++ )
  {
#ifdef TEST_ACCESS
  TRACE_ASSERT( row_idx - rowRange.first >= 0 );
  TRACE_ASSERT( rowsB[ row_idx ] - offset >= 0 );
#endif
    MATRIX::copyRow( A, B,
                     // Row to copy from
                     row_idx - rowRange.first,
                     // Row to copy to
                     rowsB[ row_idx ] - offset,
                     nCols, nCols );
  }
}

//////////////////////////////////////////////////////////////////////
// Scatters and adds rows 0, 1, ... of A to a subset of r0, r1, r2 of B
// specified with the given index range
//////////////////////////////////////////////////////////////////////
void MATRIX::scatterAddRows( const Real *A, Real *B,
                             const IntArray &rowsB, int nCols,
                             const IndexRange &rowRange,
                             // Optional offset for row indices
                             int offset,
                             Real alpha )
{
  for ( int row_idx = rowRange.first; row_idx <= rowRange.second; row_idx++ )
  {
#ifdef TEST_ACCESS
  TRACE_ASSERT( row_idx - rowRange.first >= 0 );
  TRACE_ASSERT( rowsB[ row_idx ] - offset >= 0 );
#endif
    MATRIX::addRow( A, B,
                    // Row to copy from
                    row_idx - rowRange.first,
                    // Row to copy to
                    rowsB[ row_idx ] - offset,
                    nCols, nCols,
                    // Columns to copy
                    nCols,
                    alpha );
  }
}

//////////////////////////////////////////////////////////////////////
// In place row swap
//////////////////////////////////////////////////////////////////////
void MATRIX::swapRows( Real *A, int nRows, int nCols, int row1, int row2 )
{
  for ( int col_idx = 0; col_idx < nCols; col_idx++ ) {
    Real                     tmp;
    
    tmp = MATRIX::access( A, nRows, nCols, row1, col_idx );

    MATRIX::access( A, nRows, nCols, row1, col_idx )
      = MATRIX::access( A, nRows, nCols, row2, col_idx );
    MATRIX::access( A, nRows, nCols, row2, col_idx ) = tmp;
  }
}

//////////////////////////////////////////////////////////////////////
// Copy a column (or part of it) from one matrix to another
//////////////////////////////////////////////////////////////////////
void MATRIX::copyColumn( const Real *A, Real *B,
                         int colA, int colB,
                         int nRowsA,
                         int nColsA, int nColsB,
                         int nCopyRows )
{
  if ( nCopyRows < 0 )
  {
    nCopyRows = nRowsA;
  }

#ifdef SINGLE_PRECISION
  cblas_scopy( nCopyRows,
               A + colA,
               nColsA, /* Leading dimension due to row major ordering */
               B + colB,
               nColsB );
#else
  cblas_dcopy( nCopyRows,
               A + colA,
               nColsA, /* Leading dimension due to row major ordering */
               B + colB,
               nColsB );
#endif
}

//////////////////////////////////////////////////////////////////////
// Adds a column from A to a column from B
//////////////////////////////////////////////////////////////////////
void MATRIX::addColumn( const Real *A, Real *B,
                        int colA, int colB,
                        int nRowsA,
                        int nColsA, int nColsB,
                        int nCopyRows,
                        Real alpha )
{
  if ( nCopyRows < 1 )
  {
    nCopyRows = nRowsA;
  }

  MATRIX::axpy( B + colB, // Column to copy to
                A + colA, // Column to copy from
                nCopyRows, 1, // Dimensions to copy
                alpha,
                // Leading dimensions are number of columns
                nColsB, nColsA );
}

//////////////////////////////////////////////////////////////////////
// Copies a set of columns c0, c1, c2, ... from A to columns 0, 1, 2, ... of B
//////////////////////////////////////////////////////////////////////
void MATRIX::copyColumns( const Real *A, Real *B,
                          const IntArray &colsA, int nRows,
                          int nColsA, int nColsB )
{
  for ( int col_idx = 0; col_idx < colsA.size(); col_idx++ )
  {
#ifdef TEST_ACCESS
    TRACE_ASSERT( colsA[ col_idx ] >= 0 && colsA[ col_idx ] < nColsA );
    TRACE_ASSERT( col_idx >= 0 && col_idx < nColsB );
#endif
    MATRIX::copyColumn( A, B, colsA[ col_idx ], col_idx,
                        nRows, nColsA, nColsB );
  }
}

//////////////////////////////////////////////////////////////////////
// Copies a subset of the columns c0, c1, c2, ..., specified by the
// given index range.
//////////////////////////////////////////////////////////////////////
void MATRIX::copyColumns( const Real *A, Real *B,
                          const IntArray &colsA, int nRows,
                          int nColsA, int nColsB,
                          const IndexRange &columnRange,
                          // Optional offset for column indices
                          int offset )
{
  for ( int col_idx = columnRange.first; col_idx <= columnRange.second;
        col_idx++ )
  {
#ifdef TEST_ACCESS
    TRACE_ASSERT( colsA[ col_idx ] - offset >= 0
               && colsA[ col_idx ] - offset < nColsA );
    TRACE_ASSERT( col_idx - columnRange.first >= 0
               && col_idx - columnRange.first < nColsB );
#endif
    MATRIX::copyColumn( A, B,
                        // Input column
                        colsA[ col_idx ] - offset,
                        // Output column
                        col_idx - columnRange.first,
                        nRows, nColsA, nColsB );
  }
}

//////////////////////////////////////////////////////////////////////
// Copies rows 0, 1, 2, of A in to c0, c1, c2, ... of B
//////////////////////////////////////////////////////////////////////
void MATRIX::scatterColumns( const Real *A, Real *B,
                             const IntArray &colsB, int nRows,
                             int nColsA, int nColsB )
{
  for ( int col_idx = 0; col_idx < colsB.size(); col_idx++ )
  {

#ifdef TEST_ACCESS
    TRACE_ASSERT( colsB[ col_idx ] >= 0 && colsB[ col_idx ] < nColsB );
    TRACE_ASSERT( col_idx >= 0 && col_idx < nColsA );
#endif
    MATRIX::copyColumn( A, B, col_idx, colsB[ col_idx ],
                        nRows, nColsA, nColsB );
  }
}

//////////////////////////////////////////////////////////////////////
// Scatters rows 0, 1, ... of A to a subset of c0, c1, ... of B specified
// with the given index range.
//////////////////////////////////////////////////////////////////////
void MATRIX::scatterColumns( const Real *A, Real *B,
                             const IntArray &colsB, int nRows,
                             int nColsA, int nColsB,
                             const IndexRange &columnRange,
                             // Optional offset for column indices
                             int offset )
{
  for ( int col_idx = columnRange.first; col_idx <= columnRange.second;
        col_idx++ )
  {
#ifdef TEST_ACCESS
    TRACE_ASSERT( colsB[ col_idx ] - offset >= 0
               && colsB[ col_idx ] - offset < nColsB );
    TRACE_ASSERT( col_idx - columnRange.first >= 0
               && col_idx - columnRange.first < nColsA );
#endif
    MATRIX::copyColumn( A, B,
                        // Column to copy from
                        col_idx - columnRange.first,
                        // Column to copy to
                        colsB[ col_idx ] - offset,
                        nRows, nColsA, nColsB );
  }
}

//////////////////////////////////////////////////////////////////////
// Scatters columns 0, 1, ... of A to a subset of c0, c1, ... of B specified
// with the given index range.
//////////////////////////////////////////////////////////////////////
void MATRIX::scatterAddColumns( const Real *A, Real *B,
                                const IntArray &colsB, int nRows,
                                int nColsA, int nColsB,
                                const IndexRange &columnRange,
                                // Optional offset for column indices
                                int offset,
                                Real alpha )
{
  for ( int col_idx = columnRange.first; col_idx <= columnRange.second;
        col_idx++ )
  {
#ifdef TEST_ACCESS
    TRACE_ASSERT( colsB[ col_idx ] - offset >= 0
               && colsB[ col_idx ] - offset < nColsB );
    TRACE_ASSERT( col_idx - columnRange.first >= 0
               && col_idx - columnRange.first < nColsA );
#endif
    MATRIX::addColumn( A, B,
                       // Column to copy from
                       col_idx - columnRange.first,
                       // Column to copy to
                       colsB[ col_idx ] - offset,
                       nRows, nColsA, nColsB,
                       // Rows to copy
                       nRows,
                       alpha );
  }
}

//////////////////////////////////////////////////////////////////////
// For a symmetric matrix in which the given rows/columns are filled
// In place column swap
//////////////////////////////////////////////////////////////////////
void MATRIX::swapColumns( Real *A, int nRows, int nCols, int col1, int col2 )
{
  for ( int row_idx = 0; row_idx < nRows; row_idx++ ) {
    Real                     tmp;

    tmp = MATRIX::access( A, nRows, nCols, row_idx, col1 );

    MATRIX::access( A, nRows, nCols, row_idx, col1 )
      = MATRIX::access( A, nRows, nCols, row_idx, col2 );
    MATRIX::access( A, nRows, nCols, row_idx, col2 ) = tmp;
  }
}

//////////////////////////////////////////////////////////////////////
// For a symmetric matrix in which the given rows/columns are filled
// in, scatter and add the entries to the given full matrix.
//
// Only updates the lower triangular part of B.
//////////////////////////////////////////////////////////////////////
void MATRIX::scatterAddSymmetricMatrix( const Real *A, Real *B,
                                        const IntArray &colsB, int nRows,
                                        const IndexRange &columnRange,
                                        // Optional offset
                                        int offset,
                                        Real alpha )
{
#ifdef DO_FLOP_COUNT
  MATRIX_FLOP_COUNT += range_size( columnRange ) * range_size( columnRange );
#endif

  for ( int row_idx = columnRange.first; row_idx <= columnRange.second;
        row_idx++ )
  {
    for ( int col_idx = columnRange.first; col_idx <= row_idx; col_idx++ )
    {
      MATRIX::access( B, nRows, nRows,
                      // Indices in output matrix
                      colsB[ row_idx ] - offset,
                      colsB[ col_idx ] - offset )
        += alpha * MATRIX::access( A,
                                   // Size of input matrix
                                   range_size( columnRange ),
                                   range_size( columnRange ),
                                   // Indices in input matrix
                                   row_idx - columnRange.first,
                                   col_idx - columnRange.first );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Matrix vector multiply for a matrix in compressed column form
// (ie. only a subset of the columns are stored and the others are
// assumed to be 0).
//////////////////////////////////////////////////////////////////////
void MATRIX::compressedColumnMult( const Real *A, const Real *b, Real *c,
                                   const IntArray &colsA,
                                   int nRows, int nCols,
                                   bool transpose,
                                   Real alpha, Real beta )
{
  if ( transpose )
  {
    // Scale c by beta
    MATRIX::scale( c, nCols, 1, beta );

#ifdef DO_FLOP_COUNT
  // Only count the flops we are directly responsible for
  MATRIX_FLOP_COUNT += 2 * colsA.size();
#endif

    for ( int col_idx = 0; col_idx < colsA.size(); col_idx++ )
    {
      c[ colsA[ col_idx ] ]
        += alpha * VECTOR::dot( nRows, A + col_idx /* start of column */, b,
                                colsA.size() /* leading dimension */, 1 );
    }
  }
  else
  {
    // Scale c by beta
    MATRIX::scale( c, nRows, 1, beta );

    for ( int col_idx = 0; col_idx < colsA.size(); col_idx++ )
    {
      MATRIX::axpy( c, A + col_idx, /* Align with column start */
                    nRows, 1, alpha * b[ colsA[ col_idx ] ],
                    1, colsA.size() /* Leading dimension */ );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Adds a scalar to the diagonal of a square matrix
//
// eg. A(i,i) += alpha, for i = 1:nRows
//////////////////////////////////////////////////////////////////////
void MATRIX::addToDiagonal( Real *A, int nRows, Real alpha )
{
#ifdef DO_FLOP_COUNT
  // Only count the flops we are directly responsible for
  MATRIX_FLOP_COUNT += nRows;
#endif

  for ( int i = 0; i < nRows; i++ )
  {
    MATRIX::access( A, nRows, nRows, i, i ) += alpha;
  }
}

//////////////////////////////////////////////////////////////////////
// Add one matrix to another (A = A + alpha * B
//////////////////////////////////////////////////////////////////////
void MATRIX::axpy( Real *A, const Real *B, int rows, int cols,
                   Real alpha, int incA, int incB )
{
#ifdef DO_FLOP_COUNT
  MATRIX_FLOP_COUNT += 2 * rows * cols;
#endif

#if 0
#ifdef SINGLE_PRECISION
	cblas_saxpy(rows * cols, alpha, B, incB, A, incA);
#else
	cblas_daxpy(rows * cols, alpha, B, incB, A, incA);
#endif
#endif
  for ( int i = 0; i < rows * cols; i++ ) {
    A[ i * incA ] += alpha * B[ i * incB ];
  }
}

//////////////////////////////////////////////////////////////////////
// Matrix-matrix multiplication (C = beta * C + alpha * A * B)
//////////////////////////////////////////////////////////////////////
void MATRIX::gemm( const Real *A, const Real *B, Real *C,
                   int rowsA, int colsA, int rowsB, int colsB,
                   bool transposeA, bool transposeB,
                   Real alpha, Real beta,
                   int ldaA, int ldaB, int ldaC )
{
  int requiredCols = transposeB ? rowsB : colsB;

  if ( ldaA < 0 ) {
    ldaA = colsA;
  }
  if ( ldaB < 0 ) {
    ldaB = colsB;
  }
  if ( ldaC < 0 ) {
    ldaC = requiredCols;
  }

  // FIXME:
  TRACE_ASSERT( ldaA >= colsA );
  TRACE_ASSERT( ldaB >= colsB );
  TRACE_ASSERT( ldaC >= requiredCols );

#ifdef DO_FLOP_COUNT
  MATRIX_FLOP_COUNT += 2 * ( transposeA ? colsA : rowsA )
                         * ( transposeB ? rowsB : colsB )
                         * ( transposeA ? rowsA : colsA );
#endif

#ifdef SINGLE_PRECISION
	cblas_sgemm(
			CblasRowMajor,
			transposeA ? CblasTrans : CblasNoTrans,
			transposeB ? CblasTrans : CblasNoTrans,
			transposeA ? colsA : rowsA, 
      transposeB ? rowsB : colsB,
      transposeA ? rowsA : colsA, 
			alpha,
			A, ldaA,
			B, ldaB,
			beta,
			C, ldaC);
#else
	cblas_dgemm(
			CblasRowMajor,
			transposeA ? CblasTrans : CblasNoTrans,
			transposeB ? CblasTrans : CblasNoTrans,
			transposeA ? colsA : rowsA, 
      transposeB ? rowsB : colsB,
      transposeA ? rowsA : colsA, 
			alpha,
			A, ldaA,
			B, ldaB,
			beta,
			C, ldaC);
#endif
}

//////////////////////////////////////////////////////////////////////
// Matrix-vector multiplication (C := beta * C + alpha * A * b)
//////////////////////////////////////////////////////////////////////
void MATRIX::gemv( const Real *A, const Real *b, Real *c,
                   int rowsA, int colsA,
                   bool transposeA,
                   Real alpha, Real beta )
{
#ifdef DO_FLOP_COUNT
  MATRIX_FLOP_COUNT += 2 * rowsA * colsA;
#endif

#ifdef SINGLE_PRECISION
  cblas_sgemv(
      CblasRowMajor,
      transposeA ? CblasTrans : CblasNoTrans,
      rowsA,
      colsA,
      alpha,
      A, colsA,
      b, 1,
      beta,
      c, 1);
#else
  cblas_dgemv(
      CblasRowMajor,
      transposeA ? CblasTrans : CblasNoTrans,
      rowsA,
      colsA,
      alpha,
      A, colsA,
      b, 1,
      beta,
      c, 1);
#endif
#if 0
  if (transposeA) {
    for (int j = 0; j < colsA; j++) {
      c[j] *= beta;
    }

    for (int i = 0; i < rowsA; i++) {
      const Real *base = A + i * colsA;
      Real        coef = b[i];

      for (int j = 0; j < colsA; j++) {
        c[j] += alpha * base[j] * coef;
      }
    }
  } else {
    for (int i = 0; i < rowsA; i++) {
      const Real *base = A + i * colsA;

      c[i] *= beta;

      for (int j = 0; j < colsA; j++) {
        c[i] += alpha * base[j] * b[j];
      }
    }
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Symmetrix matrix-matrix update (C := beta * C + alpha * A * A')
// Updates the lower-triangular part of C
//
// If trans == true then C := beta * C + alpha * A' * A
//////////////////////////////////////////////////////////////////////
void MATRIX::syrk( const Real *A, Real *C,
                   int rowsC, int k, /* A is either n x k or k x n */
                   bool transpose, Real alpha, Real beta,
                   int ldaA, int ldaC )
{
  int requiredCols = transpose ? rowsC : k;

  if ( ldaA < 0 ) {
    ldaA = requiredCols;
  } if ( ldaC < 0 ) {
    ldaC = rowsC;
  }

#ifdef DO_FLOP_COUNT
  MATRIX_FLOP_COUNT += k * rowsC * ( rowsC + 1 );
#endif

#ifdef SINGLE_PRECISION
  cblas_ssyrk(
      CblasRowMajor,
      CblasLower, /* Work with lower triangular part of C */
      transpose ? CblasTrans : CblasNoTrans,
      rowsC, k,
      alpha,
      A, ldaA,
      beta,
      C, ldaC);
#else
  cblas_dsyrk(
      CblasRowMajor,
      CblasLower, /* Work with lower triangular part of C */
      transpose ? CblasTrans : CblasNoTrans,
      rowsC, k,
      alpha,
      A, ldaA,
      beta,
      C, ldaC);
#endif
}

//////////////////////////////////////////////////////////////////////
// Get transpose (B = A^T)
//////////////////////////////////////////////////////////////////////
void MATRIX::transpose( Real *A, const Real *B, int rows, int cols )
{
  for (int y = 0; y < cols; y++)
    for (int x = 0; x < rows; x++)
      access( A, cols, rows, y, x ) = access( B, rows, cols, x, y );
}

//////////////////////////////////////////////////////////////////////
// Get transpose (B = A^T)
//////////////////////////////////////////////////////////////////////
void MATRIX::transposeBLAS( Real *A, const Real *B, int rows, int cols,
                            int ldaA, int ldaB )
{
  if ( ldaA < 0 )
  {
    ldaA = rows;
  }

  if ( ldaB < 0 )
  {
    ldaB = cols;
  }

#ifdef SINGLE_PRECISION
  for ( int x = 0; x < rows; x++ )
  {
    cblas_scopy( cols, B + x * ldaB, 1, A + x, ldaA );
  }
#else
  for ( int x = 0; x < rows; x++ )
  {
    cblas_dcopy( cols, B + x * ldaB, 1, A + x, ldaA );
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// In place transpose for square matrices
//////////////////////////////////////////////////////////////////////
void MATRIX::transpose( Real *A, int rows )
{
  for ( int i = 1; i < rows; i++ )
  {
    for ( int j = 0; j < i; j++ )
    {
      Real temp = access( A, rows, rows, i, j );
      access( A, rows, rows, i, j ) = access( A, rows, rows, j, i );
      access( A, rows, rows, j, i ) = temp;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Compute eigenvalues/vectors.  We require everything here,
// including workspaces, etc.
// workspace should have size (7 * rows)
// vectorWorkspace should have size (rows * rows)
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem( Real *A, int rows,
                          Real *eigenvalues, Real *eigenvectors,
                          Real *workspace, Real *vectorWorkspace )
{
  cerr << "MATRIX::eigensystem not implemented" << endl;

#if 0
  int rowsize = rows;
  int worksize = 5 * rows;

  Real *work = workspace;
  Real *valuesReal = workspace + worksize;
  Real *valuesImag = valuesReal + rows;
  Real *vectors = vectorWorkspace;

  // the actual LAPACK call
  int error;
#ifdef SINGLE_PRECISION
  sgeev("V","N", &rowsize, A, &rowsize, 
        valuesReal, valuesImag, 
        vectors, &rowsize, NULL, &rowsize,
        work, &worksize, &error);
#else
  dgeev("V","N", &rowsize, A, &rowsize, 
        valuesReal, valuesImag, 
        vectors, &rowsize, NULL, &rowsize,
        work, &worksize, &error);
#endif

  if (error != 0)
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " eigenvalue solver bombed!" << endl;
  }

  // copy out results
  for (int x = 0; x < rows; x++)
    eigenvalues[x] = valuesReal[x];

  for (int x = 0; x < rows; x++)
    for (int y = 0; y < rows; y++)
      MATRIX::access( eigenvectors, rows, rows, x, y ) = vectors[x + y * rows];
#endif
}

//////////////////////////////////////////////////////////////////////
// LU factorization - requires the matrix, and a workspace for
// pivoting information
//////////////////////////////////////////////////////////////////////
int MATRIX::LU( Real *A, int *pivotData, int rows )
{
  int info = 0;

#ifdef SINGLE_PRECISION
  info = LAPACKE_sgetrf( LAPACK_ROW_MAJOR, rows, rows, A, rows, pivotData );
#else
  info = LAPACKE_dgetrf( LAPACK_ROW_MAJOR, rows, rows, A, rows, pivotData );
#endif

  return info;
}

//////////////////////////////////////////////////////////////////////
// Solve via an LU factorization
//////////////////////////////////////////////////////////////////////
int MATRIX::LUsolve( const Real *A, const int *pivotData,
                     Real *B, int rowsB, int colsB,
                     bool transposeA,
                     int ldaA, int ldaB )
{
  int info = 0;

  if ( ldaB < colsB ) {
    ldaB = colsB;
  }

  if ( ldaA < rowsB ) {
    ldaA = rowsB;
  }

#ifdef SINGLE_PRECISION
  info = LAPACKE_sgetrs( LAPACK_ROW_MAJOR,
                         transposeA ? 'T' : 'N',
                         rowsB, colsB,
                         A, ldaA, pivotData,
                         B, ldaB );
#else
  info = LAPACKE_dgetrs( LAPACK_ROW_MAJOR,
                         transposeA ? 'T' : 'N',
                         rowsB, colsB,
                         A, ldaA, pivotData,
                         B, ldaB );
#endif
}

//////////////////////////////////////////////////////////////////////
// Computes the Cholesky factor of the matrix stored in A.  A is
// overwritten
//////////////////////////////////////////////////////////////////////
int MATRIX::cholesky( Real *A, int rows )
{
  int info = 0;

#if 0
#ifdef SINGLE_PRECISION
  spotrf( "U", &rows, A, &rows, &info );
#else
  dpotrf( "U", &rows, A, &rows, &info );
#endif
#endif

#ifdef DO_FLOP_COUNT
  MATRIX_FLOP_COUNT += (size_t)( (Real)rows * ( 1.0 / 6.0 +
                                 (Real)rows * ( 1.0 / 2.0 + 
                                 (Real)rows * ( 1.0 / 3.0 ) ) ) );
#endif

#ifdef SINGLE_PRECISION
  info = LAPACKE_spotrf( LAPACK_ROW_MAJOR, 'L', rows, A, rows );
#else
  info = LAPACKE_dpotrf( LAPACK_ROW_MAJOR, 'L', rows, A, rows );
#endif

  return info;
}

//////////////////////////////////////////////////////////////////////
// Computes a block LDL^{T} factorization using the Bunch-Kaufman method.
//
// Doesn't actual extract the factors from the matrix.... yet
//////////////////////////////////////////////////////////////////////
int MATRIX::LDL( Real *A, int rows, IntArray &pivotData )
{
  int info = 0;

  TRACE_ASSERT( pivotData.size() == rows, "Invalid pivotData size" );

#ifdef SINGLE_PRECISION
  info = LAPACKE_ssytrf( LAPACK_ROW_MAJOR, 'L', rows, A, rows,
                         pivotData.data() );
#else
  info = LAPACKE_dsytrf( LAPACK_ROW_MAJOR, 'L', rows, A, rows,
                         pivotData.data() );
#endif

  return info;
}

#if 0
//////////////////////////////////////////////////////////////////////
// Runs the LAPACK conversion routine for a LDL^T factorization
//////////////////////////////////////////////////////////////////////
int MATRIX::convertLDL( Real *A, int rows, const IntArray &permutation,
                        FloatArray &workspace )
{
  int info = 0;

  TRACE_ASSERT( permutation.size() == rows, "Invalid permutation size" );

  if ( workspace.size() != permutation.size() ) {
    workspace.resize( permutation.size() );
  }

  // Convert the column major ordering, since there is no CLAPACK
  // implementation of the conversion routine, apparently
  MATRIX::transpose( A, rows );

#ifdef SINGLE_PRECISION

#else
  dsyconv( "L", "C", &rows, A, &rows, permutation.data(), workspace.data(),
           &info );
#endif

  // Convert back to row major ordering
  MATRIX::transpose( A, rows );

  return info;
}
#endif

//////////////////////////////////////////////////////////////////////
// Builds a full permutation matrix (possibly overkill, but let's
// do it anyways)
//////////////////////////////////////////////////////////////////////
void MATRIX::buildLDLPermutation( const IntArray &pivotData,
                                  IntArray &permutation,
                                  Real *P )
{
  int                        nRows = pivotData.size();

  if ( permutation.size() != pivotData.size() )
  {
    permutation.resize( pivotData.size() );
  }

  // Initialize the identity permutation
  for ( int i = 0; i < permutation.size(); i++ )
  {
    permutation[ i ] = i;
  }

#if 0
  // Process swaps from the pivot data
  for ( int i = 0; i < pivotData.size(); i++ )
  {
    int                      swapInfo = pivotData[ i ];

    if ( swapInfo < 0 )
    {
      // Must have a 2x2 diagonal block - skip to the next entry
      i += 1;

      swapInfo *= -1;
      swapInfo -= 1;

      // Swap these permutation entries
      int                    oldValue = permutation.at( i );

      permutation.at( i ) = permutation.at( swapInfo );
      permutation.at( swapInfo ) = oldValue;
    }
    else if ( swapInfo > 0 )
    {
      swapInfo -= 1;

      int                    oldValue = permutation.at( i );

      permutation.at( i ) = permutation.at( swapInfo );
      permutation.at( swapInfo ) = oldValue;
    }
    else
    {
      TRACE_ASSERT( NULL, "Should never get here" );
    }
  }
#endif

  // Process in reverse order
  for ( int i = pivotData.size() - 1; i >= 0; i-- )
  {
    int                      swapInfo = pivotData[ i ];

    if ( swapInfo < 0 )
    {
      // 2x2 diagonal block
      swapInfo *= -1;
      swapInfo -= 1;

      // Swap these permutation entries
      int                    oldValue = permutation.at( i );

      permutation.at( i ) = permutation.at( swapInfo );
      permutation.at( swapInfo ) = oldValue;

      // To accomodate the 2x2 block
      i -= 1;
    }
    else if ( swapInfo > 0 )
    {
      swapInfo -= 1;

      int                    oldValue = permutation.at( i );

      permutation.at( i ) = permutation.at( swapInfo );
      permutation.at( swapInfo ) = oldValue;
    }
    else
    {
      TRACE_ASSERT( NULL, "Should never get here" );
    }
  }

  if ( P != NULL )
  {
    MATRIX::clear( P, nRows, nRows );

#if 0
    for ( int i = 0; i < nRows; i++ )
    {
      MATRIX::access( P, nRows, nRows,
                      i, permutation[ i ] ) = 1.0;
                      //permutation[ i ], i ) = 1.0;
    }
#endif

    // Try contstructing the permutation explicitly
    for ( int i = 0; i < nRows; i++ ) {
      MATRIX::access( P, nRows, nRows, i, i ) = 1.0;
    }

    for ( int i = pivotData.size() - 1; i >= 0; i-- )
    {
      int                      swapInfo = pivotData[ i ];

      if ( swapInfo < 0 )
      {
        // 2x2 diagonal block
        swapInfo *= -1;
        swapInfo -= 1;

        // Swap these two rows
        for ( int col_idx = 0; col_idx < nRows; col_idx++ ) {
          Real     tmp = MATRIX::access( P, nRows, nRows, i, col_idx );

          MATRIX::access( P, nRows, nRows, i, col_idx )
            = MATRIX::access( P, nRows, nRows, swapInfo, col_idx );
          MATRIX::access( P, nRows, nRows, swapInfo, col_idx ) = tmp;
        }

        // To accomodate the 2x2 block
        i -= 1;
      }
      else if ( swapInfo > 0 )
      {
        swapInfo -= 1;

        // Swap these two rows
        for ( int col_idx = 0; col_idx < nRows; col_idx++ ) {
          Real     tmp = MATRIX::access( P, nRows, nRows, i, col_idx );

          MATRIX::access( P, nRows, nRows, i, col_idx )
            = MATRIX::access( P, nRows, nRows, swapInfo, col_idx );
          MATRIX::access( P, nRows, nRows, swapInfo, col_idx ) = tmp;
        }
      }
      else
      {
        TRACE_ASSERT( NULL, "Should never get here" );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Extracts a "psychologically" lower triangular matrix from a LDL^T
// factorization; that is, a matrix which is a product of lower triangular
// and permutation matrices.
//
// Also, builds a permutation vector and its inverse along the way.
//////////////////////////////////////////////////////////////////////
void MATRIX::buildLDLTriangle( const IntArray &pivotData,
                               const Real *A,
                               IntArray &permutation,
                               IntArray &inversePermutation,
                               Real *L )
{
  // We need a workspace for extracting and permuting columns from
  // the LDL^T output
  MATRIX                     workspace( pivotData.size(), 2 );
  Real                      *workspaceData = workspace.data();

  int                        nRows = pivotData.size();
  int                        copyColumns;
  int                        swapInfo;

  if ( permutation.size() != nRows ) {
    permutation.resize( nRows );
  }

  if ( inversePermutation.size() != nRows ) {
    inversePermutation.resize( nRows );
  }

  // Initialize the permutation data
  for ( int i = 0; i < permutation.size(); i++ )
  {
    permutation[ i ] = i;
    inversePermutation[ i ] = i;
  }

  // Walk through each 1x1 or 2x2 block in the matrix
  for ( int block_start = 0; block_start < nRows; block_start++ )
  {
    swapInfo = pivotData[ block_start ];

    if ( swapInfo < 0 )
    {
      TRACE_ASSERT( swapInfo == pivotData.at( block_start + 1 ) );

      // 2x2 block
      //
      // Copy the lower triangular block to our workspace
      MATRIX::clear( workspaceData, nRows, 2 );
      MATRIX::access( workspaceData, nRows, 2, block_start, 0 ) = 1.0;
      MATRIX::access( workspaceData, nRows, 2, block_start + 1, 1 ) = 1.0;

      for ( int row = block_start + 2; row < nRows; row++ )
      {
        for ( int col = 0; col < 2; col++ )
        {
          MATRIX::access( workspaceData, nRows, 2, row, col )
            = MATRIX::access( A, nRows, nRows, row, block_start + col );
        }
      }

      // Figure out how to permute this
      swapInfo *= -1;
      swapInfo -= 1;

      // Adjust the permutation if we need to perform a swap
      if ( swapInfo != block_start + 1 )
      {
        int                  oldRow = inversePermutation[ block_start + 1 ];
        int                  newRow = inversePermutation[ swapInfo ];

        TRACE_ASSERT( swapInfo > block_start + 1 );

        //inversePermutation[ block_start + 1 ] = swapInfo;
        inversePermutation[ block_start + 1 ] = newRow;
        inversePermutation[ swapInfo ] = oldRow;

        invertIntArray( inversePermutation, permutation );
#if 0
        int                  oldRow = permutation[ block_start + 1 ];
        int                  newRow = permutation[ swapInfo ];

        permutation[ block_start + 1 ] = newRow;
        permutation[ swapInfo ] = oldRow;
#endif
      }

      copyColumns = 2;
    }
    else if ( swapInfo > 0 )
    {
      // 1x1 block
      //
      // Copy the lower triangular block to our workspace
      MATRIX::clear( workspaceData, nRows, 1 );
      MATRIX::access( workspaceData, nRows, 1, block_start, 0 ) = 1.0;

      for ( int row = block_start + 1; row < nRows; row++ )
      {
        MATRIX::access( workspaceData, nRows, 1, row, 0 )
          = MATRIX::access( A, nRows, nRows, row, block_start );
      }

      // Figure out how to permute this
      swapInfo -= 1;

      // Adjust the permutation if we need to perform a swap
      if ( swapInfo != block_start )
      {
        int                  oldRow = inversePermutation[ block_start ];
        int                  newRow = inversePermutation[ swapInfo ];

        TRACE_ASSERT( swapInfo > block_start );

        //inversePermutation[ block_start ] = swapInfo;
        inversePermutation[ block_start ] = newRow;
        inversePermutation[ swapInfo ] = oldRow;

        invertIntArray( inversePermutation, permutation );
#if 0
        int                  oldRow = permutation[ block_start ];
        int                  newRow = permutation[ swapInfo ];

        permutation[ block_start ] = newRow;
        permutation[ swapInfo ] = oldRow;
#endif
      }

      copyColumns = 1;
    }
    else
    {
      TRACE_ASSERT( NULL, "Should never get here" );
    }

    // Copy and permute rows from the workspace
    for ( int row = 0; row < nRows; row++ )
    {
      int                    row_idx = permutation[ row ];

      for ( int col = 0; col < copyColumns; col++ )
      {
        MATRIX::access( L, nRows, nRows, row, block_start + col )
          = MATRIX::access( workspaceData, nRows, copyColumns, row_idx, col );
#if 0
        MATRIX::access( L, nRows, nRows, row_idx, block_start + col )
          = MATRIX::access( workspaceData, nRows, copyColumns, row, col );
#endif
      }
    }

    if ( copyColumns == 2 ) {
      block_start += 1;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Extracts the lower triangular system from the LDL^T factorization
//////////////////////////////////////////////////////////////////////
void MATRIX::convertLDLTriangle( const IntArray &pivotData,
                                 const Real *A,
                                 Real *L )
{
  int                        rows = pivotData.size();

  MATRIX::clear( L, pivotData.size(), pivotData.size() );

  // Copy the lower triangular part of A
  for ( int i = 0; i < pivotData.size(); i++ )
  for ( int j = 0; j <= i; j++ )
  {
    MATRIX::access( L, pivotData.size(), pivotData.size(), i, j )
      = MATRIX::access( A, pivotData.size(), pivotData.size(), i, j );
  }

  // Set unit diagonal entries
  for ( int i = 0; i < pivotData.size(); i++ ) {
    MATRIX::access( L, pivotData.size(), pivotData.size(), i, i ) = 1.0;
  }

  // Zero the off-diagonal entries corresponding to 2x2 blocks
  for ( int i = 0; i < pivotData.size(); i++ )
  {
    if ( pivotData[ i ] < 0 ) {
      MATRIX::access( L, pivotData.size(), pivotData.size(), i + 1, i ) = 0.0;

      i += 1;
    }
  }

#if 0
  // Go through each row and perform swaps as necessary
  for ( int i = 0; i < pivotData.size(); i++ )
  {
    int                      pivot = pivotData[ i ];

    if ( pivot > 0 )
    {
      pivot = pivot - 1;

      for ( int j = 0; j < rows; j++ )
      //for ( int j = i; j < rows; j++ )
      //for ( int j = 0; j < i; j++ )
      {
        Real temp = MATRIX::access( L, rows, rows, pivot, j );

        MATRIX::access( L, rows, rows, pivot, j )
          = MATRIX::access( L, rows, rows, i, j );
        MATRIX::access( L, rows, rows, i, j ) = temp;
      }
    }
    else if ( pivot < 0 )
    {
      pivot = -1 * pivot;
      pivot = pivot - 1;

      for ( int j = 0; j < rows; j++ )
      //for ( int j = i; j < rows; j++ )
      //for ( int j = 0; j < i; j++ )
      {
        Real temp = MATRIX::access( L, rows, rows, pivot, j );

        MATRIX::access( L, rows, rows, pivot, j )
          = MATRIX::access( L, rows, rows, i + 1, j );
        MATRIX::access( L, rows, rows, i + 1, j ) = temp;
      }

      i = i + 1;
    }
    else
    {
      TRACE_ASSERT( NULL, "Should never get here" );
    }
  }
#endif
}

#if 0
//////////////////////////////////////////////////////////////////////
// A hacky and inefficient version of the above, which just does
// a bunch of matrix multiplication to form the psychologically lower
// triangular matrix.
//////////////////////////////////////////////////////////////////////
void MATRIX::convertLDLTHack( const IntArray &pivotData,
                              const Real *A,
                              MATRIX &L )
{
  int                        nRows = pivotData.size();
  int                        pivot;
  int                        swapRow1, swapRow2;

  L = MATRIX::ident( pivotData.size() );

  for ( int block_start = 0; block_start < pivotData.size(); block_start++ )
  {
    // Generate an identity matrix to start with, then fill in
    // the lower triangular values from L_{block_start}
    MATRIX                   L_i = MATRIX::ident( nRows );

    printf( "Hacky LDLT conversion: column %d of %d\n", block_start + 1,
            (int)pivotData.size() );

    pivot = pivotData[ block_start ];
    if ( pivot > 0 )
    {
      pivot -= 1;

      // We have a 1x1 block so just fill in this
      for ( int row = block_start + 1; row < nRows; row++ )
      {
        L_i( row, block_start )
          = MATRIX::access( A, nRows, nRows, row, block_start );
      }

      swapRow1 = block_start;
      swapRow2 = pivot;
    }
    else if ( pivot < 0 )
    {
      pivot *= -1;
      pivot -= 1;

      // We have a 2x2 block, so we need to copy 2 columns
      for ( int row = block_start + 2; row < nRows; row++ )
      for ( int col = block_start; col < block_start + 2; col++ )
      {
        L_i( row, col ) = MATRIX::access( A, nRows, nRows, row, col );
      }

      swapRow1 = block_start + 1;
#if 0
      swapRow1 = block_start;
#endif
      swapRow2 = pivot;

      block_start += 1;
    }
    else
    {
      TRACE_ASSERT( NULL, "Should never get here" );
    }

    // Swap rows, if necessary
    if ( swapRow1 != swapRow2 )
    {
#if 0
      // FIXME
      if ( pivotData[ block_start ] < 0 )
      {
        for ( int col = 0; col < nRows; col++ )
        {
          Real                 temp = L_i( swapRow1 - 1, col );

          L_i( swapRow1 - 1, col ) = L_i( swapRow2, col );
          L_i( swapRow2, col ) = temp;
        }
      }
#endif
      
      for ( int col = 0; col < nRows; col++ )
      {
        Real                 temp = L_i( swapRow1, col );

        L_i( swapRow1, col ) = L_i( swapRow2, col );
        L_i( swapRow2, col ) = temp;
      }
#if 0
      for ( int row = 0; row < nRows; row++ )
      {
        Real                 temp = L_i( row, swapRow1 );

        L_i( row, swapRow1 ) = L_i( row, swapRow2 );
        L_i( row, swapRow2 ) = temp;
      }
#endif
    }

    // Accumulate via multiplication
    L = L * L_i;
#if 0
    L = L_i * L;
#endif
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Extracts the block diagonal from a LDL^T factorization
//////////////////////////////////////////////////////////////////////
void MATRIX::extractLDLDiagonal( const IntArray &pivotData,
                                 const Real *A,
                                 Real *D )
{
  MATRIX::clear( D, pivotData.size(), pivotData.size() );

  for ( int block_start = 0; block_start < pivotData.size(); block_start++ )
  {
    if ( pivotData[ block_start ] < 0 )
    {
      // Copy a 2x2 block
      for ( int row_idx = block_start; row_idx < block_start + 2; row_idx++ )
      //for ( int col_idx = block_start; col_idx <= row_idx; col_idx++ )
      for ( int col_idx = block_start; col_idx <= row_idx; col_idx++ )
      {
        MATRIX::access( D, pivotData.size(), pivotData.size(),
                        row_idx, col_idx )
          = MATRIX::access( A, pivotData.size(), pivotData.size(),
                            row_idx, col_idx );
      }

      // Symmetrize
      MATRIX::access( D, pivotData.size(), pivotData.size(),
                      block_start, block_start + 1 )
       = MATRIX::access( D, pivotData.size(), pivotData.size(),
                         block_start + 1, block_start );

       block_start += 1;
    }
    else if ( pivotData[ block_start ] > 0 )
    {
      // Copy a single element
      MATRIX::access( D, pivotData.size(), pivotData.size(),
                      block_start, block_start )
        = MATRIX::access( A, pivotData.size(), pivotData.size(),
                          block_start, block_start );
    }
    else
    {
      TRACE_ASSERT( NULL, "Should never get here" );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Given a pivot array returned by LDL, apply the permutation implied
// by this array.  The permutation is applied to the matrix in place.
//////////////////////////////////////////////////////////////////////
void MATRIX::applyLDLPermutation( const IntArray &pivotData,
                                  Real *A, int nRows, int nCols,
                                  bool left, bool transpose )
{
  int                        swapInfo;
  int                        block_start, block_end, block_increment;

  if ( left )
  {
    int                      swapRow;

    // We will be swapping rows
    TRACE_ASSERT( nRows == pivotData.size(), "Invalid matrix size" );

    block_start = transpose ? 0 : pivotData.size() - 1;
    block_end = transpose ? pivotData.size() : -1;
    block_increment = transpose ? 1 : -1;

    for ( int pivot_block = block_start; pivot_block != block_end;
          pivot_block += block_increment )
    {
      swapInfo = pivotData[ pivot_block ];

      if ( swapInfo < 0 ) {
        swapInfo *= -1;
        swapInfo -= 1;

        // If we are moving in increasing order (ie. transposed) then
        // we need to swap the next row (pivot_block + 1)
        swapRow = transpose ? pivot_block + 1 : pivot_block;

        pivot_block += block_increment;
      }
      else {
        swapInfo -= 1;
        swapRow = pivot_block;
      }

      MATRIX::swapRows( A, nRows, nCols, swapRow, swapInfo );
    }
  }
  else
  {
    int                      swapColumn;

    // We will be swapping columns
    TRACE_ASSERT( nCols == pivotData.size(), "Invalid matrix size" );

    block_start = transpose ? pivotData.size() - 1 : 0;
    block_end = transpose ? -1 : pivotData.size();
    block_increment = transpose ? -1 : 1;

    for ( int pivot_block = block_start; pivot_block != block_end;
          pivot_block += block_increment )
    {
      swapInfo = pivotData[ pivot_block ];

      if ( swapInfo < 0 ) {
        swapInfo *= -1;
        swapInfo -= 1;

        // If we are moving in increasing order (ie. not transposed) then
        // we need to swap the next row (pivot_block + 1)
        swapColumn = transpose ? pivot_block : pivot_block + 1;

        pivot_block += block_increment;
      }
      else {
        swapInfo -= 1;
        swapColumn = pivot_block;
      }

      MATRIX::swapColumns( A, nRows, nCols, swapColumn, swapInfo );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Modified Cholesky factorization based on a Bunch-Kaufman factorization.
//
// Overwrites A with the modified Cholesky factor and fills in U and V
// with the low-rank update necessary to make A positive definite.
//////////////////////////////////////////////////////////////////////
int MATRIX::modifiedCholesky( Real *A, MATRIX &U, MATRIX &V,
                              int nRows )
{
  // Start by attemping a direct Cholesky factorization of A.
  // If it succeeds, then we can proceed.
  //
  // We need to copy A to a workspace, in case we need to recover
  // it after a failed Cholesky factorization.
  int                        info;

  int                        modificationSize = 0;

  MATRIX                     Acopy( nRows, nRows, A );

  // For now, resize the low rank modification matrices to make them
  // as big as the original matrix
  U.resizeAndWipe( nRows, nRows );
  V.resizeAndWipe( nRows, nRows );

  // FIXME: debugging
  MATRIX::write( A, nRows, nRows, "modCholInput.matrix" );

  info = MATRIX::cholesky( A, nRows );

  if ( info == 0 ) {
    return 0;
  } else if ( info < 0 ) {
    TRACE_ASSERT( NULL, "Invalid input to Cholesky solver" );
  }

  // The Cholesky factrorization failed, so copy the matrix
  // back and perform a Bunch-Kaufman factorization
  {
    IntArray                 pivotData( nRows );
    IntArray                 permutation( nRows );
    IntArray                 inversePermutation( nRows );

    // Workspaces for eigenvalue computations
    Real                     eigenvalues[ 2 ];
    Real                     eigenvectors[ 4 ];
    Real                     eigenSolverInput[ 4 ];

    Real                     Anorm1;

    int                      pivot;

    MATRIX                   L( nRows, nRows );

#ifdef USE_MKL
#ifdef SINGLE_PRECISION
    Real                     threshold = slamch( "e" );
#else
    Real                     threshold = dlamch( "e" );
#endif
#else
    // The interface for these functions seems to be different in the
    // standard LAPACKE implementation
#ifdef SINGLE_PRECISION
    Real                     threshold = LAPACKE_slamch( 'e' );
#else
    Real                     threshold = LAPACKE_dlamch( 'e' );
#endif
#endif

    MATRIX::copy( A, Acopy.data(), nRows, nRows );

    // Build a threshold based on the matrix norm
    Anorm1 = MATRIX::norm1( A, nRows, nRows );
    threshold /= 2.0;
    threshold = sqrt( threshold ) * Anorm1;

    MATRIX::LDL( A, nRows, pivotData );

    // Extract the "psychologically lower triangular" matrix L from this
    // factorization
    L.clear();
    MATRIX::buildLDLTriangle( pivotData, A, permutation, inversePermutation,
                              L.data() );

    // Step through the entries of the block diagonal matrix D
    for ( int block_start = 0; block_start < pivotData.size(); block_start++ )
    {
      pivot = pivotData[ block_start ];

      if ( pivot > 0 ) {
        // This is a 1x1 diagonal block
        Real                 d = MATRIX::access( A, nRows, nRows,
                                                 block_start, block_start );

        if ( d >= threshold ) {
          continue;
        }

        // This diagonal entry is too small, so we need to make a rank-1
        // modification to the original matrix to correct for this
        Real                 modification = max( threshold, abs( d ) );

        printf( "Bumping single entry %d from %f to %f\n",
                block_start, d, modification );
        
        modification = modification - d;

        // Append to the low-rank modification
        MATRIX::copyColumn( L.data(), U.data(), block_start, modificationSize,
                            nRows, nRows, nRows );

        // V should be zeroed initially, so this is okay
        //
        // Use addColumn here since we want to scale by modification
        MATRIX::addColumn( L.data(), V.data(), block_start, modificationSize,
                           nRows, nRows, nRows, nRows, modification );

        modificationSize += 1;
      }
      else if ( pivot < 0 ) {
        // This is a 2x2 diagonal block
        //
        // Copy the diagonal block and take it's eigen decomposition
        MATRIX::access( eigenSolverInput, 2, 2, 0, 0 )
          = MATRIX::access( A, nRows, nRows, block_start, block_start );
        MATRIX::access( eigenSolverInput, 2, 2, 1, 1 )
          = MATRIX::access( A, nRows, nRows, block_start + 1, block_start + 1 );
        MATRIX::access( eigenSolverInput, 2, 2, 1, 0 )
          = MATRIX::access( A, nRows, nRows, block_start + 1, block_start );
        MATRIX::access( eigenSolverInput, 2, 2, 0, 1 )
          = MATRIX::access( eigenSolverInput, 2, 2, 1, 0 );

#if 0
        // FIXME: debugging
        MATRIX::write( eigenSolverInput, 2, 2, "eigenInput.matrix" );
#endif

        MATRIX::eigensystem2x2( eigenSolverInput, eigenvalues, eigenvectors );

#if 0
        // FIXME: debugging
        MATRIX::write( eigenvalues, 2, 1, "eigenvalues.matrix" );
        MATRIX::write( eigenvectors, 2, 2, "eigenvectors.matrix" );
#endif

        if ( eigenvalues[ 0 ] >= threshold && eigenvalues[ 1 ] >= threshold ) {
          // Make sure we bump the column index by 2 for these cases
          block_start += 1;
          continue;
        }

        Real                 modification;

        // Bump the eigenvalues, then form the new diagonal block due to
        // this modification
        modification = max( threshold, abs( eigenvalues[ 0 ] ) );
        printf( "Bumping eigenvalue %d from %f to %f\n",
                block_start, eigenvalues[ 0 ], modification );
        eigenvalues[ 0 ] = modification - eigenvalues[ 0 ];
        modification = max( threshold, abs( eigenvalues[ 1 ] ) );
        printf( "Bumping eigenvalue %d from %f to %f\n",
                block_start + 1, eigenvalues[ 1 ], modification );
        eigenvalues[ 1 ] = modification - eigenvalues[ 1 ];

        TRACE_ASSERT( eigenvalues[ 0 ] >= 0.0 );
        TRACE_ASSERT( eigenvalues[ 1 ] >= 0.0 );

        printf( "Eigensolver input:\n" );
        printf( "[ %f, %f ]\n[ %f, %f ]\n",
                MATRIX::access( eigenSolverInput, 2, 2, 0, 0 ),
                MATRIX::access( eigenSolverInput, 2, 2, 0, 1 ),
                MATRIX::access( eigenSolverInput, 2, 2, 1, 0 ),
                MATRIX::access( eigenSolverInput, 2, 2, 1, 1 ) );

        MATRIX::access( eigenvectors, 2, 2, 0, 0 ) *= sqrt( eigenvalues[ 0 ] );
        MATRIX::access( eigenvectors, 2, 2, 1, 0 ) *= sqrt( eigenvalues[ 0 ] );
        MATRIX::access( eigenvectors, 2, 2, 0, 1 ) *= sqrt( eigenvalues[ 1 ] );
        MATRIX::access( eigenvectors, 2, 2, 1, 1 ) *= sqrt( eigenvalues[ 1 ] );

        // Form the symmetric modification
        MATRIX::syrk( eigenvectors, eigenSolverInput, 2, 2 );
        // Need to adjust upper triangular part, since syrk doesn't
        // change this
        MATRIX::access( eigenSolverInput, 2, 2, 0, 1 )
          = MATRIX::access( eigenSolverInput, 2, 2, 1, 0 );

        printf( "After reconstruction:\n" );
        printf( "[ %f, %f ]\n[ %f, %f ]\n",
                MATRIX::access( eigenSolverInput, 2, 2, 0, 0 ),
                MATRIX::access( eigenSolverInput, 2, 2, 0, 1 ),
                MATRIX::access( eigenSolverInput, 2, 2, 1, 0 ),
                MATRIX::access( eigenSolverInput, 2, 2, 1, 1 ) );

        // Append the low-rank modification
        MATRIX::copyColumn( L.data(), U.data(), block_start, modificationSize,
                            nRows, nRows, nRows );
        MATRIX::copyColumn( L.data(), U.data(), block_start + 1,
                            modificationSize + 1, nRows, nRows, nRows );

        // Multiply the columns by the square modification matrix to form V
        MATRIX::gemm( U.data() + modificationSize, /* Align with column */
                      eigenSolverInput,
                      V.data() + modificationSize, /* Align with column */
                      nRows, 2, 2, 2, /* Matrix sizes */
                      false, false, /* No transposition */
                      1.0, 0.0, /* Overwrite result */
                      nRows, 2, nRows /* Leading dimensions */ );

        modificationSize += 2;
        block_start += 1;
      }
      else {
        TRACE_ASSERT( NULL, "Should never get here" );
      }
    }

    // Restore the original matrix, add the low rank modification,
    // and take its Cholesky factor.
    MATRIX::copy( A, Acopy.data(), nRows, nRows );
#if 0
    MATRIX::write( A, nRows, nRows, "preAdd.matrix" );
#endif
    MATRIX::gemm( U.data(), V.data(), A,
                  nRows, modificationSize, /* Size of U's contents */
                  nRows, modificationSize, /* Size of V's contents */
                  false, true, /* Transpose V */
                  1.0, 1.0, /* Add to A */
                  nRows, nRows /* Leading dimensions for U and V are nRows */ );

    MATRIX::write( A, nRows, nRows, "preCholesky.matrix" );

    info = MATRIX::cholesky( A, nRows );

    if ( true || info != 0 ) {
      U.write( "Utest.matrix" );
      V.write( "Vtest.matrix" );
    }
    TRACE_ASSERT( info == 0, "Modified Cholesky failed" );
  }

  abort();

  return modificationSize;
}

//////////////////////////////////////////////////////////////////////
// Given a triangular matrix, and another input matrix, solve
// the associated triangular system.
//
// By default, solves L * X = B, where L is non-unit lower-triangular
//////////////////////////////////////////////////////////////////////
void MATRIX::triangularSolve( const Real *L, Real *B,
                              int rowsB, int colsB,
                              bool leftSide,
                              bool lower,
                              bool transpose,
                              bool unitDiag,
                              Real alpha,
                              int ldaB )
{
  CBLAS_SIDE side = leftSide ? CblasLeft : CblasRight;
  CBLAS_UPLO uplo = lower ? CblasLower : CblasUpper;
  CBLAS_TRANSPOSE transL = transpose ? CblasTrans : CblasNoTrans;
  CBLAS_DIAG diag = unitDiag ? CblasUnit : CblasNonUnit;

  if ( ldaB < 0 ) {
    ldaB = colsB;
  }

#ifdef DO_FLOP_COUNT
  MATRIX_FLOP_COUNT += leftSide ? rowsB * rowsB * colsB
                                : rowsB * colsB * colsB;
#endif

#ifdef SINGLE_PRECISION
  cblas_strsm(
      CblasRowMajor,
      side,
      uplo,
      transL,
      diag,
      rowsB, colsB,
      alpha,
      L,
      leftSide ? rowsB : colsB,
      B,
      ldaB);
#else
  cblas_dtrsm(
      CblasRowMajor,
      side,
      uplo,
      transL,
      diag,
      rowsB, colsB,
      alpha,
      L,
      leftSide ? rowsB : colsB,
      B,
      ldaB);
#endif
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but for a vector
//////////////////////////////////////////////////////////////////////
void MATRIX::triangularSolve( const Real *L, Real *b,
                              int N,
                              bool lower,
                              bool transpose,
                              bool unitDiag,
                              Real alpha,
                              int ldaB )
{
  CBLAS_UPLO uplo = lower ? CblasLower : CblasUpper;
  CBLAS_TRANSPOSE transL = transpose ? CblasTrans : CblasNoTrans;
  CBLAS_DIAG diag = unitDiag ? CblasUnit : CblasNonUnit;

  if ( ldaB < 1 ) {
    ldaB = 1;
  }

#ifdef DO_FLOP_COUNT
  MATRIX_FLOP_COUNT += N * N;
#endif

#ifdef SINGLE_PRECISION
  cblas_strsv(
      CblasRowMajor,
      uplo,
      transL,
      diag,
      N,
      L, N,
      b, ldaB);
#else
  cblas_dtrsv(
      CblasRowMajor,
      uplo,
      transL,
      diag,
      N,
      L, N,
      b, ldaB);
#endif
}

//////////////////////////////////////////////////////////////////////
// 2x2 eigensolver
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem2x2( const Real *matrix, Real *eigenvalues,
                             Real *eigenvectors )
{
  Real A = access( matrix, 2, 2, 0, 0 );
  Real B = access( matrix, 2, 2, 0, 1 );
  Real C = access( matrix, 2, 2, 1, 0 );
  Real D = access( matrix, 2, 2, 1, 1 );

  if (B * C <= 0.1e-20) {
    eigenvalues[0] = A; 
    access( eigenvectors, 2, 2, 0, 0 ) = 1; 
    access( eigenvectors, 2, 2, 1, 0 ) = 0;
    eigenvalues[1] = D; 
    access( eigenvectors, 2, 2, 0, 1 ) = 0; 
    access( eigenvectors, 2, 2, 1, 1 ) = 1;
    return;
  }

  Real tr = A + D;
  Real det = A * D - B * C;
  Real S = sqrt(tr * tr * 0.25 - det );
  eigenvalues[0] = tr * 0.5 + S;
  eigenvalues[1] = tr * 0.5 - S;

  Real temp = (A - D) * (A - D) * 0.25 + B * C;
  Real SS = (temp < 0.0) ? 0.0 : sqrt(temp);
  if (A - D < 0.0) {
    access( eigenvectors, 2, 2, 0, 0 ) = C;
    access( eigenvectors, 2, 2, 1, 0 ) = -(A - D) * 0.5 + SS;
    access( eigenvectors, 2, 2, 0, 1 ) = (A - D) * 0.5 - SS;
    access( eigenvectors, 2, 2, 1, 1 ) = B;
  } 
  else {
    access( eigenvectors, 2, 2, 0, 1 ) = C;
    access( eigenvectors, 2, 2, 1, 1 ) = -(A - D) * 0.5 - SS;
    access( eigenvectors, 2, 2, 0, 0 ) = (A - D) * 0.5 + SS;
    access( eigenvectors, 2, 2, 1, 0 ) = B;
  }

  Real n1 = sqrt(access( eigenvectors, 2, 2, 0,0) * access( eigenvectors, 2, 2, 0,0) +
                 access( eigenvectors, 2, 2, 1,0) * access( eigenvectors, 2, 2, 1,0));
  Real inv = 1.0 / n1;
  access( eigenvectors, 2, 2, 0, 0) *= inv; 
  access( eigenvectors, 2, 2, 1, 0) *= inv;
  Real n2 = sqrt(access( eigenvectors, 2, 2, 0,1) * access( eigenvectors, 2, 2, 0,1) +
                 access( eigenvectors, 2, 2, 1,1) * access( eigenvectors, 2, 2, 1,1));
  inv = 1.0 / n2;
  access( eigenvectors, 2, 2, 0, 1) *= inv; 
  access( eigenvectors, 2, 2, 1, 1) *= inv;
}

//////////////////////////////////////////////////////////////////////
// 3x3 low precision eigensolver
//////////////////////////////////////////////////////////////////////
void MATRIX::eigensystem3x3( Real *A, Real *eigenvalues, Real *eigenvectors )
{
	register int i, j, k, l;

  //float TOL = 1e-8f;
  //int MAX_SWEEPS = 500;
  float TOL = 1e-3f;
  int MAX_SWEEPS = 5;
  unsigned int n = 3;

  float a[9];
  float d[3];
  float v[9];
  i = 0;
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++, i++)
    {
      a[i] = access( A, 3, 3, x, y);
      v[i] = (x == y) ? 1.0 : 0.0;
    }
  
  
	float onorm, dnorm;
	float b, dma, q, t, c, s;
	float atemp, vtemp, dtemp;

	// Set v to the identity matrix, set the vector d to contain the
	// diagonal elements of the matrix a
	d[0] = a[0];
	d[1] = a[4];
	d[2] = a[8];

	for (l = 1; l <= MAX_SWEEPS; l++)
	{
		// Set dnorm to be the maximum norm of the diagonal elements, set
		// onorm to the maximum norm of the off-diagonal elements
		
		dnorm = (float)fabs(d[0]) + (float)fabs(d[1]) + (float)fabs(d[2]);
		onorm = (float)fabs(a[1]) + (float)fabs(a[2]) + (float)fabs(a[5]);
		// Normal end point of this algorithm.
		if((onorm/dnorm) <= TOL)
			goto Exit_now;

		for (j = 1; j < static_cast<int>(n); j++)
		{
			for (i = 0; i <= j - 1; i++)
			{

				b = a[n*i+j];
				if(fabs(b) > 0.0f)
				{
					dma = d[j] - d[i];
					if((fabs(dma) + fabs(b)) <= fabs(dma))
						t = b / dma;
					else
					{
						q = 0.5f * dma / b;
						t = 1.0f/((float)fabs(q) + (float)sqrt(1.0f+q*q));
						if (q < 0.0)
							t = -t;
					}

					c = 1.0f/(float)sqrt(t*t + 1.0f);
					s = t * c;
					a[n*i+j] = 0.0f;

					for (k = 0; k <= i-1; k++)
					{
						atemp = c * a[n*k+i] - s * a[n*k+j];
						a[n*k+j] = s * a[n*k+i] + c * a[n*k+j];
						a[n*k+i] = atemp;
					}

					for (k = i+1; k <= j-1; k++)
					{
						atemp = c * a[n*i+k] - s * a[n*k+j];
						a[n*k+j] = s * a[n*i+k] + c * a[n*k+j];
						a[n*i+k] = atemp;
					}

					for (k = j+1; k < static_cast<int>(n); k++)
					{
						atemp = c * a[n*i+k] - s * a[n*j+k];
						a[n*j+k] = s * a[n*i+k] + c * a[n*j+k];
						a[n*i+k] = atemp;
					}

					for (k = 0; k < static_cast<int>(n); k++)
					{
						vtemp = c * v[n*k+i] - s * v[n*k+j];
						v[n*k+j] = s * v[n*k+i] + c * v[n*k+j];
						v[n*k+i] = vtemp;
					}

					dtemp = c*c*d[i] + s*s*d[j] - 2.0f*c*s*b;
					d[j] = s*s*d[i] + c*c*d[j] + 2.0f*c*s*b;
					d[i] = dtemp;
				} /* end if */
			} /* end for i */
		} /* end for j */
	} /* end for l */

Exit_now:
  for (int x = 0; x < 3; x++)
    eigenvalues[x] = d[x];

  i = 0;
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++, i++)
      access( eigenvectors, 3, 3, x, y) = v[i];

	return;
}

//////////////////////////////////////////////////////////////////////
// Scales a matrix.  A *= alpha
//////////////////////////////////////////////////////////////////////
void MATRIX::scale( Real *A, int rows, int cols, Real alpha )
{
#ifdef DO_FLOP_COUNT
  MATRIX_FLOP_COUNT += rows * cols;
#endif

#ifdef SINGLE_PRECISION
  cblas_sscal(cols * rows, alpha, A, 1)
#else
  cblas_dscal(cols * rows, alpha, A, 1);
#endif
}

//////////////////////////////////////////////////////////////////////
// Scaling with a leading dimension
//////////////////////////////////////////////////////////////////////
void MATRIX::scaleMatrix( Real *A, int rows, int cols, int lda, Real alpha )
{
  lda = max( lda, cols );

#ifdef DO_FLOP_COUNT
  MATRIX_FLOP_COUNT += rows * cols;
#endif

  for ( int i = 0; i < rows; i++ )
  {
#ifdef SINGLE_PRECISION
    cblas_sscal( cols, alpha, A + i * lda, 1 );
#else
    cblas_dscal( cols, alpha, A + i * lda, 1 );
#endif
  }
}

//////////////////////////////////////////////////////////////////////
// Check to see if any entries in the matrix are NaNs
//////////////////////////////////////////////////////////////////////
bool MATRIX::is_nan( const Real *A, int rows, int cols, int lda )
{
  if ( lda < 0 )
  {
    lda = cols;
  }

  for ( int i = 0; i < rows; i++ )
  for ( int j = 0; j < cols; j++ )
  {
    if ( std::isnan( A[ i * lda + j ] ) || std::isinf( A[ i * lda + j ] ) )
    {
      return true;
    }
  }

  return false;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void MATRIX::write( const Real *A, int rows, int cols, const char *fileName )
{
  MATRIX                     tmp( rows, cols, A );

  tmp.write( fileName );
}

//////////////////////////////////////////////////////////////////////
// Matrix 1 norm (ie. maximum absolute column sum)
//////////////////////////////////////////////////////////////////////
Real MATRIX::norm1( const Real *A, int rows, int cols, int lda )
{
  Real                       maxSum = 0.0;
  Real                       colSum;

  lda = max( cols, lda );

  for ( int col_idx = 0; col_idx < cols; col_idx++ ) {
    colSum = VECTOR::absSum( rows, A + col_idx /* Start of column */, lda );

    maxSum = max( maxSum, colSum );
  }

  return maxSum;
}

#if 0
void XERBLA (const char * Name, const int * Num, const int Len) {
  printf("nUser xerbla is called :%s:%dn",Name,*Num); fflush(NULL);
  return;
}

void cblas_xerbla (const char * Name, const int Num) {
  printf("nUser cblas_xerbla is called :%s:%dn",Name,Num); fflush(NULL);
  return;
}
#endif
