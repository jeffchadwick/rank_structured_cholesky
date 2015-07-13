// MATRIX.h: interface for the MATRIX class.
//
//////////////////////////////////////////////////////////////////////

#ifndef MATRIX_H
#define MATRIX_H

#include <SETTINGS.h>
#include <TYPES.h>

#include <map>
#include <iostream>
#include <cstdio>
#include "VECTOR.h"
#include "VEC3.h"
#include "MATRIX3.h"

#include <util/trace.h>

#include <memory.h>

#define DO_FLOP_COUNT 1

#if 0
#define FLOP_COUNT_START { \
  MATRIX_FLOP_COUNT = 0; \
}

#define FLOP_COUNT_END { \
  MATRIX_FLOP_COUNT \
}
#endif

using namespace std;

//////////////////////////////////////////////////////////////////////
// An arbitrary dimension matrix class
//////////////////////////////////////////////////////////////////////
class MATRIX {

public:
#ifdef DO_FLOP_COUNT
  static size_t              MATRIX_FLOP_COUNT;
#endif

  MATRIX();
  MATRIX(int rows, int cols);
  MATRIX(int rows, int cols, const Real* data);
  MATRIX(int rows, int cols, int lda, const Real *data);
  MATRIX(const char* filename);
  MATRIX(const MATRIX& m);
  // matrix with "vec" along the diagonal
  MATRIX(VECTOR& vec);

  // make a copy of a 3x3 matrix
  MATRIX(MATRIX3& matrix3);
  virtual ~MATRIX();

  inline Real& operator()(int row, int col) {
    return _matrix[row * _cols + col];
  };

  inline Real operator()(int row, int col) const {
    return _matrix[row * _cols + col];
  };
 
  int& rows() { return _rows; };
  int& cols() { return _cols; };

  int rows() const { return _rows; };
  int cols() const { return _cols; };

  // return a pointer to the beginning of a row
  Real* row(int index) { return &_matrix[index * _cols]; };

  // wipe the whole matrix
  void clear();

  // write the matrix to a binary file
  // everything is always written as a double
  void write(const char* filename);

  // read from a binary file
  // everything is always read in as a double, then
  // converted if necessary
  void read(const char* filename);

  // resize the matrix and wipe to zero
  void resizeAndWipe(int rows, int cols);

  // overload operators
  MATRIX& operator=(const MATRIX m);
  MATRIX& operator-=(const MATRIX& m);
  MATRIX& operator+=(const MATRIX& m);
  MATRIX& operator*=(const Real& alpha);

  // return the matrix diagonal
  VECTOR diagonal();

  // return the transpose of the current matrix
  MATRIX transpose();

  // raw data pointer
  Real* data() { return _matrix; };

	// 
	Real sum()
	{
		Real sum = 0.0;
		for( int i = 0; i < _rows*_cols; i++ )
		{
			sum += _matrix[i];
		}
		return sum;
	}

  // stomp the current matrix with the given matrix 
  // starting at row number "row". It is your responsibility
  // to ensure that you don't fall off the end of this matrix.
  void setSubmatrix(MATRIX& matrix, int row);

  // Stomp the current matrix with the given vector,
  // starting at (row,col). So the first elemenet of the vector will
  // be copied into (row,col), the next into (row+1, col), etc.
  void setVector( VECTOR& vector, int row, int col );

  void copyRowFrom( MATRIX& src, int srcRow, int row );
  void copyRowFrom( VECTOR& src, int srcRow, int row );
  void copyRowFrom( VEC3F& src, int srcRow, int row );

  // BLAS axpy operation: B += alpha * A, where B is this matrix
  //
  // Note that axpy actually applies to vectors, but in this
  // case we can just treat the matrix as a vector and multiply
  // all its elements by alpha
  void axpy(Real alpha, MATRIX& A);

  // same as axpy above, but matrix contents are stomped as well
  void clearingAxpy(Real alpha, MATRIX& A);

  // BLAS gemm operation: C += alpha * A * B
  // where C is this matrix
  void gemm(Real alpha, MATRIX& A, MATRIX& B);
  void gemm(MATRIX& A, MATRIX& B) { gemm(1.0f, A, B); };

  // same as gemm above, but matrix contents are stomped as well
  void clearingGemm(Real alpha, MATRIX& A, MATRIX& B,
										bool transposeA = false, bool tranposeB = false);
  void clearingGemm(MATRIX& A, MATRIX& B,
										bool transposeA = false, bool transposeB = false)
	{
		clearingGemm(1.0f, A, B, transposeA, transposeB);
	}

  // BLAS gemv operation: y = alpha * A * x
  // where A is this matrix
  VECTOR gemv(VEC3F& x);
  VECTOR gemv(Real alpha, VEC3F& x);
  //Untested -- don't uncomment until a test case comes up
  //VECTOR gemv(VECTOR& x);
  
  // solve the linear system Ax = b, return x in the passed in b
  void solve(VECTOR& b);
	void solveSPD(VECTOR& b); // Same thing for symmetric positive-definite systems

  // Note: This version involves taking tranposes (eg. dynamic
  // allocation).
  void solveSPD(MATRIX &B);

	// Inverts this matrix in place
	void invert();

	// Returns the inverse of this matrix in A
	void inverse( MATRIX &A );

  // multiply in place
  // this * x = y
  // Do *NOT* let x == y!
  void multiplyInplace(VECTOR& x, VECTOR& y,
											 Real alpha = 1.0, bool add = false);
  void multiplyInplace(Real *x, Real *y,
											 Real alpha = 1.0, bool add = false);

  void subMatrixMultiplyInplace( VECTOR& x, VECTOR& prod, int subRows,
																 int subCols, bool transpose );

  // Assumes matrix is symmetric
  void uppertriMultiplyInplace( VECTOR& x, VECTOR& prod );

  // solve for eigenvalues
  void eigensystem(VECTOR& eigenvalues, MATRIX& eigenvectors);
 
  // copy this matrix to MATRIX3 type
  void copiesInto(MATRIX3& matrix3);

	// Returns the Frobenius-norm of the difference between this and B
	Real differenceFrobeniusSq( MATRIX& B );

  static MATRIX ident( int nRows )
  {
    MATRIX                   result( nRows, nRows );

    for ( int i = 0; i < nRows; i++ )
    {
      result( i, i ) = 1.0;
    }

    return result;
  }

  //////////////////////////////////////////////////////////////////////
  // Static routines for handling array-based matrices.
  // Use with care - nothing here does any bound checking.
  // We assume row-major format in all cases.
  //////////////////////////////////////////////////////////////////////

  // Zero out a matrix
  static inline void clear( Real *A, int rows, int cols )
  {
    memset( (void *)A, 0, rows * cols * sizeof( Real ) );
  }

  // Zero a matrix with a given leading dimension
  static inline void clear( Real *A, int rows, int cols, int lda )
  {
    for ( int row_idx = 0; row_idx < rows; row_idx++ )
    {
      clear( A + row_idx * lda, 1, cols );
    }
  }

  static inline void ident( Real *A, int rows )
  {
    clear( A, rows, rows );

    for ( int row_idx = 0; row_idx < rows; row_idx++ )
    {
      access( A, rows, rows, row_idx, row_idx ) = 1.0;
    }
  }

  // Generate a diagonal matrix
  static inline void diagonal( Real *A, Real *D, int rows )
  {
    for ( int i = 0; i < rows; i++ )
    {
      access( A, rows, rows, i, i ) = D[ i ];
    }
  }

  // Copy one matrix in to another (A <- B)
  static void copy( Real *A, const Real *B, int rows, int cols );

  // Copy a 3x3 matrix
  static void copy( Real *A, const MATRIX3 &B );

  // Copy operation in which the leading dimension of the matrix
  // changes
  static void copy( Real *A, const Real *B, int rows, int cols,
                    int ldaA, int ldaB );

  // Copies a row (or part of it) from one matrix to another
  static void copyRow( const Real *A, Real *B,
                       int rowA, int rowB,
                       int nColsA, int nColsB,
                       int nCopyCols = -1 );

  // Adds one row in A to another in B
  static void addRow( const Real *A, Real *B,
                      int rowA, int rowB,
                      int nColsA, int nColsB,
                      int nCopyCols = -1,
                      Real alpha = 1.0 );

  // Copies a set of rows r0, r1, r2, ... from A to rows 0, 1, 2, ... of B
  static void copyRows( const Real *A, Real *B,
                        const IntArray &rowsA, int nCols );

  // Copies a subset of the rows r0, r1, r2, specified by the given
  // index range.
  static void copyRows( const Real *A, Real *B,
                        const IntArray &rowsA, int nCols,
                        const IndexRange &rowRange,
                        // Optional offset for row indices
                        int offset = 0 );

  // Adds rows r0, r1, r2, ... of A to rows 0, 1, 2, ... of B
  static void copyAddRows( const Real *A, Real *B,
                           const IntArray &rowsA, int nCols,
                           const IndexRange &rowRange,
                           // Optional offset for row indices
                           int offset,
                           Real alpha = 1.0 );

  static void copyAddRows( const Real *A, Real *B,
                           const IntArray &rowsA, int nCols,
                           Real alpha )
  {
    copyAddRows( A, B, rowsA, nCols, IndexRange( 0, rowsA.size() - 1 ),
                 0, alpha );
  }

  // Copies rows 0, 1, 2, of A in to r0, r1, r2, of B
  static void scatterRows( const Real *A, Real *B,
                           const IntArray &rowsB, int nCols );

  // Scatters rows 0, 1, ... of A to a subset of r0, r1, r2 of B specified
  // with the given index range.
  static void scatterRows( const Real *A, Real *B,
                           const IntArray &rowsB, int nCols,
                           const IndexRange &rowRange,
                           // Optional offset for row indices
                           int offset );

  // Scatters and adds rows 0, 1, ... of A to a subset of r0, r1, r2 of B
  // specified with the given index range
  static void scatterAddRows( const Real *A, Real *B,
                              const IntArray &rowsB, int nCols,
                              const IndexRange &rowRange,
                              // Optional offset for row indices
                              int offset,
                              Real alpha = 1.0 );

  // In place row swap
  static void swapRows( Real *A, int nRows, int nCols, int row1, int row2 );

  // Copy a column (or part of it) from one matrix to another
  static void copyColumn( const Real *A, Real *B,
                          int colA, int colB,
                          int nRowsA,
                          int nColsA, int nColsB,
                          int nCopyRows = -1 );

  // Adds a column from A to a column from B
  static void addColumn( const Real *A, Real *B,
                         int colA, int colB,
                         int nRowsA,
                         int nColsA, int nColsB,
                         int nCopyRows = -1,
                         Real alpha = 1.0 );

  // Copies a set of columns c0, c1, c2, ... from A to columns 0, 1, 2, ... of B
  static void copyColumns( const Real *A, Real *B,
                           const IntArray &colsA, int nRows,
                           int nColsA, int nColsB );

  // Copies a subset of the columns c0, c1, c2, ..., specified by the
  // given index range.
  static void copyColumns( const Real *A, Real *B,
                           const IntArray &colsA, int nRows,
                           int nColsA, int nColsB,
                           const IndexRange &columnRange,
                           // Optional offset for column indices
                           int offset = 0 );

  // Copies rows 0, 1, 2, of A in to c0, c1, c2, ... of B
  static void scatterColumns( const Real *A, Real *B,
                              const IntArray &colsB, int nRows,
                              int nColsA, int nColsB );

  // Scatters columns 0, 1, ... of A to a subset of c0, c1, ... of B specified
  // with the given index range.
  static void scatterColumns( const Real *A, Real *B,
                              const IntArray &colsB, int nRows,
                              int nColsA, int nColsB,
                              const IndexRange &columnRange,
                              // Optional offset for column indices
                              int offset );

  // Scatters columns 0, 1, ... of A to a subset of c0, c1, ... of B specified
  // with the given index range.
  static void scatterAddColumns( const Real *A, Real *B,
                                 const IntArray &colsB, int nRows,
                                 int nColsA, int nColsB,
                                 const IndexRange &columnRange,
                                 // Optional offset for column indices
                                 int offset,
                                 Real alpha = 1.0 );

  // For a symmetric matrix in which the given rows/columns are filled
  // In place column swap
  static void swapColumns( Real *A, int nRows, int nCols, int col1, int col2 );

  // in, scatter and add the entries to the given full matrix.
  //
  // Only updates the lower triangular part of B.
  static void scatterAddSymmetricMatrix( const Real *A, Real *B,
                                         const IntArray &colsB, int nRows,
                                         const IndexRange &columnRange,
                                         // Optional offset
                                         int offset,
                                         Real alpha = 1.0 );

  // Matrix vector multiply for a matrix in compressed column form
  // (ie. only a subset of the columns are stored and the others are
  // assumed to be 0).
  static void compressedColumnMult( const Real *A, const Real *b, Real *c,
                                    const IntArray &colsA,
                                    int nRows, int nCols,
                                    bool transpose = false,
                                    Real alpha = 1.0, Real beta = 0.0 );

#define TEST_ACCESS 1
//#undef TEST_ACCESS
  // Access element from a matrix
  static inline Real &access( Real *A, int rows, int cols, int i, int j )
  {
#ifdef TEST_ACCESS
    TRACE_ASSERT( i >= 0 && i < rows && j >= 0 && j < cols );
#endif
    return A[ i * cols + j ];
  }
  static inline Real access( const Real *A, int rows, int cols, int i, int j )
  {
#ifdef TEST_ACCESS
    TRACE_ASSERT( i >= 0 && i < rows && j >= 0 && j < cols );
#endif
    return A[ i * cols + j ];
  }

  // Adds a scalar to the diagonal of a square matrix
  //
  // eg. A(i,i) += alpha, for i = 1:nRows
  static void addToDiagonal( Real *A, int nRows, Real alpha );

  // Add one matrix to another (A := A + alpha * B)
  static void axpy( Real *A, const Real *B, int rows, int cols,
                    Real alpha = 1.0,
                    int incA = 1, int incB = 1 );

  // Matrix-matrix multiplication (C := beta * C + alpha * A * B)
  static void gemm( const Real *A, const Real *B, Real *C,
                    int rowsA, int colsA, int rowsB, int colsB,
                    bool transposeA, bool transposeB,
                    Real alpha = 1.0, Real beta = 0.0,
                    int ldaA = -1, int ldaB = -1, int ldaC = -1 );

  // Matrix-vector multiplication (C := beta * C + alpha * A * b)
  static void gemv( const Real *A, const Real *b, Real *c,
                    int rowsA, int colsA,
                    bool transposeA,
                    Real alpha = 1.0, Real beta = 0.0 );

  // Symmetrix matrix-matrix update (C := beta * C + alpha * A * A')
  // Updates the lower-triangular part of C
  //
  // If trans == true then C := beta * C + alpha * A' * A
  static void syrk( const Real *A, Real *C,
                    int rowsC, int k, /* A is either n x k or k x n */
                    bool transpose = false,
                    Real alpha = 1.0, Real beta = 0.0,
                    int ldaA = -1, int ldaC = -1 );

  // Get transpose (A = B^T)
  static void transpose( Real *A, const Real *B, int rows, int cols );

  static void transposeBLAS( Real *A, const Real *B, int rows, int cols,
                             int ldaA = -1, int ldaB = -1 );

  // In place transpose for square matrices
  static void transpose( Real *A, int rows );

  // Compute eigenvalues/vectors.  We require everything here,
  // including workspaces, etc.
  // workspace should have size (7 * rows)
  // vectorWorkspace should have size (rows * rows)
  static void eigensystem( Real *A, int rows,
                           Real *eigenvalues, Real *eigenvectors,
                           Real *workspace, Real *vectorWorkspace );

  // Performs QR factorization of the given system.  Other calls
  // will be needed following this to work the factorization,
  // since this is just a wrapper for the basic LAPACK call.
  //
  // No pivoting is used here
  static int qr( Real *A, int nRows, int nCols,
                 Real *extraData, Real *workspace,
                 int workSize );

  // Performs QR factorization with column pivoting on the given
  // system.  Other calls will be needed following this to work with
  // the factorization, since this is just a wrapper for the basic
  // LAPACK call.
  static int qrPivot( Real *A, int *pivots, int nRows, int nCols,
                      Real *extraData, Real *workspace,
                      int workSize );

  // Extracts the leading columns of a QR factor computed
  // via the qr function.  Parameters should be the same.
  static int extractQRfactor( Real *A, int nRows, int nCols,
                              Real *extraData, Real *workspace,
                              int workSize );

  // LU factorization - requires the matrix, and a workspace for
  // pivoting information
  static int LU( Real *A, int *pivotData, int rows );

  // Solve via an LU factorization
  static int LUsolve( const Real *A, const int *pivotData,
                      Real *B, int rowsB, int colsB,
                      bool transposeA = false,
                      int ldaA = -1, int ldaB = -1 );

  // Computes the Cholesky factor of the matrix stored in A.  A is
  // overwritten.
  // Computes L such that A = L * L'
  static int cholesky( Real *A, int rows );

  // Computes a block LDL^T factorization using the Bunch-Kaufman method.
  //
  // Doesn't actual extract the factors from the matrix.... yet
  static int LDL( Real *A, int rows, IntArray &pivotData );

#if 0
  // Runs the LAPACK converseion routine for a LDL^T factorization
  static int convertLDL( Real *A, int rows, const IntArray &permutation,
                         FloatArray &workspace );
#endif

  // Routines for extracting LDL^T factor data.
  //
  // Builds a full permutation matrix (possibly overkill, but let's
  // do it anyways)
  static void buildLDLPermutation( const IntArray &pivotData,
                                   IntArray &permutation,
                                   Real *P = NULL );

  // Extracts a "psychologically" lower triangular matrix from a LDL^T
  // factorization; that is, a matrix which is a product of lower triangular
  // and permutation matrices.
  //
  // Also, builds a permutation vector and its inverse along the way.
  static void buildLDLTriangle( const IntArray &pivotData,
                                const Real *A,
                                IntArray &permutation,
                                IntArray &inversePermutation,
                                Real *L );

  // Extracts the lower triangular system from the LDL^T factorization
  static void convertLDLTriangle( const IntArray &pivotData,
                                  const Real *A,
                                  Real *L );

#if 0
  // A hacky and inefficient version of the above, which just does
  // a bunch of matrix multiplication to form the psychologically lower
  // triangular matrix.
  static void convertLDLTHack( const IntArray &pivotData,
                               const Real *A,
                               MATRIX &L );
#endif

  // Extracts the block diagonal from a LDL^T factorization
  static void extractLDLDiagonal( const IntArray &pivotData,
                                  const Real *A,
                                  Real *D );

  // Given a pivot array returned by LDL, apply the permutation implied
  // by this array.  The permutation is applied to the matrix in place.
  //
  // Options provided:
  //    - apply on left or right
  //    - apply transposed or untransposed permutation
  static void applyLDLPermutation( const IntArray &pivotData,
                                   Real *A, int nRows, int nCols,
                                   bool left = true, bool transpose = false );

  // Modified Cholesky factorization based on a Bunch-Kaufman factorization.
  //
  // Overwrites A with the modified Cholesky factor and fills in U and V
  // with the low-rank update necessary to make A positive definite.
  //
  // Returns the number of columns in the low-rank modification
  static int modifiedCholesky( Real *A, MATRIX &U, MATRIX &V,
                               int nRows );

  // Given a triangular matrix, and another input matrix, solve
  // the associated triangular system.
  //
  // By default, solves L * X = B, where L is non-unit lower-triangular
  static void triangularSolve( const Real *L, Real *B,
                               int rowsB, int colsB,
                               bool leftSide = true,
                               bool lower = true,
                               bool transpose = false,
                               bool unitDiag = false,
                               Real alpha = 1.0,
                               int ldaB = -1 );

  // Same as the above, but for a vector
  static void triangularSolve( const Real *L, Real *b,
                               int N,
                               bool lower = true,
                               bool transpose = false,
                               bool unitDiag = false,
                               Real alpha = 1.0,
                               int ldaB = -1 );

  // 2x2 eigensolver
  static void eigensystem2x2( const Real *A, Real *eigenvalues,
                              Real *eigenvectors );

  // 3x3 low precision eigensolver
  static void eigensystem3x3( Real *A, Real *eigenvalues, Real *eigenvectors );

  // Scales a matrix.  A *= alpha
  static void scale( Real *A, int rows, int cols, Real alpha );

  // Scaling with a leading dimension
  static void scaleMatrix( Real *A, int rows, int cols, int lda, Real alpha );

  // Check to see if any entries in the matrix are NaNs
  static bool is_nan( const Real *A, int rows, int cols, int lda = -1 );

  static void write( const Real *A, int rows, int cols, const char *fileName );

  // Matrix 1 norm (ie. maximum absolute column sum)
  static Real norm1( const Real *A, int rows, int cols, int lda = 0 );

  static long int nnz( const Real *A, int rows, int cols ) {
    long int                 nnz = 0;

    for ( int i = 0; i < rows; i++ )
    for ( int j = 0; j < cols; j++ ) {
      if ( access( A, rows, cols, i, j ) != 0.0 ) {
        nnz += 1;
      }
    }

    return nnz;
  }

protected:
  int _rows;
  int _cols;

  Real* _matrix;
};

// overloaded operators
VECTOR operator*(MATRIX& A, VECTOR& x);
MATRIX operator*(MATRIX& A, Real alpha);
MATRIX operator*(MATRIX& A, MATRIX& B);
ostream& operator<<(ostream &out, MATRIX& matrix);

// multiply by the transpose of A
VECTOR operator^(MATRIX& A, VECTOR& x);
MATRIX operator^(MATRIX& A, MATRIX& B);

#endif
