//////////////////////////////////////////////////////////////////////
// CSparse_Interface.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef CSPARSE_INTERFACE_H
#define CSPARSE_INTERFACE_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/SPARSE_MATRIX.h>
#include <linearalgebra/VECTOR.h>

#include <SETTINGS.h>
#include <TYPES.h>

#include <cs.h>

//////////////////////////////////////////////////////////////////////
// CSparse_Interface class
//
// Static routines for interfacing with the CXSparse library
//////////////////////////////////////////////////////////////////////
class CSparse_Interface {
	public:
    enum MatrixStorage {
      FULL = 0,
      LOWER,
      UPPER
    };

    // Builds a CSparse formatted sparse matrix out of a column sparse
    // matrix in our format.  Optionally permutes the matrix.  If a
    // permutation is provided then permute the matrix.  If permutation
    // is smaller than the number of columns in the matrix then the
    // remainder of the matrix will be filled in with diagonals.
    static void CopySparseSquareSymbolicMatrix(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  cs * &Acopy,
                                  const IntArray *permutation = NULL,
                                  const IntArray *inversePermutation = NULL,
                                  MatrixStorage storageType = FULL,
                                  bool storeDiagonal = true );

    // Copies a sparse symmetric symbolic matrix
    static void CopySparseSymmetricSymbolicMatrix(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  cs * &Acopy,
                                  const IntArray *permutation = NULL,
                                  const IntArray *inversePermutation = NULL );

    // Transpose routine
    static void TransposeSparseSymbolicMatrix( const cs *A, cs * &Atranspose );

    // Symbolic matrix add
    static void AddSparseSymbolicMatrices( const cs *A, const cs *B, cs * &C );

    // Constructs the elimination tree for the given matrix
    static void ConstructEliminationTree( const cs *A, IntArray &parent );

    // Builds an elimination tree post order
    static void ConstructPostOrder( const IntArray &parent,
                                    IntArray &postOrder );

    // Given an elimination tree and its post ordering, determine column
    // counts in the Cholesky factor for the given matrix
    static void CholeskyColumnCounts( const cs *A,
                                      const IntArray &parent,
                                      const IntArray &postOrder,
                                      IntArray &columnCounts );

	private:
		CSparse_Interface();

		// Destructor
		virtual ~CSparse_Interface();

};

#endif
