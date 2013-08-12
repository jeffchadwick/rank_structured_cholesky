//////////////////////////////////////////////////////////////////////
// CSparse_Interface.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "CSparse_Interface.h"

#include <util/IO.h>

//////////////////////////////////////////////////////////////////////
// Builds a CSparse formatted sparse matrix out of a column sparse
// matrix in our format.  Optionally permutes the matrix.  If a
// permutation is provided then permute the matrix.  If permutation
// is smaller than the number of columns in the matrix then the
// remainder of the matrix will be filled in with diagonals.
//////////////////////////////////////////////////////////////////////
void CSparse_Interface::CopySparseSquareSymbolicMatrix(
                               const SPARSE_MATRIX::SparseColumnMatrix &A,
                               cs * &Acopy,
                               const IntArray *permutation,
                               const IntArray *inversePermutation,
                               MatrixStorage storageType,
                               bool storeDiagonal )
{
  int                        nzmax = 0;
  int                        entryIndex = 0;
  int                        nCols;
  int                        permuted_col_idx;
  int                        permuted_row_idx;

  IntArray                   rowCounts( A._ncol, 1 );
  
  bool                       usePermutation;

  usePermutation = permutation && inversePermutation;
  
  nCols = permutation ? permutation->size() : A._ncol;

  TRACE_ASSERT( nCols <= A._ncol );
  TRACE_ASSERT( A._nrow == A._ncol );

  if ( Acopy )
  {
    cs_spfree( Acopy );
    Acopy = NULL;
  }

  // Count non-zeros from the original matrix
  for ( int col_idx = 0; col_idx < nCols; col_idx++ )
  {
    permuted_col_idx = usePermutation ? permutation->at( col_idx ) : col_idx;

    rowCounts[ col_idx ] = 0;

    for ( int row_ptr = A._p[ permuted_col_idx ];
          row_ptr < A._p[ permuted_col_idx + 1 ]; row_ptr++ )
    {
      permuted_row_idx
        = ( usePermutation && A._i[ row_ptr ] < nCols )
                         ? inversePermutation->at( A._i[ row_ptr ] )
                         : A._i[ row_ptr ];

      // Skip the diagonal if that's what we want
      if ( !storeDiagonal && permuted_row_idx == col_idx )
        continue;

      if ( storageType == FULL 
        || ( storageType == LOWER && permuted_row_idx >= col_idx )
        || ( storageType == UPPER && permuted_row_idx <= col_idx ) )
      {
        nzmax += 1;
        rowCounts[ col_idx ] += 1;
      }
    }
  }

  cout << SDUMP( rowCounts[ 0 ] ) << endl;

  // Add diagonal entries for any column not added
  nzmax += A._ncol - nCols;

  // Build the new matrix
  Acopy = cs_spalloc( A._nrow, A._ncol, nzmax,
                      false, /* no numerical values */
                      false /* not in triplet format */ );

  Acopy->p[ 0 ] = 0;

  // Fill in the matrix
  for ( int col_idx = 0; col_idx < nCols; col_idx++ )
  {
    permuted_col_idx = usePermutation ? permutation->at( col_idx ) : col_idx;

    for ( int row_ptr = A._p[ permuted_col_idx ];
          row_ptr < A._p[ permuted_col_idx + 1 ]; row_ptr++ )
    {
      permuted_row_idx
        = ( usePermutation && A._i[ row_ptr ] < nCols )
                         ? inversePermutation->at( A._i[ row_ptr ] )
                         : A._i[ row_ptr ];

      // Skip the diagonal if that's what we want
      if ( !storeDiagonal && permuted_row_idx == col_idx )
        continue;

      if ( storageType == FULL 
        || ( storageType == LOWER && permuted_row_idx >= col_idx )
        || ( storageType == UPPER && permuted_row_idx <= col_idx ) )
      {
        Acopy->i[ entryIndex ] = permuted_row_idx;
        entryIndex++;
      }
    }

    Acopy->p[ col_idx + 1 ] = entryIndex;
  }

  // Fill in diagonal entries for any column not added
  for ( int col_idx = nCols; col_idx < A._ncol; col_idx++ )
  {
    Acopy->i[ entryIndex ] = col_idx;
    entryIndex++;
    Acopy->p[ col_idx + 1 ] = entryIndex;
  }

  printf( "Original non-zeros = %d, copies non-zeros = %d\n",
          (int)A._nzmax, Acopy->nzmax );
}

//////////////////////////////////////////////////////////////////////
// Copies a sparse symmetric symbolic matrix
//////////////////////////////////////////////////////////////////////
void CSparse_Interface::CopySparseSymmetricSymbolicMatrix(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              cs * &Acopy,
                              const IntArray *permutation,
                              const IntArray *inversePermutation )
{
  cs                        *A_lowerTri = NULL;
  cs                        *A_upperTri = NULL;

  // Get the lower triangular representation of the matrix first
  CopySparseSquareSymbolicMatrix( A, A_lowerTri,
                                  permutation, inversePermutation,
                                  // FIXME: Probably want upper
                                  LOWER /* lower triangular only */ );

  // Symbolic add to get symmetrized version
  TransposeSparseSymbolicMatrix( A_lowerTri, A_upperTri );
  AddSparseSymbolicMatrices( A_lowerTri, A_upperTri, Acopy );

  cs_spfree( A_lowerTri );
  cs_spfree( A_upperTri );
}

//////////////////////////////////////////////////////////////////////
// Transpose routine
//////////////////////////////////////////////////////////////////////
void CSparse_Interface::TransposeSparseSymbolicMatrix(
                                          const cs *A, cs * &Atranspose )
{
  if ( Atranspose )
  {
    cs_spfree( Atranspose );
    Atranspose = NULL;
  }

  Atranspose = cs_transpose( A, false /* symbolic only */ );
}

//////////////////////////////////////////////////////////////////////
// Symbolic matrix add
//////////////////////////////////////////////////////////////////////
void CSparse_Interface::AddSparseSymbolicMatrices(
                                          const cs *A, const cs *B, cs * &C )
{
  if ( C )
  {
    cs_spfree( C );
    C = NULL;
  }

  C = cs_add( A, B, 1.0, 1.0 /* constants don't matter - symbolic add */ );
}

//////////////////////////////////////////////////////////////////////
// Constructs the elimination tree for the given matrix
//////////////////////////////////////////////////////////////////////
void CSparse_Interface::ConstructEliminationTree( const cs *A,
                                                  IntArray &parent )
{
  int                       *parentArray = NULL;

  TRACE_ASSERT( A->n == A->m, "Matrix is not square" );

  parent.clear();
  parent.resize( A->n );

  parentArray = cs_etree( A, 0 /* Do 'not' compute A' * A */ );

  for ( int i = 0; i < A->n; i++ )
  {
    parent[ i ] = parentArray[ i ];
  }

  cs_free( parentArray );
}

//////////////////////////////////////////////////////////////////////
// Builds an elimination tree post order
//////////////////////////////////////////////////////////////////////
void CSparse_Interface::ConstructPostOrder( const IntArray &parent,
                                            IntArray &postOrder )
{
  int                       *postArray = NULL;

  postOrder.clear();
  postOrder.resize( parent.size() );

  postArray = cs_post( parent.data(), parent.size() );

  for ( int i = 0; i < parent.size(); i++ )
  {
    postOrder[ i ] = postArray[ i ];
  }

  cs_free( postArray );
}

//////////////////////////////////////////////////////////////////////
// Given an elimination tree and its post ordering, determine column
// counts in the Cholesky factor for the given matrix
//////////////////////////////////////////////////////////////////////
void CSparse_Interface::CholeskyColumnCounts( const cs *A,
                                              const IntArray &parent,
                                              const IntArray &postOrder,
                                              IntArray &columnCounts )
{
  int                       *columnCountArray = NULL;

  columnCounts.clear();
  columnCounts.resize( parent.size() );

  columnCountArray = cs_counts( A, parent.data(), postOrder.data(),
                                0 /* Do 'not' compute A' * A */ );

  for ( int i = 0; i < parent.size(); i++ )
  {
    columnCounts[ i ] = columnCountArray[ i ];
  }

  cs_free( columnCountArray );
}

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
CSparse_Interface::CSparse_Interface()
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
CSparse_Interface::~CSparse_Interface()
{
}

