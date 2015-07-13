//////////////////////////////////////////////////////////////////////
// SparseSolver.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef SPARSE_SOLVER_H
#define SPARSE_SOLVER_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/SPARSE_MATRIX.h>
#include <linearalgebra/VECTOR.h>

#include <ordering/NestedDissection.h>
#include <ordering/Ordering.h>

#include <solver/FactorManager.h>

#include <SETTINGS.h>
#include <TYPES.h>

//////////////////////////////////////////////////////////////////////
// SparseSolver class
//
// Comments
//////////////////////////////////////////////////////////////////////
class SparseSolver {
	public:
    // Provide the sparse matrix, and a function pointer which returns
    // a spatial position given a row/column index.
    //
    // maxDenseBlock and maxDiagonalBlock refer to the maximum sizes of
    // dense off-diagonal and diagonal block appearing in the factor
    // (blocks larger than these will be compressed)
    //
    // maxBlockSize refers to 
    //
    // NOTE: The matrix will be overwritten with its permuted version
		SparseSolver( SPARSE_MATRIX::SparseColumnMatrix &A,
                  NestedDissection::PointBuilder pointPositionQuery,
                  int maxDenseBlock, int maxDiagonalBlock,
                  ExtendedFactorType extFactorType,
                  CompressionType compressionType,
                  DecompositionBlock blockType,
                  DecompositionType decompType,
                  // Whether or not to reorder supernode diagonal blocks
                  // according to the positions from pointPositionQuery
                  // Generally this should be true, assuming that
                  // pointPositionQuery provides meaningful information
                  bool reorderSeparators = true );

		// Destructor
		virtual ~SparseSolver();

    // Symbolic factorization
    void factorSymbolic( bool rebuildOrdering = true );
    
    // Just reorders a system according to our permutation
    void reorderSystem( SPARSE_MATRIX::SparseColumnMatrix &A );

    // Numeric factorization.  Right now this is just hard-coded to do
    // fixed-rank factorization (using FactorManager::blockRank)
    void factorNumeric( const SPARSE_MATRIX::SparseColumnMatrix &A );
    void factorNumeric() {
      factorNumeric( *_A );
    }

    // Solves the given system, once the factor has been computed.
    // Overwrites the provided vector with the solution.
    void solveSystem( const SPARSE_MATRIX::SparseColumnMatrix &A,
                      VECTOR &rhs )
    {
      _manager->solveSystem( A, rhs );
    }

    const IntArray &permutation() const
    {
      return _permutation;
    }

    int extendedSystemSize() const
    {
      return _manager->factorSystemSize();
    }

    FactorManager &manager()
    {
      return *_manager;
    }

    void writeTimings()
    {
      _manager->writeTimings();
    }

    void printUsage()
    {
      _manager->printExtendedUsage();
    }

    void clear()
    {
      _manager->clear();
    }

	private:

	private:
    SPARSE_MATRIX::SparseColumnMatrix   *_A;
    NestedDissection::PointBuilder       _pointPositionQuery;

    int                                  _maxDenseBlock;
    int                                  _maxDiagonalBlock;

    // System permutation and supernode set (in the permuted matrix)
    IntArray                                   _permutation;
    vector<vector<DenseBlock> >                _diagonalBlocks;
    vector<Ordering::SupernodeSpecification>   _supernodes;

    // This performs the actual factorization
    FactorManager                       *_manager;

    ExtendedFactorType                   _extFactorType;
    CompressionType                      _compressionType;
    DecompositionBlock                   _blockType;
    DecompositionType                    _decompType;

    // Rank constants for decomposition
    Real                                 _rankConstant;
    Real                                 _diagonalRankConstant;

    Real                                 _rankConstantAlpha;

    int                                  _maxFactorAttempts;

    bool                                 _useLDL;

    // Whether or not to reorder supernode diagonal blocks
    // according to the positions from pointPositionQuery
    // Generally this should be true, assuming that
    // pointPositionQuery provides meaningful information
    bool                                 _reorderSeparators;

};

#endif
