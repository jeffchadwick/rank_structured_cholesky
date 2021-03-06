//////////////////////////////////////////////////////////////////////
// Copyright (c) 2015, Jeff Chadwick and David Bindel
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// SparseSolver.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "SparseSolver.h"

#include <rschol/solver/SolverError.h>

#include <rschol/util/STLUtil.h>

//////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////
SparseSolver::SparseSolver( SPARSE_MATRIX::SparseColumnMatrix &A,
                            NestedDissection::PointBuilder pointPositionQuery,
                            int maxDenseBlock, int maxDiagonalBlock,
                            ExtendedFactorType extFactorType,
                            CompressionType compressionType,
                            DecompositionBlock blockType,
                            DecompositionType decompType,
                            bool reorderSeparators )
  : _A( &A ),
    _pointPositionQuery( pointPositionQuery ),
    _maxDenseBlock( maxDenseBlock ),
    _maxDiagonalBlock( maxDiagonalBlock ),
    _manager( NULL ),
    _extFactorType( extFactorType ),
    _compressionType( compressionType ),
    _blockType( blockType ),
    _decompType( decompType ),
    _rankConstant( FactorManager::RANK_CONSTANT ),
    _diagonalRankConstant( FactorManager::DIAGONAL_RANK_CONSTANT ),
    _rankConstantAlpha( 1.25 /* FIXME */ ),
    _maxFactorAttempts( 10 /* FIXME */ ),
    _reorderSeparators( reorderSeparators )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
SparseSolver::~SparseSolver()
{
  printf( "Factoring used %f MB of slack memory\n", _manager->slackUsageMB() );

  delete _manager;
}

//////////////////////////////////////////////////////////////////////
// Symbolic factorization
//////////////////////////////////////////////////////////////////////
void SparseSolver::factorSymbolic( bool rebuildOrdering )
{
  delete _manager;

  if ( rebuildOrdering ) {
    _permutation.clear();
    _diagonalBlocks.clear();
    _supernodes.clear();

#if 0
    vector<vector<DenseBlock> >                diagonalBlocks;
#endif
    SPARSE_MATRIX::SparseColumnMatrix          reorderedSystem;
#if 0
    vector<Ordering::SupernodeSpecification>   supernodes;
#endif

    // Permute the matrix using nested dissection
#if 0
    NestedDissection::SupernodalNestedDissection( *_A, _pointPositionQuery,
                                                  _permutation, supernodes,
                                                  _maxDenseBlock,
                                                  _maxDiagonalBlock,
                                                  _compressionType,
                                                  diagonalBlocks );
#endif

    NestedDissection::SupernodalNestedDissection( *_A, _pointPositionQuery,
                                                  _permutation, _supernodes,
                                                  _maxDenseBlock,
                                                  _maxDiagonalBlock,
                                                  _compressionType,
                                                  _diagonalBlocks,
                                                  false, /* don't split nodes */
                                                  _reorderSeparators );

    // Reorder the system
    Ordering::reorderMatrix( *_A, reorderedSystem, _permutation,
                             true /* Only lower triangular part */ );

    // FIXME: debugging
    writeVector( _permutation, "SparseSolver_permutation.vector" );

    *_A = reorderedSystem;
  }

  // Build the factor manager and perform symbolic factorization
#if 0
  _manager = new FactorManager( *_A, supernodes, &diagonalBlocks,
                                FactorManager::ERROR_BOUND_MULTIPLIER,
                                _rankConstant, _diagonalRankConstant );
#endif
  _manager = new FactorManager( *_A, _supernodes, &_diagonalBlocks,
                                FactorManager::ERROR_BOUND_MULTIPLIER,
                                _rankConstant, _diagonalRankConstant );
  _manager->factorSymbolic( _maxDenseBlock, _compressionType,
                            // Whether or not to use interior blocks
                            true );
                            //false );
}

//////////////////////////////////////////////////////////////////////
// Just reorders a system according to our permutation
//////////////////////////////////////////////////////////////////////
void SparseSolver::reorderSystem( SPARSE_MATRIX::SparseColumnMatrix &A )
{
  SPARSE_MATRIX::SparseColumnMatrix          reorderedSystem;

  // Reorder the system
  Ordering::reorderMatrix( A, reorderedSystem, _permutation,
                           true /* Only lower triangular part */ );

  A = reorderedSystem;
}

//////////////////////////////////////////////////////////////////////
// Numeric factorization.  Right now this is just hard-coded to do
// fixed-rank factorization (using FactorManager::blockRank)
//////////////////////////////////////////////////////////////////////
void SparseSolver::factorNumeric( const SPARSE_MATRIX::SparseColumnMatrix &A )
{
  int                        factorAttempt = 0;

  while ( factorAttempt < _maxFactorAttempts ) {
    SolverErrorPtr           error;

    error = _manager->factorNumeric( A, false,
                                     /* Tolerances for adaptive factorization */
                                     0.1, 0.1,
                                     /* write size and whether to write blocks */
                                     -1, false,
                                     _extFactorType,
                                     _compressionType,
                                     _blockType,
                                     _decompType );

    if ( error != SolverErrorPtr() ) {
      if ( error->type() == SolverError::INDEFINITE_DIAGONAL ) {
        // Bump the diagonal rank constant and try again
        _diagonalRankConstant *= _rankConstantAlpha;

        printf( "Increasing rank constant to %f\n", _diagonalRankConstant );

        // Set the matrix pointer to the currently factored matrix
        // before re-running the symbolic factorization, since _A might
        // not actually exist anymore.  Also, *do not* rebuild the ordering,
        // since that part is still valid.
        _A = const_cast<SPARSE_MATRIX::SparseColumnMatrix *>(&A);

        factorSymbolic( false /* don't rebuild the ordering */ );
      } else {
        printf( "Unhandled solver error\n" );
        abort();
      }
    } else {
      break;
    }

    factorAttempt += 1;
  }
#if 0
  _manager->factorNumeric( _A, true, 0.01, 0.01 );
#endif
}
