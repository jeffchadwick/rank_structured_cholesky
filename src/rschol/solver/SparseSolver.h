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
// SparseSolver.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef SPARSE_SOLVER_H
#define SPARSE_SOLVER_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/SPARSE_MATRIX.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/ordering/NestedDissection.h>
#include <rschol/ordering/Ordering.h>

#include <rschol/solver/FactorManager.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

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
