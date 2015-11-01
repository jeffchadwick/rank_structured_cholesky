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
// FactorPreconditioner.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "FactorPreconditioner.h"

#include <rschol/ordering/NestedDissection.h>

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
FactorPreconditioner::FactorPreconditioner()
  : _useTranspose( false ),
    _solver( NULL ),
    _nCalls( 0 ),
    _symbolicFactorNeeded( true ),
    _numericFactorNeeded( true )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
FactorPreconditioner::~FactorPreconditioner()
{
  if (_solver != NULL) 
  {
    delete _solver;
  }
}

//////////////////////////////////////////////////////////////////////
// Initializes the solver and performs symbolic and numeric
// rank-structured factorization, if necessary
//////////////////////////////////////////////////////////////////////
void FactorPreconditioner::initialize(
                              const SPARSE_MATRIX::SparseColumnMatrix &matrix,
                              const Vector3Array &meshPositions,
                              int maxDenseBlock, int maxDiagonalBlock,
                              ExtendedFactorType extFactorType,
                              CompressionType compressionType,
                              DecompositionBlock blockType,
                              DecompositionType decompType,
                              bool reorderSeparators )
{
  if ( !_symbolicFactorNeeded && !_numericFactorNeeded ) {
    return;
  }

#if 0
  _currentSystem = &matrix;
#endif
  _currentSystem = matrix;

  if ( _symbolicFactorNeeded ) {
    // Clear the solver and rebuild it with a new symbolic factor
    delete _solver;
    _solver = NULL;

    printf( "FactorPreconditioner::initialize: initializing DoF positions\n" );

    NestedDissection::PointBuilder pointPositionQuery
                                = boost::bind( BuildMeshPoint, _1,
                                               boost::ref( meshPositions ) );

    printf( "FactorPreconditioner::initialize: building solver\n" );
    _solver = new SparseSolver( _currentSystem, pointPositionQuery,
                                maxDenseBlock, maxDiagonalBlock,
                                extFactorType, compressionType,
                                blockType, decompType,
                                reorderSeparators );

    // Run the factoriatization
    printf( "FactorPreconditioner::initialize: symbolic factorization\n" );
    _solver->factorSymbolic();

    // Make our workspace big enough to accomodate the extended system
    _rhsWorkspace.resizeAndWipe( _solver->extendedSystemSize() );

    _symbolicFactorNeeded = false;
  } else {
#if 0
    _solver->reorderSystem( *_currentSystem );
#endif
    _solver->reorderSystem( _currentSystem );
  }

  if ( _numericFactorNeeded ) {
    printf( "FactorPreconditioner::initialize: numeric factorization\n" );
#if 0
    _solver->factorNumeric( *_currentSystem );
#endif
    _solver->factorNumeric( _currentSystem );
    _numericFactorNeeded = false;
    printf( "FactorPreconditioner::initialize: done\n" );
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void FactorPreconditioner::solve( VECTOR &x, VECTOR &b )
{
  _rhsWorkspace.clear();

  const IntArray            &permutation = _solver->permutation();

  // Copy the source vector, applying the permutation when copying
  for ( int dof_idx = 0; dof_idx < b.size(); dof_idx++ ) {
    _rhsWorkspace( dof_idx ) = b[ permutation[ dof_idx ] ];
  }

  // Solve the system
#if 0
  _solver->solveSystem( *_currentSystem, _rhsWorkspace );
#endif
  _solver->solveSystem( _currentSystem, _rhsWorkspace );

  if ( x.size() != b.size() )
  {
    x.resizeAndWipe( b.size() );
  }

  // Copy back to the output vector, being careful to apply the permutation
  // in reverse
  for ( int dof_idx = 0; dof_idx < x.size(); dof_idx++ ) {
    x( permutation[ dof_idx ] ) = _rhsWorkspace( dof_idx );
  }
}
