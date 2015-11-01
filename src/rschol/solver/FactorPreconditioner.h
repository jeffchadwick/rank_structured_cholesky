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
// FactorPreconditioner.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef FACTOR_PRECONDITIONER_H
#define FACTOR_PRECONDITIONER_H

#include <rschol/linearalgebra/PRECONDITIONER.h>
#include <rschol/linearalgebra/SPARSE_MATRIX.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/mesh/Mesh.h>

#include <rschol/solver/SparseSolver.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <iostream>

//////////////////////////////////////////////////////////////////////
// FactorPreconditioner class
//
// Preconditioner based on a rank-structured approximate factor
//////////////////////////////////////////////////////////////////////
class FactorPreconditioner : public PRECONDITIONER {
	public:
		FactorPreconditioner();

		// Destructor
		virtual ~FactorPreconditioner();

    void clearSystem()
    {
      delete _solver;
      _solver = NULL;

      _symbolicFactorNeeded = true;
      _numericFactorNeeded = true;
    }

    void clearFactor()
    {
      if ( !_solver ) {
        _symbolicFactorNeeded = true;
      }

      _numericFactorNeeded = true;
    }

    void setUseTranspose( bool useTranspose )
    {
      _useTranspose = useTranspose;
    }

    // Initializes the solver and performs symbolic and numeric
    // rank-structured factorization, if necessary
    void initialize( const SPARSE_MATRIX::SparseColumnMatrix &matrix,
                     const Vector3Array &meshPositions,
                     int maxDenseBlock, int maxDiagonalBlock,
                     ExtendedFactorType extFactorType,
                     CompressionType compressionType,
                     DecompositionBlock blockType,
                     DecompositionType decompType,
                     // Whether or not to reorder supernode diagonal blocks
                     // according to the positions from meshPositions
                     // Generally this should be true, assuming that
                     // meshPositions provides meaningful information
                     bool reorderSeparators );

    // This solver must be initialized outside of the PCG solver
    virtual void init() {}
    virtual void init( SPARSE_MATRIX &matrix )
    {
      std::cerr << "FactorPreconditioner::init(SPARSE_MATRIX &)"
                   " not implemented" << std::endl;
    }

    virtual void solve( VECTOR &x, VECTOR &b );

	protected:

  private:
    static Mesh::MeshPoint *BuildMeshPoint( int index,
                                            const Vector3Array &positions )
    {
      return new Mesh::MeshPoint( index, positions[ index ] );
    }

	private:
    bool                               _useTranspose;

    // The actual sparse matrix used by the preconditioner
    SPARSE_MATRIX::SparseColumnMatrix  _currentSystem;

    SparseSolver                      *_solver;

    VECTOR                             _rhsWorkspace;

    int                                _nCalls;

    bool                               _symbolicFactorNeeded;
    bool                               _numericFactorNeeded;

};

#endif
