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
// CSparse_Interface.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef CSPARSE_INTERFACE_H
#define CSPARSE_INTERFACE_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/SPARSE_MATRIX.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

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
