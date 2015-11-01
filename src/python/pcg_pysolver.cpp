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
// testPCG.cpp: Solves the given input system using PCG with our
//              rank-structured solver
//
//////////////////////////////////////////////////////////////////////

#include <rschol/linearalgebra/PCG_SOLVER.h>
#include <rschol/linearalgebra/SPARSE_MATRIX.h>

#include <rschol/solver/FactorPreconditioner.h>

#include <rschol/util/STLUtil.h>
#include <rschol/util/timer.h>


void pcg_pysolve(SPARSE_MATRIX::SparseColumnMatrix& AsparseCol,
                 VECTOR& rhs,
                 VECTOR& solution,
                 Vector3Array& meshPositions,
                 bool positionsFound = false)
{
  FactorPreconditioner                 factorPrecond;

  Timer                                solveTimer;
  solveTimer.tick();

  factorPrecond.initialize( AsparseCol, meshPositions,
                            200, 400,
                            //200, 400000,
                            //200000, 4000000,
                            CHOLESKY,
                            IN_PLACE,
                            DECOMPOSE_TRANSPOSE,
                            DECOMPOSE_FACTOR,
                            // Only attempt to reorder separators if we
                            // found DoF positions
                            positionsFound );

  // Relative residual error of 1e-6 and max iterations of 10000
  PCG_SOLVER                 solver( AsparseCol, 1e-6, 10000 );

  solver.solvePCG( solution, rhs, &factorPrecond );

  solveTimer.tock();

  printf( "\nPCG solution took %f seconds\n", solveTimer.getTotalSecs() );
}


void pcg_pysolve(int nrow, const int* jc, const int* ir, const double* pr,
                 const double* rhs, double* solution)
{
    SPARSE_MATRIX::SparseColumnMatrix A;
    Vector3Array meshPositions;
    A._nrow = nrow;
    A._ncol = nrow;
    A._nzmax = jc[nrow];
    A.allocate();
    std::copy(jc, jc+nrow+1, A._p);
    std::copy(ir, ir+A._nzmax, A._i);
    std::copy(pr, pr+A._nzmax, A._x);
    meshPositions.clear();
    meshPositions.resize(nrow, VEC3F( 0.0, 0.0, 0.0 ) );
    VECTOR v(nrow, (double*) rhs);
    VECTOR s(nrow);
    pcg_pysolve(A, v, s, meshPositions);
    std::copy(s.data(), s.data()+nrow, solution);
}
