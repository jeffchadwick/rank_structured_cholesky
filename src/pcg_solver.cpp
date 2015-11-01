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

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
  if ( argc < 2 ) {
    cerr << "Usage " << argv[ 0 ] << " <matrix file>" << endl;
    return 1;
  }

  SPARSE_MATRIX::SparseColumnMatrix    AsparseCol;
  VECTOR                               rhs;
  VECTOR                               solution;
  Vector3Array                         meshPositions;
  FactorPreconditioner                 factorPrecond;

  Timer                                solveTimer;

  bool                                 positionsFound = false;

  {
    char buf[ 1024 ];
    sprintf( buf, "%s_system.bcsm.gz", argv[ 1 ] );
    SPARSE_MATRIX::readFromBinaryGZ( AsparseCol, buf );
  }

  {
    char buf[ 1024 ];
    sprintf( buf, "%s_rhs.vector", argv[ 1 ] );

    rhs.read( buf );
    solution.resizeAndWipe( rhs.size() );
  }

  {
    char buf[ 1024 ];
    sprintf( buf, "%s_meshPositions.dat", argv[ 1 ] );

    positionsFound = readVector( meshPositions, buf );

    if ( !positionsFound ) {
      // If we did not find anything here, just populate this array
      // with a bunch of dummy positions
      meshPositions.clear();
      meshPositions.resize( AsparseCol._nrow, VEC3F( 0.0, 0.0, 0.0 ) );

      cerr << "Warning: no system degree of freedom positions found - "
              "using default separator ordering" << endl;
    }
  }

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

  return 0;
}
