//////////////////////////////////////////////////////////////////////
// testPCG.cpp: Solves the given input system using PCG with our
//              rank-structured solver
//
//////////////////////////////////////////////////////////////////////

#include <linearalgebra/PCG_SOLVER.h>
#include <linearalgebra/SPARSE_MATRIX.h>

#include <solver/FactorPreconditioner.h>

#include <util/STLUtil.h>
#include <util/timer.h>

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
