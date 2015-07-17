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
