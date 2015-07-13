// PCG_SOLVER.h: interface for the PCG_SOLVER class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PCG_SOLVER_H
#define PCG_SOLVER_H

#include <rschol/linearalgebra/SPARSE_MATRIX.h>
#include <rschol/linearalgebra/PRECONDITIONER.h>

//////////////////////////////////////////////////////////////////////
// A sparse matrix class that can solve using Pardiso
//////////////////////////////////////////////////////////////////////
class PCG_SOLVER {

public:
  PCG_SOLVER(SPARSE_MATRIX::SparseColumnMatrix &A,
             Real eps = 1e-4, int iterations = 1000);
  ~PCG_SOLVER();

  // solve an SPD system using PCG
  void solveCG(VECTOR& x, VECTOR& b);
  //void solveICCG(VECTOR& x, VECTOR& b);
  void solvePCG(VECTOR& x, VECTOR& b, PRECONDITIONER* M);

  int& maxIterations() { return _iterations; };
  Real& eps() { return _eps; };
  Real meanIterations() { return (Real)_totalIterations / _totalSolves; };
  Real meanResidual()   { return _totalResidual / _totalSolves; };
  VECTOR& residual()       { return _residual; };
  VECTOR& direction()      { return _direction; };
  VECTOR& q()              { return _q; };

  map<string, double>& timingBreakdown() { return _timingBreakdown; };

  double& totalTime() { return _totalTime; };

protected:
  // do allocations for symmetric matrix
  void initCG();
  void initPCG();

protected:
  // The actual matrix to solve with
  SPARSE_MATRIX::SparseColumnMatrix &_A;

  // has a factorization been done before?
  bool _firstSolve;

  // store the size of the matrix that the solver expects
  int _matrixSize;

  // CG variables
  VECTOR _residual;
  VECTOR _direction;
  VECTOR _q;
  int _iterations;
  Real _eps;

  // PCG variables
  VECTOR _s;
 
  // keep track of some stats
  int _totalSolves;
  int _totalIterations;
  Real _totalResidual;

#if 0
  // do a matrix-vector multiply that OpenMP can handle
  void parallelMultiply();
#endif
  //vector<pair<int,int> > _pairs;
  vector<int> _rowIndices;
  vector<int> _columnIndices;
  vector<Real> _entries;
  //vector<VECTOR> _qs;
  VECTOR* _qs;
  int _totalCores;

#if 0
  // precache as much as possible before going into the main PCG loop
  void initParallelMultiply();
#endif

  map<string, double> _timingBreakdown;
  double _totalTime;
};

#endif
