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



// PCG_SOLVER.h: interface for the PCG_SOLVER class.
//
//////////////////////////////////////////////////////////////////////

#include "PCG_SOLVER.h"

#include <rschol/util/timer.h>

#ifdef USING_MKL
#include <mkl_pardiso.h>
#endif
#ifdef USE_OMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor for the solver
//////////////////////////////////////////////////////////////////////
PCG_SOLVER::PCG_SOLVER(SPARSE_MATRIX::SparseColumnMatrix &A,
                       Real eps, int iterations)
  : _A(A),
    _iterations(iterations),
    _eps(eps),
    _firstSolve(true),
    _totalIterations(0),
    _totalSolves(0),
    _totalResidual(0.0),
    _totalTime(0.0),
    _qs(NULL),
    _totalCores(0)
{
}

PCG_SOLVER::~PCG_SOLVER()
{
  if (_qs) delete[] _qs;
}

//////////////////////////////////////////////////////////////////////
// solve an SPD system using CG
//////////////////////////////////////////////////////////////////////
void PCG_SOLVER::solveCG(VECTOR& x, VECTOR& b)
{
#if 0
  Timer total;
#endif
  if (_firstSolve)
  {
    initCG();
    _firstSolve = false;
  }

	// r = b - Ax
  //_residual = b - (*this) * x;
  Real temp = b.maxValue();
#if 0
  _residual = (*this) * x;
#endif
  _residual = _A * x;
  _residual = b - _residual;

  temp = _residual.maxValue();

	// d = r
  _direction.copyInplace(_residual);
  
	// deltaNew = transpose(r) * r
  Real deltaNew = _residual ^ _residual;
  
	// delta0 = deltaNew
  Real delta0 = deltaNew;

#if 0
#ifdef USE_OMP
  // precache as much as possible before going into the main PCG loop
  initParallelMultiply();
#endif
#endif

	// While deltaNew > (eps^2) * delta0
	Real maxR = 2.0f * _eps;
  int i = 0;
#ifndef _WIN32  
  cout << " Residual: ";
#endif
	while ((i < _iterations) && (maxR > _eps))
	{
#if 0
		// q = Ad
#ifdef USE_OMP
    parallelMultiply();
#else
    Timer multiply;
    _q = (*this) * _direction;
    _timingBreakdown["Matrix multiply"] += multiply.timing();
#endif
#endif
    _q = _A * _direction;

		// alpha = deltaNew / (transpose(d) * q)
#if 0
    Timer alphaTimer;
#endif
    Real alpha = _direction ^ _q;
		if (fabs(alpha) > 0.0) {
      alpha = deltaNew / alpha;
    }
#if 0
    _timingBreakdown["alpha compute"] += alphaTimer.timing();
#endif

		// x = x + alpha * d
#if 0
    Timer xTimer;
#endif
    x.axpy(alpha, _direction);
#if 0
    _timingBreakdown["x update"] += xTimer.timing();
#endif

		// r = r - alpha * q
#if 0
    Timer residualTimer;
#endif
    _residual.axpy(-alpha, _q);
#if 0
    _timingBreakdown["residual update"] += residualTimer.timing();
#endif

		// deltaOld = deltaNew
		Real deltaOld = deltaNew;

		// deltaNew = transpose(r) * r
#if 0
    Timer deltaTimer;
#endif
		deltaNew = _residual ^ _residual;
#if 0
    _timingBreakdown["delta update"] += deltaTimer.timing();
#endif

		// beta = deltaNew / deltaOld
		Real beta = deltaNew / deltaOld;

		// d = r + beta * d
#if 0
    Timer directionTimer;
#endif
    _direction *= beta;
    _direction += _residual;
#if 0
    _timingBreakdown["direction update"] += directionTimer.timing();
#endif

    // maxR = max(r);
#if 0
    Timer maxRTimer;
#endif
    //maxR = _residual.maxValue();
    //maxR = max( _residual.maxValue(), abs( _residual.minValue() ) );
		maxR = _residual.norm2() / b.norm2();
#if 0
    _timingBreakdown["maxR update"] += maxRTimer.timing();
#endif

		// i = i + 1
		i++;

#ifndef _WIN32  
    if (i % 10 == 0)
      cout << maxR << " ";
#endif
  }
#ifndef _WIN32  
  cout << endl;
#endif
  _totalIterations += i;
  _totalSolves++;
  _totalResidual += maxR;

  cout << " PCG iterations: " << i << endl;
#if 0
  _totalTime += total.timing();
  cout << " PCG timing: " << _totalTime << endl;
#endif
}

//////////////////////////////////////////////////////////////////////
// solve an SPD system using PCG
//////////////////////////////////////////////////////////////////////
void PCG_SOLVER::solvePCG(VECTOR& x, VECTOR& b, PRECONDITIONER* M)
{
#if 0
  Timer total;
#endif

  if (_firstSolve)
  {
    initCG();
    initPCG();
    _firstSolve = false;
  }

#if 0
  Timer Msetup;
#endif
  // initialize the preconditioner
  M->init();
#if 0
  _timingBreakdown["Precond. Setup"] += Msetup.timing();
#endif

	// r = b - Ax
  //_residual = b - (*this) * x;
#if 0
  _residual = (*this) * x;
#endif
  _residual = _A * x;
  _residual = b - _residual;

  // d = (M^-1) * r
  M->solve(_direction, _residual);

	// deltaNew = transpose(r) * d
  Real deltaNew = _residual ^ _direction;
  
	// delta0 = deltaNew
  Real delta0 = deltaNew;

#if 0
#ifdef USE_OMP
  // precache as much as possible before going into the main PCG loop
  initParallelMultiply();
#endif
#endif

	// While deltaNew > (eps^2) * delta0
	Real maxR = 2.0f * _eps;
  int i = 0;
  //cout << " Residual: ";
	while ((i < _iterations) && (maxR > _eps))
	{
#if 0
#ifdef USE_OMP
    parallelMultiply();
#else
    Timer multiply;
    _q = (*this) * _direction;
    _timingBreakdown["Matrix multiply"] += multiply.timing();
#endif
#endif
    _q = _A * _direction;

		// alpha = deltaNew / (transpose(d) * q)
    Real alpha = _direction ^ _q;
		if (fabs(alpha) > 0.0)
      alpha = deltaNew / alpha;
    else
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " PCG BREAKDOWN" << endl;
    }

		// x = x + alpha * d
    x.axpy(alpha, _direction);

		// r = r - alpha * q
    _residual.axpy(-alpha, _q);

		// s = (M^-1) * r
#if 0
    Timer Msolve;
#endif
    M->solve(_s, _residual);
#if 0
    _timingBreakdown["Precond. solve"] += Msolve.timing();
#endif

		// deltaOld = deltaNew
		Real deltaOld = deltaNew;

		// deltaNew = transpose(r) * s
		deltaNew = _residual ^ _s;

		// beta = deltaNew / deltaOld
		Real beta = deltaNew / deltaOld;

		// d = r + beta * d
    _direction *= beta;
    _direction += _s;

    // maxR = max(r);
#if 0
    maxR = _residual.maxValue();
#endif
		// FIXME
		maxR = _residual.norm2() / b.norm2();

		// i = i + 1
		i++;

    if (i % 10 == 0) {
      //cout << "Iteration " << i << ": " << maxR << " " << flush;
      printf( "Iteration %5d: %e\n", i, maxR );
    }
  }
  cout << endl;

  _totalIterations += i;
  _totalSolves++;
  _totalResidual += maxR;
  cout << " total iterations: " << i << endl;
  //_residual = b - (*this) * x;
#if 0
  _residual = (*this) * x;
#endif
  _residual = _A * x;
  _residual = b - _residual;
#if 0
  cout << " residual: " << _residual.normInf() << endl;
#endif
	// FIXME
  cout << " residual: " << (_residual.norm2() / b.norm2()) << endl;

#if 0
  _totalTime += total.timing();
#endif
}

//////////////////////////////////////////////////////////////////////
// initialize CG solver for first solve
//////////////////////////////////////////////////////////////////////
void PCG_SOLVER::initCG()
{
  int size = _A._nrow;
  _residual.resizeAndWipe(size);
  _direction.resizeAndWipe(size);
  _q.resizeAndWipe(size);
}

//////////////////////////////////////////////////////////////////////
// initialize PCG solver for first solve
//////////////////////////////////////////////////////////////////////
void PCG_SOLVER::initPCG()
{
  int size = _A._nrow;
  _s.resizeAndWipe(size);
}

#if 0
//////////////////////////////////////////////////////////////////////
// precache as much as possible before going into the main PCG loop
//////////////////////////////////////////////////////////////////////
void PCG_SOLVER::initParallelMultiply()
{
  // if this is the first call, cache the pair and Real locations
  if (_entries.size() == 0)
  {
    Timer cacheTimer;
    map<pair<int,int>, Real>::iterator iter;
    for (iter = _matrix.begin(); iter != _matrix.end(); iter++)
    {
      int row = iter->first.first;
      int col = iter->first.second;
      //_pairs.push_back(iter->first);
      _rowIndices.push_back(row);
      _columnIndices.push_back(col);
      _entries.push_back(iter->second);
    }
#ifdef USE_OMP
    //_qs.resize(omp_get_max_threads());
    _totalCores = omp_get_max_threads();
#else
    //_qs.resize(1);
    _totalCores = 1;
#endif
    _qs = new VECTOR[_totalCores];
    //for (int x = 0; x < _qs.size(); x++)
    for (int x = 0; x < _totalCores; x++)
      _qs[x].resizeAndWipe(_q.size());
    _timingBreakdown["Pair caching"] += cacheTimer.timing();
  }
  else
  {
    Timer cacheTimer;
    int i = 0;
    map<pair<int,int>, Real>::iterator iter;
    for (iter = _matrix.begin(); iter != _matrix.end(); iter++, i++)
      _entries[i] = iter->second;
    _timingBreakdown["Entry caching"] += cacheTimer.timing();
  }
}
#endif

#if 0
//////////////////////////////////////////////////////////////////////
// do a parallel multiply that OpenMP can handle
//////////////////////////////////////////////////////////////////////
void PCG_SOLVER::parallelMultiply()
{
  Timer multiply;
#ifdef USE_OMP
#pragma omp parallel
#endif
  { 
#ifdef USE_OMP
    const int id  = omp_get_thread_num();
#else
    const int id  = 0;
#endif
    VECTOR& q = _qs[id];
    q.clear();
#pragma omp for  schedule(static)
    for (int x = 0; x < _entries.size(); x++)
    {
      //const pair<int,int> index = _pairs[x];
      //const Real value = *(_entries[x]);
      //int i = index.first;
      //int j = index.second;
      //_qs[id](i) += _direction[j] * value;
      //_qs[id](index.first) += _direction[index.second] * _entries[x];
      //_qs[id](_rowIndices[x]) += _direction[_columnIndices[x]] * _entries[x];
      q[_rowIndices[x]] += _direction[_columnIndices[x]] * _entries[x];
    }
  }
  _timingBreakdown["Matrix multiply"] += multiply.timing();

  // combine qs
  Timer combineTimer;
  _q.clear();
  for (int x = 0; x < _totalCores; x++)
    _q += _qs[x];
  _timingBreakdown["Final combine"] += combineTimer.timing();
}
#endif

