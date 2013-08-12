// PRECONDITIONER.h: interface for the PRECONDITIONER class.
//
//////////////////////////////////////////////////////////////////////

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

#include "SPARSE_MATRIX.h"

//////////////////////////////////////////////////////////////////////
// Preconditioner for SPARSE_PCG_MATRIX
//////////////////////////////////////////////////////////////////////
class PRECONDITIONER {

public:
  PRECONDITIONER()  {};
  virtual ~PRECONDITIONER() {};

  virtual void init() = 0;
	virtual void init(SPARSE_MATRIX &matrix) = 0;
  virtual void solve(VECTOR& x, VECTOR& b) = 0;
};

#endif
