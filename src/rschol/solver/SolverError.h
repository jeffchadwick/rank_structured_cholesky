//////////////////////////////////////////////////////////////////////
// SolverError.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef SOLVER_ERROR_H
#define SOLVER_ERROR_H

#include <boost/shared_ptr.hpp>

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

//////////////////////////////////////////////////////////////////////
// SolverError class
//
// Comments
//////////////////////////////////////////////////////////////////////
class SolverError {
  public:
    SolverError()
    {
    }

    // Destructor
    virtual ~SolverError()
    {
    }

    enum ErrorType {
      INDEFINITE_DIAGONAL = 0,
      NUM_TYPES
    };

    virtual ErrorType type() const = 0;

  protected:

  private:
};

//////////////////////////////////////////////////////////////////////
// Error produced when an indefinite diagonal block is encountered
//////////////////////////////////////////////////////////////////////
class IndefiniteDiagonalError : public SolverError {
  public:
    IndefiniteDiagonalError()
    {
    }

    virtual ~IndefiniteDiagonalError()
    {
    }

    ErrorType type() const
    {
      return INDEFINITE_DIAGONAL;
    }
};

typedef boost::shared_ptr<SolverError>   SolverErrorPtr;

#endif
