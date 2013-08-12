//////////////////////////////////////////////////////////////////////
// FactorPreconditioner.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef FACTOR_PRECONDITIONER_H
#define FACTOR_PRECONDITIONER_H

#include <linearalgebra/PRECONDITIONER.h>
#include <linearalgebra/SPARSE_MATRIX.h>
#include <linearalgebra/VECTOR.h>

#include <mesh/Mesh.h>

#include <solver/SparseSolver.h>

#include <SETTINGS.h>
#include <TYPES.h>

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
