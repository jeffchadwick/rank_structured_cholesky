//////////////////////////////////////////////////////////////////////
// CHOLMOD_Supernodal.h: Interface for the CHOLMOD Supernodal class
//
//////////////////////////////////////////////////////////////////////

// (The MY_* business is so we don't interfere with CHOLMOD)

#ifndef MY_CHOLMOD_SUPERNODAL_H
#define MY_CHOLMOD_SUPERNODAL_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <rschol/util/timer.h>

#include <cholmod.h>

//////////////////////////////////////////////////////////////////////
// CHOLMOD_Supernodal class
//
// Wrapper for a supernodal CHOLMOD factorization.  This does
// things like figure out the structure of the compressed
// supernodal factorization, assign sparse right hand sides to
// appropriate supernodes, etc.
//////////////////////////////////////////////////////////////////////
class CHOLMOD_Supernodal {
	public:
    // Static instantiator, which is handy so that we can
    // make sure the factor satisfies certain conditions.
    static CHOLMOD_Supernodal *Build( cholmod_common &common,
                                      cholmod_factor *factor );

		// Destructor
		virtual ~CHOLMOD_Supernodal();

    // Return the compressed factorization
    const int *compressedP() const { return _p; }
    const int *compressedI() const { return _i; }
    int compressedN() const { return _factor->nsuper; }

    // Writes some debugging output for the compressed factorization
    void debugPrintCompressed();

    // Constructs the list of supernodes associated with
    // the non-zero entries of a sparse vector.  The provided
    // list should be at least as large as the most dense
    // column of the input matrix.
    // The output list may have duplicate supernodal indices,
    // though this should not be an issue if we do DFS with
    // markers.
    int getAssociatedSupernodes( const cholmod_sparse *B,
                                 int k, /* column to consider */
                                 int *nodes ) const;

    // Solves the given sparse right hand side system using this
    // factorization.  This only handles unpermuted lower triangular
    // solves.
    //    E is workspace - should have size at least _factor->maxesize
    //    marker_list is just a copy of _p for marking purposes
    //    k is the column of B to use as the RHS
    //    nodes is used for storage when finding supernodes associated
    //      with a particular RHS
    //    xi is the set of supernode indices that will be
    //      filled in to the solution
    //    x is the dense solution vector.  Does not need
    //      to be initialized.
    void spsolve( const cholmod_sparse *B,
                  cholmod_dense *E,
                  int *marker_list,
                  int k, int *nodes, int *xi, double *x );

    // Constructs an appropriately sized solverworkspace
    // for the given number of right hand sides
    cholmod_dense *buildSolverWorkspace( int nrhs )
    {
#if 0
      return cholmod_allocate_dense( nrhs, _factor->maxesize, nrhs,
                                     _factor->xtype, &_common );
#endif
      return cholmod_zeros( nrhs, _factor->maxesize, _factor->xtype, &_common );
    }

    Timer                      _reachTimer;
    Timer                      _solveTimer;

	protected:

  private:
    // Constructor:
    // We need the factorization itself, and whatever cholmod
    // workspace it is associated with.
		CHOLMOD_Supernodal( cholmod_common &common, cholmod_factor *factor );

    // Builds a map of matrix indices to supernodes
    void buildIndexToSupernode();

    // Compresses a indices in to supernodes and produces a
    // symbolic "factorization" which encodes the graph connecting
    // supernodes.
    void buildCompressedSymbolicFactor();

	private:
    cholmod_common            &_common;
    cholmod_factor            *_factor;

    // Maps matrix indices to supernodes
    IntArray                   _indexToSupernode;

    // Storage for compressed "factorization"
    int                       *_p;  // Column pointers
    int                       *_i;  // Row indices

};

#endif
