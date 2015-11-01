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
// CHOLMOD_Environment.h: Interface for the CHOLMOD_Environment class
//
//////////////////////////////////////////////////////////////////////

#ifndef CHOLMOD_ENVIRONMENT_H
#define CHOLMOD_ENVIRONMENT_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/SPARSE_MATRIX.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <rschol/library_interface/CHOLMOD_Supernodal.h>

#include <rschol/util/IO.h>
#include <rschol/util/timer.h>

#include <cholmod.h>

//////////////////////////////////////////////////////////////////////
// CHOLMOD_Environment class
//
// Provides functions for interacting with the CHOLMOD library.
// We will need one of these objects wherever we want to factorize
// a matrix.
//////////////////////////////////////////////////////////////////////
class CHOLMOD_Environment {
	public:
		CHOLMOD_Environment( int ordering = CHOLMOD_METIS );

		// Destructor
		virtual ~CHOLMOD_Environment();

    // Sets a new matrix to factorize
    void setMatrix( const SPARSE_MATRIX &M );
    void setMatrix( const SPARSE_MATRIX::SparseColumnMatrix &M );

    // Sets a matrix from a stored sparse matrix file
    //
    // format = 0 (SPARSE_MATRIX object written to file)
    // format = 1 (SPARSE_MATRIX::SparseColumnCopy object written to file)
    void setMatrix( const char *filename, int format );

    // Factorizes the current matrix.
    // If numeric is false, we will only perform a symbolic
    // factorization.
    void computeCholeskyFactor( bool numeric = true, bool symbolic = true,
                                int *userPerm = NULL );

    // Puts the factor in the given sparse matrix
    void getFactor( SPARSE_MATRIX &G );

    // Returns the actual CHOLMOD factor
    cholmod_factor *factor()
    {
      return _factor;
    }

    cholmod_common *common()
    {
      return &_common;
    }

    // Creates an empty factor
    void generateBlankFactor();

    // Wrapper for the CHOLMOD factor interface.  Basically just
    // a layer of abstraction in the name of good software
    // engineering, blah blah blah.
    //
    // This has the same fields as a supernodal CHOLMOD factor
    struct FactorWrapper {
      FactorWrapper()
        : _perm( NULL ),
          _x( NULL ),
          _super( NULL ),
          _pi( NULL ),
          _px( NULL ),
          _s( NULL )
      {
      }

      FactorWrapper( const FactorWrapper &wrapper )
        : _perm( NULL ),
          _x( NULL ),
          _super( NULL ),
          _pi( NULL ),
          _px( NULL ),
          _s( NULL )
      {
        copy( wrapper );
      }

      // FIXME: need copy constructor and = operator
      // to avoid doing something bad here

      virtual ~FactorWrapper()
      {
        clear();
      }

      FactorWrapper &operator=( const FactorWrapper &wrapper )
      {
        copy( wrapper );
        return *this;
      }

      void clear()
      {
        free( _perm );
        free( _x );
        free( _super );
        free( _pi );
        free( _px );
        free( _s );

        _perm = NULL;
        _x = NULL;
        _super = NULL;
        _pi = NULL;
        _px = NULL;
        _s = NULL;
      }

      void copy( const FactorWrapper &wrapper )
      {
        clear();

        FactorWrapper &factor = *this;

        factor._n = wrapper._n;
        factor._nsuper = wrapper._nsuper;
        factor._ssize = wrapper._ssize;
        factor._xsize = wrapper._xsize;
        factor._maxcsize = wrapper._maxcsize;
        factor._maxesize = wrapper._maxesize;

        // Construct copies of the fields in _factor
        if ( wrapper._perm != NULL )
        {
          factor._perm = (int *)malloc( factor._n * sizeof( int ) );
          factor._x = (double *)malloc( factor._xsize * sizeof( double ) );
          factor._super = (int *)malloc( (factor._nsuper + 1) * sizeof( int ) );
          factor._pi = (int *)malloc( (factor._nsuper + 1) * sizeof( int ) );
          factor._px = (int *)malloc( (factor._nsuper + 1) * sizeof( int ) );
          factor._s = (int *)malloc( factor._ssize * sizeof( int ) );

          memcpy( (void *)factor._perm, (void *)wrapper._perm,
                  factor._n * sizeof( int ) );

          if ( wrapper._x != NULL )
          {
            memcpy( (void *)factor._x, (void *)wrapper._x,
                    factor._xsize * sizeof( double ) );
          }
          else
          {
            factor._x = NULL;
          }

          memcpy( (void *)factor._super, (void *)wrapper._super,
                  (factor._nsuper + 1) * sizeof( int ) );
          memcpy( (void *)factor._pi, (void *)wrapper._pi,
                  (factor._nsuper + 1) * sizeof( int ) );
          memcpy( (void *)factor._px, (void *)wrapper._px,
                  (factor._nsuper + 1) * sizeof( int ) );
          memcpy( (void *)factor._s, (void *)wrapper._s,
                  factor._ssize * sizeof( int ) );
        }
      }

      size_t       _n;

      // The permutation applied to the matrix
      int         *_perm;

      // Actual numerical values for the factor
      double      *_x;

      // Supernodal stuff
      size_t       _nsuper;      // Number of supernodes
      size_t       _ssize;       // Size of the integer part of supernodes
      size_t       _xsize;       // Size of the numerical part
      size_t       _maxcsize;    // Size of largest update matrix (?)
      size_t       _maxesize;    // Max # of rows in a supernode, excl. tri. part

      int         *_super;       // size nsuper+1, first column in each supernode
      int         *_pi;          // size nsuper+1, pointers to integer patterns
      int         *_px;          // size nsuper+1, pointers to real parts
      int         *_s;           // size ssize, integer part of supernodes
    };

    // A separator node in a CHOLMOD nested dissection tree
    struct Separator {
      int          _size;

      // Range in the permuted matrix
      IndexRange   _columnRange;

      int          _numChildren;
    };

    typedef std::vector<Separator> SeparatorList;

    // Builds an ordered list of separators resulting from CHOLMOD'S
    // nexted dissection routine
    void nestedDissection( SeparatorList &separators, IntArray &permutation,
                           IntArray &separatorParents,
                           IntArray &separatorLevels,
                           IntArray &separatorMap );

    // Returns a wrapper for a CHOLMOD supernodal factorization
    bool getSupernodalFactor( FactorWrapper &factor );

    // Solves the factor system using the given right hand side.
    // If permuted == true, then this solves the system
    // (P^T * L) X = Y; that is, the permuted version of the factor system.
    void solveFactorSystem( const SPARSE_MATRIX &rhs, SPARSE_MATRIX &solution,
                            bool permuted = false );

    // Solves the matrix system (after factorization).
    void solveSystem( const VECTOR &rhs, VECTOR &solution );

    // Multiply the two given sparse matrices together.  The only
    // reason this isn't static is because we need something to
    // provide a cholmod_common environment.
    // C = A * B (or, optionally use transposes)
    void sparseMultiply( const SPARSE_MATRIX &A, bool transposeA,
                         const SPARSE_MATRIX &B, bool transposeB,
                         SPARSE_MATRIX &C );

    // Writes the permutation used by the system to the given file
    void writePermutation( const char *filename );

    // Extracts and AMD ordering for the current system
    void approximateMinimumDegree( IntArray &permutation );

    int systemSize() const { return _matrix->nrow; }

    Real cholmodMemoryUsage() const
    {
      return (Real)(_common.memory_usage) / 1048576.0;
    }

    Real factorTime() const { return _factorTimer.getTotalSecs(); }
    Real solveTime() const { return _solveTimer.getTotalSecs(); }
    Real testSolveTime() const { return _testSolveTimer.getTotalSecs(); }

    static void clearCHOLMODFactor( cholmod_factor *factor )
    {
      free( factor->Perm );
      free( factor->ColCount );
      free( factor->p );
      free( factor->i );
      free( factor->x );
      free( factor->z );
      free( factor->nz );
      free( factor->next );
      free( factor->prev );
      free( factor->super );
      free( factor->pi );
      free( factor->px );
      free( factor->s );

      factor->Perm = NULL;
      factor->ColCount = NULL;
      factor->p = NULL;
      factor->i = NULL;
      factor->x = NULL;
      factor->z = NULL;
      factor->nz = NULL;
      factor->next = NULL;
      factor->prev = NULL;
      factor->super = NULL;
      factor->pi = NULL;
      factor->px = NULL;
      factor->s = NULL;
    }

public:
    // Wrapper for either cholmod_sparse, or cholmod_factor
    struct CholmodWrapper {
      CholmodWrapper( const cholmod_factor *G )
        : n( G->n ),
          p( G->p ),
          i( G->i ),
          x( G->x )
      {
      }
      
      CholmodWrapper( const cholmod_sparse *A )
        : n( A->nrow ),
          p( A->p ),
          i( A->i ),
          x( A->x )
      {
      }

      // This version uses a symbolic supernodal factorization,
      // where indices are associated with supernodes.  This
      // is useful for finding the reach of a node in the
      // graph of supernodes.
      CholmodWrapper( const CHOLMOD_Supernodal *G )
        : n( G->compressedN() ),
          p( (void *)G->compressedP() ),
          i( (void *)G->compressedI() ),
          x( NULL ) // Symbolic only
      {
      }

      size_t n;

      const void *p;
      const void *i;
      const void *x;
    };

public:
    // Solves a sparse triangular system with a sparse right hand side.
    // xi is assumed to have size 2 * n, and x to have size n.
    // This is based entirely on cs_spsolve from [Davis, 2006], pg. 34
    static int                 spsolve( const CholmodWrapper *G,
                                        const cholmod_sparse *B,
                                        int *marker_list,
                                        int k, // The column of B
                                        int *xi,
                                        double *x,
                                        const int *pinv,
                                        int lo );

    // Computes the reach of a set of nodes, defined to be the
    // non-zero nodes in a given column of a sparse matrix
    // This is based entirely on the cs_reach algorithm
    // from Timothy Davis' "Direct Methods for Sparse Linear Systems"
    // (page 33)
    static int                 node_reach( const CholmodWrapper *G,
                                           const cholmod_sparse *B,
                                           int *marker_list,
                                           int k, // The column of B
                                           int *xi,
                                           const int *pinv );

    // Similar to the above, but with an explicitly provided
    // list of starting nodes
    static int                 node_reach( const CholmodWrapper *G,
                                           int *nodes, int numNodes,
                                           int *marker_list,
                                           int *xi,
                                           const int *pinv );

	protected:

  private:
    // Cleanup
    void clear();

    // Converts our sparse matrix format to CHOLMOD's.
    // returns copies of the matrix in both sparse and triplet
    // form (for easy modification later on).
    void                       allocateSparse( const SPARSE_MATRIX &M,
                                               cholmod_sparse * &A,
                                               cholmod_triplet * &T,
                                               bool symmetric = false,
                                               bool transpose = false );

    // Reads in a sparse matrix from a file and puts it in CHOLMOD's
    // internal format
    //
    // format = 0 (SPARSE_MATRIX object written to file)
    // format = 1 (SPARSE_MATRIX::SparseColumnCopy object written to file)
    void                       allocateSparse( const char *filename,
                                               int format,
                                               cholmod_sparse * &A,
                                               cholmod_triplet * &T,
                                               bool symmetric = false,
                                               bool transpose = false );

    // Different versions of the above for different formats
    void                       allocateSparseStandard(
                                                const char *filename,
                                                cholmod_sparse * &A,
                                                cholmod_triplet * &T,
                                                bool symmetric = false,
                                                bool transpose = false );

    void                       allocateSparseColumnCopy(
                                    const SPARSE_MATRIX::SparseColumnMatrix &S,
                                    cholmod_sparse * &A,
                                    cholmod_triplet * &T,
                                    bool symmetric = false,
                                    bool transpose = false );
                                                

    // Cleans up a sparse matrix
    void                       freeSparse( cholmod_sparse * &A );

    // Cleans up a sparse matrix in triplet form
    void                       freeTriplet( cholmod_triplet * &T );
    
    // Cleans up a sparse matrix factor
    void                       freeFactor( cholmod_factor * &G );

    // Given a nested dissection ordering, determines the level of
    // each separator in the nested dissection tree
    static void                getSeparatorLevels(
                                          const IntArray &separatorParents,
                                          IntArray &separatorLevels );

    // Recursive helper function for the above
    static void                getSeparatorLevels(
                                          int sep_idx,
                                          const IntArray &separatorParents,
                                          IntArray &separatorLevels );

    // Error handler for CHOLMOD
    static void                handleError( int status, const char *file,
                                            int line, const char *message );

    // Converts the CHOLMOD triplet format to our SPARSE_MATRIX
    // format
    static void                tripletToSparse( const cholmod_triplet *T,
                                                SPARSE_MATRIX &A );

    // Computes the reach of a node in the graph of the given
    // factor.  This is based entirely on the cs_dfs algorithm
    // from Timothy Davis' "Direct Methods for Sparse Linear Systems"
    // (page 33)
    static int                 dfs_node_reach( int j,
                                               const CholmodWrapper *G,
                                               int *marker_list,
                                               int top,
                                               int *xi,
                                               int *pstack,
                                               const int *pinv );

	private:
    // Workspace, parameters, etc.
    cholmod_common             _common;

    // Input matrix in sparse form
    cholmod_sparse            *_matrix;
    cholmod_triplet           *_matrixTriplet;

    // Matrix Cholesky factor
    cholmod_factor            *_factor;
    cholmod_sparse            *_factorSparse;
    cholmod_triplet           *_factorTriplet;

    Timer                      _factorTimer;
    Timer                      _solveTimer;
    Timer                      _testSolveTimer;

};

//////////////////////////////////////////////////////////////////////
// CHOLMOD_Environment class
//
// Provides functions for interacting with the CHOLMOD library.
// We will need one of these objects wherever we want to factorize
// a matrix.
//////////////////////////////////////////////////////////////////////
class CHOLMOD_Environment_Long {
	public:
		CHOLMOD_Environment_Long( int ordering = CHOLMOD_METIS );

		// Destructor
		virtual ~CHOLMOD_Environment_Long();

    // Sets a new matrix to factorize
    void setMatrix( const SPARSE_MATRIX &M );

    // Sets a matrix from a stored sparse matrix file
    //
    // format = 0 (SPARSE_MATRIX object written to file)
    // format = 1 (SPARSE_MATRIX::SparseColumnCopy object written to file)
    void setMatrix( const char *filename, int format );

    // Factorizes the current matrix.
    // If numeric is false, we will only perform a symbolic
    // factorization.
    void computeCholeskyFactor( bool numeric = true, bool symbolic = true,
                                long *userPerm = NULL );

    // Puts the factor in the given sparse matrix
    void getFactor( SPARSE_MATRIX &G );

    // Returns the actual CHOLMOD factor
    cholmod_factor *factor()
    {
      return _factor;
    }

    cholmod_common *common()
    {
      return &_common;
    }

    // Solves the factor system using the given right hand side.
    // If permuted == true, then this solves the system
    // (P^T * L) X = Y; that is, the permuted version of the factor system.
    void solveFactorSystem( const SPARSE_MATRIX &rhs, SPARSE_MATRIX &solution,
                            bool permuted = false );

    // Solves the matrix system (after factorization).
    void solveSystem( const VECTOR &rhs, VECTOR &solution );

    int systemSize() const { return _matrix->nrow; }

    Real cholmodMemoryUsage() const
    {
      return (Real)(_common.memory_usage) / 1048576.0;
    }

    Real factorTime() const { return _factorTimer.getTotalSecs(); }
    Real solveTime() const { return _solveTimer.getTotalSecs(); }
    Real testSolveTime() const { return _testSolveTimer.getTotalSecs(); }

    static void clearCHOLMODFactor( cholmod_factor *factor )
    {
      free( factor->Perm );
      free( factor->ColCount );
      free( factor->p );
      free( factor->i );
      free( factor->x );
      free( factor->z );
      free( factor->nz );
      free( factor->next );
      free( factor->prev );
      free( factor->super );
      free( factor->pi );
      free( factor->px );
      free( factor->s );

      factor->Perm = NULL;
      factor->ColCount = NULL;
      factor->p = NULL;
      factor->i = NULL;
      factor->x = NULL;
      factor->z = NULL;
      factor->nz = NULL;
      factor->next = NULL;
      factor->prev = NULL;
      factor->super = NULL;
      factor->pi = NULL;
      factor->px = NULL;
      factor->s = NULL;
    }

	protected:

  private:
    // Cleanup
    void clear();

    // Converts our sparse matrix format to CHOLMOD's.
    // returns copies of the matrix in both sparse and triplet
    // form (for easy modification later on).
    void                       allocateSparse( const SPARSE_MATRIX &M,
                                               cholmod_sparse * &A,
                                               cholmod_triplet * &T,
                                               bool symmetric = false,
                                               bool transpose = false );

    // Reads in a sparse matrix from a file and puts it in CHOLMOD's
    // internal format
    //
    // format = 0 (SPARSE_MATRIX object written to file)
    // format = 1 (SPARSE_MATRIX::SparseColumnCopy object written to file)
    void                       allocateSparse( const char *filename,
                                               int format,
                                               cholmod_sparse * &A,
                                               cholmod_triplet * &T,
                                               bool symmetric = false,
                                               bool transpose = false );

    // Different versions of the above for different formats
    void                       allocateSparseStandard(
                                                const char *filename,
                                                cholmod_sparse * &A,
                                                cholmod_triplet * &T,
                                                bool symmetric = false,
                                                bool transpose = false );

    void                       allocateSparseColumnCopy(
                                                const char *filename,
                                                cholmod_sparse * &A,
                                                cholmod_triplet * &T,
                                                bool symmetric = false,
                                                bool transpose = false );
                                                

    // Cleans up a sparse matrix
    void                       freeSparse( cholmod_sparse * &A );

    // Cleans up a sparse matrix in triplet form
    void                       freeTriplet( cholmod_triplet * &T );
    
    // Cleans up a sparse matrix factor
    void                       freeFactor( cholmod_factor * &G );

    // Error handler for CHOLMOD
    static void                handleError( int status, const char *file,
                                            int line, const char *message );

    // Converts the CHOLMOD triplet format to our SPARSE_MATRIX
    // format
    static void                tripletToSparse( const cholmod_triplet *T,
                                                SPARSE_MATRIX &A );

	private:
    // Workspace, parameters, etc.
    cholmod_common             _common;

    // Input matrix in sparse form
    cholmod_sparse            *_matrix;
    cholmod_triplet           *_matrixTriplet;

    // Matrix Cholesky factor
    cholmod_factor            *_factor;
    cholmod_sparse            *_factorSparse;
    cholmod_triplet           *_factorTriplet;

    Timer                      _factorTimer;
    Timer                      _solveTimer;
    Timer                      _testSolveTimer;

};

#endif
