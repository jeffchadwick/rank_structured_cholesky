//////////////////////////////////////////////////////////////////////
// FactorManager.h: Interface for the FactorManager class
//
//////////////////////////////////////////////////////////////////////

#ifndef FACTOR_MANAGER_H
#define FACTOR_MANAGER_H

#include <datastructure/WorkspaceManager.h>

#include <linearalgebra/DenseBlock.h>
#include <linearalgebra/MATRIX.h>
#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/SPARSE_MATRIX.h>
#include <linearalgebra/VECTOR.h>

#include <node/Supernode.h>

#include <ordering/NestedDissection.h>
#include <ordering/Ordering.h>

#include <solver/SolverError.h>

#include <util/IO.h>
#include <util/timer.h>

#include <SETTINGS.h>
#include <TYPES.h>

#include "config.h"

//////////////////////////////////////////////////////////////////////
// FactorManager class
//
// Manages a set of supernodes
//////////////////////////////////////////////////////////////////////
class FactorManager {
	public:
    // Constructor: Assumes that an ordering for the matrix
    // has been generated and is stored in the nodes array
    //
    // Optionally provide a list of diagonal blocks for each node in the
    // given ordering.  This is used for compression of blocks in a node's
    // main diagonal.
		FactorManager(
        const SPARSE_MATRIX::SparseColumnMatrix &A,
        //const std::vector<const NestedDissection::DissectionNode *> &nodes,
        const std::vector<Ordering::SupernodeSpecification> &nodes,
        const std::vector<std::vector<DenseBlock> > *diagonalBlocks = NULL,
        Real errorBoundMultiplier = ERROR_BOUND_MULTIPLIER,
        Real rankConstant = -1, Real diagonalRankConstant = -1 );

		// Destructor
		virtual ~FactorManager();

    void clear();

    // Builds a symbolic factor for this manager's system,
    // optionally extending the system with slack variables.
    void factorSymbolic( int maxBlockSize = -1,
                         CompressionType compressType = EXTENDED_VARIABLE,
                         bool useInteriorBlocks = false );

    // Computes the numerical factor for our system, assuming that
    // the symbolic factor has already been computed.
    // The pattern of A is assumed to match the pattern of the
    // stored matrix.  No error checking is done for this, however.
    //
    // If writeSize is positive, write all Schur complement blocks
    // exceeding that size to disk.
    //
    // If writeFactors is true, then write the low rank decompositions
    // for blocks to disk.
    //
    // compressedInteractions == true implies that we will be storing
    // low rank decompositions in place, rather than introducing
    // extended nodes
    //
    // decomposeTranspose == true implies that we will form an orthonormal
    // basis for 
    SolverErrorPtr factorNumeric(
                        const SPARSE_MATRIX::SparseColumnMatrix &A,
                        bool adaptive, Real tolerance = 0.01,
                        Real diagonalTolerance = 0.01,
                        int writeSize = -1, bool writeFactors = false,
                        ExtendedFactorType extFactorType = CHOLESKY,
                        CompressionType compressType = EXTENDED_VARIABLE,
                        DecompositionBlock blockType = DECOMPOSE_STANDARD,
                        DecompositionType decompType = DECOMPOSE_SCHUR );

    // Copies the off-diagonal slack variable data to a sparse matrix
    void copySlackOffDiagonal( SPARSE_MATRIX &S );

    // Check for any strange stuff happening in the factor
    bool checkFactor( int maxNode = -1 );

    // Counts the columns in the given list of nodes
    int nodeColumns( const IntArray &nodeList );

    void writeTimings();

    // Computes numerical solution to system using the computed factor
    inline void solveSystem( const SPARSE_MATRIX::SparseColumnMatrix &A,
                             VECTOR &x, int startNode = 0 )
    {
#if 0
      applyLDLbackwardPermutation( x, startNode );
#endif
      if ( _useInteriorBlocks ) {
        forwardSolve( A, x );
      } else {
        forwardSolve( x, startNode );
      }
#if 0
      LDLdiagonalSolve( x, startNode );
#endif
      if ( _useInteriorBlocks ) {
        backwardSolve( A, x );
      } else {
        backwardSolve( x, startNode );
      }
#if 0
      applyLDLforwardPermutation( x, startNode );
#endif
    }

    inline void solveSystem( VECTOR &x, int startNode = 0 )
    {
      solveSystem( _A, x, startNode );
    }

    int factorSystemSize() const
    {
      return _factor.back().endColumn() + 1;
    }

    int numExtendedNodes() const
    {
      return _numExtendedNodes;
    }

    int numNodes() const
    {
      return _factor.size();
    }

    long int slackUsage() const
    {
      return ( _slackDataSz - _availableSlackData );
    }

    Real slackUsageMB() const
    {
      return (Real)slackUsage() * 8.0 / 1024.0 / 1024.0;
    }

    const std::vector<InteriorBlock> &interiorBlocks() const
    {
      return _interiorBlocks;
    }

    int nodeColumns( int node_idx ) const
    {
      return _factor.at( node_idx ).numColumns();
    }

    //////////////////////////////////////////////////////////////////////
    // Diagnostic functions
    //////////////////////////////////////////////////////////////////////

    // Once all standard nodes have been handled, call this to
    // build a symbolic representation of the Schur complement formed
    // by eliminating just the standard nodes in the system.
    //
    // Write this sparse matrix to the given file.
    void writeExtendedSymbolicSchurComplement( const char *filename );

    // Writes the current contents of the extended variable block
    // to a sparse matrix
    void writeExtendedBlock( const char *filename );

    // For an extended block factored via an LDL factorization, writes
    // the diagonal and permutation matrices to disk
    void writeExtendedLDLData( const char *filename );

    // Prints out usage statistics for both off diagonal extended
    // interactions occuring in the main part of the factorization,
    // and those occuring in the extended variable space
    void printExtendedUsage();

    // Estimates the storage used to just store the diagonal
    // parts of the factor between nodes with sparse blocks
    void printUncompressedDiagonalStorage();

    void printInteriorBlocks() const;

    // Determines the storage used to just store the diagonal parts
    // of interior blocks
    void printInteriorBlockDiagonalUsage() const;

    // Extract a sub-matrix of the factor
    void extractSubMatrix( const IntArray &nodeList, MATRIX &subMatrix ) const;

    // Extract a sub-matrix of the factor in sparse format
    void extractSubMatrix( const IntArray &nodeList,
                           SPARSE_MATRIX &subMatrix ) const;

    // Extracts an interior block interaction matrix
    void extractInteriorBlockInteraction( int interior_block_idx,
                                          int interaction_idx,
                                          MATRIX &subMatrix );

  public:
    static const Real RANK_CONSTANT;
    static const Real DIAGONAL_RANK_CONSTANT;

	protected:

	private:
    // Helper function for factorNumeric
    //
    // Initializes slack variable data for all extended nodes
    void initExtendedNodes( int start_node );

    // Helper function for factorNumeric
    //
    // Builds the schur complement in extended variable space
    void buildExtendedSchurComplement( int numStandardNodes );

    // Initializes workspace data according to the size of the system
    // to be factored.
    //
    // TODO: Doesn't handle extended systems yet
    void initWorkspace();

    // Initialize workspaces for forming low-rank decompositions.
    //
    // TODO: For now, consider the fixed rank problem; that is,
    //       we will determine the desired rank for a block using
    //       the blockRank function (instead of actually trying to
    //       approximate to some tolerance)
    void initDecompositionWorkspace();

    // Walks through the structure of a numerical factorization to
    // count the number of descendents each node expects to encounter
    void countDescendents( IntArray &descendentCounts ) const;

    // Counts the number of low rank block stored in the main diagonal
    // of each node
    void countLowRankDiagonalBlocks( IntArray &blockCounts ) const;

    // Estimates the total extra data needed for decomposing
    // low-rank blocks in the factor
    long int estimateExtraStorage();

    // Sets size estimates on all extended nodes
    void estimateSlackVariableSizes();
    void clearSlackVariableEstimates();

    // Clears all workspaces, preparing them for numerical factorization
    void clearWorkspace();

    // Gets ranks for each low rank block in the given node
    void setBlockRanks( const Supernode &node );

    // Applies diagonal updates to the given node using all current
    // interior block descendents (assumes that findInteriorBlockDescendents
    // has been called)
    void applyInteriorBlockDiagonalUpdate(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  Supernode &node );

    // Applies diagonal updates to the given node using an interior block
    void applyInteriorBlockDiagonalUpdate(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  Supernode &node, int block_idx );

    // Same as the function above, but applies this for a specific block
    // in this node's diagonal.  Also, places the result in a designated
    // workspace, rather than the node's own data.
    void constructInteriorBlockDiagonalUpdate(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  Supernode &node,
                                  const DenseBlock &factorBlock,
                                  // Interior block index
                                  int block_idx,
                                  Real *schurComplement );

    // Finds the ancestor interaction between this interior block and
    // the given node
    int findInteriorBlockAncestorInteraction( const InteriorBlock &block,
                                              const Supernode &node );

    int findInteriorBlockAncestorInteraction( int block_idx, int node_idx ) {
      return findInteriorBlockAncestorInteraction( _interiorBlocks[ block_idx ],
                                                   _factor[ node_idx ] );
    }

    // The number of non-zeros stored in an interior block interaction,
    // assuming it is explicitly constructed.
    // If a rowRange is provided, then we can fill in an array of
    // local row ranges that need to be considered in each node interaction.
    size_t getInteriorBlockInteractionSize(
                          const InteriorBlock &block,
                          const InteriorBlock::Interaction &blockInteraction,
                          const IndexRange &rowRange = IndexRange(-1, -1),
                          RangeArray *nodeRowRanges = NULL );

    size_t getInteriorBlockInteractionSize(
                          int block_idx, int interaction_idx,
                          const IndexRange &rowRange = IndexRange(-1,-1),
                          RangeArray *nodeRowRanges = NULL )
    {
      return getInteriorBlockInteractionSize(
              _interiorBlocks[ block_idx ],
              _interiorBlocks[ block_idx ]._interactions[ interaction_idx ],
              rowRange, nodeRowRanges );
    }

    // Initializes a workspace for explicitly expanding the contents
    // of an interior block interaction.  Optionally provide local
    // node row ranges if we are restricting ourselves to a particular
    // block row of this interaction.
    void initializeInteriorBlockInteractionWorkspace(
                          const InteriorBlock &block,
                          const InteriorBlock::Interaction &blockInteraction,
                          size_t interactionSize,
                          Real *workspaceBase,
                          vector<Real *> &workspacePointers,
                          IntArray &interactionMap,
                          RangeArray *nodeRanges = NULL );

    // Copies sparse matrix entries from the given 'interaction' of 'blockNode'
    // (which interacts with 'node') over the given row range in to the
    // provided workspace.
    //
    // fullRowRange is the range of rows we wish to consider from the
    // sparse matrix (relative to the starting column of node).
    //
    // localRowRange is the range of rows from 'interaction' that we
    // need to worry about.
    void copyInteractionSparseMatrix(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              const Supernode &node,
                              const Supernode &blockNode,
                              const SupernodeInteraction &interaction,
                              Real *workspace,
                              IndexRange fullRowRange = IndexRange(-1,-1),
                              IndexRange localRowRange = IndexRange(-1,-1) );

    // Pushes forward fill-in when explicitly constructing the contents of
    // an interior block interaction.
    //
    // Specifically, pushes fill from 'blockNode' to all ancestors.
    // The data for each node is assumed to be stored in the provided
    // workspace pointers.  'node' is the node that this particular
    // block interaction is interacting with
    //
    // Optionally, fullRowRange specifies the block row of this interior
    // block interaction that we are trying to form.
    //
    // If fullRowRange is specified then you should also specify a list of local
    // row ranges for each sub-interaction in the interior block interaction.
    void pushInteriorBlockInteractionFill(
                              const InteriorBlock &block,
                              const InteriorBlock::Interaction &blockInteraction,
                              const Supernode &node,
                              const Supernode &blockNode,
                              const SupernodeInteraction &nodeInteraction,
                              const vector<Real *> &workspaces,
                              const IntArray &interactionMap,
                              Real *multWorkspace,
                              IndexRange fullRowRange = IndexRange(-1,-1),
                              RangeArray *localRowRanges = NULL );

    // Constructs a diagonal update matrix from 'descendent' to 'node' over
    // the given sub-block of 'node's diagonal and subtracts it from the
    // given Schur complement.
    //
    // This will be used for explicit diagonal block formation/compression
    void constructNodeDiagonalUpdate( const Supernode &node,
                                      const Supernode &descendent,
                                      const DenseBlock &factorBlock,
                                      Real *schurComplement );

#if 0
    // Initializes workspaces for interaction multiplication.
    // This basically just means stacking a bunch of random matrices
    // on top of each other in a work space.
    void initMultiplyWorkspace( const Supernode &node,
                                int fixedRank = 0 );
#endif

    // Multiplies data in _decompMultTransWorkspace on the left
    // with each low rank interaction matrix for the current node.
    // Puts the result in _decompWorkspace.
    //
    // (ie. half a step of power iteraction)
    //
    // NOTE: inputData may change if we are using the factor, since
    // we will need to apply a triangular solve to it.
    void interactionMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              bool useFactor,
                              Real *inputData,
                              Real *outputData,
                              int fixedRank );

    // Same as the above, but multiplies the entire off-diagonal block
    // all at once
    void interactionMultiply_allRows(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              bool useFactor,
                              Real *inputData,
                              Real *outputData,
                              int fixedRank );

    // Processes normal descendent multiplications when forming a low-rank
    // interaction
    void interactionMultiply_nodeDescendents(
                                    const Supernode &node,
                                    const SupernodeInteraction &interaction,
                                    const PairArray &descendents,
                                    PairArray &dataSizes,
                                    int rank,
                                    const Real *multInput,
                                    Real *outputData );

    // Overloaded version of the above function which is assumes that
    // we are simultaneously decomposing the entire off-diagonal block
    // of node
    void interactionMultiply_nodeDescendents(
                                    const Supernode &node,
                                    PairArray &dataSizes,
                                    int rank,
                                    const Real *multInput,
                                    Real *outputData,
                                    bool transpose );

    // Processes interior block descendent multiplications when forming
    // a low-rank interaction
    void interactionMultiply_blockDescendents(
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const Supernode &node,
                                    const SupernodeInteraction &interaction,
                                    int rank,
                                    const Real *multInput,
                                    Real *outputData,
                                    IntArray &nodeList,
                                    bool transpose );

    // Same as the above, but assumes we are forming an off-diagonal
    // decomposition for node all at once
    void interactionMultiply_blockDescendents(
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const Supernode &node,
                                    int rank,
                                    const Real *multInput,
                                    Real *outputData,
                                    IntArray &nodeList,
                                    bool transpose );

#if 0
    // Version of interaction multiply that uses the interior block
    // formulation
    //
    // Assumes that findInteriorBlockDescendents has been called
    void interiorBlockInteractionMultiply(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              const Real *inputData = NULL,
                              Real *outputData = NULL,
                              int fixedRank = 0 );
#endif

    // Similar to the above function, but does the multiplication with
    // a sub-matrix of the main diagonal (indexed by block_idx)
    void diagonalBlockMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node,
                                bool useFactor,
                                int block_idx,
                                int rank,
                                const Real *inputData,
                                Real *outputData,
                                bool write = false );

    // Processes normal descendent multiplications when forming a low-rank
    // diagonal block
    void diagonalBlockMultiply_nodeDescendents( const Supernode &node,
                                                const DenseBlock &block,
                                                int block_idx,
                                                PairArray &dataSizes,
                                                int rank,
                                                const Real *multInput,
                                                Real *outputData );

    // Processes interior block descendent multiplications when forming a
    // low-rank diagonal block
    void diagonalBlockMultiply_blockDescendents(
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const Supernode &node,
                                    const DenseBlock &block,
                                    int block_idx,
                                    int rank,
                                    const Real *multInput,
                                    Real *outputData,
                                    IntArray &nodeList,
                                    bool transpose );

    // Multiplies data in _decompWorkspace on the left with the
    // transpose of each low rank itneraction matrix for the current node.
    // Puts the result in _decompMultTransWorkspace.
    //
    // (half a step of power iteration)
    void interactionTransMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                   Supernode &node,
                                   bool useFactor,
                                   const Real *inputData,
                                   Real *outputData,
                                   int fixedRank );

    // Same as the above, but multiplies with the entire off-diagonal
    // all at once
    void interactionTransMultiply_allRows(
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    Supernode &node,
                                    bool useFactor,
                                    const Real *inputData,
                                    Real *outputData,
                                    int fixedRank );

    // Processes normal descendent multiplications when forming a low-rank
    // interaction
    void interactionTransMultiply_nodeDescendents(
                                    const Supernode &node,
                                    const SupernodeInteraction &interaction,
                                    const PairArray &descendents,
                                    PairArray &dataSizes,
                                    int rank,
                                    const Real *multInput,
                                    Real *outputData );

#if 0
    // Overloaded version of the above function which is assumes that
    // we are simultaneously decomposing the entire off-diagonal block
    // of node
    void interactionTransMultiply_nodeDescendents(
                                    const Supernode &node,
                                    PairArray &dataSizes,
                                    int rank,
                                    const Real *multInput,
                                    Real *outputData,
                                    bool transpose );
#endif

#if 0
    // Processes interior block descendent multiplications when forming
    // a low-rank interaction
    void interactionTransMultiply_blockDescendents(
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const Supernode &node,
                                    const SupernodeInteraction &interaction,
                                    int rank,
                                    const Real *multInput,
                                    Real *outputData,
                                    IntArray &nodeList );
#endif

    // Similar to the above function, but does the multiplication with
    // a sub-matrix of the main diagonal (indexed by block_idx)
    void diagonalBlockTransMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                     Supernode &node,
                                     bool useFactor,
                                     int block_idx,
                                     int rank,
                                     const Real *inputData,
                                     Real *outputData,
                                     bool write = false );

    // Processes normal descendent multiplications when forming a low-rank
    // diagonal block
    void diagonalBlockTransMultiply_nodeDescendents( const Supernode &node,
                                                     const DenseBlock &block,
                                                     int block_idx,
                                                     PairArray &dataSizes,
                                                     int rank,
                                                     const Real *multInput,
                                                     Real *outputData );

    // Forms the second part of a low rank interaction by multiplying
    // the interaction on the left with its basis transposed.
    //
    // transposeBasis: Whether or not to transpose the basis when
    //                 multiplying (depends on how it is stored)
    void fullBasisProject( const SPARSE_MATRIX::SparseColumnMatrix &A,
                           Supernode &node,
                           bool writeFactors,
                           int fixedRank,
                           bool transposeBasis,
                           bool transpose,
                           bool useFactor,
                           OffDiagonalCompressionType basisType );

    // Version of fullBasisProject for decomposing lower triangular
    // blocks directly
    void fullBasisProject_standard( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    Supernode &node,
                                    bool writeFactors,
                                    int fixedRank,
                                    bool transposeBasis,
                                    bool useFactor );

    // Version of the above to use when decomposing all rows at once
    void fullBasisProject_standard_allRows(
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    Supernode &node,
                                    bool writeFactors,
                                    int fixedRank,
                                    bool transposeBasis,
                                    bool useFactor );

    // Version of fullBasisProject for decomposing transposes of lower
    // triangular blocks (ie. decomposing upper triangular blocks)
    void fullBasisProject_transpose( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                     Supernode &node,
                                     bool writeFactors,
                                     int fixedRank,
                                     bool transposeBasis,
                                     bool useFactor );

    // Version to use when decomposing all rows at once
    void fullBasisProject_transpose_allRows(
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    Supernode &node,
                                    bool writeFactors,
                                    int fixedRank,
                                    bool transposeBasis,
                                    bool useFactor );

    // Reprojection procedure to guarantee that the modified system
    // remains positive definite.
    void reprojectDecomposition( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                 Supernode &node );

    // Same as the above function, b
    void diagonalBasisProject( const SPARSE_MATRIX::SparseColumnMatrix &A,
                               Supernode &node, int block_idx,
                               bool writeFactors,
                               bool transpose,
                               bool useFactor,
                               int fixedRank = 0 );

    // Version of diagonalBasisProject for decomposing lower triangular
    // blocks directly
    void diagonalBasisProject_standard(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node, int block_idx,
                              bool writeFactors,
                              bool useFactor,
                              int fixedRank = 0 );

    // Version of diagonalBasisProject for decomposing transposes of
    // lower triangular blocks (ie. decomposing upper triangular blocks)
    void diagonalBasisProject_transpose(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node, int block_idx,
                              bool writeFactors,
                              bool useFactor,
                              int fixedRank = 0 );

    void activateAllBlocks()
    {
      for ( int block_idx = 0; block_idx < _lowRankBlocks.size(); block_idx++ )
      {
        _lowRankBlockActive[ block_idx ] = true;
      }
    }

    // Performs a forward solve (L x = b) using the computed
    // numerical factor
    void forwardSolve( VECTOR &x, int startNode = 0 );

    // Forward solve function when we are using interior blocks
    void forwardSolve( const SPARSE_MATRIX::SparseColumnMatrix &A,
                       VECTOR &x );

    // Same as the above, but a backward solve (L' x = b)
    void backwardSolve( VECTOR &x, int startNode = 0 );

    // Backward solve function when we are using interior blocks
    void backwardSolve( const SPARSE_MATRIX::SparseColumnMatrix &A,
                        VECTOR &x );

    // Helper functions for running solves with interior blocks
    void interiorBlockForwardSolve(
                        const SPARSE_MATRIX::SparseColumnMatrix &A,
                        VECTOR &x, int block_idx,
                        Real *workspace );

    // Helper function for backward solves with an interior block
    void interiorBlockBackwardSolve(
                        const SPARSE_MATRIX::SparseColumnMatrix &A,
                        VECTOR &x, int block_idx,
                        Real *workspace );

    //////////////////////////////////////////////////////////////////////
    // Additional solver routines for systems which have nodes factored
    // via an LDL factorization
    //////////////////////////////////////////////////////////////////////

    void LDLdiagonalSolve( VECTOR &x, int startNode = 0 );

    void applyLDLforwardPermutation( VECTOR &x, int startNode = 0 );

    void applyLDLbackwardPermutation( VECTOR &x, int startNode = 0 );

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    // Solves a sub-system specified by the given node list.
    // Multiple right hand sides are supported.
    //
    // This is here to support the "interior block" machinery
    void forwardSubSolve( const IntArray &nodeList, int nRHS,
                          Real *rhs, Real *workspace = NULL );

    // With a node range, rather than a list
    void forwardSubSolve( IndexRange nodeRange, int nRHS,
                          Real *rhs, Real *workspace = NULL );

    // Same as the above, but for backward solves
    void backwardSubSolve( const IntArray &nodeList, int nRHS,
                           Real *rhs, Real *workspace = NULL );

    // With a node range rather than a list
    void backwardSubSolve( IndexRange nodeRange, int nRHS,
                           Real *rhs, Real *workspace = NULL );

#if 0
    // Returns ranks for a set of low rank blocks
    void blockRanks( const Supernode &node, const IntArray &lowRankBlocks,
                     IntArray &ranks );
#endif

    // Calculates slack variable memory usage in both the
    // off-diagonal part, and the Schur complement for the
    // slack variables (prior to elimination).
    void calculateSlackUsage( int maxBlockSize ) const;

    // Given a set of blocks to be decomposed for this node,
    // figure out what the blocks should actually be equal to
    // and write the results to disk.
    void writeLowRankBlocks( const SPARSE_MATRIX::SparseColumnMatrix &A,
                             Supernode &node );

    // Writes the given low rank diagonal block to disk
    void writeLowRankDiagonalBlock( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    Supernode &node, int block_idx );

    // Breaks down standard node memory usage in terms of which
    // nodes interactions actually act on
    void writeStandardMemoryUsage( int maxBlockSize );

    // Sets the block size for the current decomposition iteration
    //
    // basisType determines whether we are decomposing interactions
    // independently, or all using the same basis
    void setDecompositionBlockSizes( const Supernode &node, int iteration,
                                     bool adaptive,
                                     OffDiagonalCompressionType basisType );

    // Size of a block to use during low rank decomposition, given the
    // number of non-zero columns in the block, and the current decomposition
    // iteration
    int decompositionBlockSize( int nonZeroColumns, int iteration,
                                bool hasDescendents = true,
                                bool isDiagonalBlock = false );

    // For fixed (non-adaptive) rank decomposition
    int fixedDecompositionBlockSize( int nonZeroColumns,
                                     int nRows, int nCols,
                                     bool hasDescendents = true,
                                     bool isDiagonalBlock = false );

    int blockRank( int nRows, int nCols ) const;
    int diagonalBlockRank( int nRows, int nCols ) const;

    int maxBlockRank( int nRows, int nCols );
    int diagonalMaxBlockRank( int nRows, int nCols );
    
    int fixedBlockRank( int nRows, int nCols );
    
    //////////////////////////////////////////////////////////////////////
    // Adaptive low rank decomposition code
    //////////////////////////////////////////////////////////////////////

    // Error estimator based on power iteration.  For each
    //
    // From [Liberty et al. 2007] "Randomized algorithms for..."
    // boundMultiplier = 10.0 in this paper.
    void blockErrorEstimate( const SPARSE_MATRIX::SparseColumnMatrix &A,
                             Supernode &node, int powerIterations,
                             Real boundMultiplier,
                             bool transpose,
                             bool useFactor );

    // Error estimator used when compressing blocks in a node's diagonal
    void diagonalBlockErrorEstimate(
                             const SPARSE_MATRIX::SparseColumnMatrix &A,
                             Supernode &node, int block_idx,
                             int powerIterations, Real boundMultiplier,
                             bool transpose, bool useFactor );

    // Checks error for each block currently being decomposed in node.
    // Sets the active flag (in _lowRankBlockActive) to false for any
    // block A satisfying the error bound
    //    || A - Q Q' A ||
    //   ------------------ < tolerance
    //        || A ||
    //
    // Also, for those nodes, copy the contents of _decompWorkspace
    // to _basisWorkspace to prepare for the next iteration
    //
    // Returns false if no blocks are active anymore
    bool updateActiveBlocks( const SPARSE_MATRIX::SparseColumnMatrix &A,
                             Supernode &node,
                             int blockSize, Real tolerance,
                             bool transpose,
                             bool useFactor );

    // Similar to the above, but only acts on a single low rank block
    // from the node's diagonal
    bool updateActiveDiagonalBlock( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    Supernode &node, int block_idx,
                                    int blockSize, Real tolerance,
                                    bool transpose,
                                    bool useFactor );

    // Takes random Gaussian distributed values from _testMatrix
    // and forms a set of normalized vectors (according to column
    // sizes in low rank blocks) in _vectorWork1
    void buildRandomTestVectors( Supernode &node, int nColsTotal,
                                 bool singleVector = false );

    // Assigns workspace data to each of a node's low rank blocks.
    // This is used to incrementally build a basis for the column
    // space of the block.
    //
    // This treats the entries of _blockRanks as the maximum
    // size of any block.
    //
    // transpose == true if we are decomposing the transpose of an
    // interaction, rather than the interaction itself
    void assignBasisWorkspace( Supernode &node,
                               bool transpose,
                               OffDiagonalCompressionType basisType );

    // Generates a block of random numbers and multiplies each
    // low rank block by this block (possibly using power iteration).
    //
    // That is, forms the matrix (A * A')^q * A * G for each block to
    // be decomposed, where q is the number of power iterations
    //
    // blockSize: The number of columns in G
    void blockRandomMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node, 
                              int blockSize, int powerIterations,
                              bool transpose,
                              bool useFactor );

    // Similar to the above, but used to decompose the full off-diagonal
    // block of a node all at once
    void blockRandomMultiply_allRows(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              int blockSize, int powerIterations,
                              bool transpose, bool useFactor );

    // Generates a block of random numbers and multiplies a single low
    // rank block from the given node's diagonal by that matrix
    //
    // ie. forms the matrix (A * A')^q * A * G for the block to be
    // decomposed, where q is the number of power iterations.
    //
    // blockSize: The number of columns in G
    void diagonalBlockRandomMultiply(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node, int block_idx,
                              int blockSize, int powerIterations,
                              bool transpose, bool useFactor );

    // Counts the number of active low rank rows (that is, rows that
    // are part of blocks that still need columns added to their
    // bases).
    int countActiveLowRankRows( Supernode &node );

    // Initialization step for adaptive low rank decomposition
    //
    // transpose == true if we want to form the decomposition of block
    // transposes, rather than the blocks themselves
    void initOffDiagonalDecomposition(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  Supernode &node,
                                  int powerIterations,
                                  bool adaptive,
                                  bool transpose,
                                  bool useFactor,
                                  OffDiagonalCompressionType basisType );

    // Helper function for initOffDiagonalDecomposition
    void copyAllBases( const Supernode &node,
                       const Real *baseInputData,
                       Real *baseOutputData,
                       bool transpose,
                       OffDiagonalCompressionType basisType );

    // Initialization for adaptive decomposition of a diagonal low rank block
    void initDiagonalDecomposition( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    Supernode &node, int block_idx,
                                    int powerIterations,
                                    int &blockSize,
                                    bool adaptive,
                                    bool transpose,
                                    bool useFactor );

    // Adaptive decomposition iteration
    //
    // transpose == true if we want to decompose the transposes of
    // lower triangular blocks, rather than the blocks themselves
    //
    // useFactor == true if we are trying to directly decompose
    // factor blocks, rather than schur complements (which means that
    // the factor must be included when doing multiplication)
    void offDiagonalDecomposition( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                   Supernode &node,
                                   int powerIterations,
                                   Real tolerance,
                                   bool adaptive,
                                   bool transpose,
                                   bool useFactor,
                                   OffDiagonalCompressionType basisType );

    // Fixed-rank off-diagonal decomposition
    //
    // transpose == true if we want to decompose the transposes of
    // lower triangular blocks, rather than the blocks themselves
    //
    // useFactor == true if we are trying to directly decompose
    // factor blocks, rather than schur complements (which means that
    // the factor must be included when doing multiplication)
    void fixedOffDiagonalDecomposition(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  Supernode &node,
                                  int powerIterations,
                                  bool transpose,
                                  bool useFactor,
                                  OffDiagonalCompressionType basisType );

    // Adaptive off diagonal decomposition
    //
    // transpose == true if we want to decompose the transposes of
    // lower triangular blocks, rather than the blocks themselves
    //
    // useFactor == true if we are trying to directly decompose
    // factor blocks, rather than schur complements (which means that
    // the factor must be included when doing multiplication)
    void adaptiveOffDiagonalDecomposition(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  Supernode &node,
                                  int powerIterations,
                                  Real tolerance,
                                  bool transpose,
                                  bool useFactor,
                                  OffDiagonalCompressionType basisType );

    // Adaptive decomposition for low rank blocks on the main diagonal of
    // a node
    void diagonalDecomposition( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node,
                                int powerIterations,
                                Real tolerance,
                                bool adaptive,
                                bool transpose,
                                bool useFactor );

    // Factorization of a node's diagonal with in-place block compression
    SolverErrorPtr factorCompressedDiagonal(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  Supernode &node,
                                  int powerIterations,
                                  Real tolerance,
                                  bool adaptive,
                                  bool transpose,
                                  bool useFactor );

    // Adaptive decomposition for a single low rank block in the node's
    // main diagonal
    void diagonalDecomposition( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Supernode &node,
                                int block_idx,
                                int powerIterations,
                                Real tolerance,
                                bool adaptive,
                                bool transpose,
                                bool useFactor );

    // Fixed-rank diagonal decomposition
    void fixedDiagonalDecomposition( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                     Supernode &node,
                                     int block_idx,
                                     int powerIterations,
                                     int blockSize,
                                     bool transpose,
                                     bool useFactor );

    // Adaptive diagonal decomposition
    void adaptiveDiagonalDecomposition(
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    Supernode &node,
                                    int block_idx,
                                    int powerIterations,
                                    Real tolerance,
                                    int blockSize,
                                    bool transpose,
                                    bool useFactor );

    // Decomposes off-diagonal block block_idx, *AND* all of its children
    // by ecplicitly forming the 2x2 diagonal block for which block[ block_idx ]
    // is the off-diagonal part.  Next, this function directly forms low-rank
    // approximations for pieces of this block.
    SolverErrorPtr explicitDiagonalDecomposition(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              Supernode &node,
                              int block_idx,
                              int powerIterations,
                              bool transpose,
                              bool useFactor );

    // For each low rank block i to be decomposed, using the basis Q
    // stored in _basisStorage[ i ], project the contents of _basisWorkspace
    // in to the space orthogonal to Q.
    void basisOrthogonalProjection( Supernode &node, int blockSize,
                                    bool transpose );

    // Similar to the above, but only operates on a basis for a single
    // low rank block from the node's diagonal
    void diagonalBasisOrthogonalProjection( Supernode &node, int block_idx,
                                            int blockSize,
                                            bool transpose );

    // Append bases stored in _basisStorage based on the contents
    // of _basisWorkspace
    void appendOffDiagonalBases(
                      Supernode &node, int blockSize, int iteration,
                      const SPARSE_MATRIX::SparseColumnMatrix *A,
                      bool transpose = false );

    // Appends basis data when building a single off-diagonal decomposition
    void appendOffDiagonalBasis(
                      Supernode &node, int blockSize, int iteration,
                      const SPARSE_MATRIX::SparseColumnMatrix *A,
                      bool transpose = false );

    // Appends to a single basis for a block on the given node's diagonal
    bool appendDiagonalBasis( Supernode &node, int block_idx, int blockSize,
                              bool transpose );

    // Copies finalized adaptive bases to _decompWorkspace in preparation
    // for the second pass of low-rank decomposition
    void copyFinalBases( bool singleBasis, bool transpose );

    //////////////////////////////////////////////////////////////////////
    // For sparsification of the main diagonal block in a supernode
    //////////////////////////////////////////////////////////////////////

    // For each descendent of the node currently being handled identify
    // row/column ranges of the descendent that must be considered to
    // form each low-rank block in the current node's main diagonal.
    void findDiagonalBlockRanges( const Supernode &currentNode );

    // Gets maximum row/column counts, ranks, block sizes, etc.
    // for all low rank diagonal blocks in the given node
    void diagonalBlockSizes( const Supernode &node,
                             int &maxRows, int &maxColumns, int &maxRank );

    // Yes, this is annoying
    void cacheBlockRanks();

    void initTimers();

    // Assuming that sizes have been set in the extended part of the
    // factorization, compute a new ordering for this part of the
    // factorization using a weighted minimum degree ordering
    void buildExtendedPermutation( IntArray &permutation,
                                   IntArray &inversePermuation );

    // Once all standard nodes have been handled, call this to
    // build a symbolic representation of the Schur complement formed
    // by eliminating just the standard nodes in the system.
    //
    // Write this sparse matrix to the given file.
    void buildExtendedSchurComplement(
                          std::vector<std::set<int> > &extendedInteractions );

    //////////////////////////////////////////////////////////////////////
    // Code for using interior blocks
    //////////////////////////////////////////////////////////////////////

    // Identify "interior" block ranges; that is, ranges of supernodes
    // whose off-diagonal will not be explicitly stored
    void identifyInteriorBlocks();

    // Builds interaction sets for interior blocks
    void buildInteriorBlockInteractions();

    // Helper function for the above
    void buildSingleInteriorBlockInteractions( InteriorBlock &block );

    // Debugging routine to make sure that the interior blocks are
    // working the way they should
    void verifyInteriorBlocks() const;

#if 0
    // Node list is the list of supernodes to use in this multiplication,
    // relative to the entire matrix.
    //
    // subNodeList is the list of nodes (of the same size) relative to
    // the particular interior block
    //
    // Performs the multiplication necessary to form a low rank block
    // depending on the given interior block, without explicitly forming
    // the interior block.
    //
    // Optionally restrict to a row and column range (relative to the
    // dimensions of the product matrix)
    void interiorBlockSchurMultiply(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const IntArray &nodeList,
                                int interior_block_idx,
                                int interaction_idx1, int interaction_idx2,
                                const Real *G,
                                int nRowsG, int nColsG,
                                Real *multWorkspaceInitial,
                                Real *multWorkspaceFinal,
                                IndexRange rowRange = IndexRange ( -1, -1 ),
                                IndexRange columnRange = IndexRange( -1, -1 ),
                                bool transpose = false, bool left = false );
#endif

    // Node list is the list of supernodes to use in this multiplication,
    // relative to the entire matrix.
    //
    // subNodeList is the list of nodes (of the same size) relative to
    // the particular interior block
    //
    // Performs the multiplication necessary to form a low rank block
    // depending on the given interior block, without explicitly forming
    // the interior block.
    //
    // inverseRowList maps rows in the "full" matrix that we are forming to
    // those in a matrix in which some of these rows may be excluded.
    // We also specify the actual number of rows that we wish to consider.
    //
    // Optionally restrict to a row and column range (relative to the
    // dimensions of the product matrix).  Note that if meaningful row/column
    // ranges are provided, this will ignore the provided inverse row list and
    // assume that all rows in the desired range are required.
#if 0
    void interiorBlockSchurMultiply(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const IntArray &nodeList,
                                int interior_block_idx,
                                int interaction_idx1, int interaction_idx2,
                                const Real *G,
                                int nRowsG, int nColsG,
                                const IntArray &inverseRowList, int nRows,
                                Real *multWorkspaceInitial,
                                Real *outputData,
                                IndexRange rowRange = IndexRange ( -1, -1 ),
                                IndexRange columnRange = IndexRange( -1, -1 ),
                                bool transpose = false, bool left = false );
#endif

    // New version to work on
    void interiorBlockSchurMultiply(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                IntArray &nodeList,
                                int interior_block_idx,
                                int interaction_idx1, int interaction_idx2,
                                const Real *G,
                                int nRowsG, int nColsG,
                                const IntArray &inverseRowList, int nRows,
                                Real *multWorkspaceInitial,
                                Real *multWorkspaceIntermediate,
                                Real *outputData,
                                IndexRange rowRange = IndexRange ( -1, -1 ),
                                IndexRange columnRange = IndexRange( -1, -1 ),
                                bool transpose = false, bool left = false );

    // Same as the function above, but assumes that we are multiplying
    // with all interactions in the given block below ancestor_idx.
    //
    // Note: we don't need any row/column range stuff here, because
    // this will not be used for decomposing diagonal blocks
    void interiorBlockSchurMultiply(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                IntArray &nodeList,
                                int interior_block_idx,
                                int ancestor_idx,
                                const Real *G,
                                int nRowsG, int nColsG,
                                const IntArray &inverseRowList, int nRows,
                                Real *multWorkspaceInitial,
                                Real *multWorkspaceIntermediate,
                                Real *outputData,
                                int sparseRowStartOverride,
                                bool transpose = false );

    // Helper function for the above
    //
    //
    void transformInteriorBlockWorkspace( int interior_block_idx,
                                          int interaction_idx1,
                                          int interaction_idx2,
                                          int nCols,
                                          const Real *inputWorkspace,
                                          Real *outputWorkspace );

    // Different options for which sparse rows to use when doing an interior
    // block interaction multiplication
    enum InteriorBlockMultType {
      MULT_FULL = 0,
      MULT_ROW_SET,
      MULT_ROW_RANGE
    };

    // Performs multiplication with a single row block from an interior block
    // in the factor
    void interiorBlockInteractionMultiply(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const IntArray &nodeList,
                                int interior_block_idx, int interaction_idx,
                                const Real *G,
                                int nRowsG, int nColsG,
                                InteriorBlockMultType multType,
                                // In case we only want to use a certain
                                // subset of rows when doing sparse multiplication
                                // (ie. if multType == MULT_ROW_SET)
                                const IntArray &inverseRowList, int nRows,
                                Real *multWorkspaceInitial,
                                Real *solveWorkspace,
                                Real *outputData,
                                // Factor to multiply by when adding to
                                // the output.  Note; alpha is ignored for
                                // transposed multiplications
                                Real alpha,
                                // Optional row range in the sparse matrix
                                // (used if multType == MULT_ROW_RANGE)
                                IndexRange rowRange,
                                bool transpose,
                                int sparseRowStartOverride = -1 );

#if 0
    // Similar to the above function, but multiplies with all interactions
    // *BELOW* ancestor_idx.
    //
    // Note: We don't need the rowRange parameter here because this
    // will not be used for diagonal block decomposition
    void interiorBlockInteractionMultiply_allRows(
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const IntArray &nodeList,
                                int interior_block_idx, int ancestor_idx,
                                const Real *G,
                                int nRowsG, int nColsG,
                                InteriorBlockMultType multType,
                                // In case we only want to use a certain
                                // subset of rows when doing sparse multiplication
                                // (ie. if multType == MULT_ROW_SET)
                                const IntArray &inverseRowList, int nRows,
                                Real *multWorkspaceInitial,
                                Real *solveWorkspace,
                                Real *outputData,
                                // Factor to multiply by when adding to
                                // the output.  Note; alpha is ignored for
                                // transposed multiplications
                                Real alpha,
                                bool transpose = false );
#endif

    // Builds a node list for the function above.
    void buildInteriorBlockNodeList( int interior_block_idx,
                                     int interaction_idx1,
                                     int interaction_idx2,
                                     IntArray &nodeList );

    // Builds a node list for interior block multiplies assuming
    // that we are actually multiplying with everything "below"
    // ancestor_idx
    void buildInteriorBlockNodeList( int interior_block_idx,
                                     int ancestor_idx,
                                     IntArray &nodeList );

    // Multiplies the transpose of the sparse matrix over the given
    // row range, with the set of node indices specifying the columns to use
    void sparseSubMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                            const IntArray &nodeList,
                            int sparseRowStart, int numSparseRows,
                            const Real *G, int nColsG,
                            Real *outputMatrix,
                            Real alpha );

    // Only multiplies using a subset of the matrix rows using the provided
    // inverse row map
    void sparseSubMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                            const IntArray &nodeList,
                            int sparseRowStart, int numSparseRows,
                            const IntArray &inverseRowList,
                            const Real *G, int nColsG,
                            Real *outputMatrix,
                            Real alpha );

#if 0
    // Same as the above, but G is applied on the left
    void sparseSubLeftMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const IntArray &nodeList,
                                int fullColumnCount,
                                int sparseRowStart, int numSparseRows,
                                const Real *G, int nRowsG,
                                Real *outputMatrix );
#endif

    // Multiplies the sparse matrix over the given row range, with the
    // set of node indices specifying the columns to use
    //
    // This is a helper function for interiorBlockSchurMultiply
    void sparseSubTransposeMultiply(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            const IntArray &nodeList,
                            int sparseRowStart, int numSparseRows,
                            const Real *G, int nColsG,
                            Real *outputMatrix,
                            Real alpha );

    // Only multiplies using a subset of the matrix rows using the
    // provided inverse row map
    void sparseSubTransposeMultiply(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            const IntArray &nodeList,
                            int sparseRowStart, int numSparseRows,
                            const IntArray &inverseRowList,
                            const Real *G, int nColsG,
                            Real *outputMatrix,
                            Real alpha );


#if 0
    // Same as the above, but G is applied on the left
    void sparseSubTransposeLeftMultiply(
                            const SPARSE_MATRIX::SparseColumnMatrix &A,
                            const IntArray &nodeList,
                            int fullColumnCount,
                            int sparseRowStart, int numSparseRows,
                            const Real *G, int nRowsG,
                            Real *outputMatrix );
#endif

    // Given the list of current descendents for a node, constructs
    // a list of interior blocks to which those descendents belong
    void findInteriorBlockDescendents();
    void clearInteriorBlockDescendents();

    inline bool isInInteriorBlock( int node_idx ) {
      return _useInteriorBlocks && _interiorBlockMap[ node_idx ] != EMPTY;
    }

  private:
    static const int                               EMPTY;
    static const int                               OVER_SAMPLING;
    static const int                               POWER_ITERATIONS;
    //static const int                               DIAGONAL_POWER_ITERATIONS;

    const SPARSE_MATRIX::SparseColumnMatrix       &_A;

#if 0
    const std::vector<const NestedDissection::DissectionNode *>
                                                  &_ordering;
#endif
    const std::vector<Ordering::SupernodeSpecification>
                                                  &_ordering;

    const std::vector<std::vector<DenseBlock> >   *_diagonalBlocks;

    std::vector<Supernode>                         _factor;

    // Actual numerical factor data
    Real                                          *_factorData;

    long int                                       _factorDataSz;

    //////////////////////////////////////////////////////////////////////
    // Workspaces for various parts of the factorization
    //////////////////////////////////////////////////////////////////////
    Real                                          *_multWorkspace;
    Real                                          *_updateWorkspace;

    // Workspace for system solving
    Real                                          *_solveWorkspace;

    int                                            _multWorkspaceSz;
    int                                            _updateWorkspaceSz;
    int                                            _solveWorkspaceSz;

    WorkspaceManager<Real>                         _realWorkspaceManager;

    //////////////////////////////////////////////////////////////////////
    // Maps needed in the factorization
    //////////////////////////////////////////////////////////////////////
    IntArray                                       _scatteredMap;
    IntArray                                       _lowRankScatteredMap;

    IntArray                                       _relativeMap;
    IntArray                                       _lowRankRelativeMap;
    
    std::vector<IndexPair>                         _interactionIndices;
    std::vector<IndexPair>                         _lowRankInteractionIndices;
    std::vector<IndexPair>                         _extendedInteractionIndices;

    //////////////////////////////////////////////////////////////////////
    // Variables for tracking descendent relationships during
    // numerical factorization
    //////////////////////////////////////////////////////////////////////

    // Equivalent to Lpos in CHOLMOD (see t_cholmod_super_numeric)
    IntArray                                       _nextInteraction;
    IntArray                                       _nextInteractionCache;

    // Equivalent to the corresponding fields in CHOLMOD
    IntArray                                       _head;
    IntArray                                       _next;

    //////////////////////////////////////////////////////////////////////
    // Data and workspaces needed for forming low rank factorizations
    //////////////////////////////////////////////////////////////////////

    int                                            _numExtendedNodes;

    // Actual extra data
    Real                                          *_slackData;

    // Workspace for system solves with extended data
    Real                                          *_extendedSolveWorkspace;
    int                                            _extendedSolveWorkspaceSz;

    long int                                       _slackDataSz;
    long int                                       _availableSlackData;
    long int                                       _slackDataOffset;

    // List of interactions for a given node which are to be
    // replaced with low rank blocks
    IntArray                                       _lowRankBlocks;
    IntArray                                       _lowRankBaseRows;
    IntArray                                       _blockRanks;
    IntArray                                       _blockRanksCache;
    BoolArray                                      _lowRankBlockActive;

    // For a given interaction to be sparsified, list of
    // descendents with contributions to this interaction, and
    // the index of the relevant interaction in the descendent.
    std::vector<std::vector<IndexPair> >           _lowRankDescendents;

    // Workspace which stores the starting row in the original
    // matrix for a set of low rank interactions
    std::vector<IntArray>                          _lowRankStartRows;

    // Gaussian test matrix used for forming low-rank decompositions
    // We also need to be able to extract row sub-matrices from
    // the test matrix.
#if 0
    MATRIX                                         _testMatrix;
    Real                                          *_subMatrix;
#endif

    int                                            _overSampling;

    // TODO: This could probably be the same as the existing
    // _multWorkspace, but for now we'll keep it separate.
    Real                                          *_decompWorkspace;
    int                                            _decompWorkspaceSz;

#if 0
    Real                                          *_decompMultWorkspace;
    int                                            _decompMultWorkspaceSz;
#endif

    Real                                          *_decompMultTransWorkspace;
    int                                            _decompMultTransWorkspaceSz;

#if 0
    // We need a workspace for copied data in during low-rank
    // decomposition
    Real                                          *_copyWorkspace;
    long int                                       _copyWorkspaceSz;
#endif

    // Additional data array needed to store QR factorization information
    // (tau in the parameters to geqrf)
    Real                                          *_qrExtraData;
    int                                            _qrExtraDataSz;

    // Workspace for QR factorization
    Real                                          *_qrWorkspace;
    int                                            _qrWorkspaceSz;

    //////////////////////////////////////////////////////////////////////
    // For tolerance-adapted low-rank decomposition
    //////////////////////////////////////////////////////////////////////

    static const int                               MAX_BLOCK_SIZE = 64;
    static const int                               ERROR_POWER_ITERATIONS = 6;

    public:
    // DSB:
#ifdef HAS_CXX11_CONSTEXPR
    static constexpr Real                          ERROR_BOUND_MULTIPLIER = 1.0;
#else
    static const Real                              ERROR_BOUND_MULTIPLIER = 1.0;
#endif

    private:

    Real                                           _errorBoundMultiplier;

    // Another multiplication workspace
#if 0
    Real                                          *_blockWorkspace;
    int                                            _blockWorkspaceSz;
#endif

    // Vector workspaces used for error estimation (more generally, for
    // matrix vector multiplication)
    Real                                          *_vectorWork1;
    Real                                          *_vectorWork2;
    int                                            _vectorWorkSz;
    int                                            _vectorWork1Sz;
    int                                            _vectorWork2Sz;

    // Workspace for forming block basis columns
    Real                                          *_basisWorkspace;
    int                                            _basisWorkspaceSz;

#if 0
    // For expanding compressed interactions
    Real                                          *_expansionWorkspace;
    long int                                       _expansionWorkspaceSz;
#endif

    //////////////////////////////////////////////////////////////////////
    // Used to manage storage of low rank decompositions for a
    // given node as they are being built (since we don't put
    // them in to the actual factorization until we know how big
    // they are going to be)
    //////////////////////////////////////////////////////////////////////

    struct BlockStorageData {
      BlockStorageData()
        : _data( NULL ),
          _nRows( 0 ),
          _nCols( 0 ),
          _maxSize( 0 ),
          _full( false )
      {
      }

      Real                                        *_data;

      int                                          _nRows;
      int                                          _nCols;
      int                                          _maxSize;

      bool                                         _full;
    };

    struct BlockStorageManager {
      BlockStorageManager()
        : _data( NULL ),
          _dataSz( 0 )
      {
      }

      ~BlockStorageManager()
      {
        delete[] _data;
      }

      Real                                        *_data;
      int                                          _dataSz;

      std::vector<BlockStorageData>                _dataPtrs;
    };

    // Temporary storage for orthogonal bases Q while they are
    // being built
    BlockStorageManager                            _basisStorage;

    // Current error measures for each block
    std::vector<VEC3F>                             _blockErrors;

    //////////////////////////////////////////////////////////////////////
    // Stuff for sparsification of diagonal blocks
    //////////////////////////////////////////////////////////////////////

    // List of all descendent indices for the current node being
    // factored
    IntArray                                       _currentDescendents;

    // For each of these descendents, we will construct a list of
    // row/column ranges to consider when forming low-rank decompositions
    // for each block in the current node's main diagonal
    std::vector<std::vector<DenseBlock> >          _diagonalBlockRanges;

    // Maps all row indices from the current supernode to their
    // corresponding index in to the row list of a descendent's interaction
    IntArray                                       _inverseRowMap;

    // Which diagonal low rank blocks are currently active
    BoolArray                                      _diagonalLowRankBlockActive;

    //////////////////////////////////////////////////////////////////////
    // Variables for uncompressed supernode's implicitly (ie. without
    // storing their interactions with compressed supernodes)
    //////////////////////////////////////////////////////////////////////

    bool                                           _useInteriorBlocks;

    // The blocks themselves; identified by supernode index ranges
    std::vector<InteriorBlock>                     _interiorBlocks;

    // Inverse map
    IntArray                                       _interiorBlockMap;

    // Current interior block descendents
    IntArray                                       _currentBlockDescendents;
    BoolArray                                      _blockDescendentFound;

    //////////////////////////////////////////////////////////////////////
    // Variables here are used to support storing low-rank decompositions
    // of entire off-diagonal blocks for each node, rather than storing
    // per-interaction decompositions
    //////////////////////////////////////////////////////////////////////

    OffDiagonalCompressionType                     _offDiagCompressionType;

    //////////////////////////////////////////////////////////////////////
    // Function pointers for estimating rank and constants used for
    // this purpose
    //////////////////////////////////////////////////////////////////////
    Supernode::RankEstimator                       _offDiagonalRankEstimator;
    Supernode::RankEstimator                       _diagonalRankEstimator;

    Real                                           _rankConstant;
    Real                                           _diagonalRankConstant;

    //////////////////////////////////////////////////////////////////////
    // Timing data
    //////////////////////////////////////////////////////////////////////

    enum TimerID {
      BASE_MATRIX_COPY = 0,
      BLOCK_INIT,
      INTERACTION_INTERSECTION,
      RELATIVE_MAP,
      STANDARD_UPDATE,
      EXTENDED_UPDATE_STANDARD,
      EXTENDED_UPDATE_EXTENDED,
      LOW_RANK_DESC_LIST,
      ORTHOGONAL_PROJECTION,
      QR_FACTORIZATION,
      BASIS_COPY,
      RANDOM_MULTIPLY_GENERATE,
      RANDOM_MULTIPLY_INIT,
      RANDOM_MULTIPLY_APPLY,
      UPDATE_ACTIVE_BLOCKS,
      ERROR_ESTIMATE,
      BASIS_PROJECT,
      ASSIGN_BASIS,
      INTERACTION_EXPANSION,
      FIND_BLOCK_RANGES,
      DIAGONAL_ORTHOGONAL_PROJECTION,
      DIAGONAL_QR_FACTORIZATION,
      DIAGONAL_BASIS_COPY,
      DIAGONAL_RANDOM_MULTIPLY,
      DIAGONAL_UPDATE_ACTIVE_BLOCK,
      DIAGONAL_ERROR_ESTIMATE,
      DIAGONAL_BASIS_PROJECT,
      DIAGONAL_ASSIGN_BASIS,
      ADD_DIAGONAL_CONTRIBUTIONS,
      FACTOR_NODE,
      INIT_EXTENDED,
      EXTENDED_SCHUR_FORMATION,
      EXTENDED_SCHUR_INVERSION,
      EXTENDED_SCHUR_MULTIPLY,
      INTERIOR_BLOCK_PREAMBLE,
      INTERIOR_BLOCK_NODELIST,
      INTERIOR_BLOCK_WORKSPACE,
      INTERIOR_BLOCK_MULT_SETUP,
      INTERIOR_BLOCK_MULT_SPARSE,
      INTERIOR_BLOCK_MULT_WORKSPACE,
      INTERIOR_BLOCK_MULT_FORWARD_SOLVE,
      INTERIOR_BLOCK_MULT_BACKWARD_SOLVE,
      INTERIOR_BLOCK_MULT_TRISOLVE,
      INTERIOR_BLOCK_MULT_PROPAGATE,
      INTERIOR_BLOCK_COPY,
      NUM_TIMERS
    };

    Timer                                          _totalTime;
    std::vector<Timer>                             _timers;

    static const char                             *TIMER_NAMES[];

    // FIXME: debugging
    friend int main( int argc, char **argv );

};

#endif
