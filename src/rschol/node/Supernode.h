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
// Supernode.h: Interface for the Supernode class
//
//////////////////////////////////////////////////////////////////////

#ifndef SUPERNODE_H
#define SUPERNODE_H

#include <rschol/datastructure/WorkspaceManager.h>

#include <rschol/linearalgebra/DenseBlock.h>
#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/SPARSE_MATRIX.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <rschol/library_interface/CHOLMOD_Environment.h>

#include <rschol/node/DiagonalBlock.h>

#include <rschol/ordering/NestedDissection.h>
#include <rschol/ordering/Ordering.h>

#include <rschol/util/timer.h>

#include <boost/function.hpp>

#include <iostream>
#include <vector>

//////////////////////////////////////////////////////////////////////
// SupernodeType enum
//
// We need this to distinguish between standard supernodes appearing
// in the unmodified factorizations, and additional supernodes
// introduced by us
//
// We also allow for "compressed nodes"; that is, nodes in which
// low rank decompositions are stored in place so that no extended
// variables need to be introduced.
//////////////////////////////////////////////////////////////////////
enum SupernodeType {
  STANDARD_NODE = 0,
  EXTENDED_NODE,
  COMPRESSED_NODE,
  // For representing implicit interactions
  IMPLICIT_NODE
};

//////////////////////////////////////////////////////////////////////
// SupernodeInteraction struct
//
// Helper class summarizing interactions between supernodes in a
// sparse Cholesky factorization
//////////////////////////////////////////////////////////////////////
struct SupernodeInteraction {
  SupernodeInteraction()
    : _nodeID( -1 ),
      _type( STANDARD_NODE ),
      _active( true ),
      _numExtendedRows( 0 ),
      _extendedDataOffset( -1 ),
      _extendedOffset( -1 ),
      _extendedInteraction( -1 ),
      _firstInteraction( false ),
      _compressed( false ),
      //_compressedColumnList( NULL ),
      _compressedColumnRange( -2, -1 ),
      _lowRankDiagonalInteraction( false ),
      _initialLowRankDiagonalIndex( -1 ),
      _diagonalContributor( false ),
      _nonZeroColumns( 0 ),
      _diagonalContributionScale( 0.0 )
  {
  }

  SupernodeInteraction( int nodeID, SupernodeType type,
                        bool firstInteraction = false,
                        int extendedInteraction = -1,
                        bool compressed = false,
                        IntArray *compressedColumnList = NULL,
                        IndexRange *compressedColumnRange = NULL,
                        bool diagonalContributor = false,
                        bool lowRankDiagonalInteraction = false,
                        int initialLowRankDiagonalIndex = -1 )
    : _nodeID( nodeID ),
      _type( type ),
      _active( true ),
      _numExtendedRows( 0 ),
      _extendedDataOffset( -1 ),
      _extendedOffset( -1 ),
      _extendedInteraction( extendedInteraction ),
      _firstInteraction( firstInteraction ),
      _compressed( compressed ),
      _diagonalContributor( diagonalContributor ),
      _lowRankDiagonalInteraction( lowRankDiagonalInteraction ),
      _initialLowRankDiagonalIndex( initialLowRankDiagonalIndex ),
      _nonZeroColumns( 0 ),
      _diagonalContributionScale( 0.0 ),
      _compressedV( NULL ),
      _compressedU( NULL )
  {
    if ( compressedColumnList )
    {
      _compressedColumnList = *compressedColumnList;
    }

    if ( compressedColumnRange )
    {
      _compressedColumnRange = *compressedColumnRange;
    }
  }

  inline bool hasParent() const
  {
    return ( _type == EXTENDED_NODE && _extendedInteraction >= 0 );
  }

  static const int     NO_INTERACTION_PARENT = -10;

  int                  _nodeID;
  SupernodeType        _type;

  // The rows over which the owner of this object
  // interacts with node _nodeID.
  // This list is relative to supernode _nodeID.
  IntArray             _rowList;

  // Integer offsets in to a node's data array
  IntArray             _dataRowList;

  // Used to indicate when a particular interaction has been
  // replaced with slack variables.  Interactions which have
  // been replaced are flagged as inactive.
  bool                 _active;

  // For extended interactions, we only need to store
  // the number of rows associated with the interaction
  // since it will be full.  We also store the offset
  // for this row data in to some array.
  int                  _numExtendedRows;
  long int             _extendedDataOffset;

  // If this interaction is inactive, we also store an offset
  // from the current node to the ID of the extended node
  // introduced to compensate for this.
  int                  _extendedOffset;

  // If this interaction is inactive, store the index of the
  // interaction introduced to compensate for it.
  // 
  // If this interaction is active, and refers to an extended
  // node, and _firstInteraction == true, then store the index
  // of the original interaction that produced this one.
  int                  _extendedInteraction;

  // Stores a list of nodes occuring after _nodeID as well
  // as offsets in to their interaction lists.  This tracks
  // any additional interactions imposed by this one, and
  // allows us to quickly allocate data arrays for these
  // interactions.
  PairArray            _forwardInteractions;

  // Whether or not this node is the first node to
  // interact with another node (only used for extended interactions)
  bool                 _firstInteraction;

  // Some extended interactions will be stored in "scattered" column
  // form.  That is, is an interaction with 4 columns has only 2
  // columns filled in, we will only store those 2 columns.
  // In this case, we will use _compressedColumnList to track precisely which
  // columns are filled in.
  bool                 _compressed;
  IntArray             _compressedColumnList;

  // Also store a column range.  This is only used for interactions resulting
  // from low rank blocks on the diagonal (a bit of a hack).
  IndexRange           _compressedColumnRange;

  // Whether or not this interaction results from a low rank block
  // in the diagonal
  bool                 _lowRankDiagonalInteraction;

  // If this resulted from a low rank block on the diagonal, store the
  // initial index in to this node's diagonal block arrays assigned
  // to this interaction (in order to accomodate symbolic reordering
  // later on)
  int                  _initialLowRankDiagonalIndex;

  // Whether or not the outer product of this interaction needs to
  // be added to its parent node's diagonal
  // 
  // FIXME: Currently blocks resulting from the diagonal are not
  //        listed here.  Should they be?
  bool                 _diagonalContributor;

  // Number of columns from the original system containing non-zero entries
  int                  _nonZeroColumns;

  // Scaling factor when treating this interaction as a diagonal contribution
  Real                 _diagonalContributionScale;

  //////////////////////////////////////////////////////////////////////
  // These store the low rank decomposition of an interaction stored
  // in compressed form (ie. compressed without introducing a slack
  // variable)
  //////////////////////////////////////////////////////////////////////
  Real                *_compressedV;
  Real                *_compressedU;

};

//////////////////////////////////////////////////////////////////////
// InteriorBlock structure
//
// Models a range of supernodes in a factorization, and the
// interaction of this range with other blocks
//////////////////////////////////////////////////////////////////////
struct InteriorBlock {
  typedef std::vector<IndexPair> Interaction;

  InteriorBlock()
    : _nodeRange( -1, -1 )
  {
  }

  InteriorBlock( int start_idx, int end_idx )
    : _nodeRange( start_idx, end_idx )
  {
  }

  inline int numNodes() const {
    return range_size( _nodeRange );
  }

  inline int localNodeIdx( int node_idx ) const {
    return node_idx - _nodeRange.first;
  }

  IndexRange                 _nodeRange;

  // Interactions, and a map from node indices to their place
  // in the interaction list
  std::vector<Interaction>   _interactions;
  IntArray                   _interactionMap;

  // Row ranges for sparse matrix entries in this block's
  // off-diagonal (useful for performing multiplication later
  // on)
  IntArray                   _firstOffDiagonalRows;

};

static bool operator<( const SupernodeInteraction &a,
                       const SupernodeInteraction &b )
{
  return a._nodeID < b._nodeID;
}

static bool operator==( const SupernodeInteraction &a,
                        const SupernodeInteraction &b )
{
  return a._nodeID == b._nodeID;
}

// Structure for representing a Bunch-Kaufman factorization, rather
// than a pure Cholesky factorization
struct LDL_Data {
  // Pivot data returned by a LDL^{T} factorization (presumably
  // Bunch-Kaufman)
  IntArray               _pivotData;

  typedef Real                               Block_1x1;
  typedef Real                               Inverse_1x1;

  // 2x2 block definition - stores an array of integers to allow
  // for pivot information from an LU factorization
  typedef struct Block_2x2 {
    Real data[ 4 ];
  } Block_2x2;

  typedef struct Inverse_2x2 {
    Real data[ 4 ];
    int  pivotData[ 2 ];
  } Inverse_2x2;

  // Diagonal block data, and their inverses
  std::vector<Block_1x1>     _diagonalEntries;
  std::vector<Inverse_1x1>   _diagonalInverses;

  std::vector<Block_2x2>     _diagonalBlockEntries;
  std::vector<Inverse_2x2>   _diagonalBlockInverses;
};

//////////////////////////////////////////////////////////////////////
// Supernode class
//
// Models a supernode appearing in a sparse supernodal Cholesky
// factorization
//////////////////////////////////////////////////////////////////////
class Supernode {
	public:
    // FIXME
    static int numLargeBlocks;
    static long int multWorkspaceSz;
    static long int workspaceSz;

		Supernode( SupernodeType type = STANDARD_NODE );

		// Destructor
		virtual ~Supernode();

    // Returns the fixed floating point memory requirement
    // for this node.  That is, floating point requirements
    // based on the node's diagonal entry, and interactions
    // with other standard nodes.
    int numFixedEntries() const;

    // Flags an interaction as implicit
    void flagImplicitInteraction( int interaction_idx );

    // Builds a map from matrix indices to row indices in this
    // node's fixed data array.  This is necessary so that we
    // can copy entries from the original matrix in to the node.
    //
    // The input array is assumed to have the correct size.
    //
    // Also, assuming that we stack all low rank blocks together,
    // eg, U = [ U1; U2; ... ]
    // build s scattered map of matrix indices in to rows of U
    void constructScatteredMap( IntArray &scatteredMap,
                                IntArray &lowRankScatteredMap,
                                const std::vector<Supernode> &nodes );

    // Clears any changes made by this node to the scattered map
    void clearScatteredMap( IntArray &scatteredMap,
                            IntArray &lowRankScatteredMap,
                            const std::vector<Supernode> &nodes );

    // Copies entries from the original matrix to the appropriate
    // location in this node's data array (does nothing for an
    // extended node).
    //
    // Assumes the scattered map has been constructed for this node
    void copyMatrixData( const IntArray &scatteredMap,
                         const IntArray &lowRankScatteredMap,
                         const SPARSE_MATRIX::SparseColumnMatrix &A,
                         Real *extraData );

    // Copies a column subset of matrix data for the given interaction
    // between this node and ancestor.
    // 
    // Data is scattered according to the row list in this interaction
    void copyMatrixDataBlockColumn( int interaction_idx,
                                    const Supernode &ancestor,
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const IndexRange &columnRange,
                                    Real *workspace );

    // Copies matrix data for the node in _nodeDiagonal rooted at
    // block_idx in to the provided matrix
    void copyDiagonalBlockMatrixData(
                          int block_idx,
                          const SPARSE_MATRIX::SparseColumnMatrix &A,
                          Real *workspace ) const;

    // Finds a list of mutual off-diagonal blocks shared by this
    // node and a given descendent.
    //
    // We must be provided with the index in to d._offDiagonal
    // corresponding to "this" node.
    //
    // lowRankIndices will return interaction index pairs for
    // interactions from descendent that should be accumulated into
    // a low rank block of *this* node.
    //
    // extendedIndices returns interaction index pairs for
    // interactions from a descendent that are accumulated in to
    // an "extended" interaction.
    void findInteractionIntersection( const Supernode &descendent,
                                      int start_idx,
                                      PairArray &indices,
                                      PairArray &lowRankIndices,
                                      PairArray &extendedIndices );

    // Suppose we put the set of rows from descendent in to a matrix
    // [ L1; L2 ] where L1 corresponds to rows lying in this node's
    // column range.  Here, we build a map from rows (columns) in this
    // matrix to rows in the _data member of this node.  This is done
    // so that when we form the product
    //  / L1 \ * L1'
    //  \ L2 /
    // we know where to add components to this node's _data.
    //
    // indices is computed by findInteractionIntersection and stores
    // the intersection of the descendent._offDiagonal and this->_offDiagonal.
    void buildRelativeMap( const Supernode &descendent,
                           int start_idx,
                           const PairArray &indices,
                           IntArray &relativeMap );

    // Suppose we write the set of low rank interaction blocks
    // in this node as
    //      / U1 \
    //      | U2 |
    //  U = | .  |
    //      | .  |
    //      \ . /
    // and we wish to form the product U * G for some matrix G.
    // This function builds a relative map from data rows in descendent
    // to rows in U.
    void buildLowRankRelativeMap( const Supernode &descendent,
                                  const PairArray &lowRankIndices,
                                  const IntArray &baseRows,
                                  IntArray &relativeMap );

    // Copies data relative to the interaction between descendent and
    // soem later node in to the provided workspace so that we can do a
    // single dense multiplication
    //
    // Builds a matrix of the form [ L1; L2 ] and returns the
    // number of rows/columns associated with each piece
    void copyInteractionData( const Supernode &descendent,
                              int start_idx,
                              const PairArray &indices,
                              Real *workspace,
                              int &workCols,
                              int &workRows1, int &workRows2 ) const;

    // Given a matrix of the form [ L1; L2 ], computes the update matrix
    // [ L1; L2 ] * L1' and puts it in update.  This update matrix will
    // later be added to this node's interactions.
    //
    // This is not static because we may need to look at the structure
    // of this node's diagonal do decide how to build the matrix.
    void buildUpdateMatrix( const Real *workspace,
                            int workCols, int workRows1, int workRows2,
                            const IntArray &relativeMap,
                            Real *update,
                            int descendent_idx ) const;

    // Let L1 be the interaction in this indexed by ancestor_idx.
    // Let L2 be the interaction in this indexed by ancestor_interaction_idx.
    // This function forms a block column of the product L2 * L1' and
    // subtracts this from the given workspace matrix.
    //
    // Rows of L2 * L1' are scattered in to workspace according to the
    // provided workspace row list.
    //
    // This function is used in the process of accumulating diagonal
    // contributions as a result of introducing slack variables
    void subtractInteractionUpdateMatrix( int ancestor_idx,
                                          int ancestor_interaction_idx,
                                          const IndexRange &columnRange,
                                          const IntArray &workspaceRows,
                                          Real *multWorkspace,
                                          Real *workspace ) const;

    // Subtracts the content of the given update matrix from this
    // node's data according to the given relative map
    //
    // workRows1, workRows2 are the same variables used by
    // copyInteractionData and buildUpdateMatrix
    void applyUpdate( const IntArray &relativeMap,
                      int workRows1, int workRows2,
                      const Real *update );

    // Subtracts the outer product of the given matrix from this
    // node's diagonal.
    //
    // update's number of rows should match this node's column count
    void applyDiagonalUpdate( const Real *update, int nCols );

    // Similar to the above, but uses a row list rather than a full
    // matrix
    void applyDiagonalUpdate( const Real *update, 
                              Real *multWorkspace, int nCols,
                              const IntArray &rowList );

    // This version forms a diagonal update for a provided diagonal block
    // in this node's main diagonal.
    //
    // Forms a diagonal update and puts it in to the provided schur
    // complement matrix.
    void applyDiagonalUpdate( const DenseBlock &factorBlock,
                              const Real *update,
                              Real *multWorkspace, int nCols,
                              const IntArray &rowList,
                              Real *schurComplement,
                              // Optional range for rowList
                              IndexRange rowRange = IndexRange(-1,-1) ) const;

    // Update this node's diagonal with the contribution from a
    // descenent, but only if the relevent interaction in that descendent
    // is stored in compressed form
    void compressedDiagonalUpdate( const Supernode &descendent,
                                   int ancestor_idx );

    // A compressed diagonal update which explicitly subtracts its result
    // in to the given schur complement matrix.
    //
    // This is similar to applyDiagonalUpdate above, but update has
    // been replaced by V * U' = update
    //
    // multWorkspace should be of size at least:
    //    rank * ( rank + nCols + range_size( rowRange ) )
    //  if rowRange is provided or:
    //    rank * ( rank + nCols + rowList.size() )
    //  if not
    void compressedDiagonalUpdate(
                            const DenseBlock &factorBlock,
                            const Real *V, const Real *U, int rank,
                            Real *multWorkspace, int nCols,
                            const IntArray &rowList,
                            Real *schurComplement,
                            // Optional range for rowList
                            IndexRange rowRange = IndexRange(-1,-1) ) const;

    // Forms update matrices for each extended node interaction between
    // the given descendent node and this.  Subtract these from the
    // relevant extended interactions in this node.
    //
    // We require a workspace for multiplication, and the extra data
    // array in which extended interactions are stored.
    //
    // start_idx is the same input used to buildRelativeMap, etc.
    //
    // Since some of the extended interactions for an extended node
    // may be stored in scattered form, we use expansionWorkspace
    // to form the full matrix representation for these interactions
    // (zeros included).  This workspace must be large enough to
    // accomodate the descendent._offDiagonal[ start_idx ] and any
    // subsequent interaction in descendent.
    void applyExtendedUpdate( const Supernode &descendent,
                              const PairArray &extendedIndices,
                              Real *multWorkspace,
                              Real *extraData,
                              int start_idx,
                              Real *expansionWorkspace = NULL,
                              Timer *standardTimer = NULL,
                              Timer *extendedTimer = NULL );

    void applyExtendedUpdate( const Supernode &descendent,
                              const PairArray &extendedIndices,
                              Real *extraData,
                              int start_idx,
                              WorkspaceManager<Real> &workspaceManager,
                              Timer *standardTimer = NULL,
                              Timer *extendedTimer = NULL );

    // Returns a list of interactions to be sparsified in this node.
    //
    // Also, assuming we stack the contents of these blocks on top
    // of each other, store the starting row in this structure for
    // each of these blocks.
    void getLowRankBlocks( IntArray &lowRankBlocks,
                           IntArray &baseRows ) const;

    // For a given descendent, given the low rank interaction intersection
    // between these two nodes, and a list of all of the low rank blocks
    // for this node, identify which low rank blocks the descendent
    // interacts with.  Append the descendent's ID to a list of IDs for
    // each low rank block.
    //
    // Inputs:
    //    descendent: The descendent supernode
    //    lowRankIndices: As built by findInteractionIntersection
    //    lowRankBlocks: As built by getLowRankBlocks
    //    lowRankDescendents: List of IndexPair lists to append to.
    //                        Should be at least as large as lowRankBlocks
    void appendLowRankDescendentList(
                    const Supernode &descendent,
                    const PairArray &lowRankIndices,
                    const IntArray &lowRankBlocks,
                    std::vector<PairArray > &lowRankDescendents );

    void clearLowRankDescendentList(
                    int numEntries,
                    std::vector<PairArray> &lowRankDescendents );

#if 0
    // Given a set of interactions (U1, U2, ...) from *this* node,
    // and the interaction U_a with an ancestor node, form the
    // matrix product [ U1 U2 ... ]' * U_a' * G (where G is a given
    // matrix) and accumulate it in to the given workspace.
    void interactionMultiply( const PairArray &lowRankIndices,
                              int ancestor_idx,
                              const Supernode &ancestor,
                              const Real *G, int nColsG,
                              Real *multWorkspaceInitial,
                              Real *multWorkspaceFinal );
#endif

    // Given the interactions below ancestor_idx in this node (call these
    // (U1, U2, ...) and the interaction U_a with the ancestor node,
    // form the matrix product [ U1; U2; ... ] * U_a' * G (where G is a
    // given matrix) and place the result in the given matrix (doesn't
    // do any row selection or scattering so the input and output will
    // have to be processed)
    void interactionMultiply( int ancestor_idx,
                              const Supernode &ancestor,
                              const Real *G, int nColsG,
                              Real *multWorkspaceInitial,
                              Real *multWorkspaceFinal );

    // Row counting-routine used by the above function
    int countInteractionRows( int startInteraction, int endInteraction );

    // Similar to the above, but only handles a single interaction
    //
    // Overwrites multWorkspaceFinal
    void interactionMultiply( int interaction_idx,
                              int ancestor_idx,
                              const Supernode &ancestor,
                              const Real *G, int nColsG,
                              Real *multWorkspaceInitial,
                              Real *multWorkspaceFinal );

    // Assuming that U_a is the intersaction between this node and
    // the given ancestor node, forms a subset of the matrix U_a * U_a'
    // with the specified row and column range.  These ranges given
    // in terms of the actual set of rows stored in U_a.
    // That is, if rowRange = [ 2, 3 ] and this node stores rows
    // r0, r1, r2, r3, r4, ..., then we are considering the row subset
    // [ r2, r3 ] of U_a
    void interactionMultiply( int ancestor_idx,
                              const Supernode &ancestor,
                              const DenseBlock &ancestorBlock,
                              const Real *G, int nColsG,
                              const IndexRange &rowRange,
                              const IndexRange &columnRange,
                              Real *multWorkspaceInitial,
                              Real *multWorkspaceFinal );

#if 0
    // Given a particular interaction U1 from *this* node, and
    // the interaction U_a with an ancestor node, form the matrix
    // product Q' * U1 * U_a' (where Q is a given matrix) and
    // accumulate it in to the given workspace.
    //
    // We also require a workspace (subMatrix) to extract the necessary
    // rows from Q before multiplying (since our stored version of
    // U1 may not contain all rows)
    //
    // Q may optionally be stored directly as its transpose, in
    // which case transposeQ should be set to false
    void interactionTransMultiply( int interaction_idx,
                                   int ancestor_idx,
                                   const IntArray &rowsQ,
                                   const Supernode &ancestor,
                                   const Real *Q, int nColsQ,
                                   Real *subMatrix,
                                   Real *multWorkspaceInitial,
                                   Real *multWorkspaceFinal );

    // Similar to the above, but only multiplies Q' on the right by a
    // subset of the matrix U_a * U_a' (only used for forming pieces of
    // the diagonal block of ancestor)
    //
    // The ranges rowRange and columnRange
    // in terms of the actual set of rows stored in U_a.
    // That is, if rowRange = [ 2, 3 ] and this node stores rows
    // r0, r1, r2, r3, r4, ..., then we are considering the row subset
    // [ r2, r3 ] of U_a
    void interactionTransMultiply( int ancestor_idx,
                                   const Supernode &ancestor,
                                   const DenseBlock &ancestorBlock,
                                   const Real *Q, int nColsQ,
                                   const IndexRange &rowRange,
                                   const IndexRange &columnRange,
                                   Real *subMatrix,
                                   Real *multWorkspaceInitial,
                                   Real *multWorkspaceFinal );
#endif

    // Given the interactions below ancestor_idx in this node (call these
    // (U1, U2, ...) and the interaction U_a with the ancestor node,
    // form the matrix product U_a * [ U1; U2; ... ]' * G (where G is a
    // given matrix) and place the result in the given matrix (doesn't
    // do any row selection or scattering so the input and output will
    // have to be processed)
    void interactionTransMultiply_left( int ancestor_idx,
                                        const Supernode &ancestor,
                                        const Real *G, int nColsG,
                                        Real *multWorkspaceInitial,
                                        Real *multWorkspaceFinal );

    // Similar to the above, but computes (U1 * U_a')' * G for some
    // matrix G.  Requires the row set from the interaction we are
    // trying to build
    void interactionTransMultiply_left( int interaction_idx,
                                        int ancestor_idx,
                                        const IntArray &rowsG,
                                        const Supernode &ancestor,
                                        const Real *G, int nColsG,
                                        Real *subMatrix,
                                        Real *multWorkspaceInitial,
                                        Real *multWorkspaceFinal );

    // Similar to the above, but only multiplies G by a subset of
    // the matrix (U_a * U_a')' (only used for forming piees of the
    // diagonal block of ancestor)
    //
    // The ranges rowRange and columnRange
    // in terms of the actual set of rows stored in U_a.
    // That is, if rowRange = [ 2, 3 ] and this node stores rows
    // r0, r1, r2, r3, r4, ..., then we are considering the row subset
    // [ r2, r3 ] of U_a
    void interactionTransMultiply_left( int ancestor_idx,
                                        const Supernode &ancestor,
                                        const DenseBlock &ancestorBlock,
                                        const Real *G, int nColsG,
                                        const IndexRange &rowRange,
                                        const IndexRange &columnRange,
                                        Real *subMatrix,
                                        Real *multWorkspaceInitial,
                                        Real *multWorkspaceFinal );

    // For in-place decomposition of diagonal blocks:
    // Multiplies "previous" compressed off-diagonal block's from this
    // node's diagonal by the given matrix.
    //
    // This is part of forming the low-rank decomposition of the
    // off-diagonal block indexed by block_idx
    void multiplyPreviousInPlaceBlocks( int block_idx,
                                        const Real *G,
                                        int rowsG, int colsG,
                                        Real *output,
                                        bool left, bool transpose );

    // For in-place decomposition of diagonal blocks:
    // Multiplies (and overwrites) the input matrix with the inverse
    // of the lower-triangular matrix located "above" the off-diagonal
    // block indexed by block_idx
    void inPlaceOffDiagonalSolve( int block_idx,
                                  Real *G, int nRHS,
                                  bool transpose );

#if 0
    // For the set of interactions (U1, U2, ...) from this node that
    // are to be decomposed in to low-rank representations, take the
    // entries from the original matrix that belong in U1, U2, ... etc.
    // and multiply them by the matrix G.
    void baseSystemMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                             const IntArray &lowRankScatteredMap,
                             const IntArray &lowRankBlocks,
                             const Real *G, int nColsG,
                             Real *multWorkspaceFinal );
#endif

    // For a specific interaction in this node, take the entries from
    // the original matrix that belong in its row/column space and
    // multiply them by the matrix G.
    void baseSystemMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                             const std::vector<Supernode> &nodes,
                             int interaction_idx,
                             const Real *G, int nColsG,
                             Real *multWorkspaceFinal );

    // Multiplies by the entries of A corresponding to rows in this node's
    // full off-diagonal (rather than just a single interaction which is
    // handled by the function above)
    void baseSystemMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                             const std::vector<Supernode> &nodes,
                             const Real *G, int nColsG,
                             Real *multWorkspaceFinal );

    // Multiply by the base system of a sub-matrix in this node's main
    // diagonal
    void baseSystemMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                             const Real *G, int nColsG,
                             const IndexRange &rowRange,
                             const IndexRange &columnRange,
                             Real *multWorkspaceFinal );

    // For a single interaction from this node (U1) which is to be
    // decomposed in to a low-rank representation, take the entries
    // of the original matrix that belong in the U1 block and form
    // the product Q' * U1.
    //
    // Q can optionally be stored directly as Q', in which case
    // transposeQ should be set to false
    void baseSystemTransMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  int interaction_idx,
                                  const std::vector<Supernode> &nodes,
                                  const Real *Q, int nColsQ,
                                  Real *multWorkspaceFinal,
                                  bool transposeQ = true );

    // Multiply a matrix on the right by the base system of a sub-matrix
    // in this node's main diagonal.
    //
    // That is, form the product Q' * U1 (as in the above function)
    void baseSystemTransMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  const Real *Q, int nColsQ,
                                  const IndexRange &rowRange,
                                  const IndexRange &columnRange,
                                  Real *multWorkspaceFinal,
                                  bool transposeQ = true );

    // Similar to the above, but multiplies the transpose of the interaction,
    // with some other matrix.  (ie. U1' * Q)
    void baseSystemTransMultiply_left(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  int interaction_idx,
                                  const std::vector<Supernode> &nodes,
                                  const Real *G, int nColsG,
                                  Real *multWorkspaceFinal,
                                  bool clearOutput );

    // Multiplies by the system A in the row space given by this node's
    // full off-diagonal block (as opposed to just a single interaction block)
    void baseSystemTransMultiply_left(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  const std::vector<Supernode> &nodes,
                                  const Real *G, int nColsG,
                                  Real *multWorkspaceFinal );

    // Does the same as the function above, but uses the base system
    // of a sub-matrix in this node's main diagonal
    void baseSystemTransMultiply_left(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                                  const Real *G, int nColsG,
                                  const IndexRange &rowRange,
                                  const IndexRange &columnRange,
                                  Real *multWorkspaceFinal );

    // For a given low rank block (indexed by block_idx) in the diagonal,
    // multiply the matrix G with the sum of all diagonal contributions
    // resulting from slack variable introduction and accumulates the
    // result in multWorkspaceFinal
    //
    // We need to call expandCompressedInteractions before this in order
    // to fill in expansionWorkspace
    void diagonalContributionMult( int block_idx,
                                   const Real *G, int nColsG,
                                   const Real *extraData,
                                   Real *expansionWorkspace,
                                   Real *multWorkspaceInitial,
                                   Real *multWorkspaceFinal,
                                   bool transpose = false,
                                   bool left = true );

    void diagonalContributionMult( int block_idx,
                                   const Real *G, int nColsG,
                                   const Real *extraData,
                                   WorkspaceManager<Real> &workspaceManager,
                                   Real *multWorkspaceFinal,
                                   bool transpose = false,
                                   bool left = true );

    // Given a workspace of computed multiplications of low-rank
    // blocks with a random matrix with nColsG columns, compute
    // a QR factorization for each interaction.
    // We also require some extra data arrays for the QR factorization.
    void lowRankQR( const IntArray &lowRankBlocks,
                    const BoolArray &lowRankBlockActive,
                    Real *decompWorkspace, const IntArray &nBlockCols,
                    Real *qrExtraData, Real *qrWorkspace, int workSize,
                    bool transpose );

    // QR factorization like the function above but assuming that we
    // are decomposing the full off-diagonal of this node all at once
    void lowRankQR( Real *decompWorkspace, int rank,
                    Real *qrExtraData, Real *qrWorkspace, int workSize,
                    bool transpose );

    // Same as the above, but only performst the QR factorization for
    // a single low rank block in the node's diagonal
    void lowRankQR( int block_idx, Real *decompWorkspace, int nColsG,
                    Real *qrExtraData, Real *qrWorkspace, int workSize,
                    bool transpose );

    // Suppose that the given interaction in this node has been
    // decomposed into V * U'.  This function assigns extra data
    // pointers to interactions in a factorization and assigns
    // the data in V and U to the correct locations.
    //
    // copyWorkspace is used to scatter the rows of V in to full
    // matrix form (since V may have some zero rows)
    //
    // scatterRatio is used to determine whether or not we will store
    // V in full matrix form or in scattered form.  If the scattered
    // size of V is more than scatterRatio times the current size of
    // V, it will be stored in scattered form.
    void assignExtraData( int interaction_idx,
                          std::vector<Supernode> &nodes,
                          const Real *V, const Real *Utrans,
                          int rank, Real *extraData,
                          long int &offset, long int &remainingSize,
                          Real *copyWorkspace );

    // Assigns a single low-rank decomposition for all of this node's
    // off-diagonal content
    void assignExtraData( const Real *V, const Real *Utrans,
                          int rank, Real *extraData,
                          long int &offset, long int &remainingSize );

    // Same as the above, but for a low rank block in the diagonal
    void assignExtraDiagonalData( int block_idx,
                                  std::vector<Supernode> &nodes,
                                  const Real *V, const Real *Utrans,
                                  int rank, Real *extraData,
                                  long int &offset, long int &remainingSize );

    // Initialization for an extended node - to be called after all nodes
    // have their size set
    void initializeExtendedNode( const std::vector<Supernode> &nodes,
                                 int node_idx,
                                 Real *extraData,
                                 long int &offset, long int &remainingSize );

    // Computes the Cholesky factor piece for this node.
    // This assumes that the node has already been set up (eg. via
    // update matrices, etc.)
    void factor( Real *extraData = NULL, int writeSize = -1,
                 WorkspaceManager<Real> *workspaceManager = NULL );

    // Solves for the off-diagonal contents of the Cholesky factor for
    // this node.
    void offDiagonalSolve( bool solveCompressedBlocks,
                           Real *extraData = NULL, int writeSize = -1,
                           WorkspaceManager<Real> *workspaceManager = NULL );

    // Computes dense Cholesky factor(s) for the diagonal block
    void factorDiagonal( Real *extraData,
                         WorkspaceManager<Real> *workspaceManager = NULL );

    // Solves a system using this node's diagonal (assumed to be already
    // factored)
    //
    // Specify whether to transpose the system, and whether the solve
    // should be a left or right solve
    void diagonalSolve( Real *rhs, int nRHS, Real *extraData,
                        bool transpose = false,
                        bool left = true,
                        int startBlock = -1, int endBlock = -1 ) const;

    // Same as the above, but for left-side solves with vector data
    void diagonalSolve( Real *rhs, Real *extraData,
                        bool transpose = false,
                        int startBlock = -1, int endBlock = -1 ) const;

    // Returns the number of rows in this node's off diagonal
    // which will be replaced with low-rank decompositions
    int countLowRankRows() const;

    // Returns the full number of rows in this node's off diagonal
    // which will be replace with low-rank decompositions, assuming
    // we include zero rows.
    int countLowRankRowsFull( const std::vector<Supernode> &nodes ) const;

    // The number of uncompressed entries in the standard part of the
    // system stored by this node
    size_t uncompressedEntries() const;

    // Helps to estimate the additional compressed storage costs this
    // node will introduce, based on the given rank (estimate) for
    // its low rank blocks.
    //
    // In particular, we need the size of an expansion workspace necessary
    // to expand the set of all compressed interactions for any node
    // in the factor.
    void estimateCompressedExpansionStorage(
                                const std::vector<Supernode> &nodes,
                                int rank, IntArray &expansionStorage ) const;

    // Sets which blocks to keep in the main diagonal, and
    // which to treat as low rank
    void setDiagonalBlocks( const std::vector<DenseBlock> &diagonalBlocks,
                            bool compressInPlace = false,
                            int explicitBlockThreshold = -1 );

    // Expands each compressed interaction stored in this node in
    // to the given workspace
    //
    // If diagonalContributions is true, then we only expand the
    // compressed interactions which are needed to contribute to
    // this node's diagonal block.
    void expandCompressedInteractions(
                       const Real *extraData, Real *expansionWorkspace,
                       bool diagonalContributions = false,
                       int maxSz = -1 ) const;

    // For an interaction stored in place in compressed format V * U'
    // where V has orthonormal columns, convert these interactions
    // so that the U now has orthonormal columns (and V is arbitrary).
    void convertCompressedOrthoInteraction( int interaction_idx );

    // Converts all interactions using the function above
    void convertAllCompressedOrthoInteractions();

#if 0
    // For each interaction in this node which makes a contribution
    // to the diagonal, form the given submatrx (specified by block),
    // multiply it by the given matrix G and add to the workspace.
    //
    // expandCompressedInteractions is assumed to have been called
    // previously, with results stored in expansionWorkspace.
    void diagonalContributionMultiply( const Real *extraData,
                                       const Real *expansionWorkspace,
                                       const Real *G, int nColsG,
                                       const DenseBlock &block,
                                       Real *multWorkspaceInitial,
                                       Real *multWorkspaceFinal );
#endif

    // Adds all diagonal contributions to this node resulting from the
    // addition of slack variables.  We require both the extraData
    // array, as well as an "expansion workspace" in which all compressed
    // interactions have been expanded.
    //
    // That is, we expect expandCompressedInteractions to have been called
    void addLowRankDiagonalContributions( const Real *extraData,
                                          Real *expansionWorkspace );

    void addLowRankDiagonalContributions(
                                    const Real *extraData,
                                    WorkspaceManager<Real> &workspaceManager );

#if 0
    // Adds all low rank diagonal contributions by computing the
    // "exact" version of the contribution
    //
    // Each low rank block S must be added in the form S' * S to this
    // node's diagonal with appropriate scaling.
    //
    // Contributions introduced by earlier nodes result in the addition of
    // a scaled identity matrix to the diagonal.
    //
    // Workspaces, descendent lists, and the original matrix are provided
    // so that block columns of each low rank block S can be formed
    void addLowRankDiagonalContributions(
                const SPARSE_MATRIX::SparseColumnMatrix &A,
                const std::vector<Supernode> &nodes,
                const IntArray &lowRankBlocks,
                const std::vector<std::vector<IndexPair> > &lowRankDescendents,
                const IntArray &ancestorInteractions,
                Real *multWorkspace, Real *workspace );
#endif

    // Adds this node's contribution to the extended variable schur
    // complement.
    //
    // Need an expansion workspace to expand any compressed column
    // interactions.
    //
    // multWorkspace is used to extract the needed columns (in the case
    // of compressed column interactions) prior to multiplication
    void addExtendedSchurComplementContribution(
                const std::vector<Supernode> &nodes,
                Real *expansionWorkspace, Real *multWorkspace,
                Real *extraData,
                Timer *invertTimer = NULL,
                Timer *multiplyTimer = NULL );

    void addExtendedSchurComplementContribution(
                const std::vector<Supernode> &nodes,
                WorkspaceManager<Real> &workspaceManager,
                Real *extraData,
                Timer *invertTimer = NULL,
                Timer *multiplyTimer = NULL );

    // Runs a forward solve by solving for the components of
    // x governed by this node's indices, then making appropriate
    // updates to x to account for off-diagonal parts of this node
    void forwardSolve( Real *x,
                       Real *workspace, Real *extendedWorkspace,
                       const std::vector<Supernode> &nodes,
                       Real *extraData,
                       Real *expansionWorkspace );

    inline void forwardSolve( VECTOR &x,
                              Real *workspace, Real *extendedWorkspace,
                              const std::vector<Supernode> &nodes,
                              Real *extraData,
                              Real *expansionWorkspace )
    {
      forwardSolve( x.data(), workspace, extendedWorkspace,
                    nodes, extraData, expansionWorkspace );
    }

    // Same as the above, but solves L' x = b instead
    void backwardSolve( Real *x,
                        Real *workspace, Real *extendedWorkspace,
                        const std::vector<Supernode> &nodes,
                        Real *extraData,
                        Real *expansionWorkspace );

    inline void backwardSolve( VECTOR &x,
                               Real *workspace, Real *extendedWorkspace,
                               const std::vector<Supernode> &nodes,
                               Real *extraData,
                               Real *expansionWorkspace )
    {
      backwardSolve( x.data(), workspace, extendedWorkspace,
                     nodes, extraData, expansionWorkspace );
    }

    //////////////////////////////////////////////////////////////////////
    // Additional solve routines for LDL factorization
    //////////////////////////////////////////////////////////////////////

    // Applies inverse of the LDL factorization block-diagonal, assuming
    // that this factorization exists.  Does nothing for Cholesky-factored
    // nodes.
    void LDLdiagonalSolve( Real *x );

    inline void LDLdiagonalSolve( VECTOR &x )
    {
      LDLdiagonalSolve( x.data() );
    }

    // Applies the permutation matrix from this node's LDL factorization,
    // if it exists.  Does nothing for Chlolesky-factored nodes.
    void applyLDLforwardPermutation( Real *x );

    inline void applyLDLforwardPermutation( VECTOR &x )
    {
      applyLDLforwardPermutation( x.data() );
    }

    // Applies the inverse of the permutation matrix from this node's LDL
    // factorization, if it exists.  Does notthing for Cholesky-factored
    // nodes.
    void applyLDLbackwardPermutation( Real *x );

    inline void applyLDLbackwardPermutation( VECTOR &x )
    {
      applyLDLbackwardPermutation( x.data() );
    }

    //////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////

    // Solves for the components of x governed by this node's indices
    // and makes appropriate updates to x.  This performs a "sub-solve"
    // over only the specified set of nodes.
    //
    // This is used to support the "interior block" machinery.
    void forwardSubSolve( Real *X, int nRHS,
                          Real *workspace,
                          const IntArray &nodeList,
                          const std::vector<Supernode> &nodes,
                          Timer *solveTimer = NULL,
                          Timer *multTimer = NULL );

    // Version with a node range, rather than a list
    void forwardSubSolve( Real *X, int nRHS,
                          Real *workspace,
                          IndexRange nodeRange,
                          const std::vector<Supernode> &nodes );

    // Same as the above, but for backward solves
    void backwardSubSolve( Real *X, int nRHS,
                           Real *workspace,
                           const IntArray &nodeList,
                           const std::vector<Supernode> &nodes,
                           Timer *solveTimer = NULL,
                           Timer *multTimer = NULL );

    // Version with a node range, rather than a node list
    void backwardSubSolve( Real *X, int nRHS,
                           Real *workspace,
                           IndexRange nodeRange,
                           const std::vector<Supernode> &nodes );

    // Returns the index of the first extended node interaction
    // in this node's off-diagonal
    int firstExtendedInteraction() const;

    // Puts the entries from this supernode lying in the given
    // sub matrix range in to S
    //
    // Note: no bounds checking done here
    void copySubMatrix( const IndexRange &rowRange,
                        const IndexRange &columnRange,
                        const Real *extraData,
                        const std::vector<Supernode> &nodes,
                        SPARSE_MATRIX &S );

    // Counts non-zero columns in each interaction
    void countNonZeroInteractionColumns(
                                  const std::vector<Supernode> &nodes,
                                  const SPARSE_MATRIX::SparseColumnMatrix &A );

    // Counts non-zero columns in each low rank diagonal block
    void countNonZeroDiagonalBlockColumns(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A );

    void printInteractions() const
    {
      for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
            interaction_idx++ )
      {
        printf( "Node %d, interaction %d with node %d:\n", _nodeID,
                interaction_idx, _offDiagonal[ interaction_idx ]._nodeID );
        for ( int row_idx = 0;
              row_idx < _offDiagonal[ interaction_idx ]._rowList.size();
              row_idx++ )
        {
          printf( "\t\tRow %04d: %04d\n", row_idx,
                  _offDiagonal[ interaction_idx ]._rowList[ row_idx ] );
        }
        printf( "\n" );
      }
    }

    void writeDiagonal( const char *prefix ) const;

    void writeOffDiagonal( const char *prefix ) const;

    // Adds to a dense sub-matrix containing only the nodes specified
    // in the given list.  This is here for debugging purposes
    void extractSubMatrix( const IntArray &nodeList,
                           const std::vector<Supernode> &nodes,
                           Real *subMatrix, int nRows ) const;

#if 0
    // Adds to a sparse sub-matrix containing only the nodes
    // specified in the given list.  This is here for debugging purposes
    void extractSubMatrix( const IntArray &nodeList,
                           const std::vector<Supernode> &nodes,
                           SPARSE_MATRIX &subMatrix, int nRows ) const;
#endif

    // Writes an extended supernode
    void writeSupernode( const char *prefix, const Real *extraData ) const;

    // Another function necessary for LDL factorization.  Off-diagonal
    // interactions initially store permuted versions of what they should.
    // This function "de-permutes" the content's of descedent's interaction
    // with this node.
    void permuteDescendentInteraction( const Supernode &descendent,
                                       int start_idx,
                                       Real *extraData );

    void initDiagonalBlock()
    {
      if ( _nodeDiagonal.numBlocks() > 0 ) {
        _nodeDiagonal.initData( *this, _diagonalOffsets );
      }
    }

    inline int numColumns() const
    {
      return _columnRange.second - _columnRange.first + 1;
    }

    inline int &numExtendedColumns()
    {
      return _numColumns;
    }

    inline int numExtendedColumns() const
    {
      return _numColumns;
    }

    inline int numRows() const
    {
      return _numRows;
    }

    inline const std::vector<SupernodeInteraction> &offDiagonal() const
    {
      return _offDiagonal;
    }

    inline SupernodeType type() const
    {
      return _type;
    }

    inline int nodeID() const
    {
      return _nodeID;
    }

    inline const IntArray &interactionRowList( int interaction_idx ) const
    {
      return _offDiagonal[ interaction_idx ]._rowList;
    }

    inline int interactionRows( int interaction_idx ) const
    {
      return _offDiagonal[ interaction_idx ]._rowList.size();
    }

    inline void setColumnRange( int start, int end )
    {
      _columnRange.first = start;
      _columnRange.second = end;
    }

    inline int startColumn() const
    {
      return _columnRange.first;
    }

    inline int endColumn() const
    {
      return _columnRange.second;
    }

    inline const vector<DenseBlock> &diagonalBlocks() const
    {
      return _diagonalBlocks;
    }

    inline const vector<DenseBlock> &diagonalLowRankBlocks() const
    {
      return _diagonalLowRankBlocks;
    }

    inline int diagonalLowRankBlockNonZeroColumns( int block_idx ) const
    {
      return _diagonalLowRankBlockNonZeroColumns[ block_idx ];
    }

    inline bool lowRankDiagonal() const
    {
      return ( _diagonalBlocks.size() > 1 );
    }

    inline bool inPlaceDiagonal() const
    {
      return ( _nodeDiagonal.numDiagonalBlocks() > 0 );
    }

    inline int totalRows() const
    {
      return lowRankDiagonal() ? _firstOffDiagonalRow + _numRows
                               : _numRows + numColumns();
    }

    inline Real *offDiagonalData()
    {
      return ( _data + _firstOffDiagonalRow * numColumns() );
    }

    inline Real *data()
    {
      return _data;
    }

    inline int lowRankDiagonalInteraction( int block_idx ) const
    {
      return _compressedDiagonalInteractions[ block_idx ];
    }

    inline Real *interactionMatrix( int interaction_idx )
    {
      Real                  *matrixData = _data;
      int                    nCols = numColumns();

      matrixData += _offDiagonal[ interaction_idx ]._dataRowList[ 0 ] * nCols;

      return matrixData;
    }

    inline const Real *interactionMatrix( int interaction_idx ) const
    {
      Real                  *matrixData = _data;
      int                    nCols = numColumns();

      matrixData += _offDiagonal[ interaction_idx ]._dataRowList[ 0 ] * nCols;

      return matrixData;
    }

    inline int compressOffDiagonal() const
    {
      return _compressOffDiagonal;
    }

    inline int extendedDataOffset() const
    {
      return _extendedDataOffset;
    }

    inline const LDL_Data &ldlData() const
    {
      return _LDL_data;
    }

    inline void setLDL( bool useLDL )
    {
      _useLDL = useLDL;
    }

    inline DiagonalBlock &nodeDiagonal()
    {
      return _nodeDiagonal;
    }

  public:
    typedef boost::function<int (int, int)> RankEstimator;

    // Given a list of supernodes indexed by node IDs, build
    // off diagonal interactions based on non-zero indices
    // in the given list of super nodes.
    //
    // compressInPlace: whether to do in place compression, rather
    //                  than introducing slack variables
    static void BuildOffDiagonalInteractions(
                                std::vector<Supernode> &nodes,
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                bool compressInPlace = false );

    // Takes an existing list of supernodes, and introduces slack
    // variable nodes in order to eliminate any off-diagonal blocks
    // in the factorization beyond a given size.
    //
    // This version of the function introduces slack variables
    // "in line"; that is, slack variables are introduced immediately
    // after the supernode in which low rank interactions were discovered
    static void ExtendSystem( std::vector<Supernode> &nodes,
                              std::vector<Supernode> &newSystem,
                              int maxBlockSize );

    // This function also extends a system, but introduces all slack
    // variables at the end of the system.
    static void ExtendSystemAppend( std::vector<Supernode> &nodes,
                                    std::vector<Supernode> &newSystem,
                                    int maxBlockSize,
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    Real compressionRatio = -1.0,
                                    RankEstimator *estimator = NULL,
                                    bool compressInPlace = false );

    // Performs symbolic reordering of the extended variable space to
    // reduce fill-in.  This modifies the given interaction sets
    static void ReorderExtendedVariables(
                std::vector<Supernode> &nodes,
                std::vector<std::set<SupernodeInteraction> > &interactionSets,
                int oldSystemSize );

    // For each node, examine the row set of each of its interactions
    // in light of the fact that some of its interactions will be
    // removed and replaced with low rank decomposition.  Adjust the
    // symbolic factor to make use of these new, smaller row sets.
    static void ComputeSparsifiedFillIn(
                              std::vector<Supernode> &nodes,
                              int maxBlockSize,
                              const SPARSE_MATRIX::SparseColumnMatrix &A );

    // Given a post-ordering of a nested dissection tree for a
    // matrix, this function builds a symbolic factor for a system
    // extended from A to introduce slack variables.
    static void FactorSymbolic(
          const std::vector<const NestedDissection::DissectionNode *> &nodes,
          const SPARSE_MATRIX::SparseColumnMatrix &A,
          int maxBlockSize,
          std::vector<Supernode> &factor );

    // Given a post-ordering of a nested dissection tree for a
    // matrix, this function builds the basic symbolic factor,
    // in which for now nodes just have the correct column ranges.
    static void BuildBasicFactor(
          const std::vector<const NestedDissection::DissectionNode *> &nodes,
          const SPARSE_MATRIX::SparseColumnMatrix &A,
          std::vector<Supernode> &initialNodes );

    // Same as the above, but just takes a list of supernode specifications
    static void BuildBasicFactor(
          const std::vector<Ordering::SupernodeSpecification> &nodes,
          const SPARSE_MATRIX::SparseColumnMatrix &A,
          std::vector<Supernode> &initialNodes );

    // Counts the number of non-zero entries stored by a list
    // of supernodes
    static long int CountNonZeros( const std::vector<Supernode> &nodes );

    // Allocates all required fixed data for standard nodes, and their
    // interactions with other standard nodes.
    static Real *AllocateFixedData( std::vector<Supernode> &nodes,
                                    long int &dataSize );

#if 0
    // Subtracts the product formed above from another workspace,
    // using relativeMap to map rows from one workspace to the other.
    // This is used when building a matrix to be factored in to
    // a low-rank decomposition.
    static void ScatterLowRankUpdate( const Real *multWorkspaceFinal,
                                      Real *decompWorkspace,
                                      const IntArray &relativeMap,
                                      int nCols );
#endif

    // Similar to the above, but only handles a single interaction.
    // That is, we scatter from one row list to another.
    static void ScatterLowRankUpdate( const Real *multWorkspaceFinal,
                                      Real *decompWorkspace,
                                      const SupernodeInteraction &interaction,
                                      const IntArray &fullRows,
                                      int nCols );

    // Another version of scatterLowRankUpdate.  This scatters a
    // subset of one row list in to a full row space.
    //
    // Used in decomposition of low rank blocks in a node's diagonal
    static void ScatterLowRankUpdate( const Real *multWorkspaceFinal,
                                      Real *decompWorkspace,
                                      const SupernodeInteraction &interaction,
                                      const IndexRange &rowRange,
                                      int nCols,
                                      // Offset to subtract from row indices
                                      int offset );

    // Scatters data from multiple interactions in descendent to
    // the interaction set of ancestor, where ancestor_idx
    // specifies the interaction between descendent and ancestor
    static void ScatterMultiLowRankUpdate( const Real *inputWorkspace,
                                           Real *outputWorkspace,
                                           const Supernode &descendent,
                                           int ancestor_idx,
                                           const Supernode &ancestor,
                                           int nCols );

    // Similar to the above, once we have formed the right half of
    // a low rank decomposition Q * R', we must scatter the columns
    // of the formed matrix R' according to the contents of the
    // interaction between this node and ancestor.
    static void ScatterLowRankTransUpdate( int ancestor_idx,
                                           const Supernode &descendent,
                                           const Supernode &ancestor,
                                           const Real *multWorkspaceFinal,
                                           Real *decompWorkspace,
                                           int nRows );

    // Another version of scatterLowRankTransUpdate.  This scatters a
    // subset of one column list in to a full column space
    //
    // Used in decomposition of low rank blocks in a node's diagonal
    static void ScatterLowRankTransUpdate(
                                    const Real *multWorkspaceFinal,
                                    Real *decompWorkspace,
                                    const SupernodeInteraction &interaction,
                                    const IndexRange &columnRange,
                                    int nRows,
                                    int nColsInput, int nColsOutput,
                                    // Offset to subtract from column indices
                                    int offset );

    // Similar to the above, after a left transpose interaction multiply,
    // we must scatter the rows of multWorkspaceFinal in to the
    // appropriate locations of decompWorkspace
    static void ScatterLowRankTransUpdate_row( int ancestor_idx,
                                               const Supernode &descendent,
                                               const Supernode &ancestor,
                                               const Real *multWorkspaceFinal,
                                               Real *decompWorkspace,
                                               int nCols );

    // Builds a submatrix from fullMatrix using the rows specified
    // in the given interaction (which should be associated with the
    // given node)
    static inline void BuildInteractionSubMatrix(
                                 const SupernodeInteraction &interaction,
                                 const Real *fullMatrix, Real *subMatrix,
                                 int nCols )
    {
      MATRIX::copyRows( fullMatrix, subMatrix, interaction._rowList, nCols );
    }

    // This version also takes a range of indices to consider from the
    // row list
    static inline void BuildInteractionSubMatrix(
                                 const SupernodeInteraction &interaction,
                                 const Real *fullMatrix, Real *subMatrix,
                                 int nCols,
                                 const IndexRange &rowRange,
                                 int offset )
    {
      MATRIX::copyRows( fullMatrix, subMatrix, interaction._rowList,
                        nCols, rowRange, offset );
    }

    // Similar to the above, but now we assume that fullMatrix is itself
    // only defined over some set of rows, of which interaction stores
    // a subset
    static inline void BuildInteractionSubMatrix(
                                 const SupernodeInteraction &interaction,
                                 const IntArray &fullRows,
                                 const Real *fullMatrix, Real *subMatrix,
                                 int nCols );

    // Builds an interaction submatrix in order to multiply by multiple
    // interactions from descendent, assuming that the row space of
    // the input matrix is determined by the row space of ancestor
    static void BuildMultiInteractionSubMatrix(
                                const Supernode &descendent,
                                int ancestor_idx,
                                const Supernode &ancestor,
                                const Real *fullMatrix, Real *subMatrix,
                                int nCols );

    // Reverse of the above.  Scatters rows indexed by interaction._rowList
    // to fullRows in fullMatrix.
    static inline void scatterInteractionMatrix(
                                 const SupernodeInteraction &interaction,
                                 const IntArray &fullRows,
                                 const Real *subMatrix, Real *fullMatrix,
                                 int nCols );

	protected:

  private:
    void factorCholesky( Real *extraData, int writeSize,
                         WorkspaceManager<Real> *workspaceManager );

    // Factors this node's diagonal block using an LDLT factorization,
    // storing the result in _LDL_data
    void factorLDL( Real *extraData,
                    WorkspaceManager<Real> &workspaceManager );

    // Off diagonal solve for Cholesky-factored node
    void choleskyOffDiagonalSolve( bool solveCompressedBlocks,
                                   Real *extraData, int writeSize,
                                   WorkspaceManager<Real> *workspaceManager );

    // Off-diagonal solve for LDL-factored node
    void LDLOffDiagonalSolve( Real *extraData,
                              WorkspaceManager<Real> &workspaceManager );

    // Multiplication by this node's LDL diagonal matrix
    void applyLDLDiagonal( const Real *input, int inputRows, int inputCols,
                           Real *output,
                           bool transpose = false ) const;

#if 0
    // Multplication by this node's LDL diagonal inverse matrix
    //
    // Unfortunately, due to the structure of LU solves, we can't
    // provide a transposed version of this, so applying to a transposed
    // matrix will require some pre-processing
    void applyLDLDiagonalInverse( const Real *input,
                                  int inputRows, int inputCols,
                                  Real *output ) const;
#endif
    // Multplication by this node's LDL diagonal inverse matrix
    //
    // Unfortunately, due to the structure of LU solves, we can't
    // provide a transposed version of this, so applying to a transposed
    // matrix will require some pre-processing
    //
    // Note that this is applied in place (ie. it overwrites the input)
    void applyLDLDiagonalInverse( Real *input,
                                  int inputRows, int inputCols ) const;

    // Helper function for diagonalSolve; performs a solve with this node's
    // main diagonal assuming it has been compressed with extended variables
    void diagonalSolve_extendedVariable( Real *rhs, int nRHS, Real *extraData,
                                         bool transpose, bool left,
                                         int startBlock, int endBlock ) const;

    // Helper function for diagonalSolve; performs a solve with this node's
    // main diagonal assuming it has been compressed in place
    void diagonalSolve_inPlace( Real *rhs, int nRHS,
                                bool transpose, bool left ) const;

    // Version of the above for left solves on vector data
    void diagonalSolve_inPlace( Real *rhs, bool transpose ) const;

    // Assigns an entry to the data array.
    // Parameters are relative row and column indices, and the
    // data to assign
    inline void assign( int row_idx, int col_idx, Real data, bool add = false )
    {
      // FIXME: diagonal sparsification

      // Row-major format.
      // TODO: Check that this is actually what we want
      if ( add )
      {
        _data[ row_idx * numColumns() + col_idx ] += data;
      }
      else
      {
        _data[ row_idx * numColumns() + col_idx ] = data;
      }
    }

    // Sets the _dataRowList array in each of this node's
    // interactions
    void setDataOffsets();

    // Version of copyMatrixData called for standard nodes
    void copyMatrixData_standard( const IntArray &scatteredMap,
                                  const IntArray &lowRankScatteredMap,
                                  const SPARSE_MATRIX::SparseColumnMatrix &A );

    // Version of copyMatrixData called for extended nodes
    void copyMatrixData_extended( Real *extraData );

    // For nodes with a sparsified diagonal
    void copyMatrixData_sparseDiagonal(
                                  const SPARSE_MATRIX::SparseColumnMatrix &A );

    // Version of applyExtendedUpdate to be applied to standard nodes
    void applyExtendedUpdate_standard( const Supernode &descendent,
                                       const PairArray &extendedIndices,
                                       Real *multWorkspace,
                                       Real *extraData,
                                       int start_idx,
                                       Real *expansionWorkspace );

    void applyExtendedUpdate_standard(
                                  const Supernode &descendent,
                                  const PairArray &extendedIndices,
                                  Real *extraData,
                                  int start_idx,
                                  WorkspaceManager<Real> &workspaceManager );

    // Version of applyExtendedUpdate to be applied to extended nodes
    void applyExtendedUpdate_extended( const Supernode &descendent,
                                       const PairArray &extendedIndices,
                                       Real *extraData,
                                       int start_idx,
                                       Real *expansionWorkspace,
                                       Timer *standardTimer,
                                       Timer *extendedTimer );

    // Extended update for a descendent whose factor is stored as
    // an LDL^{T} factorization
    void applyExtendedUpdate_extendedLDL(
                                      const Supernode &descendent,
                                      const PairArray &extendedIndices,
                                      Real *extraData,
                                      int start_idx,
                                      WorkspaceManager<Real> &workspaceManager,
                                      Timer *standardTimer,
                                      Timer *extendedTimer );

    // Identifies the row set associated with the original system in each
    // of this node's interactions
    void buildBaseSystemRowSets(
                           const std::vector<Supernode> &nodes,
                           std::vector<std::set<int> > &rowSets,
                           const SPARSE_MATRIX::SparseColumnMatrix &A ) const;

    // Based on this node's diagonal set, figure out where the
    // interaction data should be placed
    void setDiagonalSize();

    // Recursively build the set of dense blocks to be decomposed
    // in the diagonal, based on the set of diagonal blocks
    void buildLowRankDiagonalBlocks( int start_idx, int end_idx,
                                     int parent = -1 );

#if 0
    // For a given piece of a low rank decomposition, add its outer
    // product to the diagonal of this node
    //
    // eg. diagonal += U' * U
    void addLowRankDiagonalContribution( const Real *U, int rank,
                                         bool transpose = false );
#endif

    // For a given piece of a low rank decomposition, add its outer
    // product to the diagonal of this node
    //
    // eg. diagonal += U' * U
    void addLowRankDiagonalContribution(
                                  const Real *U, int rank,
                                  Real *multWorkspace,
                                  const IntArray *compressedColumnList = NULL,
                                  bool transpose = false );

    void addLowRankDiagonalContribution(
                                  const Real *U, int rank,
                                  WorkspaceManager<Real> &workspaceManager,
                                  const IntArray *compressedColumnList = NULL,
                                  bool transpose = false );

    // Versions of the above for compressed matrices
    void addLowRankDiagonalContribution_compressed(
                                  const Real *U, int rank,
                                  Real *multWorkspace,
                                  const IntArray &compressedColumnList,
                                  bool transpose = false );

    void addLowRankDiagonalContribution_compressed(
                                  const Real *U, int rank,
                                  WorkspaceManager<Real> &workspaceManager,
                                  const IntArray &compressedColumnList,
                                  bool transpose = false );

    // Version for uncompressed matrices
    void addLowRankDiagonalContribution_uncompressed(
                                  const Real *U, int rank,
                                  bool transpose = false );

    // For a given interaction in this node's off-diagonal, apply the
    // inverse of this node's diagonal block to the right hand side
    // of the interaction block.  (ie. form S * L^{-T} * L^{-1} where
    // S is the interaction block)
    //
    // If the interaction is compressed, expand it before applying the
    // inverse.
    //
    // If it is compressed because it comes from a low-rank diagonal block,
    // then only form the necessary part of the matrix.
    //
    // This is a helper function for addExtendedSchurComplementContribution
    void applyInverseToInteraction( int interaction_idx,
                                    Real *extraData,
                                    Real *expansionWorkspace );

    // Multiplies the matrix stored for the given interaction with the
    // transpose of the provided matrix, which may only be defined over
    // a particular column range (stored in expansionWorkspace).
    //
    // Subtract the result from the provided output matrix
    //
    // This is a helper function for addExtendedSchurComplementContribution
    void multiplyInteraction( int interaction_idx,
                              const Real *expansionWorkspace,
                              int numRows, const IndexRange &columnRange,
                              Real *extraData, Real *multWorkspace,
                              Real *outputMatrix );

    // Version of the above function where the interaction arises from
    // a low rank diagonal block
    void multiplyInteractionLowRankDiagonal(
                              const SupernodeInteraction &interaction,
                              const Real *expansionWorkspace,
                              int numRows, const IndexRange &columnRange,
                              Real *extraData, Real *multWorkspace,
                              Real *outputMatrix );

    // Version of the above function where the interaction is stored
    // in compressed column form
    void multiplyInteractionCompressed(
                              const SupernodeInteraction &interaction,
                              const Real *expansionWorkspace,
                              int numRows, const IndexRange &columnRange,
                              Real *extraData, Real *multWorkspace,
                              Real *outputMatrix );

    // Version of the above for a "regular" uncompressed interaction
    void multiplyInteractionUncompressed(
                              const SupernodeInteraction &interaction,
                              const Real *expansionWorkspace,
                              int numRows, const IndexRange &columnRange,
                              Real *extraData, Real *multWorkspace,
                              Real *outputMatrix );

    // Builds a map from indices in this node's column range to the
    // diagonal block handling those indices
    void buildDiagonalMap();

#if 0
    // Helper function for diagonalContributionMult
    //
    // Does multiplication with contributions resulting from the
    // compression of off diagonal interactions
    void diagonalContributionMult_offDiagonalContributions(
                                               int block_idx,
                                               const Real *G, int nColsG,
                                               const Real *extraData,
                                               const Real *expansionWorkspace,
                                               Real *multWorkspaceInitial,
                                               Real *multWorkspaceFinal,
                                               bool transpose,
                                               bool left );
#endif

    // Helper function for diagonalContributionMult
    //
    // Does multiplication with contributions resulting from the
    // compression of off diagonal interactions
    void diagonalContributionMult_offDiagonalContributions(
                                               int block_idx,
                                               const Real *G, int nColsG,
                                               const Real *extraData,
                                               Real *multWorkspaceInitial,
                                               Real *multWorkspaceFinal,
                                               Real *copyWorkspace,
                                               bool transpose,
                                               bool left );

    void diagonalContributionMult_offDiagonalContributions(
                                      int block_idx,
                                      const Real *G, int nColsG,
                                      const Real *extraData,
                                      WorkspaceManager<Real> &workspaceManager,
                                      Real *multWorkspaceFinal,
                                      bool transpose,
                                      bool left );

    // Helper function for diagonalContributionMult
    //
    // Does multiplication with contributions resulting from the
    // compression of blocks from this node's diagonal
    void diagonalContributionMult_diagonalContributions(
                                               int block_idx,
                                               const Real *G, int nColsG,
                                               const Real *extraData,
                                               Real *multWorkspaceInitial,
                                               Real *multWorkspaceFinal,
                                               bool transpose,
                                               bool left );

    void diagonalContributionMult_diagonalContributions(
                                      int block_idx,
                                      const Real *G, int nColsG,
                                      const Real *extraData,
                                      WorkspaceManager<Real> &workspaceManager,
                                      Real *multWorkspaceFinal,
                                      bool transpose,
                                      bool left );

#if 0
    // Helper function for addLowRankDiagonalContributions
    //
    // Adds contributions resulting from compression of off-diagonal
    // interactinos
    void addLowRankDiagonalContributions_offDiagonalContributions(
                                            const Real *extraData,
                                            const Real *expansionWorkspace );
#endif

    // Helper function for addLowRankDiagonalContributions
    //
    // Adds contributions resulting from compression of off-diagonal
    // interactinos
    void addLowRankDiagonalContributions_offDiagonalContributions(
                                            const Real *extraData,
                                            Real *expansionWorkspace );

    void addLowRankDiagonalContributions_offDiagonalContributions(
                                    const Real *extraData,
                                    WorkspaceManager<Real> &workspaceManager );

    // Helper function for addLowRankDiagonalContributions
    //
    // Adds contributions resulting from compression of blocks from
    // this node's diagonal
    void addLowRankDiagonalContributions_diagonalContributions(
                                            const Real *extraData );

    // Function which assigns low rank decomposition information to
    // a new slack variable node introduced to compensate for the
    // given interaction
    void assignNewExtendedNode( int interaction_idx,
                                vector<Supernode> &nodes,
                                const Real *V, const Real *Utrans,
                                int rank, Real *extraData,
                                long int &offset, long int &remainingSize,
                                Real *copyWorkspace );

    // Alternative to the function above, in which we simply assign the
    // low-rank decomposition directly to a compressed interaction
    // (ie. in place compression)
    void assignCompressedInteraction( int interaction_idx,
                                      const Real *V, const Real *Utrans,
                                      int rank, Real *extraData,
                                      long int &offset,
                                      long int &remainingSize );

    // Function which assigns low rank diagonal decomposition information
    // to a new slack variable node introduced to compensate for the
    // removal of a block in the diagonal
    void assignNewExtendedNode_diagonal( int block_idx,
                                         std::vector<Supernode> &nodes,
                                         const Real *V, const Real *Utrans,
                                         int rank, Real *extraData,
                                         long int &offset,
                                         long int &remainingSize );

    // Alternative to the function above, in which we simply assign the
    // low rank-decomposition of a block in the diagonal directly to a
    // compressed interaction (ie. in place compression)
    void assignCompressedInteraction_diagonal( int block_idx,
                                               const Real *V,
                                               const Real *Utrans,
                                               int rank, Real *extraData,
                                               long int &offset,
                                               long int &remainingSize );

    //////////////////////////////////////////////////////////////////////
    // Functions for performing operations on interactions stored
    // "in place"
    //////////////////////////////////////////////////////////////////////

    // Multiplies the given matrix with an interaction stored in
    // compressed form.  Left and right multiplication, as well as
    // transposed multiplication are supported.
    //
    // In certain cases, we need to know whether or not columns/rows
    // of G need to be extracted to match the row pattern of this
    // interaction.  By default, we assume that they do.
    //
    // Optional parameters are necessary if we only wish to multiply
    // by a certain subset of this interaction (useful when forming a
    // low-rank block in a node's diagonal, for instance)
    void compressedInteractionMultiply(
                                  int interaction_idx,
                                  const Real *G,
                                  int rowsG, int colsG,
                                  Real *output,
                                  bool left, bool transpose,
                                  bool extractSubMatrix = true,
                                  IndexRange rowRange = IndexRange( -1, -1 ),
                                  int rowOffset = 0 ) const;

    // Same as the above, but actually multiplies with the full compressed
    // interaction matrix stored by this node, only for interactions following
    // the provided ancestor_idx
    void compressedInteractionMultiply_fullBlock(
                                  int ancestor_idx,
                                  const Real *G,
                                  int rowsG, int colsG,
                                  Real *output,
                                  bool left, bool transpose,
                                  bool extractSubMatrix = true );

    // Inefficient compressedInteractionMultiply used for debugging
    // purposes
    void compressedInteractionMultiplyDebug(
                                  int interaction_idx,
                                  const Real *G,
                                  int rowsG, int colsG,
                                  Real *output,
                                  bool left, bool transpose,
                                  bool extractSubMatrix = true,
                                  IndexRange rowRange = IndexRange( -1, -1 ),
                                  int rowOffset = 0 ) const;

#if 0
    // Adds the exact low rank diagonal contribution for a given low
    // rank block from this node
    //
    //  lowRankBlocks: contains the list of all sparsified block indices
    //                 in this node
    //  block_idx: Index in to lowRankBlocks of the particular block to
    //             be handled here
    //  lowRankDescendents: Of the same size of lowRankBlocks, the list
    //                      of descendents affecting each low rank block
    //  ancestorInteractions: If a node nodeID is a descendent of this node,
    //                        then ancestorInteractions[ nodeID ] gives the
    //                        index of the interaction between nodeID and
    //                        this node, in nodeID's off diagonal.
    //  scale: Scaling to apply to diagonal contribution
    void addLowRankDiagonalContribution(
                const SPARSE_MATRIX::SparseColumnMatrix &A,
                const std::vector<Supernode> &nodes,
                const IntArray &lowRankBlocks, int block_idx,
                const std::vector<std::vector<IndexPair> > &lowRankDescendents,
                const IntArray &ancestorInteractions,
                Real scale, Real *multWorkspace, Real *workspace );

    // Simpler version of the above for cases in which we just need to
    // add a scaled identity matrix to the diagonal
    void addLowRankDiagonalContribution( Real scale );
#endif

#if 0
    // Forms a block column of an interaction in this node, given
    // a descendent list
    //
    // See parameters in version of addLowRankDiagonalContribution above
    void formInteractionBlockColumn(
                const SPARSE_MATRIX::SparseColumnMatrix &A,
                const std::vector<Supernode> &nodes,
                int interaction_idx,
                const std::vector<IndexPair> &descendents,
                const IntArray &ancestorInteractions,
                const IndexRange &columnRange,
                Real *multWorkspace, Real *workspace );
#endif

  private:
    // For each node in a list, identifies which entries in the
    // original matrix contribute off-diagonal entries in the
    // supernode and adds the corresponding row indices to a
    // row set for that node.
    static void findSystemRows( const std::vector<Supernode> &nodes,
                                std::vector<std::set<int> > &rowSets,
                                const SPARSE_MATRIX::SparseColumnMatrix &A );

    // Based on the current set of rows contained in the off-diagonal
    // of each node in the nodes list, determine how additional rows
    // are introduced by interacting supernodes.
    static void propagateFillIn( const std::vector<Supernode> &nodes,
                                 std::vector<std::set<int> > &rowSets,
                                 const IntArray &nodeMap,
                                 bool compressInPlace );

    // Builds an map from node indices to super node indices
    static void buildSupernodeMap( const std::vector<Supernode> &nodes,
                                   IntArray &nodeMap );

    // Based on a set of off-diagonal rows for each supernode, build
    // all off diagonal interactions in a node list
    static void buildSupernodeInteractions(
                                std::vector<Supernode> &nodes,
                                const std::vector<std::set<int> > &rowSets,
                                IntArray &nodeMap );

    // Checks the given supernode to see whether or not parts of
    // its off-diagonal need to be sparsified and introduces slack
    // variable supernodes if so.
    static void addNode( Supernode &node,
                         const std::vector<Supernode> &nodes,
                         std::vector<Supernode> &newSystem,
                         int maxBlockSize );

    // Similar to the function above, this function checks the given
    // supernode for off-diagonal blocks larger than the given size,
    // and introduces slack variables in the place of these nodes.
    // However, in this case, all slack variables are appended to the
    // end of the system.
    static void addNodeAppend(
                int node_idx,
                std::vector<Supernode> &nodes,
                int maxBlockSize,
                std::vector<std::set<SupernodeInteraction> > &interactionSets,
                Real compressionRatio,
                RankEstimator *estimator,
                bool compressInPlace );

    // Adds new nodes/interactions in order to compress off-diagonal blocks
    // in the main diagonal of this node
    static void compressDiagonalSymbolic(
                int node_idx,
                std::vector<Supernode> &nodes,
                std::vector<std::set<SupernodeInteraction> > &interactionSets,
                RankEstimator *estimator );

    // Having added extended nodes to the interaction set for a
    // particular node, propagate the fill-in resulting from
    // these nodes to all nodes later in the factorization.
    static void propagateExtendedFillIn(
                const Supernode &node,
                std::vector<Supernode> &nodes,
                std::vector<std::set<SupernodeInteraction> > &interactionSets );

    // Adds all extended interactions from the given set to the system.
    //
    // Also sets up forward interactions amongst extended node interactions.
    // This helps us allocate storage space during numerical factorization.
    static void addExtendedInteractions(
                std::vector<Supernode> &nodes,
                std::vector<std::set<SupernodeInteraction> > &interactionSets,
                int originalSize );

    // Helper function for factorLDL which extracts diagonal data from
    // the factors and places the relevant blocks and their inverses in
    // to the provided datastructure
    static void ExtractLDLDiagonal( const Real *A, int rows,
                                    LDL_Data &factorData );

	private:
    static const int                 EMPTY;

    SupernodeType                    _type;
    int                              _nodeID;

    // Column indices from the original matrix.
    // These are not used if this is an extended node.
    IndexRange                       _columnRange;

    // For extended nodes, store a size estimate (used for symbolic
    // reordering of extended nodes)
    int                              _sizeEstimate;

    // Models interactions below the diagonal with
    // other supernodes
    vector<SupernodeInteraction>     _offDiagonal;

    int                              _extensionStart;

    // Total off-diagonal interaction rows
    int                              _numRows;

    // Number of columns - only used for extended nodes
    int                              _numColumns;

    // Pointer to this node's data
    Real                            *_data;

    bool                             _compressOffDiagonal;

    // Offset in to an extra data array to store the diagonal
    // portion of this node, if it is an extended node.
    long int                         _extendedDataOffset;

    //////////////////////////////////////////////////////////////////////
    // Optional LDL^{T} factorization data for extended nodes
    //////////////////////////////////////////////////////////////////////
    bool                             _useLDL;

    LDL_Data                         _LDL_data;

    //////////////////////////////////////////////////////////////////////
    // Storage information for sparsification of the main diagonal block
    //////////////////////////////////////////////////////////////////////

    // List of blocks on the diagonal to keep around
    vector<DenseBlock>               _diagonalBlocks;

    // Data offsets for each diagonal block
    vector<long int>                 _diagonalOffsets;

    // Maps rows/columns in this node's diagonal back to which diagonal
    // block they are contained in
    IntArray                         _diagonalMap;

    // List of off-diagonal blocks to be sparsified
    vector<DenseBlock>               _diagonalLowRankBlocks;

    // Nonzero columns in the original system in each low rank block
    // on the diagonal
    IntArray                         _diagonalLowRankBlockNonZeroColumns;

    // Interaction indices for slack variables arising from diagonal
    // block compression
    IntArray                         _compressedDiagonalInteractions;

    // In order to compress off diagonal blocks, we have to account
    // for contributions made to the diagonal as a result of the
    // introduction of slack variables.  These contributions may include
    // previously compressed diagonal block contributions.
    //
    // For each block, we store a list of other low rank blocks in the
    // diagonal that can contribute to it, as well as row and column
    // ranges relative to that block that we must consider.
    std::vector<std::vector<std::pair<int, DenseBlock> > >
                                     _previousDiagonalContributions;

    // Also store the "parent" of each off diagonal block
    IntArray                         _diagonalLowRankBlockParents;

    // If we are not storing the diagonal in "normal" form,
    // then we need to know at which local row we start storing
    // off-diagonal stuff
    int                              _firstOffDiagonalRow;

    // List of interactions whose outer products have to be
    // added to the diagonal
    IntArray                         _diagonalContributions;

    // List of interactions and column ranges that have to
    // be added to the diagonal.  This list only contains contributions
    // from interactions in the main diagonal that have been compressed
    std::vector<std::vector<std::pair<int, IndexRange> > >
                                     _mainDiagonalContributions;

    // If we are doing in-place compression, we will use the DiagonalBlock
    // data structure to represent the block
    DiagonalBlock                    _nodeDiagonal;

    //////////////////////////////////////////////////////////////////////
    // These variables are used if we want to compress this node's
    // diagonal "in place"; that is, without introducing slack variables
    //////////////////////////////////////////////////////////////////////
    bool                             _inPlaceDiagonalCompression;

    // Sizes of the compressed off diagonal blocks referred to by
    // _diagonalLowRankBlocks
    IntArray                         _diagonalLowRankBlockSizes;

    // Offsets in to the extended data array for each compressed
    // diagonal block
    LongIntArray                     _diagonalLowRankBlockOffsets;

    //////////////////////////////////////////////////////////////////////
    // These variables are used if we want to represent all of a node's
    // off-diagonal data using a single low-rank decomposition, rather
    // than per-interaction bases.  That is, the off diagonal block for
    // this node is approximated by (_compressedV * _compressedU')
    //////////////////////////////////////////////////////////////////////
    Real                            *_compressedV;
    Real                            *_compressedU;

    int                              _compressedRank;

    // FIXME: debugging
    friend int main( int argc, char **argv );

};

#endif
