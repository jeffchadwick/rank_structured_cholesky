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
// Supernode.cpp: Implementation of the Supernode class
//
//////////////////////////////////////////////////////////////////////

#include "Supernode.h"

#include <set>

#include <rschol/ordering/MinimumDegree.h>
#include <rschol/ordering/Ordering.h>

#include <rschol/util/IO.h>
#include <rschol/util/StatsCounter.h>
#include <rschol/util/STLUtil.h>

// FIXME:
//#include <mkl_types.h>
//#include <mkl_cblas.h>

#ifdef USE_MKL
#include <mkl_types.h>
#include <mkl_cblas.h>
#else
// Just use standard blas libraries
extern "C" {
#include <cblas.h>
}
#endif

using namespace std;

// FIXME
int Supernode::numLargeBlocks = 0;
long int Supernode::multWorkspaceSz = 0;
long int Supernode::workspaceSz = 0;

const int Supernode::EMPTY = -1;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
Supernode::Supernode( SupernodeType type )
  : _type( type ),
    _nodeID( -1 ),
    _columnRange( -1, -2 ), // So that numColumns() == 0
    _sizeEstimate( 0 ),
    _data( NULL ),
    _numRows( 0 ),
    _numColumns( 0 ),
    _compressOffDiagonal( false ),
    _useLDL( false ),
    _firstOffDiagonalRow( -1 ),
    _inPlaceDiagonalCompression( false ),
    _compressedV( NULL ),
    _compressedU( NULL ),
    _compressedRank( -1 )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
Supernode::~Supernode()
{
}

//////////////////////////////////////////////////////////////////////
// Returns the fixed floating point memory requirement
// for this node.  That is, floating point requirements
// based on the node's diagonal entry, and interactions
// with other standard nodes.
//////////////////////////////////////////////////////////////////////
int Supernode::numFixedEntries() const
{
  int                            numColumns;
  int                            numFixedRows = 0;

  int                            diagonalEntries = 0;
  int                            blockSz;

  numColumns = _columnRange.second - _columnRange.first + 1;

  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    // Only count active interactions with standard nodes (those
    // not referring to slack variables)
    if ( _offDiagonal[ interaction_idx ]._type == STANDARD_NODE
      && _offDiagonal[ interaction_idx ]._active )
    {
      numFixedRows += _offDiagonal[ interaction_idx ]._rowList.size();
    }
  }

  diagonalEntries = numColumns * _firstOffDiagonalRow;

  return ( numColumns * numFixedRows + diagonalEntries );
}

//////////////////////////////////////////////////////////////////////
// Flags an interaction as implicit
//////////////////////////////////////////////////////////////////////
void Supernode::flagImplicitInteraction( int interaction_idx )
{
  TRACE_ASSERT( _offDiagonal[ interaction_idx ]._type == STANDARD_NODE
             && _offDiagonal[ interaction_idx ]._active );

  _offDiagonal[ interaction_idx ]._type = IMPLICIT_NODE;

  _numRows -= _offDiagonal[ interaction_idx ]._rowList.size();
}

//////////////////////////////////////////////////////////////////////
// Builds a map from matrix indices to row indices in this
// node's fixed data array.  This is necessary so that we
// can copy entries from the original matrix in to the node.
//
// The input array is assumed to have the correct size
//////////////////////////////////////////////////////////////////////
void Supernode::constructScatteredMap( IntArray &scatteredMap,
                                       IntArray &lowRankScatteredMap,
                                       const vector<Supernode> &nodes )
{
  int                            dataRow = 0;
  int                            lowRankRow = 0;
  int                            startColumn;

  // Start by adding map entries for this node's diagonal block.
  // This is only done for unsparsified diagonal blocks.
  if ( _type == STANDARD_NODE )
  {
    if ( !lowRankDiagonal() )
    {
      for ( int row_idx = 0; row_idx < numColumns(); row_idx++ )
      {
        scatteredMap[ _columnRange.first + row_idx ] = dataRow;

        dataRow++;
      }
    }
    else
    {
      dataRow = _firstOffDiagonalRow;
    }
  }

  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    if ( _offDiagonal[ interaction_idx ]._type == STANDARD_NODE
      && _offDiagonal[ interaction_idx ]._active )
    {
      // Add entries to the scattered map for these indices
      const Supernode &node = nodes[ _offDiagonal[ interaction_idx ]._nodeID ];
      const IntArray &rowList = _offDiagonal[ interaction_idx ]._rowList;

      startColumn = node._columnRange.first;

      // rowList stored relative offsets
      for ( int row_idx = 0; row_idx < rowList.size(); row_idx++ )
      {
        scatteredMap[ startColumn + rowList[ row_idx ] ] = dataRow;

        dataRow++;
      }
    }
    else if ( (_offDiagonal[ interaction_idx ]._type == STANDARD_NODE
               && !_offDiagonal[ interaction_idx ]._active)
           || _offDiagonal[ interaction_idx ]._type == COMPRESSED_NODE )
    {
      // Add entries to the low-rank scattered map for these indices
      const Supernode &node = nodes[ _offDiagonal[ interaction_idx ]._nodeID ];
      const IntArray &rowList = _offDiagonal[ interaction_idx ]._rowList;

      startColumn = node._columnRange.first;

      // rowList stores relative offsets
      for ( int row_idx = 0; row_idx < rowList.size(); row_idx++ )
      {
        lowRankScatteredMap[ startColumn + rowList[ row_idx ] ] = lowRankRow;

        lowRankRow++;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Clears any changes made by this node to the scattered map
//////////////////////////////////////////////////////////////////////
void Supernode::clearScatteredMap( IntArray &scatteredMap,
                                   IntArray &lowRankScatteredMap,
                                   const vector<Supernode> &nodes )
{
  int                            startColumn;

  // Start by adding map entries for this node's diagonal block.
  if ( _type == STANDARD_NODE )
  {
    for ( int row_idx = 0; row_idx < numColumns(); row_idx++ )
    {
      scatteredMap[ _columnRange.first + row_idx ] = -1;
    }
  }

  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    // Add entries to the relative map for these indices
    if ( _offDiagonal[ interaction_idx ]._type == STANDARD_NODE
      && _offDiagonal[ interaction_idx ]._active )
    {
      const Supernode &node = nodes[ _offDiagonal[ interaction_idx ]._nodeID ];
      const IntArray &rowList = _offDiagonal[ interaction_idx ]._rowList;

      startColumn = node._columnRange.first;

      // rowList stored relative offsets
      for ( int row_idx = 0; row_idx < rowList.size(); row_idx++ )
      {
        scatteredMap[ startColumn + rowList[ row_idx ] ] = -1;
      }
    }
    else if ( (_offDiagonal[ interaction_idx ]._type == STANDARD_NODE
               && !_offDiagonal[ interaction_idx ]._active)
           || _offDiagonal[ interaction_idx ]._type == COMPRESSED_NODE )
    {
      // Add entries to the low-rank scattered map for these indices
      const Supernode &node = nodes[ _offDiagonal[ interaction_idx ]._nodeID ];
      const IntArray &rowList = _offDiagonal[ interaction_idx ]._rowList;

      startColumn = node._columnRange.first;

      // rowList stores relative offsets
      for ( int row_idx = 0; row_idx < rowList.size(); row_idx++ )
      {
        lowRankScatteredMap[ startColumn + rowList[ row_idx ] ] = -1;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Copies entries from the original matrix to the appropriate
// location in this node's data array (does nothing for an
// extended node).
//
// Assumes the scattered map has been constructed for this node
//////////////////////////////////////////////////////////////////////
void Supernode::copyMatrixData( const IntArray &scatteredMap,
                                const IntArray &lowRankScatteredMap,
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Real *extraData )
{
  if ( _type == STANDARD_NODE )
  {
    copyMatrixData_standard( scatteredMap, lowRankScatteredMap, A );
  }
#if 0
  else
  {
    copyMatrixData_extended( extraData );
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Copies a column subset of matrix data for the given interaction
// between this node and ancestor.
// 
// Data is scattered according to the row list in this interaction
//////////////////////////////////////////////////////////////////////
void Supernode::copyMatrixDataBlockColumn(
                                int interaction_idx,
                                const Supernode &ancestor,
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                const IndexRange &columnRange,
                                Real *workspace )
{
  int                        startColumn;
  int                        endColumn;

  int                        row_idx;

  const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

  int                        startRow;
  int                        endRow;

  int                        nCols = range_size( columnRange );

  const IntArray            &rowList = interaction._rowList;

  TRACE_ASSERT( interaction._nodeID == ancestor.nodeID(), "Node ID mismatch" );

  startRow = ancestor._columnRange.first;
  endRow = ancestor._columnRange.second;

  startColumn = _columnRange.first + columnRange.first;
  endColumn = _columnRange.first + columnRange.second;

  // Start by clearing the output matrix
  MATRIX::clear( workspace, interaction._rowList.size(), nCols );

  // FIXME: debugging
  if ( interaction._rowList.size() * nCols > workspaceSz )
  {
    TRACE_ASSERT( NULL, "Workspace too small" );
  }

  for ( int col_idx = startColumn; col_idx <= endColumn; col_idx++ )
  {
    int                      data_row_idx = 0;

    for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
          row_ptr++ )
    {
      row_idx = A._i[ row_ptr ];

      if ( row_idx < startRow || row_idx > endRow )
      {
        continue;
      }
#if 0
      else if ( row_idx > endRow )
      {
        break;
      }
#endif
      
      row_idx -= startRow;

      // Find the corresponding row in the interaction row list
      while ( data_row_idx < rowList.size()
           && rowList[ data_row_idx ] != row_idx )
      {
        data_row_idx++;
      }

      TRACE_ASSERT( data_row_idx < rowList.size(), "No matching row found" );

      // Assign the sparse matrix entry to the desired workspace location
      MATRIX::access( workspace,
                      rowList.size(), nCols, /* Output dimensions */
                      data_row_idx, col_idx - startColumn ) = A._x[ row_ptr ];
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Copies matrix data for the node in _nodeDiagonal rooted at
// block_idx in to the provided matrix
//////////////////////////////////////////////////////////////////////
void Supernode::copyDiagonalBlockMatrixData(
                                int block_idx,
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                Real *workspace ) const
{
  IndexRange                 rowRange;
  const DenseBlock          &block = _nodeDiagonal.fullBlock( block_idx );
  
  rowRange = IndexRange( _columnRange.first + block._rowRange.first,
                         _columnRange.first + block._rowRange.second );

  for ( int col_idx = rowRange.first; col_idx <= rowRange.second; col_idx++ )
  for ( int row_ptr = A._p[col_idx]; row_ptr < A._p[col_idx + 1]; row_ptr++ ) {
    int                      row_idx = A._i[ row_ptr ];

    if ( in_range( rowRange, row_idx ) ) {
      MATRIX::access( workspace, block.numRows(), block.numColumns(),
                      row_idx - rowRange.first,
                      col_idx - rowRange.first ) = A._x[ row_ptr ];
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Finds a list of mutual off-diagonal blocks shared by this
// node and a given descendent.
//
// We must be provided with the index in to d._offDiagonal
// corresponding to "this" node
//
// lowRankIndices will return interaction index pairs for
// interactions from descendent that should be accumulated into
// a low rank block of *this* node.
//////////////////////////////////////////////////////////////////////
void Supernode::findInteractionIntersection( const Supernode &descendent,
                                             int start_idx,
                                             PairArray &indices,
                                             PairArray &lowRankIndices,
                                             PairArray &extendedIndices )
{
  int                            this_idx = 0;
  int                            descendent_idx;
  int                            old_idx = -1;

  indices.clear();
  lowRankIndices.clear();
  extendedIndices.clear();

  const SupernodeInteraction    &baseInteraction
                                  = descendent._offDiagonal[ start_idx ];

  // Iterate over each list of interactions and find pairwise,
  // shared interactions.
  //
  // We start at start_idx + 1, because the interaction at start_idx
  // corresponds to this node specifically
  for ( descendent_idx = start_idx + 1;
        descendent_idx < descendent._offDiagonal.size();
        descendent_idx++ )
  {
    const SupernodeInteraction  &interaction
                                  = descendent._offDiagonal[ descendent_idx ];

    // If this interaction and the base interaction are extended interactions
    // that have have compressed column ranges that do not intersect, thenv
    // they will not contribute to fill-in
    if ( baseInteraction._lowRankDiagonalInteraction
      && interaction._lowRankDiagonalInteraction
      && !range_overlap( baseInteraction._compressedColumnRange,
                         interaction._compressedColumnRange ) )
    {
      continue;
    }

    // Skip inactive interactions
    if ( !interaction._active )
    {
      continue;
    }

    // If this interaction is active, add it to the list
    while ( this_idx < _offDiagonal.size()
         && _offDiagonal[ this_idx ]._nodeID != interaction._nodeID )
    {
      this_idx++;
    }

    TRACE_ASSERT( this_idx != old_idx, "Repeated interaction!" );

    // Make sure we actually find this interaction.
    TRACE_ASSERT( this_idx < _offDiagonal.size(),
                  "Missing interaction in current node" );

    if ( interaction._type == COMPRESSED_NODE ) {
      TRACE_ASSERT( _offDiagonal[ this_idx ]._type == COMPRESSED_NODE,
                    "Invalid interaction type" );
    }

    // If our interaction is active, we add it to the normal
    // list, otherwise, to the low rank list
    if ( !_offDiagonal[ this_idx ]._active
      || _offDiagonal[ this_idx ]._type == COMPRESSED_NODE )
    {
      lowRankIndices.push_back( IndexPair( descendent_idx, this_idx ) );
    }
    else if ( _offDiagonal[ this_idx ]._active )
    {
      // Should keep track of both standard nodes and compressed
      // nodes here
      if ( _offDiagonal[ this_idx ]._type == STANDARD_NODE
        || _offDiagonal[ this_idx ]._type == COMPRESSED_NODE )
      {
        indices.push_back( IndexPair( descendent_idx, this_idx ) );
      }
      else if ( _offDiagonal[ this_idx ]._type == EXTENDED_NODE )
      {
        TRACE_ASSERT( interaction._type == EXTENDED_NODE,
                      "Interaction type mismatch" );

        extendedIndices.push_back( IndexPair( descendent_idx, this_idx ) );
      }
      // Skip implicit interactions
    }
    else
    {
      TRACE_ASSERT( NULL, "Invalid interaction type" );
    }

    old_idx = this_idx;
  }
}

//////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////
void Supernode::buildRelativeMap( const Supernode &descendent,
                                  int start_idx,
                                  const PairArray &indices,
                                  IntArray &relativeMap )
{
  relativeMap.clear();

  // Start by looking at the interaction in descdendent with
  // this node (which must exist)
  TRACE_ASSERT( descendent._offDiagonal[ start_idx ]._nodeID == _nodeID,
                "Node ID mismatch" )

  {
    const SupernodeInteraction  &interaction
                                  = descendent._offDiagonal[ start_idx ];

    const IntArray              &rowList = interaction._rowList;

    // We assume that in this node, the diagonal block is full
    for ( int row_idx = 0; row_idx < rowList.size(); row_idx++ )
    {
      relativeMap.push_back( rowList[ row_idx ] );
    }
  }

  // Iterate over each interaction in the intersection and figure
  // out how rows in descendent map to rows in this
  for ( int interaction_idx = 0; interaction_idx < indices.size();
        interaction_idx++ )
  {
    const IndexPair             &p = indices[ interaction_idx ];

    const SupernodeInteraction  &interaction
                                  = descendent._offDiagonal[ p.first ];

    const IntArray              &rowList = interaction._rowList;

    const SupernodeInteraction  &thisInteraction = _offDiagonal[ p.second ];

    const IntArray              &thisRowList = thisInteraction._rowList;
    const IntArray              &thisDataRowList = thisInteraction._dataRowList;

    int                          thisRowIdx = 0;

    TRACE_ASSERT( interaction._nodeID == thisInteraction._nodeID,
                  "Interaction node ID mismatch" );

    for ( int row_idx = 0; row_idx < rowList.size(); row_idx++ )
    {
      // Find the corresponding row in the interaction from
      // this node (must exist)
      while ( thisRowList[ thisRowIdx ] != rowList[ row_idx ]
           && thisRowIdx < thisRowList.size() )
      {
        thisRowIdx++;
      }

      TRACE_ASSERT( thisRowIdx < thisRowList.size(),
                    "Missing row interaction" );

      relativeMap.push_back( thisDataRowList[ thisRowIdx ] );
    }
  }
}

//////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////
void Supernode::buildLowRankRelativeMap(
                              const Supernode &descendent,
                              const PairArray &lowRankIndices,
                              const IntArray &baseRows,
                              IntArray &relativeMap )
{
  int                            thisRowIdx;
  int                            baseRow;

  relativeMap.clear();

  // Iterate over each interaction in the intersection and figure
  // out how rows in descendent map to rows in this
  for ( int interaction_idx = 0; interaction_idx < lowRankIndices.size();
        interaction_idx++ )
  {
    const IndexPair             &p = lowRankIndices[ interaction_idx ];

    const SupernodeInteraction  &interaction
                                  = descendent._offDiagonal[ p.first ];

    const IntArray              &rowList = interaction._rowList;

    const SupernodeInteraction  &thisInteraction = _offDiagonal[ p.second ];

    const IntArray              &thisRowList = thisInteraction._rowList;

    baseRow = baseRows[ p.second ];

    thisRowIdx = 0;

    TRACE_ASSERT( !thisInteraction._active,
                  "Interaction is not low-rank" );

    for ( int row_idx = 0; row_idx < rowList.size(); row_idx++ )
    {
      // Find the corresponding row in the interaction from
      // this node (must exist)
      while ( thisRowList[ thisRowIdx ] != rowList[ row_idx ]
           && thisRowIdx < thisRowList.size() )
      {
        thisRowIdx++;
      }

      TRACE_ASSERT( thisRowIdx < thisRowList.size(),
                    "Missing row interaction" );

      relativeMap.push_back( baseRow + thisRowIdx );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Copies data relative to the interaction between this and
// descendent in to the provided workspace so that we can do a
// single dense multiplication
//////////////////////////////////////////////////////////////////////
void Supernode::copyInteractionData( const Supernode &descendent,
                                     int start_idx,
                                     const PairArray &indices,
                                     Real *workspace,
                                     int &workCols,
                                     int &workRows1, int &workRows2 ) const
{
  int                            nCols = numColumns();
  int                            offset = 0;
  const Real                    *descData = _data;

  workCols = nCols;

  // Start by copying all rows from the interaction between
  // descendent and this
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ start_idx ];

    if ( interaction._type == STANDARD_NODE && interaction._active ) {
      memcpy( (void *)( workspace + offset ),
              (const void *)( descData + interaction._dataRowList[0] * nCols ),
              sizeof( Real ) * interaction._rowList.size() * nCols );

      offset += interaction._rowList.size() * nCols;

      workRows1 = interaction._rowList.size();
    } else {
      workRows1 = 0;
    }
  }

  workRows2 = 0;

  for ( int interaction_idx = 0; interaction_idx < indices.size();
        interaction_idx++ )
  {
    const IndexPair             &p = indices[ interaction_idx ];

    const SupernodeInteraction  &interaction
                                  = _offDiagonal[ p.first ];

    TRACE_ASSERT( interaction._type == STANDARD_NODE && interaction._active );

    memcpy( (void *)( workspace + offset ),
            (const void *)( descData + interaction._dataRowList[ 0 ] * nCols ),
            sizeof( Real ) * interaction._rowList.size() * nCols );

    offset += interaction._rowList.size() * nCols;

    workRows2 += interaction._rowList.size();
  }
}

static bool doneWrite = false;

//////////////////////////////////////////////////////////////////////
// Given a matrix of the form [ L1; L2 ], computes the update matrix
// [ L1; L2 ] * L1' and puts it in update
//
// C should have size at least ( workRows1 + workRows2 ) * workRows1
//////////////////////////////////////////////////////////////////////
void Supernode::buildUpdateMatrix( const Real *workspace,
                                   int workCols, int workRows1, int workRows2,
                                   const IntArray &relativeMap,
                                   Real *update,
                                   int descendent_idx ) const
{
  if ( workRows1 == 0 ) {
    // Nothing to do here
    return;
  }

  if ( !lowRankDiagonal() )
  {
    // Start with the symmetric multiply C := L1 * L1'
    MATRIX::syrk( workspace, update, workRows1, workCols );

#if 0
    // FIXME
    if ( numColumns() >= 1400 && !doneWrite )
    {
      char buf[ 1024 ];
      sprintf( buf, "super_numeric/node_%d_descendent_%d_update.matrix",
              _nodeID, descendent_idx );
      MATRIX tmp( workRows1, workRows1, update );
      tmp.write( buf );

      doneWrite = true;
    }
#endif

    // Next, append C := [ C; L2 * L1' ] using a standard matrix-matrix mult
    MATRIX::gemm( workspace + workRows1 * workCols, workspace,
                  update + workRows1 * workRows1,
                  workRows2, workCols, workRows1, workCols,
                  false, /* do not transpose L2 */
                  true /* transpose L1 */ );
  }
  else
  {
    // Just build the off diagonal part
    MATRIX::gemm( workspace + workRows1 * workCols, workspace,
                  update,
                  workRows2, workCols, workRows1, workCols,
                  false, /* do not transpose L2 */
                  true /* transpose L1 */ );

    // Put the diagonal update information at the end
    update += workRows2 * workRows1;

    // Find the parts of the matrix that need to be multiplied
    // together to form diagonal blocks
    int                     block_idx = 0;
    int                     start_row_idx;
    int                     dataRowIdx;
    int                     nRows;
    int                     nCols = numColumns();

    for ( int row_idx = 0; row_idx < workRows1; )
    {
      dataRowIdx = relativeMap[ row_idx ];
      TRACE_ASSERT( dataRowIdx < nCols, "Should never be here" );

      while ( block_idx < _diagonalBlocks.size()
           && !_diagonalBlocks[ block_idx ].containsRow( dataRowIdx ) )
      {
        block_idx++;
      }

      TRACE_ASSERT( block_idx < _diagonalBlocks.size(),
                    "No diagonal block found" );

      start_row_idx = row_idx;

      // Advance to find the end row
      while ( row_idx < workRows1
        && _diagonalBlocks[ block_idx ].containsRow( relativeMap[ row_idx ] ) )
#if 0
      while ( _diagonalBlocks[ block_idx ].containsRow( dataRowIdx )
           && row_idx < workRows1 )
#endif
      {
        row_idx++;
#if 0
        dataRowIdx = relativeMap[ row_idx ];
#endif
      }

      // We now have a block range
      nRows = row_idx - start_row_idx;

      MATRIX::syrk( workspace + start_row_idx * workCols,
                    update, nRows, workCols );

      update += nRows * nRows;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Let L1 be the interaction in this indexed by ancestor_idx.
// Let L2 be the interaction in this indexed by ancestor_interaction_idx.
// This function forms a block column of the product L2 * L1' and
// subtracts this from the given workspace matrix.
//
// Rows of L2 * L1' are scattered in to workspace according to the
// provided workspace row list.
//////////////////////////////////////////////////////////////////////
void Supernode::subtractInteractionUpdateMatrix( int ancestor_idx,
                                                 int ancestor_interaction_idx,
                                                 const IndexRange &columnRange,
                                                 const IntArray &workspaceRows,
                                                 Real *multWorkspace,
                                                 Real *workspace ) const
{
  const SupernodeInteraction &interactionL1 = _offDiagonal[ ancestor_idx ];
  const SupernodeInteraction &interactionL2
                                  = _offDiagonal[ ancestor_interaction_idx ];

  const Real                *dataL1;
  const Real                *dataL2;

  IndexRange                 columnRangeL1;
  int                        numColumnsL1;

  int                        nCols = numColumns();
  int                        nColsRange = range_size( columnRange );
  
  const IntArray            &rowsL1 = interactionL1._rowList;
  const IntArray            &rowsL2 = interactionL2._rowList;

  int                        nRowsL2 = rowsL2.size();

  dataL1 = interactionMatrix( ancestor_idx );
  dataL2 = interactionMatrix( ancestor_interaction_idx );

  // Figure out how many rows (columns after transposing) that we need
  // to form the product for the given column subset
  findEntryRangeIntersection( rowsL1, columnRange, columnRangeL1 );

  TRACE_ASSERT(
       (  ( columnRangeL1.first == 0 && rowsL1[ 0 ] >= columnRange.first )
       || ( rowsL1[ columnRangeL1.first ] >= columnRange.first
         && rowsL1[ columnRangeL1.first - 1 ] < columnRange.first ) )
    && (  ( columnRangeL1.second == rowsL1.size() - 1
         && rowsL1[ columnRangeL1.second ] <= columnRange.second )
       || ( rowsL1[ columnRangeL1.second ] <= columnRange.second
         && rowsL1[ columnRangeL1.second + 1 ] > columnRange.second ) ),
    "Invalid column range" );

  if ( columnRangeL1.first == -1 )
  {
    return;
  }

  numColumnsL1 = range_size( columnRangeL1 );

  // Advance to the portion of L1 we are interested in
  dataL1 += columnRangeL1.first * nCols;

  TRACE_ASSERT( numColumnsL1 > 0, "No columns to work with" );

#if 0
  TRACE_ASSERT(
    columnRangeL1.first == 0 && columnRangeL1.second == rowsL1.size() - 1,
    "Invalid column range" );
#endif

  // FIXME: debugging
  if ( nRowsL2 * numColumnsL1 > multWorkspaceSz )
  {
    TRACE_ASSERT( NULL, "Workspace too small" );
  }

  // Multiply L2 * L1( columnRangeL1, : )'
  MATRIX::gemm( dataL2, dataL1, multWorkspace,
                nRowsL2, nCols, /* Dimensions of L2 */
                numColumnsL1, nCols, /* Dimensions of L1 subset */
                false, true /* Transpose L1 but not L2 */ );

  // Scatter back to the desired row/column space in workspace,
  // according to workspaceRows
  int                        full_row_idx = 0;
  int                        full_col_idx;

  for ( int row_idx = 0; row_idx < rowsL2.size(); row_idx++ )
  {
    // Find the corresponding row in the output matrix
    while ( full_row_idx < workspaceRows.size()
         && workspaceRows[ full_row_idx ] != rowsL2[ row_idx ] )
    {
      full_row_idx++;
    }

    TRACE_ASSERT( full_row_idx < workspaceRows.size(),
                  "No matching row index found" );

    for ( int col_idx = columnRangeL1.first; col_idx <= columnRangeL1.second;
          col_idx++ )
    {
      full_col_idx = rowsL1[ col_idx ] - columnRange.first;
      
      // Extract the desired entry
      MATRIX::access( workspace,
                      workspaceRows.size(), nColsRange, /* Output dimensions */
                      full_row_idx, full_col_idx )
        -= MATRIX::access( multWorkspace,
                           rowsL2.size(), numColumnsL1, /* Input dimensions */
                           row_idx, col_idx - columnRangeL1.first );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Subtracts the content of the given update matrix from this
// node's data according to the given relative map
//
// workRows1, workRows2 are the same variables used by
// copyInteractionData and buildUpdateMatrix
//////////////////////////////////////////////////////////////////////
void Supernode::applyUpdate( const IntArray &relativeMap,
                             int workRows1, int workRows2,
                             const Real *update )
{
  int                            dataRowIdx;
  int                            dataColIdx;
  int                            nCols = numColumns();
  int                            nRows;
  int                            workRows;
  int                            modifier;
  int                            startRow;

  Real                          *baseData;

  if ( workRows1 == 0 ) {
    // Nothing to do here
    return;
  }

  baseData = lowRankDiagonal() ? _data + _firstOffDiagonalRow * nCols : _data;
  nRows = lowRankDiagonal() ? _numRows : totalRows();

  // This compensates for the fact that we no longer store the diagonal
  // update (see below)
#if 0
  modifier = lowRankDiagonal() ? nCols : 0;
#endif
  modifier = lowRankDiagonal() ? _firstOffDiagonalRow : 0;
  startRow = lowRankDiagonal() ? workRows1 : 0;

  // We don't form the diagonal part if we are sparsifying it.
  workRows = lowRankDiagonal() ? workRows2 : workRows1 + workRows2;

#if 0
  if ( lowRankDiagonal() )
  {
    TRACE_ASSERT( relativeMap.size() == workRows2,
                  "Mismatch between relative map and update matrix" );
  }
  else
  {
#endif
    TRACE_ASSERT( relativeMap.size() == workRows1 + workRows2,
                  "Mismatch between relative map and update matrix" );
#if 0
  }
#endif

  for ( int row_idx = startRow; row_idx < relativeMap.size(); row_idx++ )
  {
    dataRowIdx = relativeMap[ row_idx ];

    // This compensates for the fact that we no longer store the diagonal
    // update
    dataRowIdx -= modifier;
    if ( dataRowIdx < 0 )
    {
      TRACE_ASSERT( NULL, "Should never get here" );
      continue;
    }

#if 0
    // FIXME
    if ( lowRankDiagonal() )
    {
      TRACE_ASSERT( dataRowIdx >= _firstOffDiagonalRow, "Something's wrong" );
    }
#endif

    for ( int col_idx = 0; col_idx < workRows1; col_idx++ )
    {
      dataColIdx = relativeMap[ col_idx ];

      MATRIX::access( baseData, nRows, nCols, dataRowIdx, dataColIdx )
        -= MATRIX::access( update, workRows, workRows1,
                           row_idx - startRow, col_idx );
    }
  }

  if ( lowRankDiagonal() )
  {
    // Add the diagonal block stuff.  This follows the same
    // structure as buildUpdateMatrix
    update += workRows2 * workRows1;

    baseData = _data;

    // Find the parts of the matrix that need to be multiplied
    // together to form diagonal blocks
    int                     block_idx = 0;
    int                     start_row_idx;
    int                     dataRowIdx;
    int                     nRowsFull;

    for ( int row_idx = 0; row_idx < workRows1; )
    {
      dataRowIdx = relativeMap[ row_idx ];
      TRACE_ASSERT( dataRowIdx < nCols, "Should never be here" );

      while ( block_idx < _diagonalBlocks.size()
           && !_diagonalBlocks[ block_idx ].containsRow( dataRowIdx ) )
      {
        nRowsFull = _diagonalBlocks[ block_idx ].numRows();
        baseData += nRowsFull * nRowsFull;

        block_idx++;
      }

      TRACE_ASSERT( block_idx < _diagonalBlocks.size(),
                    "No diagonal block found" );

      start_row_idx = row_idx;

      // Advance to find the end row
      while ( row_idx < workRows1
        && _diagonalBlocks[ block_idx ].containsRow( relativeMap[ row_idx ] ) )
#if 0
      while ( _diagonalBlocks[ block_idx ].containsRow( dataRowIdx )
           && row_idx < workRows1 )
#endif
      {
        row_idx++;
#if 0
        dataRowIdx = relativeMap[ row_idx ];
#endif
      }

      // We now have a block range
      nRows = row_idx - start_row_idx;
      nRowsFull = _diagonalBlocks[ block_idx ].numRows();

      // Fill in the block
      for ( int i = start_row_idx; i < row_idx; i++ )
      {
        // Get the row index and align with the block starting point
        dataRowIdx = relativeMap[ i ];
        dataRowIdx -= _diagonalBlocks[ block_idx ]._rowRange.first;

        TRACE_ASSERT( dataRowIdx >= 0 && dataRowIdx < nRowsFull,
                      "Row index out of range" );

        TRACE_ASSERT( i - start_row_idx >= 0 && i - start_row_idx < nRows,
                      "Input row index out of range" );

        for ( int j = start_row_idx; j < row_idx; j++ )
        {
          // Get the column index and align with the block starting point
          dataColIdx = relativeMap[ j ];
          dataColIdx -= _diagonalBlocks[ block_idx ]._columnRange.first;

          TRACE_ASSERT( dataColIdx >= 0 && dataColIdx < nRowsFull,
                        "Column index out of range" );

          TRACE_ASSERT( j - start_row_idx >= 0 && j - start_row_idx < nRows,
                        "Input column index out of range" );

          MATRIX::access(
            baseData, nRowsFull, nRowsFull, dataRowIdx, dataColIdx )
              -= MATRIX::access( update, nRows, nRows,
                                 i - start_row_idx,
                                 j - start_row_idx );
        }
      }

      update += nRows * nRows;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Subtracts the outer product of the given matrix from this
// node's diagonal.
//
// update's number of rows should match this node's column count
//////////////////////////////////////////////////////////////////////
void Supernode::applyDiagonalUpdate( const Real *update, int nCols )
{
  const Real                *baseInputData = update;
  Real                      *baseData = _data;

  if ( lowRankDiagonal() )
  {
    for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ ) {
      int                    blockRows = _diagonalBlocks[ block_idx ].numRows();

      MATRIX::syrk( baseInputData, baseData,
                    blockRows, nCols,
                    false /* don't transpose */,
                    -1.0, 1.0 /* subtract */ );

      baseInputData += blockRows * nCols;
      baseData += blockRows * blockRows;
    }
  } else {
    // Just in case _diagonalBlocks.size() == 0 (in which case the above
    // branch won't work
    MATRIX::syrk( baseInputData, baseData,
                  numColumns(), nCols,
                  false /* don't transpose */,
                  -1.0, 1.0 /* subtract */ );
  }
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but not responsible for any multiplication.
// Also, uses a row list
//////////////////////////////////////////////////////////////////////
void Supernode::applyDiagonalUpdate( const Real *update,
                                     Real *multWorkspace, int nCols,
                                     const IntArray &rowList )
{
  const Real                *inputData = update;

  if ( lowRankDiagonal() ) {
    int                      i = 0;

    while ( i < rowList.size() ) {
      int                    start_row_idx = i;
      int                    row_idx = rowList[ start_row_idx ];
      int                    block_idx = _diagonalMap[ row_idx ];
      int                    nRows;
      int                    nBlockRows;

      int                    blockRow;
      Real                  *blockData;
      
      blockRow = _diagonalBlocks[ block_idx ]._rowRange.first;
      blockData = _data + _diagonalOffsets[ block_idx ];
      nBlockRows = _diagonalBlocks[ block_idx ].numRows();

      while ( i < rowList.size() && _diagonalMap[ row_idx ] == block_idx ) {
        i++;
        row_idx = rowList[ i ];
      }

      inputData = update + start_row_idx * nCols;

      nRows = i - start_row_idx;

      // Overwrite multWorkspace with the desired symmetric update
      MATRIX::syrk( inputData, multWorkspace,
                    nRows, nCols );

      // Subtract this update
      for ( int j = start_row_idx; j < i; j++ )
      for ( int k = start_row_idx; k <= j; k++ ) {
        MATRIX::access( blockData, nBlockRows, nBlockRows,
                        rowList[ j ] - blockRow, rowList[ k ] - blockRow )
          -= MATRIX::access( multWorkspace, nRows, nRows,
                             j - start_row_idx, k - start_row_idx );
      }
    }
  } else {
    int                      nRows = rowList.size();
    int                      nBlockRows = numColumns();

    Real                    *blockData = _data;

    MATRIX::syrk( update, multWorkspace, nRows, nCols );

    for ( int j = 0; j < nRows; j++ )
    for ( int k = 0; k < nRows; k++ ) {
      MATRIX::access( blockData, nBlockRows, nBlockRows,
                      rowList[ j ], rowList[ k ] )
        -= MATRIX::access( multWorkspace, nRows, nRows, j, k );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// This version forms a diagonal update for a provided diagonal block
// in this node's main diagonal.
//
// Forms a diagonal update and puts it in to the provided schur
// complement matrix.
//////////////////////////////////////////////////////////////////////
void Supernode::applyDiagonalUpdate( const DenseBlock &factorBlock,
                                     const Real *update,
                                     Real *multWorkspace, int nCols,
                                     const IntArray &rowList,
                                     Real *schurComplement,
                                     // Optional range for rowList
                                     IndexRange rowRange ) const
{
  int                        nRows;

  // Must be a diagonal block
  TRACE_ASSERT( factorBlock.isDiagonal() );

  if ( rowRange.first < 0 || rowRange.second < rowRange.first
        || rowRange.second >= rowList.size() ) {
    rowRange.first = 0;
    rowRange.second = rowList.size() - 1;
  }

  nRows = range_size( rowRange );

  // Build the update matrix
  MATRIX::syrk( update, multWorkspace, nRows, nCols );

  // Subtract from the matrix
  for ( int j = rowRange.first; j <= rowRange.second; j++ )
  for ( int k = rowRange.first; k <= rowRange.second; k++ ) {
    MATRIX::access( schurComplement,
                    factorBlock.numRows(), factorBlock.numColumns(),
                    // Indices in to this matrix - we need to make
                    // them relative to the block we are trying to form
                    rowList[ j ] - factorBlock._rowRange.first,
                    rowList[ k ] - factorBlock._columnRange.first )
      -= MATRIX::access( multWorkspace, nRows, nRows,
                         j - rowRange.first, k - rowRange.first );
  }
}

//////////////////////////////////////////////////////////////////////
// Update this node's diagonal with the contribution from a
// descendent, but only if the relevent interaction in that descendent
// is stored in compressed form
//////////////////////////////////////////////////////////////////////
void Supernode::compressedDiagonalUpdate( const Supernode &descendent,
                                          int ancestor_idx )
{
  const SupernodeInteraction  &interaction
                                  = descendent.offDiagonal()[ ancestor_idx ];

  const IntArray              &rowList = interaction._rowList;

  int                          block_idx = 0;
  int                          start_row_idx = 0;
  int                          row_idx;

  int                          nCols = numColumns();
  int                          descendentCols = descendent.numColumns();

  int                          nRows;
  int                          nRowsFull;
  int                          blockStartRow;

  const Real                  *V = interaction._compressedV;
  const Real                  *U = interaction._compressedU;
  const Real                  *Vsub;

  Real                        *blockData;

  TRACE_ASSERT( interaction._nodeID == _nodeID );

  if ( interaction._type != COMPRESSED_NODE ) {
    return;
  }

  // This will store the product U' * U, since we will need this
  // numerous times during this computation
  MATRIX                       innerProduct( interaction._numExtendedRows,
                                             interaction._numExtendedRows );

  MATRIX::gemm( U, U, innerProduct.data(),
                // Size of U
                descendentCols, interaction._numExtendedRows,
                descendentCols, interaction._numExtendedRows,
                // Transpose, no transpose
                true, false );

  if ( lowRankDiagonal() )
  {
    // Update all blocks affected by this interaction
    while ( start_row_idx < rowList.size() )
    {
      block_idx = _diagonalMap.at( rowList[ start_row_idx ] );
      row_idx = start_row_idx;

      const DenseBlock          &block = _diagonalBlocks[ block_idx ];

      blockStartRow = block._rowRange.first;

      // Determine the row range from this interaction that affects
      // this block
      while ( row_idx < rowList.size()
           //&& block.containsRow( rowList[ row_idx ] - blockStartRow ) )
           && block.containsRow( rowList[ row_idx ] ) )
      {
        row_idx += 1;
      }

      nRows = row_idx - start_row_idx;
      nRowsFull = block.numRows();

      TRACE_ASSERT( nRows > 0 && nRowsFull > 0 );

      // Get the piece of V that we are interested in
      Vsub = V + start_row_idx * interaction._numExtendedRows;

      // Workspaces for forming the two matrix products we need here
      MATRIX                     workspace1( interaction._numExtendedRows,
                                             nRows );
      MATRIX                     workspace2( nRows, nRows );

      MATRIX::gemm( innerProduct.data(), Vsub, workspace1.data(),
                    // Size of the inner product matrix
                    interaction._numExtendedRows, interaction._numExtendedRows,
                    // Size of Vsub
                    nRows, interaction._numExtendedRows,
                    // Transpose Vsub only
                    false, true );

      MATRIX::gemm( Vsub, workspace1.data(), workspace2.data(),
                    // Size of Vsub
                    nRows, interaction._numExtendedRows,
                    // Size of the workspace
                    interaction._numExtendedRows, nRows,
                    // No transposition
                    false, false );

      // Scatter/subtract this matrix to the desired diagonal block
      blockData = _data + _diagonalOffsets[ block_idx ];

      for ( int i = 0; i < nRows; i++ )
      for ( int j = 0; j < nRows; j++ )
      {
        MATRIX::access( blockData, nRowsFull, nRowsFull,
                        rowList[ start_row_idx + i ] - blockStartRow,
                        rowList[ start_row_idx + j ] - blockStartRow )
          -= MATRIX::access( workspace2.data(), nRows, nRows, i, j );
      }

      start_row_idx = row_idx;
    }
  }
  else
  {
    // Single block update
    nRows = interaction._rowList.size();
    nRowsFull = nCols;

    // Workspaces for forming the two matrix products we need here
    MATRIX                     workspace1( interaction._numExtendedRows,
                                           nRows );
    MATRIX                     workspace2( nRows, nRows );

    MATRIX::gemm( innerProduct.data(), V, workspace1.data(),
                  // Size of the inner product matrix
                  interaction._numExtendedRows, interaction._numExtendedRows,
                  // Size of V
                  nRows, interaction._numExtendedRows,
                  // Transpose V
                  false, true );

    MATRIX::gemm( V, workspace1.data(), workspace2.data(),
                  // Size of V
                  nRows, interaction._numExtendedRows,
                  // Size of the workspace
                  interaction._numExtendedRows, nRows,
                  // No transposition
                  false, false );

    // Scatter/subtract this matrix to the desired diagonal block
    blockData = _data;

    for ( int i = 0; i < nRows; i++ )
    for ( int j = 0; j < nRows; j++ )
    {
      MATRIX::access( blockData, nRowsFull, nRowsFull,
                      rowList[ i ], rowList[ j ] )
        -= MATRIX::access( workspace2.data(), nRows, nRows, i, j );
    }

    start_row_idx = row_idx;
  }
}

//////////////////////////////////////////////////////////////////////
// A compressed diagonal update which explicitly subtracts its result
// in to the given schur complement matrix.
//
// This is similar to applyDiagonalUpdate above, but update has
// been replaced by V * U' = update
//
// multWorkspace should be of size at least:
//    rank * ( rank + 2 * range_size( rowRange ) )
//  if rowRange is provided or:
//    rank * ( rank + 2 * rowList.size() )
//  if not
//////////////////////////////////////////////////////////////////////
void Supernode::compressedDiagonalUpdate( const DenseBlock &factorBlock,
                                          const Real *V, const Real *U, int rank,
                                          Real *multWorkspace, int nCols,
                                          const IntArray &rowList,
                                          Real *schurComplement,
                                          // Optional range for rowList
                                          IndexRange rowRange ) const
{
  int                        nRows;
  Real                      *workspace1;
  Real                      *workspace2;
  Real                      *workspace3;

  // Must be a diagonal block
  TRACE_ASSERT( factorBlock.isDiagonal() );

  if ( rowRange.first < 0 || rowRange.second < rowRange.first
        || rowRange.second >= rowList.size() ) {
    rowRange.first = 0;
    rowRange.second = rowList.size() - 1;
  }

  nRows = range_size( rowRange );

  // Intermediate matrices needed during this calculation
  //    First one of size ( rank * rank )
  workspace1 = multWorkspace;

  //    Now one of size ( rank * nRows )
  workspace2 = workspace1 + rank * rank;

  //    Now one of size ( rank * nRows )
  workspace3 = workspace2 + rank * nRows;

  // Skip to the part of V that we care about
  V += rowRange.first * rank;

  // The update we want has the form ( V * U' * U * V' ), so presumably the
  // most efficient way of forming this is to start by forming U' * U, since
  // this should be small
  //
  // We don't use syrk here, since we need the whole matrix
  MATRIX::gemm( U, U, workspace1,
                // Dimensions of U (twice)
                nCols, rank, nCols, rank,
                // Transpose the first one
                true, false );

  // Next, form (U' * U) * V'
  MATRIX::gemm( workspace1, V, workspace2,
                // Input dimensions
                rank, rank, nRows, rank,
                // Transpose the second input
                false, true );

  // Form V * ((U' * U) * V')  (no transposition needed here)
  MATRIX::gemm( V, workspace2, workspace3,
                // Input dimensions
                nRows, rank, rank, nRows,
                // No transposition
                false, false );

  // Scatter the update to the schur complement block
  for ( int row_idx = rowRange.first; row_idx <= rowRange.second; row_idx++ )
  for ( int col_idx = rowRange.first; col_idx <= rowRange.second; col_idx++ ) {
    MATRIX::access( schurComplement,
                    factorBlock.numRows(), factorBlock.numColumns(),
                    rowList[ row_idx ] - factorBlock._rowRange.first,
                    rowList[ col_idx ] - factorBlock._columnRange.first )
      -= MATRIX::access( workspace3, nRows, nRows,
                         row_idx - rowRange.first, col_idx - rowRange.first );
  }
}

//////////////////////////////////////////////////////////////////////
// Forms update matrices for each extended ndoe interaction between
// the given descendent node and this.  Subtract these from the
// relevant extended interactions in this node.
//
// We require a workspace for multiplication, and the extra data
// array in which extended interactions are stored.
//
// start_idx is the same input used to buildRelativeMap, etc.
//////////////////////////////////////////////////////////////////////
void Supernode::applyExtendedUpdate( const Supernode &descendent,
                                     const PairArray &extendedIndices,
                                     Real *multWorkspace,
                                     Real *extraData,
                                     int start_idx,
                                     Real *expansionWorkspace,
                                     Timer *standardTimer,
                                     Timer *extendedTimer )
{
  if ( _type == STANDARD_NODE )
  {
#ifdef DO_TIMING
    if ( standardTimer )
    {
      standardTimer->tick();
    }
#endif
    applyExtendedUpdate_standard( descendent, extendedIndices,
                                  multWorkspace, extraData, start_idx,
                                  expansionWorkspace );
#ifdef DO_TIMING
    if ( standardTimer )
    {
      standardTimer->tock();
    }
#endif
  }
  else
  {
#ifdef DO_TIMING
    if ( extendedTimer ) {
      extendedTimer->tick();
    }
#endif
    applyExtendedUpdate_extended( descendent, extendedIndices,
                                  extraData, start_idx, NULL,
                                  standardTimer, extendedTimer );
#ifdef DO_TIMING
    if ( extendedTimer ) {
      extendedTimer->tock();
    }
#endif
  }
}

//////////////////////////////////////////////////////////////////////
// Forms update matrices for each extended ndoe interaction between
// the given descendent node and this.  Subtract these from the
// relevant extended interactions in this node.
//
// We require a workspace for multiplication, and the extra data
// array in which extended interactions are stored.
//
// start_idx is the same input used to buildRelativeMap, etc.
//////////////////////////////////////////////////////////////////////
void Supernode::applyExtendedUpdate( const Supernode &descendent,
                                     const PairArray &extendedIndices,
                                     Real *extraData,
                                     int start_idx,
                                     WorkspaceManager<Real> &workspaceManager,
                                     Timer *standardTimer,
                                     Timer *extendedTimer )
{
  if ( _type == STANDARD_NODE )
  {
#ifdef DO_TIMING
    if ( standardTimer ) {
      standardTimer->tick();
    }
#endif
    applyExtendedUpdate_standard( descendent, extendedIndices,
                                  extraData, start_idx,
                                  workspaceManager );
#ifdef DO_TIMING
    if ( standardTimer ) {
      standardTimer->tock();
    }
#endif
  }
  else
  {
    if ( _useLDL ) {
#ifdef DO_TIMING
      if ( extendedTimer ) {
        extendedTimer->tick();
      }
#endif
      applyExtendedUpdate_extendedLDL( descendent, extendedIndices,
                                       extraData, start_idx, workspaceManager,
                                       standardTimer, extendedTimer );
#ifdef DO_TIMING
      if ( extendedTimer ) {
        extendedTimer->tock();
      }
#endif
    }
    else {
#ifdef DO_TIMING
      if ( extendedTimer ) {
        extendedTimer->tick();
      }
#endif
      applyExtendedUpdate_extended( descendent, extendedIndices,
                                    extraData, start_idx, NULL,
                                    standardTimer, extendedTimer );
#ifdef DO_TIMING
      if ( extendedTimer ) {
        extendedTimer->tock();
      }
#endif
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Returns a list of interactions to be sparsified in this node
//////////////////////////////////////////////////////////////////////
void Supernode::getLowRankBlocks( IntArray &lowRankBlocks,
                                  IntArray &baseRows ) const
{
  int                            baseRow = 0;

  lowRankBlocks.clear();
  baseRows.clear();

  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    if ( (interaction._type == STANDARD_NODE && !interaction._active)
      || interaction._type == COMPRESSED_NODE )
    {
      lowRankBlocks.push_back( interaction_idx );
      baseRows.push_back( baseRow );

      baseRow += interaction._rowList.size();
    }
  }
}

//////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////
void Supernode::appendLowRankDescendentList(
                  const Supernode &descendent,
                  const PairArray &lowRankIndices,
                  const IntArray &lowRankBlocks,
                  vector<PairArray > &lowRankDescendents )
{
  int                            block_idx;

  block_idx = 0;

  for ( int interaction_idx = 0; interaction_idx < lowRankIndices.size();
        interaction_idx++ )
  {
    const IndexPair             &p = lowRankIndices[ interaction_idx ];

    // Find the corresponding entry in the full list of
    // low rank blocks for this node.
    while ( lowRankBlocks[ block_idx ] != p.second
         && block_idx < lowRankBlocks.size() )
    {
      block_idx++;
    }

    TRACE_ASSERT( block_idx < lowRankBlocks.size(),
                  "No matching interaction found" );
    
    lowRankDescendents[ block_idx ].push_back( IndexPair( descendent._nodeID,
                                                          p.first ) );
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void Supernode::clearLowRankDescendentList(
                  int numEntries,
                  vector<PairArray> &lowRankDescendents )
{
  for ( int i = 0; i < numEntries; i++ )
  {
    lowRankDescendents[ i ].clear();
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Given a set of interactions (U1, U2, ...) from *this* node,
// and the interaction U_a with an ancestor node, form the
// matrix product [ U1; U2; ... ] * U_a' * G (where G is a given
// matrix) and accumulate it in to the given workspace.
//////////////////////////////////////////////////////////////////////
void Supernode::interactionMultiply( const PairArray &lowRankIndices,
                                     int ancestor_idx,
                                     const Supernode &ancestor,
                                     const Real *G, int nColsG,
                                     Real *multWorkspaceInitial,
                                     Real *multWorkspaceFinal )
{
  Real                          *baseData;
  Real                          *baseOutputData;
  int                            nCols = numColumns();
  int                            nRows;

  TRACE_ASSERT( _offDiagonal[ ancestor_idx ]._nodeID == ancestor._nodeID,
                "Node ID mismatch" );

  // Start by forming the product
  //    multWorkspaceInitial = U_a' * G
  baseData = _data + _offDiagonal[ ancestor_idx ]._dataRowList[ 0 ] * nCols;
  nRows = _offDiagonal[ ancestor_idx ]._rowList.size();

  MATRIX::gemm( baseData, G, multWorkspaceInitial,
                nRows, nCols, nRows, nColsG,
                true, /* Transpose U_a */
                false /* Don't transpose G */ );

  baseOutputData = multWorkspaceFinal;

  for ( int interaction_idx = 0; interaction_idx < lowRankIndices.size();
        interaction_idx++ )
  {
    const IndexPair             &p = lowRankIndices[ interaction_idx ];

    const SupernodeInteraction  &interaction = _offDiagonal[ p.first ];

    baseData = _data + interaction._dataRowList[ 0 ] * nCols;
    nRows = interaction._rowList.size();

    MATRIX::gemm( baseData, multWorkspaceInitial, baseOutputData,
                  nRows, nCols, nCols, nColsG,
                  false, false /* No transposition */ );

    baseOutputData += nRows * nColsG;
  }
}
#endif

// FIXME: debugging for compressed interaction multiplies
#define DEBUG_COMPRESSED_MULT 1
#undef DEBUG_COMPRESSED_MULT

//////////////////////////////////////////////////////////////////////
// Given the interactions below ancestor_idx in this node (call these
// (U1, U2, ...) and the interaction U_a with the ancestor node,
// form the matrix product [ U1 U2 ... ]' * U_a' * G (where G is a
// given matrix) and place the result in the given matrix (doesn't
// do any row scattering so the output will have to be processed)
//////////////////////////////////////////////////////////////////////
void Supernode::interactionMultiply( int ancestor_idx,
                                     const Supernode &ancestor,
                                     const Real *G, int nColsG,
                                     Real *multWorkspaceInitial,
                                     Real *multWorkspaceFinal )
{
  int                              nCols = numColumns();

  if ( _offDiagonal[ ancestor_idx ]._type == STANDARD_NODE )
  {
    Real                    *baseData;
    Real                    *baseOutputData;
    int                      nRows = 0;

    // Start by multiplying U_a' * G
    baseData = _data + _offDiagonal[ ancestor_idx ]._dataRowList[ 0 ] * nCols;
    nRows = _offDiagonal[ ancestor_idx ]._rowList.size();

    MATRIX::gemm( baseData, G, multWorkspaceInitial,
                  nRows, nCols, nRows, nColsG,
                  true, /* Transpose U_a */
                  false /* Don't transpose G */ );

    baseOutputData = multWorkspaceFinal;

    // Count rows from interactions below ancestor_idx in this node
    nRows = countInteractionRows( ancestor_idx + 1, _offDiagonal.size() - 1 );

    // Move to the part of _compressedV that we want to multiply with
    baseData
      = _data + _offDiagonal[ ancestor_idx + 1 ]._dataRowList[ 0 ] * nCols;

    MATRIX::gemm( baseData, multWorkspaceInitial, baseOutputData,
                  nRows, nCols, nCols, nColsG,
                  false, false /* No transposition */ );
  }
  else if ( _offDiagonal[ ancestor_idx ]._type == COMPRESSED_NODE )
  {
    // The ancestor interaction
    const SupernodeInteraction  &interaction1 = _offDiagonal[ ancestor_idx ];

    // Transposed multiplication by the ancestor interaction
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( ancestor_idx, G,
                                   // Size of G (it has already been
                                   // compressed)
                                   interaction1._rowList.size(), nColsG,
                                   multWorkspaceInitial,
                                   // Left, transposed multiplication
                                   true, true,
                                   // No need to extract a submatrix of
                                   // G (should already have been done)
                                   false );

    // Untransposed multiplication with the off-diagonal block below
    // the ancestor interaction
    compressedInteractionMultiply_fullBlock( ancestor_idx, multWorkspaceInitial,
                                             nCols, nColsG,
                                             multWorkspaceFinal,
                                             // Left, untransposed
                                             true, false,
                                             // No need to extract a submatrix
                                             // of G (shouldn't matter for this
                                             // type of multiplication anyways)
                                             false );
  }
  else
  {
    TRACE_ASSERT( NULL, "Invalid interaction type" );
  }
}

//////////////////////////////////////////////////////////////////////
// Row counting-routine used by the above function
//////////////////////////////////////////////////////////////////////
int Supernode::countInteractionRows( int startInteraction, int endInteraction )
{
  int                        nRows = 0;

  for ( int interaction_idx = startInteraction;
        interaction_idx <= endInteraction; interaction_idx++ )
  {
    nRows += _offDiagonal[ interaction_idx ]._rowList.size();
  }

  return nRows;
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but only handles a single interaction
//////////////////////////////////////////////////////////////////////
void Supernode::interactionMultiply( int interaction_idx,
                                     int ancestor_idx,
                                     const Supernode &ancestor,
                                     const Real *G, int nColsG,
                                     Real *multWorkspaceInitial,
                                     Real *multWorkspaceFinal )
{
  int                              nCols = numColumns();

  if ( _offDiagonal[ ancestor_idx ]._type == STANDARD_NODE )
  {
    Real                          *baseData;
    Real                          *baseOutputData;
    int                            nRows;

    TRACE_ASSERT( _offDiagonal[ ancestor_idx ]._nodeID == ancestor._nodeID,
                  "Node ID mismatch" );

    // Start by forming the product
    //    multWorkspaceInitial = U_a' * G
    baseData = _data + _offDiagonal[ ancestor_idx ]._dataRowList[ 0 ] * nCols;
    nRows = _offDiagonal[ ancestor_idx ]._rowList.size();

    MATRIX::gemm( baseData, G, multWorkspaceInitial,
                  nRows, nCols, nRows, nColsG,
                  true, /* Transpose U_a */
                  false /* Don't transpose G */ );

    baseOutputData = multWorkspaceFinal;

    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    baseData = _data + interaction._dataRowList[ 0 ] * nCols;
    nRows = interaction._rowList.size();

    MATRIX::gemm( baseData, multWorkspaceInitial, baseOutputData,
                  nRows, nCols, nCols, nColsG,
                  false, false /* No transposition */ );
  }
  else if ( _offDiagonal[ ancestor_idx ]._type == COMPRESSED_NODE )
  {
    // The two interactions from this node involved in forming the
    // Schur complement product
    const SupernodeInteraction  &interaction1 = _offDiagonal[ ancestor_idx ];
    const SupernodeInteraction  &interaction2 = _offDiagonal[ interaction_idx ];

    // For now anyways, assume that they are both stored in compressed form
    //
    // TODO: Should this be fixed later?
    TRACE_ASSERT( interaction2._type == COMPRESSED_NODE );

    // Transposed multiplication by interaction 1
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( ancestor_idx, G,
                                   // Size of G (it has already been
                                   // compressed)
                                   interaction1._rowList.size(), nColsG,
                                   multWorkspaceInitial,
                                   // Left, transposed multiplication
                                   true, true,
                                   // No need to extract a submatrix of
                                   // G (should already have been done)
                                   false );

    // Untransposed multiplication by interaction 2
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( interaction_idx, multWorkspaceInitial,
                                   // Size of matrix formed in the last
                                   // multiplication
                                   nCols, nColsG,
                                   multWorkspaceFinal,
                                   // Left, untransposed multiplication
                                   true, false,
                                   // No need to extract a submatrix here
                                   // (shouldn't matter for this type of
                                   // multiplication anyways)
                                   false );
  }
  else
  {
    TRACE_ASSERT( NULL, "Invalid interaction type" );
  }
}

//////////////////////////////////////////////////////////////////////
// Assuming that U_a is the intersaction between this node and
// the given ancestor node, forms a subset of the matrix U_a * U_a'
// with the specified row and column range.  These ranges given
// in terms of the actual set of rows stored in U_a.
// That is, if rowRange = [ 2, 3 ] and this node stores rows
// r0, r1, r2, r3, r4, ..., then we are considering the row subset
// [ r2, r3 ] of U_a
//////////////////////////////////////////////////////////////////////
void Supernode::interactionMultiply( int ancestor_idx,
                                     const Supernode &ancestor,
                                     const DenseBlock &ancestorBlock,
                                     const Real *G, int nColsG,
                                     const IndexRange &rowRange,
                                     const IndexRange &columnRange,
                                     Real *multWorkspaceInitial,
                                     Real *multWorkspaceFinal )
{
  int                              nCols;
  int                              nCols_full = numColumns();
  int                              nRows;

  nRows = range_size( rowRange );
  nCols = range_size( columnRange );

  if ( _offDiagonal[ ancestor_idx ]._type == STANDARD_NODE )
  {
    Real                          *baseData;
    Real                          *baseOutputData;

    const SupernodeInteraction    &interaction = _offDiagonal[ ancestor_idx ];

    TRACE_ASSERT( _offDiagonal[ ancestor_idx ]._nodeID == ancestor._nodeID,
                  "Node ID mismatch" );

    // Start by forming the product
    //    multWorkspaceInitial = U_a' * G
    //
    // First we need to get the proper sub-matrix from this interaction's data
    baseData = _data + interaction._dataRowList[ columnRange.first ] * nCols_full;

    // Multiply by G - we have already extracted the necessary submatrix
    // with nCols rows from the full version of G, so this is safe
    MATRIX::gemm( baseData, G, multWorkspaceInitial,
                  nCols, nCols_full, nCols, nColsG,
                  true, /* Transpose U_a */
                  false /* Don't transpose G */ );

    // Next, multiply by the appropriate sub-matrix of U_a
    baseOutputData = multWorkspaceFinal;

    baseData = _data + interaction._dataRowList[ rowRange.first ] * nCols_full;

    MATRIX::gemm( baseData, multWorkspaceInitial, baseOutputData,
                  nRows, nCols_full, nCols_full, nColsG,
                  false, false /* No transposition */ );
  }
  else if ( _offDiagonal[ ancestor_idx ]._type == COMPRESSED_NODE )
  {
    // The two interactions from this node involved in forming the
    // Schur complement product
    const SupernodeInteraction  &interaction1 = _offDiagonal[ ancestor_idx ];

    // Multiply by the columnRange subset of the interaction
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( ancestor_idx, G,
                                   // Size of G (it has already been
                                   // compressed)
                                   nCols, nColsG,
                                   multWorkspaceInitial,
                                   // Left, transposed multiplication
                                   true, true,
                                   // No need to extract a submatrix of G
                                   false,
                                   // Provide the row range (row offset not
                                   // needed here, since no extraction is done)
                                   columnRange );

    // Multiply by the rowRange subset of the interaction
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( ancestor_idx, multWorkspaceInitial,
                                   // Size of matrix formed in the last
                                   // multiplication
                                   nCols_full, nColsG,
                                   multWorkspaceFinal,
                                   // Left, untransposed multiplication
                                   true, false,
                                   // No need to extract a submatrix here
                                   // (shouldn't matter for this type of
                                   // multiplication anyways)
                                   false,
                                   // Provide the row range (row offset not
                                   // needed here, since no extraction
                                   rowRange );
  }
  else
  {
    TRACE_ASSERT( NULL, "Invalid interaction type" );
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Given a particular interaction U1 from *this* node, and
// the interaction U_a with an ancestor node, form the matrix
// product Q' * U1 * U_a' (where Q is a given matrix) and
// accumulate it in to the given workspace.
//
// We also require a workspace (subMatrix) to extract the necessary
// rows from Q before multiplying (since our stored version of
// U1 may not contain all rows)
//////////////////////////////////////////////////////////////////////
void Supernode::interactionTransMultiply( int interaction_idx,
                                          int ancestor_idx,
                                          const IntArray &rowsQ,
                                          const Supernode &ancestor,
                                          const Real *Q, int nColsQ,
                                          Real *subMatrix,
                                          Real *multWorkspaceInitial,
                                          Real *multWorkspaceFinal )
{
  int                              nCols = numColumns();
  int                              nRows;

  // Copy rows from Q to submatrix
  BuildInteractionSubMatrix( _offDiagonal[ interaction_idx ], rowsQ,
                             Q, subMatrix, nColsQ );

  nRows = _offDiagonal[ interaction_idx ]._rowList.size();

  if ( _offDiagonal[ ancestor_idx ]._type == STANDARD_NODE )
  {
    Real                          *baseData;

    TRACE_ASSERT( _offDiagonal[ ancestor_idx ]._nodeID == ancestor._nodeID,
                  "Node ID mismatch" );

    baseData = _data + _offDiagonal[ interaction_idx ]._dataRowList[ 0 ] * nCols;

    // Form X = Q' * U1
    MATRIX::gemm( subMatrix, baseData, multWorkspaceInitial,
                  nRows, nColsQ, nRows, nCols,
                  true /* Transpose Q */, false /* Don't transpose U1 */ );

    // Form X * U_a'
    baseData = _data + _offDiagonal[ ancestor_idx ]._dataRowList[ 0 ] * nCols;
    nRows = _offDiagonal[ ancestor_idx ]._rowList.size();

    MATRIX::gemm( multWorkspaceInitial, baseData, multWorkspaceFinal,
                  nColsQ, nCols, nRows, nCols,
                  false /* Don't transpose X */, true /* Transpose U_a */ );
  }
  else if ( _offDiagonal[ ancestor_idx ]._type == COMPRESSED_NODE )
  {
    // The two interactions from this node involved in forming the
    // Schur complement product
    const SupernodeInteraction  &interaction1 = _offDiagonal[ ancestor_idx ];
    const SupernodeInteraction  &interaction2 = _offDiagonal[ interaction_idx ];

    // For now anyways, assume that they are both stored in compressed form
    //
    // TODO: Should this be fixed later?
    TRACE_ASSERT( interaction2._type == COMPRESSED_NODE );

    // Transpose the submatrix to prepare it for compressed multiplication
    MATRIX                       Qtrans( nColsQ, nRows );

    MATRIX::transposeBLAS( Qtrans.data(), subMatrix, nRows, nColsQ );
    subMatrix = Qtrans.data();

    // Form Q' * U1
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( interaction_idx, subMatrix,
                                   // Size of subMatrix (it has already been
                                   // compressed and should have it's
                                   // column number should match the number
                                   // of rows in interaction 1)
                                   nColsQ, nRows,
                                   multWorkspaceInitial,
                                   // Right, untransposed multplication
                                   false, false,
                                   // No need to extract a submatrix
                                   false );

    // Form (Q' * U1) * U_a'
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( ancestor_idx, multWorkspaceInitial,
                                   // Size of the matrix formed in the last
                                   // multiplication
                                   nColsQ, nCols,
                                   multWorkspaceFinal,
                                   // Right, transposed multiplication
                                   false, true,
                                   // No need to extract a submatrix
                                   false );
  }
  else
  {
    TRACE_ASSERT( NULL, "Invalid interaction type" );
  }
}
#endif

#if 0
//////////////////////////////////////////////////////////////////////
// Similar to the above, but only multiplies Q' on the right by a
// subset of the matrix U_a * U_a' (only used for forming pieces of
// the diagonal block of ancestor)
//
// The ranges rowRange and columnRange
// in terms of the actual set of rows stored in U_a.
// That is, if rowRange = [ 2, 3 ] and this node stores rows
// r0, r1, r2, r3, r4, ..., then we are considering the row subset
// [ r2, r3 ] of U_a
//////////////////////////////////////////////////////////////////////
void Supernode::interactionTransMultiply( int ancestor_idx,
                                          const Supernode &ancestor,
                                          const DenseBlock &ancestorBlock,
                                          const Real *Q, int nColsQ,
                                          const IndexRange &rowRange,
                                          const IndexRange &columnRange,
                                          Real *subMatrix,
                                          Real *multWorkspaceInitial,
                                          Real *multWorkspaceFinal )
{
  int                              nCols;
  int                              nCols_full = numColumns();
  int                              nRows;

  const SupernodeInteraction      &interaction = _offDiagonal[ ancestor_idx ];

  nRows = range_size( rowRange );
  nCols = range_size( columnRange );

  // Copy rows from Q to submatrix
  BuildInteractionSubMatrix( interaction, Q, subMatrix, nColsQ,
                             rowRange, ancestorBlock._rowRange.first );

  if ( _offDiagonal[ ancestor_idx ]._type == STANDARD_NODE )
  {
    Real                          *baseData;

    TRACE_ASSERT( _offDiagonal[ ancestor_idx ]._nodeID == ancestor._nodeID,
                  "Node ID mismatch" );

    baseData = _data + interaction._dataRowList[ rowRange.first ] * nCols_full;

    // Form X = Q' * U_a
    MATRIX::gemm( subMatrix, baseData, multWorkspaceInitial,
                  nRows, nColsQ, nRows, nCols_full,
                  true /* Transpose Q */, false /* Don't transpose U1 */ );

    // Form X * U_a'
    baseData = _data + interaction._dataRowList[ columnRange.first ] * nCols_full;
    //nRows = _offDiagonal[ ancestor_idx ]._rowList.size();

    MATRIX::gemm( multWorkspaceInitial, baseData, multWorkspaceFinal,
                  nColsQ, nCols_full, nCols, nCols_full,
                  false /* Don't transpose X */, true /* Transpose U_a */ );
  }
  else if ( _offDiagonal[ ancestor_idx ]._type == COMPRESSED_NODE )
  {
    // Transpose the submatrix to prepare it for compressed multiplication
    MATRIX                       Qtrans( nColsQ, nRows );

    MATRIX::transposeBLAS( Qtrans.data(), subMatrix, nRows, nColsQ );
    subMatrix = Qtrans.data();

    // Form Q' * U_a using rowRange
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( ancestor_idx, subMatrix,
                                   // Size of subMatrix (it has already
                                   // been compressed and should match the
                                   // number 
                                   nColsQ, nRows,
                                   multWorkspaceInitial,
                                   // Right, untransposed multiplication
                                   false, false,
                                   // No need to extract a submatrix
                                   false );

    // Form (Q' * U_a) * U_a'
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( ancestor_idx, multWorkspaceInitial,
                                   // Size of the matrix formed in the last
                                   // multiplication
                                   nColsQ, nCols_full,
                                   multWorkspaceFinal,
                                   // Right, transposed multiplication
                                   false, true,
                                   // No need to extract a submatrix
                                   false );
  }
  else
  {
    TRACE_ASSERT( NULL, "Invalid interaction type" );
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Given the interactions below ancestor_idx in this node (call these
// (U1, U2, ...) and the interaction U_a with the ancestor node,
// form the matrix product U_a * [ U1; U2; ... ]' * G (where G is a
// given matrix) and place the result in the given matrix (doesn't
// do any row selection or scattering so the input and output will
// have to be processed)
//////////////////////////////////////////////////////////////////////
void Supernode::interactionTransMultiply_left( int ancestor_idx,
                                               const Supernode &ancestor,
                                               const Real *G, int nColsG,
                                               Real *multWorkspaceInitial,
                                               Real *multWorkspaceFinal )
{
  int                        nCols = numColumns();

  if ( _offDiagonal[ ancestor_idx ]._type == STANDARD_NODE )
  {
    Real                    *baseData;
    Real                    *baseOutputData;
    int                      nRows = 0;

    // Start by multiplying [ U1; U2; ... ]' * G

    // Count rows from interactions below ancestor_idx in this node
    nRows = countInteractionRows( ancestor_idx + 1, _offDiagonal.size() - 1 );
#if 0
    cout << "Counted " << SDUMP( nRows ) << endl;
#endif

    // Move to the part of _compressedV that we want to multiply with
    baseData
      = _data + _offDiagonal[ ancestor_idx + 1 ]._dataRowList[ 0 ] * nCols;

    MATRIX::gemm( baseData, G, multWorkspaceInitial,
                  nRows, nCols, nRows, nColsG,
                  true, false /* Transpose [ U1; U2; ... ] */ );

    // Multiply by U_a
    baseData = _data + _offDiagonal[ ancestor_idx ]._dataRowList[ 0 ] * nCols;
    nRows = _offDiagonal[ ancestor_idx ]._rowList.size();

    MATRIX::gemm( baseData, multWorkspaceInitial, multWorkspaceFinal,
                  nRows, nCols, nCols, nColsG,
                  false, false /* No transposition */ );
  }
  else if ( _offDiagonal[ ancestor_idx ]._type == COMPRESSED_NODE )
  {
    // The ancestor interaction
    const SupernodeInteraction  &interaction1 = _offDiagonal[ ancestor_idx ];

    int                          nRows;

    // Count rows from interactions below ancestor_idx in this node
    nRows = countInteractionRows( ancestor_idx + 1, _offDiagonal.size() - 1 );
#if 0
    cout << "Counted " << SDUMP( nRows ) << endl;
#endif

    // Transposed multiplication by [ U1; U2; ... ]
    compressedInteractionMultiply_fullBlock( ancestor_idx, G,
                                             // Size of G (assuming it has
                                             // already been compressed)
                                             nRows, nColsG,
                                             multWorkspaceInitial,
                                             // Left, transposed multiplication
                                             true, true,
                                             // No need to extract a submatrix
                                             // of G (should already be done)
                                             false );

    // Untransposed multiplication with U_a
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( ancestor_idx, multWorkspaceInitial,
                                   // Size of the matrix formed in the last
                                   // multiplication
                                   nCols, nColsG,
                                   multWorkspaceFinal,
                                   // Left, untransposed multiplication
                                   true, false,
                                   // No need to extract a submatrix
                                   false );
  }
  else
  {
    TRACE_ASSERT( NULL, "Invalid interaction type" );
  }
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but computes (U1 * U_a')' * G for some
// matrix G.  Requires the row set from the interaction we are
// trying to build
//////////////////////////////////////////////////////////////////////
void Supernode::interactionTransMultiply_left( int interaction_idx,
                                               int ancestor_idx,
                                               const IntArray &rowsG,
                                               const Supernode &ancestor,
                                               const Real *G, int nColsG,
                                               Real *subMatrix,
                                               Real *multWorkspaceInitial,
                                               Real *multWorkspaceFinal )
{
  int                            nCols = numColumns();
  int                            nRows;

  nRows = _offDiagonal[ interaction_idx ]._rowList.size();

  // Copy rows from G to submatrix
  BuildInteractionSubMatrix( _offDiagonal[ interaction_idx ], rowsG,
                             G, subMatrix, nColsG );

  if ( _offDiagonal[ ancestor_idx ]._type == STANDARD_NODE )
  {
    Real                          *baseData;

    TRACE_ASSERT( _offDiagonal[ ancestor_idx ]._nodeID == ancestor._nodeID,
                  "Node ID mismatch" );

    baseData = _data + _offDiagonal[ interaction_idx ]._dataRowList[ 0 ] * nCols;

    // Form X = U1' * G
    MATRIX::gemm( baseData, subMatrix, multWorkspaceInitial,
                  nRows, nCols, nRows, nColsG,
                  true /* Tranpose U1 */, false /* Don't transpose G */ );

    // Form U_a * X
    baseData = _data + _offDiagonal[ ancestor_idx ]._dataRowList[ 0 ] * nCols;
    nRows = _offDiagonal[ ancestor_idx ]._rowList.size();

    MATRIX::gemm( baseData, multWorkspaceInitial, multWorkspaceFinal,
                  nRows, nCols, nCols, nColsG,
                  false /* Don't transpose U_a */,
                  false /* Don't transpose X */ );
  }
  else if ( _offDiagonal[ ancestor_idx ]._type == COMPRESSED_NODE )
  {
    // The two interactions from this node involved in forming the
    // Schur complement product
    const SupernodeInteraction  &interaction1 = _offDiagonal[ ancestor_idx ];
    const SupernodeInteraction  &interaction2 = _offDiagonal[ interaction_idx ];

    // For now anyways, assume that they are both stored in compressed form
    //
    // TODO: Should this be fixed later?
    TRACE_ASSERT( interaction2._type == COMPRESSED_NODE );

    // Form X = U1' * G
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( interaction_idx, subMatrix,
                                   // Size of subMatrix (it has already been
                                   // compressed and should have it's
                                   // column number should match the number
                                   // of rows in interaction 1)
                                   nRows, nColsG,
                                   multWorkspaceInitial,
                                   // Left, transposed multplication
                                   true, true,
                                   // No need to extract a submatrix
                                   false );

    // Form U_a * X
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( ancestor_idx, multWorkspaceInitial,
                                   // Size of the matrix formed in the last
                                   // multiplication
                                   nCols, nColsG,
                                   multWorkspaceFinal,
                                   // Left, untransposed multiplication
                                   true, false,
                                   // No need to extract a submatrix
                                   false );
  }
  else
  {
    TRACE_ASSERT( NULL, "Invalid interaction type" );
  }
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but only multiplies G by a subset of
// the matrix (U_a * U_a')' (only used for forming piees of the
// diagonal block of ancestor)
//
// The ranges rowRange and columnRange
// in terms of the actual set of rows stored in U_a.
// That is, if rowRange = [ 2, 3 ] and this node stores rows
// r0, r1, r2, r3, r4, ..., then we are considering the row subset
// [ r2, r3 ] of U_a
//////////////////////////////////////////////////////////////////////
void Supernode::interactionTransMultiply_left( int ancestor_idx,
                                               const Supernode &ancestor,
                                               const DenseBlock &ancestorBlock,
                                               const Real *G, int nColsG,
                                               const IndexRange &rowRange,
                                               const IndexRange &columnRange,
                                               Real *subMatrix,
                                               Real *multWorkspaceInitial,
                                               Real *multWorkspaceFinal )
{
  int                            nCols;
  int                            nCols_full = numColumns();
  int                            nRows;

  const SupernodeInteraction    &interaction = _offDiagonal[ ancestor_idx ];

  nRows = range_size( rowRange );
  nCols = range_size( columnRange );

  // Copy rows from G to submatrix
  BuildInteractionSubMatrix( interaction, G, subMatrix, nColsG,
                             rowRange,
                             // Use an offset, since G is only as big as
                             // it needs to be (it doesn't fill the full
                             // row space of the interaction)
                             ancestorBlock._rowRange.first );

  if ( interaction._type == STANDARD_NODE )
  {
    Real                          *baseData;

    TRACE_ASSERT( _offDiagonal[ ancestor_idx ]._nodeID == ancestor._nodeID,
                  "Node ID mismatch" );

    // Start by forming 
    //    X = U_a( rowRange, : )' * G
    baseData = _data + interaction._dataRowList[ rowRange.first ] * nCols_full;

    MATRIX::gemm( baseData, subMatrix, multWorkspaceInitial,
                  nRows, nCols_full, nRows, nColsG,
                  true /* Tranpose U1 */, false /* Don't transpose G */ );

    // Form U_a( columnRange, : ) * X
    baseData = _data + interaction._dataRowList[ columnRange.first ] * nCols_full;

    MATRIX::gemm( baseData, multWorkspaceInitial, multWorkspaceFinal,
                  nCols, nCols_full, nCols_full, nColsG,
                  false /* Don't transpose U_a */,
                  false /* Don't transpose X */ );
  }
  else if ( interaction._type == COMPRESSED_NODE )
  {
    // Start by forming 
    //    X = U_a( rowRange, : )' * G
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( ancestor_idx, subMatrix,
                                   // Size of subMatrix (it has already
                                   // been compressed and should match the
                                   // number 
                                   nRows, nColsG,
                                   multWorkspaceInitial,
                                   // Left, transposed multiplication
                                   true, true,
                                   // No need to extract a submatrix
                                   false,
                                   // Provide the row range (row offset not
                                   // needed here, since no extraction is done)
                                   rowRange );

    // Form U_a( columnRange, : ) * X
#ifdef DEBUG_COMPRESSED_MULT
    compressedInteractionMultiplyDebug
#else
    compressedInteractionMultiply
#endif
                                 ( ancestor_idx, multWorkspaceInitial,
                                   // Size of the matrix formed in the
                                   // last multiplication
                                   nCols_full, nColsG,
                                   multWorkspaceFinal,
                                   // Left, untransposed multiplication
                                   true, false,
                                   // No need to extract a submatrix
                                   false,
                                   // Provide the row range (row offset not
                                   // needed here, since no extraction)
                                   columnRange );
  }
  else
  {
    TRACE_ASSERT( NULL, "Invalid interaction type" );
  }
}

//////////////////////////////////////////////////////////////////////
// For in-place decomposition of diagonal blocks:
// Multiplies "previous" compressed off-diagonal block's from this
// node's diagonal by the given matrix
//
// This is part of forming the low-rank decomposition of the
// off-diagonal block indexed by block_idx
//////////////////////////////////////////////////////////////////////
void Supernode::multiplyPreviousInPlaceBlocks( int block_idx,
                                               const Real *G,
                                               int rowsG, int colsG,
                                               Real *output,
                                               bool left, bool transpose )
{
  TRACE_ASSERT( inPlaceDiagonal() );

  _nodeDiagonal.offDiagonalMultiply( block_idx, G, rowsG, colsG, output,
                                     left, transpose );
}

//////////////////////////////////////////////////////////////////////
// For in-place decomposition of diagonal blocks:
// Multiplies (and overwrites) the input matrix with the inverse
// of the lower-triangular matrix located "above" the off-diagonal
// block indexed by block_idx
//////////////////////////////////////////////////////////////////////
void Supernode::inPlaceOffDiagonalSolve( int block_idx,
                                         Real *G, int nRHS,
                                         bool transpose )
{
  TRACE_ASSERT( inPlaceDiagonal() );

  _nodeDiagonal.offDiagonalSolve( block_idx, G, nRHS, transpose );
}

#if 0
//////////////////////////////////////////////////////////////////////
// For the set of interactions (U1, U2, ...) from this node that
// are to be decomposed in to low-rank representations, take the
// entries from the original matrix that belong in U1, U2, ... etc.
// and multiply them by the matrix G.
//////////////////////////////////////////////////////////////////////
void Supernode::baseSystemMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const IntArray &lowRankScatteredMap,
                                    const IntArray &lowRankBlocks,
                                    const Real *G, int nColsG,
                                    Real *multWorkspaceFinal )
{
  int                            row_idx;
  int                            workspace_row_idx;
  int                            interaction_idx;
  int                            totalRows = 0;

  // Clear a large enough chunk of the workspace
  for ( int block_idx = 0; block_idx < lowRankBlocks.size();
        block_idx++ )
  {
    interaction_idx = lowRankBlocks[ block_idx ];

    totalRows += _offDiagonal[ interaction_idx ]._rowList.size();
  }

  MATRIX::clear( multWorkspaceFinal, totalRows, nColsG );

  // For each column in the matrix, walk through its row entries.
  // If the row entry has a valid index in lowRankScatteredMap,
  // multiply by the corresponding row in G and add to the workspace
  // accordingly.
  for ( int col_idx = _columnRange.first; col_idx <= _columnRange.second;
        col_idx++ )
  for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
        row_ptr++ )
  {
    row_idx = A._i[ row_ptr ];
    workspace_row_idx = lowRankScatteredMap[ row_idx ];

    if ( workspace_row_idx >= 0 )
    {
      MATRIX::axpy( multWorkspaceFinal + workspace_row_idx * nColsG,
                    // Align with correct row
                    G + ( col_idx - _columnRange.first ) * nColsG,
                    1 /* one row */, nColsG,
                    A._x[ row_ptr ] /* Multiply with matrix entry */ );
    }
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// For a specific interaction in this node, take the entries from
// the original matrix that belong in its row/column space and
// multiply them by the matrix G.
//////////////////////////////////////////////////////////////////////
void Supernode::baseSystemMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const vector<Supernode> &nodes,
                                    int interaction_idx,
                                    const Real *G, int nColsG,
                                    Real *multWorkspaceFinal )
{
  int                            row_idx;
  int                            workspace_row_idx;
  int                            totalRows = 0;
  int                            this_row_idx;
  int                            startRow, endRow;

  SupernodeInteraction          &interaction = _offDiagonal[ interaction_idx ];

  const Supernode               &ancestor = nodes[ interaction._nodeID ];

  startRow = ancestor._columnRange.first;
  endRow = ancestor._columnRange.second;

  MATRIX::clear( multWorkspaceFinal, interaction._rowList.size(), nColsG );

  // For each column in the matrix, walk through its row entries.
  // If the row entry has a valid index in lowRankScatteredMap,
  // multiply by the corresponding row in G and add to the workspace
  // accordingly.
  for ( int col_idx = _columnRange.first; col_idx <= _columnRange.second;
        col_idx++ )
  {
    this_row_idx = 0;

    for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
          row_ptr++ )
    {
      row_idx = A._i[ row_ptr ];

      if ( row_idx < startRow || row_idx > endRow )
      {
        continue;
      }

      // Find the corresponding interaction row
      while ( this_row_idx < interaction._rowList.size()
           && interaction._rowList[ this_row_idx ] + startRow != row_idx )
      {
        this_row_idx++;
      }

      TRACE_ASSERT( this_row_idx < interaction._rowList.size(),
                    "No corresponding interaction row found" );

      MATRIX::axpy( multWorkspaceFinal + this_row_idx * nColsG,
                    // Align with correct row
                    G + ( col_idx - _columnRange.first ) * nColsG,
                    1 /* one row */, nColsG,
                    A._x[ row_ptr ] /* Multiply with matrix entry */ );

#if 0
      if ( workspace_row_idx >= 0 )
      {
        MATRIX::axpy( multWorkspaceFinal + workspace_row_idx * nColsG,
                      // Align with correct row
                      G + ( col_idx - _columnRange.first ) * nColsG,
                      1 /* one row */, nColsG,
                      A._x[ row_ptr ] /* Multiply with matrix entry */ );
      }
#endif
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Multiplies by the entries of A corresponding to rows in this node's
// full off-diagonal (rather than just a single interaction which is
// handled by the function above)
//////////////////////////////////////////////////////////////////////
void Supernode::baseSystemMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const std::vector<Supernode> &nodes,
                                    const Real *G, int nColsG,
                                    Real *multWorkspaceFinal )
{
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    baseSystemMultiply( A, nodes, interaction_idx, G, nColsG,
                        multWorkspaceFinal );

    multWorkspaceFinal += interaction._rowList.size() * nColsG;
  }
}

//////////////////////////////////////////////////////////////////////
// Multiply by the base system of a sub-matrix in this node's main
// diagonal
//////////////////////////////////////////////////////////////////////
void Supernode::baseSystemMultiply( const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    const Real *G, int nColsG,
                                    const IndexRange &rowRange,
                                    const IndexRange &columnRange,
                                    Real *multWorkspaceFinal )
{
  int                            row_idx;
  int                            workspace_row_idx;
  int                            totalRows = 0;
  int                            startRow, endRow;
  int                            startCol, endCol;
  int                            nRows;

  startRow = _columnRange.first + rowRange.first;
  endRow = _columnRange.first + rowRange.second;

  startCol = _columnRange.first + columnRange.first;
  endCol = _columnRange.first + columnRange.second;

  nRows = range_size( rowRange );

  MATRIX::clear( multWorkspaceFinal, nRows, nColsG );

  // For each column in the matrix, walk through its row entries.
  // If the row entry has a valid index in lowRankScatteredMap,
  // multiply by the corresponding row in G and add to the workspace
  // accordingly.
  for ( int col_idx = startCol; col_idx <= endCol; col_idx++ )
  {
    for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
          row_ptr++ )
    {
      row_idx = A._i[ row_ptr ];

      if ( row_idx < startRow || row_idx > endRow )
      {
        continue;
      }

      // Put the row index in to the range of the output matrix
      row_idx -= startRow;

      MATRIX::axpy( multWorkspaceFinal + row_idx * nColsG,
                    // Align with correct row
                    G + ( col_idx - startCol ) * nColsG,
                    1 /* one row */, nColsG,
                    A._x[ row_ptr ] /* Multiply with matrix entry */ );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// For a single interaction from this node (U1) which is to be
// decomposed in to a low-rank representation, take the entries
// of the original matrix that belong in the U1 block and form
// the product Q' * U1.
//////////////////////////////////////////////////////////////////////
void Supernode::baseSystemTransMultiply(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              int interaction_idx,
                              const vector<Supernode> &nodes,
                              const Real *Q, int nColsQ,
                              Real *multWorkspaceFinal,
                              bool transposeQ )
{
  int                            row_idx, q_row_idx;
  int                            nCols = numColumns();
  const SupernodeInteraction    &interaction = _offDiagonal[ interaction_idx ];
  const IntArray                &rowList = interaction._rowList;
  int                            inputOffset, outputOffset;
  int                            qLDA;

  // Get the first row of the block U1
  int                            ancestor_idx = interaction._nodeID;
  const Supernode               &ancestor = nodes[ ancestor_idx ];
  int                            rowStart = ancestor._columnRange.first;
  int                            rowEnd = ancestor._columnRange.second;

  // If we need to transpose Q (ie. Q is stored in its non-tranposed form),
  // then we will need to extract rows from Q, in which case the leading
  // dimension is 1.  Otherwise, extracting rows implies extracting a
  // column from Qtrans, in which case the leading dimension is the
  // number of rows in Q.
  qLDA = transposeQ ? 1 : (int)rowList.size();

  // Clear the workspace
  MATRIX::clear( multWorkspaceFinal, nColsQ, nCols );

  // Form the result column by column by multiplying Q' with
  // each column independently
  for ( int col_idx = _columnRange.first; col_idx <= _columnRange.second;
        col_idx++ )
  {
    q_row_idx = 0;

    for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
          row_ptr++ )
    {
      row_idx = A._i[ row_ptr ];

      // Skip rows outside of the desired column range
      if ( row_idx < rowStart || row_idx > rowEnd )
      {
        continue;
      }

      // Find the corresponding row in Q
      while ( q_row_idx < rowList.size()
           && rowList[ q_row_idx ] + rowStart != row_idx )
      {
        q_row_idx++;
      }

      TRACE_ASSERT( q_row_idx < rowList.size(),
                    "No corresponding row found in the interaction" );

      // Just move to the start of a column
      outputOffset = col_idx - _columnRange.first;

      if ( transposeQ )
      {
        // Move to the start of a row
        inputOffset = q_row_idx * nColsQ;
      }
      else
      {
        // Move to the start of a column, since Q' is stored
        inputOffset = q_row_idx;
      }

      // Take a row from Q, and add it as a column in the output
      MATRIX::axpy( multWorkspaceFinal + outputOffset,
                    Q + inputOffset,
                    nColsQ, 1,
                    A._x[ row_ptr ] /* Scale by matrix entry */,
                    nCols /* Copy to column of output */,
                    qLDA /* Copy row of Q */ );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Multiply a matrix on the right by the base system of a sub-matrix
// in this node's main diagonal.
//
// That is, form the product Q' * U1 (as in the above function)
//////////////////////////////////////////////////////////////////////
void Supernode::baseSystemTransMultiply(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              const Real *Q, int nColsQ,
                              const IndexRange &rowRange,
                              const IndexRange &columnRange,
                              Real *multWorkspaceFinal,
                              bool transposeQ )
{
  int                            row_idx;
  int                            nCols;
  int                            nRows;
  //const SupernodeInteraction    &interaction = _offDiagonal[ interaction_idx ];
  //const IntArray                &rowList = interaction._rowList;
  int                            inputOffset, outputOffset;
  int                            qLDA;

  // Get the first row of the block U1
  //int                            ancestor_idx = interaction._nodeID;
  //const Supernode               &ancestor = nodes[ ancestor_idx ];
  int                            rowStart, rowEnd;
  int                            colStart, colEnd;

  rowStart = _columnRange.first + rowRange.first;
  rowEnd = _columnRange.first + rowRange.second;

  colStart = _columnRange.first + columnRange.first;
  colEnd = _columnRange.first + columnRange.second;

  nRows = range_size( rowRange );
  nCols = range_size( columnRange );

  // If we need to transpose Q (ie. Q is stored in its non-tranposed form),
  // then we will need to extract rows from Q, in which case the leading
  // dimension is 1.  Otherwise, extracting rows implies extracting a
  // column from Qtrans, in which case the leading dimension is the
  // number of rows in Q.
  qLDA = transposeQ ? 1 : nRows;

  // Clear the workspace
  MATRIX::clear( multWorkspaceFinal, nColsQ, nCols );

  // Form the result column by column by multiplying Q' with
  // each column independently
  for ( int col_idx = colStart; col_idx <= colEnd; col_idx++ )
  {
    for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
          row_ptr++ )
    {
      row_idx = A._i[ row_ptr ];

      // Skip rows outside of the desired column range
      if ( row_idx < rowStart || row_idx > rowEnd )
      {
        continue;
      }

      // Just move to the start of a column
      outputOffset = col_idx - colStart;

      // Make the row index relative to the row space we are working in
      row_idx -= rowStart;

      if ( transposeQ )
      {
        // Move to the start of a row
        inputOffset = row_idx * nColsQ;
      }
      else
      {
        // Move to the start of a column, since Q' is stored
        inputOffset = row_idx;
      }

      // Take a row from Q, and add it as a column in the output
      MATRIX::axpy( multWorkspaceFinal + outputOffset,
                    Q + inputOffset,
                    nColsQ, 1,
                    A._x[ row_ptr ] /* Scale by matrix entry */,
                    nCols /* Copy to column of output */,
                    qLDA /* Copy row of Q */ );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but multiplies the transpose of the interaction,
// with some other matrix.  (ie. U1' * Q)
//////////////////////////////////////////////////////////////////////
void Supernode::baseSystemTransMultiply_left(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              int interaction_idx,
                              const vector<Supernode> &nodes,
                              const Real *G, int nColsG,
                              Real *multWorkspaceFinal,
                              bool clearOutput )
{
  int                            nRows;
  int                            nCols;
  int                            inputOffset, outputOffset;
  int                            row_idx, g_row_idx;

  SupernodeInteraction          &interaction = _offDiagonal[ interaction_idx ];

  // Get the first row of the block U1
  int                            ancestor_idx = interaction._nodeID;
  const Supernode               &ancestor = nodes[ ancestor_idx ];
  int                            rowStart = ancestor._columnRange.first;
  int                            rowEnd = ancestor._columnRange.second;

  nRows = interaction._rowList.size();
  nCols = numColumns();

  // Clear the workspace
  if ( clearOutput ) {
    MATRIX::clear( multWorkspaceFinal, nCols, nColsG );
  }

  // Form the result column by column by multiplying Q' with
  // each column independently
  for ( int col_idx = _columnRange.first; col_idx <= _columnRange.second;
        col_idx++ )
  {
    g_row_idx = 0;

    for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
          row_ptr++ )
    {
      row_idx = A._i[ row_ptr ];

      // Skip rows outside of the desired column range
      if ( row_idx < rowStart || row_idx > rowEnd )
      {
        continue;
      }

      // Find the corresponding row index in to G
      while ( g_row_idx < interaction._rowList.size()
          && interaction._rowList[ g_row_idx ] + rowStart != row_idx )
      {
        g_row_idx++;
      }

      TRACE_ASSERT( g_row_idx < interaction._rowList.size(),
                    "No corresponding interaction row" );

      inputOffset = g_row_idx * nColsG;
      outputOffset = ( col_idx - _columnRange.first ) * nColsG;

      // Scaled row add
      MATRIX::axpy( multWorkspaceFinal + outputOffset, G + inputOffset,
                    1, nColsG, A._x[ row_ptr ] );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Multiplies by the system A in the row space given by this node's
// full off-diagonal block (as opposed to just a single interaction block)
//////////////////////////////////////////////////////////////////////
void Supernode::baseSystemTransMultiply_left(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              const std::vector<Supernode> &nodes,
                              const Real *G, int nColsG,
                              Real *multWorkspaceFinal )
{
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ ) 
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    baseSystemTransMultiply_left( A, interaction_idx, nodes, G, nColsG,
                                  multWorkspaceFinal,
                                  // Only clear on first pass
                                  interaction_idx == 0 );

    G += interaction._rowList.size() * nColsG;
  }
}

//////////////////////////////////////////////////////////////////////
// Does the same as the function above, but uses the base system
// of a sub-matrix in this node's main diagonal
//////////////////////////////////////////////////////////////////////
void Supernode::baseSystemTransMultiply_left(
                              const SPARSE_MATRIX::SparseColumnMatrix &A,
                              const Real *G, int nColsG,
                              const IndexRange &rowRange,
                              const IndexRange &columnRange,
                              Real *multWorkspaceFinal )
{
  int                            nRows;
  int                            nCols;
  int                            inputOffset, outputOffset;
  int                            row_idx;

  // Get the first row of the block U1
  //int                            ancestor_idx = interaction._nodeID;
  //const Supernode               &ancestor = nodes[ ancestor_idx ];
  int                            rowStart, rowEnd;
  int                            colStart, colEnd;

  rowStart = _columnRange.first + rowRange.first;
  rowEnd = _columnRange.first + rowRange.second;

  colStart = _columnRange.first + columnRange.first;
  colEnd = _columnRange.first + columnRange.second;

  nRows = range_size( rowRange );
  nCols = range_size( columnRange );

  // Clear the workspace
  MATRIX::clear( multWorkspaceFinal, nCols, nColsG );

  // Form the result column by column by multiplying Q' with
  // each column independently
  for ( int col_idx = colStart; col_idx <= colEnd; col_idx++ )
  {
    for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
          row_ptr++ )
    {
      row_idx = A._i[ row_ptr ];

      // Skip rows outside of the desired column range
      if ( row_idx < rowStart || row_idx > rowEnd )
      {
        continue;
      }

      // Put the row index in the row space we are working in
      row_idx -= rowStart;

      inputOffset = row_idx * nColsG;
      outputOffset = ( col_idx - colStart ) * nColsG;

      // Scaled row add
      MATRIX::axpy( multWorkspaceFinal + outputOffset, G + inputOffset,
                    1, nColsG, A._x[ row_ptr ] );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// For a given low rank block (indexed by block_idx) in the diagonal,
// multiply the matrix G with the sum of all diagonal contributions
// resulting from slack variable introduction and accumulates the
// result in multWorkspaceFinal
//
// We need to call expandCompressedInteractions before this in order
// to fill in expansionWorkspace
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalContributionMult( int block_idx,
                                          const Real *G, int nColsG,
                                          const Real *extraData,
                                          Real *expansionWorkspace,
                                          Real *multWorkspaceInitial,
                                          Real *multWorkspaceFinal,
                                          bool transpose,
                                          bool left )
{
  diagonalContributionMult_offDiagonalContributions(
                                  block_idx, G, nColsG, extraData,
                                  multWorkspaceInitial, multWorkspaceFinal,
                                  expansionWorkspace,
                                  transpose, left );

  diagonalContributionMult_diagonalContributions(
                                  block_idx, G, nColsG, extraData,
                                  multWorkspaceInitial, multWorkspaceFinal,
                                  transpose, left );
}

//////////////////////////////////////////////////////////////////////
// For a given low rank block (indexed by block_idx) in the diagonal,
// multiply the matrix G with the sum of all diagonal contributions
// resulting from slack variable introduction and accumulates the
// result in multWorkspaceFinal
//
// We need to call expandCompressedInteractions before this in order
// to fill in expansionWorkspace
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalContributionMult(
                                      int block_idx,
                                      const Real *G, int nColsG,
                                      const Real *extraData,
                                      WorkspaceManager<Real> &workspaceManager,
                                      Real *multWorkspaceFinal,
                                      bool transpose,
                                      bool left )
{
  diagonalContributionMult_offDiagonalContributions(
                                  block_idx, G, nColsG, extraData,
                                  workspaceManager, multWorkspaceFinal,
                                  transpose, left );

  diagonalContributionMult_diagonalContributions(
                                  block_idx, G, nColsG, extraData,
                                  workspaceManager, multWorkspaceFinal,
                                  transpose, left );
}

//////////////////////////////////////////////////////////////////////
// Given a workspace of computed multiplications of low-rank
// blocks with a random matrix with nColsG columns, compute
// a QR factorization for each interaction.
//////////////////////////////////////////////////////////////////////
void Supernode::lowRankQR( const IntArray &lowRankBlocks,
                           const BoolArray &lowRankBlockActive,
                           Real *decompWorkspace, const IntArray &nBlockCols,
                           Real *qrExtraData, Real *qrWorkspace, int workSize,
                           bool transpose )
{
  int                            nRows;
  int                            info;
  Real                          *baseData = decompWorkspace;
  int                            interaction_idx;
  int                            nColsG;
  int                            nCols = numColumns();

  for ( int block_idx = 0; block_idx < lowRankBlocks.size(); block_idx++ )
  {
    if ( !lowRankBlockActive[ block_idx ] )
    {
      continue;
    }

    nColsG = nBlockCols[ block_idx ];

    interaction_idx = lowRankBlocks[ block_idx ];

    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    if ( (interaction._type == STANDARD_NODE && !interaction._active)
      || interaction._type == COMPRESSED_NODE )
    {
      // We are replacing this node with slack variables, so we
      // must have already formed the multiple of its contents
      // and a random matrix and placed this in to decomp workspace
      nRows = transpose ? nCols : interaction._rowList.size();

      info = MATRIX::qr( baseData, nRows, nColsG,
                         qrExtraData, qrWorkspace, workSize );

      if ( info != 0 )
      {
        printf( "QR parameter %d had an illegal value\n", -info );
        cout << SDUMP( nRows ) << endl;
        cout << SDUMP( nColsG ) << endl;
        cout << SDUMP( baseData ) << endl;
        cout << SDUMP( _nodeID ) << endl;

        MATRIX tmp( nRows, nColsG, baseData );
        tmp.write( "badMatrix.matrix" );

        TRACE_ASSERT( info == 0, "QR factorization failed" );
      }

      // Get the actual factor from this
      info = MATRIX::extractQRfactor( baseData, nRows, nColsG,
                                      qrExtraData, qrWorkspace, workSize );

      if ( info != 0 )
      {
        printf( "QR extraction parameter %d had an illegal value\n", -info );
        cout << SDUMP( nRows ) << endl;
        cout << SDUMP( nColsG ) << endl;
        cout << SDUMP( baseData ) << endl;
        cout << SDUMP( numColumns() ) << endl;

        MATRIX tmp( nRows, nColsG, baseData );
        tmp.write( "badMatrixExtract.matrix" );

        TRACE_ASSERT( info == 0, "Q extraction failed" );
      }

      baseData += nRows * nColsG;
    }
    else {
      TRACE_ASSERT( NULL, "Invalid interaction type" );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// QR factorization like the function above but assuming that we
// are decomposing the full off-diagonal of this node all at once
//////////////////////////////////////////////////////////////////////
void Supernode::lowRankQR( Real *decompWorkspace, int rank,
                           Real *qrExtraData, Real *qrWorkspace, int workSize,
                           bool transpose )
{
  int                        nRows;
  int                        info;
  Real                      *baseData = decompWorkspace;
  int                        nCols = numColumns();

  nRows = transpose ? nCols : countLowRankRows();

  info = MATRIX::qr( baseData, nRows, rank,
                     qrExtraData, qrWorkspace, workSize );

  if ( info != 0 ) {
    printf( "Full decomp. QR parameter %d had an illegal value\n", -info );
    cout << SDUMP( nRows ) << endl;
    cout << SDUMP( rank ) << endl;
    cout << SDUMP( baseData ) << endl;
    cout << SDUMP( _nodeID ) << endl;

    MATRIX tmp( nRows, rank, baseData );
    tmp.write( "badMatrix.matrix" );

    printf( "Enter a character to continue: " );
    char waitChar;
    cin >> waitChar;
    printf( "\n" );

    TRACE_ASSERT( info == 0, "QR factorization failed" );
  }

  // Pull out the actual factor
  info = MATRIX::extractQRfactor( baseData, nRows, rank,
                                  qrExtraData, qrWorkspace, workSize );

#if 0
  // FIXME: debugging
  MATRIX::write( baseData, nRows, rank, "postQR.matrix" );
#endif

  if ( info != 0 ) {
    printf( "QR extraction parameter %d had an illegal value\n", -info );
    cout << SDUMP( nRows ) << endl;
    cout << SDUMP( rank ) << endl;
    cout << SDUMP( baseData ) << endl;
    cout << SDUMP( numColumns() ) << endl;

    MATRIX tmp( nRows, rank, baseData );
    tmp.write( "badMatrixExtract.matrix" );

    TRACE_ASSERT( info == 0, "Q extraction failed" );
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but only performst the QR factorization for
// a single low rank block in the node's diagonal
//////////////////////////////////////////////////////////////////////
void Supernode::lowRankQR( int block_idx, Real *decompWorkspace, int nColsG,
                           Real *qrExtraData, Real *qrWorkspace,
                           int workSize, bool transpose )
{
  const DenseBlock          &block = _diagonalLowRankBlocks[ block_idx ];

  int                        nRows = transpose ? block.numColumns()
                                               : block.numRows();
  int                        info;

  info = MATRIX::qr( decompWorkspace, nRows, nColsG,
                     qrExtraData, qrWorkspace, workSize );

  if ( info != 0 )
  {
    printf( "Diagonal block QR parameter %d had an illegal value\n", -info );
    cout << SDUMP( nRows ) << endl;
    cout << SDUMP( nColsG ) << endl;
    cout << SDUMP( decompWorkspace ) << endl;
    cout << SDUMP( _nodeID ) << endl;

    MATRIX tmp( nRows, nColsG, decompWorkspace );
    tmp.write( "badMatrix.matrix" );

    printf( "Enter a character to continue: " );
    char waitChar;
    cin >> waitChar;
    printf( "\n" );

    TRACE_ASSERT( info == 0, "QR factorization failed" );
  }

  // Get the actual factor
  info = MATRIX::extractQRfactor( decompWorkspace, nRows, nColsG,
                                  qrExtraData, qrWorkspace, workSize );

  TRACE_ASSERT( info == 0, "Q extraction failed" );
}

//////////////////////////////////////////////////////////////////////
// Suppose that the given interaction in this node has been
// decomposed into V * U'.  This function assigns extra data
// pointers to interactions in a factorization and assigns
// the data in V and U to the correct locations.
//////////////////////////////////////////////////////////////////////
void Supernode::assignExtraData( int interaction_idx,
                                 vector<Supernode> &nodes,
                                 const Real *V, const Real *Utrans,
                                 int rank, Real *extraData,
                                 long int &offset, long int &remainingSize,
                                 Real *copyWorkspace )
{
  if ( _offDiagonal[ interaction_idx ]._type == STANDARD_NODE ) {
    assignNewExtendedNode( interaction_idx, nodes, V, Utrans, rank,
                           extraData, offset, remainingSize,
                           copyWorkspace );
  }
  else if ( _offDiagonal[ interaction_idx ]._type == COMPRESSED_NODE ) {
    assignCompressedInteraction( interaction_idx, V, Utrans, rank,
                                 extraData, offset, remainingSize );
  }
  else {
    TRACE_ASSERT( NULL, "Invalid interaction type" );
  }
}

//////////////////////////////////////////////////////////////////////
// Assigns a single low-rank decomposition for all of this node's
// off-diagonal content
//////////////////////////////////////////////////////////////////////
void Supernode::assignExtraData( const Real *V, const Real *Utrans,
                                 int rank, Real *extraData,
                                 long int &offset, long int &remainingSize )
{
  int                        nCols = numColumns();
  int                        nRows;

  Real                      *baseData;
  
  nRows = countInteractionRows( 0, _offDiagonal.size() - 1 );

  _compressedRank = rank;
  _compressedV = extraData + offset;
  offset += nRows * rank;
  _compressedU = extraData + offset;
  offset += nCols * rank;

  remainingSize -= rank * ( nRows + nCols );

  // Copy
  MATRIX::copy( _compressedV, V, nRows, rank );

  MATRIX::transposeBLAS( _compressedU, Utrans, rank, nCols );

  baseData = _compressedV;

  // Also, assign things to each individual interaction, since this
  // will be convenient in other parts of our code
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    SupernodeInteraction    &interaction = _offDiagonal[ interaction_idx ];

    TRACE_ASSERT( interaction._type == COMPRESSED_NODE,
                  "Invalid interaction type" );

    interaction._compressedV = baseData;

    // All interactions have the same column basis
    interaction._compressedU = _compressedU;

#if 0
    // FIXME: debugging
    {
      char buf[1024];

      sprintf( buf, "super_numeric/node_%d_interaction_%d_V.matrix",
               _nodeID, interaction_idx );

      MATRIX::write( interaction._compressedV, interaction._rowList.size(),
                     rank, buf );

      sprintf( buf, "super_numeric/node_%d_interaction_%d_U.matrix",
               _nodeID, interaction_idx );

      MATRIX::write( interaction._compressedU, nCols, rank, buf );
    }
#endif

    interaction._numExtendedRows = rank;
    //interaction._numExtendedRows = 0;

    baseData += interaction._rowList.size() * rank;
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but for a low rank block in the diagonal
//////////////////////////////////////////////////////////////////////
void Supernode::assignExtraDiagonalData(
                              int block_idx,
                              std::vector<Supernode> &nodes,
                              const Real *V, const Real *Utrans,
                              int rank, Real *extraData,
                              long int &offset, long int &remainingSize )
{
  if ( inPlaceDiagonal() ) {
    assignCompressedInteraction_diagonal( block_idx, V, Utrans,
                                          rank, extraData,
                                          offset, remainingSize );
  } else {
    assignNewExtendedNode_diagonal( block_idx, nodes, V, Utrans,
                                    rank, extraData,
                                    offset, remainingSize );
  }
}

//////////////////////////////////////////////////////////////////////
// Initialization for an extended node - to be called after all nodes
// have their size set
//////////////////////////////////////////////////////////////////////
void Supernode::initializeExtendedNode( const vector<Supernode> &nodes,
                                        int node_idx,
                                        Real *extraData,
                                        long int &offset,
                                        long int &remainingSize )
{
  int                        nCols = numExtendedColumns();
  int                        blockSize;
  long int                   totalSize = 0;

  TRACE_ASSERT( _type == EXTENDED_NODE,
                "Extended initialization routine called on standard node" );
  TRACE_ASSERT( nCols > 0, "Node column size not set" );
  TRACE_ASSERT( &nodes[ node_idx ] == this, "Wrong node index" );
  TRACE_ASSERT( node_idx > 0, "Not an extended node" );

  // Set the column range for this node
  const Supernode           &lastNode = nodes[ node_idx - 1 ];

  TRACE_ASSERT( lastNode.numColumns() > 0,
                "Previous node's column range is not initialized" );

  _columnRange.first = lastNode._columnRange.second + 1;
  _columnRange.second = lastNode._columnRange.second + nCols;

  // Set data offset for this node's diagonal block
  _extendedDataOffset = offset;
  offset += nCols * nCols;
  remainingSize -= nCols * nCols;

  totalSize = nCols * nCols;

  // Set data offsets for all interactions
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    SupernodeInteraction    &interaction = _offDiagonal[ interaction_idx ];
    const Supernode         &ancestor = nodes[ interaction._nodeID ];

    TRACE_ASSERT( interaction._type == EXTENDED_NODE );

#if 0
    // FIXME: debugging
    if ( interaction._numExtendedRows != ancestor.numExtendedColumns() )
    {
      cerr << SDUMP( interaction._numExtendedRows ) << endl;
      cerr << SDUMP( ancestor.numExtendedColumns() ) << endl;
    }
#endif

    TRACE_ASSERT( interaction._numExtendedRows == ancestor.numExtendedColumns(),
                  "Extended variable size mismatch" );

    blockSize = nCols * interaction._numExtendedRows;

    interaction._extendedDataOffset = offset;
    offset += blockSize;
    remainingSize -= blockSize;
    totalSize += blockSize;
  }

  MATRIX::clear( extraData + _extendedDataOffset, totalSize, 1 );
  
  // Put an identity matrix in the diagonal block
  copyMatrixData_extended( extraData );
}

//////////////////////////////////////////////////////////////////////
// Computes the Cholesky factor piece for this node.
// This assumes that the node has already been set up (eg. via
// update matrices, etc.)
//////////////////////////////////////////////////////////////////////
void Supernode::factor( Real *extraData, int writeSize,
                        WorkspaceManager<Real> *workspaceManager )
{
  if ( _useLDL ) {
    TRACE_ASSERT( workspaceManager != NULL );

    factorLDL( extraData, *workspaceManager );
  }
  else {
    factorCholesky( extraData, writeSize, workspaceManager );
  }
}

//////////////////////////////////////////////////////////////////////
// Solves for the off-diagonal contents of the Cholesky factor for
// this node.
//////////////////////////////////////////////////////////////////////
void Supernode::offDiagonalSolve( bool solveCompressedBlocks,
                                  Real *extraData, int writeSize,
                                  WorkspaceManager<Real> *workspaceManager )
{
  if ( _useLDL ) {
    TRACE_ASSERT( workspaceManager != NULL );

    LDLOffDiagonalSolve( extraData, *workspaceManager );
  }
  else {
    choleskyOffDiagonalSolve( solveCompressedBlocks, extraData,
                              writeSize, workspaceManager );
  }
}

//////////////////////////////////////////////////////////////////////
// Computes dense Cholesky factor(s) for the diagonal block
//////////////////////////////////////////////////////////////////////
void Supernode::factorDiagonal( Real *extraData,
                                WorkspaceManager<Real> *workspaceManager )
{
  Real                          *diagonalData;
  Real                          *diagonalWorkspace = NULL;
  int                            nCols;
  int                            info;

  diagonalData = ( _type == STANDARD_NODE ) ? _data
                                            : extraData + _extendedDataOffset;

  if ( !lowRankDiagonal() )
  {
    // We store the full diagonal
    nCols = numColumns();

    if ( _type == EXTENDED_NODE && workspaceManager )
    {
      RealWorkspace          workspace( *workspaceManager,
                                        IndexPair( nCols, nCols ) );

      diagonalWorkspace = workspace.workspaceData( 0 );

      // Make a copy of the original matrix
      MATRIX::copy( diagonalWorkspace, diagonalData, nCols, nCols );
    }

    if ( _type == EXTENDED_NODE ) {
      MATRIX U, V;
      int modificationSize;
      modificationSize = MATRIX::modifiedCholesky( diagonalData, U, V, nCols );

      if ( modificationSize > 0 ) {
        printf( "Adding rank %d modification to diagonal block\n",
                modificationSize );
        printf( "Block column [ %d, %d ]\n",
                _columnRange.first, _columnRange.second );
      }
    }
    else {
      info = MATRIX::cholesky( diagonalData, nCols );

      if ( info < 0 ) {
        TRACE_ASSERT( NULL,
                      "ERROR: Invalid input to Cholesky solver" );
      }
      else if ( info > 0 ) {
        MATRIX::write( diagonalData, nCols, nCols,
                       "super_numeric/badDiagonal.matrix" );

        printf( "Indefinite matrix found in diagonal of node %d\n", _nodeID );

        TRACE_ASSERT( NULL, "Matrix is not positive definite" );
        abort();
      }
    }
  }
  else
  {
    // Factor each diagonal block independently
    for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ )
    {
      nCols = range_size( _diagonalBlocks[ block_idx ]._columnRange );

      info = MATRIX::cholesky( diagonalData, nCols );

      if ( info < 0 )
      {
        TRACE_ASSERT( NULL,
                      "ERROR: Invalid input to Cholesky solver" );
      }
      else if ( info > 0 )
      {
        cerr << endl << "Block error in node " << _nodeID << endl;
        TRACE_ASSERT( NULL,
                      "ERROR: Matrix not positive definite" );
      }

      diagonalData += nCols * nCols;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Solves a system using this node's diagonal (assumed to be already
// factored)
//
// Specify whether to transpose the system, and whether the solve
// should be a left or right solve
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalSolve( Real *rhs, int nRHS, Real *extraData,
                               bool transpose, bool left,
                               int startBlock, int endBlock ) const
{
  if ( inPlaceDiagonal() ) {
    diagonalSolve_inPlace( rhs, nRHS, transpose, left );
  } else {
    diagonalSolve_extendedVariable( rhs, nRHS, extraData,
                                    transpose, left,
                                    startBlock, endBlock );
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but for left-side solves with vector data
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalSolve( Real *rhs, Real *extraData,
                               bool transpose,
                               int startBlock, int endBlock ) const
{
  if ( inPlaceDiagonal() ) {
    diagonalSolve_inPlace( rhs, transpose );
  } else {
    // No special version of the function here
    diagonalSolve_extendedVariable( rhs, 1, extraData,
                                    transpose, true /* left */,
                                    startBlock, endBlock );
  }
}

//////////////////////////////////////////////////////////////////////
// Returns the number of rows in this node's off diagonal
// which will be replaced with low-rank decompositions
//////////////////////////////////////////////////////////////////////
int Supernode::countLowRankRows() const
{
  int                            numLowRankRows = 0;

  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    if ( !interaction._active || interaction._type == COMPRESSED_NODE )
    {
      numLowRankRows += interaction._rowList.size();
    }
  }

  return numLowRankRows;
}

//////////////////////////////////////////////////////////////////////
// Returns the full number of rows in this node's off diagonal
// which will be replace with low-rank decompositions, assuming
// we include zero rows
//////////////////////////////////////////////////////////////////////
int Supernode::countLowRankRowsFull( const vector<Supernode> &nodes ) const
{
  int                            numLowRankRows = 0;

  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    if ( !interaction._active )
    {
      const Supernode           &ancestor = nodes[ interaction._nodeID ];

      numLowRankRows += ancestor.numColumns();
    }
  }

  return numLowRankRows;
}

//////////////////////////////////////////////////////////////////////
// The number of uncompressed entries in the standard part of the
// system stored by this node
//////////////////////////////////////////////////////////////////////
size_t Supernode::uncompressedEntries() const
{
  size_t                     entries = 0;
  int                        nCols = numColumns();

  if ( _diagonalBlocks.size() > 0 )
  {
    for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ )
    {
      const DenseBlock      &block = _diagonalBlocks[ block_idx ];

      entries += block.numRows() * block.numColumns();
    }
  }
  else
  {
    entries += nCols * nCols;
  }

  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    if ( interaction._type == EXTENDED_NODE || !interaction._active
      || interaction._type == COMPRESSED_NODE )
    {
      continue;
    }

    entries += nCols * interaction._rowList.size();
  }

  return entries;
}

//////////////////////////////////////////////////////////////////////
// Helps to estimate the additional compressed storage costs this
// node will introduce, based on the given rank (estimate) for
// its low rank blocks.
//
// In particular, we need the size of an expansion workspace necessary
// to expand the set of all compressed interactions for any node
// in the factor.
//////////////////////////////////////////////////////////////////////
void Supernode::estimateCompressedExpansionStorage(
                            const std::vector<Supernode> &nodes,
                            int rank, IntArray &expansionStorage ) const
{
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    if ( interaction._active )
    {
      continue;
    }

    const SupernodeInteraction  &extendedInteraction
                          = _offDiagonal[ interaction._extendedInteraction ];

    const PairArray &forwardInteractions
                          = extendedInteraction._forwardInteractions;

    for ( int next_interaction_idx = 0;
          next_interaction_idx < forwardInteractions.size();
          next_interaction_idx++ )
    {
      const IndexPair &indices = forwardInteractions[ next_interaction_idx ];

      const Supernode &nextNode = nodes[ indices.first ];

      const SupernodeInteraction &nextInteraction
                                = nextNode._offDiagonal[ indices.second ];

      // If this node is compressed, add the cost of expanding it to
      // nextNode
      if ( nextInteraction._compressed )
      {
        // Estimated cost of expanding this interaction
        expansionStorage[ indices.first ] += rank * nextNode.numColumns();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Sets which blocks to keep in the main diagonal, and
// which to treat as low rank
//////////////////////////////////////////////////////////////////////
void Supernode::setDiagonalBlocks( const vector<DenseBlock> &diagonalBlocks,
                                   bool compressInPlace,
                                   int explicitBlockThreshold )
{
  long int                   offset = 0;

  if ( diagonalBlocks.size() <= 1 && !compressInPlace ) {
    return;
  } else if ( diagonalBlocks.size() == 0 && compressInPlace ) {
    _diagonalBlocks.push_back( DenseBlock( IndexRange( 0, numColumns() - 1 ),
                                           IndexRange( 0, numColumns() - 1 ) ) );
  } else {
    _diagonalBlocks = diagonalBlocks;
  }

  // Set offsets for each block
  _diagonalOffsets.resize( _diagonalBlocks.size() );
  for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ ) {
    int                      nCols = _diagonalBlocks[ block_idx ].numColumns();

    _diagonalOffsets[ block_idx ] = offset;

    offset += nCols * nCols;
  }

  // Build a map that allows us to map indices in this node back to
  // the appropriate diagonal block.
  buildDiagonalMap();

  setDiagonalSize();

  // Set up extended variable information for compressed blocks
  _mainDiagonalContributions.resize( _diagonalBlocks.size() );

  if ( !compressInPlace ) {
    // FIXME: print out diagonal blocks
    printf( "Node %d has %d diagonal blocks\n",
            _nodeID, (int)_diagonalBlocks.size() );

    buildLowRankDiagonalBlocks( 0, diagonalBlocks.size() );
  } else {
    // Initialize the compressed diagonal block
    _nodeDiagonal.init( _diagonalBlocks, _diagonalOffsets,
                        explicitBlockThreshold );

    // Initialize the list of low-rank blocks on the diagonal
    _diagonalLowRankBlocks.clear();
    for ( int block_idx = 0; block_idx < _nodeDiagonal.numOffDiagonalBlocks();
          block_idx++ )
    {
      _diagonalLowRankBlocks.push_back(
        _nodeDiagonal.offDiagonalBlock( block_idx ) );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Expands each compressed interaction stored in this node in
// to the given workspace
//////////////////////////////////////////////////////////////////////
void Supernode::expandCompressedInteractions(
                               const Real *extraData, Real *expansionWorkspace,
                               bool diagonalContributions,
                               int maxSz ) const
{
  const Real                *baseData;
  int                        nCols = numColumns();

  int                        numInteractions;
  int                        interaction_idx;

  // FIXME
  long int totalSz = 0;

  numInteractions = diagonalContributions ? _diagonalContributions.size()
                                          : _offDiagonal.size();

  for ( int i = 0; i < numInteractions; i++ )
  {
    interaction_idx = diagonalContributions ? _diagonalContributions[ i ] : i;

    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    // DSB: Is ! right?
    if ( !interaction._type == EXTENDED_NODE || !interaction._compressed )
    {
      continue;
    }

    baseData = extraData + interaction._extendedDataOffset;

    // FIXME
    totalSz += interaction._numExtendedRows * nCols;
    if ( maxSz > 0 && ( totalSz <= 0 || totalSz > maxSz ) )
    {
      cout << SDUMP( totalSz ) << "    " << SDUMP( maxSz ) << endl;
      TRACE_ASSERT( NULL, "Out of space in expansion workspace" );
    }

    MATRIX::clear( expansionWorkspace, interaction._numExtendedRows, nCols );
    MATRIX::scatterColumns( baseData, expansionWorkspace,
                            interaction._compressedColumnList,
                            interaction._numExtendedRows,
                            interaction._compressedColumnList.size(),
                            nCols );

    expansionWorkspace += interaction._numExtendedRows * nCols;
  }
}

//////////////////////////////////////////////////////////////////////
// For an interaction stored in place in compressed format V * U'
// where V has orthonormal columns, convert these interactions
// so that the U now has orthonormal columns (and V is arbitrary).
//////////////////////////////////////////////////////////////////////
void Supernode::convertCompressedOrthoInteraction( int interaction_idx )
{
  SupernodeInteraction      &interaction = _offDiagonal[ interaction_idx ];

  Real                      *V = interaction._compressedV;
  Real                      *U = interaction._compressedU;

  int                        nCols = numColumns();

  cout << "Converting compressed interaction" << endl;

  printf( "V (old) norm: %f    U (old) norm: %f\n",
          VECTOR::norm2( V, interaction._rowList.size()
                            * interaction._numExtendedRows ),
          VECTOR::norm2( U, nCols * interaction._numExtendedRows ) );
  
  TRACE_ASSERT( interaction._type == COMPRESSED_NODE,
                "Invalid interaction type" );

  // Some workspaces for conversion
  MATRIX                     Ucopy( nCols, interaction._numExtendedRows, U );
  MATRIX                     Vcopy( interaction._rowList.size(),
                                    interaction._numExtendedRows, V );
  MATRIX                     Rtrans( interaction._numExtendedRows,
                                     interaction._numExtendedRows );
  MATRIX                     qrData( min( Ucopy.rows(), Ucopy.cols() ), 1 );
                                            
  // Get a QR factorization for the U matrix
  MATRIX::qr( Ucopy.data(), Ucopy.rows(), Ucopy.cols(), qrData.data(),
              NULL, -1 );
  MATRIX::extractQRfactor( Ucopy.data(), Ucopy.rows(), Ucopy.cols(),
                           qrData.data(), NULL, -1 );

  // At this point, U = Ucopy * R, where R = Ucopy' * U.
  // We want to form V <-- V * R' = V * U' * Ucopy
  MATRIX::gemm( U, Ucopy.data(), Rtrans.data(),
                // Dimensions of U and Ucopy
                nCols, interaction._numExtendedRows,
                nCols, interaction._numExtendedRows,
                // Transpose U, but not Ucopy
                true, false );

  // Multiply by V to form the new (non-orthonormal) V
  MATRIX::gemm( Vcopy.data(), Rtrans.data(), V,
                // Dimensions of V and Rtrans
                Vcopy.rows(), Vcopy.cols(), Rtrans.rows(), Rtrans.cols(),
                // No transposition
                false, false );

  // Replace U
  MATRIX::copy( U, Ucopy.data(), Ucopy.rows(), Ucopy.cols() );

  printf( "V (new) norm: %f    U (new) norm: %f\n\n",
          VECTOR::norm2( V, interaction._rowList.size()
                            * interaction._numExtendedRows ),
          VECTOR::norm2( U, nCols * interaction._numExtendedRows ) );
}

//////////////////////////////////////////////////////////////////////
// Converts all interactions using the function above
//////////////////////////////////////////////////////////////////////
void Supernode::convertAllCompressedOrthoInteractions()
{
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    if ( _offDiagonal[ interaction_idx ]._type == COMPRESSED_NODE ) {
      convertCompressedOrthoInteraction( interaction_idx );
    }
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// For each interaction in this node which makes a contribution
// to the diagonal, form the given submatrx (specified by block),
// multiply it by the given matrix G and add to the workspace.
//
// expandCompressedInteractions is assumed to have been called
// previously, with results stored in expansionWorkspace.
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalContributionMultiply( const Real *extraData,
                                              const Real *expansionWorkspace,
                                              const Real *G, int nColsG,
                                              const DenseBlock &block,
                                              Real *multWorkspaceInitial,
                                              Real *multWorkspaceFinal )
{
  int                        interaction_idx;
  const Real                *baseData;
  int                        nCols = numColumns();

  int                        startRow, startCol;
  int                        nBlockRows, nBlockCols;

  // Get the row/column range for the block we're working with
  startRow = block._rowRange.first;
  startCol = block._columnRange.first;

  nBlockRows = range_size( block._rowRange );
  nBlockCols = range_size( block._columnRange );

  for ( int i = 0; i < _diagonalContributions.size(); i++ )
  {
    interaction_idx = _diagonalContributions[ i ];

    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    if ( interaction._compressed )
    {
      baseData = expansionWorkspace;
      expansionWorkspace += interaction._numExtendedRows * nCols;
    }
    else
    {
      baseData = extraData + interaction._extendedDataOffset;
    }

    // We only want a row subset of the matrix.  However, note
    // that the transpose of the matrix is stored, so we need
    // to align with a column to get the desired row.
    //
    // If V is the stored version of the matrix, we want to form the
    // product V(:, rowRange)' * V(:, colRange) * G
    MATRIX::gemm( baseData + startCol, /* Align with correct column */
                  G, multWorkspaceInitial,
                  interaction._numExtendedRows, /* Rows in V */
                  nBlockCols, /* Column size for the block */
                  nBlockCols, nColsG, /* Dimensions of G */
                  false, false, /* Do not transpose */
                  1.0, 0.0, /* Overwrite with result */
                  nCols /* Leading dimension for submatrix of V */ );

    MATRIX::gemm( baseData + startRow, /* Align with correct column */
                  multWorkspaceInitial, multWorkspaceFinal,
                  interaction._numExtendedRows, /* Rows in V */
                  nBlockRows, /* Row size for the block */
                  interaction._numExtendedRows, /* Rows in prior result */
                  nColsG, /* Columns in prior result */
                  true /* Transpose V submatrix */, false, /* No transpose */
                  1.0, 1.0, /* Add to result */
                  nCols /* Leading dimension for submatrix of V */ );
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Adds all diagonal contributions to this node resulting from the
// addition of slack variables.  We require both the extraData
// array, as well as an "expansion workspace" in which all compressed
// interactions have been expanded.
//
// That is, we expect expandCompressedInteractions to have been called
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContributions(
                                      const Real *extraData,
                                      Real *expansionWorkspace )
{
  addLowRankDiagonalContributions_offDiagonalContributions( extraData,
                                                            expansionWorkspace );

  addLowRankDiagonalContributions_diagonalContributions( extraData );
}

//////////////////////////////////////////////////////////////////////
// Adds all diagonal contributions to this node resulting from the
// addition of slack variables.  We require both the extraData
// array, as well as an "expansion workspace" in which all compressed
// interactions have been expanded.
//
// That is, we expect expandCompressedInteractions to have been called
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContributions(
                                const Real *extraData,
                                WorkspaceManager<Real> &workspaceManager )
{
  addLowRankDiagonalContributions_offDiagonalContributions( extraData,
                                                            workspaceManager );

  addLowRankDiagonalContributions_diagonalContributions( extraData );
}

#if 0
//////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContributions(
                      const SPARSE_MATRIX::SparseColumnMatrix &A,
                      const vector<Supernode> &nodes,
                      const IntArray &lowRankBlocks,
                      const vector<vector<IndexPair> > &lowRankDescendents,
                      const IntArray &ancestorInteractions,
                      Real *multWorkspace, Real *workspace )
{
  int                        lowrank_block_idx = 0;

  Real                       contributionScale;

  for ( int diagonal_contrib_idx = 0;
        diagonal_contrib_idx < _diagonalContributions.size();
        diagonal_contrib_idx++ )
  {
    const SupernodeInteraction &contribInteraction
          = _offDiagonal[ _diagonalContributions[ diagonal_contrib_idx ] ];

    contributionScale = contribInteraction._diagonalContributionScale;
    // Since we are squaring the matrix
    contributionScale *= contributionScale;

    TRACE_ASSERT( contributionScale > 0.0,
                  "Contribution scale has not been set" );

    // Two cases here:
    //    1) This contribution is the result of a sparsified block from
    //       earlier in the factorization, in which case we just add a
    //       scaled identity matrix to the diagonal.
    //
    //    2) This contribution is the result of some sparsified block S
    //       in this node, in which case we add S' * S to the diagonal
    //       (with some scaling factor)
    if ( !contribInteraction.hasParent() )
    {
      addLowRankDiagonalContribution( contributionScale );
    }
    else
    {
      // Find the low rank block producing this interaction
      TRACE_ASSERT(
                lowRankBlocks[ lowrank_block_idx ]
                  == contribInteraction._extendedInteraction,
                  "No corresponding low rank interaction found" );

      addLowRankDiagonalContribution( A, nodes,
                                      lowRankBlocks, lowrank_block_idx,
                                      lowRankDescendents, ancestorInteractions,
                                      contributionScale,
                                      multWorkspace, workspace );

      lowrank_block_idx++;
    }
  }

  TRACE_ASSERT( lowrank_block_idx == lowRankBlocks.size(),
                "Missed some low rank block" );
}
#endif

//////////////////////////////////////////////////////////////////////
// Adds this node's contribution to the extended variable schur
// complement.
//
// Need an expansion workspace to expand any compressed column
// interactions.
//
// multWorkspace is used to extract the needed columns (in the case
// of compressed column interactions) prior to multiplication
//////////////////////////////////////////////////////////////////////
void Supernode::addExtendedSchurComplementContribution(
                          const vector<Supernode> &nodes,
                          Real *expansionWorkspace, Real *multWorkspace,
                          Real *extraData,
                          Timer *invertTimer,
                          Timer *multiplyTimer )
{
  int                        nCols = numColumns();
  Real                      *baseData;

  int                        numMultiplyColumns;

  IndexRange                 fullRange( 0, nCols - 1 );

  // Iterate over every extended interaction, and push out necessary
  // Schur complement contributions
  for ( int interaction_idx = firstExtendedInteraction();
        interaction_idx < _offDiagonal.size(); interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // Subtract from
    const Supernode             &schurNode = nodes[ interaction._nodeID ];
    int                          schur_interaction_idx = 0;

    TRACE_ASSERT( schurNode.type() == EXTENDED_NODE );

    const IndexRange            &interactionRange
                                  = interaction._lowRankDiagonalInteraction ?
                                        interaction._compressedColumnRange
                                      : fullRange;

    TRACE_ASSERT( interaction._numExtendedRows == schurNode.numColumns() );

    // Apply the diagonal matrix inverse to this interaction and put
    // the result in expansionWorkspace
#ifdef DO_TIMING
    if ( invertTimer )
    {
      invertTimer->tick();
    }
#endif
    applyInverseToInteraction( interaction_idx, extraData, expansionWorkspace );
#ifdef DO_TIMING
    if ( invertTimer )
    {
      invertTimer->tock();
    }
#endif

    // Subtract from the diagonal of schurNode
#ifdef DO_TIMING
    if ( multiplyTimer )
    {
      multiplyTimer->tick();
    }
#endif
    multiplyInteraction( interaction_idx, expansionWorkspace,
                         interaction._numExtendedRows,
                         interactionRange,
                         extraData, multWorkspace,
                         extraData + schurNode._extendedDataOffset );
#ifdef DO_TIMING
    if ( multiplyTimer )
    {
      multiplyTimer->tock();
    }
#endif

    // Add to the remaining interactions in schurNode
    for ( int next_interaction_idx = interaction_idx + 1;
          next_interaction_idx < _offDiagonal.size(); next_interaction_idx++ )
    {
      const SupernodeInteraction 
                      &nextInteraction = _offDiagonal[ next_interaction_idx ];

      // If both interaction and nextInteraction are the result of compression
      // blocks in the diagonal, check for overlap in their column ranges
      // before propogating
      if ( interaction._lowRankDiagonalInteraction
        && nextInteraction._lowRankDiagonalInteraction
        && !range_overlap( interaction._compressedColumnRange,
                           nextInteraction._compressedColumnRange ) )
      {
        continue;
      }

      // Find the interaction we are filling in
      while ( schur_interaction_idx < schurNode._offDiagonal.size()
           && ( schurNode._offDiagonal[ schur_interaction_idx ]._nodeID
             != nextInteraction._nodeID ) )
      {
        schur_interaction_idx++;
      }

      TRACE_ASSERT( schur_interaction_idx < schurNode._offDiagonal.size(),
                    "Schur complement interaction not found!" );

      const SupernodeInteraction  &schurInteraction
                            = schurNode._offDiagonal[ schur_interaction_idx ];

      TRACE_ASSERT( schurInteraction._nodeID == nextInteraction._nodeID );

      // Subtract from the off-diagonal of schurNode
#ifdef DO_TIMING
      if ( multiplyTimer )
      {
        multiplyTimer->tick();
      }
#endif
      multiplyInteraction( next_interaction_idx, expansionWorkspace,
                           interaction._numExtendedRows,
                           interactionRange,
                           extraData, multWorkspace,
                           extraData + schurInteraction._extendedDataOffset );
#ifdef DO_TIMING
      if ( multiplyTimer )
      {
        multiplyTimer->tock();
      }
#endif
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Adds this node's contribution to the extended variable schur
// complement.
//
// Need an expansion workspace to expand any compressed column
// interactions.
//
// multWorkspace is used to extract the needed columns (in the case
// of compressed column interactions) prior to multiplication
//////////////////////////////////////////////////////////////////////
void Supernode::addExtendedSchurComplementContribution(
                          const std::vector<Supernode> &nodes,
                          WorkspaceManager<Real> &workspaceManager,
                          Real *extraData,
                          Timer *invertTimer,
                          Timer *multiplyTimer )
{
  int                        nCols = numColumns();
  Real                      *baseData;

  int                        numMultiplyColumns;

  IndexRange                 fullRange( 0, nCols - 1 );

  // Workspace stuff
  PairArray                  dataSizes( 2 );

  Real                      *expansionWorkspace;
  Real                      *multWorkspace;

  // Iterate over every extended interaction, and push out necessary
  // Schur complement contributions
  for ( int interaction_idx = firstExtendedInteraction();
        interaction_idx < _offDiagonal.size(); interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // Subtract from
    const Supernode             &schurNode = nodes[ interaction._nodeID ];
    int                          schur_interaction_idx = 0;

    TRACE_ASSERT( schurNode.type() == EXTENDED_NODE );

    const IndexRange            &interactionRange
                                  = interaction._lowRankDiagonalInteraction ?
                                        interaction._compressedColumnRange
                                      : fullRange;

    TRACE_ASSERT( interaction._numExtendedRows == schurNode.numColumns() );

    // Set workspace sizes
    if ( interaction._lowRankDiagonalInteraction )
    {
      dataSizes[ 0 ]
        = IndexPair( interaction._numExtendedRows,
                     range_size( interaction._compressedColumnRange ) );
    }
    else
    {
      dataSizes[ 0 ] = IndexPair( interaction._numExtendedRows, nCols );
    }
    dataSizes[ 1 ] = dataSizes[ 0 ];

    RealWorkspace            workspace( workspaceManager, dataSizes );

    expansionWorkspace = workspace.workspaceData( 0 );
    multWorkspace = workspace.workspaceData( 1 );

    // Apply the diagonal matrix inverse to this interaction and put
    // the result in expansionWorkspace
#ifdef DO_TIMING
    if ( invertTimer )
    {
      invertTimer->tick();
    }
#endif
    applyInverseToInteraction( interaction_idx, extraData, expansionWorkspace );
#ifdef DO_TIMING
    if ( invertTimer )
    {
      invertTimer->tock();
    }
#endif

    // Subtract from the diagonal of schurNode
#ifdef DO_TIMING
    if ( multiplyTimer )
    {
      multiplyTimer->tick();
    }
#endif
    multiplyInteraction( interaction_idx, expansionWorkspace,
                         interaction._numExtendedRows,
                         interactionRange,
                         extraData, multWorkspace,
                         extraData + schurNode._extendedDataOffset );
#ifdef DO_TIMING
    if ( multiplyTimer )
    {
      multiplyTimer->tock();
    }
#endif

    // Add to the remaining interactions in schurNode
    for ( int next_interaction_idx = interaction_idx + 1;
          next_interaction_idx < _offDiagonal.size(); next_interaction_idx++ )
    {
      const SupernodeInteraction 
                      &nextInteraction = _offDiagonal[ next_interaction_idx ];

      // If both interaction and nextInteraction are the result of compression
      // blocks in the diagonal, check for overlap in their column ranges
      // before propogating
      if ( interaction._lowRankDiagonalInteraction
        && nextInteraction._lowRankDiagonalInteraction
        && !range_overlap( interaction._compressedColumnRange,
                           nextInteraction._compressedColumnRange ) )
      {
        continue;
      }

      // Find the interaction we are filling in
      while ( schur_interaction_idx < schurNode._offDiagonal.size()
           && ( schurNode._offDiagonal[ schur_interaction_idx ]._nodeID
             != nextInteraction._nodeID ) )
      {
        schur_interaction_idx++;
      }

      TRACE_ASSERT( schur_interaction_idx < schurNode._offDiagonal.size(),
                    "Schur complement interaction not found!" );

      const SupernodeInteraction  &schurInteraction
                            = schurNode._offDiagonal[ schur_interaction_idx ];

      TRACE_ASSERT( schurInteraction._nodeID == nextInteraction._nodeID );

      // Subtract from the off-diagonal of schurNode
#ifdef DO_TIMING
      if ( multiplyTimer )
      {
        multiplyTimer->tick();
      }
#endif
      multiplyInteraction( next_interaction_idx, expansionWorkspace,
                           interaction._numExtendedRows,
                           interactionRange,
                           extraData, multWorkspace,
                           extraData + schurInteraction._extendedDataOffset );
#ifdef DO_TIMING
      if ( multiplyTimer )
      {
        multiplyTimer->tock();
      }
#endif
    }
  }
}

#define INTERIOR_SOLVE_TIMING 1
//#undef INTERIOR_SOLVE_TIMING

//////////////////////////////////////////////////////////////////////
// Runs a forward solve by solving for the components of
// x governed by this node's indices, then making appropriate
// updates to x to account for off-diagonal parts of this node
//////////////////////////////////////////////////////////////////////
void Supernode::forwardSolve( Real *x,
                              Real *workspace, Real *extendedWorkspace,
                              const vector<Supernode> &nodes,
                              Real *extraData,
                              Real *expansionWorkspace )
{
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_START( "Forward solve (other):" );
#endif

  int                            systemRow;
  int                            work_idx = 0;
  int                            extended_work_idx = 0;
  int                            nCols = numColumns();
  int                            nRows;
  Real                          *baseData;

  Real                          *diagonalData;

  bool                           foundCompressed = false;

  int                            startCol = _columnRange.first;

  diagonalData = ( _type == STANDARD_NODE ) ? _data
                                            : extraData + _extendedDataOffset;

#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Forward solve (other):" );
  TIMING_START( "Forward solve copy" );
#endif
  // Gather components of x corresponding to off-diagonal stuff
  // in to the workspace
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // If the interaction is with an active, standard node, copy data
    // to the work space row by row.
    //
    // If it is for an extended node, we can copy the full block of
    // data corresponding to this extended node.
    if ( (interaction._type == STANDARD_NODE && interaction._active)
      || (interaction._type == COMPRESSED_NODE && _compressedRank > 0) ) {
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Forward solve (other):" );
      TIMING_START( "Forward solve standard copy" );
#endif
#endif
      for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ ) {
        systemRow = nodes[ interaction._nodeID ]._columnRange.first
                  + interaction._rowList[ row_idx ];

        TRACE_ASSERT( systemRow >= 0, "Something is very wrong" );

        workspace[ work_idx ] = x[ systemRow ];

        work_idx++;
      }
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Forward solve standard copy" );
      TIMING_START( "Forward solve (other):" );
#endif
#endif
    } else if ( interaction._type == EXTENDED_NODE ) {
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Forward solve (other):" );
      TIMING_START( "Forward solve extended copy" );
#endif
#endif
      nRows = nodes[ interaction._nodeID ].numColumns();

      systemRow = nodes[ interaction._nodeID ]._columnRange.first;

      MATRIX::copy( extendedWorkspace + extended_work_idx,
                    x + systemRow, nRows, 1 );

      extended_work_idx += nRows;
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Forward solve extended copy" );
      TIMING_START( "Forward solve (other):" );
#endif
#endif
    }

    // FIXME: treat all standard node interactions as compressed
    //if ( interaction._compressed )
    if ( _type == STANDARD_NODE && interaction._type == EXTENDED_NODE ) {
      foundCompressed = true;
    }
  }
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Forward solve copy" );
  TIMING_START( "Forward solve (other):" );
#endif

  // FIXME: diagonal sparsification

  // Do the initial triangular solve
  //    Solve L1 * x1 = x1  (yes... bad notation)
#if 0
  MATRIX::triangularSolve( diagonalData, x + startCol, nCols );
#endif
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Forward solve (other):" );
  if ( _diagonalBlocks.size() > 1 ) {
    FLOP_COUNT_START( "Forward compressed diagonal solve" );
    TIMING_START( "Forward compressed diagonal solve" );
  } else {
    FLOP_COUNT_START( "Forward diagonal solve" );
    TIMING_START( "Forward diagonal solve" );
  }
#endif

#if 0
  diagonalSolve( x + startCol, 1, extraData );
#endif
  diagonalSolve( x + startCol, extraData );

#ifdef INTERIOR_SOLVE_TIMING
  if ( _diagonalBlocks.size() > 1 ) {
    TIMING_STOP( "Forward compressed diagonal solve" );
    FLOP_COUNT_END( "Forward compressed diagonal solve" );
  } else {
    TIMING_STOP( "Forward diagonal solve" );
    FLOP_COUNT_END( "Forward diagonal solve" );
  }
    TIMING_START( "Forward solve (other):" );
#endif

  // If we have any compressed nodes, we also need to form
  // the solution to L1^T * y = x1
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Forward solve (other):" );
  FLOP_COUNT_START( "Forward expansion diagonal solve" );
  TIMING_START( "Forward expansion diagonal solve" );
#endif
  if ( foundCompressed ) {
    MATRIX::copy( expansionWorkspace, x + startCol, nCols, 1 );
#if 0
    diagonalSolve( expansionWorkspace, 1, extraData, true /* transpose */ );
#endif
    diagonalSolve( expansionWorkspace, extraData, true /* transpose */ );
  }
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Forward expansion diagonal solve" );
  FLOP_COUNT_END( "Forward expansion diagonal solve" );
  TIMING_START( "Forward solve (other):" );
#endif

#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Forward solve (other):" );
  FLOP_COUNT_START( "Forward standard propagate" );
  TIMING_START( "Forward standard propagate" );
#endif
  if ( _numRows > 0 ) {
    // Accumulate changes in to the workspace
    //    workspace = workspace - L2 * x1
    MATRIX::gemv( offDiagonalData(), x + startCol, workspace,
                  _numRows, nCols, false,
                  -1.0, /* alpha */
                  1.0 /* beta */ );
  } else if ( _compressedRank > 0 ) {
    VECTOR                   multWorkspace( _compressedRank );

    MATRIX::gemv( _compressedU, x + startCol, multWorkspace.data(),
                  nCols, _compressedRank, true );
    MATRIX::gemv( _compressedV, multWorkspace.data(), workspace,
                  work_idx, _compressedRank, false,
                  -1.0, 1.0 /* subtract */ );
  }
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Forward standard propagate" );
  FLOP_COUNT_END( "Forward standard propagate" );
  TIMING_START( "Forward solve (other):" );
#endif

  // Accumulate changes to the extended workspace
  extended_work_idx = 0;

#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Forward solve (other):" );
  FLOP_COUNT_START( "Forward std.-ext. propagate" );
  TIMING_START( "Forward std.-ext. propagate" );
#endif
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    if ( interaction._type != EXTENDED_NODE ) {
      continue;
    }

    nRows = nodes[ interaction._nodeID ].numColumns();

    baseData = extraData + interaction._extendedDataOffset;

    if ( interaction._compressed ) {
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Forward solve (other):" );
      FLOP_COUNT_START( "Forward std.-ext. propagate (compressed)" );
      TIMING_START( "Forward std.-ext. propagate (compressed)" );
#endif
#endif

      // baseData stores the Schur complement in column-compressed form, rather
      // than the factor itself.  We need an additional triangular solve here
      // to do the multiplication.  The result of this triangular solve
      // is stored in expansionWorkspace
      MATRIX::compressedColumnMult( baseData, expansionWorkspace,
                                    extendedWorkspace + extended_work_idx,
                                    interaction._compressedColumnList,
                                    nRows, nCols,
                                    false, /* don't transpose */
                                    -1.0, /* alpha */
                                    1.0 /* beta */ );

#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Forward std.-ext. propagate (compressed)" );
      FLOP_COUNT_END( "Forward std.-ext. propagate (compressed)" );
      TIMING_START( "Forward solve (other):" );
#endif
#endif
    } // FIXME: treat all interactions in a standard node as compressed
    else if ( _type == STANDARD_NODE ) {
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Forward solve (other):" );
      FLOP_COUNT_START( "Forward std.-ext. propagate (uncompressed)" );
      TIMING_START( "Forward std.-ext. propagate (uncompressed)" );
#endif
#endif

      MATRIX::gemv( baseData, expansionWorkspace,
                    extendedWorkspace + extended_work_idx,
                    nRows, nCols, false,
                    -1.0, /* alpha */
                    1.0 /* beta */ );

#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Forward std.-ext. propagate (uncompressed)" );
      FLOP_COUNT_END( "Forward std.-ext. propagate (uncompressed)" );
      TIMING_START( "Forward solve (other):" );
#endif
#endif
    } else {
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Forward solve (other):" );
      FLOP_COUNT_START( "Forward ext.-ext. propagate" );
      TIMING_START( "Forward ext.-ext. propagate" );
#endif
#endif

      // Accumulate changes in the workspace
      //    workspace = workspace - L2 * x1
      MATRIX::gemv( baseData, x + startCol,
                    extendedWorkspace + extended_work_idx,
                    nRows, nCols, false,
                    -1.0, /* alpha */
                    1.0 /* beta */ );

#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Forward ext.-ext. propagate" );
      FLOP_COUNT_END( "Forward ext.-ext. propagate" );
      TIMING_START( "Forward solve (other):" );
#endif
#endif
    }

    extended_work_idx += nRows;
  }
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Forward std.-ext. propagate" );
  FLOP_COUNT_END( "Forward std.-ext. propagate" );
  TIMING_START( "Forward solve (other):" );
#endif

  // Scatter workspace back to x
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Forward solve (other):" );
  TIMING_START( "Forward scatter back" );
#endif
  work_idx = 0;
  extended_work_idx = 0;
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    if ( (interaction._type == STANDARD_NODE && interaction._active)
         || (interaction._type == COMPRESSED_NODE && _compressedRank > 0) ) {
      for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ ) {
        systemRow = nodes[ interaction._nodeID ]._columnRange.first
                  + interaction._rowList[ row_idx ];

        x[ systemRow ] = workspace[ work_idx ];

        work_idx++;
      }
    } else if ( interaction._type == EXTENDED_NODE ) {
      nRows = nodes[ interaction._nodeID ].numColumns();

      systemRow = nodes[ interaction._nodeID ]._columnRange.first;

      MATRIX::copy( x + systemRow, extendedWorkspace + extended_work_idx,
                    nRows, 1 );

      extended_work_idx += nRows;
    } else if ( interaction._type == COMPRESSED_NODE ) {
      // Directly apply updates due to compressed interactions here
      const Real            *U = interaction._compressedU;
      const Real            *V = interaction._compressedV;

      // Workspaces for vector multipliction
      VECTOR                 workVector1( interaction._numExtendedRows );
      VECTOR                 workVector2( interaction._rowList.size() );

      MATRIX::gemv( U, x + startCol, workVector1.data(),
                    // Dimensions of U
                    nCols, interaction._numExtendedRows,
                    // Transpose
                    true );

      MATRIX::gemv( V, workVector1.data(), workVector2.data(),
                    // Dimensions of V
                    interaction._rowList.size(), interaction._numExtendedRows,
                    // Don't transpose
                    false );

      // Subtract from the output vector
      for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ ) {
        systemRow = nodes[ interaction._nodeID ]._columnRange.first;
        systemRow += interaction._rowList[ row_idx ];

        x[ systemRow ] -= workVector2( row_idx );
      }
    }
  }
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Forward scatter back" );
#endif
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but solves L' x = b instead
//////////////////////////////////////////////////////////////////////
void Supernode::backwardSolve( Real *x,
                               Real *workspace, Real *extendedWorkspace,
                               const vector<Supernode> &nodes,
                               Real *extraData,
                               Real *expansionWorkspace )
{
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_START( "Backward solve (other)" );
#endif

  int                            systemRow;
  int                            work_idx = 0;
  int                            extended_work_idx = 0;
  int                            nCols = numColumns();
  int                            nRows;
  Real                          *baseData;

  Real                          *diagonalData;

  bool                           foundCompressed = false;
  bool                           foundExtended = false;

  int                            startCol = _columnRange.first;

  diagonalData = ( _type == STANDARD_NODE ) ? _data
                                            : extraData + _extendedDataOffset;

  // Gather components of x corresponding to off-diagonal stuff
  // in to the workspace
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Backward solve (other)" );
  TIMING_START( "Backward solve copy" );
#endif
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    if ( (interaction._type == STANDARD_NODE && interaction._active )
         || (interaction._type == COMPRESSED_NODE && _compressedRank > 0) ) {
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_START( "Backward solve standard copy" );
#endif
#endif
      for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ ) {
        systemRow = nodes[ interaction._nodeID ]._columnRange.first
                  + interaction._rowList[ row_idx ];

        TRACE_ASSERT( systemRow >= 0, "Something is very wrong" );

        workspace[ work_idx ] = x[ systemRow ];

        work_idx++;
      }
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Backward solve standard copy" );
#endif
#endif
    } else if ( interaction._type == EXTENDED_NODE ) {
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_START( "Backward solve extended copy" );
#endif
#endif
      nRows = nodes[ interaction._nodeID ].numColumns();

      systemRow = nodes[ interaction._nodeID ]._columnRange.first;

      MATRIX::copy( extendedWorkspace + extended_work_idx,
                    x + systemRow, nRows, 1 );

      extended_work_idx += nRows;

      foundExtended = true;
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      TIMING_STOP( "Backward solve extended copy" );
#endif
#endif
    } else if ( interaction._type == COMPRESSED_NODE ) {
      // Directly apply updates from compressed interactions here
      const Real            *U = interaction._compressedU;
      const Real            *V = interaction._compressedV;

      // Workspaces for vector multipliction
      VECTOR                 workVector1( interaction._rowList.size() );
      VECTOR                 workVector2( interaction._numExtendedRows );

      // Copy the necessary rows from the input vector
      for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ )
      {
        systemRow = nodes[ interaction._nodeID ]._columnRange.first;
        systemRow += interaction._rowList[ row_idx ];

        workVector1( row_idx ) = x[ systemRow ];
      }

      MATRIX::gemv( V, workVector1.data(), workVector2.data(),
                    // Dimensions of V
                    interaction._rowList.size(), interaction._numExtendedRows,
                    // Transpose
                    true );

      // Subtract from the input vector
      MATRIX::gemv( U, workVector2.data(), x + startCol,
                    // Dimensions of U
                    nCols, interaction._numExtendedRows,
                    // Don't transpose
                    false,
                    // Subtract
                    -1.0, 1.0 );
    }

    if ( interaction._compressed ) {
      foundCompressed = true;
    }
  }
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Backward solve copy" );
#endif

#ifdef INTERIOR_SOLVE_TIMING
  FLOP_COUNT_START( "Backward standard propagate" );
  TIMING_START( "Backward standard propagate" );
#endif
  if ( _numRows > 0 ) {
    // Accumulate previous changes to x
    //    x1 = x1 - L2' * workspace
    MATRIX::gemv( offDiagonalData(), workspace, x + startCol,
                  _numRows, nCols, true, /* transpose */
                  -1.0, /* alpha */
                  1.0 /* beta */ );
  } else if ( _compressedRank > 0 ) {
    VECTOR                   multWorkspace( _compressedRank );

    MATRIX::gemv( _compressedV, workspace, multWorkspace.data(),
                  work_idx, _compressedRank, true );
    MATRIX::gemv( _compressedU, multWorkspace.data(), x + startCol,
                  nCols, _compressedRank, false,
                  -1.0, 1.0 /* subtract */ );
  }
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Backward standard propagate" );
  FLOP_COUNT_END( "Backward standard propagate" );
  TIMING_START( "Backward solve (other)" );
#endif

  // Accumulate changes in the extended workspace
  extended_work_idx = 0;

  if ( foundExtended && _type == STANDARD_NODE )
  {
    MATRIX::clear( expansionWorkspace, nCols, 1 );
  }

#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Backward solve (other)" );
  FLOP_COUNT_START( "Backward solve std.-ext propagate" );
  TIMING_START( "Backward solve std.-ext propagate" );
#endif
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    if ( interaction._type != EXTENDED_NODE ) {
      continue;
    }

    nRows = nodes[ interaction._nodeID ].numColumns();

    baseData = extraData + interaction._extendedDataOffset;

    if ( interaction._compressed )
    {
      // baseData stores the Schur complement in column-compressed form, rather
      // than the factor itself.  We need an additional triangular solve here
      // to do the multiplication.
#if 0
      MATRIX::compressedColumnMult( baseData,
                                    extendedWorkspace + extended_work_idx,
                                    expansionWorkspace,
                                    interaction._compressedColumnList,
                                    nRows, nCols,
                                    true /* transpose */ );

      diagonalSolve( expansionWorkspace, 1, extraData,
                     false /* don't transpose */ );

      // Subtract from x + startCol
      MATRIX::axpy( x + startCol, expansionWorkspace, nCols, 1, -1.0 );
#endif
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      FLOP_COUNT_START( "Backward std.-ext. propagate (compressed)" );
      TIMING_START( "Backward std.-ext. propagate (compressed)" );
#endif
#endif

      MATRIX::compressedColumnMult( baseData,
                                    extendedWorkspace + extended_work_idx,
                                    expansionWorkspace,
                                    interaction._compressedColumnList,
                                    nRows, nCols,
                                    true, /* transpose */
                                    1.0, 1.0 /* add */ );

#if 0
#ifdef INTERIOR_SOLVE_TIMING
      FLOP_COUNT_END( "Backward std.-ext. propagate (compressed)" );
      TIMING_STOP( "Backward std.-ext. propagate (compressed)" );
#endif
#endif
    }
    // FIXME: treat all interactions as compressed in a standard node for now
    else if ( _type == STANDARD_NODE )
    {
#if 0
      MATRIX::gemv( baseData, extendedWorkspace + extended_work_idx,
                    expansionWorkspace,
                    nRows, nCols, true, /* transpose */
                    1.0, /* alpha */
                    0.0 /* beta */ );

      diagonalSolve( expansionWorkspace, 1, extraData,
                     false /* don't transpose */ );

      // Subtract from x + startCol
      MATRIX::axpy( x + startCol, expansionWorkspace, nCols, 1, -1.0 );
#endif
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      FLOP_COUNT_START( "Backward std.-ext. propagate (uncompressed)" );
      TIMING_START( "Backward std.-ext. propagate (uncompressed)" );
#endif
#endif

      MATRIX::gemv( baseData, extendedWorkspace + extended_work_idx,
                    expansionWorkspace,
                    nRows, nCols, true, /* transpose */
                    1.0, 1.0 /* add */ );

#if 0
#ifdef INTERIOR_SOLVE_TIMING
      FLOP_COUNT_END( "Backward std.-ext. propagate (uncompressed)" );
      TIMING_STOP( "Backward std.-ext. propagate (uncompressed)" );
#endif
#endif
    }
    else
    {
#if 0
#ifdef INTERIOR_SOLVE_TIMING
      FLOP_COUNT_START( "Backward ext.-ext. propagate" );
      TIMING_START( "Backward ext.-ext. propagate" );
#endif
#endif

      // Accumulate changes in the workspace
      //    workspace = workspace - L2 * x1
      MATRIX::gemv( baseData, extendedWorkspace + extended_work_idx,
                    x + startCol,
                    nRows, nCols, true, /* transpose */
                    -1.0, /* alpha */
                    1.0 /* beta */ );

#if 0
#ifdef INTERIOR_SOLVE_TIMING
      FLOP_COUNT_END( "Backward ext.-ext. propagate" );
      TIMING_STOP( "Backward ext.-ext. propagate" );
#endif
#endif
    }

    extended_work_idx += nRows;
  }
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Backward solve std.-ext propagate" );
  FLOP_COUNT_END( "Backward solve std.-ext propagate" );
#endif

  // Solve for the compressed part
#ifdef INTERIOR_SOLVE_TIMING
  FLOP_COUNT_START( "Backward expansion diagonal solve" );
  TIMING_START( "Backward expansion diagonal solve" );
#endif
  if ( foundExtended && _type == STANDARD_NODE ) {
#if 0
    diagonalSolve( expansionWorkspace, 1, extraData, false /* no transpose */ );
#endif
    diagonalSolve( expansionWorkspace, extraData, false /* no transpose */ );

    // Subtract
    MATRIX::axpy( x + startCol, expansionWorkspace, nCols, 1, -1.0 );
  }
#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Backward expansion diagonal solve" );
  FLOP_COUNT_END( "Backward expansion diagonal solve" );
#endif

  // Do the final triangular solve
  //    Solve L1' x1 = x1
#if 0
  MATRIX::triangularSolve( diagonalData, x + startCol, nCols,
                           true, /* Store lower triangular part */
                           true /* Transpose */ );
#endif
#ifdef INTERIOR_SOLVE_TIMING
  FLOP_COUNT_START( "Backward diagonal solve" );
  TIMING_START( "Backward diagonal solve" );
#endif

#if 0
  diagonalSolve( x + startCol, 1, extraData, true /* transpose */ );
#endif
  diagonalSolve( x + startCol, extraData, true /* transpose */ );

#ifdef INTERIOR_SOLVE_TIMING
  TIMING_STOP( "Backward diagonal solve" );
  FLOP_COUNT_END( "Backward diagonal solve" );
#endif
}

//////////////////////////////////////////////////////////////////////
// Applies inverse of the LDL factorization block-diagonal, assuming
// that this factorization exists.  Does nothing for Cholesky-factored
// nodes.
//////////////////////////////////////////////////////////////////////
void Supernode::LDLdiagonalSolve( Real *x )
{
  if ( !_useLDL )
    return;

  int                            nCols = numColumns();
  int                            startCol = _columnRange.first;

  // Move to the part of x that we care about
  x = x + startCol;

  applyLDLDiagonalInverse( x, nCols, 1 );
}

//////////////////////////////////////////////////////////////////////
// Applies the permutation matrix from this node's LDL factorization,
// if it exists.  Does nothing for Chlolesky-factored nodes.
//////////////////////////////////////////////////////////////////////
void Supernode::applyLDLforwardPermutation( Real *x )
{
  if ( !_useLDL )
    return;

  int                            nCols = numColumns();
  int                            startCol = _columnRange.first;

  // Move to the part of x that we care about
  x = x + startCol;

  // Left side application, forward (non-transposed) permutation
  MATRIX::applyLDLPermutation( _LDL_data._pivotData, x, nCols, 1 );
}

//////////////////////////////////////////////////////////////////////
// Applies the inverse of the permutation matrix from this node's LDL
// factorization, if it exists.  Does notthing for Cholesky-factored
// nodes.
//////////////////////////////////////////////////////////////////////
void Supernode::applyLDLbackwardPermutation( Real *x )
{
  if ( !_useLDL )
    return;

  int                            nCols = numColumns();
  int                            startCol = _columnRange.first;

  // Move to the part of x that we care about
  x = x + startCol;

  // Left side application, backward (transposed) permutation
  MATRIX::applyLDLPermutation( _LDL_data._pivotData, x, nCols, 1,
                               true, true );
}

//////////////////////////////////////////////////////////////////////
// Solves for the components of x governed by this node's indices
// and makes appropriate updates to x.  This performs a "sub-solve"
// over only the specified set of nodes.
//
// This is used to support the "interior block" machinery.
//////////////////////////////////////////////////////////////////////
void Supernode::forwardSubSolve( Real *X, int nRHS,
                                 Real *workspace,
                                 const IntArray &nodeList,
                                 const std::vector<Supernode> &nodes,
                                 Timer *solveTimer, Timer *multTimer )
{
  TIMING_START( "forwardSubSolve (nodeList): preamble" );

#if 0
  int                        systemRow = _firstOffDiagonalRow;
#endif
  int                        startRow = 0;
  int                        currentRow;
  int                        block_node_idx = 0;
  int                        nCols = numColumns();

  Real                      *baseData;
  Real                      *baseOffDiagonalData;

  TRACE_ASSERT( _type == STANDARD_NODE );

  // Find this node's position in the sub-list of nodes
  for ( ; block_node_idx < nodeList.size(); block_node_idx++ ) {
    if ( nodeList[ block_node_idx ] == _nodeID ) {
      // We've found this node's position in the list
      break;
    }

    // Navigate to the correct starting column for this node
    startRow += nodes[ nodeList[ block_node_idx ] ].numColumns();
  }

  TRACE_ASSERT( block_node_idx < nodeList.size(), "Node index not found" );
  TRACE_ASSERT( startRow >= 0 );

  // Navigate to the row of X corresponding to this node's indices
  // and do a diagonal solve
  baseData = X + startRow * nRHS;

  TIMING_STOP( "forwardSubSolve (nodeList): preamble" );

  TIMING_START( "forwardSubSolve (nodeList): diagonalSolve" );
#ifdef DO_TIMING
  if ( solveTimer ) {
    solveTimer->tick();
  }
#endif
  diagonalSolve( baseData, nRHS, NULL /* No extra data needed */ );
#ifdef DO_TIMING
  if ( solveTimer ) {
    solveTimer->tock();
  }
#endif
  TIMING_STOP( "forwardSubSolve (nodeList): diagonalSolve" );

  TIMING_START( "forwardSubSolve (nodeList): propagate" );

  block_node_idx += 1;
  currentRow = startRow + nCols;

  // Gather components from the right hand side corresponding
  // to off-diagonal parts in to the workspace
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // Skip until we're at the next node index in our sub-list
    while ( block_node_idx < nodeList.size()
         && nodeList[ block_node_idx ] < interaction._nodeID )
    {
      // Skip to the next matrix block
      currentRow += nodes.at( nodeList[ block_node_idx ] ).numColumns();

      block_node_idx++;
    }

    if ( block_node_idx >= nodeList.size() ) {
      // We've encountered all of the relevant interactions
      break;
    }

    // Skip until we're at the next node index in our sub-list
    if ( interaction._nodeID != nodeList[ block_node_idx ] ) {
      TRACE_ASSERT( nodeList[ block_node_idx ] > interaction._nodeID );
      continue;
    }

    TRACE_ASSERT( interaction._type == STANDARD_NODE && interaction._active  );

#if 0
    for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ ) {
      systemRow = currentRow + interaction._rowList[ row_idx ];

      TRACE_ASSERT( systemRow >= 0, "Something is very wrong" );


    }
#endif
#ifdef DO_TIMING
    if ( multTimer ) {
      multTimer->tick();
    }
#endif

    // Copy the relevant row data to the workspace row by row
    MATRIX::copyRows( X, workspace, interaction._rowList, nRHS,
                      // Copy all rows for this interaction
                      IndexRange( 0, interaction._rowList.size() - 1 ),
                      // Negative offset, since we want to add
                      // currentRow to every row we take from X
                      -1 * currentRow );

    // Multiply by the solution for this block, and subtract from
    // the workspace
    baseOffDiagonalData = interactionMatrix( interaction_idx );
    MATRIX::gemm( baseOffDiagonalData, baseData, workspace,
                  interaction._rowList.size(), /* Number of rows */
                  nCols,
                  // RHS dimensions
                  nCols, nRHS,
                  false, false, /* No transpose */
                  -1.0, /* alpha */
                  1.0 /* beta */ );

    // Scatter the workspace back to the relevant part of X
    MATRIX::scatterRows( workspace, X, interaction._rowList, nRHS,
                         // Scatter all the rows
                         IndexRange( 0, interaction._rowList.size() - 1 ),
                         // Negative offset, since we want to add
                         // currentRow to every row we are pushing back
                         // to X
                         -1 * currentRow );

#ifdef DO_TIMING
    if ( multTimer ) {
      multTimer->tock();
    }
#endif
  }

  TIMING_STOP( "forwardSubSolve (nodeList): propagate" );
}

//////////////////////////////////////////////////////////////////////
// Version with a node range, rather than a list
//////////////////////////////////////////////////////////////////////
void Supernode::forwardSubSolve( Real *X, int nRHS,
                                 Real *workspace,
                                 IndexRange nodeRange,
                                 const std::vector<Supernode> &nodes )
{
  int                        work_idx = 0;
  int                        startRow = 0;
  int                        currentRow;
  int                        systemRow;
  int                        block_node_idx = 0;
  int                        nCols = numColumns();

  Real                      *baseData;
  Real                      *baseOffDiagonalData;

  TRACE_ASSERT( _type == STANDARD_NODE );
  TRACE_ASSERT( _nodeID >= nodeRange.first && _nodeID <= nodeRange.second );

  // Find a start row
  startRow = startColumn();
  startRow -= nodes[ nodeRange.first ].startColumn();

  // Navigate to the row of X corresponding to this node's indices
  // and do a diagonal solve
  baseData = X + startRow * nRHS;

  FLOP_COUNT_START( "forwardSubSolve (nodeRange): diagonal solve" );
  TIMING_START( "forwardSubSolve (nodeRange): diagonal solve" );
#if 0
  diagonalSolve( baseData, nRHS, NULL /* No extra data needed */ );
#endif
  diagonalSolve( baseData, NULL /* No extra data needed */ );
#if 0
  cblas_dtrsv( CblasRowMajor, CblasLower, CblasNoTrans,
               CblasNonUnit, nCols, _data, nCols, baseData, 1 );
#endif
#if 0
  cblas_dtrsm( CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans,
               CblasNonUnit, nCols, 1, 1.0, _data, nCols, baseData, 1 );
#endif
  TIMING_STOP( "forwardSubSolve (nodeRange): diagonal solve" );
  FLOP_COUNT_END( "forwardSubSolve (nodeRange): diagonal solve" );

  startRow = nodes[ nodeRange.first ].startColumn();

#if 0
  // Gather components of x corresponding to off-diagonal stuff in the
  // workspace
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // If we've reached the end of the node range, we're done
    if ( interaction._nodeID > nodeRange.second ) {
      break;
    }

    for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ ) {
      systemRow = nodes[ interaction._nodeID ]._columnRange.first
                + interaction._rowList[ row_idx ] - startRow;

      TRACE_ASSERT( systemRow >= 0, "Something is very wrong" );

      workspace[ work_idx ] = X[ systemRow ];

      work_idx++;
    }
  }

  FLOP_COUNT_START( "forwardSubSolve (nodeRange): propagate" );
  TIMING_START( "forwardSubSolve (nodeRange): propagate" );
  if ( work_idx > 0 ) {
    // Accumulate changes in to the workspace
    MATRIX::gemv( offDiagonalData(), baseData, workspace,
                  // Num. rows to consider == work_idx
                  work_idx, nCols, false,
                  -1.0, /* alpha */
                  1.0 /* beta */ );
  }
  TIMING_STOP( "forwardSubSolve (nodeRange): propagate" );
  FLOP_COUNT_END( "forwardSubSolve (nodeRange): propagate" );

  // Copy back to the solution vector
  work_idx = 0;
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // If we've reached the end of the node range, we're done
    if ( interaction._nodeID > nodeRange.second ) {
      break;
    }

    for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ ) {
      systemRow = nodes[ interaction._nodeID ]._columnRange.first
                + interaction._rowList[ row_idx ] - startRow;

      TRACE_ASSERT( systemRow >= 0, "Something is very wrong" );

      X[ systemRow ] = workspace[ work_idx ];

      work_idx++;
    }
  }
#endif

  // Count rows
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // If we've reached the end of the node range, we're done
    if ( interaction._nodeID > nodeRange.second ) {
      break;
    }

    work_idx += interaction._rowList.size();
  }

  FLOP_COUNT_START( "forwardSubSolve (nodeRange): propagate" );
  TIMING_START( "forwardSubSolve (nodeRange): propagate" );
  if ( work_idx > 0 ) {
    // Accumulate changes in to the workspace
    MATRIX::gemv( offDiagonalData(), baseData, workspace,
                  // Num. rows to consider == work_idx
                  work_idx, nCols, false,
                  -1.0, /* alpha */
                  0.0 /* beta */ );
  }
  TIMING_STOP( "forwardSubSolve (nodeRange): propagate" );
  FLOP_COUNT_END( "forwardSubSolve (nodeRange): propagate" );

  // Add to the solution vector
  work_idx = 0;
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // If we've reached the end of the node range, we're done
    if ( interaction._nodeID > nodeRange.second ) {
      break;
    }

    for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ ) {
      systemRow = nodes[ interaction._nodeID ]._columnRange.first
                + interaction._rowList[ row_idx ] - startRow;

      TRACE_ASSERT( systemRow >= 0, "Something is very wrong" );

      X[ systemRow ] += workspace[ work_idx ];

      work_idx++;
    }
  }

#if 0
  // Gather components from the right hand side corresponding
  // to off-diagonal parts in to the workspace
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // If we've reached the end of the node range, we're done
    if ( interaction._nodeID > nodeRange.second ) {
      break;
    }

    TRACE_ASSERT( interaction._type == STANDARD_NODE && interaction._active  );

    // Figure out the base row for this interaction in the
    // right hand side
    currentRow = nodes[ interaction._nodeID ].startColumn();
    currentRow -= nodes[ nodeRange.first ] .startColumn();

    // Copy the relevant row data to the workspace row by row
    MATRIX::copyRows( X, workspace, interaction._rowList, nRHS,
                      // Copy all rows for this interaction
                      IndexRange( 0, interaction._rowList.size() - 1 ),
                      // Negative offset, since we want to add
                      // currentRow to every row we take from X
                      -1 * currentRow );

    // Multiply by the solution for this block, and subtract from
    // the workspace
    baseOffDiagonalData = interactionMatrix( interaction_idx );
    MATRIX::gemm( baseOffDiagonalData, baseData, workspace,
                  interaction._rowList.size(), /* Number of rows */
                  nCols,
                  // RHS dimensions
                  nCols, nRHS,
                  false, false, /* No transpose */
                  -1.0, /* alpha */
                  1.0 /* beta */ );

    // Scatter the workspace back to the relevant part of X
    MATRIX::scatterRows( workspace, X, interaction._rowList, nRHS,
                         // Scatter all the rows
                         IndexRange( 0, interaction._rowList.size() - 1 ),
                         // Negative offset, since we want to add
                         // currentRow to every row we are pushing back
                         // to X
                         -1 * currentRow );
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but for backward solves
//////////////////////////////////////////////////////////////////////
void Supernode::backwardSubSolve( Real *X, int nRHS,
                                  Real *workspace,
                                  const IntArray &nodeList,
                                  const std::vector<Supernode> &nodes,
                                  Timer *solveTimer, Timer *multTimer )
{
  int                        startRow = 0;
  int                        currentRow;
  int                        block_node_idx = 0;
  int                        nCols = numColumns();

  Real                      *baseData;
  Real                      *baseOffDiagonalData;

  TRACE_ASSERT( _type == STANDARD_NODE );

  // Find this node's position in the sub-list of nodes
  for ( ; block_node_idx < nodeList.size(); block_node_idx++ ) {
    if ( nodeList[ block_node_idx ] == _nodeID ) {
      // We've found this node's position in the list
      break;
    }

    // Navigate to the correct starting column for this node
    startRow += nodes[ nodeList[ block_node_idx ] ].numColumns();
  }

  TRACE_ASSERT( block_node_idx < nodeList.size(), "Node index not found" );
  TRACE_ASSERT( startRow >= 0 );

  // Navigate to the row of X corresponding to this node's indices
  baseData = X + startRow * nRHS;

  block_node_idx += 1;
  currentRow = startRow + nCols;

  // Gather components from the right hand side corresponding
  // to off-diagonal parts in to the workspace
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // Skip until we're at the next node index in our sub-list
    while ( block_node_idx < nodeList.size()
         && nodeList[ block_node_idx ] < interaction._nodeID )
    {
      // Skip to the next matrix block
      currentRow += nodes.at( nodeList[ block_node_idx ] ).numColumns();

      block_node_idx++;
    }

    if ( block_node_idx >= nodeList.size() ) {
      // We've encountered all of the relevant interactions
      break;
    }

    // Skip until we're at the next node index in our sub-list
    if ( interaction._nodeID != nodeList[ block_node_idx ] ) {
      TRACE_ASSERT( nodeList[ block_node_idx ] > interaction._nodeID );
      continue;
    }

    TRACE_ASSERT( interaction._type == STANDARD_NODE && interaction._active  );

#ifdef DO_TIMING
    if ( multTimer ) {
      multTimer->tick();
    }
#endif

    // Copy the relevant row data to the workspace row by row
    MATRIX::copyRows( X, workspace, interaction._rowList, nRHS,
                      // Copy all rows for this interaction
                      IndexRange( 0, interaction._rowList.size() - 1 ),
                      // Negative offset, since we want to add
                      // currentRow to every row we take from X
                      -1 * currentRow );

    // Accumulate previous changes to x
    //      x1 = x1 - L2' * workspace
    baseOffDiagonalData = interactionMatrix( interaction_idx );
    MATRIX::gemm( baseOffDiagonalData, workspace, baseData,
                  interaction._rowList.size(), /* Number of rows */
                  nCols,
                  // Dimensions of workspace
                  interaction._rowList.size(), nRHS,
                  true, false, /* Transpose off-diagonal part */
                  -1.0, /* alpha */
                  1.0 /* beta */ );

#ifdef DO_TIMING
    if ( multTimer ) {
      multTimer->tock();
    }
#endif
  }

#ifdef DO_TIMING
  if ( solveTimer ) {
    solveTimer->tick();
  }
#endif

  // Perform the final triangular solve
  diagonalSolve( baseData, nRHS, NULL /* No extra data */,
                 true /* transpose */ );

#ifdef DO_TIMING
  if ( solveTimer ) {
    solveTimer->tock();
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Version with a node range, rather than a node list
//////////////////////////////////////////////////////////////////////
void Supernode::backwardSubSolve( Real *X, int nRHS,
                                  Real *workspace,
                                  IndexRange nodeRange,
                                  const std::vector<Supernode> &nodes )
{
  int                        startRow = 0;
  int                        currentRow;
  int                        block_node_idx = 0;
  int                        nCols = numColumns();
  int                        work_idx = 0;
  int                        systemRow;

  Real                      *baseData;
  Real                      *baseOffDiagonalData;

  TRACE_ASSERT( _type == STANDARD_NODE );
  TRACE_ASSERT( _nodeID >= nodeRange.first && _nodeID <= nodeRange.second );

#if 0
  // Find this node's position in the sub-list of nodes
  for ( ; block_node_idx < nodeList.size(); block_node_idx++ )
  {
    if ( nodeList[ block_node_idx ] == _nodeID )
    {
      // We've found this node's position in the list
      break;
    }

    // Navigate to the correct starting column for this node
    startRow += nodes[ nodeList[ block_node_idx ] ].numColumns();
  }

  TRACE_ASSERT( block_node_idx < nodeList.size(), "Node index not found" );
  TRACE_ASSERT( startRow >= 0 );
#endif

  // Find a start row
  startRow = startColumn();
  startRow -= nodes[ nodeRange.first ].startColumn();

  // Navigate to the row of X corresponding to this node's indices
  baseData = X + startRow * nRHS;

#if 0
  block_node_idx += 1;
  currentRow = startRow + nCols;
#endif

#if 0
  // Gather components from the right hand side corresponding
  // to off-diagonal parts in to the workspace
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

#if 0
    // Skip until we're at the next node index in our sub-list
    while ( block_node_idx < nodeList.size()
         && nodeList[ block_node_idx ] < interaction._nodeID )
    {
      // Skip to the next matrix block
      currentRow += nodes.at( nodeList[ block_node_idx ] ).numColumns();

      block_node_idx++;
    }

    if ( block_node_idx >= nodeList.size() )
    {
      // We've encountered all of the relevant interactions
      break;
    }

    // Skip until we're at the next node index in our sub-list
    if ( interaction._nodeID != nodeList[ block_node_idx ] )
    {
      TRACE_ASSERT( nodeList[ block_node_idx ] > interaction._nodeID );
      continue;
    }
#endif

    // If we've passed the end of the node range, we're done
    if ( interaction._nodeID > nodeRange.second ) {
      break;
    }

    TRACE_ASSERT( interaction._type == STANDARD_NODE && interaction._active  );

    // Figure out the base row for this interaction in the
    // right hand side
    currentRow = nodes[ interaction._nodeID ].startColumn();
    currentRow -= nodes[ nodeRange.first ] .startColumn();

    // Copy the relevant row data to the workspace row by row
    MATRIX::copyRows( X, workspace, interaction._rowList, nRHS,
                      // Copy all rows for this interaction
                      IndexRange( 0, interaction._rowList.size() - 1 ),
                      // Negative offset, since we want to add
                      // currentRow to every row we take from X
                      -1 * currentRow );

    // Accumulate previous changes to x
    //      x1 = x1 - L2' * workspace
    baseOffDiagonalData = interactionMatrix( interaction_idx );
#if 0
    MATRIX::gemm( baseOffDiagonalData, workspace, baseData,
                  interaction._rowList.size(), /* Number of rows */
                  nCols,
                  // Dimensions of workspace
                  interaction._rowList.size(), nRHS,
                  true, false, /* Transpose off-diagonal part */
                  -1.0, /* alpha */
                  1.0 /* beta */ );
#endif
    MATRIX::gemv( baseOffDiagonalData, workspace, baseData,
                  interaction._rowList.size(), /* Number of rows */
                  nCols,
                  true, /* Transpose off-diagonal part */
                  -1.0, /* alpha */
                  1.0 /* beta */ );
  }
#endif

  startRow = nodes[ nodeRange.first ].startColumn();

  // Gather components of x in to a workspace for multiplication
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // If we've reached the end of the node range, we're done
    if ( interaction._nodeID > nodeRange.second ) {
      break;
    }

    for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ ) {
      systemRow = nodes[ interaction._nodeID ]._columnRange.first
                + interaction._rowList[ row_idx ] - startRow;

      TRACE_ASSERT( systemRow >= 0, "Something is very wrong" );

      workspace[ work_idx ] = X[ systemRow ];

      work_idx++;
    }
  }

  // Multiply and subtract from X for this node's column range
  if ( work_idx > 0 ) {
    MATRIX::gemv( offDiagonalData(), workspace, baseData,
                  // Num. rows to consider == work_idx
                  work_idx, nCols, true, /* Transpose */
                  -1.0, /* alpha */
                  1.0 /* beta */ );
  }

  // Perform the final triangular solve
#if 0
  diagonalSolve( baseData, nRHS, NULL /* No extra data */,
                 true /* transpose */ );
#endif
  diagonalSolve( baseData, NULL /* No extra data */, true /* transpose */ );
}

//////////////////////////////////////////////////////////////////////
// Returns the index of the first extended node interaction
// in this node's off-diagonal
//////////////////////////////////////////////////////////////////////
int Supernode::firstExtendedInteraction() const
{
  int                            start_interaction = 0;

  while ( start_interaction < _offDiagonal.size()
       && _offDiagonal[ start_interaction ]._type != EXTENDED_NODE )
  {
    start_interaction++;
  }

  return start_interaction;
}

//////////////////////////////////////////////////////////////////////
// Puts the entries from this supernode lying in the given
// sub matrix range in to S
//
// Note: no bounds checking done here
//////////////////////////////////////////////////////////////////////
void Supernode::copySubMatrix( const IndexRange &rowRange,
                               const IndexRange &columnRange,
                               const Real *extraData,
                               const vector<Supernode> &nodes,
                               SPARSE_MATRIX &S )
{
  int                            nCols = numColumns();
  int                            nRows;
  int                            full_row_idx, full_col_idx;
  const Real                    *baseData;

  TRACE_ASSERT( NULL, "This is broken, and needs to be fixed to accomodate"
                      " diagonal sparsification" );

  if ( _columnRange.first > columnRange.second
    || _columnRange.second < columnRange.first )
  {
    return;
  }

  // Copy the diagonal block if this supernode's column range
  // intersects with the desired row range
  if ( _columnRange.first <= rowRange.second
    && _columnRange.second >= rowRange.first )
  {
    cout << "Copying diagonal data" << endl;

    baseData = ( _type == STANDARD_NODE ) ? _data
                                          : extraData + _extendedDataOffset;

    for ( int row_idx = 0; row_idx < nCols; row_idx++ )
    {
      full_row_idx = _columnRange.first + row_idx;

      if ( full_row_idx < rowRange.first
        || full_row_idx > rowRange.second )
      {
        continue;
      }

      for ( int col_idx = 0; col_idx < nCols; col_idx++ )
      {
        full_col_idx = _columnRange.first + col_idx;

        if ( full_col_idx < columnRange.first
          || full_col_idx > columnRange.second )
        {
          continue;
        }

        // FIXME: diagonal sparsification
        
        S( full_row_idx - rowRange.first, full_col_idx - columnRange.first )
          = MATRIX::access( baseData, nCols, nCols, row_idx, col_idx );
      }
    }
  }

  // Copy off-diagonal stuff
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];
    const Supernode             &ancestor = nodes[ interaction._nodeID ];

    if ( ancestor._columnRange.first > rowRange.second
      || ancestor._columnRange.second < rowRange.first )
    {
      continue;
    }

    baseData = ( interaction._type == STANDARD_NODE ) ?
                    _data + nCols * interaction._dataRowList[ 0 ]
                  : extraData + interaction._extendedDataOffset;

    nRows = ( interaction._type == STANDARD_NODE ) ?
                    interaction._rowList.size()
                  : interaction._numExtendedRows;

    for ( int row_idx = 0; row_idx < nRows; row_idx++ )
    {
      full_row_idx = ( interaction._type == STANDARD_NODE ) ?
                ancestor._columnRange.first + interaction._rowList[ row_idx ]
              : ancestor._columnRange.first + row_idx;

      if ( full_row_idx < rowRange.first
        || full_row_idx > rowRange.second )
      {
        continue;
      }

      for ( int col_idx = 0; col_idx < nCols; col_idx++ )
      {
        full_col_idx = _columnRange.first + col_idx;

        if ( full_col_idx < columnRange.first
          || full_col_idx > columnRange.second )
        {
          continue;
        }

        S( full_row_idx - rowRange.first, full_col_idx - columnRange.first )
          = MATRIX::access( baseData, nRows, nCols, row_idx, col_idx );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Counts non-zero columns in each interaction
//////////////////////////////////////////////////////////////////////
void Supernode::countNonZeroInteractionColumns(
                                  const vector<Supernode> &nodes,
                                  const SPARSE_MATRIX::SparseColumnMatrix &A )
{
  int                        rowStart;
  int                        rowEnd;
  int                        row_idx;

  // Nothing to do for extended nodes
  if ( _type == EXTENDED_NODE )
  {
    return;
  }

  // Count the non-zero columns in A for each interaction
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    SupernodeInteraction    &interaction = _offDiagonal[ interaction_idx ];

    const Supernode         &ancestor = nodes[ interaction._nodeID ];

    interaction._nonZeroColumns = 0;

    if ( interaction._type == EXTENDED_NODE )
    {
      continue;
    }

    rowStart = ancestor._columnRange.first;
    rowEnd = ancestor._columnRange.second;

    // Check each column from the original matrix for non-zero entries
    // in the row range of this interaction
    for ( int col_idx = _columnRange.first; col_idx <= _columnRange.second;
          col_idx++ )
    {
      for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
            row_ptr++ )
      {
        row_idx = A._i[ row_ptr ];

        if ( row_idx >= rowStart )
        {
          if ( row_idx <= rowEnd )
          {
            // Found a non-zero entry in this column of the interaction
            interaction._nonZeroColumns++;
          }

          // We can stop as soon as we have reached the start
          // row for this interaction, since rows are stored in
          // sorted order
          break;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Counts non-zero columns in each low rank diagonal block
//////////////////////////////////////////////////////////////////////
void Supernode::countNonZeroDiagonalBlockColumns(
                              const SPARSE_MATRIX::SparseColumnMatrix &A )
{
  int                        baseColumn = _columnRange.first;
  IndexRange                 columnRange;
  int                        row_idx;

  _diagonalLowRankBlockNonZeroColumns.clear();

  for ( int block_idx = 0; block_idx < _diagonalLowRankBlocks.size();
        block_idx++ )
  {
    const DenseBlock        &block = _diagonalLowRankBlocks[ block_idx ];
    int                      nonZeroColumns = 0;

    columnRange.first = baseColumn + block._columnRange.first;
    columnRange.second = baseColumn + block._columnRange.second;

    for ( int col_idx = columnRange.first; col_idx <= columnRange.second;
          col_idx++ )
    {
      for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
            row_ptr++ )
      {
        row_idx = A._i[ row_ptr ];

        if ( in_range( block._rowRange, row_idx - baseColumn ) )
        {
          nonZeroColumns++;
          break;
        }
      }
    }

    _diagonalLowRankBlockNonZeroColumns.push_back( nonZeroColumns );
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void Supernode::writeDiagonal( const char *prefix ) const
{
  const Real                *baseData = _data;
  int                        nCols;
  char                       buf[ 1024 ];

  if ( !lowRankDiagonal() )
  {
    nCols = numColumns();

    MATRIX                   tmp( nCols, nCols, baseData );

    sprintf( buf, "super_numeric/node_%d_%sdiagonal.matrix", _nodeID, prefix );

    tmp.write( buf );
  }
  else
  {
    for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ )
    {
      const DenseBlock      &block = _diagonalBlocks[ block_idx ];

      nCols = block.numColumns();

      MATRIX                tmp( nCols, nCols, baseData );

      sprintf( buf, "super_numeric/node_%d_%sdiagonal_block_%d.matrix",
               _nodeID, prefix, block_idx );

      tmp.write( buf );

      baseData += nCols * nCols;
    }
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void Supernode::writeOffDiagonal( const char *prefix ) const
{
  int                        nRows;
  int                        nCols = numColumns();

  const Real                *baseData;
  char                       buf[ 1024 ];

  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    if ( !interaction._active || interaction._type == EXTENDED_NODE )
    {
      continue;
    }

    nRows = interaction._rowList.size();

    baseData = interactionMatrix( interaction_idx );

    sprintf( buf, "%s_node_%d_interaction_%d.matrix",
             prefix, _nodeID, interaction_idx );

    MATRIX tmp( nRows, nCols, baseData );

    tmp.write( buf );
  }
}

//////////////////////////////////////////////////////////////////////
// Adds to a dense sub-matrix containing only the nodes specified
// in the given list.  This is here for debugging purposes
//////////////////////////////////////////////////////////////////////
void Supernode::extractSubMatrix( const IntArray &nodeList,
                                  const vector<Supernode> &nodes,
                                  Real *subMatrix, int nRows ) const
{
  int                        startRow = 0;
  int                        block_node_idx = 0;
  int                        nCols = numColumns();
  int                        subMatrixLDA = nRows;

  Real                      *baseData;
  const Real                *baseOffDiagonalData;

  TRACE_ASSERT( _type == STANDARD_NODE && _diagonalLowRankBlocks.size() == 0 );

  // Find this node's position in the sub-list of nodes
  for ( ; block_node_idx < nodeList.size(); block_node_idx++ )
  {
    if ( nodeList[ block_node_idx ] == _nodeID )
    {
      // We've found this node's position in the list
      break;
    }

    // Navigate to the correct starting column for this node
    startRow += nodes[ nodeList[ block_node_idx ] ].numColumns();
  }

  TRACE_ASSERT( block_node_idx < nodeList.size(), "Node index not found" );
  TRACE_ASSERT( startRow >= 0 );

  // Copy this node's diagonal block
  baseData = subMatrix + startRow * subMatrixLDA + startRow;

  for ( int row_idx = 0; row_idx < nCols; row_idx++ )
  {
    MATRIX::copyRow( _data, baseData, row_idx, row_idx, nCols, subMatrixLDA,
                      nCols /* Number of columns to actually copy */ );
  }

  // Move to the start of the off-diagonal part
  baseData += nCols * subMatrixLDA;

  // Go to the next node in our list
  block_node_idx++;

  // Copy relevant interactions
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // Skip until we're at the next node index in our sub-list
    while ( block_node_idx < nodeList.size()
         && nodeList[ block_node_idx ] < interaction._nodeID )
    {
      // Skip to the next matrix block
      baseData
        += nodes.at( nodeList[ block_node_idx ] ).numColumns() * subMatrixLDA;

      block_node_idx++;
    }

    if ( block_node_idx >= nodeList.size() )
    {
      // We've encountered all of the relevant interactions
      break;
    }

    // Check to see if this interaction is in the node list we're
    // concerned with
    if ( interaction._nodeID != nodeList[ block_node_idx ] )
    {
      TRACE_ASSERT( nodeList[ block_node_idx ] > interaction._nodeID );
      continue;
    }

    TRACE_ASSERT( interaction._type == STANDARD_NODE && interaction._active  );

    baseOffDiagonalData = interactionMatrix( interaction_idx );
    
    for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ )
    {
      MATRIX::copyRow( baseOffDiagonalData, baseData, row_idx,
                       interaction._rowList[ row_idx ],
                       nCols, subMatrixLDA, nCols );
    }
#if 0
    baseData += nodes[ interaction._nodeID ].numColumns() * subMatrixLDA;
#endif
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Adds to a sparse sub-matrix containing only the nodes
// specified in the given list.  This is here for debugging purposes
//////////////////////////////////////////////////////////////////////
void Supernode::extractSubMatrix( const IntArray &nodeList,
                                  const vector<Supernode> &nodes,
                                  SPARSE_MATRIX &subMatrix, int nRows ) const
{
  int                        startRow = 0;
  int                        block_node_idx = 0;
  int                        nCols = numColumns();
  int                        subMatrixLDA = nRows;

  Real                      *baseData;
  const Real                *baseOffDiagonalData;

  TRACE_ASSERT( _type == STANDARD_NODE && _diagonalLowRankBlocks.size() == 0 );

  // Find this node's position in the sub-list of nodes
  for ( ; block_node_idx < nodeList.size(); block_node_idx++ )
  {
    if ( nodeList[ block_node_idx ] == _nodeID )
    {
      // We've found this node's position in the list
      break;
    }

    // Navigate to the correct starting column for this node
    startRow += nodes[ nodeList[ block_node_idx ] ].numColumns();
  }

  TRACE_ASSERT( block_node_idx < nodeList.size(), "Node index not found" );
  TRACE_ASSERT( startRow >= 0 );

  for ( int row_idx = 0; row_idx < nCols; row_idx++ )
  for ( int col_idx = row_idx; col_idx < nCols; col_idx++ )
  {
    subMatrix( row_idx + startRow, col_idx + startRow )
      = MATRIX::access( _data, nCols, nCols, row_idx, col_idx );
  }

  startRow += nCols;

  // Move to the start of the off-diagonal part
  baseData += nCols * subMatrixLDA;

  // Copy relevant interactions
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

    // Skip until we're at the next node index in our sub-list
    if ( interaction._nodeID != block_node_idx )
    {
      continue;
    }

    TRACE_ASSERT( interaction._type == STANDARD_NODE && interaction._active  );

    baseOffDiagonalData = interactionMatrix( interaction_idx );
    
    for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ )
    {
      MATRIX::copyRow( baseOffDiagonalData, baseData, row_idx,
                       interaction._rowList[ row_idx ],
                       nCols, subMatrixLDA, nCols );
    }

    baseData += nodes[ interaction._nodeID ].numColumns() * subMatrixLDA;
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Writes an extended supernode
//////////////////////////////////////////////////////////////////////
void Supernode::writeSupernode( const char *prefix, const Real *extraData ) const
{
  const Real                *baseData = extraData + _extendedDataOffset;
  char                       buf[ 1024 ];

  int                        nRows = numColumns();

  TRACE_ASSERT( _type == EXTENDED_NODE );

  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    nRows += _offDiagonal[ interaction_idx ]._numExtendedRows;
  }

  sprintf( buf, "%s_extended_node_%06d_data.matrix", prefix, _nodeID );

  MATRIX tmp( nRows, numColumns(), baseData );

  tmp.write( buf );
}

//////////////////////////////////////////////////////////////////////
// Another function necessary for LDL factorization.  Off-diagonal
// interactions initially store permuted versions of what they should.
// This function "de-permutes" the content's of descedent's interaction
// with this node.
//////////////////////////////////////////////////////////////////////
void Supernode::permuteDescendentInteraction( const Supernode &descendent,
                                              int start_idx,
                                              Real *extraData )
{
  int                        nCols = numColumns();
  int                        descendentCols = descendent.numColumns();

  Real                      *baseData;

  const SupernodeInteraction &interaction
                                  = descendent.offDiagonal()[ start_idx ];

  TRACE_ASSERT( interaction._nodeID == _nodeID );
  TRACE_ASSERT( _LDL_data._pivotData.size() == nCols );

  baseData = extraData + interaction._extendedDataOffset;

  // Apply our permutation
  MATRIX::applyLDLPermutation( _LDL_data._pivotData, baseData,
                               nCols, descendentCols,
                               /* Left side, transpose permutation */
                               true, true );
}

//////////////////////////////////////////////////////////////////////
// Given a list of supernodes indexed by node IDs, build
// off diagonal interactions based on non-zero indices
// in the given list of super nodes.
//////////////////////////////////////////////////////////////////////
void Supernode::BuildOffDiagonalInteractions(
                                vector<Supernode> &nodes,
                                const SPARSE_MATRIX::SparseColumnMatrix &A,
                                bool compressInPlace )
{
  // Keeps track of the set of rows stored in each supernode's
  // off-diagonal.
  vector<set<int> >              rowSets( nodes.size() );
  IntArray                       nodeMap;

  buildSupernodeMap( nodes, nodeMap );

  // Start by initializing each node according to what rows
  // have entries in the original matrix.
  findSystemRows( nodes, rowSets, A );

  // Figure out all additional off-diagonal rows resulting
  // from fill-in
  propagateFillIn( nodes, rowSets, nodeMap, compressInPlace );

  // Fill in the actual interaction data
  buildSupernodeInteractions( nodes, rowSets, nodeMap );
}

//////////////////////////////////////////////////////////////////////
// Takes an existing list of supernodes, and introduces slack
// variable nodes in order to eliminate any off-diagonal blocks
// in the factorization beyond a given size.
//////////////////////////////////////////////////////////////////////
void Supernode::ExtendSystem( std::vector<Supernode> &nodes,
                              vector<Supernode> &newSystem,
                              int maxBlockSize )
{
  TRACE_ASSERT( NULL, "Needs to be updated" );

  // Map from old node IDs to new IDs - needed because we are
  // introducing new supernodes inline
  IntArray               nodeIDMap( nodes.size() );

  newSystem.clear();

  // Check each node to see if we need to introduce slack
  // variables
  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  {
    nodeIDMap.push_back( newSystem.size() );

    addNode( nodes[ node_idx ], nodes, newSystem, maxBlockSize );
  }

  // Use the node ID map to correct node-node interactions
  for ( int node_idx = 0; node_idx < newSystem.size(); node_idx++ )
  {
    Supernode           &node = newSystem[ node_idx ];

    for ( int block_idx = 0; block_idx < node._offDiagonal.size(); block_idx++ )
    {
      SupernodeInteraction &interaction = node._offDiagonal[ block_idx ];

      // For all interactions referring to non-extended nodes,
      // we need to correct the index that they refer to.
      //
      // Interactions with extended nodes should already refer
      // to the correct index.
      if ( interaction._type == STANDARD_NODE )
      {
        interaction._nodeID = nodeIDMap[ interaction._nodeID ];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// This function also extends a system, but introduces all slack
// variables at the end of the system.
//////////////////////////////////////////////////////////////////////
void Supernode::ExtendSystemAppend( vector<Supernode> &nodes,
                                    vector<Supernode> &newSystem,
                                    int maxBlockSize,
                                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                                    Real compressionRatio,
                                    RankEstimator *estimator,
                                    bool compressInPlace )
{
  vector<set<SupernodeInteraction> > interactionSets( nodes.size() );
  int                                oldSize = nodes.size();

  printf( "Extending system with %d nodes.\n", (int)nodes.size() );

  newSystem = nodes;

#if 0
  if ( maxBlockSize <= 0 )
    return;
#endif

#if 0
  // At this point, we need to figure out which rows we
  // can remove from interactions, based on the removal
  // of low rank interactions
  ComputeSparsifiedFillIn( newSystem, maxBlockSize, A );
#endif

  // Check each node to see if we need to introduce slack
  // variables
  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  {
    addNodeAppend( node_idx, newSystem, maxBlockSize,
                   interactionSets, compressionRatio, estimator,
                   compressInPlace );
  }

  printf( "Added %d nodes.\n", (int)newSystem.size() - (int)oldSize );
  printf( "New system has %d nodes.\n", (int)newSystem.size() );

  // Do variable reordering, provided that we have estimates for
  // extended variable sizes
  if ( estimator )
  {
    ReorderExtendedVariables( newSystem, interactionSets, oldSize );
  }

  // Add interaction sets for each extended node, and propagate fill in
  interactionSets.resize( newSystem.size() );

  TRACE_ASSERT( newSystem.size() >= oldSize );
  if ( newSystem.size() == oldSize )
  {
    return;
  }

  TRACE_ASSERT( interactionSets[ oldSize ].size() == 0,
                "Invalid interaction set" );

  for ( int node_idx = 0; node_idx < newSystem.size(); node_idx++ )
  {
    propagateExtendedFillIn( newSystem[ node_idx ], newSystem,
                             interactionSets );
  }

  // Now, add all new interactions to the system.  This will
  // also set up some local variables in the interactions themselves.
  addExtendedInteractions( newSystem, interactionSets, oldSize );
}

//////////////////////////////////////////////////////////////////////
// Performs symbolic reordering of the extended variable space to
// reduce fill-in.  This modifies the given interaction sets
//////////////////////////////////////////////////////////////////////
void Supernode::ReorderExtendedVariables(
                      vector<Supernode> &nodes,
                      vector<set<SupernodeInteraction> > &interactionSets,
                      int oldSystemSize )
{
  vector<set<int> >          schurInteractions;
#if 0
  vector<set<int> >          offDiagonalInteractions( oldSystemSize );
#endif
  vector<IntArray>           offDiagonalInteractions( oldSystemSize );
  vector<vector<IndexRange> > offDiagonalColumnRanges( oldSystemSize );
  IntArray                   blockSizes( nodes.size() - oldSystemSize );
  IntArray                   permutation;
  IntArray                   inversePermutation;

  TRACE_ASSERT( interactionSets.size() == oldSystemSize,
                "System must be reordered before extended interaction sets"
                " are built" );

  // Set up block sizes using the estimates for each extended node
  for ( int node_idx = oldSystemSize; node_idx < nodes.size(); node_idx++ )
  {
    int                      sizeEstimate = nodes[ node_idx ]._sizeEstimate;

    TRACE_ASSERT( sizeEstimate > 0 );

    blockSizes.at( node_idx - oldSystemSize ) = sizeEstimate;

#if 0
    printf( "Block size for extended node %d = %d\n",
            node_idx - oldSystemSize, sizeEstimate );
#endif

    // Reset the size estimate to zero - we don't need it anymore
    nodes[ node_idx ]._sizeEstimate = 0;
  }

  // Build index lists for the off-diagonal interactions already
  // present in the standard node set
  for ( int node_idx = 0; node_idx < oldSystemSize; node_idx++ )
  {
    set<SupernodeInteraction> &interactions = interactionSets[ node_idx ];
#if 0
    set<int> &interactionIndices = offDiagonalInteractions[ node_idx ];

    for ( set<SupernodeInteraction>::iterator iter = interactions.begin();
          iter != interactions.end(); iter++ )
    {
      interactionIndices.insert( iter->_nodeID );
    }
#endif
    IntArray &interactionIndices = offDiagonalInteractions[ node_idx ];
    vector<IndexRange> &interactionColumnRanges
                                    = offDiagonalColumnRanges[ node_idx ];

    for ( set<SupernodeInteraction>::iterator iter = interactions.begin();
          iter != interactions.end(); iter++ )
    {
      interactionIndices.push_back( iter->_nodeID );

      if ( iter->_lowRankDiagonalInteraction )
      {
        interactionColumnRanges.push_back( iter->_compressedColumnRange );
      }
      else
      {
        interactionColumnRanges.push_back( IndexRange( 0, 0 ) );
      }
    }
  }

  // Build the schur complement interaction set
  Ordering::BuildBlockSchurComplement( offDiagonalInteractions,
                                       schurInteractions,
                                       oldSystemSize, /* start index */
                                       nodes.size(), /* end index */
                                       &offDiagonalColumnRanges );

  MinimumDegree::BlockMinimumDegree( schurInteractions, blockSizes,
                                     permutation, inversePermutation );

  // Diagnostic: print out usage
  MinimumDegree::PrintBlockOrderingStats( schurInteractions, blockSizes,
                                          permutation, inversePermutation );

  // Use the permutation to renumber all interactions in the extended set
  for ( int node_idx = 0; node_idx < oldSystemSize; node_idx++ )
  {
    set<SupernodeInteraction> &oldInteractions = interactionSets[ node_idx ];
    set<SupernodeInteraction>  newInteractions;

    Supernode                 &node = nodes[ node_idx ];

    for ( set<SupernodeInteraction>::iterator iter = oldInteractions.begin();
          iter != oldInteractions.end(); iter++ )
    {
      int                      oldNodeID = iter->_nodeID;
      int                      newNodeID;

      // Express this relative to the start of the extended variable space
      oldNodeID -= oldSystemSize;

      // Look up the new relative index in the permutation
      newNodeID = inversePermutation.at( oldNodeID );

      // Put this back in to full variable space
      newNodeID += oldSystemSize;

      SupernodeInteraction   interactionCopy = *iter;

      interactionCopy._nodeID = newNodeID;

      newInteractions.insert( interactionCopy );
    }

    // Now that we've built a new interaction set, we can just replace
    // the old one.
    oldInteractions = newInteractions;

    // Apply permutation to any extended offsets stored in the
    // standard part of the factorization
    for ( int interaction_idx = 0; interaction_idx < node._offDiagonal.size();
          interaction_idx++ )
    {
      SupernodeInteraction  &interaction = node._offDiagonal[ interaction_idx ];

      if ( interaction._active )
      {
        // Nothing to be done here
        continue;
      }

      TRACE_ASSERT( interaction._extendedOffset > 0,
                    "No extended offset found in inactive interaction" );

      int                    oldNodeID;
      int                    newNodeID;
      
      oldNodeID = node_idx + interaction._extendedOffset;

      // Express this relative to the start of the extended variable space
      // and permute
      oldNodeID -= oldSystemSize;
      newNodeID = inversePermutation.at( oldNodeID );
      newNodeID += oldSystemSize;

      interaction._extendedOffset = newNodeID - node_idx;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// For each node, examine the row set of each of its interactions
// in light of the fact that some of its interactions will be
// removed and replaced with low rank decomposition.  Adjust the
// symbolic factor to make use of these new, smaller row sets.
//////////////////////////////////////////////////////////////////////
void Supernode::ComputeSparsifiedFillIn(
                              vector<Supernode> &nodes,
                              int maxBlockSize,
                              const SPARSE_MATRIX::SparseColumnMatrix &A )
{
  // Linked lists stored during numerical factorization
  IntArray                   nextInteraction( nodes.size(), EMPTY );
  IntArray                   head( nodes.size(), EMPTY );
  IntArray                   next( nodes.size(), EMPTY );

  int                        descendent_idx;
  int                        descendent_idx_next;
  int                        start_idx;
  int                        descendentAncestor;
  int                        parentNode;

  int                        curr_interaction_idx;
  int                        desc_interaction_idx;

  int                        nRows;

  PairArray                  interactionIndices;
  PairArray                  scratch1, scratch2;

  vector<set<int> >          rowSets;

  // This mirrors the numerical factorization with nodes removed
  // (since they will be treated as low rank blocks)
  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  {
    if ( node_idx % 1000 == 0 )
    {
      printf( "computeSparsifiedFillIn: node %d of %d\r",
              node_idx + 1, (int)nodes.size() );
    }

    Supernode               &node = nodes[ node_idx ];

    // Figure out which rows belong in this node's interactions
    // based just on the original matrix
    node.buildBaseSystemRowSets( nodes, rowSets, A );

    TRACE_ASSERT( rowSets.size() == node._offDiagonal.size(),
                  "Wrong row set size" );

    // Look through this node's descendents to see which other rows we need
    // to include.  Note that the list of descendents is built in a way
    // which excludes low rank contributions
    for ( descendent_idx = head[ node_idx ]; descendent_idx != EMPTY;
          descendent_idx = descendent_idx_next )
    {
      Supernode             &descendent = nodes[ descendent_idx ];

      descendent_idx_next = next[ descendent_idx ];

      // Figure out symbolically what the Schur complement
      // should look like based on the descendent list
      start_idx = nextInteraction[ descendent_idx ];

      node.findInteractionIntersection( descendent, start_idx,
                                        interactionIndices,
                                        scratch1, scratch2 );

      // For each interaction affected by this descendent, include any
      // rows that the descendent fills in to this node.
      for ( int interaction_idx = 0; interaction_idx < interactionIndices.size();
            interaction_idx++ )
      {
        desc_interaction_idx = interactionIndices[ interaction_idx ].first;
        curr_interaction_idx = interactionIndices[ interaction_idx ].second;

        set<int>             &rowSet = rowSets[ curr_interaction_idx ];

        SupernodeInteraction &descInteraction
                          = descendent._offDiagonal[ desc_interaction_idx ];

        SupernodeInteraction &currInteraction
                          = node._offDiagonal[ curr_interaction_idx ];

        TRACE_ASSERT( currInteraction._nodeID == descInteraction._nodeID,
                      "Interaction node ID mismatch" );

        for ( int row_idx = 0; row_idx < descInteraction._rowList.size();
              row_idx++ )
        {
          rowSet.insert( descInteraction._rowList[ row_idx ] );
        }
      }

      start_idx += 1;
      nextInteraction[ descendent_idx ] = start_idx;

      // Find the next node this descendent is going to interact with
      while ( start_idx < descendent._offDiagonal.size()
           && !descendent._offDiagonal[ start_idx ]._active )
      {
        // Skip inactive nodes, since they will not make any further
        // contributions to the factorization.
        start_idx += 1;
        nextInteraction[ descendent_idx ] = start_idx;
      }

      // Adjust linked list entries
      if ( start_idx < descendent._offDiagonal.size() )
      {
        descendentAncestor = descendent._offDiagonal[ start_idx ]._nodeID;

        TRACE_ASSERT( descendentAncestor > node_idx
                   && descendentAncestor < nodes.size(),
                   "Invalid descendent ancestor" );

        next[ descendent_idx ] = head[ descendentAncestor ];
        head[ descendentAncestor ] = descendent_idx;
      }
    }

    // Build a new off diagonal for this node
    vector<SupernodeInteraction>   oldInteractions = node._offDiagonal;

    node._offDiagonal.clear();

    nRows = 0;

    // Adjust the row set for each interaction in this node
    for ( int interaction_idx = 0; interaction_idx < oldInteractions.size();
          interaction_idx++ )
    {
      set<int>              &rowSet = rowSets[ interaction_idx ];

      // Only include this interaction if it has a non-empty row set
      if ( rowSet.size() == 0 )
      {
        continue;
      }

      nRows += rowSet.size();

      node._offDiagonal.push_back(
          SupernodeInteraction( oldInteractions[ interaction_idx ]._nodeID,
                                STANDARD_NODE ) );

      // Put the adjusted row list in to the new interaction
      SupernodeInteraction  &newInteraction = node._offDiagonal.back();

      newInteraction._rowList.reserve( rowSet.size() );
      newInteraction._rowList.clear();

      for ( set<int>::iterator i = rowSet.begin(); i != rowSet.end(); i++ )
      {
        newInteraction._rowList.push_back( *i );
      }
    }

    node._numRows = nRows;

    if ( maxBlockSize >= 0 && node.numColumns() > maxBlockSize )
    {
      // Temporarily flag all interactions exceeding maxBlockSize
      // as inactive
      for ( int interaction_idx = 0; interaction_idx < node._offDiagonal.size();
            interaction_idx++ )
      {
        node._offDiagonal[ interaction_idx ]._active = false;
      }
    }

    // Set up the linked list for this node
    start_idx = 0;

    while ( start_idx < node._offDiagonal.size()
         && !node._offDiagonal[ start_idx ]._active )
    {
      // Skip inactive nodes, since they will not make any further
      // contributions to the factorization.
      start_idx++;
    }

    // Update descendent relationships
    if ( start_idx < node._offDiagonal.size() )
    {
      nextInteraction[ node_idx ] = start_idx;

      parentNode = node._offDiagonal[ start_idx ]._nodeID;

      next[ node_idx ] = head[ parentNode ];
      head[ parentNode ] = node_idx;
    }

    head[ node_idx ] = EMPTY;
  }

  // Reset all interactions to be active
  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  for ( int interaction_idx = 0;
        interaction_idx < nodes[ node_idx ]._offDiagonal.size();
        interaction_idx++ )
  {
    nodes[ node_idx ]._offDiagonal[ interaction_idx ]._active = true;
  }

  printf( "\n" );
}

//////////////////////////////////////////////////////////////////////
// Given a post-ordering of a nested dissection tree for a
// matrix, this function builds a symbolic factor for a system
// extended from A to introduce slack variables.
//
// A is assumed to have been reordered according to the
// permutation given by the dissection node list.
//////////////////////////////////////////////////////////////////////
void Supernode::FactorSymbolic(
      const vector<const NestedDissection::DissectionNode *> &nodes,
      const SPARSE_MATRIX::SparseColumnMatrix &A,
      int maxBlockSize,
      vector<Supernode> &factor )
{
  TRACE_ASSERT( NULL, "Needs to be updated" );

  vector<Supernode>            initialNodes;

  factor.clear();

  BuildBasicFactor( nodes, A, initialNodes );

  BuildOffDiagonalInteractions( initialNodes, A );

  ExtendSystem( initialNodes, factor, maxBlockSize );
}

//////////////////////////////////////////////////////////////////////
// Given a post-ordering of a nested dissection tree for a
// matrix, this function builds the basic symbolic factor,
// in which for now nodes just have the correct column ranges.
//////////////////////////////////////////////////////////////////////
void Supernode::BuildBasicFactor(
      const vector<const NestedDissection::DissectionNode *> &nodes,
      const SPARSE_MATRIX::SparseColumnMatrix &A,
      vector<Supernode> &initialNodes )
{
  int                          lastIndex = 0;

  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  {
    initialNodes.push_back( Supernode() );

    Supernode                 &node = initialNodes.back();

    node._nodeID = node_idx;

    node._columnRange.first = lastIndex;
    node._columnRange.second = lastIndex + nodes[ node_idx ]->size() - 1;

    IndexRange                 columnRange = nodes[ node_idx ]->reorderedRange();

    node._firstOffDiagonalRow = range_size( node._columnRange );

    lastIndex = node._columnRange.second + 1;
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but just takes a list of supernode specifications
//////////////////////////////////////////////////////////////////////
void Supernode::BuildBasicFactor(
                    const vector<Ordering::SupernodeSpecification> &nodes,
                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                    vector<Supernode> &initialNodes )
{
  int                          lastIndex = 0;

  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  {
    initialNodes.push_back( Supernode() );

    Supernode                 &node = initialNodes.back();

    node._nodeID = node_idx;

    node._columnRange = nodes[ node_idx ]._columnRange;
    node._compressOffDiagonal = nodes[ node_idx ]._compressOffDiagonal;
    node._firstOffDiagonalRow = range_size( node._columnRange );

    lastIndex = node._columnRange.second + 1;
  }
}

//////////////////////////////////////////////////////////////////////
// Counts the number of non-zero entries stored by a list
// of supernodes
//////////////////////////////////////////////////////////////////////
long int Supernode::CountNonZeros( const std::vector<Supernode> &nodes )
{
  TRACE_ASSERT( NULL, "Needs to be updated" );

  long int                     nnz = 0;
  int                          nrows, ncols;

  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  {
    ncols = nodes[ node_idx ].numColumns();
    nrows = 0;

    for ( int block_idx = 0; block_idx < nodes[ node_idx ]._offDiagonal.size();
          block_idx++ )
    {
      if ( nodes[ node_idx ]._offDiagonal[ block_idx ]._active )
      {
        nrows += nodes[ node_idx ]._offDiagonal[ block_idx ]._rowList.size();
      }
    }

    // FIXME: diagonal sparsification

    nnz += ncols * ncols + ncols * nrows;
  }

  return nnz;
}

//////////////////////////////////////////////////////////////////////
// Allocates all required fixed data for standard nodes, and their
// interactions with other standard nodes.
//////////////////////////////////////////////////////////////////////
Real *Supernode::AllocateFixedData( std::vector<Supernode> &nodes,
                                    long int &dataSize )
{
  Real                         dataSizeMB;
  Real                        *data;
  Real                        *dataCopy;

  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  {
    dataSize += nodes[ node_idx ].numFixedEntries();
  }

  dataSizeMB = (Real)( dataSize * 8 ) / 1024.0 / 1024.0;

  printf( "Allocating %ld nonzeros, %f MB of data\n", dataSize, dataSizeMB );

  data = new Real[ dataSize ];
  dataCopy = data;

  memset( (void *)data, 0, dataSize * sizeof( Real ) );

  // Assign data pointers to the node list
  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  {
    // This only assigns data pointers to standard nodes.
    // Extended nodes are handled during factorization.
    if ( nodes[ node_idx ]._type == STANDARD_NODE )
    {
      nodes[ node_idx ]._data = dataCopy;

      dataCopy += nodes[ node_idx ].numFixedEntries();

      // FIXME: diagonal sparsification

      nodes[ node_idx ].setDataOffsets();
      nodes[ node_idx ].initDiagonalBlock();
    }
  }

  return data;
}

#if 0
//////////////////////////////////////////////////////////////////////
// Subtracts the product formed above from another workspace,
// using relativeMap to map rows from one workspace to the other.
// This is used when building a matrix to be factored in to
// a low-rank decomposition.
//////////////////////////////////////////////////////////////////////
void Supernode::ScatterLowRankUpdate( const Real *multWorkspaceFinal,
                                      Real *decompWorkspace,
                                      const IntArray &relativeMap,
                                      int nCols )
{
  for ( int row_idx = 0; row_idx < relativeMap.size(); row_idx++ )
  {
    // A := A - B (one row at a time)
    MATRIX::axpy( decompWorkspace + relativeMap[ row_idx ] * nCols,
                  multWorkspaceFinal + row_idx * nCols,
                  1, nCols, -1.0 );
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Similar to the above, but only handles a single interaction.
// That is, we scatter from one row list to another.
//////////////////////////////////////////////////////////////////////
void Supernode::ScatterLowRankUpdate( const Real *multWorkspaceFinal,
                                      Real *decompWorkspace,
                                      const SupernodeInteraction &interaction,
                                      const IntArray &fullRows,
                                      int nCols )
{
  int                          full_row_idx = 0;

  // Have to copy one row at a time in this case
  for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ )
  {
    // Find the appropriate row in the full matrix
    while ( full_row_idx < fullRows.size()
         && fullRows[ full_row_idx ] != interaction._rowList[ row_idx ] )
    {
      full_row_idx++;
    }

    TRACE_ASSERT( full_row_idx < fullRows.size(),
                  "No corresponding row found in the original interaction" );

    // Subtract the contribution
    MATRIX::axpy( decompWorkspace + full_row_idx * nCols,
                  multWorkspaceFinal + row_idx * nCols,
                  1, nCols, -1.0 );
  }
}

//////////////////////////////////////////////////////////////////////
// Another version of scatterLowRankUpdate.  This scatters a
// subset of one row list in to a full row space.
//////////////////////////////////////////////////////////////////////
void Supernode::ScatterLowRankUpdate( const Real *multWorkspaceFinal,
                                      Real *decompWorkspace,
                                      const SupernodeInteraction &interaction,
                                      const IndexRange &rowRange,
                                      int nCols,
                                      // Offset to subtract from row indices
                                      int offset )
{
  int                        full_row_idx;

  // Have to copy one row at a time
  for ( int row_idx = rowRange.first; row_idx <= rowRange.second; row_idx++ )
  {
    // Find the full row index, using the given offset
    full_row_idx = interaction._rowList[ row_idx ] - offset;

    // Subtract the contribution
    MATRIX::axpy( decompWorkspace + full_row_idx * nCols,
                  multWorkspaceFinal + ( row_idx - rowRange.first ) * nCols,
                  1, nCols, -1.0 );
  }
}

//////////////////////////////////////////////////////////////////////
// Scatters data from multiple interactions in descendent to
// the interaction set of ancestor, where ancestor_idx
// specifies the interaction between descendent and ancestor
//////////////////////////////////////////////////////////////////////
void Supernode::ScatterMultiLowRankUpdate( const Real *inputWorkspace,
                                           Real *outputWorkspace,
                                           const Supernode &descendent,
                                           int ancestor_idx,
                                           const Supernode &ancestor,
                                           int nCols )
{
  int                        interaction_idx = 0;
  int                        nRows;

  for ( int desc_interaction_idx = ancestor_idx + 1;
        desc_interaction_idx < descendent.offDiagonal().size();
        desc_interaction_idx++ )
  {
    // Every interaction in descendent must have a corresponding
    // interaction in ancestor
    while ( interaction_idx < ancestor.offDiagonal().size()
         && descendent.offDiagonal()[ desc_interaction_idx ]._nodeID
          != ancestor.offDiagonal()[ interaction_idx ]._nodeID )
    {
      // Move to the next position in the output workspace
      nRows = ancestor.offDiagonal()[ interaction_idx ]._rowList.size();
      outputWorkspace += nRows * nCols;

      interaction_idx++;
    }

    TRACE_ASSERT( interaction_idx < ancestor.offDiagonal().size(),
                  "No corresponding interaction found in ancestor" );

    // Scatter data
    ScatterLowRankUpdate( inputWorkspace, outputWorkspace,
                          descendent.offDiagonal()[ desc_interaction_idx ],
                          ancestor.offDiagonal()[ interaction_idx ]._rowList,
                          nCols );

    // Move to the next position in the input workspace
    nRows = descendent.offDiagonal()[ desc_interaction_idx ]._rowList.size();
    inputWorkspace += nRows * nCols;
  }
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, once we have formed the right half of
// a low rank decomposition Q * R', we must scatter the columns
// of the formed matrix R' according to the contents of the
// interaction between this node and ancestor.
//////////////////////////////////////////////////////////////////////
void Supernode::ScatterLowRankTransUpdate( int ancestor_idx,
                                           const Supernode &descendent,
                                           const Supernode &ancestor,
                                           const Real *multWorkspaceFinal,
                                           Real *decompWorkspace,
                                           int nRows )
{
  const IntArray &rowList = descendent._offDiagonal[ ancestor_idx ]._rowList;

  int                          nColsInput = rowList.size();
  int                          nColsOutput = ancestor.numColumns();

  for ( int col_idx = 0; col_idx < rowList.size(); col_idx++ )
  {
    // A := A - B (one column at a time)
    MATRIX::axpy( // Offset to start of columns
                  decompWorkspace + rowList[ col_idx ],
                  multWorkspaceFinal + col_idx,
                  nRows, 1, -1.0,
                  // Leading dimension is number of columns, since
                  // we are copying columns instead of rows
                  nColsOutput, nColsInput );
  }
}

//////////////////////////////////////////////////////////////////////
// Another version of scatterLowRankTransUpdate.  This scatters a
// subset of one column list in to a full column space
//
// Used in decomposition of low rank blocks in a node's diagonal
//////////////////////////////////////////////////////////////////////
void Supernode::ScatterLowRankTransUpdate(
                                const Real *multWorkspaceFinal,
                                Real *decompWorkspace,
                                const SupernodeInteraction &interaction,
                                const IndexRange &columnRange,
                                int nRows,
                                int nColsInput, int nColsOutput,
                                // Offset to subtract from column indices
                                int offset )
{
  const IntArray            &rowList = interaction._rowList;

  for ( int col_idx = columnRange.first; col_idx <= columnRange.second;
        col_idx++ )
  {
    // A := A - B (one column at a time)
    MATRIX::axpy( // Offset to start of columns
                  decompWorkspace + ( rowList[ col_idx ] - offset ),
                  multWorkspaceFinal + ( col_idx - columnRange.first ),
                  nRows, 1, -1.0,
                  // Leading dimension is number of columns, since
                  // we are copying columns instead of rows
                  nColsOutput, nColsInput );
  }
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, after a left transpose interaction multiply,
// we must scatter the rows of multWorkspaceFinal in to the
// appropriate locations of decompWorkspace
//////////////////////////////////////////////////////////////////////
void Supernode::ScatterLowRankTransUpdate_row( int ancestor_idx,
                                               const Supernode &descendent,
                                               const Supernode &ancestor,
                                               const Real *multWorkspaceFinal,
                                               Real *decompWorkspace,
                                               int nCols )
{
  const IntArray &rowList = descendent._offDiagonal[ ancestor_idx ]._rowList;

  int                          nColsInput = rowList.size();
  int                          nColsOutput = ancestor.numColumns();

  // Scatter in to the full column space of ancestor
  for ( int row_idx = 0; row_idx < rowList.size(); row_idx++ )
  {
    // A := A - B (one column at a time)
    MATRIX::axpy( // Offset to start of rows
                  decompWorkspace + rowList[ row_idx ] * nCols,
                  multWorkspaceFinal + row_idx * nCols,
                  1, nCols, -1.0,
                  // Leading dimension is 1, since we are copying rows
                  1, 1 );
  }
}

//////////////////////////////////////////////////////////////////////
// Similar to the above, but now we assume that fullMatrix is itself
// only defined over some set of rows, of which interaction stores
// a subset
//////////////////////////////////////////////////////////////////////
void Supernode::BuildInteractionSubMatrix(
                             const SupernodeInteraction &interaction,
                             const IntArray &fullRows,
                             const Real *fullMatrix, Real *subMatrix,
                             int nCols )
{
  int                          full_row_idx = 0;

  // Have to copy one row at a time in this case
  for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ )
  {
    // Find the appropriate row in the full matrix
    while ( full_row_idx < fullRows.size()
         && fullRows[ full_row_idx ] != interaction._rowList[ row_idx ] )
    {
      full_row_idx++;
    }

    TRACE_ASSERT( full_row_idx < fullRows.size(),
                  "No corresponding row found in the original interaction" );

    // Copy the row
    MATRIX::copyRow( fullMatrix, subMatrix, full_row_idx, row_idx,
                     nCols, nCols );
  }
}

//////////////////////////////////////////////////////////////////////
// Builds an interaction submatrix in order to multiply by multiple
// interactions from descendent, assuming that the row space of
// the input matrix is determined by the row space of ancestor
//////////////////////////////////////////////////////////////////////
void Supernode::BuildMultiInteractionSubMatrix(
                            const Supernode &descendent,
                            int ancestor_idx,
                            const Supernode &ancestor,
                            const Real *fullMatrix, Real *subMatrix,
                            int nCols )
{
  int                        interaction_idx = 0;
  int                        nRows;

  int                        totalRows = 0;

  for ( int desc_interaction_idx = ancestor_idx + 1;
        desc_interaction_idx < descendent.offDiagonal().size();
        desc_interaction_idx++ )
  {
    // Every interaction in descendent must have a corresponding
    // interaction in ancestor
    while ( interaction_idx < ancestor.offDiagonal().size()
         && descendent.offDiagonal()[ desc_interaction_idx ]._nodeID
          != ancestor.offDiagonal()[ interaction_idx ]._nodeID )
    {
      // Move to the next position in the input workspace
      nRows = ancestor.offDiagonal()[ interaction_idx ]._rowList.size();
      fullMatrix += nRows * nCols;

      interaction_idx++;
    }

    TRACE_ASSERT( interaction_idx < ancestor.offDiagonal().size(),
                  "No corresponding interaction found in ancestor" );

#if 0
    // Scatter data
    ScatterLowRankUpdate( inputWorkspace, outputWorkspace,
                          descendent.offDiagonal()[ desc_interaction_idx ],
                          ancestor.offDiagonal()[ interaction_idx ]._rowList,
                          nCols );
#endif
    // Grab the necessary rows from the input data and put them in the
    // sub-matrix
    BuildInteractionSubMatrix(
                          descendent.offDiagonal()[ desc_interaction_idx ],
                          ancestor.offDiagonal()[ interaction_idx ]._rowList,
                          fullMatrix, subMatrix, nCols );

    // Move to the next position in the sub-matrix
    nRows = descendent.offDiagonal()[ desc_interaction_idx ]._rowList.size();
    subMatrix += nRows * nCols;

    totalRows += nRows;
  }

#if 0
  cout << "Found " << SDUMP( totalRows ) << endl;
#endif
}

//////////////////////////////////////////////////////////////////////
// Reverse of the above.  Scatters rows indexed by interaction._rowList
// to fullRows in fullMatrix.
//////////////////////////////////////////////////////////////////////
void Supernode::scatterInteractionMatrix(
                             const SupernodeInteraction &interaction,
                             const IntArray &fullRows,
                             const Real *subMatrix, Real *fullMatrix,
                             int nCols )
{
  int                          full_row_idx = 0;

  // Have to copy one row at a time in this case
  for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ )
  {
    // Find the appropriate row in the full matrix
    while ( full_row_idx < fullRows.size()
         && fullRows[ full_row_idx ] != interaction._rowList[ row_idx ] )
    {
      full_row_idx++;
    }

    TRACE_ASSERT( full_row_idx < fullRows.size(),
                  "No corresponding row found in the original interaction" );

    // Copy the row
    MATRIX::copyRow( subMatrix, fullMatrix, row_idx, full_row_idx,
                     nCols, nCols );
  }
}

//////////////////////////////////////////////////////////////////////
// Sets the _dataRowList array in each of this node's
// interactions
//////////////////////////////////////////////////////////////////////
void Supernode::setDataOffsets()
{
  // FIXME: diagonal sparsification

  int                          dataOffset = _firstOffDiagonalRow;

  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    SupernodeInteraction      &interaction = _offDiagonal[ interaction_idx ];

    if ( interaction._type == STANDARD_NODE && interaction._active )
    {
      interaction._dataRowList.resize( interaction._rowList.size() );

      for ( int row_idx = 0; row_idx < interaction._rowList.size(); row_idx++ )
      {
        interaction._dataRowList[ row_idx ] = dataOffset;

        dataOffset++;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void Supernode::factorCholesky( Real *extraData, int writeSize,
                                WorkspaceManager<Real> *workspaceManager )
{
  int                            nCols = numColumns();
  int                            nRows = _numRows + nCols;
  Real                          *diagonalData;

  diagonalData = ( _type == STANDARD_NODE ) ? _data
                                            : extraData + _extendedDataOffset;

  // Write all blocks larger than writeSize
  if ( writeSize > 0 )
  {
    for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
          interaction_idx++ )
    {
      const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

      if ( interaction._type == STANDARD_NODE && interaction._active
        //&& ( interaction._rowList.size() * nCols ) >= writeSize )
        && ( nCols ) >= writeSize )
      {
        MATRIX tmp( interaction._rowList.size(), nCols,
                    _data + nCols * interaction._dataRowList[ 0 ] );

        char buf[ 1024 ];

        sprintf( buf, "super_numeric/node_%09d_interaction_%02d.matrix",
                 _nodeID, interaction_idx );

        tmp.write( buf );
      }
    }
  }

  // Computes dense Cholesky factor(s) for the diagonal block(s)
  factorDiagonal( extraData, workspaceManager );
}

//////////////////////////////////////////////////////////////////////
// Factors this node's diagonal block using an LDLT factorization,
// storing the result in _LDL_data
//////////////////////////////////////////////////////////////////////
void Supernode::factorLDL( Real *extraData,
                           WorkspaceManager<Real> &workspaceManager )
{
  TRACE_ASSERT( _type == EXTENDED_NODE );

  int                        nCols = numColumns();
  Real                      *diagonalData = extraData + _extendedDataOffset;

  // Workspaces for extracting matrices from the Bunch-Kaufman factorization
  IntArray                   permutation;
  IntArray                   inversePermutation;

  // Initialize this nodes LDL data structure and compute the
  // Bunch-Kaufman factorization of it's diagonal block
  _LDL_data._pivotData.resize( nCols );

  MATRIX::LDL( diagonalData, nCols, _LDL_data._pivotData );

  // Pull out the 1x1 and 2x2 diagonal blocks
  ExtractLDLDiagonal( diagonalData, nCols, _LDL_data );

  // Extract a psychologically lower triangular matrix from the
  // LDL matrix, and overwrite this node's diagonal data with it
  {
    RealWorkspace            workspace( workspaceManager,
                                        IndexPair( nCols, nCols ) );

    IntArray                 permutation;
    IntArray                 inversePermutation;

    MATRIX::buildLDLTriangle( _LDL_data._pivotData, diagonalData,
                              permutation, inversePermutation,
                              workspace.workspaceData( 0 ) );

    MATRIX::copy( diagonalData, workspace.workspaceData( 0 ), nCols, nCols );
  }

  // Permute the psychologically lower triangular matrix using the
  // inverse permutation, which will put it in lower-triangular form
  MATRIX::applyLDLPermutation( _LDL_data._pivotData, diagonalData,
                               nCols, nCols,
                               /* Left side, transposed multiplication */
                               true, true );
}

//////////////////////////////////////////////////////////////////////
// Off diagonal solve for Cholesky-factored node
//////////////////////////////////////////////////////////////////////
void Supernode::choleskyOffDiagonalSolve(
                              bool solveCompressedBlocks,
                              Real *extraData, int writeSize,
                              WorkspaceManager<Real> *workspaceManager )
{
  Real                      *baseData;

  if ( _numRows > 0 )
  {
    TRACE_ASSERT( _type == STANDARD_NODE,
                  "Extended nodes cannot have standard off diagonal rows" );

    TRACE_ASSERT( _firstOffDiagonalRow > 0,
                  "First diagonal row not specified" );

#if 0
    if ( _nodeID == 9225 ) {
      printf( "Writing off diagonal data, firstRow = %d, nCols = %d\n",
              _firstOffDiagonalRow, numColumns() );
      MATRIX::write( offDiagonalData(), _numRows, numColumns(),
                     "super_numeric/realInput.matrix" );
    }
#endif

    diagonalSolve( offDiagonalData(), _numRows, extraData,
                   true /* transpose */, false /* right side solve */ );

#if 0
    // FIXME: debugging
    if ( _nodeID < 253 ) {
      for ( int i = 0; i < _offDiagonal.size(); i++ ) {
        if ( _offDiagonal[ i ]._nodeID == 253
          && _offDiagonal[ i ]._type == STANDARD_NODE )
        {
          Real *data = interactionMatrix( i );

          char buf[ 1024 ];
          sprintf( buf, "test/node_%d_interaction.matrix", _nodeID );

          MATRIX::write( data, _offDiagonal[ i ]._rowList.size(), numColumns(),
                         buf );
        }
      }
    } else {
      abort();
    }
#endif

#if 0
    // FIXME: debugging
    if ( _nodeID == 9225 ) {
      MATRIX::write( offDiagonalData(), _numRows, numColumns(),
                     "super_numeric/realData.matrix" );
      MATRIX::write( _data, numColumns(), numColumns(),
                     "super_numeric/realDiagonal.matrix" );
    }
#endif
  }

  if ( !extraData )
  {
    return;
  }

  // If we have extra data, we can also perform triangular solves
  // FIXME: Technically we could start at the first extended node,
  //        assuming they all come at the end.
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    // Skip standard nodes, and extended nodes stored in compressed
    // form
    if ( _offDiagonal[ interaction_idx ]._type == STANDARD_NODE
      || _offDiagonal[ interaction_idx ]._compressed )
    {
      continue;
    }

    // FIXME: treat all extended interactions as compressed
    if ( _type == STANDARD_NODE
      && _offDiagonal[ interaction_idx ]._type == EXTENDED_NODE )
    {
      continue;
    }

    if ( _offDiagonal[ interaction_idx ]._type == EXTENDED_NODE ) {
      baseData
        = extraData + _offDiagonal[ interaction_idx ]._extendedDataOffset;

      TRACE_ASSERT( _offDiagonal[ interaction_idx ]._extendedDataOffset >= 0,
                    "Invalid extended data offset" );

      diagonalSolve( baseData,
                     _offDiagonal[ interaction_idx ]._numExtendedRows,
                     extraData,
                     true /* transpose */, false /* right side */ );
    }
    else if ( solveCompressedBlocks
           && _offDiagonal[ interaction_idx ]._type == COMPRESSED_NODE ) {
      // Do a diagonal solve on the "U" part of the factorization.
      // Here, we do a left, non-transposed solve, since U is not stored
      // transposed.
      SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

      diagonalSolve( interaction._compressedU,
                     interaction._numExtendedRows,
                     extraData,
                     false /* don't transpose */, true /* left side */ );
    }
    // FIXME: debugging
    else if ( !solveCompressedBlocks
           && _offDiagonal[ interaction_idx ]._type == COMPRESSED_NODE ) {
#if 0
      // Just write stuff to disk
      SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

      char buf[ 1024 ];
      sprintf( buf, "super_numeric/node_%d_interaction_%d_U.matrix",
               _nodeID, interaction_idx );

      MATRIX::write( interaction._compressedU, numColumns(),
                     interaction._numExtendedRows, buf );
#endif
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Off-diagonal solve for LDL-factored node
//////////////////////////////////////////////////////////////////////
void Supernode::LDLOffDiagonalSolve(
                              Real *extraData,
                              WorkspaceManager<Real> &workspaceManager )
{
  int                        nCols = numColumns();

  const Real                *diagonalData = extraData + _extendedDataOffset;
  Real                      *baseData;

  // Form LDL off-diagonal segments
  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    // Skip standard nodes, and extended nodes stored in compressed
    // form
    if ( _offDiagonal[ interaction_idx ]._type == STANDARD_NODE
      || _offDiagonal[ interaction_idx ]._compressed )
    {
      continue;
    }

    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    baseData = extraData + interaction._extendedDataOffset;

    TRACE_ASSERT( interaction._extendedDataOffset >= 0,
                  "Invalid extended data offset" );

    // Apply this node's LDL permutation to the right hand side of
    // the interaction data
    //
    // If U is the schur complement in this interaction block and this
    // node has factorization P * L * D * L' * P', then
    // what we want for now is:
    //      U <-- U * P * L^{-T} * D^{-1}
    MATRIX::applyLDLPermutation( _LDL_data._pivotData, baseData,
                                 interaction._numExtendedRows, nCols,
                                 /* Right side, not transposed */
                                 false, false );
    
    MATRIX::triangularSolve( diagonalData, baseData,
                             interaction._numExtendedRows, nCols,
                             false, /* Right side */
                             true, /* Lower triangular matrix */
                             true, /* Transpose */
                             true /* Unit diagonal */ );

    // Here, we need workspaces for transposition and diagonal
    // inverse application
    {
      RealWorkspace          workspace(
                                workspaceManager,
                                IndexPair( nCols,
                                           interaction._numExtendedRows ) );

      MATRIX::transposeBLAS( workspace.workspaceData( 0 ), baseData,
                             interaction._numExtendedRows, nCols );

      applyLDLDiagonalInverse( workspace.workspaceData( 0 ),
                               nCols, interaction._numExtendedRows );

      // Copy back to node storage
      MATRIX::transposeBLAS( baseData, workspace.workspaceData( 0 ),
                             nCols, interaction._numExtendedRows );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Multiplication by this node's LDL diagonal matrix
//////////////////////////////////////////////////////////////////////
void Supernode::applyLDLDiagonal( const Real *input,
                                  int inputRows, int inputCols,
                                  Real *output,
                                  bool transpose ) const
{
  int                        nCols = numColumns();
  const IntArray            &pivotData = _LDL_data._pivotData;

  const Real                *blockData;
  const Real                *baseInputData;  
  Real                      *baseOutputData;
  int                        blockSize;

  int                        blockIndex_1x1 = 0;
  int                        blockIndex_2x2 = 0;

  TRACE_ASSERT(    ( transpose && nCols == inputCols )
                || ( !transpose && nCols == inputRows ) );
  TRACE_ASSERT( pivotData.size() == nCols );

  for ( int block_start = 0; block_start < pivotData.size(); block_start++ )
  {
    // If transposing, need to align with start of input column.
    // Otherwise, align with start of input row.
    baseInputData = transpose ? ( input + block_start )
                              : ( input + block_start * inputCols );

    baseOutputData = transpose ? ( output + block_start * inputRows )
                               : ( output + block_start * inputCols );

    if ( pivotData[ block_start ] < 0 ) {
      const LDL_Data::Block_2x2 &block
                      = _LDL_data._diagonalBlockEntries[ blockIndex_2x2 ];

      blockData = block.data;
      blockSize = 2;

      blockIndex_2x2 += 1;
      block_start += 1;
    }
    else {
      const LDL_Data::Block_1x1 &block
                      = _LDL_data._diagonalEntries[ blockIndex_1x1 ];
  
      blockData = &block;
      blockSize = 1;

      blockIndex_1x1 += 1;
    }

    // Multiply with the diagonal block
    MATRIX::gemm( blockData, baseInputData, baseOutputData,
                  blockSize, blockSize, /* Diagonal block dimensions */
                  /* Input block column dimensions */
                  transpose ? inputRows : blockSize,
                  transpose ? blockSize : inputCols,
                  /* Whether or not to do transposition */
                  false, transpose,
                  /* Overwrite */
                  1.0, 0.0,
                  /* Leading dimensions */
                  blockSize, inputCols );
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Multplication by this node's LDL diagonal inverse matrix
//////////////////////////////////////////////////////////////////////
void Supernode::applyLDLDiagonalInverse( const Real *input,
                                         int inputRows, int inputCols,
                                         Real *output ) const
{
  int                        nCols = numColumns();
  const IntArray            &pivotData = _LDL_data._pivotData;

  const Real                *blockData;
  const Real                *baseInputData;  
  Real                      *baseOutputData;

  int                        blockIndex_1x1 = 0;
  int                        blockIndex_2x2 = 0;

  TRACE_ASSERT( nCols == inputRows );
  TRACE_ASSERT( pivotData.size() == nCols );

  for ( int block_start = 0; block_start < pivotData.size(); block_start++ )
  {
    // If transposing, need to align with start of input column.
    // Otherwise, align with start of input row.
    baseInputData = input + block_start * inputCols;

    baseOutputData = output + block_start * inputCols;

    if ( pivotData[ block_start ] < 0 ) {
      const LDL_Data::Inverse_2x2 &block
                      = _LDL_data._diagonalBlockInverses[ blockIndex_2x2 ];

      blockData = block.data;

      blockIndex_2x2 += 1;
      block_start += 1;

      // Apply the inverse of the 2x2 system
      MATRIX::copy( baseOutputData, baseInputData, 2, inputCols );
      MATRIX::LUsolve( blockData, block.pivotData,
                       baseOutputData, 2, inputCols );
    }
    else {
      const LDL_Data::Inverse_1x1 &block
                      = _LDL_data._diagonalInverses[ blockIndex_1x1 ];
  
      blockData = &block;

      blockIndex_1x1 += 1;

      // Multiplication by the scalar inverse
      MATRIX::gemm( blockData, baseInputData, baseOutputData,
                    1, 1, /* Diagonal block dimensions */
                    /* Input block column dimensions */
                    1, inputCols,
                    false, false /* No transposition */ );
    }
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Multplication by this node's LDL diagonal inverse matrix
//
// Unfortunately, due to the structure of LU solves, we can't
// provide a transposed version of this, so applying to a transposed
// matrix will require some pre-processing
//
// Note that this is applied in place (ie. it overwrites the input)
//////////////////////////////////////////////////////////////////////
void Supernode::applyLDLDiagonalInverse( Real *input,
                                         int inputRows, int inputCols ) const
{
  int                        nCols = numColumns();
  const IntArray            &pivotData = _LDL_data._pivotData;

  const Real                *blockData;
  Real                      *baseData;

  int                        blockIndex_1x1 = 0;
  int                        blockIndex_2x2 = 0;

  TRACE_ASSERT( nCols == inputRows );
  TRACE_ASSERT( pivotData.size() == nCols );

  for ( int block_start = 0; block_start < pivotData.size(); block_start++ )
  {
    // Align with start of input row.
    baseData = input + block_start * inputCols;

    if ( pivotData[ block_start ] < 0 ) {
      const LDL_Data::Inverse_2x2 &block
                      = _LDL_data._diagonalBlockInverses[ blockIndex_2x2 ];

      blockData = block.data;

      blockIndex_2x2 += 1;
      block_start += 1;

      // Apply the inverse of the 2x2 system
      MATRIX::LUsolve( blockData, block.pivotData,
                       baseData, 2, inputCols );
    }
    else {
      const LDL_Data::Inverse_1x1 &block
                      = _LDL_data._diagonalInverses[ blockIndex_1x1 ];
  
      blockData = &block;

      blockIndex_1x1 += 1;

      MATRIX::scale( baseData, 1, inputCols, blockData[ 0 ] );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for diagonalSolve; performs a solve with this node's
// main diagonal assuming it has been compressed with extended variables
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalSolve_extendedVariable( Real *rhs, int nRHS,
                                                Real *extraData,
                                                bool transpose, bool left,
                                                int startBlock,
                                                int endBlock ) const
{
  Real                          *diagonalData;
  int                            nCols = numColumns();
  int                            rhsSize = nCols;

  // Check to see if we are only doing a solve with a subset of this
  // node's diagonal blocks
  if ( startBlock >= 0 && endBlock >= 0 )
  {
    TRACE_ASSERT( lowRankDiagonal() );
    TRACE_ASSERT( endBlock - startBlock >= 0 );
    TRACE_ASSERT( endBlock < _diagonalBlocks.size() );

    // Figure out the size of the right hand side, based on the
    // number of blocks being considered
    rhsSize = 0;

    for ( int block_idx = startBlock; block_idx <= endBlock; block_idx++ )
    {
      rhsSize += range_size( _diagonalBlocks[ block_idx ]._rowRange );
    }
  }
  else
  {
    startBlock = 0;
    endBlock = _diagonalBlocks.size() - 1;
  }

  if ( lowRankDiagonal() )
  {
    TRACE_ASSERT( _type == STANDARD_NODE,
                  "Only standard nodes can have sparsified diagonals" );

    Real                        *rhsData = rhs;
    int                          nRows;
    
    diagonalData = _data;
    for ( int block_idx = 0; block_idx < startBlock; block_idx++ )
    {
      const DenseBlock &block = _diagonalBlocks[ block_idx ];

      diagonalData += block.numRows() * block.numRows();
    }

    //for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ )
    for ( int block_idx = startBlock; block_idx <= endBlock; block_idx++ )
    {
      const DenseBlock &block = _diagonalBlocks[ block_idx ];

      nRows = range_size( block._rowRange );

      if ( left )
      {
        // We are applying the inverse to a block row
        MATRIX::triangularSolve( diagonalData, rhsData, nRows, nRHS,
                                 left, true /* lower triangular */, transpose );

        // Align with next block row
        rhsData += nRows * nRHS;
      }
      else
      {
        // We are applying the inverse to a block column
        MATRIX::triangularSolve( diagonalData, rhsData, nRHS, nRows,
                                 left, true /* lower triangular */, transpose,
                                 false /* Non-unit diagonal */,
                                 1.0 /* Alpha */,
                                 rhsSize /* LDA necessary for block column */ );

        // Align with next block column
        rhsData += nRows;
      }

      diagonalData += nRows * nRows;
    }
  }
  else
  {
    diagonalData = ( _type == STANDARD_NODE ) ? _data
                                              : extraData + _extendedDataOffset;

    MATRIX::triangularSolve( diagonalData, rhs,
                             left ? nCols : nRHS,
                             left ? nRHS : nCols,
                             left, true /* lower triangular */, transpose );
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for diagonalSolve; performs a solve with this node's
// main diagonal assuming it has been compressed in place
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalSolve_inPlace( Real *rhs, int nRHS,
                                       bool transpose, bool left ) const
{
  _nodeDiagonal.fullDiagonalSolve( rhs, nRHS, transpose, left );
}

//////////////////////////////////////////////////////////////////////
// Version of the above for left solves on vector data
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalSolve_inPlace( Real *rhs, bool transpose ) const
{
  _nodeDiagonal.fullDiagonalSolve( rhs, transpose );
}

//////////////////////////////////////////////////////////////////////
// Version of copyMatrixData called for standard nodes
//////////////////////////////////////////////////////////////////////
void Supernode::copyMatrixData_standard(
                              const IntArray &scatteredMap,
                              const IntArray &lowRankScatteredMap,
                              const SPARSE_MATRIX::SparseColumnMatrix &A )
{
  int                            columnOffset;
  int                            real_idx;

  int                            baseRow;

  // So we can skip the diagonal if it is being sparsified
  baseRow = lowRankDiagonal() ? _columnRange.second : -1;

  // Copy lower-triangular matrix data
  for ( int col_idx = _columnRange.first; col_idx <= _columnRange.second;
        col_idx++ )
  {
    for ( int row_idx = A._p[ col_idx ]; row_idx < A._p[ col_idx + 1 ];
          row_idx++ )
    {
      // Only copy lower-triangular part.
      if ( A._i[ row_idx ] < col_idx || A._i[ row_idx ] <= baseRow )
        continue;

      // We can also skip anything that belongs in a low rank block
      if ( lowRankScatteredMap[ A._i[ row_idx ] ] >= 0 )
        continue;

#if 0
      TRACE_ASSERT( scatteredMap[ A._i[ row_idx ] ] >= 0,
                    "ERROR: Invalid row index!!!!" );
#endif
      // This matrix entry belong to an implicit interaction
      if ( scatteredMap[ A._i[ row_idx ] ] < 0 ) {
        continue;
      }

      assign( scatteredMap[ A._i[ row_idx ] ], col_idx - _columnRange.first,
              A._x[ row_idx ], true /* add */ );
    }
  }

  // Handle the diagonal here, if were are introducing low rank blocks
  if ( lowRankDiagonal() )
  {
    copyMatrixData_sparseDiagonal( A );
  }
}

//////////////////////////////////////////////////////////////////////
// Version of copyMatrixData called for extended nodes
//////////////////////////////////////////////////////////////////////
void Supernode::copyMatrixData_extended( Real *extraData )
{
  // Here, we just need to assign the diagonal block to the
  // identity matrix
  int                            nCols = numColumns();
  Real                          *baseData = extraData + _extendedDataOffset;

  for ( int i = 0; i < nCols; i++ )
  {
    MATRIX::access( baseData, nCols, nCols, i, i ) = 1.0;
  }
}

//////////////////////////////////////////////////////////////////////
// For nodes with a sparsified diagonal
//////////////////////////////////////////////////////////////////////
void Supernode::copyMatrixData_sparseDiagonal(
                              const SPARSE_MATRIX::SparseColumnMatrix &A )
{
  Real                    *baseData = _data;

  int                      colStart, colEnd;
  int                      rowStart, rowEnd;

  int                      nRows, nCols;

  int                      dataRowIdx;

  // Copy matrix data in to the individual diagonal blocks
  for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ )
  {
    const DenseBlock      &block = _diagonalBlocks[ block_idx ];

    colStart = _columnRange.first + block._columnRange.first;
    colEnd = _columnRange.first + block._columnRange.second;

    rowStart = _columnRange.first + block._rowRange.first;
    rowEnd = _columnRange.first + block._rowRange.second;

    nRows = block.numRows();
    nCols = block.numColumns();

    for ( int col_idx = colStart; col_idx <= colEnd; col_idx++ )
    {
      for ( int row_idx = A._p[ col_idx ]; row_idx < A._p[ col_idx + 1 ];
            row_idx++ )
      {
        dataRowIdx = A._i[ row_idx ];

        if ( dataRowIdx < rowStart || dataRowIdx > rowEnd
          || dataRowIdx < col_idx )
        {
          continue;
        }

        MATRIX::access( baseData, nRows, nCols,
                        A._i[ row_idx ] - rowStart, /* row index */
                        col_idx - colStart /* col index */ ) = A._x[ row_idx ];
      }
    }

    baseData += nRows * nCols;
  }
}

//////////////////////////////////////////////////////////////////////
// Version of applyExtendedUpdate to be applied to standard nodes
//////////////////////////////////////////////////////////////////////
void Supernode::applyExtendedUpdate_standard( const Supernode &descendent,
                                              const PairArray &extendedIndices,
                                              Real *multWorkspace,
                                              Real *extraData,
                                              int start_idx,
                                              Real *expansionWorkspace )
{
  Real                          *baseData;
  const Real                    *descendentData;
  const Real                    *baseDescendentData;
  int                            descendentCols = descendent.numColumns();

  // Handle each off-diagonal interaction
  const SupernodeInteraction    &interaction
                      = descendent._offDiagonal[ start_idx ];

  TRACE_ASSERT( interaction._nodeID == _nodeID, "Invalid start index" );

  // Get the data for the interaction between descendent and this node
  descendentData
    = descendent._data + interaction._dataRowList[ 0 ] * descendentCols;

  for ( int interaction_idx = 0; interaction_idx < extendedIndices.size();
        interaction_idx++ )
  {
    const IndexPair             &p = extendedIndices[ interaction_idx ];

    const SupernodeInteraction  &descInteraction
                        = descendent._offDiagonal[ p.first ];

    const SupernodeInteraction  &thisInteraction = _offDiagonal[ p.second ];

    TRACE_ASSERT( !thisInteraction._compressed,
                  "Can't accumulate to a scattered interaction" );

    baseData = extraData + thisInteraction._extendedDataOffset;

    baseDescendentData = extraData + descInteraction._extendedDataOffset;

    // If this node is in scattered form, we must form a full
    // representation of the matrix
    if ( descInteraction._compressed )
    {
      TRACE_ASSERT( descInteraction._compressedColumnList.size() > 0,
                    "Empty scattered column list" );

      MATRIX::clear( expansionWorkspace, descInteraction._numExtendedRows,
                     descendentCols );

      // Scatter stored interaction data
      MATRIX::scatterColumns( baseDescendentData, expansionWorkspace,
                              descInteraction._compressedColumnList,
                              // Number of rows in this interaction
                              descInteraction._numExtendedRows,
                              descInteraction._compressedColumnList.size(),
                              descendentCols );

      // Since the interaction was stored in scattered format, we
      // only stored the Schur complement, not the factor sub-matrix.
      // Apply the inverse transpose of the diagonal factor to the right
      // of this data to form the factor contribution.
      descendent.diagonalSolve( expansionWorkspace,
                                descInteraction._numExtendedRows, extraData,
                                true, /* Transpose */
                                false /* Right side solve */ );

      baseDescendentData = expansionWorkspace;
    }

    // interaction may not store all of the rows covered by this
    // node, so we need to store the product in a workspace, then
    // scatter it back to thisInteraction's extra data.
    MATRIX::gemm( baseDescendentData,
                  descendentData, multWorkspace,
                  descInteraction._numExtendedRows, descendentCols,
                  interaction._rowList.size(), descendentCols,
                  false /* Don't transpose the first matrix */,
                  true /* Transpose the second matrix */ );

    // Scatter the product back to the result matrix
    ScatterLowRankTransUpdate( start_idx, descendent, *this,
                               multWorkspace, baseData,
                               descInteraction._numExtendedRows );
  }
}

//////////////////////////////////////////////////////////////////////
// Version of applyExtendedUpdate to be applied to standard nodes
//////////////////////////////////////////////////////////////////////
void Supernode::applyExtendedUpdate_standard(
                                    const Supernode &descendent,
                                    const PairArray &extendedIndices,
                                    Real *extraData,
                                    int start_idx,
                                    WorkspaceManager<Real> &workspaceManager )
{
  Real                          *baseData;
  const Real                    *descendentData;
  const Real                    *baseDescendentData;
  int                            descendentCols = descendent.numColumns();

  // For specifying data sizes
  PairArray                      dataSizes( 2 );

  Real                          *expansionWorkspace;
  Real                          *multWorkspace;

  if ( extendedIndices.size() == 0 ) {
    return;
  }

  // Handle each off-diagonal interaction
  const SupernodeInteraction    &interaction
                      = descendent._offDiagonal[ start_idx ];

  TRACE_ASSERT( interaction._nodeID == _nodeID, "Invalid start index" );

  // Get the data for the interaction between descendent and this node
  descendentData
    = descendent._data + interaction._dataRowList[ 0 ] * descendentCols;

  for ( int interaction_idx = 0; interaction_idx < extendedIndices.size();
        interaction_idx++ )
  {
    const IndexPair             &p = extendedIndices[ interaction_idx ];

    const SupernodeInteraction  &descInteraction
                        = descendent._offDiagonal[ p.first ];

    const SupernodeInteraction  &thisInteraction = _offDiagonal[ p.second ];

    TRACE_ASSERT( !thisInteraction._compressed,
                  "Can't accumulate to a scattered interaction" );

    baseData = extraData + thisInteraction._extendedDataOffset;

    baseDescendentData = extraData + descInteraction._extendedDataOffset;

    // Create workspaces
    if ( descInteraction._compressed )
    {
      // Set workspace data sizes
      dataSizes[ 0 ] = IndexPair( descInteraction._numExtendedRows,
                                 descendentCols );
    }
    else
    {
      dataSizes[ 0 ] = IndexPair( 0, 0 );
    }
    dataSizes[ 1 ] = IndexPair( descInteraction._numExtendedRows,
                               interaction._rowList.size() );

    RealWorkspace            workspace( workspaceManager, dataSizes );

    // If this node is in scattered form, we must form a full
    // representation of the matrix
    if ( descInteraction._compressed )
    {
      TRACE_ASSERT( descInteraction._compressedColumnList.size() > 0,
                    "Empty scattered column list" );

      expansionWorkspace = workspace.workspaceData( 0 );

      MATRIX::clear( expansionWorkspace, descInteraction._numExtendedRows,
                     descendentCols );

      // Scatter stored interaction data
      MATRIX::scatterColumns( baseDescendentData, expansionWorkspace,
                              descInteraction._compressedColumnList,
                              // Number of rows in this interaction
                              descInteraction._numExtendedRows,
                              descInteraction._compressedColumnList.size(),
                              descendentCols );

      // Since the interaction was stored in scattered format, we
      // only stored the Schur complement, not the factor sub-matrix.
      // Apply the inverse transpose of the diagonal factor to the right
      // of this data to form the factor contribution.
      descendent.diagonalSolve( expansionWorkspace,
                                descInteraction._numExtendedRows, extraData,
                                true, /* Transpose */
                                false /* Right side solve */ );

      baseDescendentData = expansionWorkspace;
    }

    multWorkspace = workspace.workspaceData( 1 );

    // interaction may not store all of the rows covered by this
    // node, so we need to store the product in a workspace, then
    // scatter it back to thisInteraction's extra data.
    MATRIX::gemm( baseDescendentData,
                  descendentData, multWorkspace,
                  descInteraction._numExtendedRows, descendentCols,
                  interaction._rowList.size(), descendentCols,
                  false /* Don't transpose the first matrix */,
                  true /* Transpose the second matrix */ );

    // Scatter the product back to the result matrix
    ScatterLowRankTransUpdate( start_idx, descendent, *this,
                               multWorkspace, baseData,
                               descInteraction._numExtendedRows );
  }
}

//////////////////////////////////////////////////////////////////////
// Version of applyExtendedUpdate to be applied to extended nodes
//////////////////////////////////////////////////////////////////////
void Supernode::applyExtendedUpdate_extended( const Supernode &descendent,
                                              const PairArray &extendedIndices,
                                              Real *extraData,
                                              int start_idx,
                                              Real *expansionWorkspace,
                                              Timer *standardTimer,
                                              Timer *extendedTimer )
{
  Real                          *baseData;
  const Real                    *descendentData;
  const Real                    *baseDescendentData;
  int                            nCols = numColumns();
  int                            descendentCols = descendent.numColumns();

#if 0
  if ( descendent._type == STANDARD_NODE && extendedTimer )
  {
    extendedTimer->tick();
  }
#endif

  // FIXME: skip standard nodes
  if ( descendent.type() == STANDARD_NODE )
  {
#if 0
    if ( extendedTimer )
    {
      extendedTimer->tock();
    }
#endif
    return;
  }

  // Form the update matrix for the diagonal block
  {
    const SupernodeInteraction  &interaction
                        = descendent._offDiagonal[ start_idx ];

    TRACE_ASSERT( interaction._nodeID == _nodeID, "Invalid start index" );
    TRACE_ASSERT( interaction._numExtendedRows == nCols,
                  "Invalid interaction size" );

    baseData = extraData + _extendedDataOffset;

    descendentData = extraData + interaction._extendedDataOffset;

    // FIXME: Treat everything as compressed for now
    if ( interaction._compressed || descendent._type == STANDARD_NODE )
    {
#if 0
      if ( descendent._type == STANDARD_NODE && extendedTimer )
      {
        extendedTimer->tock();
      }
      if ( descendent._type == STANDARD_NODE && standardTimer )
      {
        standardTimer->tick();
      }
#endif

      // FIXME: debugging
      TRACE_ASSERT( NULL, "Shouldn't be here!!!!" );

      // Scatter stored interaction data
      //
      // FIXME: "scatter" even for uncompressed interactions
      if ( interaction._compressed )
      {
        TRACE_ASSERT( interaction._compressedColumnList.size() > 0,
                      "Empty scattered column list" );

        MATRIX::clear( expansionWorkspace, nCols, descendentCols );

        MATRIX::scatterColumns( descendentData, expansionWorkspace,
                                interaction._compressedColumnList,
                                nCols, // Number of rows in this interaction
                                interaction._compressedColumnList.size(),
                                descendentCols );
      }
      else
      {
        MATRIX::copy( expansionWorkspace, descendentData,
                      nCols, descendentCols );
      }

      // Since the interaction was stored in scattered format, we
      // only stored the Schur complement, not the factor sub-matrix.
      // Apply the inverse transpose of the diagonal factor to the right
      // of this data to form the factor contribution.
      descendent.diagonalSolve( expansionWorkspace, nCols, extraData,
                                true, /* Transpose */
                                false /* Right side solve */ );

      descendentData = expansionWorkspace;

      expansionWorkspace += nCols * descendentCols;

#if 0
      if ( descendent._type == STANDARD_NODE && standardTimer )
      {
        standardTimer->tock();
      }
      if ( descendent._type == STANDARD_NODE && extendedTimer )
      {
        extendedTimer->tick();
      }
#endif
    }

    MATRIX::syrk( descendentData, baseData, nCols, descendentCols,
                  false /* Do not transpose */,
                  -1.0 /* Subtract from the total */, 1.0 /* beta = 1.0 */ );
  }

  // Handle each off diagonal interaction
  for ( int interaction_idx = 0; interaction_idx < extendedIndices.size();
        interaction_idx++ )
  {
    const IndexPair             &p = extendedIndices[ interaction_idx ];

    const SupernodeInteraction  &interaction
                        = descendent._offDiagonal[ p.first ];

    const SupernodeInteraction  &thisInteraction = _offDiagonal[ p.second ];

    TRACE_ASSERT( !thisInteraction._compressed,
                  "Can't accumulate to a scattered interaction" );

    baseData = extraData + thisInteraction._extendedDataOffset;
    baseDescendentData = extraData + interaction._extendedDataOffset;

    // If this node is in scattered form, we must form a full
    // representation of the matrix
    //
    // FIXME: treat everything as compressed for now
    if ( interaction._compressed || descendent._type == STANDARD_NODE )
    {
      TRACE_ASSERT( NULL, "Should never get here" );
#if 0
      if ( descendent._type == STANDARD_NODE && extendedTimer )
      {
        extendedTimer->tock();
      }
      if ( descendent._type == STANDARD_NODE && standardTimer )
      {
        standardTimer->tick();
      }
#endif

      // Scatter stored interaction data
      // FIXME: scatter even for uncompressed interactions
      if ( interaction._compressed )
      {
        TRACE_ASSERT( interaction._compressedColumnList.size() > 0,
                      "Empty scattered column list" );

        TRACE_ASSERT( interaction._compressedColumnList.size() <= descendentCols,
                      "Invalid compressed column list" );

        MATRIX::clear( expansionWorkspace, interaction._numExtendedRows,
                       descendentCols );

        MATRIX::scatterColumns( baseDescendentData, expansionWorkspace,
                                interaction._compressedColumnList,
                                interaction._numExtendedRows, // Num rows in
                                                              // this interaction
                                interaction._compressedColumnList.size(),
                                descendentCols );
      }
      else
      {
        MATRIX::copy( expansionWorkspace, baseDescendentData,
                      interaction._numExtendedRows, descendentCols );
      }

      // Since the interaction was stored in scattered format, we
      // only stored the Schur complement, not the factor sub-matrix.
      // Apply the inverse transpose of the diagonal factor to the right
      // of this data to form the factor contribution.
      descendent.diagonalSolve( expansionWorkspace,
                                interaction._numExtendedRows, extraData,
                                true, /* Transpose */
                                false /* Right side solve */ );

      baseDescendentData = expansionWorkspace;

#if 0
      if ( descendent._type == STANDARD_NODE && standardTimer )
      {
        standardTimer->tock();
      }
      if ( descendent._type == STANDARD_NODE && extendedTimer )
      {
        extendedTimer->tick();
      }
#endif
    }

    MATRIX::gemm( baseDescendentData, descendentData, baseData,
                  interaction._numExtendedRows, descendentCols,
                  nCols, descendentCols,
                  false /* Don't transpose the first matrix */,
                  true /* Transpose the second matrix */,
                  -1.0, 1.0 /* Subtract from the total */ );
  }

#if 0
  if ( descendent._type == STANDARD_NODE && extendedTimer )
  {
    extendedTimer->tock();
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Extended update for a descendent whose factor is stored as
// an LDL^{T} factorization
//////////////////////////////////////////////////////////////////////
void Supernode::applyExtendedUpdate_extendedLDL(
                                  const Supernode &descendent,
                                  const PairArray &extendedIndices,
                                  Real *extraData,
                                  int start_idx,
                                  WorkspaceManager<Real> &workspaceManager,
                                  Timer *standardTimer,
                                  Timer *extendedTimer )
{
  Real                      *baseData;
  const Real                *descendentData;
  int                        nCols = numColumns();
  int                        descendentCols = descendent.numColumns();

  if ( descendent.type() == STANDARD_NODE ) {
    return;
  }

  {
    const SupernodeInteraction  &interaction
                        = descendent._offDiagonal[ start_idx ];

    TRACE_ASSERT( interaction._nodeID == _nodeID, "Invalid start index" );
    TRACE_ASSERT( interaction._numExtendedRows == nCols,
                  "Invalie interaction size" );

    descendentData = extraData + interaction._extendedDataOffset;
  }

  // Suppose the factorization for descendent is given by
  //  [ P L D L' P' ]
  //  [      .      ]
  //  [      .      ]
  //  [     Z_s     ]
  //  [     Z_l     ]
  // where Z_s is the interaction between descendent and this.
  // The update we want is:
  //    [ Z_s ] * D * Z_s'
  //    [ Z_l ]
  //
  // Start by forming D * Z_s', since we will need this throughout
  RealWorkspace              workspace( workspaceManager,
                                        IndexPair( descendentCols,
                                                   nCols ) );

  descendent.applyLDLDiagonal( descendentData,
                               nCols, descendentCols,
                               workspace.workspaceData( 0 ),
                               true /* transpose Z_s */ );

  // Update this node's diagonal
  baseData = extraData + _extendedDataOffset;
  MATRIX::gemm( descendentData, workspace.workspaceData( 0 ), baseData,
                nCols, descendentCols,
                descendentCols, nCols,
                /* No transposition necesssary; should have been handled
                   during diagonal application */
                false, false,
                -1.0, 1.0 /* Subtract from the total */ );

  // Handle each off diagonal interaction
  for ( int interaction_idx = 0; interaction_idx < extendedIndices.size();
        interaction_idx++ )
  {
    const IndexPair             &p = extendedIndices[ interaction_idx ];

    const SupernodeInteraction  &interaction
                        = descendent._offDiagonal[ p.first ];

    const SupernodeInteraction  &thisInteraction = _offDiagonal[ p.second ];

    TRACE_ASSERT( !thisInteraction._compressed,
                  "Can't accumulate to a scattered interaction" );

    baseData = extraData + thisInteraction._extendedDataOffset;
    descendentData = extraData + interaction._extendedDataOffset;

    MATRIX::gemm( descendentData, workspace.workspaceData( 0 ), baseData,
                  interaction._numExtendedRows, descendentCols,
                  descendentCols, nCols,
                  /* No transposition necesssary; should have been handled
                     during diagonal application */
                  false, false,
                  -1.0, 1.0 /* Subtract from the total */ );
  }
}

//////////////////////////////////////////////////////////////////////
// Identifies the row set associated with the original system in each
// of this node's interactions
//////////////////////////////////////////////////////////////////////
void Supernode::buildBaseSystemRowSets(
                           const vector<Supernode> &nodes,
                           vector<set<int> > &rowSets,
                           const SPARSE_MATRIX::SparseColumnMatrix &A ) const
{
  rowSets.resize( _offDiagonal.size() );

  for ( int interaction_idx = 0; interaction_idx < _offDiagonal.size();
        interaction_idx++ )
  {
    const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];
    const Supernode             &ancestor = nodes[ interaction._nodeID ];

    set<int>                    &rowSet = rowSets[ interaction_idx ];

    rowSet.clear();

    for ( int col_idx = _columnRange.first; col_idx <= _columnRange.second;
          col_idx++ )
    {
      for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
            row_ptr++ )
      {
        // Check to see if the matrix has any entries in this interaction's
        // submatrix
        if ( A._i[ row_ptr ] >= ancestor._columnRange.first
          && A._i[ row_ptr ] <= ancestor._columnRange.second )
        {
          // Add the row index relative to the interaction
          rowSet.insert( A._i[ row_ptr ] - ancestor._columnRange.first );
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Based on this node's diagonal set, figure out where the
// interaction data should be placed
//////////////////////////////////////////////////////////////////////
void Supernode::setDiagonalSize()
{
  long int                   numEntries = 0;
  int                        nRows;
  int                        nCols = numColumns();
  int                        padding;

  // Count the number of entries on the diagonal
  for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ )
  {
    const DenseBlock        &block = _diagonalBlocks[ block_idx ];

    TRACE_ASSERT( block._rowRange == block._columnRange,
                  "Diagonal blocks must be square" );

    nRows = range_size( block._rowRange );

    numEntries += nRows * nRows;
  }

  // Pad this with extra space so that the remaining interaction data
  // is row-aligned.
  if ( numEntries % nCols > 0 ) {
    padding = nCols - ( numEntries % nCols );
  } else {
    padding = 0;
  }
  numEntries += padding;

  TRACE_ASSERT( numEntries % nCols == 0, "Not row-aligned" );

  _firstOffDiagonalRow = numEntries / nCols;

#if 0
  printf( "Node %d: first off diagonal row = %d\n",
          _nodeID, _firstOffDiagonalRow );
#endif
}

//////////////////////////////////////////////////////////////////////
// Recursively build the set of dense blocks to be decomposed
// in the diagonal, based on the set of diagonal blocks
//////////////////////////////////////////////////////////////////////
void Supernode::buildLowRankDiagonalBlocks( int start_idx, int end_idx,
                                            int parent )
{
  int                        rangeSize = end_idx - start_idx;
  int                        mid_idx_low, mid_idx_high;
  IndexRange                 rowRange, columnRange;
  int                        next_parent;

  // Base case: we can only build interactions between two or
  // more blocks
  if ( rangeSize <= 1 )
  {
    return;
  }

  TRACE_ASSERT( rangeSize % 2 == 0, "Number of diagonal blocks must be"
                                    " a power of 2" );

  mid_idx_high = rangeSize / 2 + start_idx;
  mid_idx_low = mid_idx_high - 1;

  rowRange.first = _diagonalBlocks[ mid_idx_high ]._rowRange.first;
  rowRange.second = _diagonalBlocks[ end_idx - 1 ]._rowRange.second;
  
  columnRange.first = _diagonalBlocks[ start_idx ]._columnRange.first;
  columnRange.second = _diagonalBlocks[ mid_idx_low ]._columnRange.second;

  _diagonalLowRankBlockParents.push_back( parent );

  parent = _diagonalLowRankBlocks.size();

  _diagonalLowRankBlocks.push_back( DenseBlock( rowRange, columnRange ) );

  // Work out previous diagonal contributions for each parent
  _previousDiagonalContributions.push_back( vector<pair<int, DenseBlock> >() );

  vector<pair<int, DenseBlock> > &previousContributions
                                  = _previousDiagonalContributions.back();

  const DenseBlock          &thisBlock = _diagonalLowRankBlocks.back();

  // Append to the list of contributions for each diagonal block
  for ( int diag_block_idx = start_idx; diag_block_idx < end_idx;
        diag_block_idx++ )
  {
    const DenseBlock        &diagBlock = _diagonalBlocks[ diag_block_idx ];

    int                      rowStart;
    int                      rowEnd;

    rowStart = diagBlock._rowRange.first - columnRange.first;
    rowEnd = diagBlock._rowRange.second - columnRange.first;

    vector<pair<int, IndexRange> > &mainDiagonalContributions
                              = _mainDiagonalContributions[ diag_block_idx ];

    mainDiagonalContributions.push_back(
      pair<int, IndexRange>( parent, IndexRange( rowStart, rowEnd ) ) );
  }

  // Walk through the ancestor branch and determine which subblock of
  // each ancestor is needed to form this block
  for ( next_parent = _diagonalLowRankBlockParents.back(); next_parent >= 0;
        next_parent = _diagonalLowRankBlockParents[ next_parent ] )
  {
    const DenseBlock        &ancestor = _diagonalLowRankBlocks[ next_parent ];

    int                      rowStart, rowEnd;
    int                      colStart, colEnd;

    if ( thisBlock._rowRange.first < ancestor._rowRange.first )
    {
      // This is a bit tricky to explain...
      //
      // We are extracting from the "column range" of the ancestor, which
      // comes first in the interaction
      rowStart = range_size( ancestor._columnRange );
      //rowStart -= ancestor._rowRange.first - thisBlock._rowRange.first + 1;
      rowStart -= ancestor._rowRange.first - thisBlock._rowRange.first;

      colStart = thisBlock._columnRange.first - ancestor._columnRange.first;
    }
    else
    {
      // This is a bit tricky to explain...
      //
      // We are extracting from the "row range" of the ancestor, which
      // comes second in the interaction
      rowStart = thisBlock._rowRange.first - ancestor._rowRange.first;
      rowStart += range_size( ancestor._columnRange );

      colStart = thisBlock._columnRange.first - ancestor._columnRange.first;
    }

    rowEnd = rowStart + range_size( thisBlock._rowRange ) - 1;
    colEnd = colStart + range_size( thisBlock._columnRange ) - 1;

    previousContributions.push_back(
                    pair<int, DenseBlock>( next_parent,
                                           DenseBlock( rowStart, rowEnd,
                                                       colStart, colEnd ) ) );
  }

  buildLowRankDiagonalBlocks( start_idx, mid_idx_high, parent );
  buildLowRankDiagonalBlocks( mid_idx_high, end_idx, parent );
}

#if 0
//////////////////////////////////////////////////////////////////////
// For a given piece of a low rank decomposition, add its outer
// product to the diagonal of this node
//
// eg. diagonal += U' * U
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContribution( const Real *U, int rank,
                                                bool transpose )
{
  if ( !lowRankDiagonal() )
  {
    // Add to the full diagonal
    MATRIX::syrk( U, _data, numColumns(), rank,
                  transpose, 1.0, 1.0 /* Add to the diagonal */ );
  }
  else
  {
    Real                    *baseData = _data;
    int                      nRows;
    int                      nCols = numColumns();

    for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ )
    {
      const DenseBlock      &block = _diagonalBlocks[ block_idx ];

      nRows = range_size( block._rowRange );

      MATRIX::syrk( U, baseData, nRows, rank,
                    transpose, 1.0, 1.0, /* Add to the diagonal */
                    // We need to explicitly set the leading dimension,
                    // since we are only using submatrices here
                    transpose ? nCols : rank );

      // Align with the next block, whose position will depend on
      // whether or not we are storing the transpose
      if ( transpose )
      {
        U += nRows;
      }
      else
      {
        U += nRows * rank;
      }

      baseData += nRows * nRows;
    }
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// For a given piece of a low rank decomposition, add its outer
// product to the diagonal of this node
//
// eg. diagonal += U' * U
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContribution(
                              const Real *U, int rank,
                              Real *multWorkspace,
                              const IntArray *compressedColumnList,
                              bool transpose )
{
  if ( compressedColumnList )
  {
    addLowRankDiagonalContribution_compressed( U, rank, multWorkspace,
                                               *compressedColumnList,
                                               transpose );
  }
  else
  {
    addLowRankDiagonalContribution_uncompressed( U, rank, transpose );
  }
}

//////////////////////////////////////////////////////////////////////
// For a given piece of a low rank decomposition, add its outer
// product to the diagonal of this node
//
// eg. diagonal += U' * U
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContribution(
                              const Real *U, int rank,
                              WorkspaceManager<Real> &workspaceManager,
                              const IntArray *compressedColumnList,
                              bool transpose )
{
  if ( compressedColumnList )
  {
    addLowRankDiagonalContribution_compressed( U, rank, workspaceManager,
                                               *compressedColumnList,
                                               transpose );
  }
  else
  {
    addLowRankDiagonalContribution_uncompressed( U, rank, transpose );
  }
}

//////////////////////////////////////////////////////////////////////
// Versions of the above for compressed matrices
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContribution_compressed(
                              const Real *U, int rank,
                              Real *multWorkspace,
                              const IntArray &compressedColumnList,
                              bool transpose )
{
  if ( !lowRankDiagonal() )
  {
    int                    nCols = compressedColumnList.size();
    IndexRange             fullRange( 0, compressedColumnList.size() - 1 );

    // Form the compressed product in a workspace
    MATRIX::syrk( U, multWorkspace, nCols, rank, transpose );

    // Scatter to the diagonal
    MATRIX::scatterAddSymmetricMatrix( multWorkspace, _data,
                                       // Rows/columns to scatter to
                                       compressedColumnList,
                                       // Full size of matrix
                                       numColumns(),
                                       // Use full range of
                                       // compressedColumnList
                                       fullRange,
                                       // No offset, and scale by 1.0
                                       0, 1.0 );
  }
  else
  {
    IndexRange               interactionRange( 0, 0 );
    int                      current_block_idx;
    int                      diag_block_idx = 0;

    Real                    *baseData = _data;
    const Real              *interactionData = U;

    while ( interactionRange.first < compressedColumnList.size() )
    {
      // Find the block associated with the current sub-range
      current_block_idx = _diagonalMap[ compressedColumnList[
                                                interactionRange.first ] ];

      // Move the data pointer up to the current block
      while ( diag_block_idx < current_block_idx )
      {
        int                  nRows;
        
        nRows = _diagonalBlocks[ diag_block_idx ].numRows();
        baseData += nRows * nRows;

        diag_block_idx++;
      }

      const DenseBlock      &block = _diagonalBlocks[ current_block_idx ];

      // Figure out how big this range is
      while ( interactionRange.second < compressedColumnList.size()
           && _diagonalMap[ compressedColumnList[ interactionRange.second ] ]
              == current_block_idx )
      {
        interactionRange.second += 1;
      }
      interactionRange.second -= 1;
      TRACE_ASSERT(
        _diagonalMap[ compressedColumnList[ interactionRange.second ] ]
          == current_block_idx );

      int                    nCols;

      nCols = range_size( interactionRange );

      // Build the update matrix
      MATRIX::syrk( interactionData, multWorkspace, nCols, rank,
                    transpose, 1.0, 0.0, /* Just get the matrix */
                    // We need to explicitly set the leading dimension
                    // since we are only using submatrices here
                    transpose ? compressedColumnList.size() : rank );

      // Scatter to the diagonal block
      MATRIX::scatterAddSymmetricMatrix( multWorkspace, baseData,
                                         // Rows/columns to scatter to
                                         compressedColumnList,
                                         // Full size of diagonal block
                                         block.numColumns(),
                                         // Sub-range of compressedColumnList
                                         // to use
                                         interactionRange,
                                         // Offset by start of diagonal block
                                         block._rowRange.first,
                                         // Scale by 1.0
                                         1.0 );

      // Move to the next block
      if ( transpose )
      {
        interactionData += nCols;
      }
      else
      {
        interactionData += nCols * rank;
      }

      interactionRange.first = interactionRange.second + 1;
      interactionRange.second = interactionRange.first;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Versions of the above for compressed matrices
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContribution_compressed(
                                    const Real *U, int rank,
                                    WorkspaceManager<Real> &workspaceManager,
                                    const IntArray &compressedColumnList,
                                    bool transpose )
{
  if ( !lowRankDiagonal() )
  {
    int                    nCols = compressedColumnList.size();
    IndexRange             fullRange( 0, compressedColumnList.size() - 1 );

    PairArray              dataSizes( 1 );
    Real                  *multWorkspace;

    // Get a workspace for multiplication
    dataSizes[ 0 ] = IndexPair( nCols, nCols );

    RealWorkspace          workspace( workspaceManager, dataSizes );

    multWorkspace = workspace.workspaceData( 0 );

    // Form the compressed product in a workspace
    MATRIX::syrk( U, multWorkspace, nCols, rank, transpose );

    // Scatter to the diagonal
    MATRIX::scatterAddSymmetricMatrix( multWorkspace, _data,
                                       // Rows/columns to scatter to
                                       compressedColumnList,
                                       // Full size of matrix
                                       numColumns(),
                                       // Use full range of
                                       // compressedColumnList
                                       fullRange,
                                       // No offset, and scale by 1.0
                                       0, 1.0 );
  }
  else
  {
    IndexRange               interactionRange( 0, 0 );
    int                      current_block_idx;
    int                      diag_block_idx = 0;

    Real                    *baseData = _data;
    const Real              *interactionData = U;

    PairArray                dataSizes( 1 );
    Real                    *multWorkspace;

    while ( interactionRange.first < compressedColumnList.size() )
    {
      // Find the block associated with the current sub-range
      current_block_idx = _diagonalMap[ compressedColumnList[
                                                interactionRange.first ] ];

      // Move the data pointer up to the current block
      while ( diag_block_idx < current_block_idx )
      {
        int                  nRows;
        
        nRows = _diagonalBlocks[ diag_block_idx ].numRows();
        baseData += nRows * nRows;

        diag_block_idx++;
      }

      const DenseBlock      &block = _diagonalBlocks[ current_block_idx ];

      // Figure out how big this range is
      while ( interactionRange.second < compressedColumnList.size()
           && _diagonalMap[ compressedColumnList[ interactionRange.second ] ]
              == current_block_idx )
      {
        interactionRange.second += 1;
      }
      interactionRange.second -= 1;
      TRACE_ASSERT(
        _diagonalMap[ compressedColumnList[ interactionRange.second ] ]
          == current_block_idx );

      int                    nCols;

      nCols = range_size( interactionRange );

      // Get a workspace for multiplication
      dataSizes[ 0 ] = IndexPair( nCols, nCols );

      RealWorkspace          workspace( workspaceManager, dataSizes );

      multWorkspace = workspace.workspaceData( 0 );

      // Build the update matrix
      MATRIX::syrk( interactionData, multWorkspace, nCols, rank,
                    transpose, 1.0, 0.0, /* Just get the matrix */
                    // We need to explicitly set the leading dimension
                    // since we are only using submatrices here
                    transpose ? compressedColumnList.size() : rank );

      // Scatter to the diagonal block
      MATRIX::scatterAddSymmetricMatrix( multWorkspace, baseData,
                                         // Rows/columns to scatter to
                                         compressedColumnList,
                                         // Full size of diagonal block
                                         block.numColumns(),
                                         // Sub-range of compressedColumnList
                                         // to use
                                         interactionRange,
                                         // Offset by start of diagonal block
                                         block._rowRange.first,
                                         // Scale by 1.0
                                         1.0 );

      // Move to the next block
      if ( transpose )
      {
        interactionData += nCols;
      }
      else
      {
        interactionData += nCols * rank;
      }

      interactionRange.first = interactionRange.second + 1;
      interactionRange.second = interactionRange.first;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Version for uncompressed matrices
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContribution_uncompressed(
                              const Real *U, int rank,
                              bool transpose )
{
  if ( !lowRankDiagonal() )
  {
    // Add to the full diagonal
    MATRIX::syrk( U, _data, numColumns(), rank,
                  transpose, 1.0, 1.0 /* Add to the diagonal */ );
  }
  else
  {
    Real                    *baseData = _data;
    int                      nRows;
    int                      nCols = numColumns();

    for ( int block_idx = 0; block_idx < _diagonalBlocks.size(); block_idx++ )
    {
      const DenseBlock      &block = _diagonalBlocks[ block_idx ];

      nRows = range_size( block._rowRange );

      MATRIX::syrk( U, baseData, nRows, rank,
                    transpose, 1.0, 1.0, /* Add to the diagonal */
                    // We need to explicitly set the leading dimension,
                    // since we are only using submatrices here
                    transpose ? nCols : rank );

      // Align with the next block, whose position will depend on
      // whether or not we are storing the transpose
      if ( transpose )
      {
        U += nRows;
      }
      else
      {
        U += nRows * rank;
      }

      baseData += nRows * nRows;
    }
  }

}

//////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////
void Supernode::applyInverseToInteraction( int interaction_idx,
                                           Real *extraData,
                                           Real *expansionWorkspace )
{
  const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

  int                          nCols = numColumns();

  int                          startBlock, endBlock;

  if ( interaction._lowRankDiagonalInteraction )
  {
    int                        numCompressedCols;
    
    TRACE_ASSERT( lowRankDiagonal(),
                  "Node does not have low rank diagonal blocks" );

    // In this case, we only want to apply inverse to the desired range
    startBlock = _diagonalMap[ interaction._compressedColumnRange.first ];
    endBlock = _diagonalMap[ interaction._compressedColumnRange.second ];

    TRACE_ASSERT( endBlock - startBlock >= 0 );

    numCompressedCols = range_size( interaction._compressedColumnRange );

    // Copy to the expansion workspace - no decompression required
    MATRIX::copy( expansionWorkspace,
                  // Original interaction data
                  extraData + interaction._extendedDataOffset,
                  // Matrix size
                  interaction._numExtendedRows, numCompressedCols );
  }
  else
  {
    startBlock = -1;
    endBlock = -1;

    if ( interaction._compressed )
    {
      // Expand the compressed interaction
      MATRIX::clear( expansionWorkspace, interaction._numExtendedRows, nCols );

      // Uncompress the interaction data
      MATRIX::scatterColumns( // Original matrix data
                              extraData + interaction._extendedDataOffset,
                              expansionWorkspace,
                              interaction._compressedColumnList,
                              interaction._numExtendedRows,
                              // Number of columns in the input
                              interaction._compressedColumnList.size(),
                              // Number of columns in the output
                              nCols );
    }
    else
    {
      MATRIX::copy( expansionWorkspace,
                    // Original interaction data
                    extraData + interaction._extendedDataOffset,
                    // Matrix size
                    interaction._numExtendedRows, nCols );
    }
  }

  // We now have the fully uncompressed interaction, so we can just
  // apply the inverse directly
  //
  // We only solve over the selected diagonal block range (if this node
  // has low rank diagonal blocks)
  //
  // Two successive triangular solves:
  //    ie. if S is the block stored in this interaction, and L is the
  //    diagonal factor block for this node, then we form
  //      S * L^{-T} * L^{-1}
  diagonalSolve( expansionWorkspace, interaction._numExtendedRows,
                 extraData,
                 true, /* transposed solve */
                 false, /* right side solve */
                 startBlock, endBlock );

  diagonalSolve( expansionWorkspace, interaction._numExtendedRows,
                 extraData,
                 false, /* not transposed */
                 false, /* right side solve */
                 startBlock, endBlock );

}

//////////////////////////////////////////////////////////////////////
// Multiplies the matrix stored for the given interaction with the
// transpose of the provided matrix, which may only be defined over
// a particular column range (stored in expansionWorkspace).
//
// Subtract the result from the provided output matrix
//
// This is a helper function for addExtendedSchurComplementContribution
//////////////////////////////////////////////////////////////////////
void Supernode::multiplyInteraction( int interaction_idx,
                                     const Real *expansionWorkspace,
                                     int numRows, const IndexRange &columnRange,
                                     Real *extraData, Real *multWorkspace,
                                     Real *outputMatrix )
{
  const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

  const Real                  *interactionData;
  const Real                  *multiplyData;
  int                          numMultiplyColumns;
  int                          numOriginalColumns;
  int                          nCols = numColumns();

  int                          interactionLDA;

  if ( interaction._lowRankDiagonalInteraction )
  {
    multiplyInteractionLowRankDiagonal( interaction, expansionWorkspace,
                                        numRows, columnRange, extraData,
                                        multWorkspace, outputMatrix );
  }
  else if ( interaction._compressed )
  {
    multiplyInteractionCompressed( interaction, expansionWorkspace,
                                   numRows, columnRange, extraData,
                                   multWorkspace, outputMatrix );
  }
  else
  {
    multiplyInteractionUncompressed( interaction, expansionWorkspace,
                                     numRows, columnRange, extraData,
                                     multWorkspace, outputMatrix );
  }
}

//////////////////////////////////////////////////////////////////////
// Version of the above function where the interaction arises from
// a low rank diagonal block
//////////////////////////////////////////////////////////////////////
void Supernode::multiplyInteractionLowRankDiagonal(
                          const SupernodeInteraction &interaction,
                          const Real *expansionWorkspace,
                          int numRows, const IndexRange &columnRange,
                          Real *extraData, Real *multWorkspace,
                          Real *outputMatrix )
{
  const Real                  *interactionData;
  const Real                  *multiplyData;
  int                          numMultiplyColumns;
  int                          numOriginalColumns;
  int                          nCols = numColumns();

  int                          interactionLDA;

  IndexRange                   interactionRange;
  int                          multiplyLDA;

  interactionRange = range_intersection( interaction._compressedColumnRange,
                                         columnRange );

  numMultiplyColumns = range_size( interactionRange );

  if ( numMultiplyColumns < 1 )
  {
    return;
  }

  TRACE_ASSERT(
    interactionRange.first >= interaction._compressedColumnRange.first );
  TRACE_ASSERT(
    interactionRange.second <= interaction._compressedColumnRange.second );
  TRACE_ASSERT( interactionRange.first >= columnRange.first );
  TRACE_ASSERT( interactionRange.second <= columnRange.second );

  interactionData = extraData + interaction._extendedDataOffset;

  // Offset by the start of the column range
  interactionData += ( interactionRange.first
                       - interaction._compressedColumnRange.first );

  multiplyData = expansionWorkspace;

  // Offset by the start of the column range
  multiplyData += ( interactionRange.first - columnRange.first );

  interactionLDA = range_size( interaction._compressedColumnRange );
  multiplyLDA = range_size( columnRange );

  // Multiply
  MATRIX::gemm( interactionData, multiplyData, outputMatrix,
                interaction._numExtendedRows, numMultiplyColumns,
                numRows, numMultiplyColumns,
                // Transpose just the second matrix
                false, true,
                // Subtract
                -1.0, 1.0,
                // Leading dimension of interactionData
                interactionLDA,
                // Leading dimension of multiplyData
                multiplyLDA );
}

//////////////////////////////////////////////////////////////////////
// Version of the above function where the interaction is stored
// in compressed column form
//////////////////////////////////////////////////////////////////////
void Supernode::multiplyInteractionCompressed(
                          const SupernodeInteraction &interaction,
                          const Real *expansionWorkspace,
                          int numRows, const IndexRange &columnRange,
                          Real *extraData, Real *multWorkspace,
                          Real *outputMatrix )
{
  const Real                  *interactionData;
  const Real                  *multiplyData;
  int                          numMultiplyColumns;
  int                          numOriginalColumns;
  int                          nCols = numColumns();

  int                          interactionLDA;

  IndexRange                   interactionRange;

  // Figure out over which range the compressed column list for this
  // interaction intersects the given column range.
  findEntryRangeIntersection( interaction._compressedColumnList,
                              columnRange, interactionRange );

  if ( interactionRange.first < 0 || interactionRange.second < 0 )
  {
    // Nothing to do here
    return;
  }

  numMultiplyColumns = range_size( interactionRange );

  numOriginalColumns = range_size( columnRange );

  TRACE_ASSERT( numMultiplyColumns <= numOriginalColumns );

  // Extract the desired column range
  MATRIX::copyColumns( expansionWorkspace, multWorkspace,
                       // Columns to copy
                       interaction._compressedColumnList,
                       numRows,
                       // Number of columns in expansion workspace
                       numOriginalColumns,
                       // Number of columns to be copied
                       numMultiplyColumns,
                       // Index range to copy over
                       interactionRange,
                       // Offset everything by the start of columnRange
                       columnRange.first );

  multiplyData = multWorkspace;

  interactionData = extraData + interaction._extendedDataOffset;

  // Offset to the first relavant column
  interactionData += interactionRange.first;

  interactionLDA = interaction._compressedColumnList.size();

  MATRIX::gemm( interactionData, multiplyData, outputMatrix,
                interaction._numExtendedRows, numMultiplyColumns,
                numRows, numMultiplyColumns,
                // Transpose just the second matrix
                false, true,
                // Subtract
                -1.0, 1.0,
                // Leading dimension of interactionData.
                // All other leading dimensions are standard.
                interactionLDA );
}

//////////////////////////////////////////////////////////////////////
// Version of the above for a "regular" uncompressed interaction
//////////////////////////////////////////////////////////////////////
void Supernode::multiplyInteractionUncompressed(
                          const SupernodeInteraction &interaction,
                          const Real *expansionWorkspace,
                          int numRows, const IndexRange &columnRange,
                          Real *extraData, Real *multWorkspace,
                          Real *outputMatrix )
{
  const Real                  *interactionData;
  const Real                  *multiplyData;
  int                          numMultiplyColumns;
  int                          numOriginalColumns;
  int                          nCols = numColumns();

  int                          interactionLDA;

  // We are storing the full, uncompressed interaction, so multiply
  // according to columnRange
  interactionData = extraData + interaction._extendedDataOffset;

  numMultiplyColumns = range_size( columnRange );

  // Align with the first column in columnRange
  interactionData += columnRange.first;

  multiplyData = expansionWorkspace;

  interactionLDA = nCols;

  MATRIX::gemm( interactionData, multiplyData, outputMatrix,
                interaction._numExtendedRows, numMultiplyColumns,
                numRows, numMultiplyColumns,
                // Transpose just the second matrix
                false, true,
                // Subtract
                -1.0, 1.0,
                // Leading dimension of interactionData.
                // All other leading dimensions are standard.
                interactionLDA );
}

//////////////////////////////////////////////////////////////////////
// Builds a map from indices in this node's column range to the
// diagonal block handling those indices
//////////////////////////////////////////////////////////////////////
void Supernode::buildDiagonalMap()
{
  _diagonalMap.clear();
  _diagonalMap.resize( numColumns(), -1 );

  for ( int diag_block_idx = 0; diag_block_idx < _diagonalBlocks.size();
        diag_block_idx++ )
  {
    const DenseBlock        &block = _diagonalBlocks[ diag_block_idx ];

    for ( int i = block._rowRange.first; i <= block._rowRange.second; i++ )
    {
      TRACE_ASSERT( _diagonalMap[ i ] == -1 );

      _diagonalMap[ i ] = diag_block_idx;
    }
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Helper function for diagonalContributionMult
//
// Does multiplication with contributions resulting from the
// compression of off diagonal interactions
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalContributionMult_offDiagonalContributions(
                                           int block_idx,
                                           const Real *G, int nColsG,
                                           const Real *extraData,
                                           const Real *expansionWorkspace,
                                           Real *multWorkspaceInitial,
                                           Real *multWorkspaceFinal,
                                           bool transpose,
                                           bool left )
{
  int                        interaction_idx;
  int                        nCols = numColumns();
  const Real                *baseData;
  const Real                *rowData;
  const Real                *columnData;

  const DenseBlock          &block = _diagonalLowRankBlocks[ block_idx ];

  int                        nRowsBlock;
  int                        nColsBlock;

  nRowsBlock = transpose ? range_size( block._columnRange )
                         : range_size( block._rowRange );

  nColsBlock = transpose ? range_size( block._rowRange )
                         : range_size( block._columnRange );

  // Start with all contributions resulting from compression of
  // off-diagonal blocks
  //
  // FIXME: Currently this list does not include diagonal blocks,
  //        but if this changes then this code will also have to change
  for ( int contribution_idx = 0;
        contribution_idx < _diagonalContributions.size(); contribution_idx++ )
  {
    interaction_idx = _diagonalContributions[ contribution_idx ];

    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    if ( interaction._compressed )
    {
      baseData = expansionWorkspace;

      expansionWorkspace += nCols * interaction._numExtendedRows;
    }
    else
    {
      baseData = extraData + interaction._extendedDataOffset;
    }

    // Align with the correct row/column to form the particular block we need
    rowData = transpose ? baseData + block._columnRange.first
                        : baseData + block._rowRange.first;

    columnData = transpose ? baseData + block._rowRange.first
                           : baseData + block._columnRange.first;

    if ( left )
    {
      // Left Multiplication of the given matrix
      //

      // Let U be the current piece of interaction data we are working
      // with.  Start by forming Z = U( columnRange, : )' * G.
      // No need to transpose for this call, since U is stored in transposed
      // form.
      //
      // Z = U( rowRange, : )' * G if transpose == true
      MATRIX::gemm( columnData, G, multWorkspaceInitial,
                    interaction._numExtendedRows, /* Rows in U' */
                    nColsBlock, /* Columns in the subset of U' */
                    nColsBlock, nColsG, /* Size of G */
                    false, false, /* No transposition */
                    1.0, 0.0, /* Overwrite with the result */
                    // Explicitly set leading dimension because we
                    // are extracting a subset of U'
                    //
                    // Other leading dimensions can be left as is
                    nCols );

      // Form the product U( rowRange, : ) * Z, and add the result to
      // the output matrix
      //
      // Compute U( columnRange, : ) * Z if transpose == true
      MATRIX::gemm( rowData, multWorkspaceInitial, multWorkspaceFinal,
                    // Size of U subset (transpose is stored)
                    interaction._numExtendedRows, nRowsBlock,
                    interaction._numExtendedRows, nColsG, /* Size of Z */
                    true, false, /* Transpose U' to get U */
                    1.0, 1.0, /* Add to result */
                    // Explicitly set leading dimension, as above
                    nCols );
    }
    else
    {
      // Right multiplication of given matrix
      //

      // Form the product Z = G' * U( rowRange, : )
      //
      // Computes G' * U( columnRange, : ) if transpose == true
      MATRIX::gemm( G, rowData, multWorkspaceInitial,
                    nRowsBlock, nColsG, /* Size of G */
                    // Size of U subset (transpose is stored)
                    interaction._numExtendedRows, nRowsBlock,
                    true, true, /* Transpose both matrices */
                    1.0, 0.0, /* Overwrite with result */
                    nColsG, /* Leading dimension of G */
                    nCols /* Leading dimension of U */ );

      // Form the product Z * U( columnRange, : )' and add the result to
      // the output matrix
      //
      // Compute Z * U( rowRange, : )' if transpose == true
      MATRIX::gemm( multWorkspaceInitial, columnData, multWorkspaceFinal,
                    nColsG, interaction._numExtendedRows, /* Size of Z */
                    // Size of U subset (transpose is stored)
                    interaction._numExtendedRows, nColsBlock,
                    false, false, /* No transpose */
                    1.0, 1.0, /* Add to result */
                    interaction._numExtendedRows, /* Leading dimension of Z */
                    nCols /* Leading dimensino of U */ );
    }
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Helper function for diagonalContributionMult
//
// Does multiplication with contributions resulting from the
// compression of off diagonal interactions
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalContributionMult_offDiagonalContributions(
                                           int block_idx,
                                           const Real *G, int nColsG,
                                           const Real *extraData,
                                           Real *multWorkspaceInitial,
                                           Real *multWorkspaceFinal,
                                           Real *copyWorkspace,
                                           bool transpose,
                                           bool left )
{
  int                        interaction_idx;
  int                        nColsFull = numColumns();
  const Real                *baseData;
  const Real                *rowData;
  const Real                *columnData;
  const Real                *subMatrix;

  if ( _diagonalContributions.empty() ) {
    return;
  }

  const DenseBlock          &block = _diagonalLowRankBlocks[ block_idx ];

  int                        nRowsBlock;
  int                        nColsBlock;

  int                        numInteractionColumns;

  IndexRange                 interactionColumnRange;
  IndexRange                 interactionRowRange;

  const IndexRange          &rowRange = transpose ? block._columnRange
                                                  : block._rowRange;
  const IndexRange          &columnRange = transpose ? block._rowRange
                                                     : block._columnRange;

  // Start with all contributions resulting from compression of
  // off-diagonal blocks
  //
  // FIXME: Currently this list does not include diagonal blocks,
  //        but if this changes then this code will also have to change
  for ( int contribution_idx = 0;
        contribution_idx < _diagonalContributions.size(); contribution_idx++ )
  {
    interaction_idx = _diagonalContributions[ contribution_idx ];

    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    baseData = extraData + interaction._extendedDataOffset;

    if ( interaction._compressed )
    {
      // Figure out the range over which the compressed column list
      // interacts with the desired row and column ranges for the
      // block.
      findEntryRangeIntersection( interaction._compressedColumnList,
                                  columnRange, interactionColumnRange );
      findEntryRangeIntersection( interaction._compressedColumnList,
                                  rowRange, interactionRowRange );

      numInteractionColumns = interaction._compressedColumnList.size();
    }
    else
    {
      interactionColumnRange = columnRange;
      interactionRowRange = rowRange;

      numInteractionColumns = nColsFull;
    }

    if ( interactionColumnRange.first < 0 || interactionColumnRange.second < 0
      || interactionRowRange.first < 0 || interactionRowRange.second < 0 )
    {
      // Nothing to do for this interaction: it's row or column range doesn't
      // interact with the block we are trying to form.
      continue;
    }

    // Get the desired row range from G if we are doing left multiplication
    // and the column range if we are doing right multiplication.
    //
    // Note that we use G' when multiplying on the right, so we will
    // always technically extract rows from G
    if ( interaction._compressed )
    {
      MATRIX::copyRows( G, copyWorkspace,
                        // Rows to copy
                        interaction._compressedColumnList, nColsG,
                        // Range of rows from compressedColumnList
                        // to copy
                        left ? interactionColumnRange : interactionRowRange,
                        // Offset by the first column in this block
                        left ? columnRange.first : rowRange.first );

      subMatrix = copyWorkspace;
    }
    else
    {
      // We use the whole G matrix in this case
      subMatrix = G;
    }

    // Align with the first row we want
    rowData = baseData + interactionRowRange.first;
    columnData = baseData + interactionColumnRange.first;

    nRowsBlock = range_size( interactionRowRange );
    nColsBlock = range_size( interactionColumnRange );

    if ( left )
    {
      // Left Multiplication of the given matrix
      //

      // Let U be the current piece of interaction data we are working
      // with.  Start by forming Z = U( columnRange, : )' * G.
      // No need to transpose for this call, since U is stored in transposed
      // form.
      //
      // Z = U( rowRange, : )' * G if transpose == true
      MATRIX::gemm( columnData, subMatrix, multWorkspaceInitial,
                    interaction._numExtendedRows, /* Rows in U' */
                    nColsBlock, /* Columns in the subset of U' */
                    nColsBlock, nColsG, /* Size of G sub-matrix */
                    false, false, /* No transposition since U is
                                     stored transposed */
                    1.0, 0.0, /* Overwrite with the result */
                    // Explicitly set leading dimension because we
                    // are extracting a subset of U'
                    //
                    // Other leading dimensions can be left as is
                    numInteractionColumns );

      if ( interaction._compressed )
      {
        // Form the product U( rowRange, : ) * Z, and add the result to
        // the output matrix
        //
        // Compute U( columnRange, : ) * Z if transpose == true
        MATRIX::gemm( rowData, multWorkspaceInitial, copyWorkspace,
                      // Size of U subset (transpose is stored)
                      interaction._numExtendedRows, nRowsBlock,
                      interaction._numExtendedRows, nColsG, /* Size of Z */
                      true, false, /* Transpose U' to get U,
                                      since U is stored transposed */
                      1.0, 0.0, /* Overwrite with the result */
                      // Explicitly set leading dimension, as above
                      numInteractionColumns );

        // Finally, scatter and add the contribution to the desired
        // row sub-set of the output matrix
        MATRIX::scatterAddRows( copyWorkspace, multWorkspaceFinal,
                                // Rows to copy
                                interaction._compressedColumnList, nColsG,
                                // Range of rows to copy from this list
                                interactionRowRange,
                                // Offset by the first row in the block
                                // we are forming
                                rowRange.first,
                                // Add
                                1.0 );
      }
      else
      {
        // Form the product U( rowRange, : ) * Z, and add the result to
        // the output matrix
        //
        // Compute U( columnRange, : ) * Z if transpose == true
        //
        // We can do this directly in this case, since U is dense
        MATRIX::gemm( rowData, multWorkspaceInitial, multWorkspaceFinal,
                      // Size of U subset (transpose is stored)
                      interaction._numExtendedRows, nRowsBlock,
                      interaction._numExtendedRows, nColsG, /* Size of Z */
                      true, false, /* Transpose U' to get U,
                                      since U is stored transposed */
                      1.0, 1.0, /* Add to the result */
                      // Explicitly set leading dimension, as above
                      numInteractionColumns );
      }
    }
    else
    {
      // Right multiplication of given matrix
      //

      // Form the product Z = G' * U( rowRange, : )
      //
      // Computes G' * U( columnRange, : ) if transpose == true
      MATRIX::gemm( subMatrix, rowData, multWorkspaceInitial,
                    nRowsBlock, nColsG, /* Size of G */
                    // Size of U subset (transpose is stored)
                    interaction._numExtendedRows, nRowsBlock,
                    true, true, /* Transpose both matrices */
                    1.0, 0.0, /* Overwrite with result */
                    nColsG, /* Leading dimension of G */
                    numInteractionColumns /* Leading dimension of U */ );

      if ( interaction._compressed )
      {
        // Form the product Z * U( columnRange, : )' and add the result to
        // the output matrix
        //
        // Compute Z * U( rowRange, : )' if transpose == true
        MATRIX::gemm( multWorkspaceInitial, columnData, copyWorkspace,
                      nColsG, interaction._numExtendedRows, /* Size of Z */
                      // Size of U subset (transpose is stored)
                      interaction._numExtendedRows, nColsBlock,
                      false, false, /* No transpose */
                      1.0, 0.0, /* Overwrite result */
                      interaction._numExtendedRows, /* Leading dimension of Z */
                      numInteractionColumns /* Leading dimension of U */ );

        // Finally, scatter and add the contribution to the desired
        // row sub-set of the output matrix
        MATRIX::scatterAddColumns( copyWorkspace, multWorkspaceFinal,
                                   // Rows to copy
                                   interaction._compressedColumnList, nColsG,
                                   // Number of columns in each matrix
                                   range_size( interactionColumnRange ),
                                   range_size( columnRange ),
                                   // Range of rows to copy from this list
                                   interactionColumnRange,
                                   // Offset by the first column in the block
                                   // we are forming
                                   columnRange.first,
                                   // Add
                                   1.0 );
      }
      else
      {
        // Form the product Z * U( columnRange, : )' and add the result to
        // the output matrix
        //
        // Compute Z * U( rowRange, : )' if transpose == true
        //
        // Can just add directly in this case because U is dense
        MATRIX::gemm( multWorkspaceInitial, columnData, multWorkspaceFinal,
                      nColsG, interaction._numExtendedRows, /* Size of Z */
                      // Size of U subset (transpose is stored)
                      interaction._numExtendedRows, nColsBlock,
                      false, false, /* No transpose */
                      1.0, 1.0, /* Overwrite result */
                      interaction._numExtendedRows, /* Leading dimension of Z */
                      numInteractionColumns /* Leading dimensino of U */ );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for diagonalContributionMult
//
// Does multiplication with contributions resulting from the
// compression of off diagonal interactions
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalContributionMult_offDiagonalContributions(
                                  int block_idx,
                                  const Real *G, int nColsG,
                                  const Real *extraData,
                                  WorkspaceManager<Real> &workspaceManager,
                                  Real *multWorkspaceFinal,
                                  bool transpose,
                                  bool left )
{
  int                        interaction_idx;
  int                        nColsFull = numColumns();
  const Real                *baseData;
  const Real                *rowData;
  const Real                *columnData;
  const Real                *subMatrix;

  if ( _diagonalContributions.empty() ) {
    return;
  }

  const DenseBlock          &block = _diagonalLowRankBlocks[ block_idx ];

  int                        nRowsBlock;
  int                        nColsBlock;

  int                        numInteractionColumns;

  IndexRange                 interactionColumnRange;
  IndexRange                 interactionRowRange;

  const IndexRange          &rowRange = transpose ? block._columnRange
                                                  : block._rowRange;
  const IndexRange          &columnRange = transpose ? block._rowRange
                                                     : block._columnRange;

  // Workspace data
  PairArray                  dataSizes( 2 );

  Real                      *copyWorkspace;
  Real                      *multWorkspaceInitial;

  // Start with all contributions resulting from compression of
  // off-diagonal blocks
  //
  // FIXME: Currently this list does not include diagonal blocks,
  //        but if this changes then this code will also have to change
  for ( int contribution_idx = 0;
        contribution_idx < _diagonalContributions.size(); contribution_idx++ )
  {
    interaction_idx = _diagonalContributions[ contribution_idx ];

    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    baseData = extraData + interaction._extendedDataOffset;

    if ( interaction._compressed )
    {
      // Figure out the range over which the compressed column list
      // interacts with the desired row and column ranges for the
      // block.
      findEntryRangeIntersection( interaction._compressedColumnList,
                                  columnRange, interactionColumnRange );
      findEntryRangeIntersection( interaction._compressedColumnList,
                                  rowRange, interactionRowRange );

      numInteractionColumns = interaction._compressedColumnList.size();
    }
    else
    {
      interactionColumnRange = columnRange;
      interactionRowRange = rowRange;

      numInteractionColumns = nColsFull;
    }

    if ( interactionColumnRange.first < 0 || interactionColumnRange.second < 0
      || interactionRowRange.first < 0 || interactionRowRange.second < 0 )
    {
      // Nothing to do for this interaction: it's row or column range doesn't
      // interact with the block we are trying to form.
      continue;
    }

    nRowsBlock = range_size( interactionRowRange );
    nColsBlock = range_size( interactionColumnRange );

    // Determine workspace sizes
    if ( interaction._compressed )
    {
      dataSizes[ 0 ] = IndexPair( max( nRowsBlock, nColsBlock ), nColsG );
    }
    else
    {
      dataSizes[ 0 ] = IndexPair( 0, 0 );
    }

    if ( left )
    {
      dataSizes[ 1 ] = IndexPair( interaction._numExtendedRows, nColsG );
    }
    else
    {
      dataSizes[ 1 ] = IndexPair( nColsG, interaction._numExtendedRows );
    }
    
    RealWorkspace            workspace( workspaceManager, dataSizes );

    copyWorkspace = workspace.workspaceData( 0 );
    multWorkspaceInitial = workspace.workspaceData( 1 );

    // Get the desired row range from G if we are doing left multiplication
    // and the column range if we are doing right multiplication.
    //
    // Note that we use G' when multiplying on the right, so we will
    // always technically extract rows from G
    if ( interaction._compressed )
    {
      MATRIX::copyRows( G, copyWorkspace,
                        // Rows to copy
                        interaction._compressedColumnList, nColsG,
                        // Range of rows from compressedColumnList
                        // to copy
                        left ? interactionColumnRange : interactionRowRange,
                        // Offset by the first column in this block
                        left ? columnRange.first : rowRange.first );

      subMatrix = copyWorkspace;
    }
    else
    {
      // We use the whole G matrix in this case
      subMatrix = G;
    }

    // Align with the first row we want
    rowData = baseData + interactionRowRange.first;
    columnData = baseData + interactionColumnRange.first;

    if ( left )
    {
      // Left Multiplication of the given matrix
      //

      // Let U be the current piece of interaction data we are working
      // with.  Start by forming Z = U( columnRange, : )' * G.
      // No need to transpose for this call, since U is stored in transposed
      // form.
      //
      // Z = U( rowRange, : )' * G if transpose == true
      MATRIX::gemm( columnData, subMatrix, multWorkspaceInitial,
                    interaction._numExtendedRows, /* Rows in U' */
                    nColsBlock, /* Columns in the subset of U' */
                    nColsBlock, nColsG, /* Size of G sub-matrix */
                    false, false, /* No transposition since U is
                                     stored transposed */
                    1.0, 0.0, /* Overwrite with the result */
                    // Explicitly set leading dimension because we
                    // are extracting a subset of U'
                    //
                    // Other leading dimensions can be left as is
                    numInteractionColumns );

      if ( interaction._compressed )
      {
        // Form the product U( rowRange, : ) * Z, and add the result to
        // the output matrix
        //
        // Compute U( columnRange, : ) * Z if transpose == true
        MATRIX::gemm( rowData, multWorkspaceInitial, copyWorkspace,
                      // Size of U subset (transpose is stored)
                      interaction._numExtendedRows, nRowsBlock,
                      interaction._numExtendedRows, nColsG, /* Size of Z */
                      true, false, /* Transpose U' to get U,
                                      since U is stored transposed */
                      1.0, 0.0, /* Overwrite with the result */
                      // Explicitly set leading dimension, as above
                      numInteractionColumns );

        // Finally, scatter and add the contribution to the desired
        // row sub-set of the output matrix
        MATRIX::scatterAddRows( copyWorkspace, multWorkspaceFinal,
                                // Rows to copy
                                interaction._compressedColumnList, nColsG,
                                // Range of rows to copy from this list
                                interactionRowRange,
                                // Offset by the first row in the block
                                // we are forming
                                rowRange.first,
                                // Add
                                1.0 );
      }
      else
      {
        // Form the product U( rowRange, : ) * Z, and add the result to
        // the output matrix
        //
        // Compute U( columnRange, : ) * Z if transpose == true
        //
        // We can do this directly in this case, since U is dense
        MATRIX::gemm( rowData, multWorkspaceInitial, multWorkspaceFinal,
                      // Size of U subset (transpose is stored)
                      interaction._numExtendedRows, nRowsBlock,
                      interaction._numExtendedRows, nColsG, /* Size of Z */
                      true, false, /* Transpose U' to get U,
                                      since U is stored transposed */
                      1.0, 1.0, /* Add to the result */
                      // Explicitly set leading dimension, as above
                      numInteractionColumns );
      }
    }
    else
    {
      // Right multiplication of given matrix
      //

      // Form the product Z = G' * U( rowRange, : )
      //
      // Computes G' * U( columnRange, : ) if transpose == true
      MATRIX::gemm( subMatrix, rowData, multWorkspaceInitial,
                    nRowsBlock, nColsG, /* Size of G */
                    // Size of U subset (transpose is stored)
                    interaction._numExtendedRows, nRowsBlock,
                    true, true, /* Transpose both matrices */
                    1.0, 0.0, /* Overwrite with result */
                    nColsG, /* Leading dimension of G */
                    numInteractionColumns /* Leading dimension of U */ );

      if ( interaction._compressed )
      {
        // Form the product Z * U( columnRange, : )' and add the result to
        // the output matrix
        //
        // Compute Z * U( rowRange, : )' if transpose == true
        MATRIX::gemm( multWorkspaceInitial, columnData, copyWorkspace,
                      nColsG, interaction._numExtendedRows, /* Size of Z */
                      // Size of U subset (transpose is stored)
                      interaction._numExtendedRows, nColsBlock,
                      false, false, /* No transpose */
                      1.0, 0.0, /* Overwrite result */
                      interaction._numExtendedRows, /* Leading dimension of Z */
                      numInteractionColumns /* Leading dimension of U */ );

        // Finally, scatter and add the contribution to the desired
        // row sub-set of the output matrix
        MATRIX::scatterAddColumns( copyWorkspace, multWorkspaceFinal,
                                   // Rows to copy
                                   interaction._compressedColumnList, nColsG,
                                   // Number of columns in each matrix
                                   range_size( interactionColumnRange ),
                                   range_size( columnRange ),
                                   // Range of rows to copy from this list
                                   interactionColumnRange,
                                   // Offset by the first column in the block
                                   // we are forming
                                   columnRange.first,
                                   // Add
                                   1.0 );
      }
      else
      {
        // Form the product Z * U( columnRange, : )' and add the result to
        // the output matrix
        //
        // Compute Z * U( rowRange, : )' if transpose == true
        //
        // Can just add directly in this case because U is dense
        MATRIX::gemm( multWorkspaceInitial, columnData, multWorkspaceFinal,
                      nColsG, interaction._numExtendedRows, /* Size of Z */
                      // Size of U subset (transpose is stored)
                      interaction._numExtendedRows, nColsBlock,
                      false, false, /* No transpose */
                      1.0, 1.0, /* Overwrite result */
                      interaction._numExtendedRows, /* Leading dimension of Z */
                      numInteractionColumns /* Leading dimensino of U */ );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for diagonalContributionMult
//
// Does multiplication with contributions resulting from the
// compression of blocks from this node's diagonal
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalContributionMult_diagonalContributions(
                                           int block_idx,
                                           const Real *G, int nColsG,
                                           const Real *extraData,
                                           Real *multWorkspaceInitial,
                                           Real *multWorkspaceFinal,
                                           bool transpose,
                                           bool left )
{
  int                        interaction_idx;
  int                        ancestor_block;
  const Real                *baseData;
  const Real                *rowData;
  const Real                *columnData;

  if ( _previousDiagonalContributions.empty() ) {
    return;
  }

  const DenseBlock          &block = _diagonalLowRankBlocks[ block_idx ];

  int                        nRowsBlock;
  int                        nColsBlock;

  nRowsBlock = transpose ? range_size( block._columnRange )
                         : range_size( block._rowRange );

  nColsBlock = transpose ? range_size( block._rowRange )
                         : range_size( block._columnRange );

  // Diagonal contributions from "parent" blocks in the diagonal that
  // have already been sparsified
  const vector<pair<int, DenseBlock> > &previousContributions
                                = _previousDiagonalContributions[ block_idx ];

  // Now handle previous contributions from the diagonal
  for ( int contribution_idx = 0;
        contribution_idx < previousContributions.size(); contribution_idx++ )
  {
    ancestor_block = previousContributions[ contribution_idx ].first;

    interaction_idx = _compressedDiagonalInteractions[ ancestor_block ];

    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    const DenseBlock &contributionRanges
                        = previousContributions[ contribution_idx ].second;

    baseData = extraData + interaction._extendedDataOffset;

    TRACE_ASSERT( contributionRanges._rowRange.first >= 0
               && ( contributionRanges._rowRange.second
                  < interaction._compressedColumnList.size() ),
                  "Row subset out of range" );

    TRACE_ASSERT( contributionRanges._columnRange.first >= 0
               && ( contributionRanges._columnRange.second
                  < interaction._compressedColumnList.size() ),
                  "Row subset out of range" );

    rowData = transpose ? baseData + contributionRanges._columnRange.first
                        : baseData + contributionRanges._rowRange.first;

    columnData = transpose ? baseData + contributionRanges._rowRange.first
                           : baseData + contributionRanges._columnRange.first;

    if ( left )
    {
      // Left multiplication of G
      //

      // Let U be the current piece of interaction data we are working
      // with.  Start by forming Z = U( columnRange, : )' * G.
      // No need to transpose for this call, since U is stored in transposed
      // form.
      MATRIX::gemm( columnData, G, multWorkspaceInitial,
                    interaction._numExtendedRows, /* Rows in U' */
                    nColsBlock, /* Columns in the subset of U' */
                    nColsBlock, nColsG, /* Size of G */
                    false, false, /* No transposition */
                    1.0, 0.0, /* Overwrite with the result */
                    // Explicitly set leading dimension because we
                    // are extracting a subset of U'
                    //
                    // Other leading dimensions can be left as is
                    interaction._compressedColumnList.size() );

      // Form the product U( rowRange, : ) * Z, and add the result to
      // the output matrix
      MATRIX::gemm( rowData, multWorkspaceInitial, multWorkspaceFinal,
                    // Size of U subset (transpose is stored)
                    interaction._numExtendedRows, nRowsBlock,
                    interaction._numExtendedRows, nColsG, /* Size of Z */
                    true, false, /* Transpose U' to get U */
                    1.0, 1.0, /* Add to result */
                    // Explicitly set leading dimension, as above
                    interaction._compressedColumnList.size() );
    }
    else
    {
      // Right multiplication of given matrix
      //

      // Form the product Z = G' * U( rowRange, : )
      //
      // Computes G' * U( columnRange, : ) if transpose == true
      MATRIX::gemm( G, rowData, multWorkspaceInitial,
                    nRowsBlock, nColsG, /* Size of G */
                    // Size of U subset (transpose is stored)
                    interaction._numExtendedRows, nRowsBlock,
                    true, true, /* Transpose both matrices */
                    1.0, 0.0, /* Overwrite with result */
                    nColsG, /* Leading dimension of G */
                    // Leading dimension of U
                    interaction._compressedColumnList.size() );

      // Form the product Z * U( columnRange, : )' and add the result to
      // the output matrix
      //
      // Compute Z * U( rowRange, : )' if transpose == true
      MATRIX::gemm( multWorkspaceInitial, columnData, multWorkspaceFinal,
                    nColsG, interaction._numExtendedRows, /* Size of Z */
                    // Size of U subset (transpose is stored)
                    interaction._numExtendedRows, nColsBlock,
                    false, false, /* No transpose */
                    1.0, 1.0, /* Add to result */
                    interaction._numExtendedRows, /* Leading dimension of Z */
                    // Leading dimension of U
                    interaction._compressedColumnList.size() );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for diagonalContributionMult
//
// Does multiplication with contributions resulting from the
// compression of blocks from this node's diagonal
//////////////////////////////////////////////////////////////////////
void Supernode::diagonalContributionMult_diagonalContributions(
                                  int block_idx,
                                  const Real *G, int nColsG,
                                  const Real *extraData,
                                  WorkspaceManager<Real> &workspaceManager,
                                  Real *multWorkspaceFinal,
                                  bool transpose,
                                  bool left )
{
  int                        interaction_idx;
  int                        ancestor_block;
  const Real                *baseData;
  const Real                *rowData;
  const Real                *columnData;

  const DenseBlock          &block = _diagonalLowRankBlocks[ block_idx ];

  int                        nRowsBlock;
  int                        nColsBlock;

  if ( _previousDiagonalContributions.empty() ) {
    return;
  }

  nRowsBlock = transpose ? range_size( block._columnRange )
                         : range_size( block._rowRange );

  nColsBlock = transpose ? range_size( block._rowRange )
                         : range_size( block._columnRange );

  // Diagonal contributions from "parent" blocks in the diagonal that
  // have already been sparsified
  const vector<pair<int, DenseBlock> > &previousContributions
                                = _previousDiagonalContributions[ block_idx ];

  // Workspace stuff
  PairArray                  dataSizes( 1 );

  Real                      *multWorkspaceInitial;

  // Now handle previous contributions from the diagonal
  for ( int contribution_idx = 0;
        contribution_idx < previousContributions.size(); contribution_idx++ )
  {
    ancestor_block = previousContributions[ contribution_idx ].first;

    interaction_idx = _compressedDiagonalInteractions[ ancestor_block ];

    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    const DenseBlock &contributionRanges
                        = previousContributions[ contribution_idx ].second;

    baseData = extraData + interaction._extendedDataOffset;

    TRACE_ASSERT( contributionRanges._rowRange.first >= 0
               && ( contributionRanges._rowRange.second
                  < interaction._compressedColumnList.size() ),
                  "Row subset out of range" );

    TRACE_ASSERT( contributionRanges._columnRange.first >= 0
               && ( contributionRanges._columnRange.second
                  < interaction._compressedColumnList.size() ),
                  "Row subset out of range" );

    rowData = transpose ? baseData + contributionRanges._columnRange.first
                        : baseData + contributionRanges._rowRange.first;

    columnData = transpose ? baseData + contributionRanges._rowRange.first
                           : baseData + contributionRanges._columnRange.first;

    if ( left )
    {
      // Left multiplication of G
      //

      // Generate workspaces
      dataSizes[ 0 ] = IndexPair( interaction._numExtendedRows, nColsG );
      dataSizes[ 1 ] = IndexPair( nRowsBlock, nColsG );

      RealWorkspace          workspace( workspaceManager, dataSizes );

      multWorkspaceInitial = workspace.workspaceData( 0 );

      // Let U be the current piece of interaction data we are working
      // with.  Start by forming Z = U( columnRange, : )' * G.
      // No need to transpose for this call, since U is stored in transposed
      // form.
      MATRIX::gemm( columnData, G, multWorkspaceInitial,
                    interaction._numExtendedRows, /* Rows in U' */
                    nColsBlock, /* Columns in the subset of U' */
                    nColsBlock, nColsG, /* Size of G */
                    false, false, /* No transposition */
                    1.0, 0.0, /* Overwrite with the result */
                    // Explicitly set leading dimension because we
                    // are extracting a subset of U'
                    //
                    // Other leading dimensions can be left as is
                    interaction._compressedColumnList.size() );

      // Form the product U( rowRange, : ) * Z, and add the result to
      // the output matrix
      MATRIX::gemm( rowData, multWorkspaceInitial, multWorkspaceFinal,
                    // Size of U subset (transpose is stored)
                    interaction._numExtendedRows, nRowsBlock,
                    interaction._numExtendedRows, nColsG, /* Size of Z */
                    true, false, /* Transpose U' to get U */
                    1.0, 1.0, /* Add to result */
                    // Explicitly set leading dimension, as above
                    interaction._compressedColumnList.size() );
    }
    else
    {
      // Right multiplication of given matrix
      //

      // Generate workspaces
      dataSizes[ 0 ] = IndexPair( nColsG, interaction._numExtendedRows );
      dataSizes[ 1 ] = IndexPair( nColsG, nColsBlock );

      RealWorkspace          workspace( workspaceManager, dataSizes );

      multWorkspaceInitial = workspace.workspaceData( 0 );

      // Form the product Z = G' * U( rowRange, : )
      //
      // Computes G' * U( columnRange, : ) if transpose == true
      MATRIX::gemm( G, rowData, multWorkspaceInitial,
                    nRowsBlock, nColsG, /* Size of G */
                    // Size of U subset (transpose is stored)
                    interaction._numExtendedRows, nRowsBlock,
                    true, true, /* Transpose both matrices */
                    1.0, 0.0, /* Overwrite with result */
                    nColsG, /* Leading dimension of G */
                    // Leading dimension of U
                    interaction._compressedColumnList.size() );

      // Form the product Z * U( columnRange, : )' and add the result to
      // the output matrix
      //
      // Compute Z * U( rowRange, : )' if transpose == true
      MATRIX::gemm( multWorkspaceInitial, columnData, multWorkspaceFinal,
                    nColsG, interaction._numExtendedRows, /* Size of Z */
                    // Size of U subset (transpose is stored)
                    interaction._numExtendedRows, nColsBlock,
                    false, false, /* No transpose */
                    1.0, 1.0, /* Add to result */
                    interaction._numExtendedRows, /* Leading dimension of Z */
                    // Leading dimension of U
                    interaction._compressedColumnList.size() );
    }
  }
}

#if 0
//////////////////////////////////////////////////////////////////////
// Helper function for addLowRankDiagonalContributions
//
// Adds contributions resulting from compression of off-diagonal
// interactions
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContributions_offDiagonalContributions(
                                            const Real *extraData,
                                            const Real *expansionWorkspace )
{
  int                        interaction_idx;
  const Real                *baseData;
  int                        nCols = numColumns();

  for ( int diagonal_contrib_idx = 0;
        diagonal_contrib_idx < _diagonalContributions.size();
        diagonal_contrib_idx++ )
  {
    interaction_idx = _diagonalContributions[ diagonal_contrib_idx ];

    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    if ( interaction._compressed )
    {
      // Get the interaction from the expansion workspace.
      baseData = expansionWorkspace;

      expansionWorkspace += nCols * interaction._numExtendedRows;
    }
    else
    {
      baseData = extraData + interaction._extendedDataOffset;
    }

    addLowRankDiagonalContribution( baseData, interaction._numExtendedRows,
                                    true /* Stored in transposed form */ );
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// Helper function for addLowRankDiagonalContributions
//
// Adds contributions resulting from compression of off-diagonal
// interactinos
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContributions_offDiagonalContributions(
                                        const Real *extraData,
                                        Real *expansionWorkspace )
{
  int                        interaction_idx;
  const Real                *baseData;
  int                        nCols = numColumns();

  for ( int diagonal_contrib_idx = 0;
        diagonal_contrib_idx < _diagonalContributions.size();
        diagonal_contrib_idx++ )
  {
    interaction_idx = _diagonalContributions[ diagonal_contrib_idx ];

    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    baseData = extraData + interaction._extendedDataOffset;

    addLowRankDiagonalContribution(
      baseData, interaction._numExtendedRows, expansionWorkspace,
      interaction._compressed ? &interaction._compressedColumnList : NULL,
      true /* Stored in transposed form */ );
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for addLowRankDiagonalContributions
//
// Adds contributions resulting from compression of off-diagonal
// interactinos
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContributions_offDiagonalContributions(
                                      const Real *extraData,
                                      WorkspaceManager<Real> &workspaceManager )
{
  int                        interaction_idx;
  const Real                *baseData;
  int                        nCols = numColumns();

  for ( int diagonal_contrib_idx = 0;
        diagonal_contrib_idx < _diagonalContributions.size();
        diagonal_contrib_idx++ )
  {
    interaction_idx = _diagonalContributions[ diagonal_contrib_idx ];

    const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

    baseData = extraData + interaction._extendedDataOffset;

    addLowRankDiagonalContribution(
      baseData, interaction._numExtendedRows, workspaceManager,
      interaction._compressed ? &interaction._compressedColumnList : NULL,
      true /* Stored in transposed form */ );
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for addLowRankDiagonalContributions
//
// Adds contributions resulting from compression of blocks from
// this node's diagonal
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContributions_diagonalContributions(
                                            const Real *extraData )
{
  int                        interaction_idx;
  const Real                *baseData;
  Real                      *baseOutputData = _data;

  // Add contributions resulting from compression of blocks in this
  // node's diagonal
  for ( int diag_block_idx = 0; diag_block_idx < _diagonalBlocks.size();
        diag_block_idx++ )
  {
    vector<pair<int, IndexRange> > &mainDiagonalContributions
                                = _mainDiagonalContributions[ diag_block_idx ];

    const DenseBlock        &diagBlock = _diagonalBlocks[ diag_block_idx ];
    int                      numBlockRows = diagBlock.numRows();

    for ( int contrib_idx = 0; contrib_idx < mainDiagonalContributions.size();
          contrib_idx++ )
    {
      // Get the index of the interaction corresponding to this
      // particular block
      interaction_idx = mainDiagonalContributions[ contrib_idx ].first;
      interaction_idx = _compressedDiagonalInteractions[ interaction_idx ];

      const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

      baseData = extraData + interaction._extendedDataOffset;

      // Align with the correct row (matrix is stored transposed)
      baseData += mainDiagonalContributions[ contrib_idx ].second.first;

      MATRIX::syrk( baseData, baseOutputData,
                    numBlockRows, /* Rows in the output matrix */
                    interaction._numExtendedRows,
                    true, /* Transpose the matrix */
                    1.0, 1.0, /* Add to the diagonal */
                    // Leading dimension for interaction data
                    interaction._compressedColumnList.size() );
    }

    baseOutputData += numBlockRows * numBlockRows;
  }
}

//////////////////////////////////////////////////////////////////////
// Function which assigns low rank decomposition information to
// a new slack variable node introduced to compensate for the
// given interaction
//////////////////////////////////////////////////////////////////////
void Supernode::assignNewExtendedNode( int interaction_idx,
                                       vector<Supernode> &nodes,
                                       const Real *V, const Real *Utrans,
                                       int rank, Real *extraData,
                                       long int &offset,
                                       long int &remainingSize,
                                       Real *copyWorkspace )
{
  int                            nRows;
  int                            nCols = numColumns();
  int                            newNodeID;
  int                            blockSize;
  int                            node_idx;
  const SupernodeInteraction    &interaction = _offDiagonal[ interaction_idx ];
  const Real                    *Vcompressed;

  Real                           Vnorm, Unorm, normScale;

#if 0
  // FIXME: debugging
  char buf[ 1024 ];
  sprintf( buf, "super_numeric/node_%d_block_%d_V.matrix", _nodeID,
           interaction_idx );
  MATRIX::write( V, interaction._rowList.size(), rank, buf );

  sprintf( buf, "super_numeric/node_%d_block_%d_U.matrix", _nodeID,
           interaction_idx );
  MATRIX::write( Utrans, rank, nCols, buf );
#endif

  SupernodeInteraction &newInteraction
                          = _offDiagonal[ interaction._extendedInteraction ];

  const PairArray &forwardInteractions = newInteraction._forwardInteractions;

  newNodeID = newInteraction._nodeID;

  Supernode &newNode = nodes[ newNodeID ];
  Supernode &ancestorNode = nodes[ interaction._nodeID ];

  TRACE_ASSERT( newInteraction._type == EXTENDED_NODE, "Not an extended node" );
  TRACE_ASSERT( newNode.numColumns() == 0, "Node already allocated" );
  TRACE_ASSERT( newInteraction._numExtendedRows == 0,
                "Interaction already allocated" );
  TRACE_ASSERT( ancestorNode._type == STANDARD_NODE,
                "Ancestor not a standard node" );

  // Scatter the compressed version of V to copyWorkspace
  MATRIX::clear( copyWorkspace, ancestorNode.numColumns(), rank );
  MATRIX::scatterRows( V, copyWorkspace, interaction._rowList, rank );

  // Get the frobenius (vector) norms of both U and V, so that we
  // can balance their relative sizes
  Unorm = VECTOR::norm2( Utrans, rank * nCols );
  Vnorm = VECTOR::norm2( V, interaction._rowList.size() * rank );
  
  // Figure out a scaling to apply to the two matrices
  normScale = Unorm / Vnorm;
  normScale = sqrt( normScale );
  // FIXME
  //normScale = 1.0;

  // Compressed and scattered versions of V
  Vcompressed = V;
  V = copyWorkspace;

  // Assign extra data to the new interaction
  blockSize = rank * nCols;

  newInteraction._extendedDataOffset = offset;
  newInteraction._numExtendedRows = rank;

  // Copy U' in to this data area
  MATRIX::copy( extraData + offset, Utrans, rank, nCols );

  // FIXME
  MATRIX::scale( extraData + offset, rank, nCols, 1.0 / normScale );

  // Set the scaling for later when we actually add the contribution
  newInteraction._diagonalContributionScale = 1.0 / normScale;

  offset += blockSize;
  remainingSize -= blockSize;

  // Look through each forward interaction and assign data
  TRACE_ASSERT( forwardInteractions.size() > 0, "No forward interactions" );

  for ( int forward_idx = 0; forward_idx < forwardInteractions.size();
        forward_idx++ )
  {
    node_idx = forwardInteractions[ forward_idx ].first;
    interaction_idx = forwardInteractions[ forward_idx ].second;

    Supernode &nextNode = nodes[ node_idx ];
    SupernodeInteraction &nextInteraction
                              = nextNode._offDiagonal[ interaction_idx ];

    if ( nextInteraction._compressed )
    {
      nCols = nextInteraction._compressedColumnList.size();
    }
    else
    {
      nCols = nextNode.numColumns();
    }
    blockSize = rank * nCols;

#if 0
    nextInteraction._extendedDataOffset = offset;
#endif
    nextInteraction._numExtendedRows = rank;

    if ( node_idx == interaction._nodeID )
    {
      // If this nextInteraction is in the off-diagonal of the other
      // node associated with interaction then place V' in this data block
      if ( nextInteraction._compressed )
      {
        TRACE_ASSERT( nextInteraction._compressedColumnList.size() > 0,
                      "Empty compressed column list" );

        TRACE_ASSERT(
          nextInteraction._compressedColumnList.size()
            == interaction._rowList.size(), "Row list size mismatch" );

        // Copy the compressed version to extra storage
        MATRIX::transposeBLAS( extraData + offset, Vcompressed, nCols, rank );
      }
      else
      {
        // Copy the full, scattered version of the matrix
        MATRIX::transposeBLAS( extraData + offset, V, nCols, rank );
      }

      // Set the scaling for later when we actually add the contribution
      MATRIX::scale( extraData + offset, rank, nCols, -1.0 * normScale );

      nextInteraction._diagonalContributionScale = normScale;

      nextInteraction._extendedDataOffset = offset;
      offset += blockSize;
      remainingSize -= blockSize;
    }
    else
    {
      TRACE_ASSERT( !nextInteraction._compressed,
                    "This interaction cannot be stored in compressed form" );

#if 0
      // Otherwise, just zero out the space
      MATRIX::clear( extraData + offset, rank, nCols );
#endif
    }

#if 0
    offset += blockSize;
    remainingSize -= blockSize;
#endif
  }

  // Allocate diagonal data and a column range in the new node.
  blockSize = rank * rank;

  newNode._numColumns = rank;
#if 0
  newNode._extendedDataOffset = offset;

  newNode.setColumnRange( nodes[ newNodeID - 1 ]._columnRange.second + 1,
                          nodes[ newNodeID - 1 ]._columnRange.second + rank );

  MATRIX::clear( extraData + offset, rank, rank );

  offset += blockSize;
  remainingSize -= blockSize;
#endif
}

//////////////////////////////////////////////////////////////////////
// Alternative to the function above, in which we simply assign the
// low-rank decomposition directly to a compressed interaction
// (ie. in place compression)
//////////////////////////////////////////////////////////////////////
void Supernode::assignCompressedInteraction( int interaction_idx,
                                             const Real *V, const Real *Utrans,
                                             int rank, Real *extraData,
                                             long int &offset,
                                             long int &remainingSize )
{
  SupernodeInteraction      &interaction = _offDiagonal[ interaction_idx ];
  int                        nCols = numColumns();

  interaction._numExtendedRows = rank;

  // Assign matrix data
  interaction._compressedV = extraData + offset;
  offset += interaction._rowList.size() * rank;
  interaction._compressedU = extraData + offset;
  offset += nCols * rank;

  remainingSize -= rank * ( interaction._rowList.size() + nCols );

  // Copy
  MATRIX::copy( interaction._compressedV, V,
                interaction._rowList.size(), rank );

  MATRIX::transposeBLAS( interaction._compressedU, Utrans, rank, nCols );

#if 0
  // FIXME: debugging
  char buf[ 1024 ];
  sprintf( buf, "super_numeric/node_%d_block_%d_V.matrix", _nodeID,
           interaction_idx );
  MATRIX::write( interaction._compressedV,
                 interaction._rowList.size(), rank, buf );

  sprintf( buf, "super_numeric/node_%d_block_%d_U.matrix", _nodeID,
           interaction_idx );
  MATRIX::write( interaction._compressedU, nCols, rank, buf );
#endif
}

//////////////////////////////////////////////////////////////////////
// Function which assigns low rank diagonal decomposition information
// to a new slack variable node introduced to compensate for the
// removal of a block in the diagonal
//////////////////////////////////////////////////////////////////////
void Supernode::assignNewExtendedNode_diagonal(
                                      int block_idx,
                                      std::vector<Supernode> &nodes,
                                      const Real *V, const Real *Utrans,
                                      int rank, Real *extraData,
                                      long int &offset,
                                      long int &remainingSize )
{
  int                        nRows;
  int                        nCols;
  int                        interaction_idx;
  int                        interactionSz;
  int                        newNodeID;
  int                        blockSize;
  int                        node_idx;

  Real                       Vnorm, Unorm, normScale;

  const DenseBlock          &block = _diagonalLowRankBlocks[ block_idx ];

  nRows = block.numRows();
  nCols = block.numColumns();

  Real                      *interactionData;

  interaction_idx = _compressedDiagonalInteractions[ block_idx ];

  // Get the interaction corresponding to this low rank block
  // and the node for the variables it introduces.
  SupernodeInteraction      &interaction = _offDiagonal[ interaction_idx ];

  Supernode                 &newNode = nodes[ interaction._nodeID ];

  const PairArray &forwardInteractions = interaction._forwardInteractions;

  newNodeID = newNode._nodeID;
  interactionSz = nRows + nCols;

  // FIXME: debugging
  if ( interaction._compressedColumnList.size() != interactionSz )
  {
    cout << SDUMP( interaction._compressedColumnList.size() ) << endl;
    cout << SDUMP( interactionSz ) << endl;
    cout << SDUMP( newNode.numExtendedColumns() ) << endl;
    cout << SDUMP( interaction._numExtendedRows ) << endl;
  }

  // Make sure everything has the right size
  TRACE_ASSERT( interaction._compressedColumnList.size() == interactionSz,
                "Compressed diagonal interaction has wrong size" );
  TRACE_ASSERT( interaction._compressed, "Interaction is not compressed" );

  // Assign interaction data
  blockSize = rank * interactionSz;

  interaction._extendedDataOffset = offset;
  interaction._numExtendedRows = rank;

  interactionData = extraData + offset;

  // Get the frobenius (vector) norms of both U and V, so that we
  // can balance their relative sizes
  Unorm = VECTOR::norm2( Utrans, rank * nCols );
  Vnorm = VECTOR::norm2( V, nRows * rank );
  
  // Figure out a scaling to apply to the two matrices
  normScale = Unorm / Vnorm;
  normScale = sqrt( normScale );

#if 0
  // FIXME
  if ( _diagonalLowRankBlocks.size() <= 3 )
  {
    // Only apply rescaling when system gets big enough?
    normScale = 1.0;
  }
#endif

  // Copy data, being careful to set the leading dimension
  MATRIX::copy( interactionData, Utrans, rank, nCols,
                interactionSz /* output leading dimension */,
                nCols /* input leading dimension */ );

  // FIXME
  MATRIX::scaleMatrix( interactionData, rank, nCols,
                       interactionSz /* lda */, 1.0 / normScale );

  interactionData += nCols;

  // Copy/transpose data, again being careful to set the output
  // leading dimension
  MATRIX::transposeBLAS( interactionData, V, nRows, rank, /* original size */
                         interactionSz /* output leading dimension */ );

  // FIXME: scale!!!!!!
  MATRIX::scaleMatrix( interactionData, rank, nRows,
                       interactionSz /* lda */, -1.0 * normScale );

  offset += blockSize;
  remainingSize -= blockSize;

  // Look through each forward interaction and assign data sizes
  for ( int forward_idx = 0; forward_idx < forwardInteractions.size();
        forward_idx++ )
  {
    node_idx = forwardInteractions[ forward_idx ].first;
    interaction_idx = forwardInteractions[ forward_idx ].second;

    Supernode               &nextNode = nodes[ node_idx ];
    SupernodeInteraction    &nextInteraction
                                  = nextNode._offDiagonal[ interaction_idx ];

    TRACE_ASSERT( !nextInteraction._compressed,
                  "This interaction should not be compressed" );

#if 0
    nCols = nextNode.numColumns();

    blockSize = rank * nCols;

    nextInteraction._extendedDataOffset = offset;
#endif
    nextInteraction._numExtendedRows = rank;

#if 0
    // Zero out the data
    MATRIX::clear( extraData + offset, rank, nCols );

    offset += blockSize;
    remainingSize -= blockSize;
#endif
  }

  // Allocate diagonal data and a column range for the new node
  blockSize = rank * rank;

  newNode._numColumns = rank;
#if 0
  newNode._extendedDataOffset = offset;

  newNode.setColumnRange( nodes[ newNodeID - 1 ]._columnRange.second + 1,
                          nodes[ newNodeID - 1 ]._columnRange.second + rank );

  MATRIX::clear( extraData + offset, rank, rank );

  offset += blockSize;
  remainingSize -= blockSize;
#endif
}

//////////////////////////////////////////////////////////////////////
// Alternative to the function above, in which we simply assign the
// low rank-decomposition of a block in the diagonal directly to a
// compressed interaction (ie. in place compression)
//////////////////////////////////////////////////////////////////////
void Supernode::assignCompressedInteraction_diagonal(
                                          int block_idx,
                                          const Real *V,
                                          const Real *Utrans,
                                          int rank, Real *extraData,
                                          long int &offset,
                                          long int &remainingSize )
{
  _nodeDiagonal.assignOffDiagonalBlock( block_idx, V, Utrans, rank,
                                        extraData, offset, remainingSize );
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but actually multiplies with the full compressed
// interaction matrix stored by this node, only for interactions following
// the provided ancestor_idx
//////////////////////////////////////////////////////////////////////
void Supernode::compressedInteractionMultiply_fullBlock(
                                                int ancestor_idx,
                                                const Real *G,
                                                int rowsG, int colsG,
                                                Real *output,
                                                bool left, bool transpose,
                                                bool extractSubMatrix )
{
  int                        nCols = numColumns();

  int                        nRows;

  // The part of _compressedV that we want starts here
  const Real                *V = _offDiagonal[ ancestor_idx + 1 ]._compressedV;
  const Real                *U = _compressedU;
  
  nRows = countInteractionRows( ancestor_idx + 1, _offDiagonal.size() - 1 );

  Real                      *outputWorkspace = output;

  MATRIX                     transposeWorkspace;
  MATRIX                     transposeOutputWorkspace;
  MATRIX                     copyWorkspace;
  MATRIX                     workspace;

  TIMING_START( "compressedInteractionMultiply_fullBlock: preamble" );

  TRACE_ASSERT( left, "Should never get here" );
  TRACE_ASSERT( !extractSubMatrix, "extractSubMatrix not implemented" );

  workspace.resizeAndWipe( _compressedRank, colsG );

  TIMING_STOP( "compressedInteractionMultiply_fullBlock: preamble" );

  if ( transpose ) {
    // Multiply by V'
    TIMING_START( "compressedInteractionMultiply_fullBlock: multiply V'" );
    MATRIX::gemm( V, G, workspace.data(),
                  // Dimensions of V
                  nRows, _compressedRank,
                  // Dimensions of G
                  nRows, colsG,
                  // Transpose V
                  true, false );
    TIMING_STOP( "compressedInteractionMultiply_fullBlock: multiply V'" );

    // Multiply by U, overwriting result
    TIMING_START( "compressedInteractionMultiply_fullBlock: multiply U" );
    MATRIX::gemm( U, workspace.data(), outputWorkspace,
                  // Dimensions of U
                  nCols, _compressedRank,
                  // Dimensions of workspace
                  _compressedRank, colsG,
                  // No transposition
                  false, false );
    TIMING_STOP( "compressedInteractionMultiply_fullBlock: multiply U" );
  }
  else {
    // Check to make sure G is big enough
    TRACE_ASSERT( nCols == rowsG );

    // Multiply by U'
    TIMING_START( "compressedInteractionMultiply_fullBlock: multiply U'" );
    MATRIX::gemm( U, G, workspace.data(),
                  // Dimensions of U
                  nCols, _compressedRank,
                  // Dimensions of G
                  rowsG, colsG,
                  // Transpose U
                  true, false );
    TIMING_STOP( "compressedInteractionMultiply_fullBlock: multiply U'" );

    // Multiply by V
    TIMING_START( "compressedInteractionMultiply_fullBlock: multiply V" );
    MATRIX::gemm( V, workspace.data(), outputWorkspace,
                  // Dimensions of V
                  nRows, _compressedRank,
                  // Dimensions of workspace
                  _compressedRank, colsG,
                  // No transposition
                  false, false );
    TIMING_STOP( "compressedInteractionMultiply_fullBlock: multiply V" );
  }
}

//////////////////////////////////////////////////////////////////////
// Multiplies the given matrix with an interaction stored in
// compressed form.  Left and right multiplication, as well as
// transposed multiplication are supported.
//////////////////////////////////////////////////////////////////////
void Supernode::compressedInteractionMultiply( int interaction_idx,
                                               const Real *G,
                                               int rowsG, int colsG,
                                               Real *output,
                                               bool left, bool transpose,
                                               bool extractSubMatrix,
                                               IndexRange rowRange,
                                               int rowOffset ) const
{
  int                          nCols = numColumns();

  const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

  int                          nRows = interaction._rowList.size();

  const Real                  *V = interaction._compressedV;
  const Real                  *U = interaction._compressedU;

  Real                        *outputWorkspace;

  MATRIX                       transposeWorkspace;
  MATRIX                       transposeOutputWorkspace;
  MATRIX                       copyWorkspace;
  MATRIX                       workspace;

  TIMING_START( "compressedInteractionMultiply: preamble" );

  if ( rowRange.first < 0 || rowRange.second < 0 ) {
    rowRange.first = 0;
    rowRange.second = nRows - 1;
  }
  nRows = range_size( rowRange );

  // Align V with the correct space in the row range
#if 0
  V += rowRange.first * nCols;
#endif
  V += rowRange.first * interaction._numExtendedRows;

  // If we are doing right side multiplication, we will just transpose
  // the matrix and treat this as though we were doing
  // right multiplication
  if ( !left ) {
    TRACE_ASSERT( NULL, "Should never get here" );
    transposeWorkspace.resizeAndWipe( colsG, rowsG );

    MATRIX::transposeBLAS( transposeWorkspace.data(), G, rowsG, colsG );

    G = transposeWorkspace.data();

    transposeOutputWorkspace.resizeAndWipe( transpose ? nRows : nCols,
                                            rowsG );

    outputWorkspace = transposeOutputWorkspace.data();

    transpose = !transpose;

    // Swap rows and columns
    int tmp = rowsG;
    rowsG = colsG;
    colsG = tmp;
  }
  else {
    outputWorkspace = output;
  }

  workspace.resizeAndWipe( interaction._numExtendedRows, colsG );

  TIMING_STOP( "compressedInteractionMultiply: preamble" );

  // We have now made everything look like a left multiplication
  if ( transpose ) {
    // Apply the V matrix first.  We need to start by extracting only the
    // rows of G referred to by the row set of this interaction (since these
    // are the only rows explicitly represented in V)
    if ( extractSubMatrix ) {
      copyWorkspace.resizeAndWipe( nRows, colsG );
      MATRIX::copyRows( G, copyWorkspace.data(), interaction._rowList, colsG,
                        rowRange, rowOffset );
      G = copyWorkspace.data();
    }

    // Multiply by V'
    TIMING_START( "compressedInteractionMultiply: multiply V'" );
    MATRIX::gemm( V, G, workspace.data(),
                  // Dimensions of V
                  nRows, interaction._numExtendedRows,
                  // Dimensions of G
                  nRows, colsG,
                  // Transpose V, but not copy workspace
                  true, false );
    TIMING_STOP( "compressedInteractionMultiply: multiply V'" );

    // Multiply by U, overwriting result
    TIMING_START( "compressedInteractionMultiply: multiply U" );
    MATRIX::gemm( U, workspace.data(), outputWorkspace,
                  // Dimensions of U
                  nCols, interaction._numExtendedRows,
                  // Dimensions of workspace
                  interaction._numExtendedRows, colsG,
                  // No transposition
                  false, false );
    TIMING_STOP( "compressedInteractionMultiply: multiply U" );
  }
  else {
    // U' represents all columns in the interaction, so no copying is
    // necessary here

    // Check to make sure G is big enough
    TRACE_ASSERT( nCols == rowsG );

    // Multiply by U'
    TIMING_START( "compressedInteractionMultiply: multiply U'" );
    MATRIX::gemm( U, G, workspace.data(),
                  // Dimensions of U
                  nCols, interaction._numExtendedRows,
                  // Dimensions of G
                  rowsG, colsG,
                  // Transpose U
                  true, false );
    TIMING_STOP( "compressedInteractionMultiply: multiply U'" );

    // Multiply by V
    TIMING_START( "compressedInteractionMultiply: multiply V" );
    MATRIX::gemm( V, workspace.data(), outputWorkspace,
                  // Dimensions of V
                  nRows, interaction._numExtendedRows,
                  // Dimensions of workspace
                  interaction._numExtendedRows, colsG,
                  // No transposition
                  false, false );
    TIMING_STOP( "compressedInteractionMultiply: multiply V" );
  }

  // If we just did a right multiplication, transpose back to the output
  if ( !left ) {
    TRACE_ASSERT( output != outputWorkspace );

    MATRIX::transposeBLAS( output, outputWorkspace,
                           transposeOutputWorkspace.rows(),
                           transposeOutputWorkspace.cols() );
  }

#if 0
  if ( nRows < interaction._rowList.size() ) {
    MATRIX::write( G, rowsG, colsG, "super_numeric/multInput.matrix" );
    MATRIX::write( output, transpose ? nCols : nRows, colsG,
                   "super_numeric/result.matrix" );
    abort();
  }
#endif
}

//////////////////////////////////////////////////////////////////////
// Inefficient compressedInteractionMultiply used for debugging
// purposes
//////////////////////////////////////////////////////////////////////
void Supernode::compressedInteractionMultiplyDebug(
                              int interaction_idx,
                              const Real *G,
                              int rowsG, int colsG,
                              Real *output,
                              bool left, bool transpose,
                              bool extractSubMatrix,
                              IndexRange rowRange,
                              int rowOffset ) const
{
  // Explicitly build the interaction matrix
  int                          nCols = numColumns();

  const SupernodeInteraction  &interaction = _offDiagonal[ interaction_idx ];

  int                          nRows = interaction._rowList.size();

  const Real                  *V = interaction._compressedV;
  const Real                  *U = interaction._compressedU;

  if ( rowRange.first < 0 || rowRange.second < 0 ) {
    rowRange.first = 0;
    rowRange.second = nRows - 1;
  }
  nRows = range_size( rowRange );

  V += rowRange.first * interaction._numExtendedRows;

  MATRIX                       interactionWorkspace( nRows, nCols );

  MATRIX::gemm( V, U, interactionWorkspace.data(),
                // V dimensions
                nRows, interaction._numExtendedRows,
                // U dimensions
                nCols, interaction._numExtendedRows,
                // Transpose U
                false, true );

#if 0
  // FIXME:
  // Ignore row offsets of the time being
  TRACE_ASSERT( rowRange.first == -1 && rowRange.second == -1,
                "rowRange functionality not implemented" );
#endif

  // Multiply
  if ( left && !transpose ) {
    MATRIX::gemm( interactionWorkspace.data(), G, output,
                  nRows, nCols, rowsG, colsG,
                  false, false );
  }
  else if ( left && transpose ) {
    MATRIX::gemm( interactionWorkspace.data(), G, output,
                  nRows, nCols, rowsG, colsG,
                  true, false );
  }
  else if ( !left && !transpose ) {
    MATRIX::gemm( G, interactionWorkspace.data(), output,
                  rowsG, colsG, nRows, nCols,
                  false, false );
  }
  else if ( !left && transpose ) {
    MATRIX::gemm( G, interactionWorkspace.data(), output,
                  rowsG, colsG, nRows, nCols,
                  false, true );
  }

#if 0
  if ( nRows < interaction._rowList.size() ) {
    MATRIX::write( G, rowsG, colsG, "super_numeric/multInput.matrix" );
    MATRIX::write( output, transpose ? nCols : nRows, colsG,
                   "super_numeric/result.matrix" );
    abort();
  }
#endif
}

#if 0
//////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContribution(
                      const SPARSE_MATRIX::SparseColumnMatrix &A,
                      const vector<Supernode> &nodes,
                      const IntArray &lowRankBlocks, int block_idx,
                      const vector<vector<IndexPair> > &lowRankDescendents,
                      const IntArray &ancestorInteractions,
                      Real scale, Real *multWorkspace, Real *workspace )
{
  Real                      *diagonalData = _data;
  int                        interaction_idx = lowRankBlocks[ block_idx ];
  int                        nCols = numColumns();

  const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

  int                        nRows = interaction._rowList.size();

  if ( !lowRankDiagonal() )
  {
    IndexRange               fullRange( 0, numColumns() - 1 );

    formInteractionBlockColumn( A, nodes,
                                // Actual interaction index in this node's
                                // full off-diagonal
                                lowRankBlocks[ block_idx ],
                                // List of descendents affecting this
                                // interaction
                                lowRankDescendents[ block_idx ],
                                ancestorInteractions,
                                fullRange, multWorkspace, workspace );

    // FIXME
    if ( _nodeID == 13456 )
    {
      TRACE_ASSERT( nRows * nCols <= workspaceSz, "Workspace too small" );

      MATRIX tmp( nRows, nCols, workspace );

      tmp.write( "super_numeric/testBlock.matrix" );
    }

    cout << "Adding contribution scaled by " << scale << endl;

    // Add to the diagonal
    MATRIX::syrk( workspace, diagonalData, nCols, nRows,
                  true, /* Transpose the matrix, eg. add S' * S */
                  scale, 1.0 /* Add to the existing diagonal block */ );

    // FIXME
    if ( _nodeID == 13456 )
    {
      MATRIX tmp( nCols, nCols );

      MATRIX::syrk( workspace, tmp.data(), nCols, nRows,
                    true,
                    1.0, 0.0 );

      tmp.write( "super_numeric/testProduct.matrix" );

      MATRIX::syrk( workspace, tmp.data(), nCols, nRows,
                    true,
                    scale, 0.0 );

      tmp.write( "super_numeric/testScaledProduct.matrix" );
    }
  }
  else
  {
    for ( int diag_block_idx = 0; diag_block_idx < _diagonalBlocks.size();
          diag_block_idx++ )
    {
      const IndexRange      &blockRange
                                = _diagonalBlocks[ diag_block_idx ]._rowRange;

      formInteractionBlockColumn( A, nodes,
                                  // Actual interaction index in this node's
                                  // full off-diagonal
                                  lowRankBlocks[ block_idx ],
                                  // List of descendents affecting this
                                  // interaction
                                  lowRankDescendents[ block_idx ],
                                  ancestorInteractions,
                                  blockRange, multWorkspace, workspace );
      
      nCols = range_size( blockRange );

      // Add to this diagonal block
      MATRIX::syrk( workspace, diagonalData, nCols, nRows,
                    true, /* Transpose the matrix, eg. add S' * S */
                    scale, 1.0 /* Add to the existing diagonal block */ );

      // Align data with the next diagonal block
      diagonalData += nCols * nCols;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Simpler version of the above for cases in which we just need to
// add a scaled identity matrix to the diagonal
//////////////////////////////////////////////////////////////////////
void Supernode::addLowRankDiagonalContribution( Real scale )
{
  if ( !lowRankDiagonal() )
  {
    cout << "Adding " << scale << " to the diagonal" << endl;
    MATRIX::addToDiagonal( _data, numColumns(), scale );
  }
  else
  {
    Real                    *diagonalData = _data;
    int                      nCols;

    for ( int diag_block_idx = 0; diag_block_idx < _diagonalBlocks.size();
          diag_block_idx++ )
    {
      nCols = range_size( _diagonalBlocks[ diag_block_idx ]._rowRange );

      MATRIX::addToDiagonal( diagonalData, nCols, scale );

      // Align with the next diagonal block
      diagonalData += nCols * nCols;
    }
  }
}
#endif

#if 0
//////////////////////////////////////////////////////////////////////
// Forms a block column of an interaction in this node, given
// a descendent list
//////////////////////////////////////////////////////////////////////
void Supernode::formInteractionBlockColumn(
                      const SPARSE_MATRIX::SparseColumnMatrix &A,
                      const vector<Supernode> &nodes,
                      int interaction_idx,
                      const vector<IndexPair> &descendents,
                      const IntArray &ancestorInteractions,
                      const IndexRange &columnRange,
                      Real *multWorkspace, Real *workspace )
{
  const SupernodeInteraction &interaction = _offDiagonal[ interaction_idx ];

  int                        ancestor_idx = interaction._nodeID;
  int                        descendent_idx;
  int                        descendent_interaction_idx;
  int                        descendent_ancestor_idx;
  
  const Supernode           &ancestor = nodes[ ancestor_idx ];

  // Get the base matrix data
  copyMatrixDataBlockColumn( interaction_idx, ancestor, A,
                             columnRange, workspace );

  for ( int i = 0; i < descendents.size(); i++ )
  {
    descendent_idx = descendents[ i ].first;
    descendent_interaction_idx = descendents[ i ].second;

    // Which interaction from descendent interactions with *this* node
    descendent_ancestor_idx = ancestorInteractions[ descendent_idx ];

    const Supernode         &descendent = nodes[ descendent_idx ];

    const SupernodeInteraction &descInteraction
                      = descendent._offDiagonal[ descendent_interaction_idx ];

    TRACE_ASSERT( descInteraction._nodeID == interaction._nodeID,
                  "Node ID mismatch" );

    TRACE_ASSERT(
      descendent._offDiagonal[ descendent_ancestor_idx ]._nodeID == _nodeID,
      "Ancestor node ID mismatch" );

    descendent.subtractInteractionUpdateMatrix( descendent_ancestor_idx,
                                                descendent_interaction_idx,
                                                columnRange,
                                                // Row list for this
                                                // block column
                                                interaction._rowList,
                                                multWorkspace, workspace );
  }
}
#endif

//////////////////////////////////////////////////////////////////////
// For each node in a list, identifies which entries in the
// original matrix contribute off-diagonal entries in the
// supernode and adds the corresponding row indices to a
// row set for that node.
//////////////////////////////////////////////////////////////////////
void Supernode::findSystemRows( const vector<Supernode> &nodes,
                                vector<set<int> > &rowSets,
                                const SPARSE_MATRIX::SparseColumnMatrix &A )
{
  TRACE_ASSERT( nodes.size() == rowSets.size(), "Row set size mismatch" );

  printf( "A has %d non-zero entries\n", (int)A._nzmax );

  for ( int i = 0; i < nodes.size(); i++ )
  {
    set<int>            &rowSet = rowSets[ i ];
    const Supernode     &node = nodes[ i ];

    if ( i % 1000 == 0 )
    {
      printf( "In Supernode::findSystemRows: node %06d of %06d\r",
              i + 1, (int)nodes.size() );
    }

    for ( int col_idx = node._columnRange.first;
          col_idx <= node._columnRange.second; col_idx++ )
    {
      for ( int row_ptr = A._p[ col_idx ]; row_ptr < A._p[ col_idx + 1 ];
            row_ptr++ )
      {
        if ( A._i[ row_ptr ] > node._columnRange.second )
        {
          rowSet.insert( A._i[ row_ptr ] );
        }
      }
    }
  }
  printf( "\n" );
}

//////////////////////////////////////////////////////////////////////
// Based on the current set of rows contained in the off-diagonal
// of each node in the nodes list, determine how additional rows
// are introduced by interacting supernodes.
//
// TODO: This probably isn't very fast, but I don't think it's too
// important right now
//////////////////////////////////////////////////////////////////////
void Supernode::propagateFillIn( const vector<Supernode> &nodes,
                                 vector<set<int> > &rowSets,
                                 const IntArray &nodeMap,
                                 bool compressInPlace )
{
  int                    nodeID;
  long int               totalAdded = 0;

  printf( "In propagateFillIn\n\n" );

  cout << SDUMP( rowSets.size() ) << endl;

  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  {
    set<int>            &rowSet = rowSets[ node_idx ];

    totalAdded += rowSet.size() * nodes[ node_idx ].numColumns();

    if ( node_idx % 1000 == 0 )
    {
      printf( "In Supernode::propagateFillIn: node %06d of %06d\r",
              node_idx + 1, (int)nodes.size() );
    }

    // Do not propagate fill in from compressed nodes
    //
    // BUT; if we are doing in place compression, then these interactions
    // will still affect fill in
    if ( nodes[ node_idx ]._compressOffDiagonal && !compressInPlace )
    {
      continue;
    }

    // Visit this node's off-diagonal rows in ascending order
    set<int>::iterator iter = rowSet.begin();
    while ( iter != rowSet.end() )
    {
      // Get the node ID for this row interaction
      nodeID = nodeMap[ *iter ];

      TRACE_ASSERT( nodeID > node_idx,
                    "nodeID should only refer to off-diagonal" );

      // Skip to the end of rows associated with this node,
      // since these will be included in the node's diagonal
      //for ( ; iter != rowSet.end(); iter++ )
      while ( iter != rowSet.end() )
      {
        if ( nodeMap[ *iter ] != nodeID )
        {
          break;
        }

        iter++;
      }

      // Iterate over the remaining rows in the set, and
      // add these as off-diagonal rows for nodeID
      for ( set<int>::iterator subIter = iter;
            subIter != rowSet.end(); subIter++ )
      {
        rowSets[ nodeID ].insert( *subIter );
      }

      // Iterate over the remaining rows in the set, and
      // add these as off-diagonal rows for nodeID
      //for ( set<int>::iterator subIter = iter;
    }
  }
  printf( "\n" );

  printf( "Added %ld non-zeros to the system\n", totalAdded );
}

//////////////////////////////////////////////////////////////////////
// Builds an map from node indices to super node indices
//////////////////////////////////////////////////////////////////////
void Supernode::buildSupernodeMap( const vector<Supernode> &nodes,
                                   IntArray &nodeMap )
{
  nodeMap.clear();

  nodeMap.resize( nodes.back()._columnRange.second + 1, -1 );

  for ( int i = 0; i < nodes.size(); i++ )
  {
    const Supernode     &node = nodes[ i ];

    if ( i % 1000 == 0 )
    {
      printf( "Building Supernode map for node %06d of %06d\r",
              i, (int)nodes.size() );
    }

    for ( int j = node._columnRange.first; j <= node._columnRange.second; j++ )
    {
      nodeMap.at( j ) = node._nodeID;
    }
  }
  printf( "\n" );
}

//////////////////////////////////////////////////////////////////////
// Based on a set of off-diagonal rows for each supernode, build
// all off diagonal interactions in a node list
//////////////////////////////////////////////////////////////////////
void Supernode::buildSupernodeInteractions(
                            vector<Supernode> &nodes,
                            const vector<set<int> > &rowSets,
                            IntArray &nodeMap )
{
  int            offDiagonalID;
  int            relativeRowIdx;

  // Visit each node and build its interaction list
  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  {
    if ( node_idx % 1000 == 0 )
    {
      printf( "Node %06d of %06d\r", node_idx + 1, (int)nodes.size() );
    }

    Supernode           &node = nodes[ node_idx ];
    const set<int>      &rowSet = rowSets[ node_idx ];

    node._numRows = rowSet.size();

    node._offDiagonal.clear();

    //for ( set<int>::iterator iter = rowSet.begin();
          //iter != rowSet.end(); iter++ )
    set<int>::iterator iter = rowSet.begin();
    while ( iter != rowSet.end() )
    {
      node._offDiagonal.push_back( SupernodeInteraction() );

      offDiagonalID = nodeMap[ *iter ];

      node._offDiagonal.back()._nodeID = offDiagonalID;
      node._offDiagonal.back()._type = STANDARD_NODE;

      // Go through the remaining rows from this node and
      // add them to the new interaction
      for ( ; iter != rowSet.end() && ( nodeMap[ *iter ] == offDiagonalID );
            iter++ )
      {
        // Get the row index relative to this supernode
        relativeRowIdx = *iter - nodes[ offDiagonalID ]._columnRange.first;

        node._offDiagonal.back()._rowList.push_back( relativeRowIdx );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Checks the given supernode to see whether or not parts of
// its off-diagonal need to be sparsified and introduces slack
// variable supernodes if so.
//
// TODO: Currently this only looks at off-diagonal
// blocks and doesn't try to eliminate off-diagonal blocks from
// the diagonal (if that makes sense...)
//////////////////////////////////////////////////////////////////////
void Supernode::addNode( Supernode &node,
                         const vector<Supernode> &nodes,
                         vector<Supernode> &newSystem,
                         int maxBlockSize )
{
  TRACE_ASSERT( NULL, "Needs to be updated" );

  // Temporary list for all slack variable nodes
  // introduced.
  vector<Supernode>              extendedNodes;
  IntArray                       extendedIndices;
  vector<SupernodeInteraction>   oldInteractions;
  int                            ncols, nrows;
  int                            numExtended = 0;
  Supernode                      nodeC = node;

  ncols = nodeC._columnRange.second - nodeC._columnRange.first + 1;

  // Look at each off-diagonal block in the nodeC, and flag
  // the ones we will mark as inactive
  for ( int block_idx = 0; block_idx < nodeC._offDiagonal.size(); block_idx++ )
  {
    SupernodeInteraction        &interaction = nodeC._offDiagonal[ block_idx ];

    nrows = interaction._rowList.size();

    if ( nrows * ncols > maxBlockSize )
    {
      extendedIndices.push_back( block_idx );

      // Zero out this interaction and replace it with a new set of
      // slack variables.
      interaction._active = false;

      numExtended += 1;

      interaction._extendedOffset = numExtended;
    }
  }

  // Introduce new supernodes
  extendedNodes.resize( numExtended, Supernode( EXTENDED_NODE ) );

  // Assign IDs to each extended node and have them build
  // interaction lists with each node appearing after it.
  for ( int node_idx = 0; node_idx < extendedNodes.size(); node_idx++ )
  {
    extendedNodes[ node_idx ]._nodeID = newSystem.size() + node_idx + 1;

    for ( int sub_node_idx = node_idx + 1; sub_node_idx < extendedNodes.size();
          sub_node_idx++ )
    {
      // Create an interaction with an empty row list.  This will
      // have to be determined at factorization time.
      extendedNodes[ node_idx ]._offDiagonal.push_back(
              SupernodeInteraction( newSystem.size() + sub_node_idx + 1,
                                    EXTENDED_NODE ) );
    }

    // This node has all valid interactions from the original source node.
    // It also has interactions with a subset of the nodes whose interactions
    // have been marked invalid.
    // 
    // Suppose the input node marks interactions indexed by B_1, B_2, B_3
    // as invalid.  These will introduce respective extended nodes
    // E_1, E_2, E_3.  Node E_1 will interact with the node originally
    // associated with B_1, E_2 with B_1 and B_2, and E_3 with B_1, B_2 and
    // B_3.
    for ( int block_idx = 0; block_idx < nodeC._offDiagonal.size(); block_idx++ )
    {
      SupernodeInteraction        &interaction = nodeC._offDiagonal[ block_idx ];

      if ( interaction._active )
      {
        // Copy the interaction directly
        extendedNodes[ node_idx ]._offDiagonal.push_back( interaction );
      }
      // Copy only inactive interactions occuring before the one
      // that created extendedNodes[ node_idx ]
      else if ( block_idx <= extendedIndices[ node_idx ] )
      {
        SupernodeInteraction       newInteraction = interaction;

        // This interaction will be dense, extending over all rows
        // in the affected supernode
        newInteraction._active = true;
        newInteraction._rowList.clear();

        // This interaction will be dense, so we need all
        // rows associated with it
        for ( int row_idx = 0;
              row_idx < nodes[ newInteraction._nodeID ].numColumns();
              row_idx++ )
        {
          newInteraction._rowList.push_back( row_idx );
        }

        extendedNodes[ node_idx ]._offDiagonal.push_back( newInteraction );
      }
    }
  }

  // Finally, we are ready to modify the input node itself
  oldInteractions = nodeC._offDiagonal;
  nodeC._offDiagonal.clear();

  nodeC._nodeID = newSystem.size();

  // Put interactions with the new extended nodes at the
  // beginning of the list.
  for ( int node_idx = 0; node_idx < extendedNodes.size(); node_idx++ )
  {
    // Leave the row list blank, since it won't be known until
    // we compute the numerical factor
    nodeC._offDiagonal.push_back(
      SupernodeInteraction( newSystem.size() + node_idx + 1, EXTENDED_NODE ) );
  }

  // Put the old interactions in after this
  for ( int block_idx = 0; block_idx < oldInteractions.size(); block_idx++ )
  {
    nodeC._offDiagonal.push_back( oldInteractions[ block_idx ] );
  }

  newSystem.push_back( nodeC );

  // Add the extended nodes to the new system
  for ( int node_idx = 0; node_idx < extendedNodes.size(); node_idx++ )
  {
    extendedNodes[ node_idx ]._nodeID = newSystem.size();

    newSystem.push_back( extendedNodes[ node_idx ] );
  }
}

//////////////////////////////////////////////////////////////////////
// Similar to the function above, this function checks the given
// supernode for off-diagonal blocks larger than the given size,
// and introduces slack variables in the place of these nodes.
// However, in this case, all slack variables are appended to the
// end of the system.
//////////////////////////////////////////////////////////////////////
void Supernode::addNodeAppend(
                   int node_idx,
                   vector<Supernode> &nodes,
                   int maxBlockSize,
                   vector<set<SupernodeInteraction> > &interactionSets,
                   Real compressionRatio,
                   RankEstimator *estimator,
                   bool compressInPlace )
{
  // Temporary list for all slack variable nodes
  // introduced.
  int                            ncols, nrows;
  int                            oldsize = nodes.size();
  int                            numAdded = 0;
  int                            node_interaction_idx;
  int                            nColsFull, nColsCompressed;

  //ncols = node._columnRange.second - node._columnRange.first + 1;
  ncols = nodes[ node_idx ].numColumns();

  if ( maxBlockSize > 0 && ncols > maxBlockSize )
  {
    Supernode::numLargeBlocks++;
  }

  //node._extensionStart = node._offDiagonal.size();
  nodes[ node_idx ]._extensionStart = nodes[ node_idx ]._offDiagonal.size();

  int offDiagSize = nodes[ node_idx ]._offDiagonal.size();

  // Look at each off-diagonal block in the node, and flag
  // the ones we will mark as inactive
  for ( int block_idx = 0; block_idx < nodes[ node_idx ]._offDiagonal.size();
        block_idx++ )
  {
    // Get the node each time, since the vector is changing and we
    // can't just hold on to the same reference
    Supernode                   &node = nodes[ node_idx ];
    SupernodeInteraction        &interaction = node._offDiagonal[ block_idx ];

    TRACE_ASSERT( nodes[ node_idx ]._offDiagonal.size() == offDiagSize,
                  "Uh-oh" );

    nrows = interaction._rowList.size();

    //if ( nrows * ncols > maxBlockSize )
    //if ( maxBlockSize > 0 && ncols > maxBlockSize )
    if ( node._compressOffDiagonal )
    {
      node._numRows -= interaction._rowList.size();

      // If we are doing in place compression here, then we just need
      // to flag this interaction as compressed.
      //
      // Otherwise, we need to introduce slack variables
      if ( compressInPlace )
      {
        interaction._type = COMPRESSED_NODE;
        interaction._compressedV = NULL;
        interaction._compressedU = NULL;
      }
      else
      {
        // Zero out this interaction and replace it with a new set of
        // slack variables.
        interaction._active = false;
        interaction._extendedOffset = nodes.size() - node._nodeID;

        // Add interactions with the two nodes affected by this
        interactionSets[ node._nodeID ].insert(
                  SupernodeInteraction( nodes.size(), EXTENDED_NODE,
                                        true /* First interaction */,
                                        block_idx, /* Parent interaction */
                                        false, /* not compressed */
                                        NULL, /* no compressed column list */
                                        NULL, /* no compressed column range */
                                        true /* diagonal contributor */ ) );

        const Supernode       &nextNode = nodes[ interaction._nodeID ];

        // Figure out how many additional columns will be introduced in
        // nextNode's new interaction if we store the interaction in
        // full vs. compressed form.
        nColsFull = nextNode.numColumns();
        nColsCompressed = nrows;

        if ( compressionRatio > 0.0
          && ( (Real)nColsFull / (Real)nColsCompressed ) > compressionRatio )
        {
          // Add a compressed interaction
          interactionSets[ interaction._nodeID ].insert(
                    SupernodeInteraction( nodes.size(), EXTENDED_NODE,
                                          false, /* firstInteraction */
                                          -1, /* extendedInteraction */
                                          true, /* compressed */
                                          &interaction._rowList,
                                          NULL, /* no compressed column range */
                                          true /* diagonal contributor */ ) );
        }
        else
        {
          interactionSets[ interaction._nodeID ].insert(
                    SupernodeInteraction( nodes.size(), EXTENDED_NODE,
                                          false, /* firstInteraction */
                                          -1, /* extendedInteraction */
                                          false, /* not compressed */
                                          NULL, /* no compressed column list */
                                          NULL, /* no compressed column range */
                                          true /* diagonal contributor */ ) );
        }

        nodes.push_back( Supernode( EXTENDED_NODE ) );
        nodes.back()._nodeID = nodes.size() - 1;

        if ( estimator )
        {
          // Estimate the size of this node
          nodes.back()._sizeEstimate = ( *estimator )( nrows, ncols );
        }

        numAdded++;
      }
    }
  }

  // Handle the node's main diagonal
  if ( !compressInPlace ) {
    compressDiagonalSymbolic( node_idx, nodes, interactionSets, estimator );
  }
}

//////////////////////////////////////////////////////////////////////
// Adds new nodes/interactions in order to compress off-diagonal blocks
// in the main diagonal of this node
//////////////////////////////////////////////////////////////////////
void Supernode::compressDiagonalSymbolic(
            int node_idx,
            vector<Supernode> &nodes,
            vector<set<SupernodeInteraction> > &interactionSets,
            RankEstimator *estimator )
{
  IntArray                   columnList;
  IndexRange                 columnRange;

  Supernode                 &node = nodes[ node_idx ];

  // We *don't* want a reference here, since we are modifying the
  // array and node could change.
  const vector<DenseBlock>   blocks = node._diagonalLowRankBlocks;

  // For each off-diagonal block in the main diagonal (confusing
  // terminology, I know), we add a new node, and a compressed
  // interaction with that node.
  for ( int block_idx = 0; block_idx < blocks.size(); block_idx++ )
  {
    // Build the full row and column range for this block
    const DenseBlock        &block = blocks[ block_idx ];

    columnList.clear();

    for ( int col_idx = block._columnRange.first;
          col_idx <= block._columnRange.second; col_idx++ )
    {
      columnList.push_back( col_idx );
    }

    for ( int row_idx = block._rowRange.first;
          row_idx <= block._rowRange.second; row_idx++ )
    {
      columnList.push_back( row_idx );
    }

    columnRange.first = block._columnRange.first;
    columnRange.second = block._rowRange.second;

#if 0
    printf( "Node %d: Introducing interaction to "
            "compensate for block [%d, %d] x [%d, %d]\n",
            node_idx, block._rowRange.first, block._rowRange.second,
            block._columnRange.first, block._columnRange.second );
#endif

    interactionSets[ node_idx ].insert(
              SupernodeInteraction( nodes.size(), EXTENDED_NODE,
                                    true, /* First interaction */
                                    /* No interaction as a parent */
                                    SupernodeInteraction::NO_INTERACTION_PARENT,
                                    true, /* compressed */
                                    &columnList,  // Provide a compressed
                                    &columnRange, // column list and range
                                    false, /* not a diagonal contributor */
                                    true, /* comes from a diagonal block */
                                    block_idx /* initial index in the node's
                                                 arrays */ ) );

    nodes.push_back( Supernode( EXTENDED_NODE ) );
    nodes.back()._nodeID = nodes.size() - 1;

    if ( estimator )
    {
      nodes.back()._sizeEstimate
        = ( *estimator )( block.numRows(), block.numColumns() );
    }
  }
}

static int totalInteractions = 0;

//////////////////////////////////////////////////////////////////////
// Having added extended nodes to the interaction set for a
// particular node, propagate the fill-in resulting from
// these nodes to all nodes later in the factorization.
//////////////////////////////////////////////////////////////////////
void Supernode::propagateExtendedFillIn(
                       const Supernode &node,
                       vector<Supernode> &nodes,
                       vector<set<SupernodeInteraction> > &interactionSets )
{
  set<SupernodeInteraction> &interactionSet = interactionSets[ node._nodeID ];

  set<SupernodeInteraction>::iterator interaction_iter;
  set<SupernodeInteraction>::iterator prev_iter;

  IndexRange                 baseColumnRange;
  IndexRange                 nextColumnRange;

  // Propagate all extended interactions forward through the factor
  for ( interaction_iter = interactionSet.begin();
        interaction_iter != interactionSet.end(); interaction_iter++ )
  {
    const SupernodeInteraction  &interaction = *( interaction_iter );

    TRACE_ASSERT( interaction._type == EXTENDED_NODE,
                  "Not propagating an extended node" );

    // Consider all interactions in this node up to this one
    for ( int prev_interaction_idx = 0;
          prev_interaction_idx < node._offDiagonal.size();
          prev_interaction_idx++ )
    {
      const SupernodeInteraction &prevInteraction
                      = node._offDiagonal[ prev_interaction_idx ];

      // If prevInteraction is active (standard or extended), we must
      // propogate the influence of interaction to the corresponding node
      if ( prevInteraction._active )
      {
        TRACE_ASSERT( NULL, "Should never get here: "
                            "functionality not implemented" );

        // Add an interaction to this node's set (the set stores unique
        // values, so repeated values aren't a problem).
        int oldSize = interactionSets[ prevInteraction._nodeID ].size();
        interactionSets[ prevInteraction._nodeID ].insert(
                SupernodeInteraction( interaction._nodeID, EXTENDED_NODE ) );

        if ( interactionSets[ prevInteraction._nodeID ].size() > oldSize )
        {
          printf( "Adding extended fill-in interaction to node %d\n",
                  prevInteraction._nodeID );
        }
      }
    }

    // Do the same thing for all new interactions
    for ( prev_iter = interactionSet.begin();
          prev_iter != interactionSet.end(); prev_iter++ )
    {
      // Stop once we reach the current interaction
      if ( prev_iter->_nodeID == interaction._nodeID )
      {
        break;
      }

      TRACE_ASSERT( nodes[ prev_iter->_nodeID ]._type == EXTENDED_NODE,
                    "This shouldn't happen" );
      TRACE_ASSERT( interaction._nodeID > prev_iter->_nodeID,
                    "Invalid interaction relationship" );

#if 0
      if ( interactionSets[ prev_iter->_nodeID ].find(
              SupernodeInteraction( interaction._nodeID, EXTENDED_NODE ) ) 
          != interactionSets[ prev_iter->_nodeID ].end() )
      {
        cout << "existing interaction found" << endl;
      }
#endif

      // FIXME
      if ( interaction._lowRankDiagonalInteraction )
      {
        TRACE_ASSERT( range_size( interaction._compressedColumnRange ) > 1 );
      }
      if ( prev_iter->_lowRankDiagonalInteraction )
      {
        TRACE_ASSERT( range_size( prev_iter->_compressedColumnRange ) > 1 );
      }

      // Skip this fill in operation if both interactions come from
      // compressing blocks on the diagonal AND their column ranges
      // do not overlap
      if ( interaction._lowRankDiagonalInteraction
        && prev_iter->_lowRankDiagonalInteraction
        && !range_overlap( interaction._compressedColumnRange,
                           prev_iter->_compressedColumnRange ) )
      {
        continue;
      }

      interactionSets[ prev_iter->_nodeID ].insert(
                SupernodeInteraction( interaction._nodeID, EXTENDED_NODE ) );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Adds all extended interactions from the given set to the system.
//
// Also sets up forward interactions amongst extended node interactions.
// This helps us allocate storage space during numerical factorization.
//////////////////////////////////////////////////////////////////////
void Supernode::addExtendedInteractions(
                vector<Supernode> &nodes,
                vector<set<SupernodeInteraction> > &interactionSets,
                int originalSize )
{
  // For each extended node, keep track of the indices of the
  // node/interaction that produced the extended node.
  PairArray                      firstInteractions( nodes.size() - originalSize,
                                                    IndexPair( -1, -1 ) );

  int                            extended_idx;

  set<SupernodeInteraction>::iterator interaction_iter;

  TRACE_ASSERT( nodes.size() == interactionSets.size(),
                "Mismatch between node and interaction sets" );

  for ( int node_idx = 0; node_idx < nodes.size(); node_idx++ )
  {
    Supernode                   &node = nodes[ node_idx ];
    set<SupernodeInteraction>   &interactions = interactionSets[ node_idx ];

    // We will need to build permutations of low rank diagonal block indices
    // so that for each node.  We use this to adjust internal arrays to respect
    // the possibly altered ordering of interactions resulting from low rank
    // blocks on the diagonal.
    IntArray diagonalPermutation( node._diagonalLowRankBlocks.size() );
    IntArray diagonalInversePermutation( node._diagonalLowRankBlocks.size() );
    int      currentDiagonalBlockIndex = 0;

    node._compressedDiagonalInteractions.clear();

    // These should be in sorted order of node ID
    for ( interaction_iter = interactions.begin();
          interaction_iter != interactions.end(); interaction_iter++ )
    {
      extended_idx = interaction_iter->_nodeID - originalSize;

      if ( interaction_iter->_firstInteraction )
      {
        if ( interaction_iter->_extendedInteraction
          != SupernodeInteraction::NO_INTERACTION_PARENT )
        {
          node._offDiagonal[
            interaction_iter->_extendedInteraction ]._extendedInteraction
              = node._offDiagonal.size();

          TRACE_ASSERT(
            node._offDiagonal[
              interaction_iter->_extendedInteraction ]._extendedOffset
                == ( interaction_iter->_nodeID - node._nodeID ),
            "Extended offset points to the wrong node" );
        }
        else if ( interaction_iter->_lowRankDiagonalInteraction )
        {
          TRACE_ASSERT( interaction_iter->_initialLowRankDiagonalIndex >= 0,
                        "No initial low rank diagonal block index" );

          int diag_block_idx = interaction_iter->_initialLowRankDiagonalIndex;

          diagonalPermutation.at( currentDiagonalBlockIndex ) = diag_block_idx;
          diagonalInversePermutation.at( diag_block_idx )
                                          = currentDiagonalBlockIndex;

          currentDiagonalBlockIndex++;

          // No interaction parent means that this is introduced by
          // the main diagonal, so store the interaction index in
          // the diagonal interaction list
          node._compressedDiagonalInteractions.push_back(
                                                  node._offDiagonal.size() );
        }

        firstInteractions[ extended_idx ].first = node_idx;
        firstInteractions[ extended_idx ].second = node._offDiagonal.size();
      }
      else
      {
        // We *must* have seen an interaction with this node already
        TRACE_ASSERT( firstInteractions[ extended_idx ].first >= 0,
                      "Extended node not yet visited" );

        // Find the first interaction for the node associated with this
        // interaction
        IndexPair &p = firstInteractions[ extended_idx ];

        PairArray &forwardInteractions
          = nodes[ p.first ]._offDiagonal[ p.second ]._forwardInteractions;

        // Put this new interaction in the list
        forwardInteractions.push_back( IndexPair( node_idx,
                                                  node._offDiagonal.size() ) );
      }

      if ( interaction_iter->_diagonalContributor )
      {
        node._diagonalContributions.push_back( node._offDiagonal.size() );
      }

      node._offDiagonal.push_back( *( interaction_iter ) );
    }

    TRACE_ASSERT(
            currentDiagonalBlockIndex == node._diagonalLowRankBlocks.size(),
            "We haven't built a full diagonal low rank block permuation" );

    // Now we want to un-permute the _compressedDiagonalInteractions index
    // list so that we can still compress the main diagonal as though
    // all interactions are in the "standard" order.
    IntArray                 newCompressedDiagonalInteractions;

    for ( int diag_block_idx = 0;
          diag_block_idx < node._compressedDiagonalInteractions.size();
          diag_block_idx++ )
    {
      newCompressedDiagonalInteractions.push_back(
                  node._compressedDiagonalInteractions.at(
                          diagonalInversePermutation.at( diag_block_idx ) ) );
    }
    node._compressedDiagonalInteractions = newCompressedDiagonalInteractions;
  }
}

//////////////////////////////////////////////////////////////////////
// Helper function for factorLDL which extracts diagonal data from
// the factors and places the relevant blocks and their inverses in
// to the provided datastructure
//////////////////////////////////////////////////////////////////////
void Supernode::ExtractLDLDiagonal( const Real *A, int rows,
                                    LDL_Data &factorData )
{
  int                        numSingleBlocks = 0;
  int                        numDoubleBlocks = 0;
  
  // Start by counting the number of single and double blocks
  // in the diagonal
  TRACE_ASSERT( rows == factorData._pivotData.size() );

  for ( int block_start = 0; block_start < rows; block_start++ ) {
    if ( factorData._pivotData[ block_start ] > 0 ) {
      numSingleBlocks += 1;
    }
    else if ( factorData._pivotData[ block_start ] < 0 ) {
      numDoubleBlocks += 1;

      // Advance to the next block
      block_start += 1;
    }
    else {
      TRACE_ASSERT( NULL, "Should never get here" );
    }
  }

  factorData._diagonalEntries.resize( numSingleBlocks );
  factorData._diagonalInverses.resize( numSingleBlocks );

  factorData._diagonalBlockEntries.resize( numDoubleBlocks );
  factorData._diagonalBlockInverses.resize( numDoubleBlocks );

  // Copy block data from the matrix and generate inverses
  numSingleBlocks = 0;
  numDoubleBlocks = 0;
  for ( int block_start = 0; block_start < rows; block_start++ )
  {
    if ( factorData._pivotData[ block_start ] > 0 )
    {
      factorData._diagonalEntries[ numSingleBlocks ]
        = MATRIX::access( A, rows, rows, block_start, block_start );

      factorData._diagonalInverses[ numSingleBlocks ]
        = 1.0 / factorData._diagonalEntries[ numSingleBlocks ];

      numSingleBlocks += 1;
    }
    else
    {
      Real *block = factorData._diagonalBlockEntries[ numDoubleBlocks ].data;

      for ( int row_idx = 0; row_idx < 2; row_idx++ )
      for ( int col_idx = 0; col_idx <= row_idx; col_idx++ )
      {
        MATRIX::access( block, 2, 2, row_idx, col_idx )
          = MATRIX::access( A, rows, rows,
                            block_start + row_idx, block_start + col_idx );
      }
      // Symmetrize
      MATRIX::access( block, 2, 2, 0, 1 ) = MATRIX::access( block, 2, 2, 1, 0 );

      // Form the block's inverse, via LU factorization
      MATRIX::copy(
                factorData._diagonalBlockInverses[ numDoubleBlocks ].data,
                block,
                2, 2 );
      MATRIX::LU(
                factorData._diagonalBlockInverses[ numDoubleBlocks ].data,
                factorData._diagonalBlockInverses[ numDoubleBlocks ].pivotData,
                2 );

      // Bump to the next block
      block_start += 1;

      numDoubleBlocks += 1;
    }
  }
}
