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
// NestedDissection.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "NestedDissection.h"

#include "Ordering.h"

#include <rschol/datastructure/SeparatorTree.h>

#include <rschol/geometry/FiniteDifference.h>

#include <rschol/util/IO.h>
#include <rschol/util/MERSENNETWISTER.h>
#include <rschol/util/STLUtil.h>

#ifdef USE_OMP
#include <omp.h>
#endif

#include <stack>

using namespace std;

//////////////////////////////////////////////////////////////////////
// NestedDissection::DissectionNode Implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
NestedDissection::DissectionNode::DissectionNode( int variableDimensions )
  : _left( NULL ),
    _right( NULL ),
    _variableDimensions( variableDimensions )
{
}

//////////////////////////////////////////////////////////////////////
// Desctructor
//////////////////////////////////////////////////////////////////////
NestedDissection::DissectionNode::~DissectionNode()
{
  delete _left;
  delete _right;
}

//////////////////////////////////////////////////////////////////////
// Builds a leaf node, which has no separators or children
// defined.
//////////////////////////////////////////////////////////////////////
NestedDissection::DissectionNode *
NestedDissection::DissectionNode::buildLeafNode( IntArray &indices,
                                                 int variableDimensions )
{
  DissectionNode *node = new DissectionNode( variableDimensions );

  node->_indices = indices;

  return node;
}

// Builds a permutation vector based on this nested dissection
// ordering.
void
NestedDissection::DissectionNode::buildPermutation( IntArray &permutation,
                                                    int variableDimensions )
{
  if ( _left == NULL && _right == NULL )
  {
    for ( int i = 0; i < _indices.size(); i++ )
    {
      for ( int d = 0; d < variableDimensions; d++ )
      {
        permutation.push_back( _indices[ i ] * variableDimensions + d );
      }
    }
  }
  else
  {
    if ( _left ) _left->buildPermutation( permutation, variableDimensions );
    if ( _right ) _right->buildPermutation( permutation, variableDimensions );

    for ( int i = 0; i < _separatorIndices.size(); i++ )
    {
      for ( int d = 0; d < variableDimensions; d++ )
      {
        permutation.push_back( _separatorIndices[ i ] * variableDimensions + d );
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Generates an array of all leaf nodes in this subtree
//////////////////////////////////////////////////////////////////////
void
NestedDissection::DissectionNode::getLeafNodes(
                          std::vector<const DissectionNode *> &leafNodes ) const
{
  if ( isLeaf() )
  {
    leafNodes.push_back( this );
  }
  else
  {
    if ( _left ) _left->getLeafNodes( leafNodes );
    if ( _right ) _right->getLeafNodes( leafNodes );
  }
}

//////////////////////////////////////////////////////////////////////
// Generates a post-ordering of this subtree
//////////////////////////////////////////////////////////////////////
void
NestedDissection::DissectionNode::postOrder(
                          std::vector<const DissectionNode *> &nodes ) const
{
  if ( !isLeaf() )
  {
    if ( _left ) _left->postOrder( nodes );
    if ( _right ) _right->postOrder( nodes );
  }

  nodes.push_back( this );
}

//////////////////////////////////////////////////////////////////////
// Builds an interior node, which has a separator and children
// defined.
//////////////////////////////////////////////////////////////////////
NestedDissection::DissectionNode *
NestedDissection::DissectionNode::buildInteriorNode( IntArray &separator,
                                                     DissectionNode *left,
                                                     DissectionNode *right,
                                                     int variableDimensions )
{
  DissectionNode *node = new DissectionNode( variableDimensions );

  node->_left = left;
  node->_right = right;

  TRACE_ASSERT( left != NULL || right != NULL,
                "Trying to build an interior node with NULL children" );

  node->_separatorIndices = separator;

  return node;
}

//////////////////////////////////////////////////////////////////////
// Determines the indices in a system reordered according to this
// tree structure that would be associated with this node.
//////////////////////////////////////////////////////////////////////
int NestedDissection::DissectionNode::buildReorderedIndices(
                                            int start, int variableDimensions )
{
  int totalIndices = 0;

  if ( _left == NULL && _right == NULL )
  {
    TRACE_ASSERT( _indices.size() > 0, "Empty leaf node" );

    // Base case: indices just correspond to whatever is
    // stored in this node.
    _reorderedRange.first = start;
    _reorderedRange.second = start + _indices.size() * variableDimensions - 1;

    totalIndices = _indices.size() * variableDimensions;
  }
  else
  {
    if ( _left )
    {
      totalIndices += _left->buildReorderedIndices( start,
                                                    variableDimensions );
    }
    if ( _right )
    {
      totalIndices += _right->buildReorderedIndices( start + totalIndices,
                                                     variableDimensions );
    }

    totalIndices += _separatorIndices.size() * variableDimensions;

    _reorderedRange.first = start;
    _reorderedRange.second = start + totalIndices - 1;
  }

  return totalIndices;
}

//////////////////////////////////////////////////////////////////////
// NestedDissection::DissectionTree Implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
NestedDissection::DissectionTree::DissectionTree()
  : _root( NULL )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
NestedDissection::DissectionTree::~DissectionTree()
{
}

//////////////////////////////////////////////////////////////////////
// NestedDissection Implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
NestedDissection::NestedDissection()
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
NestedDissection::~NestedDissection()
{
}

//////////////////////////////////////////////////////////////////////
// Builds a nested dissection ordering tree for a
// rectangular finite difference grid of the
// given size.
//////////////////////////////////////////////////////////////////////
NestedDissection::DissectionTree *
NestedDissection::nestedDissectionFD( Tuple3i gridDivisions,
                                      int maxBlockSize,
                                      int maxLevels,
                                      bool splitLargest,
                                      int variableDimensions,
                                      int maxSeparatorSize )
{
  Tuple3i gridStart( 0, 0, 0 );
  Tuple3i gridEnd( gridDivisions[0] - 1,
                   gridDivisions[1] - 1,
                   gridDivisions[2] - 1 );

  DissectionNode *root = nestedDissectionFD( gridDivisions, gridStart,
                                             gridEnd, maxBlockSize, maxLevels,
                                             splitLargest,
                                             variableDimensions,
                                             maxSeparatorSize );

  DissectionTree *tree = new DissectionTree();

  tree->root() = root;

  tree->buildReorderedIndices();

  return tree;
}

//////////////////////////////////////////////////////////////////////
// Converts a nested dissection node list in to a list of supernode
// specifications (directly assuming that separators correspond to
// supernodes)
//
// Assumes that reordered column ranges have been computed for the
// nodes
//
// maxBlockSize specifies the size above which off diagonals in the
// supernodes should be compressed.
//////////////////////////////////////////////////////////////////////
void NestedDissection::ConvertNodeList(
              const std::vector<const DissectionNode *> &nodeList,
              int maxBlockSize,
              std::vector<Ordering::SupernodeSpecification> &supernodes )
{
  supernodes.resize( nodeList.size() );

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    const DissectionNode    *node = nodeList[ node_idx ];

    IndexRange               columnRange = node->reorderedRange();

    if ( node->isLeaf() )
    {
      supernodes[ node_idx ]._compressOffDiagonal = false;
    }
    else
    {
      // Get the correct column range (ie. only separator columns)
      if ( node->right() )
      {
        columnRange.first = node->right()->reorderedRange().second + 1;
      }
      else
      {
        columnRange.first = node->left()->reorderedRange().second + 1;
      }


      if ( range_size( columnRange ) > maxBlockSize )
      {
        supernodes[ node_idx ]._compressOffDiagonal = true;
      }
    }

    supernodes[ node_idx ]._columnRange = columnRange;
  }
}

//////////////////////////////////////////////////////////////////////
// Performs nested dissection on an arbitraty matrix, given a way
// of converting point indices in to mesh points
//////////////////////////////////////////////////////////////////////
void NestedDissection::SupernodalNestedDissection(
                  const SPARSE_MATRIX::SparseColumnMatrix &A,
                  PointBuilder &pointPositionQuery,
                  IntArray &permutationVector,
                  std::vector<Ordering::SupernodeSpecification> &supernodes,
                  int maxBlockSize, int maxDiagonalBlockSize,
                  CompressionType compressionType,
                  std::vector<std::vector<DenseBlock> > &diagonalBlocks,
                  bool splitNodes,
                  bool reorderSeparators )
{
  // We will use CHOLMOD and CSparse for things like nested dissection,
  // eliminatation trees, etc.
  CHOLMOD_Environment::SeparatorList    separators;

  IntArray                   inversePermutationVector;
  IntArray                   separatorParents;
  IntArray                   separatorMap;
  IntArray                   separatorMapPermuted( A._nrow );
  IntArray                   separatorLevels;

  IntArray                   parent;
  IntArray                   postOrder( A._nrow );
  IntArray                   columnCounts;

  // Flat copies of a few things for CSparse
  int                       *permutation;
  int                       *inversePermutation;

  // Matrices in CSparse format (yes, I know this is messy)
  cs                        *Acopy = NULL;
  cs                        *Apermuted = NULL;

  // Get a nested dissection ordering from CHOLMOD
  {
    CHOLMOD_Environment      solver;

    solver.setMatrix( A );
    solver.nestedDissection( separators, permutationVector,
                             separatorParents, separatorLevels, separatorMap );

    // FIXME: Dump some nested dissection information to disk
    writeVector( separators, "separators.vector" );
    writeVector( separatorParents, "separatorParents.vector" );
  }

  // Permute the separator map
  //
  // Also build an identity post order for later on
  for ( int i = 0; i < permutationVector.size(); i++ )
  {
    separatorMapPermuted[ i ] = separatorMap[ permutationVector[ i ] ];
    postOrder[ i ] = i;
  }

  // Build the matrix in CSparse format
  CSparse_Interface::CopySparseSquareSymbolicMatrix( A, Acopy,
                                                     NULL, NULL,
                                                     CSparse_Interface::UPPER,
                                                     false );

  // Permutation inverse, etc.
  invertIntArray( permutationVector, inversePermutationVector );
  permutation = permutationVector.data();
  inversePermutation = inversePermutationVector.data();

  Apermuted = cs_symperm( Acopy, inversePermutation, 0 /* No values */ );

  // Build the elimination tree
  CSparse_Interface::ConstructEliminationTree( Apermuted, parent );

  // The separator tree is already post-ordered, so we can just use
  // an identity post order
  //
  // Get non-zero counts for each column
  CSparse_Interface::CholeskyColumnCounts( Apermuted, parent, postOrder,
                                           columnCounts );

  // Build our supernodes
  BuildRelaxedSupernodes( parent, columnCounts, Apermuted,
                          separators, separatorParents,
                          separatorLevels, maxBlockSize, supernodes );

  // Reorder diagonal blocks of large supernodes
  ReorderNestedSeparators( pointPositionQuery, supernodes,
                           maxDiagonalBlockSize, compressionType,
                           permutationVector,
                           diagonalBlocks,
                           1, /* variable dimension */
                           reorderSeparators );

  if ( splitNodes ) {
    SplitNodes( supernodes, diagonalBlocks );
  }

  // FIXME:
  int maxNodeSize = 0;
  for ( int i = 0; i < supernodes.size(); i++ ) {
    maxNodeSize = max( maxNodeSize, range_size( supernodes[i]._columnRange ) );
  }
  printf( "Max node size = %d\n", maxNodeSize );
#if 0
  char tmp;
  cout << "Enter any key..." << endl;
  cin >> tmp;
#endif

#if 0
  // FIXME: try splitting the nodes according to some maximum size
  SplitNodes( supernodes, diagonalBlocks, 12000 );
#endif

  cs_spfree( Acopy );
  cs_spfree( Apermuted );
}

//////////////////////////////////////////////////////////////////////
// Given a post ordering of a dissection tree from a finite difference
// grid computed
// with the function above, reorders separators in the tree according
// to a maximum block size.  This is done to expose low rank behaviour
// on the diagonals of blocks formed for these separators.
//
// This function adjusts the given permutation according to these
// reorderings.  It also provides a list of diagonal blocks for each
// node in the tree (though this list will be empty for many nodes).
//////////////////////////////////////////////////////////////////////
void NestedDissection::ReorderNestedSeparators(
                    Tuple3i gridDivisions, Real h,
                    const vector<const DissectionNode *> &nodeList,
                    int maxBlockSize,
                    IntArray &permutation,
                    vector<vector<DenseBlock> > &diagonalBlocks )
{
  const DissectionNode      *node;
  IndexRange                 nodeRange;

  vector<FiniteDifference::GridPoint *> separatorPermutation;

  vector<IndexRange>         diagonalRanges;

  typedef SeparatorTree<FiniteDifference::GridPoint> PointTree;

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    node = nodeList[ node_idx ];

    diagonalBlocks.push_back( vector<DenseBlock>() );

    if ( node->isLeaf() )
    {
      continue;
    }

    separatorPermutation.clear();

    nodeRange = node->reorderedRange();
    // We only want the indices contained in the separator
    if ( node->right() )
    {
      nodeRange.first = node->right()->reorderedRange().second + 1;
    }
    else
    {
      nodeRange.first = node->left()->reorderedRange().second + 1;
    }

    // Skip any blocks that we know are too small
    if ( range_size( nodeRange ) <= maxBlockSize )
    {
      continue;
    }

    for ( int point_idx = nodeRange.first; point_idx <= nodeRange.second;
          point_idx++ )
    {
      separatorPermutation.push_back(
        new FiniteDifference::GridPoint( permutation[ point_idx ],
                                         gridDivisions, h ) );
    }

    PointTree separatorTree( PointTree::LONGEST_AXIS,
                             // Max size of separator diagonal block
                             maxBlockSize,
                             PointTree::FIXED, 0.0 );

    separatorTree.build( separatorPermutation, false, true /* full tree */ );

    // FIXME: debugging
    if ( node_idx == nodeList.size() - 1 )
    {
      printf( "\n" );
      for ( int i = 0; i < separatorPermutation.size(); i++ )
      {
        printf( "separatorPermutation[ %d ] = %d", i,
                separatorPermutation[ i ]->index() );
        cout << "  " << separatorPermutation[ i ]->bbox().minBound() << endl;
      }
    }

    // Apply a sub-ordering to the permutation block storing indices
    // in this separator
    for ( int point_idx = nodeRange.first; point_idx <= nodeRange.second;
          point_idx++ )
    {
      permutation[ point_idx ]
            = separatorPermutation[ point_idx - nodeRange.first ]->index();
    }

    vector<DenseBlock>      &nodeDiagonalBlocks = diagonalBlocks.back();

    diagonalRanges.clear();

    // Get a list of diagonal blocks for the permutation we just
    // built
    separatorTree.leafRanges( diagonalRanges );

    for ( int diag_block_idx = 0; diag_block_idx < diagonalRanges.size();
          diag_block_idx++ )
    {
      nodeDiagonalBlocks.push_back(
                      DenseBlock( diagonalRanges[ diag_block_idx ],
                                  diagonalRanges[ diag_block_idx ] ) );
    }

    clearVectorContents( separatorPermutation );
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but uses a mesh with explicit vertex positions
//////////////////////////////////////////////////////////////////////
void NestedDissection::ReorderNestedSeparators(
                    const Mesh &mesh,
                    const vector<const DissectionNode *> &nodeList,
                    int maxBlockSize,
                    IntArray &permutation,
                    vector<vector<DenseBlock> > &diagonalBlocks,
                    int variableDimensions )
{
  const DissectionNode      *node;
  IndexRange                 nodeRange;

  vector<Mesh::MeshPoint *>  separatorPermutation;

  vector<IndexRange>         diagonalRanges;

  typedef SeparatorTree<Mesh::MeshPoint> PointTree;

  int                        separatorSize;
  int                        startPoint, endPoint;
  int                        permutedPoint;

  for ( int node_idx = 0; node_idx < nodeList.size(); node_idx++ )
  {
    node = nodeList[ node_idx ];

    diagonalBlocks.push_back( vector<DenseBlock>() );

    if ( node->isLeaf() )
    {
      continue;
    }

    separatorPermutation.clear();

    nodeRange = node->reorderedRange();
    // We only want the indices contained in the separator
    if ( node->right() )
    {
      nodeRange.first = node->right()->reorderedRange().second + 1;
    }
    else
    {
      nodeRange.first = node->left()->reorderedRange().second + 1;
    }

    separatorSize = range_size( nodeRange );

    TRACE_ASSERT( separatorSize % variableDimensions == 0,
                  "Invalid separator size" );

    separatorSize /= variableDimensions;

    // Skip any blocks that we know are too small
    if ( separatorSize <= maxBlockSize )
    {
      continue;
    }

    TRACE_ASSERT( nodeRange.first % variableDimensions == 0
               && ( nodeRange.second + 1 ) % variableDimensions == 0,
               "Invalid node range values" );

    startPoint = nodeRange.first / variableDimensions;
    endPoint = ( ( nodeRange.second + 1 ) / variableDimensions ) - 1;

    for ( int point_idx = startPoint; point_idx <= endPoint; point_idx++ )
    {
      permutedPoint = permutation[ point_idx * variableDimensions ];

      TRACE_ASSERT( permutedPoint % variableDimensions == 0,
                    "Invalid permuted point" );

      permutedPoint /= variableDimensions;

      separatorPermutation.push_back(
        new Mesh::MeshPoint( permutedPoint,
                             mesh.restPose()[ permutedPoint ] ) );
    }

    PointTree separatorTree( PointTree::LONGEST_AXIS,
                             // Max size of separator diagonal block
                             maxBlockSize,
                             PointTree::FIXED, 0.0 );

    separatorTree.build( separatorPermutation, false, true /* full tree */ );

    // Apply a sub-ordering to the permutation block storing indices
    // in this separator.
    //
    // This is a permutation of nodes, not matrix indices
    for ( int point_idx = startPoint; point_idx <= endPoint; point_idx++ )
    {
#if 0
      // TODO: This is wrong.  Fix it!
      permutedPoint = point_idx - startPoint;
      permutedPoint *= variableDimensions;
      permutedPoint = separatorPermutation[ permutedPoint ]->index();
#endif
      permutedPoint = separatorPermutation[ point_idx - startPoint ]->index();

      for ( int d = 0; d < variableDimensions; d++ )
      {
        permutation[ point_idx * variableDimensions + d ]
            = permutedPoint * variableDimensions + d;
      }

#if 0
      permutation[ point_idx ]
            = separatorPermutation[ point_idx - nodeRange.first ]->index();
#endif
    }

    vector<DenseBlock>      &nodeDiagonalBlocks = diagonalBlocks.back();

    diagonalRanges.clear();

    // Get a list of diagonal blocks for the permutation we just
    // built
    separatorTree.leafRanges( diagonalRanges );

    for ( int diag_block_idx = 0; diag_block_idx < diagonalRanges.size();
          diag_block_idx++ )
    {
      IndexRange             diagonalRange = diagonalRanges[ diag_block_idx ];

      // Modify to respect variable dimensions
      diagonalRange.first *= variableDimensions;

      diagonalRange.second += 1;
      diagonalRange.second *= variableDimensions;
      diagonalRange.second -= 1;

      nodeDiagonalBlocks.push_back(
                      DenseBlock( diagonalRange, diagonalRange ) );
    }

    clearVectorContents( separatorPermutation );
  }
}

//////////////////////////////////////////////////////////////////////
// Same as the above, but node positions will be specified implicitly
// by the given PointBuilder function
//////////////////////////////////////////////////////////////////////
void NestedDissection::ReorderNestedSeparators(
              const PointBuilder &pointPositionQuery,
              const std::vector<Ordering::SupernodeSpecification> &supernodes,
              int maxBlockSize,
              CompressionType compressionType,
              IntArray &permutation,
              vector<vector<DenseBlock> > &diagonalBlocks,
              int variableDimensions,
              bool reorderSeparators )
{
  vector<Mesh::MeshPoint *>  separatorPermutation;

  vector<IndexRange>         diagonalRanges;

  typedef SeparatorTree<Mesh::MeshPoint> PointTree;

  int                        separatorSize;
  int                        startPoint, endPoint;
  int                        permutedPoint;

  MERSENNETWISTER            generator;

  for ( int node_idx = 0; node_idx < supernodes.size(); node_idx++ )
  {
    const Ordering::SupernodeSpecification &node = supernodes[ node_idx ];

    diagonalBlocks.push_back( vector<DenseBlock>() );

    if ( compressionType != IN_PLACE && !node._compressOffDiagonal ) {
      continue;
    }

#if 0
    // FIXME: debugging
    if ( node_idx >= supernodes.size() - 1 ) {
      continue;
    }
#endif

    separatorPermutation.clear();

    const IndexRange &nodeRange = node._columnRange;

    separatorSize = range_size( nodeRange );

    TRACE_ASSERT( separatorSize % variableDimensions == 0,
                  "Invalid separator size" );

    separatorSize /= variableDimensions;

    // Skip any blocks that we know are too small
    if ( separatorSize <= maxBlockSize )
    {
      continue;
    }

    TRACE_ASSERT( nodeRange.first % variableDimensions == 0
               && ( nodeRange.second + 1 ) % variableDimensions == 0,
               "Invalid node range values" );

    startPoint = nodeRange.first / variableDimensions;
    endPoint = ( ( nodeRange.second + 1 ) / variableDimensions ) - 1;

    for ( int point_idx = startPoint; point_idx <= endPoint; point_idx++ )
    {
      permutedPoint = permutation[ point_idx * variableDimensions ];

      TRACE_ASSERT( permutedPoint % variableDimensions == 0,
                    "Invalid permuted point" );

      permutedPoint /= variableDimensions;

      separatorPermutation.push_back( pointPositionQuery( permutedPoint ) );
    }

    PointTree separatorTree( PointTree::LONGEST_AXIS,
                             // Max size of separator diagonal block
                             maxBlockSize,
                             PointTree::FIXED, 0.0 );

    separatorTree.build( separatorPermutation, false, true /* full tree */ );

    // Build the sub-ordering to apply here
    separatorPermutation.clear();
    separatorTree.orderNodes( separatorPermutation );

#if 0
    // FIXME: try randomly shuffling the separator indices
    int nSepPoints = separatorPermutation.size();
    for ( int i = 0; i < 10 * nSepPoints; i++) {
      int s1, s2;
      Mesh::MeshPoint *tmp;

      s1 = generator.randInt( nSepPoints - 1 );
      s2 = generator.randInt( nSepPoints - 1 );

      // Swap
      tmp = separatorPermutation[ s1 ];
      separatorPermutation[ s1 ] = separatorPermutation[ s2 ];
      separatorPermutation[ s2 ] = tmp;
    }
#endif

    if ( reorderSeparators ) {
      // Apply a sub-ordering to the permutation block storing indices
      // in this separator.
      //
      // This is a permutation of nodes, not matrix indices
      for ( int point_idx = startPoint; point_idx <= endPoint; point_idx++ )
      {
        permutedPoint = separatorPermutation[ point_idx - startPoint ]->index();

        for ( int d = 0; d < variableDimensions; d++ )
        {
          permutation[ point_idx * variableDimensions + d ]
              = permutedPoint * variableDimensions + d;
        }
      }
    }

    vector<DenseBlock>      &nodeDiagonalBlocks = diagonalBlocks.back();

    diagonalRanges.clear();

    // Get a list of diagonal blocks for the permutation we just
    // built
    separatorTree.leafRanges( diagonalRanges );

    for ( int diag_block_idx = 0; diag_block_idx < diagonalRanges.size();
          diag_block_idx++ )
    {
#if 0
      // FIXME: debugging
      if ( node_idx == supernodes.size() - 1 )
      {
        BoundingBox blockBBox;
        printf( "\n" );
        for ( int i = diagonalRanges[ diag_block_idx ].first;
              i <= diagonalRanges[ diag_block_idx ].second; i++ )
        {
#if 0
          printf( "Block %d; perm[ %d ] = %d", diag_block_idx, i,
                  separatorPermutation[ i ]->index() );
          cout << "  " << separatorPermutation[ i ]->bbox().minBound() << endl;
#endif

          blockBBox += separatorPermutation[ i ]->bbox().minBound();
        }
#if 0
        printf( "\n" );
        printf( "Block %d bounds: [ %f, %f ]\n", diag_block_idx,
                blockBBox.minBound()[ 0 ], blockBBox.maxBound()[ 0 ] );
        printf( "                 [ %f, %f ]\n",
                blockBBox.minBound()[ 1 ], blockBBox.maxBound()[ 1 ] );
        printf( "                 [ %f, %f ]\n",
                blockBBox.minBound()[ 2 ], blockBBox.maxBound()[ 2 ] );
        cout << SDUMP( blockBBox.longestAxis() ) << endl;
        printf( "\n" );
#endif
      }
#endif

      IndexRange             diagonalRange = diagonalRanges[ diag_block_idx ];

      // Modify to respect variable dimensions
      diagonalRange.first *= variableDimensions;

      diagonalRange.second += 1;
      diagonalRange.second *= variableDimensions;
      diagonalRange.second -= 1;

      nodeDiagonalBlocks.push_back(
                      DenseBlock( diagonalRange, diagonalRange ) );
    }

#if 0
    // FIXME: debugging
    nodeDiagonalBlocks.clear();
#endif

    clearVectorContents( separatorPermutation );
  }
}

//////////////////////////////////////////////////////////////////////
// Modifies supernode set so that each diagonal block is explicitly
// treated as a supernode.
//////////////////////////////////////////////////////////////////////
void NestedDissection::SplitNodes(
                      vector<Ordering::SupernodeSpecification> &supernodes,
                      vector<vector<DenseBlock> > &diagonalBlocks )
{
  vector<Ordering::SupernodeSpecification> newNodes;
  vector<vector<DenseBlock> >              newBlocks;

  for ( int node_idx = 0; node_idx < supernodes.size(); node_idx++ ) {
    const Ordering::SupernodeSpecification &nodeSpec = supernodes[ node_idx ];

    vector<DenseBlock> &blocks = diagonalBlocks[ node_idx ];

    if ( blocks.size() > 1 ) {
      // This node has diagonal blocks, so split it
      for ( int block_idx = 0; block_idx < blocks.size(); block_idx++ ) {
        const DenseBlock &block = blocks[ block_idx ];
        
        Ordering::SupernodeSpecification newNodeSpec;

        int startCol = nodeSpec._columnRange.first;

        newNodeSpec._columnRange.first = startCol + block._columnRange.first;
        newNodeSpec._columnRange.second = startCol + block._columnRange.second;
        newNodeSpec._compressOffDiagonal = true;

        newNodes.push_back( newNodeSpec );
        newBlocks.push_back( vector<DenseBlock>() );
      }
    } else {
      // This node is represented as a single block, so just copy it
      newNodes.push_back( nodeSpec );
      newBlocks.push_back( vector<DenseBlock>() );
    }
  }

  supernodes = newNodes;
  diagonalBlocks = newBlocks;
}

//////////////////////////////////////////////////////////////////////
// Splits nodes according to some maximum size
//////////////////////////////////////////////////////////////////////
void NestedDissection::SplitNodes(
                      vector<Ordering::SupernodeSpecification> &supernodes,
                      vector<vector<DenseBlock> > &diagonalBlocks,
                      int maxSize )
{
  vector<Ordering::SupernodeSpecification> newNodes;
  vector<vector<DenseBlock> >              newBlocks;

  for ( int node_idx = 0; node_idx < supernodes.size(); node_idx++ ) {
    const Ordering::SupernodeSpecification &nodeSpec = supernodes[ node_idx ];

    vector<DenseBlock> &blocks = diagonalBlocks[ node_idx ];

    vector<Ordering::SupernodeSpecification>   newSplitNodes;
    vector<vector<DenseBlock> >                newSplitBlocks;

    SplitNode( nodeSpec, blocks, newSplitNodes, newSplitBlocks, maxSize );

    TRACE_ASSERT( newSplitNodes.size() == newSplitBlocks.size() );
    for ( int i = 0; i < newSplitNodes.size(); i++ ) {
      newNodes.push_back( newSplitNodes[ i ] );
      newBlocks.push_back( newSplitBlocks[ i ] );
    }
  }

  supernodes = newNodes;
  diagonalBlocks = newBlocks;
}

//////////////////////////////////////////////////////////////////////
// Splits a single node according to some maximum size constraint
//////////////////////////////////////////////////////////////////////
void NestedDissection::SplitNode(
                      const Ordering::SupernodeSpecification &nodeSpec,
                      const vector<DenseBlock> &diagonalBlocks,
                      vector<Ordering::SupernodeSpecification> &newNodes,
                      vector<vector<DenseBlock> > &newBlocks,
                      int maxSize )
{
  stack<Ordering::SupernodeSpecification>  nodesToProcess;
  stack<IndexRange>                        nodeBlockRanges;

  if ( range_size( nodeSpec._columnRange ) <= maxSize ) {
    newNodes.push_back( nodeSpec );
    newBlocks.push_back( diagonalBlocks );

    return;
  }

  printf( "Splitting node with %d columns\n",
          range_size( nodeSpec._columnRange ) );

  nodesToProcess.push( nodeSpec );
  nodeBlockRanges.push( IndexRange( 0, diagonalBlocks.size() ) );

  while ( !nodesToProcess.empty() ) {
    Ordering::SupernodeSpecification       currentSpec = nodesToProcess.top();
    IndexRange                             currentRange = nodeBlockRanges.top();

    nodesToProcess.pop();
    nodeBlockRanges.pop();

    // Base case: if the node is small enough already, take its blocks
    // and move on
    if ( range_size( currentSpec._columnRange ) <= maxSize ) {
      newNodes.push_back( currentSpec );
      newBlocks.push_back( vector<DenseBlock>() );

      vector<DenseBlock>                  &currentBlocks = newBlocks.back();

      // We need the starting row for this block range, so that we can
      // set diagonal blocks appropriately for this new node
      int firstRow = diagonalBlocks[ currentRange.first ]._rowRange.first;

      for ( int block_idx = currentRange.first; block_idx < currentRange.second;
            block_idx++ )
      {
        // Transform this block in to the new node
        DenseBlock oldBlock = diagonalBlocks[ block_idx ];

        oldBlock._rowRange.first -= firstRow;
        oldBlock._rowRange.second -= firstRow;
        oldBlock._columnRange.first -= firstRow;
        oldBlock._columnRange.second -= firstRow;

        currentBlocks.push_back( oldBlock );
      }

      continue;
    }

    // Split this node in to two nodes, each with half as many
    // diagonal blocks
    IndexRange                             range1 = currentRange;
    IndexRange                             range2 = currentRange;

    range1.second = ( currentRange.first + currentRange.second ) / 2;
    range2.first = range1.second;

    Ordering::SupernodeSpecification       spec1 = currentSpec;
    Ordering::SupernodeSpecification       spec2 = currentSpec;

    // Figure out where these nodes start and end
    spec1._columnRange.second
      //=   spec1._columnRange.first
      =   nodeSpec._columnRange.first
        + diagonalBlocks[ range1.second - 1 ]._rowRange.second;
    spec2._columnRange.first = spec1._columnRange.second + 1;

    printf( "Splitting block range from [%d, %d] to [%d, %d]"
            " and [%d, %d]\n",
            currentRange.first, currentRange.second,
            range1.first, range1.second,
            range2.first, range2.second );

    printf( "Splitting column range from [%d, %d] to [%d, %d]"
            " and [%d, %d]\n",
            currentSpec._columnRange.first, currentSpec._columnRange.second,
            spec1._columnRange.first, spec1._columnRange.second,
            spec2._columnRange.first, spec2._columnRange.second );

    cout << "Pushing nodes" << endl;

    nodesToProcess.push( spec2 );
    nodesToProcess.push( spec1 );

    nodeBlockRanges.push( range2 );
    nodeBlockRanges.push( range1 );
  }
}

//////////////////////////////////////////////////////////////////////
// Recursive finite difference nested dissection function.
//////////////////////////////////////////////////////////////////////
NestedDissection::DissectionNode *
NestedDissection::nestedDissectionFD( Tuple3i gridDivisions,
                                      Tuple3i gridStart,
                                      Tuple3i gridEnd,
                                      int maxBlockSize,
                                      int maxLevels,
                                      bool splitLargest,
                                      int variableDimensions,
                                      int maxSeparatorSize )
{
  DissectionNode *node = NULL;

  Tuple3i gridSize( gridEnd[0] - gridStart[0] + 1,
                    gridEnd[1] - gridStart[1] + 1,
                    gridEnd[2] - gridStart[2] + 1 );

  int sz = gridSize[0] * gridSize[1] * gridSize[2];

  if ( sz == 0 )
  {
    // Nothing to do here
    return NULL;
  }

  if ( sz <= maxBlockSize || maxLevels <= 0 )
  {
    // Base case:  Take all indices in this range and put
    // them in a leaf node.
    IntArray indices;

    indices.reserve( sz );

    FOR_ALL_3D_POINTS_IN_RANGE( gridStart, gridEnd, idx )
    {
      int nodeIndex = INDEX_3D_GRID( gridDivisions, idx );

      indices.push_back( nodeIndex );
    }

    node = NestedDissection::DissectionNode::buildLeafNode(
                                                indices, variableDimensions );
  }
  else
  {
    int maxDimension, separatorIndex;

    if ( splitLargest )
    {
      // Get the largest dimension
      maxDimension = gridSize.maxIndex();
    }
    else
    {
      // Try splitting along the smallest dimension
      int sortedIndices[3];
      gridSize.sortIndices( sortedIndices, true );

      for ( int index_number = 0; index_number < 3; index_number++ )
      {
        maxDimension = sortedIndices[ index_number ];

        if ( gridSize[ maxDimension ] >= 3 )
          break;
      }
    }

    // Find a separator along this dimension
    separatorIndex = gridSize[ maxDimension ] / 2;

#if 0
    cout << SDUMP( maxLevels ) << endl;
    cout << SDUMP( maxDimension ) << endl;
    cout << SDUMP( separatorIndex ) << endl;
    cout << SDUMP( gridStart ) << endl;
    cout << SDUMP( gridEnd ) << endl;
    cout << SDUMP( gridSize ) << endl;
    cout << SDUMP( gridSize[ maxDimension ] ) << endl;
    cout << endl;
#endif

    int separatorSize = sz / gridSize[ maxDimension ];
    
    // If the separator size is smaller than the maximum separator
    // size, then we can treat this as a leaf node
    if ( separatorSize < maxSeparatorSize )
    {
      IntArray indices;

      indices.reserve( sz );

      FOR_ALL_3D_POINTS_IN_RANGE( gridStart, gridEnd, idx )
      {
        int nodeIndex = INDEX_3D_GRID( gridDivisions, idx );

        indices.push_back( nodeIndex );
      }

      node = NestedDissection::DissectionNode::buildLeafNode(
                                                  indices, variableDimensions );
    }
    else
    {
      // Figure out start and end indices for the two child grids

      // Left child - everything before the separator
      Tuple3i gridStart_left = gridStart;
      Tuple3i gridEnd_left = gridEnd;

      gridEnd_left[ maxDimension ]
        = gridStart[ maxDimension ] + separatorIndex - 1;

      // Right child - everything after the separator
      Tuple3i gridStart_right = gridStart;
      Tuple3i gridEnd_right = gridEnd;

      gridStart_right[ maxDimension ]
        = gridStart[ maxDimension ] + separatorIndex + 1;

      // Build our children
      DissectionNode *left = nestedDissectionFD( gridDivisions,
                                                 gridStart_left,
                                                 gridEnd_left,
                                                 maxBlockSize,
                                                 maxLevels - 1,
                                                 splitLargest,
                                                 variableDimensions,
                                                 maxSeparatorSize );
      DissectionNode *right = nestedDissectionFD( gridDivisions,
                                                  gridStart_right,
                                                  gridEnd_right,
                                                  maxBlockSize,
                                                  maxLevels - 1,
                                                  splitLargest,
                                                  variableDimensions,
                                                  maxSeparatorSize );

      // Get the separator indices
      Tuple3i separatorStart = gridStart;
      Tuple3i separatorEnd = gridEnd;

      separatorStart[ maxDimension ]
        = gridStart[ maxDimension ] + separatorIndex;
      separatorEnd[ maxDimension ]
        = gridStart[ maxDimension ] + separatorIndex;

      IntArray separatorIndices;

      separatorIndices.reserve( separatorSize );
      
      FOR_ALL_3D_POINTS_IN_RANGE( separatorStart, separatorEnd, idx )
      {
        int nodeIndex = INDEX_3D_GRID( gridDivisions, idx );

        separatorIndices.push_back( nodeIndex );
      }

      node = NestedDissection::DissectionNode::buildInteriorNode(
                                                          separatorIndices,
                                                          left, right,
                                                          variableDimensions );
    }
  }

  return node;
}

// Builds supernodes out of an ordering (with column counts).
// This code is heavily adapted from CHOLMOD's super_symbolic function,
// which performs supernode relaxation.
//
// Specify the list of separators from which the permutation is
// built.  Also indicate the maximum block size that we will leave
// uncompressed.  Any separators larger than this size will be treated
// as supernodes themselves.
void NestedDissection::BuildRelaxedSupernodes(
            const IntArray &Parent, const IntArray &columnCounts,
            const cs *Apermuted,
            const CHOLMOD_Environment::SeparatorList &separators,
            const IntArray &separatorParents,
            const IntArray &separatorLevels,
            int maxBlockSize,
            std::vector<Ordering::SupernodeSpecification> &supernodes,
            int nrelax0, int nrelax1, int nrelax2,
            Real zrelax0, Real zrelax1, Real zrelax2 )
{
  Real xxsize ;
  int *Ap, *Ai, *Ls, *Lpi, *Lpx, *Fnz,
      *Anz, *SuperMap, *Nscol, *Zeros, *Fp, *Fj,
      *ColCount, *Lpi2, *Lsuper, *Iwork ;
  int nsuper, d, n, j, k, s, mark, parent, p, pend, k1, k2, packed, nscol,
      nsrow, ndrow1, ndrow2, stype, ssize, xsize, sparent, plast, slast,
      csize, maxcsize, ss, nscol0, nscol1, ns, nfsuper, newzeros, totzeros,
      merge, snext, esize, maxesize, Asorted ;
  size_t w ;
  int ok = true ;

  int                        currentSeparator = 0;
  int                        maxCompressLevel = 0;

  /* ---------------------------------------------------------------------- */
  /* allocate workspace */
  /* ---------------------------------------------------------------------- */

  n = Apermuted->n ;

  /* ---------------------------------------------------------------------- */
  /* get inputs */
  /* ---------------------------------------------------------------------- */

  /* A is now either A or triu(A(p,p)) for the symmetric case.  It is either
   * A or A(p,f) for the unsymmetric case (both in column form).  It can be
   * either packed or unpacked, and either sorted or unsorted.  Entries in
   * the lower triangular part may be present if A is symmetric, but these
   * are ignored. */

  Ap = Apermuted->p ;
  Ai = Apermuted->i ;
  Anz = NULL;

  zrelax0 = std::isnan(zrelax0) ? 0 : zrelax0 ;
  zrelax1 = std::isnan(zrelax1) ? 0 : zrelax1 ;
  zrelax2 = std::isnan(zrelax2) ? 0 : zrelax2 ;

  /* ---------------------------------------------------------------------- */
  /* get workspace */
  /* ---------------------------------------------------------------------- */

  /* Sparent, Snz, and Merged could be allocated later, of size nfsuper */

  IntArray                   Wi( Parent.size() );
  IntArray                   Wj( Parent.size() );
  IntArray                   Sparent( Parent.size() );
  IntArray                   Snz( Parent.size() );
  IntArray                   Merged( Parent.size() );
  IntArray                   Super( Parent.size() );

  BoolArray                  IsSeparator( Parent.size(), false );

  IntArray                   Flag( Parent.size() );
  IntArray                   Head( Parent.size() + 1 );

  int                        EMPTY = -1;
  int                        TRUE = 1;

  /* ---------------------------------------------------------------------- */
  /* find the fundamental supernodes */
  /* ---------------------------------------------------------------------- */

  /* count the number of children of each node, using Wi [ */
  for (j = 0 ; j < n ; j++)
  {
    Wi [j] = 0 ;
  }
  for (j = 0 ; j < n ; j++)
  {
    parent = Parent [j] ;
    if (parent >= 0)
    {
      Wi [parent]++ ;
    }
  }

  // Figure out how many levels of the nested dissection tree to
  // compress
  for ( int sep_idx = 0; sep_idx < separators.size(); sep_idx++ )
  {
    const CHOLMOD_Environment::Separator &separator = separators[ sep_idx ];

    if ( separator._size > maxBlockSize && separator._numChildren > 0 )
    {
      maxCompressLevel = max( maxCompressLevel, separatorLevels[ sep_idx ] );
    }
  }

  // Flag array specifying whether to compress a given node
  BoolArray    compressNode( separators.size(), false );
  for ( int sep_idx = 0; sep_idx < separators.size(); sep_idx++ )
  {
    const CHOLMOD_Environment::Separator &separator = separators[ sep_idx ];

    if ( ( separator._size > maxBlockSize && separator._numChildren > 0 ) )
    {
      int                    parent_idx = separatorParents[ sep_idx ];

      compressNode[ sep_idx ] = true;

      while ( parent_idx >= 0 )
      {
        compressNode[ parent_idx ] = true;

        parent_idx = separatorParents[ parent_idx ];
      }
    }
  }

  nfsuper = 0;
  int col_idx = 0;
  bool newSegment = true;
  for ( int sep_idx = 0; sep_idx < separators.size(); sep_idx++ )
  {
    const CHOLMOD_Environment::Separator &separator = separators[ sep_idx ];

    //if ( separator._size > maxBlockSize && separator._numChildren > 0 )
    //if ( separatorLevels[ sep_idx ] <= maxCompressLevel )
    if ( compressNode[ sep_idx ] )
    {
      // Treat this whole separator as a single supernode
      Super[ nfsuper ] = separator._columnRange.first;
      IsSeparator[ nfsuper ] = true;
      nfsuper += 1;

      newSegment = true;

#if 0
      printf( "Fundamental supernode with %d columns\n", separator._size );
#endif
    }
    else
    {
      if ( newSegment )
      {
        Super[ nfsuper ] = separator._columnRange.first;
        nfsuper += 1;
        col_idx = separator._columnRange.first + 1;
        newSegment = false;
      }

      // Do standard supernodal formation in the sub separator pieces
      for ( ; col_idx <= separator._columnRange.second; col_idx += 1 )
      {
        if ( Parent[ col_idx - 1 ] != col_idx
          || columnCounts[ col_idx - 1 ] != columnCounts[ col_idx ] + 1
          || Wi[ col_idx ] > 1 )
        {
          // Start a new super node at this column
          Super[ nfsuper ] = col_idx;
          nfsuper += 1;
        }
      }
    }
  }
  Super [nfsuper] = n ;

  printf( "Found %d fundamental supernodes\n", nfsuper );

  /* contents of Wi no longer needed for child count ] */

  Nscol = Wi.data() ; /* use Wi as size-nfsuper workspace for Nscol [ */

  /* ---------------------------------------------------------------------- */
  /* find the mapping of fundamental nodes to supernodes */
  /* ---------------------------------------------------------------------- */

  SuperMap = Wj.data() ;	/* use Wj as workspace for SuperMap [ */

  /* SuperMap [k] = s if column k is contained in supernode s */
  for (s = 0 ; s < nfsuper ; s++)
  {
    for (k = Super [s] ; k < Super [s+1] ; k++)
    {
      SuperMap [k] = s ;
    }
  }

  /* ---------------------------------------------------------------------- */
  /* construct the fundamental supernodal etree */
  /* ---------------------------------------------------------------------- */

  for (s = 0 ; s < nfsuper ; s++)
  {
    j = Super [s+1] - 1 ;	/* last node in supernode s */
    parent = Parent [j] ;	/* parent of last node */
    Sparent [s] = (parent == EMPTY) ? EMPTY : SuperMap [parent] ;
  }

  /* contents of Wj no longer needed as workspace for SuperMap ]
   * SuperMap will be recomputed below, for the relaxed supernodes. */

  Zeros = Wj.data() ;   /* use Wj for Zeros, workspace of size nfsuper [ */

  /* ---------------------------------------------------------------------- */
  /* relaxed amalgamation */
  /* ---------------------------------------------------------------------- */

  for (s = 0 ; s < nfsuper ; s++)
  {
    Merged [s] = EMPTY ;			/* s not merged into another */
    Nscol [s] = Super [s+1] - Super [s] ;	/* # of columns in s */
    Zeros [s] = 0 ;				/* # of zero entries in s */
    Snz [s] = columnCounts [Super [s]] ;  /* # of entries in leading col of s */
  }

  for (s = nfsuper-2 ; s >= 0 ; s--)
  {
    ss = Sparent [s] ;
    if (ss == EMPTY)
    {
      continue ;
    }

    // Don't merge the supernodes that we are explicitly keeping around
    if ( IsSeparator[ s ] || IsSeparator[ s + 1 ] )
    {
      continue;
    }

    /* find the current parent of s (perform path compression as needed) */
    for (ss = Sparent [s] ; Merged [ss] != EMPTY ; ss = Merged [ss]) ;
    sparent = ss ;

    /* ss is the current parent of s */
    for (ss = Sparent [s] ; Merged [ss] != EMPTY ; ss = snext)
    {
      snext = Merged [ss] ;
      Merged [ss] = sparent ;
    }

    /* if s+1 is not the current parent of s, do not merge */
    if (sparent != s+1)
    {
      continue ;
    }

    nscol0 = Nscol [s] ;	/* # of columns in s */
    nscol1 = Nscol [s+1] ;	/* # of columns in s+1 */
    ns = nscol0 + nscol1 ;

    totzeros = Zeros [s+1] ;	/* current # of zeros in s+1 */

    /* determine if supernodes s and s+1 should merge */
    if (ns <= nrelax0)
    {
      merge = TRUE ;
    }
    else
    {
      /* use double to avoid integer overflow */
      double lnz0 = Snz [s] ;	/* # entries in leading column of s */
      double lnz1 = Snz [s+1] ;	/* # entries in leading column of s+1 */
      double xnewzeros = nscol0 * (lnz1 + nscol0 - lnz0) ;

      /* use Int for the final update of Zeros [s] below */
      newzeros = nscol0 * (Snz [s+1] + nscol0 - Snz [s]) ;

      if (xnewzeros == 0)
      {
        /* no new zeros, so go ahead and merge */
        merge = TRUE ;
      }
      else
      {
        /* # of zeros if merged */
        double xtotzeros = ((double) totzeros) + xnewzeros ;

        /* xtotsize: total size of merged supernode, if merged: */
        double xns = (double) ns ;
        double xtotsize  = (xns * (xns+1) / 2) + xns * (lnz1 - nscol1) ;
        double z = xtotzeros / xtotsize ;

        int totsize ;
        totsize  = (ns * (ns+1) / 2) + ns * (Snz [s+1] - nscol1) ;

        /* use Int for the final update of Zeros [s] below */
        totzeros += newzeros ;

        /* do not merge if supernode would become too big
         * (Int overflow).  Continue computing; not (yet) an error. */
        /* fl.pt. compare, but no NaN's can occur here */
        merge = ((ns <= nrelax1 && z < zrelax0) ||
            (ns <= nrelax2 && z < zrelax1) ||
            (z < zrelax2));

      }
    }

    if (merge)
    {
      Zeros [s] = totzeros ;
      Merged [s+1] = s ;
      Snz [s] = nscol0 + Snz [s+1] ;
      Nscol [s] += Nscol [s+1] ;
    }
  }

  /* contents of Wj no longer needed for Zeros ] */
  /* contents of Wi no longer needed for Nscol ] */
  /* contents of Sparent no longer needed (recomputed below) */

  /* ---------------------------------------------------------------------- */
  /* construct the relaxed supernode list */
  /* ---------------------------------------------------------------------- */

  nsuper = 0 ;
  for (s = 0 ; s < nfsuper ; s++)
  {
    if (Merged [s] == EMPTY)
    {
      Super [nsuper] = Super [s] ;
      Snz [nsuper] = Snz [s] ;
      IsSeparator[ nsuper ] = IsSeparator[ s ];
      nsuper++ ;
    }
  }
  Super [nsuper] = n ;

  /* Merged no longer needed ] */

  /* ---------------------------------------------------------------------- */
  /* find the mapping of relaxed nodes to supernodes */
  /* ---------------------------------------------------------------------- */

  /* use Wj as workspace for SuperMap { */

  /* SuperMap [k] = s if column k is contained in supernode s */
  for (s = 0 ; s < nsuper ; s++)
  {
    for (k = Super [s] ; k < Super [s+1] ; k++)
    {
      SuperMap [k] = s ;
    }
  }

  /* ---------------------------------------------------------------------- */
  /* construct the relaxed supernodal etree */
  /* ---------------------------------------------------------------------- */

  for (s = 0 ; s < nsuper ; s++)
  {
    j = Super [s+1] - 1 ;	/* last node in supernode s */
    parent = Parent [j] ;	/* parent of last node */
    Sparent [s] = (parent == EMPTY) ? EMPTY : SuperMap [parent] ;
  }

  /* ---------------------------------------------------------------------- */
  /* determine the size of L->s and L->x */
  /* ---------------------------------------------------------------------- */

  ssize = 0 ;
  xsize = 0 ;
  xxsize = 0 ;
  for (s = 0 ; s < nsuper ; s++)
  {
    nscol = Super [s+1] - Super [s] ;
#if 0
    cout << "Supernode " << s << " has " << nscol << " columns.";
    cout << " IsSeparator = " << IsSeparator[ s ] << endl;
#endif
    nsrow = Snz [s] ;
    ssize += nsrow ;
    xsize += nscol * nsrow ;
    /* also compute xsize in double to guard against Int overflow */
    xxsize += ((double) nscol) * ((double) nsrow) ;
  }
  xsize = max (1, xsize) ;
  ssize = max (1, ssize) ;

  /* ---------------------------------------------------------------------- */
  /* allocate L (all except real part L->x) */
  /* ---------------------------------------------------------------------- */

  cout << SDUMP( nsuper ) << endl;

  // Fill in the list of supernode specifications
  supernodes.clear();
  for ( int super_idx = 0; super_idx < nsuper; super_idx++ )
  {
    Ordering::SupernodeSpecification node;

    node._columnRange = IndexRange( Super[ super_idx ],
                                    Super[ super_idx + 1 ] - 1 );
    node._compressOffDiagonal = IsSeparator[ super_idx ];

    supernodes.push_back( node );
  }

  cout << SDUMP( supernodes.size() ) << endl;

#if 0
  L->ssize = ssize ;
  L->xsize = xsize ;
  L->nsuper = nsuper ;

  CHOLMOD(change_factor) (CHOLMOD_PATTERN, TRUE, TRUE, TRUE, TRUE, L, Common);

  if (Common->status < CHOLMOD_OK)
  {
    /* out of memory; L is still a valid simplicial symbolic factor */
    FREE_WORKSPACE ;
    return (FALSE) ;
  }

  DEBUG (CHOLMOD(dump_factor) (L, "L to symbolic super", Common)) ;
  ASSERT (L->is_ll && L->xtype == CHOLMOD_PATTERN && L->is_super) ;

  Lpi = L->pi ;
  Lpx = L->px ;
  Ls = L->s ;
  Ls [0] = 0 ;    /* flag for cholmod_check_factor; supernodes are defined */
  Lpx [0] = for_cholesky ? 0 : 123456 ;   /* magic number for sparse QR */
  Lsuper = L->super ;

  /* copy the list of relaxed supernodes into the final list in L */
  for (s = 0 ; s <= nsuper ; s++)
  {
    Lsuper [s] = Super [s] ;
  }

  /* Head no longer needed as workspace for fundamental Super list ) */

  Super = Lsuper ;	    /* Super is now the list of relaxed supernodes */

  /* ---------------------------------------------------------------------- */
  /* construct column pointers of relaxed supernodal pattern (L->pi) */
  /* ---------------------------------------------------------------------- */

  p = 0 ;
  for (s = 0 ; s < nsuper ; s++)
  {
    Lpi [s] = p ;
    p += Snz [s] ;
    PRINT1 (("Snz ["ID"] = "ID", Super ["ID"] = "ID"\n",
          s, Snz [s], s, Super[s])) ;
  }
  Lpi [nsuper] = p ;
  ASSERT ((Int) (L->ssize) == MAX (1,p)) ;

  /* ---------------------------------------------------------------------- */
  /* construct pointers for supernodal values (L->px) */
  /* ---------------------------------------------------------------------- */

  if (for_cholesky)
  {
    /* L->px is not needed for QR factorization (it may lead to Int
       overflow, anyway, if xsize caused Int overflow above) */
    p = 0 ;
    for (s = 0 ; s < nsuper ; s++)
    {
      nscol = Super [s+1] - Super [s] ;   /* number of columns in s */
      nsrow = Snz [s] ;           /* # of rows, incl triangular part*/
      Lpx [s] = p ;               /* pointer to numerical part of s */
      p += nscol * nsrow ;
    }
    Lpx [s] = p ;
    ASSERT ((Int) (L->xsize) == MAX (1,p)) ;
  }

  /* Snz no longer needed ] */

  /* ---------------------------------------------------------------------- */
  /* symbolic analysis to construct the relaxed supernodal pattern (L->s) */
  /* ---------------------------------------------------------------------- */

  Lpi2 = Wi ;	    /* copy Lpi into Lpi2, using Wi as workspace for Lpi2 [ */
  for (s = 0 ; s < nsuper ; s++)
  {
    Lpi2 [s] = Lpi [s] ;
  }

  Asorted = A->sorted ;

  for (s = 0 ; s < nsuper ; s++)
  {
    /* sth supernode is in columns k1 to k2-1.
     * compute nonzero pattern of L (k1:k2-1,:). */

    /* place rows k1 to k2-1 in leading column of supernode s */
    k1 = Super [s] ;
    k2 = Super [s+1] ;
    PRINT1 (("=========>>> Supernode "ID" k1 "ID" k2-1 "ID"\n",
          s, k1, k2-1)) ;
    for (k = k1 ; k < k2 ; k++)
    {
      Ls [Lpi2 [s]++] = k ;
    }

    /* compute nonzero pattern each row k1 to k2-1 */
    for (k = k1 ; k < k2 ; k++)
    {
      /* compute row k of L.  In the symmetric case, the pattern of L(k,:)
       * is the set of nodes reachable in the supernodal etree from any
       * row i in the nonzero pattern of A(0:k,k).  In the unsymmetric
       * case, the pattern of the kth column of A*A' is the set union
       * of all columns A(0:k,j) for each nonzero F(j,k). */

      /* clear the Flag array and mark the current supernode */
      /* mark = CHOLMOD(clear_flag) (Common) ; */
      CHOLMOD_CLEAR_FLAG (Common) ;
      mark = Common->mark ;
      Flag [s] = mark ;
      ASSERT (s == SuperMap [k]) ;

      /* traverse the row subtree for each nonzero in A or AA' */
      if (stype != 0)
      {
        subtree (k, k, Ap, Ai, Anz, SuperMap, Sparent, mark,
            Asorted, k1, Flag, Ls, Lpi2) ;
      }
      else
      {
        /* for each j nonzero in F (:,k) do */
        p = Fp [k] ;
        pend = (packed) ? (Fp [k+1]) : (p + Fnz [k]) ;
        for ( ; p < pend ; p++)
        {
          subtree (Fj [p], k, Ap, Ai, Anz, SuperMap, Sparent, mark,
              Asorted, k1, Flag, Ls, Lpi2) ;
        }
      }
    }
  }
#ifndef NDEBUG
  for (s = 0 ; s < nsuper ; s++)
  {
    PRINT1 (("Lpi2[s] "ID" Lpi[s+1] "ID"\n", Lpi2 [s], Lpi [s+1])) ;
    ASSERT (Lpi2 [s] == Lpi [s+1]) ;
    CHOLMOD(dump_super) (s, Super, Lpi, Ls, NULL, NULL, 0, Common) ;
  }
#endif

  /* contents of Wi no longer needed for Lpi2 ] */
  /* Sparent no longer needed ] */

  /* ---------------------------------------------------------------------- */
  /* determine the largest update matrix (L->maxcsize) */
  /* ---------------------------------------------------------------------- */

  /* maxcsize could be determined before L->s is allocated and defined, which
   * would mean that all memory requirements for both the symbolic and numeric
   * factorizations could be computed using O(nnz(A)+O(n)) space.  However, it
   * would require a lot of extra work.  The analysis phase, above, would need
   * to be duplicated, but with Ls not kept; instead, the algorithm would keep
   * track of the current s and slast for each supernode d, and update them
   * when a new row index appears in supernode d.  An alternative would be to
   * do this computation only if the allocation of L->s failed, in which case
   * the following code would be skipped.
   *
   * The csize for a supernode is the size of its largest contribution to
   * a subsequent ancestor supernode.  For example, suppose the rows of #'s
   * in the figure below correspond to the columns of a subsequent supernode,
   * and the dots are the entries in that ancestore.
   *
   *	    c
   *	    c c
   *	    c c c
   *	    x x x
   *	    x x x
   *	    # # #   .
   *	    # # #   . .
   *	    * * *   . .
   *	    * * *   . .
   *	    * * *   . .
   *	            . .
   *
   * Then for this update, the csize is 3-by-2, or 6, because there are 3
   * rows of *'s which is the number of rows in the update, and there are
   * 2 rows of #'s, which is the number columns in the update.  The csize
   * of a supernode is the largest such contribution for any ancestor
   * supernode.  maxcsize, for the whole matrix, has a rough upper bound of
   * the maximum size of any supernode.  This bound is loose, because the
   * the contribution must be less than the size of the ancestor supernodal
   * that it's updating.  maxcsize of a completely dense matrix, with one
   * supernode, is zero.
   *
   * maxesize is the column dimension for the workspace E needed for the
   * solve.  E is of size nrhs-by-maxesize, where the nrhs is the number of
   * columns in the right-hand-side.  The maxesize is the largest esize of
   * any supernode.  The esize of a supernode is the number of row indices
   * it contains, excluding the column indices of the supernode itself.
   * For the following example, esize is 4:
   *
   *	    c
   *	    c c
   *	    c c c
   *	    x x x
   *	    x x x
   *	    x x x
   *	    x x x
   *
   * maxesize can be no bigger than n.
   */

  maxcsize = 1 ;
  maxesize = 1 ;

  /* Do not need to guard csize against Int overflow since xsize is OK. */

  if (for_cholesky)
  {
    /* this is not needed for QR factorization */
    for (d = 0 ; d < nsuper ; d++)
    {
      nscol = Super [d+1] - Super [d] ;
      p = Lpi [d] + nscol ;
      plast = p ;
      pend = Lpi [d+1] ;
      esize = pend - p ;
      maxesize = MAX (maxesize, esize) ;
      slast = (p == pend) ? (EMPTY) : (SuperMap [Ls [p]]) ;
      for ( ; p <= pend ; p++)
      {
        s = (p == pend) ? (EMPTY) : (SuperMap [Ls [p]]) ;
        if (s != slast)
        {
          /* row i is the start of a new supernode */
          ndrow1 = p - plast ;
          ndrow2 = pend - plast ;
          csize = ndrow2 * ndrow1 ;
          PRINT1 (("Supernode "ID" ancestor "ID" C: "ID"-by-"ID
                "  csize "ID"\n", d, slast, ndrow1, ndrow2, csize)) ;
          maxcsize = MAX (maxcsize, csize) ;
          plast = p ;
          slast = s ;
        }
      }
    }
    PRINT1 (("max csize "ID"\n", maxcsize)) ;
  }

  /* Wj no longer needed for SuperMap } */

  L->maxcsize = maxcsize ;
  L->maxesize = maxesize ;
  L->is_super = TRUE ;
  ASSERT (L->xtype == CHOLMOD_PATTERN && L->is_ll) ;

  /* ---------------------------------------------------------------------- */
  /* supernodal symbolic factorization is complete */
  /* ---------------------------------------------------------------------- */

  FREE_WORKSPACE ;
  return (TRUE) ;
#endif
}
