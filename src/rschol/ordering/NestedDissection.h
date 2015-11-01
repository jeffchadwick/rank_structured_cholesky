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
// NestedDissection.h: Interface for the NestedDissection class
//
//////////////////////////////////////////////////////////////////////

#ifndef NESTED_DISSECTION_H
#define NESTED_DISSECTION_H

#include <rschol/library_interface/CHOLMOD_Environment.h>
#include <rschol/library_interface/CSparse_Interface.h>

#include <rschol/linearalgebra/DenseBlock.h>
#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/mesh/Mesh.h>

#include <rschol/util/IO.h>

#include <rschol/ordering/Ordering.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <boost/bind.hpp>
#include <boost/function.hpp>

//////////////////////////////////////////////////////////////////////
// NestedDissection class
//
// Provides data structures and routines for handling nested
// dissection orderings of matrices.
//////////////////////////////////////////////////////////////////////
class NestedDissection {
	public:
		NestedDissection();

		// Destructor
		virtual ~NestedDissection();

    // Tree node class for modelling nested dissection orderings
    class DissectionNode {
      public:
        // Constructor
        DissectionNode( int variableDimensions );

        // Desctructor
        virtual ~DissectionNode();

        // Builds a leaf node, which has no separators or children
        // defined.
        static DissectionNode *buildLeafNode( IntArray &indices,
                                              int variableDimensions );

        // Builds an interior node, which has a separator and children
        // defined.
        static DissectionNode *buildInteriorNode( IntArray &separator,
                                                  DissectionNode *left,
                                                  DissectionNode *right,
                                                  int variableDimensions );

        // Determines the indices in a system reordered according to this
        // tree structure that would be associated with this node.
        int buildReorderedIndices( int start, int variableDimensions );

        // Builds a permutation vector based on this nested dissection
        // ordering.
        void buildPermutation( IntArray &permutation,
                               int variableDimensions );

        // Generates an array of all leaf nodes in this subtree
        void getLeafNodes( std::vector<const DissectionNode *> &leafNodes ) const;

        // Generates a post-ordering of this subtree
        void postOrder( std::vector<const DissectionNode *> &nodes ) const;

        IndexRange reorderedRange() const { return _reorderedRange; }

        const DissectionNode *left() const { return _left; }
        const DissectionNode *right() const { return _right; }

        bool isLeaf() const { return _left == NULL && _right == NULL; }

        // For a leaf node, return the number of entries stored in
        // the leaf node, otherwise, store the number of separator
        // entries
        int size() const
        {
          if ( isLeaf() )
          {
            return _indices.size() * _variableDimensions;
          }
          else
          {
            return _separatorIndices.size() * _variableDimensions;
          }
        }

      private:
        // Children
        DissectionNode    *_left;
        DissectionNode    *_right;

        // The indices associated with this node.
        // This should only be non-empty for leaf nodes.
        IntArray           _indices;

        // Which indices are associated with the separator
        // This should only be non-empty for interior nodes.
        IntArray           _separatorIndices;

        // The range of indices associated with this node assuming
        // that the associated system has been reordered using
        // nested dissection.  Call buildReorderedIndices
        // to construct this.
        IndexRange         _reorderedRange;

        int                _variableDimensions;
    };

    // Tree class for modelling nested dissection orderings
    class DissectionTree {
      public:
        // Constructor
        DissectionTree();

        // Destructor
        virtual ~DissectionTree();

        DissectionNode * &root() { return _root; }
        const DissectionNode *root() const { return _root; }

        void buildReorderedIndices( int variableDimensions = 1 )
        {
          _root->buildReorderedIndices( 0, variableDimensions );
        }

        void buildPermutation( IntArray &permutation,
                               int variableDimensions = 1 )
        {
          _root->buildPermutation( permutation, variableDimensions );
        }

        // Generates an array of all leaf nodes in the tree
        void getLeafNodes( std::vector<const DissectionNode *> &leafNodes ) const
        {
          leafNodes.clear();
          _root->getLeafNodes( leafNodes );
        }

        // Generates a post-ordering of the tree (which corresponds
        // to the desired ordering of index blocks in the final
        // matrix)
        void postOrder( std::vector<const DissectionNode *> &nodes ) const
        {
          nodes.clear();
          _root->postOrder( nodes );
        }

      private:
        DissectionNode    *_root;
    };

    // Builds a nested dissection ordering tree for a
    // rectangular finite difference grid of the
    // given size.
    static DissectionTree *nestedDissectionFD( Tuple3i gridDivisions,
                                               int maxBlockSize,
                                               int maxLevels,
                                               bool splitLargest = true,
                                               int variableDimensions = 1,
                                               int maxSeparatorSize = 0 );

    // Converts a nested dissection node list in to a list of supernode
    // specifications (directly assuming that separators correspond to
    // supernodes)
    //
    // Assumes that reordered column ranges have been computed for the
    // nodes
    //
    // maxBlockSize specifies the size above which off diagonals in the
    // supernodes should be compressed.
    static void ConvertNodeList(
                  const std::vector<const DissectionNode *> &nodeList,
                  int maxBlockSize,
                  std::vector<Ordering::SupernodeSpecification> &supernodes );

    typedef boost::function<Mesh::MeshPoint * (int index)> PointBuilder;

    // Performs nested dissection on an arbitraty matrix, given a way
    // of converting point indices in to mesh points
    static void SupernodalNestedDissection(
                    const SPARSE_MATRIX::SparseColumnMatrix &A,
                    PointBuilder &pointPositionQuery,
                    IntArray &permutation,
                    std::vector<Ordering::SupernodeSpecification> &supernodes,
                    int maxBlockSize, int maxDiagonalBlockSize,
                    CompressionType compressionType,
                    std::vector<std::vector<DenseBlock> > &diagonalBlocks,
                    // If this is true, then rather then splitting nodes
                    // hierarchically, we simply represent each diagonal
                    // block as a separate supernode.
                    bool splitNodes = false,
                    // Whether or not to apply sub-ordering to separator
                    // diagonal blocks.  Generally this should be true,
                    // assuming that pointPositionQuery provides useful
                    // information about DoF positions
                    bool reorderSeparators = true );

    // Given a post ordering of a dissection tree from a finite difference
    // grid computed
    // with the function above, reorders separators in the tree according
    // to a maximum block size.  This is done to expose low rank behaviour
    // on the diagonals of blocks formed for these separators.
    //
    // This function adjusts the given permutation according to these
    // reorderings.  It also provides a list of diagonal blocks for each
    // node in the tree (though this list will be empty for many nodes).
    static void ReorderNestedSeparators(
            Tuple3i gridDivisions, Real h,
            const std::vector<const DissectionNode *> &nodeList,
            int maxBlockSize,
            IntArray &permutation,
            std::vector<std::vector<DenseBlock> > &diagonalBlocks );

    // Same as the above, but uses a mesh with explicit vertex positions
    static void ReorderNestedSeparators(
            const Mesh &mesh,
            const std::vector<const DissectionNode *> &nodeList,
            int maxBlockSize,
            IntArray &permutation,
            std::vector<std::vector<DenseBlock> > &diagonalBlocks,
            int variableDimensions = 1 );

    // Same as the above, but node positions will be specified implicitly
    // by the given PointBuilder function
    static void ReorderNestedSeparators(
            const PointBuilder &pointPositionQuery,
            const std::vector<Ordering::SupernodeSpecification> &supernodes,
            int maxBlockSize,
            CompressionType compressionType,
            IntArray &permutation,
            std::vector<std::vector<DenseBlock> > &diagonalBlocks,
            int variableDimensions = 1,
            // Whether or not to apply sub-ordering to separator
            // diagonal blocks.  Generally this should be true,
            // assuming that pointPositionQuery provides useful
            // information about DoF positions
            bool reorderSeparators = true );

    // Modifies supernode set so that each diagonal block is explicitly
    // treated as a supernode.
    static void SplitNodes(
            std::vector<Ordering::SupernodeSpecification> &supernodes,
            std::vector<std::vector<DenseBlock> > &diagonalBlocks );

    // Splits nodes according to some maximum size
    static void SplitNodes(
            std::vector<Ordering::SupernodeSpecification> &supernodes,
            std::vector<std::vector<DenseBlock> > &diagonalBlocks,
            int maxSize );

    // Splits a single node according to some maximum size constraint
    static void SplitNode(
            const Ordering::SupernodeSpecification &nodeSpec,
            const std::vector<DenseBlock> &diagonalBlocks,
            std::vector<Ordering::SupernodeSpecification> &newNodes,
            std::vector<std::vector<DenseBlock> > &newBlocks,
            int maxSize );

	protected:

	private:
    // Recursive finite difference nested dissection function.
    static DissectionNode *nestedDissectionFD( Tuple3i gridDivisions,
                                               Tuple3i gridStart,
                                               Tuple3i gridEnd,
                                               int maxBlockSize,
                                               int maxLevels,
                                               bool splitLargest,
                                               int variableDimensions,
                                               int maxSeparatorSize );

    // Builds supernodes out of an ordering (with column counts).
    // This code is heavily adapted from CHOLMOD's super_symbolic function,
    // which performs supernode relaxation.
    //
    // Specify the list of separators from which the permutation is
    // built.  Also indicate the maximum block size that we will leave
    // uncompressed.  Any separators larger than this size will be treated
    // as supernodes themselves.
    static void BuildRelaxedSupernodes(
                const IntArray &Parent, const IntArray &columnCounts,
                const cs *Apermuted,
                const CHOLMOD_Environment::SeparatorList &separators,
                const IntArray &separatorParents,
                const IntArray &separatorLevels,
                int maxBlockSize,
                std::vector<Ordering::SupernodeSpecification> &supernodes,
                int nrelax0 = 4, int nrelax1 = 16, int nrelax2 = 48,
                Real zrelax0 = 0.8, Real zrelax1 = 0.1, Real zrelax2 = 0.05 );

};

#endif
