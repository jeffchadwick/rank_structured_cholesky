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
// KDTree.h: Interface for the KDTree class
//
//////////////////////////////////////////////////////////////////////

#ifndef KD_TREE_H
#define KD_TREE_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <rschol/geometry/BoundingBox.h>

#include <set>
#include <vector>

//////////////////////////////////////////////////////////////////////
// KDTree class
//
// A standard KD-tree for partitioning objects (possibly volumetric
// objects like tetrahedra, etc.).
//
// All classes T used in this template must provide the function
// BoundingBox & bbox();
// which returns the bounding box for the object.
// Single points may simply return a degenerate bounding box.
//////////////////////////////////////////////////////////////////////
template <class T>
class KDTree {
	public:
		// KD Tree splitting directions.
		enum KDSplit {
			SPLIT_X = 0,
			SPLIT_Y,
			SPLIT_Z
		};

    // The type of splitting to do.  We can either cycle through
    // splitting directions, or always split along the longest
    // axis of each node.
    enum SplitType {
      CYCLICAL = 0,
      LONGEST_AXIS
    };

		class KDCompare {
			public:
				KDCompare( KDSplit split )
					: _split( split )
				{
				}

				virtual ~KDCompare() {}

				inline bool								 operator()( T *t1, T *t2 )
																	 {
																		 return ( t1->bbox().center()[_split]
																		 			< t2->bbox().center()[_split] );
																	 }

			private:
				// The index to compare on
				KDSplit									 _split;
		};

		// Tree nodes.  These do most of the work.
		class KDNode {
			public:
				KDNode( KDSplit split );

				// Destructor
				virtual ~KDNode();

				// Build a new node given a set of objects
				void										 build( std::vector<T *> &data,
                                        SplitType splitType = CYCLICAL,
                                        int maxLeafSize = 1 );

				inline BoundingBox			&bbox()
																 {
																	 return _bbox;
																 }

				// Intersect the given point with this node, and append to a
				// set of T objects (possibly) intersecting with this point
				void										 intersect( const VEC3F &x,
																						std::set<T *> &entries );

				// Same as the above, except that elements are put in to a
				// vector.  Note that this allows duplicate entries.
				void										 intersect( const VEC3F &x,
																						std::vector<T *> &entries );

        // Returns an ordering of the nodes in the tree, via a left
        // to right traversal along the bottom level.
        void                     orderNodes( std::vector<T *> &entries );

        // Returns a list of node sizes for each node in the tree.
        void                     nodeSizes( IntArray &sizes, int idx ) const;

        // Accessors
        int                      numElements() const { return _numElements; }
        const KDNode            *left() const { return _left; }
        const KDNode            *right() const { return _right; }

			private:
				// Bounds for all objects stored in this
				// node.
				BoundingBox							 _bbox;

				// This node's splitting direction
				KDSplit									 _split;

				// Children
				KDNode									*_left;
				KDNode									*_right;

				// List of data for this node - only non-empty at leaves
				std::vector<T *>				 _data;

				// How many elements are stored beneath this node
				int											 _numElements;

		};

		KDTree( SplitType splitType = CYCLICAL, int maxLeafSize = 1 );

		// Destructor
		virtual ~KDTree();

		// Build's a KD tree out of the given data
		void												 build( std::vector<T *> &data );

		// Intersect the given point with this tree, and append to a
		// set of T objects (possibly) intersecting with this point
		void												 intersect( const VEC3F &x,
																						std::set<T *> &entries );

		// Same as the above, except that elements are put in to a
		// vector.  Note that this allows duplicate entries.
		void												 intersect( const VEC3F &x,
																						std::vector<T *> &entries );

    // Returns an ordering of the nodes in the tree, via a left
    // to right traversal along the bottom level.
    void                         orderNodes( std::vector<T *> &entries );

    // Returns a list of node sizes for each node in the tree.
    void                         nodeSizes( IntArray &sizes ) const;

    // Returns the root
    const KDNode                *root() const { return _root; }

	protected:

	private:
		KDNode											*_root;
    SplitType                    _splitType;
    int                          _maxLeafSize;

};

#include <rschol/datastructure/KDTree.cpp>

#endif
