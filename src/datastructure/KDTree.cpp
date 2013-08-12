//////////////////////////////////////////////////////////////////////
// KDTree.cpp: Implementation of the KDTree class
//
//////////////////////////////////////////////////////////////////////

#include "KDTree.h"

#include <algorithm>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
template <class T>
KDTree<T>::KDTree( KDTree::SplitType splitType, int maxLeafSize )
	: _root( NULL ),
    _splitType( splitType ),
    _maxLeafSize( maxLeafSize )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
template <class T>
KDTree<T>::~KDTree()
{
	delete _root;
}

//////////////////////////////////////////////////////////////////////
// Build's a KD tree out of the given data
//////////////////////////////////////////////////////////////////////
template <class T>
void KDTree<T>::build( vector<T *> &data )
{
	if ( !_root )
	{
		delete _root;
	}

	_root = new KDTree<T>::KDNode( KDTree::SPLIT_X );

	_root->build( data, _splitType, _maxLeafSize );
}

//////////////////////////////////////////////////////////////////////
// Intersect the given point with this tree, and append to a
// set of T objects (possibly) intersecting with this point
//////////////////////////////////////////////////////////////////////
template <class T>
void KDTree<T>::intersect( const VEC3F &x, set<T *> &entries )
{
	if ( _root )
	{
		_root->intersect( x, entries );
	}
}

//////////////////////////////////////////////////////////////////////
// Same as the above, except that elements are put in to a
// vector.  Note that this allows duplicate entries.
//////////////////////////////////////////////////////////////////////
template <class T>
void KDTree<T>::intersect( const VEC3F &x, vector<T *> &entries )
{
	if ( _root )
	{
		_root->intersect( x, entries );
	}
}

//////////////////////////////////////////////////////////////////////
// Returns an ordering of the nodes in the tree, via a left
// to right traversal along the bottom level.
//////////////////////////////////////////////////////////////////////
template <class T>
void KDTree<T>::orderNodes( std::vector<T *> &entries )
{
  entries.clear();

  _root->orderNodes( entries );
}

//////////////////////////////////////////////////////////////////////
// Returns a list of node sizes for each node in the tree.
//////////////////////////////////////////////////////////////////////
template <class T>
void KDTree<T>::nodeSizes( IntArray &sizes ) const
{
  sizes.clear();
  _root->nodeSizes( sizes, 0 );
}

//////////////////////////////////////////////////////////////////////
// KDNode class definitions
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
template <class T>
KDTree<T>::KDNode::KDNode( KDSplit split )
	: _split( split ),
		_left( NULL ),
		_right( NULL ),
		_data( 0 )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
template <class T>
KDTree<T>::KDNode::~KDNode()
{
	delete _left;
	delete _right;
}

//////////////////////////////////////////////////////////////////////
// Build's a tree with the given data
//////////////////////////////////////////////////////////////////////
template <class T>
void KDTree<T>::KDNode::build( vector<T *> &data,
                               KDTree<T>::SplitType splitType,
                               int maxLeafSize )
{
	_numElements = data.size();

	// Grow the bounding box for this node to accomodate
	// its data
	for ( int i = 0; i < data.size(); i++ )
	{
		_bbox += data[i]->bbox();
	}

  // Override our split direction if we want to split
  // along the longest box axis.
  if ( splitType == LONGEST_AXIS )
  {
    switch ( _bbox.longestAxis() )
    {
      case BoundingBox::X_AXIS:
      default:
      {
        _split = SPLIT_X;
        break;
      }
      case BoundingBox::Y_AXIS:
      {
        _split = SPLIT_Y;
        break;
      }
      case BoundingBox::Z_AXIS:
      {
        _split = SPLIT_Z;
        break;
      }
    }
  }

	// If we are down to the maximum leaf size then we're done
	if ( data.size() <= maxLeafSize )
	{
		_data = data;
		return;
	}

	int										 splitPoint;
	vector<T *>						 leftData, rightData;

	// Sort the data along the split axis direction
	KDTree<T>::KDCompare	 compare( _split );

	sort( data.begin(), data.end(), compare );

	splitPoint = data.size() / 2;

	leftData.resize( splitPoint );
	rightData.resize( data.size() - splitPoint );

	// Create data arrays for our left and right children,
	// and grow our bounding box while we're at it.
	for ( int i = 0; i < splitPoint; i++ )
	{
		leftData[i] = data[i];
	}

	for ( int i = splitPoint; i < data.size(); i++ )
	{
		rightData[i - splitPoint] = data[i];
	}

	// If we haven't reduced the size of our subtrees, then
	// we stop here (don't want to recurse indefinitely, after all).
	if ( leftData.size() == data.size() )
	{
		_data = data;
		return;
	}

	// Create children for ourselves
	_left = new KDTree<T>::KDNode( (KDSplit)( ((int)_split + 1) % 3 ) );
	_right = new KDTree<T>::KDNode( (KDSplit)( ((int)_split + 1) % 3 ) );

	_left->build( leftData, splitType, maxLeafSize );
	_right->build( rightData, splitType, maxLeafSize );
}

//////////////////////////////////////////////////////////////////////
// Intersect a point with this KD node, returning a set of T pointers
// which (possibly) intersect with the point
//////////////////////////////////////////////////////////////////////
template <class T>
void KDTree<T>::KDNode::intersect( const VEC3F &x, set<T *> &entries )
{
	// If we are at a leaf node, check all entries against the
	// point and insert in to the set, if necessary.
	if ( _data.size() != 0 )
	{
		for ( int i = 0; i < _data.size(); i++ )
		{
			if ( _data[i]->bbox().isInside( x ) )
			{
				entries.insert( _data[i] );
			}
		}

		return;
	}

	// Figure out which children we have to intersect with.
	if ( _left->bbox().isInside( x ) )
	{
		_left->intersect( x, entries );
	}

	if ( _right->bbox().isInside( x ) )
	{
		_right->intersect( x, entries );
	}
}

//////////////////////////////////////////////////////////////////////
// Same as the above, except that elements are put in to a
// vector.  Note that this allows duplicate entries.
//////////////////////////////////////////////////////////////////////
template <class T>
void KDTree<T>::KDNode::intersect( const VEC3F &x, vector<T *> &entries )
{
	// If we are at a leaf node, check all entries against the
	// point and insert in to the set, if necessary.
	if ( _data.size() != 0 )
	{
		for ( int i = 0; i < _data.size(); i++ )
		{
			if ( _data[i]->bbox().isInside( x ) )
			{
				entries.push_back( _data[i] );
			}
		}

		return;
	}

	// Figure out which children we have to intersect with.
	if ( _left->bbox().isInside( x ) )
	{
		_left->intersect( x, entries );
	}

	if ( _right->bbox().isInside( x ) )
	{
		_right->intersect( x, entries );
	}
}

//////////////////////////////////////////////////////////////////////
// Returns an ordering of the nodes in the tree, via a left
// to right traversal along the bottom level.
//////////////////////////////////////////////////////////////////////
template <class T>
void KDTree<T>::KDNode::orderNodes( std::vector<T *> &entries )
{
  // If we're at a leaf node, copy its contents
  if ( _data.size() != 0 )
  {
    for ( int i = 0; i < _data.size(); i++ )
    {
      entries.push_back( _data[i] );
    }

    return;
  }

  _left->orderNodes( entries );
  _right->orderNodes( entries );
}

//////////////////////////////////////////////////////////////////////
// Returns a list of node sizes for each node in the tree.
//////////////////////////////////////////////////////////////////////
template <class T>
void KDTree<T>::KDNode::nodeSizes( IntArray &sizes, int idx ) const
{
  while ( idx >= sizes.size() )
  {
    sizes.push_back( 0 );
  }

  sizes[ idx ] = _numElements;

  if ( _left ) _left->nodeSizes( sizes, 2 * idx + 1 );
  if ( _right ) _right->nodeSizes( sizes, 2 * idx + 2 );
}
