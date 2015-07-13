//////////////////////////////////////////////////////////////////////
// SeparatorTree.cpp: Implementation of the SeparatorTree class
//
//////////////////////////////////////////////////////////////////////

#include "SeparatorTree.h"

#include <algorithm>

#include <util/IO.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
template <class T>
SeparatorTree<T>::SeparatorTree( SeparatorTree::SplitType splitType,
                                 int maxLeafSize,
                                 SeparatorTree::SeparatorType separatorType,
                                 Real separatorThickness,
                                 int maxLevels )
	: _root( NULL ),
    _splitType( splitType ),
    _maxLeafSize( maxLeafSize ),
    _maxLevels( maxLevels ),
    _separatorType( separatorType ),
    _separatorThickness( separatorThickness )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
template <class T>
SeparatorTree<T>::~SeparatorTree()
{
	delete _root;
}

//////////////////////////////////////////////////////////////////////
// Build's a KD tree out of the given data
//////////////////////////////////////////////////////////////////////
template <class T>
void SeparatorTree<T>::build( vector<T *> &data, bool storeSeparator,
                              bool fullTree )
{
	if ( !_root )
	{
		delete _root;
	}

	_root = new SeparatorTree<T>::KDNode( SeparatorTree::SPLIT_X,
                                        _separatorType, _separatorThickness );

	_root->build( data, _splitType, _maxLeafSize, _maxLevels, storeSeparator,
                fullTree );
}

//////////////////////////////////////////////////////////////////////
// Intersect the given point with this tree, and append to a
// set of T objects (possibly) intersecting with this point
//////////////////////////////////////////////////////////////////////
template <class T>
void SeparatorTree<T>::intersect( const VEC3F &x, set<T *> &entries )
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
void SeparatorTree<T>::intersect( const VEC3F &x, vector<T *> &entries )
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
void SeparatorTree<T>::orderNodes( std::vector<T *> &entries )
{
  entries.clear();

  _root->orderNodes( entries );
}

//////////////////////////////////////////////////////////////////////
// Returns a list of node sizes for each node in the tree.
//////////////////////////////////////////////////////////////////////
template <class T>
void SeparatorTree<T>::nodeSizes( IntArray &sizes,
                                  bool excludeSeparator ) const
{
  sizes.clear();
  _root->nodeSizes( sizes, 0, excludeSeparator );
}

//////////////////////////////////////////////////////////////////////
// List of ranges for each leaf in the separator tree
//////////////////////////////////////////////////////////////////////
template <class T>
void SeparatorTree<T>::leafRanges( vector<IndexRange> &sizes ) const
{
  sizes.clear();
  _root->leafRanges( sizes );
}

//////////////////////////////////////////////////////////////////////
// KDNode class definitions
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
template <class T>
SeparatorTree<T>::KDNode::KDNode( KDSplit split,
                                  SeparatorTree<T>::SeparatorType separatorType,
                                  Real separatorThickness )
	: _split( split ),
		_left( NULL ),
		_right( NULL ),
    _center( NULL ),
		_data( 0 ),
    _separatorType( separatorType ),
    _separatorThickness( separatorThickness )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
template <class T>
SeparatorTree<T>::KDNode::~KDNode()
{
	delete _left;
	delete _right;
}

//////////////////////////////////////////////////////////////////////
// Build's a tree with the given data
//////////////////////////////////////////////////////////////////////
template <class T>
int SeparatorTree<T>::KDNode::build( vector<T *> &data,
                                     SeparatorTree<T>::SplitType splitType,
                                     int maxLeafSize, int maxLevels,
                                     bool storeSeparator, bool fullTree )
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
	if ( data.size() <= maxLeafSize || maxLevels == 0 )
	{
		_data = data;
		return 0;
	}

	vector<T *>						 leftData, centerData, rightData;
  int                    leftSize, centerSize, rightSize;

	// Sort the data along the split axis direction
	SeparatorTree<T>::KDCompare	 compare( _split );

	sort( data.begin(), data.end(), compare );

#if 0
  // FIXME:
  cout << "Data after sorting" << endl;
  cout << "Split direction was " << _split << endl;
  for ( int i = 0; i < data.size(); i++ ) {
    cout << data[i]->bbox().center() << endl;
  }
  abort();
#endif

  findChildSizes( data, leftSize, centerSize, rightSize );

  if ( storeSeparator )
  {
    leftData.resize( leftSize );
    centerData.resize( centerSize );
  }
  else
  {
    leftData.resize( leftSize + centerSize );
  }
	rightData.resize( rightSize );

	// Create data arrays for our left and right children,
	// and grow our bounding box while we're at it.
  if ( storeSeparator )
  {
    for ( int i = 0; i < leftSize; i++ )
    {
      leftData[i] = data[i];
    }

    for ( int i = leftSize; i < leftSize + centerSize; i++ )
    {
      centerData[i - leftSize] = data[i];
    }
  }
  else
  {
    for ( int i = 0; i < leftSize + centerSize; i++ )
    {
      leftData[i] = data[i];
    }
  }

  for ( int i = leftSize + centerSize; i < data.size(); i++ )
  {
    rightData[ i - leftSize - centerSize ] = data[i];
  }

	// If we haven't reduced the size of our subtrees, then
	// we stop here (don't want to recurse indefinitely, after all).
	if ( leftData.size() == data.size() || rightData.size() == 0 )
	{
		_data = data;
		return 0;
	}

	// Create children for ourselves
	_left = new SeparatorTree<T>::KDNode( (KDSplit)( ((int)_split + 1) % 3 ),
                                        _separatorType, _separatorThickness );
	_right = new SeparatorTree<T>::KDNode( (KDSplit)( ((int)_split + 1) % 3 ),
                                         _separatorType, _separatorThickness );

  if ( storeSeparator )
  {
    _center = new SeparatorTree<T>::KDNode(
                                (KDSplit)( ((int)_split + 1) % 3 ),
                                _separatorType, _separatorThickness );
  }
  else
  {
    _center = NULL;
  }

  int                        myHeight = 1;
  int                        leftHeight;
  int                        rightHeight;

  if ( fullTree &&
       ( rightData.size() > maxLeafSize && leftData.size() <= maxLeafSize ) )
  {
    // This is a real hack
    leftHeight = _left->build( leftData, splitType,
                                leftData.size() - 1 /* Override size */,
                                maxLevels - 1, storeSeparator, fullTree );
  }
  else
  {
    leftHeight = _left->build( leftData, splitType, maxLeafSize,
                                maxLevels - 1, storeSeparator, fullTree );
  }

  if ( fullTree &&
       ( leftData.size() > maxLeafSize && rightData.size() <= maxLeafSize ) )
  {
    rightHeight = _right->build( rightData, splitType,
                                 rightData.size() - 1 /* Override size */,
                                 maxLevels - 1, storeSeparator, fullTree );
  }
  else
  {
    rightHeight = _right->build( rightData, splitType, maxLeafSize,
                                 maxLevels - 1, storeSeparator, fullTree );
  }

  if ( fullTree && leftHeight != rightHeight )
  {
    cout << SDUMP( leftHeight ) << SDUMP( rightHeight ) << endl;
    cout << SDUMP( leftData.size() ) << SDUMP( rightData.size() ) << endl;
    cout << SDUMP( maxLeafSize ) << endl;
  }
  TRACE_ASSERT( !fullTree || leftHeight == rightHeight,
                "Full tree child mismatch" );

  myHeight = 1 + max( leftHeight, rightHeight );

  if ( storeSeparator )
  {
    _center->build( centerData, splitType, maxLeafSize, maxLevels - 1, true,
                    false );
  }

  return myHeight;
}

//////////////////////////////////////////////////////////////////////
// Intersect a point with this KD node, returning a set of T pointers
// which (possibly) intersect with the point
//////////////////////////////////////////////////////////////////////
template <class T>
void SeparatorTree<T>::KDNode::intersect( const VEC3F &x, set<T *> &entries )
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

  if ( _center->bbox().isInside( x ) )
  {
    _center->intersect( x, entries );
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
void SeparatorTree<T>::KDNode::intersect( const VEC3F &x, vector<T *> &entries )
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

  if ( _center->bbox().isInside( x ) )
  {
    _center->intersect( x, entries );
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
void SeparatorTree<T>::KDNode::orderNodes( std::vector<T *> &entries )
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

  if ( _left ) _left->orderNodes( entries );
  if ( _right ) _right->orderNodes( entries );
  if ( _center ) _center->orderNodes( entries );
}

//////////////////////////////////////////////////////////////////////
// Returns a list of node sizes for each node in the tree.
//////////////////////////////////////////////////////////////////////
template <class T>
void SeparatorTree<T>::KDNode::nodeSizes( IntArray &sizes, int idx,
                                          bool excludeSeparator ) const
{
  while ( idx >= sizes.size() )
  {
    sizes.push_back( 0 );
  }

  sizes[ idx ] = _numElements;

  if ( excludeSeparator )
  {
    // Treat this more as a binary tree, ignoring the size
    // of the middle (separator) node
    if ( _left ) _left->nodeSizes( sizes, 2 * idx + 1, excludeSeparator );
    if ( _right ) _right->nodeSizes( sizes, 2 * idx + 2, excludeSeparator );
  }
  else
  {
    if ( _left ) _left->nodeSizes( sizes, 3 * idx + 1 );
    if ( _right ) _right->nodeSizes( sizes, 3 * idx + 2 );
    if ( _center ) _center->nodeSizes( sizes, 3 * idx + 3 );
  }
}

//////////////////////////////////////////////////////////////////////
// Returns a list of leaf node sizes for the leafs of this tree
//////////////////////////////////////////////////////////////////////
template <class T>
void SeparatorTree<T>::KDNode::leafRanges( vector<IndexRange> &ranges ) const
{
  TRACE_ASSERT( !_center, "Doesn't work if we store the separator" );

  if ( !_left )
  {
    TRACE_ASSERT( !_right, "Something is wrong" );

    if ( ranges.size() == 0 )
    {
      ranges.push_back( IndexRange( 0, _data.size() - 1 ) );
    }
    else
    {
      ranges.push_back( IndexRange( ranges.back().second + 1,
                                    ranges.back().second + _data.size() ) );
    }
  }
  else
  {
    _left->leafRanges( ranges );
    _right->leafRanges( ranges );
  }
}

//////////////////////////////////////////////////////////////////////
// Returns indices corresponding to the sizes of the
// left, right and center chunks of the given index list.
//////////////////////////////////////////////////////////////////////
template <class T>
void SeparatorTree<T>::KDNode::findChildSizes( std::vector<T *> &data,
                                               int &leftSize,
                                               int &centerSize,
                                               int &rightSize ) const
{
  int          splitPoint = data.size() / 2;
  int          idx = -1;
  Real         dist = FLT_MAX;

  Real         separatorDistance;

  switch ( _separatorType )
  {
    case FIXED:
    default:
    {
      separatorDistance = _separatorThickness;
      break;
    }
    case LINEAR_LENGTH:
    {
      separatorDistance = _bbox.axislength( _split ) * _separatorThickness;
      break;
    }
  }

  do
  {
    dist = data[ splitPoint ]->bbox().center()[ _split ]
         - data[ idx + 1 ]->bbox().center()[ _split ];
    idx++;
  } while ( dist > separatorDistance );

  leftSize = idx;

  if ( separatorDistance == 0.0 )
  {
    centerSize = 0;
  }
  else
  {
    while ( abs( dist ) <= separatorDistance )
    {
      dist = data[ splitPoint ]->bbox().center()[ _split ]
           - data[ idx + 1 ]->bbox().center()[ _split ];
      idx++;
    }

    centerSize = idx - leftSize;
  }

  // FIXME
  leftSize = data.size() / 2;
  centerSize = 0;

  rightSize = data.size() - leftSize - centerSize;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class T>
int SeparatorTree<T>::KDNode::maxDepth()
{
  int                        maxDepth = 0;

  if ( _left )
  {
    maxDepth = _left->maxDepth();
  }

  if ( _right )
  {
    maxDepth = max( maxDepth, _right->maxDepth() );
  }

  return maxDepth;
}
