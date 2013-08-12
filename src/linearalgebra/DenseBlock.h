//////////////////////////////////////////////////////////////////////
// DenseBlock.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef DENSE_BLOCK_H
#define DENSE_BLOCK_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/VECTOR.h>

#include <SETTINGS.h>
#include <TYPES.h>

//////////////////////////////////////////////////////////////////////
// DenseBlock class
//
//////////////////////////////////////////////////////////////////////
class DenseBlock {
	public:
    DenseBlock()
      : _rowRange( 0, 0 ),
        _columnRange( 0, 0 )
    {
    }

    // Destructor
    virtual ~DenseBlock()
    {
    }

    DenseBlock( const IndexRange &rowRange, const IndexRange &columnRange )
      : _rowRange( rowRange ),
        _columnRange( columnRange )
    {
    }

    DenseBlock( int startRow, int endRow, int startCol, int endCol )
      : _rowRange( startRow, endRow ),
        _columnRange( startCol, endCol )
    {
    }

    inline bool valid() const
    {
      return ( _rowRange.first >= 0
            && _rowRange.second >= _rowRange.first
            && _columnRange.first >= 0
            && _columnRange.second >= _columnRange.first );
    }

    inline int numRows() const
    {
      return range_size( _rowRange );
    }

    inline int numColumns() const
    {
      return range_size( _columnRange );
    }

    inline bool containsRow( int row ) const
    {
      return ( row >= _rowRange.first && row <= _rowRange.second );
    }

    inline bool containsColumn( int column ) const
    {
      return ( column >= _columnRange.first && column <= _columnRange.second );
    }

    inline bool isDiagonal() const
    {
      return _rowRange == _columnRange;
    }

    IndexRange             _rowRange;
    IndexRange             _columnRange;

	protected:

	private:
};

#endif
