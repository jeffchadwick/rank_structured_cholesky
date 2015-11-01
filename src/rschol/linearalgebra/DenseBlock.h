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
// DenseBlock.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef DENSE_BLOCK_H
#define DENSE_BLOCK_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

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
