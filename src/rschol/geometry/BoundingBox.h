//////////////////////////////////////////////////////////////////////
// BoundingBox.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>
#include <rschol/linearalgebra/MATRIX.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <math.h>
#include <float.h>

using namespace std;

//////////////////////////////////////////////////////////////////////
// BoundingBox class
//
// Helpful class for storing object bounds
//////////////////////////////////////////////////////////////////////
class BoundingBox {
  public:
    BoundingBox( bool is3D = true )
      : _is3D( is3D )
    {
      _minBound = FLT_MAX;
      _maxBound = -FLT_MAX;
    }

    BoundingBox( const Vector3Array &vertices, bool is3D = true )
      : _is3D( is3D )
    {
      _minBound = FLT_MAX;
      _maxBound = -FLT_MAX;

      for ( int i = 0; i < vertices.size(); i++ )
      {
        if ( vertices.at( i )[0] > xmax() )
        {
          _maxBound[0] = vertices.at( i )[0];
        }
        if ( vertices.at( i )[1] > ymax() )
        {
          _maxBound[1] = vertices.at( i )[1];
        }
        
        if ( _is3D )
        {
          if ( vertices.at( i )[2] > zmax() )
          {
            _maxBound[2] = vertices.at( i )[2];
          }
        }

        if ( vertices.at( i )[0] < xmin() )
        {
          _minBound[0] = vertices.at( i )[0];
        }
        if ( vertices.at( i )[1] < ymin() )
        {
          _minBound[1] = vertices.at( i )[1];
        }

        if ( _is3D )
        {
          if ( vertices.at( i )[2] < zmin() )
          {
            _minBound[2] = vertices.at( i )[2];
          }
        }
      }
    }

    // Helpful for representing bounding box axes.
    enum BoxAxis {
      X_AXIS = 0,
      Y_AXIS,
      Z_AXIS
    };

    BoundingBox( const VEC3F &minBound, const VEC3F &maxBound,
                 bool is3D = true )
      : _minBound( minBound ),
        _maxBound( maxBound ),
        _is3D( is3D )
    {
    }

    void             setMinBound( const VEC3F &minBound )
                     {
                       _minBound = minBound;
                     }
    void             setMaxBound( const VEC3F &maxBound )
                     {
                       _maxBound = maxBound;
                     }

    VEC3F            &minBound()
                     {
                       return _minBound;
                     }
    VEC3F            &maxBound()
                     {
                       return _maxBound;
                     }

    VEC3F             minBound() const
                     {
                       return _minBound;
                     }
    VEC3F             maxBound() const
                     {
                       return _maxBound;
                     }

    Real             xmin() const
                     {
                       return _minBound[0];
                     }
    Real             ymin() const
                     {
                       return _minBound[1];
                     }
    Real             zmin() const
                     {
                       return _minBound[2];
                     }

    Real             axismin( int i ) const
                     {
                       return _minBound[i];
                     }

    Real             xmax() const
                     {
                       return _maxBound[0];
                     }
    Real             ymax() const
                     {
                       return _maxBound[1];
                     }
    Real             zmax() const
                     {
                       return _maxBound[2];
                     }

    Real             axismax( int i ) const
                     {
                       return _maxBound[i];
                     }

    Real             axislength( int i ) const
                     {
                       return axismax( i ) - axismin( i );
                     }

    // Returns the longest axis
    BoxAxis          longestAxis() const
                     {
                       Real maxLength = -1.0;
                       BoxAxis maxAxis;
                       int numAxes = _is3D ? 3 : 2;

                       for ( int i = 0; i < numAxes; i++ )
                       {
                         if ( axislength( i ) > maxLength )
                         {
                           maxAxis = (BoxAxis)i;
                           maxLength = axislength( i );
                         }
                       }

                       return maxAxis;
                     }

    bool             isInside( const VEC3F &x ) const
                     {
                       if ( _is3D )
                       {
                         return ( x[0] >= xmin() &&
                                  x[1] >= ymin() &&
                                  x[2] >= zmin() &&
                                  x[0] <= xmax() &&
                                  x[1] <= ymax() &&
                                  x[2] <= zmax() );
                       }
                       else
                       {
                         return ( x[0] >= xmin() &&
                                  x[1] >= ymin() &&
                                  x[0] <= xmax() &&
                                  x[1] <= ymax() );
                       }
                     }

    // Returns the distance to the boundary of the box.
    // Negative if outside the box - positive if inside.
    Real             boundaryDistance( const VEC3F &x ) const
                     {
                       Real xdist = min( abs( x[0] - xmin() ),
                                          abs( x[0] - xmax() ) );
                       Real ydist = min( abs( x[1] - ymin() ),
                                          abs( x[1] - ymax() ) );
                       Real zdist = min( abs( x[2] - zmin() ),
                                          abs( x[2] - zmax() ) );

                       if ( isInside( x ) )
                         return min( xdist, min( ydist, zdist ) );
                        else
                          return -1.0 * min( xdist, min( ydist, zdist ) );
                     }

    VEC3F            center() const
                     {
                       return (Real)0.5 * ( _minBound + _maxBound );
                     }

    // Whether or not two bounding boxes intersect
    bool             intersects( BoundingBox &b ) const
                     {
                       if ( _is3D )
                       {
                         return ( ( (xmin() <= b.xmax() && xmax() >= b.xmax())
                                 || (b.xmin() <= xmax() && b.xmax() >= xmax()) )
                               && ( (ymin() <= b.ymax() && ymax() >= b.ymax())
                                 || (b.ymin() <= ymax() && b.ymax() >= ymax()) )
                               && ( (zmin() <= b.zmax() && zmax() >= b.zmax())
                                 || (b.zmin() <= zmax() && b.zmax() >= zmax()) )
                                );
                       }
                       else
                       {
                        return ( ( (xmin() <= b.xmax() && xmax() >= b.xmax())
                                 || (b.xmin() <= xmax() && b.xmax() >= xmax()) )
                               && ( (ymin() <= b.ymax() && ymax() >= b.ymax())
                                 || (b.ymin() <= ymax() && b.ymax() >= ymax()) )
                                );
                       }
                     }

    // Projects a point to the boundary of the box, if it lies outside.
    void             projectToBoundary( VEC3F &x ) const
                     {
                       x[0] = min( _maxBound[0], x[0] );
                       x[1] = min( _maxBound[1], x[1] );

                       if ( _is3D )
                       {
                         x[2] = min( _maxBound[2], x[2] );
                       }

                       x[0] = max( _minBound[0], x[0] );
                       x[1] = max( _minBound[1], x[1] );

                       if ( _is3D )
                       {
                         x[2] = max( _minBound[2], x[2] );
                       }
                     }

    // Intersects this bounding box with a vertex, making it
    // larger if necessary
    BoundingBox     &operator+=( const VEC3F &x )
                     {
                       if ( _is3D )
                       {
                         _minBound[0] = min( _minBound[0], x[0] );
                         _minBound[1] = min( _minBound[1], x[1] );
                         _minBound[2] = min( _minBound[2], x[2] );

                         _maxBound[0] = max( _maxBound[0], x[0] );
                         _maxBound[1] = max( _maxBound[1], x[1] );
                         _maxBound[2] = max( _maxBound[2], x[2] );
                       }
                       else
                       {
                         _minBound[0] = min( _minBound[0], x[0] );
                         _minBound[1] = min( _minBound[1], x[1] );

                         _maxBound[0] = max( _maxBound[0], x[0] );
                         _maxBound[1] = max( _maxBound[1], x[1] );
                       }

                       return *this;
                     }

    // Intersects this box with another one
    BoundingBox     &operator+=( const BoundingBox &b )
                     {
                       if ( _is3D )
                       {
                         _minBound[0] = min( _minBound[0], b.xmin() );
                         _minBound[1] = min( _minBound[1], b.ymin() );
                         _minBound[2] = min( _minBound[2], b.zmin() );
   
                         _maxBound[0] = max( _maxBound[0], b.xmax() );
                         _maxBound[1] = max( _maxBound[1], b.ymax() );
                         _maxBound[2] = max( _maxBound[2], b.zmax() );
                       }
                       else
                       {
                         _minBound[0] = min( _minBound[0], b.xmin() );
                         _minBound[1] = min( _minBound[1], b.ymin() );
   
                         _maxBound[0] = max( _maxBound[0], b.xmax() );
                         _maxBound[1] = max( _maxBound[1], b.ymax() );
                       }

                       return *this;
                     }

  private:
    VEC3F              _minBound;
    VEC3F              _maxBound;

    bool               _is3D;
};

// Overloaded operators
inline BoundingBox operator*( BoundingBox &b, Real alpha )
{
  VEC3F minBound = b.center() + alpha * ( b.minBound() - b.center() );
  VEC3F maxBound = b.center() + alpha * ( b.maxBound() - b.center() );

  return BoundingBox( minBound, maxBound );
}

#endif
