#ifndef _VEC3_H_
#define _VEC3_H_

//////////////////////////////////////////////////////////////////////
// Credit: stolen from libgfx by Michael Garland
//////////////////////////////////////////////////////////////////////

#include <rschol/SETTINGS.h>
#include <iostream>
#include <sstream>
#include <math.h>

template<class T>
class TVEC3 {
private:
  T element[3];

public:
  // Standard constructors
  TVEC3(T s=0) { *this = s; }
  TVEC3(T x, T y, T z) { element[0]=x; element[1]=y; element[2]=z; }

  // Copy constructors & assignment operators
  template<class U> TVEC3(const TVEC3<U>& v) { *this = v; }
  template<class U> TVEC3(const U v[3])
    { element[0]=v[0]; element[1]=v[1]; element[2]=v[2]; }
  template<class U> TVEC3& operator=(const TVEC3<U>& v)
    { element[0]=v[0];  element[1]=v[1];  element[2]=v[2];  return *this; }
  TVEC3& operator=(T s) { element[0]=element[1]=element[2]=s; return *this; }

  // Descriptive interface
  typedef T value_type;
  static int dim() { return 3; }

  // Access methods
  operator       T*()       { return element; }
  operator const T*() const { return element; }

  T& operator[](int i)       { return element[i]; }
  T  operator[](int i) const { return element[i]; }
  operator const T*()        { return element; }

  // Assignment and in-place arithmetic methods
  inline TVEC3& operator+=(const TVEC3& v);
  inline TVEC3& operator-=(const TVEC3& v);
  inline TVEC3& operator*=(T s);
  inline TVEC3& operator/=(T s);

	inline void assign(T x, T y, T z)
	{
		element[0] = x;
		element[1] = y;
		element[2] = z;
	}

  T normalize() {
		T l = norm2(*this);
		if( l!=1.0 && l!=0.0 )
		{
			l = sqrt(l);
			*this /= l;
		}
		return l;
  };
  void clear() {
    T zero = 0.0;
    element[0] = zero;
    element[1] = zero;
    element[2] = zero;
  };

	inline bool hasNaNs() const
	{
		return element[0] != element[0]
			|| element[1] != element[1]
			|| element[2] != element[2];
	}

	void set( const char* s )
	{
		std::istringstream ss( s );
		ss >> *this;
	}

  // Maximum element index
  int maxIndex() const
  {
    return ( element[0] >= element[1] && element[0] >= element[2] ) ? 0 :
         ( ( element[1] >= element[0] && element[1] >= element[2] ) ? 1 : 2 );
  }

  // Minimum element index
  int minIndex() const
  {
    return ( element[0] <= element[1] && element[0] <= element[2] ) ? 0 :
         ( ( element[1] <= element[0] && element[1] <= element[2] ) ? 1 : 2 );
  }

  // Sorts the three indices
  void sortIndices( int indices[3], bool ascend )
  {
    if ( ascend )
    {
      int minIdx = minIndex();

      switch ( minIdx )
      {
        case 0:
        {
          indices[0] = 0;

          indices[1] = ( element[1] <= element[2] ) ? 1 : 2;
          indices[2] = ( element[1] <= element[2] ) ? 2 : 1;

          break;
        }
        case 1:
        {
          indices[0] = 1;

          indices[1] = ( element[0] <= element[2] ) ? 0 : 2;
          indices[2] = ( element[0] <= element[2] ) ? 2 : 0;

          break;
        }
        case 2:
        default:
        {
          indices[0] = 2;

          indices[1] = ( element[0] <= element[1] ) ? 0 : 1;
          indices[2] = ( element[0] <= element[1] ) ? 1 : 0;

          break;
        }
      }
    }
    else
    {
      int maxIdx = minIndex();

      switch ( maxIdx )
      {
        case 0:
        {
          indices[0] = 0;

          indices[1] = ( element[1] >= element[2] ) ? 1 : 2;
          indices[2] = ( element[1] >= element[2] ) ? 2 : 1;

          break;
        }
        case 1:
        {
          indices[0] = 1;

          indices[1] = ( element[0] >= element[2] ) ? 0 : 2;
          indices[2] = ( element[0] >= element[2] ) ? 2 : 0;

          break;
        }
        case 2:
        default:
        {
          indices[0] = 2;

          indices[1] = ( element[0] >= element[1] ) ? 0 : 1;
          indices[2] = ( element[0] >= element[1] ) ? 1 : 0;

          break;
        }
      }
    }
  }
};

////////////////////////////////////////////////////////////////////////
// Method definitions
////////////////////////////////////////////////////////////////////////

template<class T> inline TVEC3<T>& TVEC3<T>::operator+=(const TVEC3<T>& v)
  { element[0] += v[0];   element[1] += v[1];   element[2] += v[2];  return *this; }

template<class T> inline TVEC3<T>& TVEC3<T>::operator-=(const TVEC3<T>& v)
  { element[0] -= v[0];   element[1] -= v[1];   element[2] -= v[2];  return *this; }

template<class T> inline TVEC3<T>& TVEC3<T>::operator*=(T s)
  { element[0] *= s;   element[1] *= s;   element[2] *= s;  return *this; }

template<class T> inline TVEC3<T>& TVEC3<T>::operator/=(T s)
  { element[0] /= s;   element[1] /= s;   element[2] /= s;  return *this; }


////////////////////////////////////////////////////////////////////////
// Operator definitions
////////////////////////////////////////////////////////////////////////

template<class T>
inline TVEC3<T> operator+(const TVEC3<T> &u, const TVEC3<T>& v)
  { return TVEC3<T>(u[0]+v[0], u[1]+v[1], u[2]+v[2]); }

template<class T>
inline TVEC3<T> operator-(const TVEC3<T> &u, const TVEC3<T>& v)
  { return TVEC3<T>(u[0]-v[0], u[1]-v[1], u[2]-v[2]); }

template<class T> inline TVEC3<T> operator-(const TVEC3<T> &v)
  { return TVEC3<T>(-v[0], -v[1], -v[2]); }

  template<class T> inline TVEC3<T> operator*(T s, const TVEC3<T> &v)
  { return TVEC3<T>(v[0]*s, v[1]*s, v[2]*s); }
  template<class T> inline TVEC3<T> operator*(const TVEC3<T> &v, T s)
  { return s*v; }

  template<class T> inline TVEC3<T> operator/(const TVEC3<T> &v, T s)
  { return TVEC3<T>(v[0]/s, v[1]/s, v[2]/s); }

template<class T> inline T operator*(const TVEC3<T> &u, const TVEC3<T>& v)
  { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; }

template<class T> inline TVEC3<T> cross(const TVEC3<T>& u, const TVEC3<T>& v)
  { return TVEC3<T>( u[1]*v[2] - v[1]*u[2], -u[0]*v[2] + v[0]*u[2], u[0]*v[1] - v[0]*u[1] ); }

template<class T>
inline TVEC3<T> operator^(const TVEC3<T>& u, const TVEC3<T>& v)
  { return cross(u, v); }


template<class T>
inline std::ostream &operator<<(std::ostream &out, const TVEC3<T>& v)
  { return out << v[0] << " " << v[1] << " " << v[2]; }

template<class T>
inline std::istream &operator>>(std::istream &in, TVEC3<T>& v)
  { return in >> v[0] >> v[1] >> v[2]; }

////////////////////////////////////////////////////////////////////////
// Misc. function definitions
////////////////////////////////////////////////////////////////////////

template<class T> inline T norm2(const TVEC3<T>& v)  { return v*v; }
template<class T> inline T norm(const TVEC3<T>& v)   { return sqrt(norm2(v)); }

template<class T> inline T infty_norm(const TVEC3<T>& v)   
{
  T result = 0;
  for (int i = 0; i < 3; i++)
  {
    T absVal = (v[i] >= 0) ? v[i] : -v[i];
    result =  (result >= absVal) ? result : absVal;
  }
  return result;
}

////////////////////////////////////////////////////////////////////////
// TODO: wtf is with max/min being taken.  should resolve this issue.
////////////////////////////////////////////////////////////////////////
template<class T> inline TVEC3<T> vec_max (const TVEC3<T>& u, const TVEC3<T>& v)   
{
  TVEC3<T> result;
  for (int i = 0; i < 3; i++)
  {
    result[i] = (u[i] > v[i]) ? u[i] : v[i];
  }
  return result;
}

template<class T> inline TVEC3<T> vec_min (const TVEC3<T>& u, const TVEC3<T>& v)   
{
  TVEC3<T> result;
  for (int i = 0; i < 3; i++)
  {
    result[i] = (u[i] < v[i]) ? u[i] : v[i];
  }
  return result;
}

template<class T> inline void unitize(TVEC3<T>& v)
{
    T l = norm2(v);
    if( l!=1.0 && l!=0.0 )  v /= sqrt(l);
}

template<class T> inline TVEC3<T> project_onto (const TVEC3<T>& vec_to_project, const TVEC3<T> vec_to_project_onto)
{
    return ( vec_to_project_onto * (vec_to_project * vec_to_project_onto / norm2 (vec_to_project_onto) ) );
}

template<class T> inline TVEC3<T> get_orthonormal_vector (TVEC3<T> v)
{
  // make v be a unit vector
  unitize (v);

  // let a be a vector not parallel to v
  TVEC3<T> a (1, 0, 0);
  if (a * v > 0.75)
  {
      a = TVEC3<T> (0, 1, 0);
  }

  // remove from a its projection onto v
  TVEC3<T> aOnV = v * (v * a);
  a -= aOnV;

  // unitize a
  unitize (a);
  return a;
}

typedef TVEC3<Real> VEC3;
typedef TVEC3<Real> VEC3F;

#endif

