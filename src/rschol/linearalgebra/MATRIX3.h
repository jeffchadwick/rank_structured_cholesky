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



#ifndef _MAT3_H_
#define _MAT3_H_

#include <rschol/linearalgebra/VEC3.h>

//////////////////////////////////////////////////////////////////////
// 3x3 Matrix class - also stolen from libgfx by Michael Garland
// (I think - TK)
//////////////////////////////////////////////////////////////////////

class MATRIX3
{
private:
  VEC3F row[3];

public:
  // Standard constructors
  MATRIX3() { *this = 0.0; }
  MATRIX3(const VEC3F& r0,const VEC3F& r1,const VEC3F& r2)
    { row[0]=r0; row[1]=r1; row[2]=r2; }
  MATRIX3(const MATRIX3& m) { *this = m; }
  MATRIX3(Real* data);

	// Constructs the standard skew-symmetric cross product matrix
	// for a given vector; that is, the matrix A such that
	// A * x = cross(w, x)
	static MATRIX3 crossProductMatrix( const VEC3F &w );

	// Construct standard rotation matrices
	static MATRIX3 rotationX( Real angle );
	static MATRIX3 rotationY( Real angle );
	static MATRIX3 rotationZ( Real angle );

  // Access methods
  Real& operator()(int i, int j)       { return row[i][j]; }
  Real  operator()(int i, int j) const { return row[i][j]; }
  VEC3F&       operator[](int i)       { return row[i]; }
  const VEC3F& operator[](int i) const { return row[i]; }
  inline VEC3F col(int i) const {return VEC3F(row[0][i],row[1][i],row[2][i]);}

  operator       Real*()       { return row[0]; }
  operator const Real*()       { return row[0]; }
  operator const Real*() const { return row[0]; }

  // Assignment methods
  inline MATRIX3& operator=(const MATRIX3& m);
  inline MATRIX3& operator=(Real s);

  inline MATRIX3& operator+=(const MATRIX3& m);
  inline MATRIX3& operator-=(const MATRIX3& m);
  inline MATRIX3& operator*=(Real s);
  inline MATRIX3& operator/=(Real s);

	// In place multiplication with a 3 vector embedded in some
	// larger vector
	inline void multiplyInPlace(Real *data);

  // Construction of standard matrices
  static MATRIX3 I();
  static MATRIX3 outer_product(const VEC3F& u, const VEC3F& v);
  static MATRIX3 outer_product(const VEC3F& v);

  MATRIX3& diag(Real d);
  MATRIX3& ident() { return diag(1.0); }

  MATRIX3 inverse() {
    MATRIX3 A = adjoint();
    Real d = A[0] * row[0];
    if (d == 0.0f)
      return ident();
    A = A.transpose();
    A /= d;
    return A;
  };
  MATRIX3 adjoint() {
    return MATRIX3(row[1]^row[2], row[2]^row[0], row[0]^row[1]);
  };
  MATRIX3 transpose() {
    return MATRIX3(this->col(0), this->col(1), this->col(2));
  };

  void clear() {
    row[0].clear(); row[1].clear(); row[2].clear();
  };
};

////////////////////////////////////////////////////////////////////////
// Method definitions
////////////////////////////////////////////////////////////////////////

inline MATRIX3& MATRIX3::operator=(const MATRIX3& m)
	{ row[0] = m[0]; row[1] = m[1]; row[2] = m[2];  return *this; }

inline MATRIX3& MATRIX3::operator=(Real s)
	{ row[0]=s;  row[1]=s;  row[2]=s;  return *this; }

inline MATRIX3& MATRIX3::operator+=(const MATRIX3& m)
	{ row[0] += m[0]; row[1] += m[1]; row[2] += m[2]; return *this; }

inline MATRIX3& MATRIX3::operator-=(const MATRIX3& m)
	{ row[0] -= m[0]; row[1] -= m[1]; row[2] -= m[2]; return *this; }

inline MATRIX3& MATRIX3::operator*=(Real s)
	{ row[0] *= s; row[1] *= s; row[2] *= s;  return *this; }

inline MATRIX3& MATRIX3::operator/=(Real s)
	{ row[0] /= s; row[1] /= s; row[2] /= s;  return *this; }

inline void MATRIX3::multiplyInPlace(Real *data)
{
	Real temp0 = row[0][0] * data[0] + row[0][1] * data[1] + row[0][2] * data[2];
	Real temp1 = row[1][0] * data[0] + row[1][1] * data[1] + row[1][2] * data[2];
	data[2] = row[2][0] * data[0] + row[2][1] * data[1] + row[2][2] * data[2];
	data[1] = temp1;
	data[0] = temp0;
}

////////////////////////////////////////////////////////////////////////
// Operator definitions
////////////////////////////////////////////////////////////////////////

inline MATRIX3 operator+(const MATRIX3& n, const MATRIX3& m)
	{ return MATRIX3(n[0]+m[0], n[1]+m[1], n[2]+m[2]); }

inline MATRIX3 operator-(const MATRIX3& n, const MATRIX3& m)
	{ return MATRIX3(n[0]-m[0], n[1]-m[1], n[2]-m[2]); }

inline MATRIX3 operator-(const MATRIX3& m)
	{ return MATRIX3(-m[0], -m[1], -m[2]); }

inline MATRIX3 operator*(Real s, const MATRIX3& m)
	{ return MATRIX3(m[0]*s, m[1]*s, m[2]*s); }
inline MATRIX3 operator*(const MATRIX3& m, Real s)
	{ return s*m; }

inline MATRIX3 operator/(const MATRIX3& m, Real s)
	{ return MATRIX3(m[0]/s, m[1]/s, m[2]/s); }

inline VEC3F operator*(const MATRIX3& m, const VEC3F& v)
	{ return VEC3F(m[0]*v, m[1]*v, m[2]*v); }

extern MATRIX3 operator*(const MATRIX3& n, const MATRIX3& m);

inline std::ostream &operator<<(std::ostream &out, const MATRIX3& M)
	{ return out << "[" << std::endl 
               << M[0] << std::endl << M[1] << std::endl << M[2]
               << std::endl << "]" << std::endl; }

inline std::istream &operator>>(std::istream &in, MATRIX3& M)
	{ return in >> M[0] >> M[1] >> M[2]; }

////////////////////////////////////////////////////////////////////////
// Misc. function definitions
////////////////////////////////////////////////////////////////////////

inline Real det(const MATRIX3& m) { return m[0] * (m[1] ^ m[2]); }

inline Real trace(const MATRIX3& m) { return m(0,0) + m(1,1) + m(2,2); }

#endif

