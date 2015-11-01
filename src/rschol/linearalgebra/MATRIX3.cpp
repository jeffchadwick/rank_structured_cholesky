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



#include "MATRIX3.h"

#include <math.h>

MATRIX3::MATRIX3(Real* data)
{
  for (int y = 0; y < 3; y++)
    for (int x = 0; x < 3; x++)
      (*this)(x,y) = data[x + y * 3];
}

MATRIX3 MATRIX3::crossProductMatrix( const VEC3F &w )
{
	return MATRIX3( VEC3F( 0.0, -w[2], w[1] ),
									VEC3F( w[2], 0.0, -w[0] ),
									VEC3F( -w[1], w[0], 0.0 ) );
}

MATRIX3 MATRIX3::rotationX( Real angle )
{
	MATRIX3 rotation = MATRIX3::I();
	Real cos_theta = cos( angle );
	Real sin_theta = sin( angle );

	rotation(1,1) = cos_theta;
	rotation(2,2) = cos_theta;
	rotation(1,2) = -sin_theta;
	rotation(2,1) = sin_theta;

	return rotation;
}

MATRIX3 MATRIX3::rotationY( Real angle )
{
	MATRIX3 rotation = MATRIX3::I();
	Real cos_theta = cos( angle );
	Real sin_theta = sin( angle );

	rotation(0,0) = cos_theta;
	rotation(2,2) = cos_theta;
	rotation(0,2) = sin_theta;
	rotation(2,0) = -sin_theta;

	return rotation;
}

MATRIX3 MATRIX3::rotationZ( Real angle )
{
	MATRIX3 rotation = MATRIX3::I();
	Real cos_theta = cos( angle );
	Real sin_theta = sin( angle );

	rotation(0,0) = cos_theta;
	rotation(1,1) = cos_theta;
	rotation(0,1) = -sin_theta;
	rotation(1,0) = sin_theta;

	return rotation;
}

MATRIX3 MATRIX3::I() { return MATRIX3(VEC3F(1,0,0), VEC3F(0,1,0), VEC3F(0,0,1)); }

MATRIX3 &MATRIX3::diag(Real d)
{
  *this = 0.0;
  row[0][0] = row[1][1] = row[2][2] = d;
  return *this;
}

MATRIX3 diag(const VEC3F& v)
{
  return MATRIX3(VEC3F(v[0],0,0),  VEC3F(0,v[1],0),  VEC3F(0,0,v[2]));
}

MATRIX3 MATRIX3::outer_product(const VEC3F& v)
{
  MATRIX3 A;
  Real x=v[0], y=v[1], z=v[2];

  A(0,0) = x*x;  A(0,1) = x*y;  A(0,2) = x*z;
  A(1,0)=A(0,1); A(1,1) = y*y;  A(1,2) = y*z;
  A(2,0)=A(0,2); A(2,1)=A(1,2); A(2,2) = z*z;

  return A;
}

MATRIX3 MATRIX3::outer_product(const VEC3F& u, const VEC3F& v)
{
  MATRIX3 A;

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      A(i, j) = u[i]*v[j];

  return A;
}

MATRIX3 operator*(const MATRIX3& n, const MATRIX3& m)
{
  MATRIX3 A;

  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      A(i,j) = n[i]*m.col(j);
  return A;
}

