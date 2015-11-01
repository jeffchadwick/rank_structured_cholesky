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

// IsoHexElement.h - Definition for the IsoHexElement class

#ifndef ISOHEXELEMENT_H
#define ISOHEXELEMENT_H

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <rschol/mesh/LinearElement.h>

class Mesh;

// IsoHexElement class
// Represents and Isoparametric, 8 node hexahedral element
class IsoHexElement : public LinearElement
{
	public:
		// Constructs a new element given a mesh and
		// the indices of the vertices included in this
		// element.  Also needed are the elastic modulus,
		// Poisson ratio and density.
		IsoHexElement( Mesh *mesh, int index,
									 Real E, Real v, Real density,
									 int v0, int v1, int v2, int v3,
									 int v4, int v5, int v6, int v7,
									 int numGaussPoints );

		virtual ~IsoHexElement()
		{
		}

		// Adds this element's faces to a face map data structure
		virtual void addFaces( FaceMap *faceMap );

		// Given a face number, provide an array of indices that can
		// be used to draw triangles representing the given face
		virtual void getFaceTriangles( int faceNum, IntArray &indices ) const;

		// Returns the strain experienced at the center of the given
		// element face
		virtual VECTOR &getFaceStrain( VECTOR &meshState, int faceNum );

		// Return strain at the center of the element
		virtual VECTOR &getStrain( VECTOR &meshState );

	protected:
		// Virtual pre-computation functions
		virtual void precomputeStiffness();
		virtual void precomputeMass();

	private:
		// Adds the stiffness integrand, evaluated at the given
		// coordinates and scaled by the given weight, to the
		// stiffness matrix
		void addStiffnessIntegrandEvaluation( VEC3F &p, Real weight );

		// Adds the mass integrand, evaluated at the given
		// coordinates and scaled by the given weight, to the
		// mass matrix
		void addMassIntegrandEvaluation( VEC3F &p, Real weight );

		// Evaluates strain at the given xi, eta, zeta coordinates
		void evaluateStrain( VEC3F &p, VECTOR &meshState );

	private:
		Real E, v, density;

		int vertices[8];

		// Store x y and z values explicitly to sync up
		// more easily with generated code.
		Real x0, x1, x2, x3, x4, x5, x6, x7;
		Real y0, y1, y2, y3, y4, y5, y6, y7;
		Real z0, z1, z2, z3, z4, z5, z6, z7;

		// The number of quadrature points to
		// use when evaluating integrals
		int _numGaussPoints;
};

#endif
