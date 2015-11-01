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

// LinearElement.h - Header for the LinearElement class

#ifndef LINEARELEMENT_H
#define LINEARELEMENT_H

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/VECTOR.h>
#include <rschol/linearalgebra/VEC3.h>

class FaceMap;

// Element class
// Abstract class used to define basic element routines for a
// linear element.
class LinearElement
{
	public:
		LinearElement( int index )
			: _index(index), _strain(6), _stress(6)
		{
		}

		virtual ~LinearElement()
		{
		}

		// Accessors
		MATRIX &localStiffness() { return _localStiffness; }
		MATRIX &localMass() { return _localMass; }
		IntArray &localMap() { return _localMap; }

		// Precompute necessary data (ie. local stiffness and mass matrices)
		// for this element, given the initial mesh state.
		void precompute();

		// Computes local internal forces based on the given
		// mesh state
		VECTOR& computeLocalInternalForces( VECTOR &meshState );

		// Computes stiffness proportional damping forces for this
		// element, given the current mesh velocity.  Beta is the
		// stiffness damping coefficient.
		VECTOR& computeStiffnessProportionalDamping( VECTOR &meshVelocity,
																								 Real beta );

		int index() { return _index; }

		// Adds this element's faces to a face map data structure
		virtual void addFaces( FaceMap *faceMap ) = 0;

		// Given a face number, provide an array of indices that can
		// be used to draw triangles representing the given face
		virtual void getFaceTriangles( int faceNum, IntArray &indices ) const = 0;

		// Returns the strain experienced at the center of the given
		// element face
		virtual VECTOR &getFaceStrain( VECTOR &meshState, int faceNum ) = 0;

		// Return strain somewhere (ie. probably the center) of the
		// element
		virtual VECTOR &getStrain( VECTOR &meshState ) = 0;

	protected:
		// Virtual pre-computation functions
		virtual void precomputeStiffness() = 0;
		virtual void precomputeMass() = 0;

		// Sets the local state, given an overall mesh state
		void setLocalState( VECTOR &meshState );

	protected:
		// The local stiffness matrix, which can be precomputed
		// due to linearity
		MATRIX _localStiffness;

		// Local mass matrix
		MATRIX _localMass;

		// Local map, which maps element node numbers to
		// global node numbers
		IntArray _localMap;

		// Storage for internal forces computed via this
		// element
		VECTOR _internalForces;

		// Used as temporary storage for the degrees of freedom
		// associated with this particular element
		VECTOR _localState;

		// The global element index
		int _index;

		// Temporary storage for stress and strain
		VECTOR _strain, _stress;
};

#endif
