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

// LinearElement.cpp - Definition for the LinearElement class

#include "LinearElement.h"

// Precompute necessary data for this element
void LinearElement::precompute()
{
	precomputeStiffness();
	precomputeMass();

	_internalForces.resizeAndWipe( _localMap.size() * 3 );
	_localState.resizeAndWipe( _localMap.size() * 3 );
}

// Compute a vector of internal forces acting on the
// nodes of this element.  This can be done via a matrix-vector
// multiply using the local stiffness matrix and the local state
// for this element.
VECTOR& LinearElement::computeLocalInternalForces( VECTOR &meshState )
{
#if 0
	int localIdx = 0;

	// Start by initializing the local state for this element
	for (int i = 0; i < _localMap.size(); i++)
	{
		int globalIdx = _localMap.at(i) * 3;

		for (int j = 0; j < 3; j++)
		{
			_localState(localIdx) = meshState(globalIdx);

			localIdx += 1;
			globalIdx += 1;
		}
	}
#endif
	setLocalState( meshState );

	// Internal forces can now be computed via a simple
	// matrix-vector multiply
	_localStiffness.multiplyInplace( _localState, _internalForces );

	return _internalForces;
}

// Compute a vector of stiffness-proportional damping forces for
// this element.
VECTOR& LinearElement::computeStiffnessProportionalDamping(
																			VECTOR &meshVelocity, Real beta )
{
#if 0
	int localIdx = 0;

	// Start by initializing the local velocities for this element
	for (int i = 0; i < _localMap.size(); i++)
	{
		int globalIdx = _localMap.at(i) * 3;

		for (int j = 0; j < 3; j++)
		{
			_localState(localIdx) = meshVelocity(globalIdx);

			localIdx += 1;
			globalIdx += 1;
		}
	}
#endif
	setLocalState( meshVelocity );

	// Damping forces can now be computed via a simple
	// matrix-vector multiply
	_localStiffness.multiplyInplace( _localState, _internalForces );

	_internalForces *= beta;

	return _internalForces;
}

// Sets the element's local state
void LinearElement::setLocalState( VECTOR &meshState )
{
	int localIdx = 0;

	// Start by initializing the local state for this element
	for (int i = 0; i < _localMap.size(); i++)
	{
		int globalIdx = _localMap.at(i) * 3;

		if ( globalIdx < meshState.size() )
		{
			for (int j = 0; j < 3; j++)
			{
				_localState(localIdx) = meshState(globalIdx);

				localIdx += 1;
				globalIdx += 1;
			}
		}
		else
		{
			// Nodes must be constrained, so displacements are zero
			// TODO - support time-dependent displacements here?
			for (int j = 0; j < 3; j++)
			{
				_localState(localIdx) = 0;

				localIdx += 1;
			}
		}
	}
}
