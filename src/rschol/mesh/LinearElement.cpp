// LinearElement.cpp - Definition for the LinearElement class
// Jeffrey Chadwick (jnc52)

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
