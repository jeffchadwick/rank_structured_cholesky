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
// ThreadSpecificData.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include <rschol/datastructure/ThreadSpecificData.h>

#ifdef USE_OMP
#include <omp.h>
#endif

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
template <class T>
ThreadSpecificData<T>::ThreadSpecificData()
	: _useDefault( false )
{
}

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
template <class T>
ThreadSpecificData<T>::ThreadSpecificData( T defaultValue )
	: _useDefault( true ),
		_defaultValue( defaultValue )
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
template <class T>
ThreadSpecificData<T>::~ThreadSpecificData()
{
}

//////////////////////////////////////////////////////////////////////
// Returns the data associated with the current thread
//////////////////////////////////////////////////////////////////////
template <class T>
T & ThreadSpecificData<T>::get()
{
#ifdef USE_OMP
	int threadNum = omp_get_thread_num();
	int numThreads = omp_get_max_threads();
#else
	int threadNum = 0;
	int numThreads = 1;
#endif

#if 0
	// Enlarge the data array if we need to.  This
	// needs to be done atomically
#pragma omp critical
	{
		if ( numThreads > _data.size() )
		{
			_data.resize( numThreads );
		}
	}

	return _data[ threadNum ];
#endif
#pragma omp critical
	{
		if ( _data.find( threadNum ) == _data.end() )
		{
			if ( _useDefault )
			{
				_data[ threadNum ] = _defaultValue;
			}
			else
			{
				_data[ threadNum ] = T();
			}
		}
	}

	return _data[ threadNum ];
}
