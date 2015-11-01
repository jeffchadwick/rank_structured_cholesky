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
// ThreadSpecificData.h: Interface for the ThreadSpecificData class
//
//////////////////////////////////////////////////////////////////////

#ifndef THREAD_SPECIFIC_DATA_H
#define THREAD_SPECIFIC_DATA_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <vector>
#include <set>

//////////////////////////////////////////////////////////////////////
// ThreadSpecificData class
//
// A generic container for data which should be replicated
// across threads.
// The template class T must provide a default constructor
// with no arguments
//////////////////////////////////////////////////////////////////////
template <class T>
class ThreadSpecificData {
	public:
		ThreadSpecificData();
		ThreadSpecificData( T defaultValue );

		// Destructor
		virtual ~ThreadSpecificData();

		// Returns the data associated with the current thread
		T &get();

		// Returns all data
		std::map<int, T> &getAll() { return _data; }

	private:
#if 0
		std::vector<T>				_data;
#endif
		std::map<int, T>			_data;

		bool									_useDefault;
		T											_defaultValue;
};

#include <rschol/datastructure/ThreadSpecificData.cpp>

#endif
