//////////////////////////////////////////////////////////////////////
// ThreadSpecificData.h: Interface for the ThreadSpecificData class
//
//////////////////////////////////////////////////////////////////////

#ifndef THREAD_SPECIFIC_DATA_H
#define THREAD_SPECIFIC_DATA_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/VECTOR.h>

#include <SETTINGS.h>
#include <TYPES.h>

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

#include "ThreadSpecificData.cpp"

#endif
