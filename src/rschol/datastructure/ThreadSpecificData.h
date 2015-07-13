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