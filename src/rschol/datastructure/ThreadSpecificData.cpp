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
