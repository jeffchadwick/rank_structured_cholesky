//////////////////////////////////////////////////////////////////////
// WorkspaceManager.cpp: Implementation of the WorkspaceManager class
//
//////////////////////////////////////////////////////////////////////

// Base size of 16 entries
template <class T>
const size_t WorkspaceManager<T>::BASE_SIZE = 16;

//////////////////////////////////////////////////////////////////////
// WorkspaceManager implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
template <class T>
WorkspaceManager<T>::WorkspaceManager( size_t baseSize )
  : _dataSize( baseSize )
{
  TRACE_ASSERT( baseSize > 0, "Invalid base size" );

  _allData = new T[ _dataSize ];
  
  memset( (void *)_allData, 0, _dataSize * sizeof( T ) );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
template <class T>
WorkspaceManager<T>::~WorkspaceManager()
{
  delete[] _allData;
}

//////////////////////////////////////////////////////////////////////
// Requests a workspace of the given size. _allData will be resized
// if necessary
//////////////////////////////////////////////////////////////////////
template <class T>
T *WorkspaceManager<T>::requestData( size_t dataSize )
{
  // Check the data size and reallocate if necessary
  if ( dataSize > _dataSize )
  {
    printf( "Warning: resizing data from %zd ", _dataSize );

    delete[] _allData;

    while ( dataSize > _dataSize )
    {
      _dataSize *= 2;
    }

    printf( "to %zd\n", _dataSize );

    _allData = new T[ _dataSize ];
  }

  return _allData;
}

//////////////////////////////////////////////////////////////////////
// Workspace implementation
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor for vector workspaces
//////////////////////////////////////////////////////////////////////
template <class T>
Workspace<T>::Workspace( WorkspaceManager<T> &manager,
                         const SizeArray &dataSizes )
  : _dataSizes( dataSizes )
{
  initDataPointers( manager );
}

//////////////////////////////////////////////////////////////////////
// Constructor for matrix workspaces
//////////////////////////////////////////////////////////////////////
template <class T>
Workspace<T>::Workspace( WorkspaceManager<T> &manager,
                         const PairArray &dataSizes )
{
  for ( int i = 0; i < dataSizes.size(); i++ )
  {
    _dataSizes.push_back( dataSizes[ i ].first * dataSizes[ i ].second );
  }

  initDataPointers( manager );
}

//////////////////////////////////////////////////////////////////////
// Constructor if we only want a single vector workspace
//////////////////////////////////////////////////////////////////////
template <class T>
Workspace<T>::Workspace( WorkspaceManager<T> &manager, size_t dataSize )
{
  _dataSizes.push_back( dataSize );

  initDataPointers( manager );
}

//////////////////////////////////////////////////////////////////////
// Constructor if we only want a single matrix workspace
//////////////////////////////////////////////////////////////////////
template <class T>
Workspace<T>::Workspace( WorkspaceManager<T> &manager, IndexPair dataSize )
{
  _dataSizes.push_back( dataSize.first * dataSize.second );

  initDataPointers( manager );
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
template <class T>
Workspace<T>::~Workspace()
{
}

//////////////////////////////////////////////////////////////////////
// Initializes data pointers based on _dataSizes
//////////////////////////////////////////////////////////////////////
template <class T>
void Workspace<T>::initDataPointers( WorkspaceManager<T> &manager )
{
  size_t                     totalSize = 0;
  T                         *baseData;

  _dataPointers.clear();
  _dataPointers.resize( _dataSizes.size() );

  for ( int i = 0; i < _dataSizes.size(); i++ )
  {
    totalSize += _dataSizes[ i ];
  }

  baseData = manager.requestData( totalSize );

  for ( int i = 0; i < _dataSizes.size(); i++ )
  {
    _dataPointers[ i ] = baseData;

    baseData += _dataSizes[ i ];
  }
}
