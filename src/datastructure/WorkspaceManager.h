//////////////////////////////////////////////////////////////////////
// WorkspaceManager.h: Interface for the WorkspaceManager class
//
//////////////////////////////////////////////////////////////////////

#ifndef WORKSPACE_MANAGER_H
#define WORKSPACE_MANAGER_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/VECTOR.h>

#include <util/trace.h>

#include <SETTINGS.h>
#include <TYPES.h>

template <class T>
class Workspace;

//////////////////////////////////////////////////////////////////////
// WorkspaceManager class
//
// 
//////////////////////////////////////////////////////////////////////
template <class T>
class WorkspaceManager {
	public:
		WorkspaceManager( size_t baseSize = BASE_SIZE );

    WorkspaceManager( const WorkspaceManager &manager )
    {
      TRACE_ASSERT( "We probably don't want to be here" );
    }

		// Destructor
		virtual ~WorkspaceManager();

    WorkspaceManager &operator=( const WorkspaceManager &manager )
    {
      TRACE_ASSERT( "We probably don't want to be here" );

      return *this;
    }

	protected:

  private:
    friend class Workspace<T>;

    // Requests a workspace of the given size. _allData will be resized
    // if necessary
    T *requestData( size_t dataSize );

	private:
    static const size_t      BASE_SIZE;

    T                       *_allData;
    size_t                   _dataSize;

};

template <class T>
class Workspace {
  public: 
    // Constructor for vector workspaces
    Workspace( WorkspaceManager<T> &manager, const SizeArray &dataSizes );

    // Constructor for matrix workspaces
    Workspace( WorkspaceManager<T> &manager, const PairArray &dataSizes );

    // Constructor if we only want a single vector workspace
    Workspace( WorkspaceManager<T> &manager, size_t dataSize );

    // Constructor if we only want a single matrix workspace
    Workspace( WorkspaceManager<T> &manager, IndexPair dataSize );

    // Desctructor
    virtual ~Workspace();

    T *workspaceData( int index )
    {
      return _dataPointers[ index ];
    }

    size_t workspaceSize( int index ) const
    {
      return _dataSizes[ index ];
    }

  private:
    // Initializes data pointers based on _dataSizes
    void initDataPointers( WorkspaceManager<T> &manager );

  private:
    SizeArray                _dataSizes;
    std::vector<T *>         _dataPointers;

};

#include <datastructure/WorkspaceManager.cpp>

typedef Workspace<Real>      RealWorkspace;
typedef Workspace<int>       IntWorkspace;

#endif
