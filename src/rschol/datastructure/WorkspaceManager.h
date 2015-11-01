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
// WorkspaceManager.h: Interface for the WorkspaceManager class
//
//////////////////////////////////////////////////////////////////////

#ifndef WORKSPACE_MANAGER_H
#define WORKSPACE_MANAGER_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/util/trace.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

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

#include <rschol/datastructure/WorkspaceManager.cpp>

typedef Workspace<Real>      RealWorkspace;
typedef Workspace<int>       IntWorkspace;

#endif
