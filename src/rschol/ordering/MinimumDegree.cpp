//////////////////////////////////////////////////////////////////////
// MinimumDegree.cpp: Implementation of the MinimumDegree class
//
//////////////////////////////////////////////////////////////////////

#include "MinimumDegree.h"
#include "Ordering.h"

#include <rschol/util/IO.h>

#include <map>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Computes a minimum degree ordering for a block matrix given the
// set of interactions between blocks, and block sizes.
//
// blockInteractions will be modified
//////////////////////////////////////////////////////////////////////
void MinimumDegree::BlockMinimumDegree(
                  const vector<set<int> > &oldInteractions,
                  const IntArray &blockSizes,
                  IntArray &permutation, IntArray &inversePermutation )
{
  // Use a multimap to hack together a priority queue.
  // The first template parameter here is the current "degree" of each node
  DegreeHeap                     degreeHeap;

  vector<set<int> >              blockInteractions( oldInteractions );

  // Track current node degrees
  vector<long int>               blockDegrees;

  int                            numNodes = blockInteractions.size();

  TRACE_ASSERT( numNodes <= blockSizes.size() );

  SymmetrizeInteractionSets( blockInteractions );

  InitializeHeap( blockInteractions, degreeHeap, blockDegrees, blockSizes );

  permutation.resize( numNodes );
  inversePermutation.resize( numNodes );

  // Build the permutation by choosing the best node to eliminate at
  // each iteration
  for ( int block_idx = 0; block_idx < numNodes; block_idx++ )
  {
    if ( block_idx % 1024 == 0 )
    {
      printf( "Completed block %d of %d\n", block_idx + 1, numNodes );
    }

    // Get the cheapest node to eliminate
    DegreeHeap::iterator iter = degreeHeap.begin();

    int                      toRemove = iter->second;

    permutation[ block_idx ] = toRemove;
    inversePermutation[ toRemove ] = block_idx;

    EliminateBlock( blockInteractions, degreeHeap, blockDegrees,
                    blockSizes, toRemove );
  }
  printf( "\n" );
}

//////////////////////////////////////////////////////////////////////
// Prints out usage statistics comparing the number of non-zeros
// in unpermuted and permuted versions of a block system
//////////////////////////////////////////////////////////////////////
void MinimumDegree::PrintBlockOrderingStats(
                  vector<set<int> > &blockInteractions,
                  const IntArray &blockSizes,
                  const IntArray &permutation,
                  const IntArray &inversePermutation )
{
  vector<set<int> >          orderedInteractions( blockInteractions.size() );
  IntArray                   orderedBlockSizes;
  int                        numNodes = blockInteractions.size();
  
  long int                   unorderedNNZ;
  long int                   orderedNNZ;

  // Set up ordered interactions list
  for ( int block_idx = 0; block_idx < numNodes; block_idx++ )
  {
    orderedBlockSizes.push_back( blockSizes.at( permutation.at( block_idx ) ) );

    set<int> &oldInteractions
                  = blockInteractions.at( permutation.at( block_idx ) );

    // Modify indices in these interactions, and remove any upper
    // triangular entries
    set<int> &newInteractions = orderedInteractions.at( block_idx );

    // Copy old interactions which are below the diagonal in
    // the new ordering
    for ( set<int>::iterator iter = oldInteractions.begin();
          iter != oldInteractions.end(); iter++ )
    {
      TRACE_ASSERT( *iter > permutation.at( block_idx ), "Bad index" );

      int newIdx = inversePermutation.at( *iter );

      if ( newIdx > block_idx )
      {
        newInteractions.insert( newIdx );
      }
      else if ( newIdx == block_idx )
      {
        TRACE_ASSERT( NULL, "Should never get here" );
      }
      else
      {
        orderedInteractions.at( newIdx ).insert( block_idx );
      }
    }
  }

  // Figure out fill in behaviour for each ordering
  unorderedNNZ = Ordering::BlockFillIn( blockInteractions, blockSizes );
  orderedNNZ = Ordering::BlockFillIn( orderedInteractions, orderedBlockSizes );

  printf( "\n\n" );
  printf( "Unordered system requires %ld non-zeros (%f MB)\n",
          unorderedNNZ, (Real)( unorderedNNZ * 8 ) / ( 1024.0 * 1024.0 ) );
  printf( "Ordered system requires %ld non-zeros (%f MB)\n",
          orderedNNZ, (Real)( orderedNNZ * 8 ) / ( 1024.0 * 1024.0 ) );
  printf( "\n\n" );
}

//////////////////////////////////////////////////////////////////////
// Given a lower triangular interaction set, compute the symmetric version
// (effectively a graph adjacency matrix)
//////////////////////////////////////////////////////////////////////
void MinimumDegree::SymmetrizeInteractionSets(
                  vector<set<int> > &blockInteractions )
{
  for ( int block_idx = 0; block_idx < blockInteractions.size(); block_idx++ )
  {
    const set<int>          &nodeInteractions = blockInteractions[ block_idx ];

    // Currently this node blockInteractions with a set of nodes of larger
    // index - add "back edges" from these nodes to this one
    for ( set<int>::iterator iter = nodeInteractions.begin();
          iter != nodeInteractions.end(); iter++ )
    {
      // Allows us to consider rectangular blocks for partial minimum
      // degree ordering
      if ( *iter >= blockInteractions.size() )
      {
        continue;
      }

      blockInteractions.at( *iter ).insert( block_idx );
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Initializes the degree heap
//////////////////////////////////////////////////////////////////////
void MinimumDegree::InitializeHeap(
                  vector<std::set<int> > &blockInteractions,
                  DegreeHeap &degreeHeap,
                  vector<long int> &blockDegrees,
                  const IntArray &blockSizes )
{
  for ( int block_idx = 0; block_idx < blockInteractions.size(); block_idx++ )
  {
    if ( block_idx % 1024 == 0 )
    {
      printf( "Insertion: finished block %d of %d\n", block_idx + 1,
              (int)blockInteractions.size() );
    }

    blockDegrees.push_back( BlockDegree( blockInteractions,
                                         blockSizes, block_idx ) );

    HeapInsert( degreeHeap, blockDegrees.back(), block_idx );
  }
  printf( "\n\n" );
}

//////////////////////////////////////////////////////////////////////
// Eliminates a node from the given system and generates the required
// fill-in.  Modifies the priority queue accordingly.
//////////////////////////////////////////////////////////////////////
void MinimumDegree::EliminateBlock(
                  vector<set<int> > &blockInteractions,
                  DegreeHeap &degreeHeap,
                  vector<long int> &blockDegrees,
                  const IntArray &blockSizes,
                  int block_idx )
{
  // Go through each other block this block interacts with
  set<int>                  &interactions = blockInteractions.at( block_idx );

  // Remove the given node from the heap
  HeapDelete( degreeHeap, blockDegrees.at( block_idx ), block_idx );

  for ( set<int>::iterator iter = interactions.begin();
        iter != interactions.end(); iter++ )
  {
    // Allows us to consider rectangular blocks for the purposes
    // of partial minimum degree ordering
    if ( *iter >= blockInteractions.size() )
    {
      continue;
    }

    set<int> &neighbourInteractions = blockInteractions.at( *iter );

    // Remove this node's index from its neighbour
    neighbourInteractions.erase( block_idx );

    // Connect all other neighbours of node block_idx
    for ( set<int>::iterator neighbourIter = interactions.begin();
          neighbourIter != interactions.end(); neighbourIter++ )
    {
      if ( *neighbourIter == *iter )
      {
        continue;
      }

      // Add this in to the current neighbour's interaction set
      // (that is; form a clique amongst all current neighbours)
      neighbourInteractions.insert( *neighbourIter );
    }

    // Remove the neighbour from the heap, update its degree, then
    // re-insert it
    HeapDelete( degreeHeap, blockDegrees.at( *iter ), *iter );

    blockDegrees.at( *iter ) = BlockDegree( blockInteractions, blockSizes,
                                            *iter );

    HeapInsert( degreeHeap, blockDegrees.at( *iter ), *iter );
  }
}

//////////////////////////////////////////////////////////////////////
// Computes a degree for a given block, given its interactions
//////////////////////////////////////////////////////////////////////
long int MinimumDegree::BlockDegree(
                  vector<set<int> > &blockInteractions,
                  const IntArray &blockSizes,
                  int block_idx )
{
  const set<int>            &interactions = blockInteractions.at( block_idx );
  long int                   degree = 0;
  int                        blockSize = blockSizes.at( block_idx );

  for ( set<int>::iterator iter = interactions.begin();
        iter != interactions.end(); iter++ )
  {
    degree += blockSizes.at( *iter );
  }

  // This is the total cost of fill-in introduced by eliminating
  // the given node
  degree *= degree;

  return degree;
}

//////////////////////////////////////////////////////////////////////
// Heap-like data structure operations implemented using a multimap
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Insert
//////////////////////////////////////////////////////////////////////
void MinimumDegree::HeapInsert(
                  DegreeHeap &degreeMap,
                  pair<long int, int> &entry )
{
  pair<DegreeHeap::iterator, DegreeHeap::iterator> rangeIter;

  DegreeHeap::iterator iter;

  rangeIter = degreeMap.equal_range( entry.first );

  for ( iter = rangeIter.first; iter != rangeIter.second; iter++ )
  {
    if ( iter->second == entry.second )
    {
      break;
    }
  }

  TRACE_ASSERT( iter == rangeIter.second, "Duplicate entry found" );

  degreeMap.insert( entry );
}

//////////////////////////////////////////////////////////////////////
// Delete
//////////////////////////////////////////////////////////////////////
void MinimumDegree::HeapDelete(
                  DegreeHeap &degreeMap,
                  pair<long int, int> &entry )
{
  pair<DegreeHeap::iterator, DegreeHeap::iterator> rangeIter;

  DegreeHeap::iterator iter;

  rangeIter = degreeMap.equal_range( entry.first );

  for ( iter = rangeIter.first; iter != rangeIter.second; iter++ )
  {
    if ( iter->second == entry.second )
    {
      degreeMap.erase( iter );
      break;
    }
  }

  TRACE_ASSERT( iter != rangeIter.second, "Entry to remove not found" );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
MinimumDegree::MinimumDegree()
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
MinimumDegree::~MinimumDegree()
{
}

