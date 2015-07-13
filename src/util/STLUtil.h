//////////////////////////////////////////////////////////////////////
// STLUtil.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef STL_UTIL_H
#define STL_UTIL_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/VECTOR.h>

#include <util/trace.h>

#include <SETTINGS.h>
#include <TYPES.h>

#include <map>
#include <string>
#include <vector>

//////////////////////////////////////////////////////////////////////
// Some random utility functions that I needed a place for
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Deletes all contents in a vector of pointers
//////////////////////////////////////////////////////////////////////
template <class T>
void clearVectorContents( std::vector<T *> &data )
{
  for ( int i = 0; i < data.size(); i++ )
  {
    delete data[ i ];
  }
}

//////////////////////////////////////////////////////////////////////
// I/O routines for STL vectors
//////////////////////////////////////////////////////////////////////
template <class T>
void writeVector( const std::vector<T> &data, FILE *f )
{
  int size = data.size();

  fwrite( (void*)&size, sizeof(int), 1, f );

  fwrite( (void*)data.data(), sizeof(T), size, f );
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class T>
void writeVector( const std::vector<T> &data, const char *fileName )
{
  FILE *f = fopen(fileName, "wb");

  if ( f == NULL )
  {
    cerr << "Error opening file " << fileName << endl;
    return;
  }

  writeVector( data, f );

  fclose(f);
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class T>
bool readVector( std::vector<T> &data, FILE *f )
{
  int size;
  size_t bytes_read;

  bytes_read = fread( (void *)&size, sizeof(int), 1, f );

  data.clear();
  data.resize( size );

  bytes_read = fread( (void *)data.data(), sizeof(T), size, f );

  return true;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
template <class T>
bool readVector( std::vector<T> &data, const char *fileName )
{
  FILE *f = fopen(fileName, "rb" );

  if ( f == NULL )
  {
    cerr << "Error opening file " << fileName << endl;
    return false;
  }

  bool result = readVector( data, f );

  fclose(f);

  return result;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
inline int firstNonNegativeEntry( const IntArray &data,
                                  int start_idx = -1, int end_idx = -1 )
{
  int                        entry = -1;

  if ( start_idx < 0 )
  {
    start_idx = 0;
  }
  if ( end_idx < 0 )
  {
    end_idx = data.size() - 1;
  }

  for ( int i = start_idx; i <= end_idx; i++ )
  {
    if ( data[ i ] >= 0 )
    {
      return data[ i ];
    }
  }

  return -1;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
inline int lastNonNegativeEntry( const IntArray &data,
                                 int start_idx = -1, int end_idx = -1 )
{
  if ( start_idx < 0 )
  {
    start_idx = 0;
  }
  if ( end_idx < 0 )
  {
    end_idx = data.size() - 1;
  }

  for ( int i = end_idx; i >= start_idx; i-- )
  {
    if ( data[ i ] >= 0 )
    {
      return data[ i ];
    }
  }

  return -1;
}

//////////////////////////////////////////////////////////////////////
// Array must have size at least 1
//////////////////////////////////////////////////////////////////////
inline int maxEntry( const IntArray &data )
{
  int                        currentMax = data[ 0 ];

  for ( int i = 1; i < data.size(); i++ )
  {
    currentMax = max( currentMax, data[ i ] );
  }

  return currentMax;
}

//////////////////////////////////////////////////////////////////////
// Finds which entries of a sorted integer array lie within
// the given range
//////////////////////////////////////////////////////////////////////
inline void findEntryRangeIntersection( const IntArray &data,
                                        const IndexRange &range,
                                        IndexRange &outputRange )
{
  int                        i = 0;

  outputRange.first = -1;
  outputRange.second = -1;

  if ( data.back() < range.first || data.front() > range.second )
  {
    return;
  }

  while ( i < data.size() && !in_range( range, data[ i ] ) )
  {
    i++;
  }

  if ( i < data.size() )
  {
    outputRange.first = i;
    outputRange.second = i;

    while ( outputRange.second < data.size()
         && in_range( range, data[ outputRange.second ] ) )
    {
      outputRange.second++;
    }

    // Since the previous while loop will "overshoot"
    outputRange.second--;
  }

  return;
}

//////////////////////////////////////////////////////////////////////
// Print memory usage data stored in a map
//////////////////////////////////////////////////////////////////////
inline void printUsageMap( std::map<std::string, long int> &usage,
                           const char *title )
{
  printf( "%s:\n", title );
  for ( std::map<std::string, long int>::iterator i = usage.begin();
        i != usage.end(); i++ )
  {
    printf( "    %40s: %f MB\n", i->first.c_str(),
            (Real)( i->second ) / ( 1024.0 * 1024.0 ) );
  }
  printf( "\n\n" );
}

//////////////////////////////////////////////////////////////////////
// Inverts an integer array
//////////////////////////////////////////////////////////////////////
inline void invertIntArray( const IntArray &input, IntArray &inverse )
{
  if ( inverse.size() != input.size() )
  {
    inverse.resize( input.size() );
  }

  for ( int i = 0; i < input.size(); i++ )
  {
    inverse.at( input.at( i ) ) = i;
  }
}

//////////////////////////////////////////////////////////////////////
// Inverts a partial integer array
//////////////////////////////////////////////////////////////////////
inline void invertPartialIntArray( const IntArray &input, IntArray &inverse,
                                   IndexRange range = IndexRange(-1,-1) )
{
  if ( range.first < 0 || range.second < range.first ) {
    range.first = 0;
    range.second = input.size() - 1;
  }

  for ( int i = range.first; i <= range.second; i++ ) {
    inverse[ input[ i ] ] = i - range.first;
  }
}

inline void clearPartialInverseArray( const IntArray &input, IntArray &inverse,
                                      int clearValue,
                                      IndexRange range = IndexRange(-1,-1) )
{
  if ( range.first < 0 || range.second < range.first ) {
    range.first = 0;
    range.second = input.size() - 1;
  }

  for ( int i = range.first; i <= range.second; i++ ) {
    inverse[ input[ i ] ] = clearValue;
  }
}

#if 0
int findEntryInRange( const IntArray &data, const IndexRange &range,
                      int start_idx, int end_idx )
{
  if ( data[ start_idx ] > range.second || data[ end_idx ] < range.first )
  {
    return -1;
  }

  int                        mid_idx;

  mid_idx = start_idx + ( end_idx - start_idx + 1 ) / 2;

  if ( in_range( range, data[ mid_idx ] ) )
  {
    return mid_idx;
  }
  else if ( mid_idx == start_idx )
  {
    return -1;
  }
  else if ( data[ mid_idx ] < 
}

inline int findEntryInRange( const IntArray &data, const IndexRange &range )
{
}
#endif

#endif
