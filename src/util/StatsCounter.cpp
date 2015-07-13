//////////////////////////////////////////////////////////////////////
// StatsCounter.cpp: Implementation
//
//////////////////////////////////////////////////////////////////////

#include "StatsCounter.h"

StatsCounter StatsCounter::ourInstance;

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void StatsCounter::writeTimings( const char *filename )
{
  FILE                      *fid = stdout;

  if ( filename != NULL ) {
    fid = fopen( filename, "w" );
  }

  fprintf( fid, "Timings:\n" );
  for ( TimerMap::iterator i = _timers.begin();
        i != _timers.end(); i++ )
  {
    fprintf( fid, "    %60s: %f seconds\n",
            i->first.c_str(), i->second.getTotalSecs() );
  }
  fprintf( fid, "\n" );

  if ( filename != NULL ) {
    fclose( fid );
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void StatsCounter::writePerCycleTimings( const char *filename )
{
  FILE                      *fid = stdout;

  if ( filename != NULL ) {
    fid = fopen( filename, "w" );
  }

  fprintf( fid, "Per-cycle timings:\n" );
  for ( TimerMap::iterator i = _timers.begin();
        i != _timers.end(); i++ )
  {
    fprintf( fid, "    %60s: %f ms\n",
            i->first.c_str(), i->second.getMsPerCycle() );
  }
  fprintf( fid, "\n" );

  if ( filename != NULL ) {
    fclose( fid );
  }
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void StatsCounter::writeFlopCounts( const char *filename )
{
  FILE                      *fid = stdout;

  if ( filename != NULL ) {
    fid = fopen( filename, "w" );
  }

  fprintf( fid, "FLOP counts:\n" );
  for ( SizeMap::iterator i = _flopCounts.begin();
        i != _flopCounts.end(); i++ )
  {
    fprintf( fid, "    %60s: %zd flops\n",
             i->first.c_str(), i->second );
  }
  fprintf( fid, "\n" );

  if ( filename != NULL ) {
    fclose( fid );
  }
}

//////////////////////////////////////////////////////////////////////
// Writes the time per flop ratio for each quantity for
// which this ratio exists
//////////////////////////////////////////////////////////////////////
void StatsCounter::writeRatios( const char *filename )
{
  FILE                      *fid = stdout;

  if ( filename != NULL ) {
    fid = fopen( filename, "w" );
  }

  fprintf( fid, "Timings:\n" );
  for ( TimerMap::iterator i = _timers.begin();
        i != _timers.end(); i++ )
  {
    Real                     ratio = i->second.getTotalSecs();
    size_t                   flopCount = _flopCounts[ i->first ];

    if ( flopCount == 0 )
    {
      continue;
    }

    ratio = ratio / (Real)flopCount;

    fprintf( fid, "    %60s: %e seconds per flop\n", i->first.c_str(), ratio );
  }
  fprintf( fid, "\n" );

  if ( filename != NULL ) {
    fclose( fid );
  }
}

//////////////////////////////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////////////////////////////
StatsCounter::StatsCounter()
{
}

//////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////
StatsCounter::~StatsCounter()
{
  writeTimings( "timing_test.txt" );
  writePerCycleTimings( "cycle_timing_test.txt" );
  writeFlopCounts( "flopcount_test.txt" );
  writeRatios( "ratio_test.txt" );
}

