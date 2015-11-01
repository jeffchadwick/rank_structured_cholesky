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

