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
// StatsCounter.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef STATS_COUNTER_H
#define STATS_COUNTER_H

#include <rschol/linearalgebra/MATRIX.h>
#include <rschol/linearalgebra/MATRIX3.h>
#include <rschol/linearalgebra/VECTOR.h>

#include <rschol/util/timer.h>

#include <rschol/SETTINGS.h>
#include <rschol/TYPES.h>

#include <map>
#include <string>

//////////////////////////////////////////////////////////////////////
// StatsCounter class
//
// Comments
//////////////////////////////////////////////////////////////////////
class StatsCounter {
	public:
#if 0
    typedef std::map<std::string, Timer>   TimerMap;
    typedef std::map<std::string, size_t>  SizeMap;
#endif
    typedef UnorderedMap<std::string, Timer>::type   TimerMap;
    typedef UnorderedMap<std::string, size_t>::type  SizeMap;

    static StatsCounter &get()
    {
      return ourInstance;
    }

    TimerMap &timers()
    {
      return _timers;
    }

    SizeMap &flopCounts()
    {
      return _flopCounts;
    }

    SizeMap &memoryUsage()
    {
      return _memoryUsage;
    }

    void writeTimings( const char *filename = NULL );
    void writePerCycleTimings( const char *filename = NULL );
    void writeFlopCounts( const char *filename = NULL );

    // Writes the time per flop ratio for each quantity for
    // which this ratio exists
    void writeRatios( const char *filename = NULL );

    void reset()
    {
      _timers.clear();
      _flopCounts.clear();
      _memoryUsage.clear();
    }

  private:
		StatsCounter();

		// Destructor
		virtual ~StatsCounter();

	private:
    static StatsCounter      ourInstance;

    TimerMap                 _timers;
    SizeMap                  _flopCounts;
    SizeMap                  _memoryUsage;

};

// Useful macros
#define DO_TIMING 1
#define DO_FLOP_COUNT 1
#undef DO_TIMING
#undef DO_FLOP_COUNT

#ifdef DO_FLOP_COUNT
#define FLOP_COUNT_START { \
  MATRIX::MATRIX_FLOP_COUNT = 0; \
}

#define FLOP_COUNT_END( counterName ) { \
  StatsCounter::get().flopCounts()[ counterName ] += MATRIX::MATRIX_FLOP_COUNT; \
}
#else
#define FLOP_COUNT_START
#define FLOP_COUNT_END( counterName ) { \
}
#endif

#ifdef DO_TIMING
#define TIMING_START( timerName ) { \
  if ( !StatsCounter::get().timers()[ timerName ].hasName() ) { \
    StatsCounter::get().timers()[ timerName ].setName( timerName ); \
  } \
  StatsCounter::get().timers()[ timerName ].tick(); \
}
#define TIMING_STOP( timerName ) { \
  if ( !StatsCounter::get().timers()[ timerName ].hasName() ) { \
    StatsCounter::get().timers()[ timerName ].setName( timerName ); \
  } \
  StatsCounter::get().timers()[ timerName ].tock(); \
}
#else
#define TIMING_START( timerName ) { \
}
#define TIMING_STOP( timerName ) { \
}
#endif

#endif
