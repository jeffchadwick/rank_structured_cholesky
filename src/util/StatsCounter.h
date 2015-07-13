//////////////////////////////////////////////////////////////////////
// StatsCounter.h: Interface
//
//////////////////////////////////////////////////////////////////////

#ifndef STATS_COUNTER_H
#define STATS_COUNTER_H

#include <linearalgebra/MATRIX.h>
#include <linearalgebra/MATRIX3.h>
#include <linearalgebra/VECTOR.h>

#include <util/timer.h>

#include <SETTINGS.h>
#include <TYPES.h>

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
