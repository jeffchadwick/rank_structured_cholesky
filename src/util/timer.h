#ifndef _TIMER_H
#define _TIMER_H

#include <iostream>
#include <string>

//----------------------------------------
//  Windows counter
//  Uses QueryPerformanceCounter for microsecond accuracy.
//----------------------------------------

#ifdef WIN32

#include <windows.h>
#include <winnt.h>

class Timer
{
	public:

		Timer( const char* name );
		Timer();

		~Timer()
		{
			if( dumpOnDestroy ) {
				dump();
			}
		}

		void tick();
		void tock();
		void pause();
		void unpause();

		//----------------------------------------
		//  Call this between tick/tocks to delimit sections of your code.
		//  Useful for seeing what sections contribute how much to the total cycle time.
		//----------------------------------------
		//void nextSection();

		double getLastCycleSecs();
		double getTotalSecs() const;

		double getMsPerCycle();

		void dump( std::ostream& os );
		void dump() { dump( std::cout ); }

		//----------------------------------------
		//  Returns the number of tick/tock cycles collected.
		//----------------------------------------
		long getNumCycles() { return numEvaluations; }

		void setDumpOnDestroy( bool val ) {
			dumpOnDestroy = val;

		}

	private:

		bool dumpOnDestroy;

		bool inTick;

		const char* _name;

		void initialize();

		LARGE_INTEGER ticksPerSecond;
		LARGE_INTEGER totalTicks;
		long numEvaluations;

		LARGE_INTEGER ti;
		LARGE_INTEGER to;
};

#endif // WIN32

//----------------------------------------
//  Unix/Apple timer
//  Copied from Jernej Barbic's performanceCounter code
//  TODO - errr not done yet heh
//----------------------------------------

#if (defined __unix__) || (defined __APPLE__) || (defined __LINUX__)

#include <stdlib.h>
#include <sys/time.h>

class Timer
{
	public:

		Timer( const char* name ) :
			_totalTime(0.0),
			numEvaluations(0),
			inTick(false),
			_name( name ),
			dumpOnDestroy( false ),
      _hasName(true)
		{
		}

		Timer() :
			_totalTime(0.0),
			numEvaluations(0),
			inTick(false),
			_name("Unnamed Timer"),
      _hasName(false)
		{
		}

		~Timer()
		{
			if( dumpOnDestroy ) {
				dump( std::cout );	
			}
		}

		void setDumpOnDestroy( bool val ) {
			dumpOnDestroy = val;
		}

    void setName( const char *name ) {
      _name = name;
      _hasName = true;
    }

    bool hasName() const {
      return _hasName;
    }

		void tick() {
			if( inTick )
      {
        cerr << "WARNING: tick was called right after another tick"
                " on timer " << _name << endl;;
      }
			inTick = true;

			numEvaluations++;
			gettimeofday(&_startTime,NULL);
		}

		double tock() {
			if( !inTick )
      {
        cerr << "WARNING: tock was called before tick was called"
                " on timer " << _name << endl;;
      }
			inTick = false;

			gettimeofday(&_endTime,NULL);

			double startSeconds = (double)_startTime.tv_sec
                              + (1.e-6)*(double)_startTime.tv_usec;
			double endSeconds = (double)_endTime.tv_sec
                              + (1.e-6)*(double)_endTime.tv_usec;
			_totalTime += endSeconds - startSeconds;
			return endSeconds - startSeconds;
		}

		void pause() {}
		void unpause() {}

		double getTotalSecs() const {
			return _totalTime;
		}

		double getMsPerCycle() {
			return (getTotalSecs()/(double)numEvaluations)*1000;
		}

		void dump( std::ostream& os )
		{
			os << "Timer: " << _name << "\t";
			os << getMsPerCycle() << " ms/cycle\t";
			os << getTotalSecs()*1000 << " total ms\t";
			os << endl;
		}
		void dump() { dump(std::cout); }

    const std::string &name() const
    {
      return _name;
    }

		long getNumCycles() { return numEvaluations; }

	private:

		bool dumpOnDestroy;

		bool inTick;
		std::string _name;
		long numEvaluations;
		long tickSec, tockSec, tickMicroSec, tockMicroSec;

    bool _hasName;

    timeval _startTime;
    timeval _endTime;
    double _totalTime;

		/*
  struct timeval tv;

  gettimeofday(&tv,NULL);

  stopCountSec = tv.tv_sec;
  stopCountMicroSec = tv.tv_usec;


inline double PerformanceCounter::GetElapsedTime()
{
  float elapsedTime = 1.0 * (stopCountSec-startCountSec) + 1E-6 * (stopCountMicroSec - startCountMicroSec);
  return elapsedTime;
		 */
		
};

#endif // (defined __unix__) || (defined __APPLE__)

typedef Timer DirectSolveTimer;

#endif // _TIMER_H
