#ifndef RUN_MAP_FACTORY_H
#define RUN_MAP_FACTORY_H

// #include "IObject.h"
// using namespace jdb;


#include <vector>
#include <map>
#include <unordered_map>
#include <string>

using namespace std;

class RunMapFactory 
{
public:
	virtual const char * classname( ) const { return "RunMapFactory"; }

	RunMapFactory( bool _dynamic = false ){
		dynamic = _dynamic;
		const vector<int> *runs = &run_15_pp_200_picodst_runs;
		const vector<int> *bruns = &run_15_pp_200_bad_runs;

		if ( false == dynamic ){
			int runIndex = 0;

			for ( int iRun : (*runs) ){
				runIndices[ iRun ] = runIndex;
				indexRuns[ runIndex ] = iRun;
				runIndex++;
			}
			for ( int iRun : (*bruns) ){
				badRuns[ iRun ] = true;
			}

		} else {
			nextRunIndex = -1; // so that it starts at 0
		}
		
	}
	~RunMapFactory(){}

	inline int indexForRun( int rn ){
		if ( dynamic ){
			if ( runIndices.count( rn ) >= 1 )  // already set, then return it
				return runIndices[ rn ];
			else { // push this one onto the map
				nextRunIndex++; // makes it start at 0
				runIndices[ rn ] = nextRunIndex;
				return nextRunIndex;
			}
		}
		return runIndices[ rn ];
	}

	inline int runForIndex( int index ){
		return indexRuns[ index ];
	}

	bool isRunBad( int run ){
		if ( badRuns.count( run ) > 0 && badRuns[ run ] )
			return true;
		return false;
	}

	static int year_run_number( int run ){
		return (run / (int)1e6 ) * (int)1e6;
	}
	static int day( int run ){
		return ( run - year_run_number( run ) ) / 1e3;
	}
	static int run_in_day( int run ){
		return ( run - year_run_number( run ) ) - day( run ) * (int)1e3;
	}


protected:
	map< int, int > runIndices;
	map< int, int > indexRuns;
	unordered_map< int, bool> badRuns;
	static const vector<int> run_15_pp_200_picodst_runs;
	static const vector<int> run_15_pp_200_bad_runs;
	

	bool dynamic;
	int nextRunIndex;

};

#endif