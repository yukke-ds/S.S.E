#include <iostream>
#include <algorithm>
#include <cfloat>
#include <cmath>

#include "search.h"
#include "timeman.h"
#include "usi.h"

TimeManagement Time; // Our global time management object

// init() is called at the beginning of the search and calculates the allowed
// thinking time out of the time control and current game ply. We support four
// different kinds of time controls, passed in 'limits':
//
//  inc == 0 && movestogo == 0 means: x basetime  [sudden death!]
//  inc == 0 && movestogo != 0 means: x moves in y minutes
//  inc >  0 && movestogo == 0 means: x basetime + z increment
//  inc >  0 && movestogo != 0 means: x moves in y minutes + z increment
void TimeManagement::init(Search::LimitsType& limits, Color us, int ply)
{
	int minThinkingTime	= Options["MinimumThinkingTime"];
	int moveOverhead	= Options["MoveOverhead"];
	int slowMover		= Options["SlowMover"];
	int npmsec			= Options["nodestime"];

	// If we have to play in 'nodes as time' mode, then convert from time
	// to nodes, and use resulting values in time management formulas.
	// WARNING: Given npms (nodes per millisecond) must be much lower then
	// the real engine speed to avoid time losses.
	if (npmsec)
	{
		if (!availableNodes) // Only once at game start
			availableNodes = npmsec * limits.time[us]; // Time is in msec

		// Convert from millisecs to nodes
		limits.time[us] = (int)availableNodes;
		limits.inc[us] *= npmsec;
		limits.byoyomi *= npmsec;
		limits.npmsec = npmsec;
	}

	int remainTime;
	if (limits.inc[us] > 0)
		remainTime = limits.time[us] + limits.inc[us] - moveOverhead;
	else
		remainTime = limits.time[us] + limits.byoyomi - moveOverhead;

	optimumTime = maximumTime = std::max(remainTime, minThinkingTime);

	if (Options["Ponder"])
		optimumTime += optimumTime / 4;

	std::cout << "info string MaximumTime:" << maximumTime 
			  << " OptimumTime:" << optimumTime << std::endl;
}