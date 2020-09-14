#ifndef TIMEMAN_H_INCLUDED
#define TIMEMAN_H_INCLUDED

#include "misc.h"
#include "search.h"
#include "thread.h"

// The TimeManagement class computes the optimal time to think depending on
// the maximum available time, the game move number and other parameters.
class TimeManagement {
public:
	void init(Search::LimitsType& limits, Color us, int ply);
	void start() { startTime = startTimeFromPonderhit = now(); }
	void start_ponderhit() { startTimeFromPonderhit = now(); }
	int optimum() const { return optimumTime; }
	int maximum() const { return maximumTime; }
	int elapsed() const { return int(Search::Limits.npmsec ? Threads.nodes_searched() : now() - startTime); }
	int elapsed_from_ponderhit() const { return int(Search::Limits.npmsec ? Threads.nodes_searched() : now() - startTimeFromPonderhit); }

	int64_t availableNodes; // When in 'nodes as time' mode

private:
	TimePoint startTime;
	TimePoint startTimeFromPonderhit;
	int optimumTime;
	int maximumTime;
};

extern TimeManagement Time;

#endif // #ifndef TIMEMAN_H_INCLUDED