#ifndef SEARCH_H_INCLUDED
#define SEARCH_H_INCLUDED

#include <atomic>
#include <vector>

#include "misc.h"
#include "movepick.h"
#include "tt.h"
#include "types.h"

class Position;

namespace Search {

// Threshold used for countermoves based pruning
const int CounterMovePruneThreshold = 0;

// Stack struct keeps track of the information we need to remember from nodes
// shallower and deeper in the tree during the search. Each search thread has
// its own array of Stack objects, indexed by the current ply.
struct Stack {
	Move* pv;
	PieceToHistory* contHistory;
	int ply;
	Move currentMove;
	Move excludedMove;
	Move killers[2];
	Value staticEval;
	int statScore;
	int moveCount;
};

// RootMove struct is used for moves at the root of the tree. For each root move
// we store a score and a PV (really a refutation in the case of moves which
// fail low). Score is normally set at -VALUE_INFINITE for all non-pv moves.
struct RootMove {

	explicit RootMove(Move m) : pv(1, m) {}
	bool extract_ponder_from_tt(Position& pos);
	bool operator==(const Move& m) const { return pv[0] == m; }
	bool operator<(const RootMove& m) const { // Sort in descending order
		return m.score != score ? m.score < score
								: m.previousScore < previousScore;
	}

	Value score = -VALUE_INFINITE;
	Value previousScore = -VALUE_INFINITE;
	int selDepth = 0;
	std::vector<Move> pv;
};

typedef std::vector<RootMove> RootMoves;

// LimitsType struct stores information sent by GUI about available time to
// search the current move, maximum depth/time, if we are in analysis mode or
// if we have to ponder while it's our opponent's turn to move.
struct LimitsType {

	LimitsType() { // Init explicitly due to broken value-initialization of non POD in MSVC
		nodes = time[WHITE] = time[BLACK] = inc[WHITE] = inc[BLACK] = 
		byoyomi = npmsec = movestogo = depth = movetime = mate = infinite = 0;
		silent = false;
	}

	bool use_time_management() const {
		return !(mate | movetime | depth | nodes | infinite);
	}

	std::vector<Move> searchmoves;	
	int time[COLOR_NB], inc[COLOR_NB], byoyomi, npmsec, movestogo, depth,
		movetime, mate, infinite;
	uint64_t nodes;
	bool silent;
};

extern LimitsType Limits;

void init();
void clear();

} // namespace Search

#endif // ifndef SEARCH_H_INCLUDED