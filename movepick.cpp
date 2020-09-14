#include <cassert>

#include "movepick.h"

namespace {

enum Stages {
	MAIN_TT, CAPTURE_INIT, GOOD_CAPTURE, REFUTATION, QUIET_INIT, QUIET, BAD_CAPTURE,
	EVASION_TT, EVASION_INIT, EVASION,
	PROBCUT_TT, PROBCUT_INIT, PROBCUT,
	QSEARCH_TT, QCAPTURE_INIT, QCAPTURE, QCHECK_INIT, QCHECK
};

const Value LVATable[PIECE_WHITE] = {
	Value(0), Value(1) /*ï‡*/, Value(2)/*çÅ*/, Value(3)/*åj*/, Value(4)/*ã‚*/, Value(7)/*äp*/, Value(8)/*îÚ*/, Value(6)/*ã‡*/,
	Value(10000)/*â§*/, Value(5)/*Ç∆*/, Value(5)/*ê¨çÅ*/, Value(5)/*ê¨åj*/, Value(5)/*ê¨ã‚*/, Value(9)/*în*/, Value(10)/*ó¥*/,Value(11)/*ê¨ã‡*/
};

Value LVA(PieceType pt) { return LVATable[pt]; }

// Helper filter used with select()
const auto Any = []() { return true; };

// partial_insertion_sort() sorts moves in descending order up to and including
// a given limit. The order of moves smaller than the limit is left unspecified.
void partial_insertion_sort(ExtMove* begin, ExtMove* end, int limit) {

	for (ExtMove *sortedEnd = begin, *p = begin + 1; p < end; ++p)
		if (p->value >= limit)
		{
			ExtMove tmp = *p, *q;
			*p = *++sortedEnd;
			for (q = sortedEnd; q != begin && *(q - 1) < tmp; --q)
				*q = *(q - 1);
			*q = tmp;
		}
}

} // namespace

// Constructors of the MovePicker class. As arguments we pass information
// to help it to return the (presumably) good moves first, to decide which
// moves to return (in the quiescence search, for instance, we only want to
// search captures, promotions, and some checks) and how important good move
// ordering is at the current node.

// MovePicker constructor for the main search
MovePicker::MovePicker(const Position& p, Move ttm, Depth d, const ButterflyHistory* mh,
					   const PieceToHistory** ch, Move cm, Move* killers)
	: pos(p), mainHistory(mh), contHistory(ch), 
	  refutations{ { killers[0], 0 },{ killers[1], 0 },{ cm, 0 } }, depth(d) {

	assert(d > DEPTH_ZERO);

	stage = pos.checkers() ? EVASION_TT : MAIN_TT;
	ttMove = ttm && pos.pseudo_legal(ttm) ? ttm : MOVE_NONE;
	stage += (ttMove == MOVE_NONE);
}

// MovePicker constructor for quiescence search
MovePicker::MovePicker(const Position& p, Move ttm, Depth d, const ButterflyHistory* mh, Square rs)
			: pos(p), mainHistory(mh), recaptureSquare(rs), depth(d) {

	assert(d <= DEPTH_ZERO);

	stage = pos.checkers() ? EVASION_TT : QSEARCH_TT;
	ttMove =   ttm
			&& pos.pseudo_legal(ttm)
			&& (depth > DEPTH_QS_RECAPTURES || to_sq(ttm) == recaptureSquare) ? ttm : MOVE_NONE;
	stage += (ttMove == MOVE_NONE);
}

// MovePicker constructor for ProbCut: we generate captures with SEE higher
// than or equal to the given threshold.
MovePicker::MovePicker(const Position& p, Move ttm, Value th)
			: pos(p), threshold(th) {

	assert(!pos.checkers());

	stage = PROBCUT_TT;
	ttMove =   ttm
			&& pos.pseudo_legal(ttm)
			&& pos.capture(ttm)
			&& pos.see_ge(ttm, threshold) ? ttm : MOVE_NONE;
	stage += (ttMove == MOVE_NONE);
}

// score() assigns a numerical value to each move in a list, used for sorting.
// Captures are ordered by Most Valuable Victim (MVV), preferring captures
// with a good history. Quiets are ordered using the histories.
template<GenType Type>
void MovePicker::score() {

	static_assert(Type == CAPTURES || Type == QUIETS || Type == EVASIONS, "Wrong type");

	for (auto& m : *this)
		if (Type == CAPTURES)
		{
			PieceType pt = type_of(pos.piece_on(from_sq(m)));
			m.value = Eval::CapturePieceValue[pos.piece_on(to_sq(m))] - LVA(pt);
		}
		else if (Type == QUIETS)
		{
			Piece movedPiece = pos.moved_piece(m);
			Square movedSq = to_sq(m);

			m.value = (*mainHistory)[pos.side_to_move()][from_to(m)]
					+ (*contHistory[0])[movedPiece][movedSq]
					+ (*contHistory[1])[movedPiece][movedSq]
					+ (*contHistory[3])[movedPiece][movedSq];
		}
		else // Type == EVASIONS
		{
			if (pos.capture(m))
				m.value = Eval::CapturePieceValue[pos.piece_on(to_sq(m))]
						  - LVA(type_of(pos.moved_piece(m)));
			else
				m.value = (*mainHistory)[pos.side_to_move()][from_to(m)] - (1 << 28);
		}
}

// MovePicker::select() returns the next move satisfying a predicate function.
// It never returns the TT move.
template<MovePicker::PickType T, typename Pred>
Move MovePicker::select(Pred filter) {

	while (cur < endMoves)
	{
		if (T == Best)
			std::swap(*cur, *std::max_element(cur, endMoves));

		move = *cur++;

		if (move != ttMove && filter())
			return move;
	}
	return move = MOVE_NONE;
}

// MovePicker::next_move() is the most important method of the MovePicker class. It
// returns a new pseudo legal move every time it is called until there are no more
// moves left, picking the move with the highest score from a list of generated moves.
Move MovePicker::next_move(bool skipQuiets) {

top:
	switch (stage) {

	case MAIN_TT:
	case EVASION_TT:
	case QSEARCH_TT:
	case PROBCUT_TT:
		++stage;
		return ttMove;

	case CAPTURE_INIT:
	case PROBCUT_INIT:
	case QCAPTURE_INIT:
		cur = endBadCaptures = moves;
		endMoves = generate<CAPTURES>(pos, cur);

		score<CAPTURES>();
		++stage;
		goto top;

	case GOOD_CAPTURE:
		if (select<Best>([&]() {
			return pos.see_ge(move, Value(-55 * (cur - 1)->value / 1024)) ?
				// Move losing capture to endBadCaptures to be tried later
				true : (*endBadCaptures++ = move, false); }))
			return move;

			// Prepare the pointers to loop over the refutations array
			cur = std::begin(refutations);
			endMoves = std::end(refutations);

			// If the countermove is the same as a killer, skip it
			if (   refutations[0].move == refutations[2].move
				|| refutations[1].move == refutations[2].move)
				--endMoves;

			++stage;
			/* fallthrough */

	case REFUTATION:
		if (select<Next>([&]() { return    move != MOVE_NONE
			&& !pos.capture_or_pawn_promotion(move)
			&& pos.pseudo_legal(move); }))
			return move;
		++stage;
		/* fallthrough */

	case QUIET_INIT:
		cur = endBadCaptures;
		endMoves = generate<QUIETS>(pos, cur);

		score<QUIETS>();
		partial_insertion_sort(cur, endMoves, -4000 * depth / ONE_PLY);
		++stage;
		/* fallthrough */

	case QUIET:
		if (!skipQuiets
			&& select<Next>([&]() {return   move != refutations[0]
										 && move != refutations[1]
										 && move != refutations[2]; }))
			return move;

		// Prepare the pointers to loop over the bad captures
		cur = moves;
		endMoves = endBadCaptures;

		++stage;
		/* fallthrough */

	case BAD_CAPTURE:
		return select<Next>(Any);

	case EVASION_INIT:
		cur = moves;
		endMoves = generate<EVASIONS>(pos, cur);

		score<EVASIONS>();
		++stage;
		/* fallthrough */

	case EVASION:
		return select<Best>(Any);

	case PROBCUT:
		return select<Best>([&]() { return pos.see_ge(move, threshold); });

	case QCAPTURE:
		if (select<Best>([&]() { return   depth > DEPTH_QS_RECAPTURES
									   || to_sq(move) == recaptureSquare; }))
			return move;

		// If we did not find any move and we do not try checks, we have finished
		if (depth != DEPTH_QS_CHECKS)
			return MOVE_NONE;

		++stage;
		/* fallthrough */

	case QCHECK_INIT:
		cur = moves;
		endMoves = generate_checks<QUIET_CHECKS>(pos, cur);

		++stage;
		/* fallthrough */

	case QCHECK:
		return select<Next>(Any);
	}

	assert(false);
	return MOVE_NONE; // Silence warning
}
