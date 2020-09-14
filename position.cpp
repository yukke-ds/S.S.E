#include <algorithm>
#include <cassert>
#include <cstring> // For std::memset, std::memcmp
#include <iomanip>
#include <sstream>
#include <fstream>

#include "misc.h"
#include "movegen.h"
#include "position.h"
#include "thread.h"
#include "tt.h"
#include "usi.h"

using std::string;

namespace Zobrist {

Key psq[COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];
Key hand[COLOR_NB][PIECE_HAND_NB];
Key side;

} // namespace Zobrist

namespace {

const string PieceToChar(" PLNSBRGK        plnsbrgk");

// min_attacker() is a helper function used by see() to locate the least
// valuable attacker for the side to move, remove the attacker we just found
// from the bitboards and scan for new X-ray attacks behind it.
template<PieceType Pt>
PieceType min_attacker(const Position& pos, Color stm, Square to, Bitboard stmAttackers,
	Bitboard& occupied, Bitboard& attackers) {

	Bitboard b = stmAttackers & pos.pieces(Pt);
	if (!b)
		return min_attacker<PieceType(Pt + 1)>(pos, stm, to, stmAttackers, occupied, attackers);

	Square from = b.pop();
	occupied ^= from;

	if (Pt == PAWN || Pt == LANCE)
		attackers |= (lance_attacks(~stm, to, occupied) & (pos.pieces(ROOK, DRAGON) | pos.pieces(stm, LANCE)));

	if (Pt == GOLD || Pt == PRO_PAWN || PRO_LANCE || Pt == PRO_KNIGHT || Pt == PRO_SILVER || Pt == HORSE || Pt == DRAGON)
		attackers |= (lance_attacks(~stm, to, occupied) & pos.pieces(stm, LANCE))
				   | (lance_attacks(stm, to, occupied) & pos.pieces(~stm, LANCE))
				   | (rook_attacks(to, occupied) & pos.pieces(ROOK, DRAGON))
				   | (bishop_attacks(to, occupied) & pos.pieces(BISHOP, HORSE));

	if (Pt == SILVER)
		attackers |= (lance_attacks(~stm, to, occupied) & pos.pieces(stm, LANCE))
				   | (rook_attacks(to, occupied) & pos.pieces(ROOK, DRAGON))
				   | (bishop_attacks(to, occupied) & pos.pieces(BISHOP, HORSE));

	if (Pt == BISHOP)
		attackers |= (bishop_attacks(to, occupied) & pos.pieces(BISHOP, HORSE));

	if (Pt == ROOK)
		attackers |= (lance_attacks(~stm, to, occupied) & pos.pieces(stm, LANCE))
				   | (lance_attacks(stm, to, occupied) & pos.pieces(~stm, LANCE))
				   | (rook_attacks(to, occupied) & pos.pieces(ROOK, DRAGON));

	attackers &= occupied; // After X-ray that may add already processed pieces

	if ( !(Pt & PIECE_PROMOTION) 
		&& Pt != GOLD
		&& (is_promotion_area(stm, to) || is_promotion_area(stm, from)))
		return promoted_type(Pt);
	else
		return Pt;
}

template<>
PieceType min_attacker<KING>(const Position&, Color, Square, Bitboard, Bitboard&, Bitboard&) {
	return KING; // No need to update bitboards: it is the last cycle
}

} // namespace

// Position::init()は、ハッシュキー関連の初期化を行う
void Position::init() {

	PRNG rng(1070372);

	for (Color c = BLACK; c <= WHITE; ++c)
		for (PieceType pt = PAWN; pt <= DRAGON; ++pt)
			for (Square s = SQ_11; s <= SQ_99; ++s)
				Zobrist::psq[c][pt][s] = rng.rand<Key>() & ~1ULL;

	for (Color c = BLACK; c <= WHITE; ++c)
		for (PieceType pt = PAWN; pt <= GOLD; ++pt)
			Zobrist::hand[c][pt] = rng.rand<Key>() & ~1ULL;

	Zobrist::side = 1ULL;
}

// Position::set()は、与えられたFEN文字列と共にpositionオブジェクトを初期化する。
// この関数はとても脆弱で、GUIからFEN文字列が正しく入力されることを前提としている。
Position& Position::set(const string& sfenStr, StateInfo* si, Thread* th) {

	unsigned char token;
	bool promotion = false;
	size_t idx;
	File f = FILE_9;
	Rank r = RANK_1;
	std::istringstream ss(sfenStr);

	std::memset(this, 0, sizeof(Position));
	std::memset(si, 0, sizeof(StateInfo));
	st = si;

	ss >> std::noskipws;

	// 1. EvalListの初期化
	// PieceListを更新する上で、どの駒がどこにあるかを設定しなければならないが、
	// それぞれの駒をどこまで使ったかのカウンタ
	PieceNumber piece_no_count[PIECE_HAND_NB] = {
		PIECE_NUMBER_ZERO, PAWN_NUMBER, LANCE_NUMBER, KNIGHT_NUMBER, 
		SILVER_NUMBER, BISHOP_NUMBER, ROOK_NUMBER, GOLD_NUMBER
	};

	evalList.clear();

	for (PieceNumber pn = PIECE_NUMBER_ZERO; pn < PIECE_NUMBER_NB; ++pn)
		evalList.put_board_piece(pn, NO_PIECE, SQUARE_ZERO);

	// 2. 駒の配置
	while ((ss >> token) && !isspace(token))
	{
		if (isdigit(token))
			f -= File(token - '0');

		else if (token == '/')
		{
			f = FILE_9;
			++r;
		}

		else if (token == '+')
			promotion = true;

		else if ((idx = PieceToChar.find(token)) != string::npos)
		{
			Piece pc = Piece(idx + (promotion ? PIECE_PROMOTION : 0));
			Square sq = f | r;
			PieceNumber pn = (idx == B_KING) ? BKING_NUMBER :
						 (idx == W_KING) ? WKING_NUMBER :
						 piece_no_count[hand_type(type_of(Piece(idx)))]++;

			put_piece(color_of(Piece(idx)), type_of(pc), sq);
			evalList.put_board_piece(pn, pc, sq);

			promotion = false;
			--f;
		}
	}

	// 3. 手番
	ss >> token;
	sideToMove = (token == 'b' ? BLACK : WHITE);
	ss >> token;

	// 4. 持ち駒
	int counts = 0;
	while ((ss >> token) && !isspace(token))
	{
		if (token == '-')
			break;

		if (isdigit(token))
			counts = (token - '0') + counts * 10;

		else if ((idx = PieceToChar.find(token)) != string::npos)
		{
			counts = std::max(counts, 1);
			put_hand(hand[color_of(Piece(idx))], type_of(Piece(idx)), counts);

			for (int i = 0; i < counts; ++i)
			{
				PieceType pt = type_of(Piece(idx));
				PieceNumber pn = piece_no_count[pt]++;
				assert(is_ok(pn));
				evalList.put_hand_piece(pn, pt, color_of(Piece(idx)), i);
			}

			counts = 0;
		}
	}

	// 5. 手数
	ss >> std::skipws >> gamePly;

	thisThread = th;
	set_state(st);

	assert(pos_is_ok());

	return *this;
}

// 局面のsfen文字列を取得する (Position::set()の逆変換)
const std::string Position::sfen() const {

	std::ostringstream ss;

	// 1. 盤面
	int emptyCnt;
	for (Rank r = RANK_1; r <= RANK_9; ++r)
	{
		for (File f = FILE_9; f >= FILE_1; --f)
		{
			// 駒のない升の数は数字で表現するので，それを数える
			for (emptyCnt = 0; f >= FILE_1 && piece_on(f | r) == NO_PIECE; --f)
				++emptyCnt;

			// 駒のない升の数を出力
			if (emptyCnt)
				ss << emptyCnt;

			// 駒があったなら，それに対応する駒文字を出力
			if (f >= FILE_1)
			{
				Piece pc = piece_on(f | r);
				if ((pc & PIECE_PROMOTION) && type_of(pc) != KING) {
					ss << '+';
					ss << PieceToChar[pc - PIECE_PROMOTION];
				}
				else
					ss << PieceToChar[pc];
			}
		}

		// 最下段以外では次の行があるのでセパレータである'/'を出力する
		if (r < RANK_9)
			ss << '/';
	}

	// 2. 手番
	ss << (sideToMove == BLACK ? " b " : " w ");

	// 3. 持ち駒
	int handCnt;
	bool found = false;
	for (Color c = BLACK; c <= WHITE; ++c) {
		for (int pn = 0; pn < 7; ++pn)
		{
			// 飛，角，金，銀，桂，香，歩の順で入れる
			const PieceType USI_Hand[7] = { ROOK, BISHOP, GOLD, SILVER, KNIGHT, LANCE, PAWN };

			PieceType pt = USI_Hand[pn];

			// その種類の手駒の枚数
			handCnt = count_of(hand[c], pt);

			// その種類の手駒を持っているか
			if (handCnt != 0)
			{
				// 手駒が1枚でも見つかった
				found = true;

				// その種類の駒の枚数．1ならば出力を省略
				if (handCnt != 1)
					ss << handCnt;

				ss << PieceToChar[make_piece(c, pt)];
			}
		}
	}

	// 手駒がない場合はハイフンを出力
	ss << (found ? " " : "- ");

	// 4. 手数
	ss << gamePly;

	return ss.str();
}

// Position::set_check_info() sets king attacks to detect if a move gives check
void Position::set_check_info(StateInfo* si) const {

	si->blockersForKing[BLACK] = slider_blockers(pieces(WHITE), king_square(BLACK), si->pinnersForKing[BLACK]);
	si->blockersForKing[WHITE] = slider_blockers(pieces(BLACK), king_square(WHITE), si->pinnersForKing[WHITE]);

	Color them = ~sideToMove;		// 相手の手番
	Square ksq = king_square(them);	// 敵玉の位置
	Bitboard occupied = pieces();	// 盤上の駒の位置が1のビットボード

	// それぞれの駒によって、相手玉に王手できるマスが1のビットボード
	si->checkSquares[PAWN] = pawn_attacks(them, ksq);
	si->checkSquares[LANCE] = lance_attacks(them, ksq, occupied);
	si->checkSquares[KNIGHT] = knight_attacks(them, ksq);
	si->checkSquares[SILVER] = silver_attacks(them, ksq);
	si->checkSquares[BISHOP] = bishop_attacks(ksq, occupied);
	si->checkSquares[ROOK] = rook_attacks(ksq, occupied);
	si->checkSquares[GOLD] = gold_attacks(them, ksq);
	si->checkSquares[KING] = ALL0BB;
	si->checkSquares[PRO_PAWN] = si->checkSquares[GOLD];
	si->checkSquares[PRO_LANCE] = si->checkSquares[GOLD];
	si->checkSquares[PRO_KNIGHT] = si->checkSquares[GOLD];
	si->checkSquares[PRO_SILVER] = si->checkSquares[GOLD];
	si->checkSquares[HORSE] = si->checkSquares[BISHOP] | king_attacks(ksq);
	si->checkSquares[DRAGON] = si->checkSquares[ROOK] | king_attacks(ksq);
}

// Position::set_state() computes the hash keys of the position, and other
// data that once computed is updated incrementally as moves are made.
// The function is only used when a new position is set up, and to verify
// the correctness of the StateInfo data when running in debug mode.
void Position::set_state(StateInfo* si) const {

	si->checkersBB = attackers_to(king_square(sideToMove)) & pieces(~sideToMove);

	si->boardKey = sideToMove == BLACK ? 0ULL : Zobrist::side;
	si->handKey = 0ULL;

	set_check_info(si);

	// 盤面のハッシュキーを計算
	for (Bitboard b = pieces(); b; )
	{
		Square s = b.pop();
		Piece pc = piece_on(s);
		si->boardKey += Zobrist::psq[color_of(pc)][type_of(pc)][s];
	}

	// 持ち駒のハッシュキーを計算
	for (Color c = BLACK; c <= WHITE; ++c)
		for (PieceType pt = PAWN; pt <= GOLD; ++pt)
			si->handKey += Zobrist::hand[c][pt] * (int64_t)count_of(hand[c], pt);

	si->hand = hand[sideToMove];
	si->material = Eval::material(*this);
	Eval::all_calculate(*this);
}

// Position::slider_blockers() returns a bitboard of all the pieces (both colors)
// that are blocking attacks on the square 's' from 'sliders'. A piece blocks a
// slider if removing that piece from the board would result in a position where
// square 's' is attacked. For example, a king-attack blocking piece can be either
// a pinned or a discovered check piece, according if its color is the opposite
// or the same of the color of the slider.
Bitboard Position::slider_blockers(Bitboard sliders, Square s, Bitboard& pinners) const {

	Color us = color_of(piece_on(s));
	Bitboard result = ALL0BB;
	pinners = ALL0BB;

	// Snipers are sliders that attack 's' when a piece is removed
	Bitboard snipers = ((pieces(ROOK, DRAGON) & rook_attacks(s, ALL0BB))
		| (pieces(BISHOP, HORSE) & bishop_attacks(s, ALL0BB))
		| (pieces(LANCE) & lance_attacks(us, s, ALL0BB))) & sliders;

	while (snipers)
	{
		// ksqとpinnerに挟まれている駒を列挙
		Square sniperSq = snipers.pop();
		Bitboard blocker = between_bb(s, sniperSq) & pieces();

		// ksqとpinnerの間に2駒以上なければ開き王手可能
		if (!more_than_one(blocker))
		{
			result |= blocker;
			if (blocker & pieces(us))
				pinners |= sniperSq;
		}
	}

	return result;
}

// Position::attackers_to() computes a bitboard of all pieces which attack a
// given square. Slider attacks use the occupied bitboard to indicate occupancy.
Bitboard Position::attackers_to(Square s, Bitboard occupied) const {

	return (pawn_attacks(BLACK, s) & pieces(WHITE, PAWN))
		| (pawn_attacks(WHITE, s) & pieces(BLACK, PAWN))
		| (lance_attacks(BLACK, s, occupied) & pieces(WHITE, LANCE))
		| (lance_attacks(WHITE, s, occupied) & pieces(BLACK, LANCE))
		| (knight_attacks(BLACK, s) & pieces(WHITE, KNIGHT))
		| (knight_attacks(WHITE, s) & pieces(BLACK, KNIGHT))
		| (silver_attacks(BLACK, s) & pieces(WHITE, SILVER))
		| (silver_attacks(WHITE, s) & pieces(BLACK, SILVER))
		| (gold_attacks(BLACK, s) & gold_pieces(WHITE))
		| (gold_attacks(WHITE, s) & gold_pieces(BLACK))
		| (bishop_attacks(s, occupied) & pieces(BISHOP, HORSE))
		| (rook_attacks(s, occupied) & pieces(ROOK, DRAGON))
		| (king_attacks(s) & (pieces(KING) | pieces(DRAGON, HORSE)));
}

// Position::legal() tests whether a pseudo-legal move is legal
bool Position::legal(Move m) const {

	assert(is_ok(m));

	Color us = sideToMove;
	Square from = from_sq(m);

	assert(piece_on(king_square(us)) == make_piece(us, KING));

	if (type_of(m) == DROP)
		return true;

	assert(color_of(moved_piece(m)) == us);

	// If the moving piece is a king, check whether the destination
	// square is attacked by the opponent.
	if (type_of(piece_on(from)) == KING)
		return !(attackers_to(to_sq(m)) & pieces(~us));

	// A non-king move is legal if and only if it is not pinned or it
	// is moving along the ray towards or away from the king.
	return !(pinned_pieces(us) & from)
		|| aligned(from, to_sq(m), king_square(us));
}

// Position::pseudo_legal() takes a random move and tests whether the move is
// pseudo legal. It is used to validate moves from TT that can be corrupted
// due to SMP concurrent access or hash position key aliasing.
bool Position::pseudo_legal(const Move m) const {

	Color us = sideToMove;
	Square to = to_sq(m);

	if (type_of(m) == DROP)
	{
		PieceType pt = drop_type(m);
		assert(pt >= PAWN && pt < KING);

		if (piece_on(to) != NO_PIECE || count_of(hand[us], pt) == 0)
			return false;

		if (checkers())
		{
			Bitboard target = checkers();
			Square checksq = target.pop();

			// double checkers
			if (target)
				return false;

			if (!(between_bb(checksq, king_square(us)) & to))
				return false;
		}

		if (pt == PAWN && !legal_drop_pawn(to))
			return false;
	}
	else {
		Square from = from_sq(m);
		Piece pc = piece_on(from);

		// If the 'from' square is not occupied by a piece belonging to the side to
		// move, the move is obviously not legal.
		if (pc == NO_PIECE || color_of(pc) != us)
			return false;

		if (!(attacks_from(pc, from) & to))
			return false;

		// The destination square cannot be occupied by a friendly piece
		if (pieces(us) & to)
			return false;

		PieceType pt = type_of(pc);
		if (type_of(m) == PROMOTION)
		{
			if (pt >= GOLD)
				return false;

			if (!(promotion_area(us) & (Bitboard(from) | Bitboard(to))))
				return false;
		}
		else {
			if (pt == PAWN || pt == LANCE)
				if (relative_rank(us, rank_of(to)) == RANK_1)
					return false;
		}

		// Evasions generator already takes care to avoid some kind of illegal moves
		// and legal() relies on this. We therefore have to take care that the same
		// kind of moves are filtered out here.
		if (checkers())
		{
			if (pt != KING)
			{
				// Double check? In this case a king move is required
				if (more_than_one(checkers()))
					return false;

				// Our move must be a blocking evasion or a capture of the checking piece
				if (!((between_bb(checkers().pop(), king_square(us)) | checkers()) & to))
					return false;
			}

			// In case of king moves under check we have to remove king so as to catch
			// invalid moves like b1a1 when opposite queen is on c1.
			else if (attackers_to(to, pieces() ^ from) & pieces(~us))
				return false;
		}
	}

	return true;
}

// Position::gives_check() tests whether a pseudo-legal move gives a check
bool Position::gives_check(Move m) const {

	assert(is_ok(m));

	Square to = to_sq(m);

	if (type_of(m) == DROP)
		return st->checkSquares[drop_type(m)] & to;

	else
	{
		Square from = from_sq(m);
		PieceType ptFrom = type_of(piece_on(from));
		PieceType ptTo = type_of(m) == PROMOTION ? promoted_type(ptFrom) : ptFrom;

		// Is there a direct check?
		if (st->checkSquares[ptTo] & to)
			return true;

		// Is there a discovered check?
		if ((discovered_check_candidates() & from)
			&& !(LineBB[from][to] & king_square(~sideToMove)))
			return true;
	}

	return false;
}

// Position::do_move() makes a move, and saves all information necessary
// to a StateInfo object. The move is assumed to be legal. Pseudo-legal
// moves should be filtered out before this function is called.
void Position::do_move(Move m, StateInfo& newSt, bool givesCheck) {

	assert(is_ok(m));
	assert(&newSt != st);

	thisThread->nodes.fetch_add(1, std::memory_order_relaxed);
	Key k = st->boardKey ^ Zobrist::side;
	Key h = st->handKey;

	// Copy some fields of the old state to our new StateInfo object except the
	// ones which are going to be recalculated from scratch anyway and then switch
	// our state pointer to point to the new (ready to be updated) state.
	std::memcpy(&newSt, st, offsetof(StateInfo, boardKey));
	newSt.previous = st;
	st = &newSt;

	// Increment ply counters.
	++gamePly;
	++st->pliesFromNull;

	st->eval.p[0][0] = VALUE_NOT_EVALUATED;

	Color us = sideToMove;
	Color them = ~us;
	Square to = to_sq(m);
	PieceType captured = type_of(piece_on(to));
	Value materialDiff;
	MovedPieceInfo& mpi = st->mpi;
	mpi.numPieceMoved = 1;

	if (type_of(m) == DROP)
	{
		PieceType pt = drop_type(m);
		PieceNumber pn = piece_no_of(us, pt);

		assert(piece_on(to) == NO_PIECE);
		assert(is_ok(pn));

		sub_hand(hand[us], pt);
		put_piece(us, pt, to);

		// Update about diff evaluation(drop)
		mpi.movedPieceNumber[0] = pn;
		mpi.changedPiece[0].oldPiece = evalList.eval_index(pn);
		evalList.put_board_piece(pn, make_piece(us, pt), to);
		mpi.changedPiece[0].newPiece = evalList.eval_index(pn);

		k += Zobrist::psq[us][pt][to];
		h -= Zobrist::hand[us][pt];

		materialDiff = VALUE_ZERO;
		captured = NO_PIECE_TYPE;

		if (givesCheck) {
			st->checkersBB = Bitboard(to);
			st->continuousCheck[us] += 2;
		}
		else {
			st->checkersBB = ALL0BB;
			st->continuousCheck[us] = 0;
		}
	}
	else
	{
		Square from = from_sq(m);
		PieceType ptFrom = type_of(piece_on(from));
		PieceType ptTo;
		
		if (type_of(m) == PROMOTION) {
			ptTo = promoted_type(ptFrom);
			materialDiff = Eval::PromotionDiff[piece_on(from)];
		}
		else {
			ptTo = ptFrom;
			materialDiff = VALUE_ZERO;
		}

		if (captured)
		{
			PieceType ptHand = hand_type(captured);
			PieceNumber pnCaptured = piece_no_of(to);

			assert(is_ok(pnCaptured));

			materialDiff += Eval::CapturePieceValue[piece_on(to)];

			// Update about diff evaluation(captured)
			mpi.numPieceMoved = 2;
			mpi.movedPieceNumber[1] = pnCaptured;
			mpi.changedPiece[1].oldPiece = evalList.eval_index(pnCaptured);
			evalList.put_hand_piece(pnCaptured, ptHand, us, count_of(hand[us], ptHand));
			mpi.changedPiece[1].newPiece = evalList.eval_index(pnCaptured);

			// Update board
			remove_piece(them, captured, to);
			add_hand(hand[us], ptHand);

			// Update board hash key and hand hash key
			k -= Zobrist::psq[them][captured][to];
			h += Zobrist::hand[us][ptHand];
		}

		PieceNumber pnFrom = piece_no_of(from);

		assert(is_ok(pnFrom));

		// Update board
		remove_piece(us, ptFrom, from);
		put_piece(us, ptTo, to);

		// Update about diff evaluation(from-to)
		mpi.movedPieceNumber[0] = pnFrom;
		mpi.changedPiece[0].oldPiece = evalList.eval_index(pnFrom);
		evalList.put_board_piece(pnFrom, make_piece(us, ptTo), to);
		mpi.changedPiece[0].newPiece = evalList.eval_index(pnFrom);

		// Update board hash key
		k += Zobrist::psq[us][ptTo][to] - Zobrist::psq[us][ptFrom][from];

		// Calculate checkers bitboard (if move gives check)
		if (givesCheck) {
			st->checkersBB = attackers_to(kingSquare[them]) & pieces(us);
			st->continuousCheck[us] += 2;
		}
		else {
			st->checkersBB = ALL0BB;
			st->continuousCheck[us] = 0;
		}
	}

	st->material = st->previous->material + (us == BLACK ? materialDiff : -materialDiff);

	// Set capture piece
	st->capturedType = captured;

	// Update the keys with the final value
	st->boardKey = k;
	st->handKey = h;

	sideToMove = ~sideToMove;

	// Set hand
	st->hand = hand[sideToMove];

	// Update king attacks used for fast check detection
	set_check_info(st);

	assert(pos_is_ok());
}

// Position::undo_move() unmakes a move. When it returns, the position should
// be restored to exactly the same state as before the move was made.
void Position::undo_move(Move m) {

	assert(is_ok(m));

	sideToMove = ~sideToMove;

	Color us = sideToMove;
	Square to = to_sq(m);
	PieceNumber pnFrom = piece_no_of(to);

	assert(is_ok(pnFrom));

	if (type_of(m) == DROP)
	{
		PieceType pt = drop_type(m);

		evalList.put_hand_piece(pnFrom, pt, us, count_of(hand[us], pt));

		remove_piece(us, pt, to);
		add_hand(hand[us], pt);
	}
	else {
		Square from = from_sq(m);
		PieceType ptTo = type_of(piece_on(to));
		PieceType ptFrom = type_of(m) == PROMOTION ? hand_type(ptTo) : ptTo;

		remove_piece(us, ptTo, to);
		put_piece(us, ptFrom, from);

		evalList.put_board_piece(pnFrom, make_piece(us, ptFrom), from);

		if (st->capturedType)
		{		
			PieceNumber pnCaptured = piece_no_of(us, hand_type(st->capturedType));

			sub_hand(hand[us], hand_type(st->capturedType));
			put_piece(~us, st->capturedType, to); // Restore the captured piece

			assert(is_ok(pnCaptured));
			evalList.put_board_piece(pnCaptured, make_piece(~us, st->capturedType), to);
		}
	}

	// Finally point our state pointer back to the previous state
	st = st->previous;
	--gamePly;

	assert(pos_is_ok());
}

// Position::do(undo)_null_move() is used to do(undo) a "null move": It flips
// the side to move without executing any move on the board.
void Position::do_null_move(StateInfo& newSt) {

	assert(!checkers());
	assert(&newSt != st);

	std::memcpy(&newSt, st, sizeof(StateInfo));
	newSt.previous = st;
	st = &newSt;

	st->boardKey ^= Zobrist::side;
	prefetch(TT.first_entry(key()));

	st->pliesFromNull = 0;

	sideToMove = ~sideToMove;

	set_check_info(st);

	assert(pos_is_ok());
}

void Position::undo_null_move() {

	assert(!checkers());

	st = st->previous;
	sideToMove = ~sideToMove;
}

// Position::key_after() computes the new hash key after the given move. Needed
// for speculative prefetch.
Key Position::key_after(Move m) const {

	Key k = st->boardKey ^ Zobrist::side;
	Key h = st->handKey;

	Square to = to_sq(m);
	Color us = sideToMove;

	if (type_of(m) == DROP)
	{
		PieceType pt = drop_type(m);

		k += Zobrist::psq[us][pt][to];
		h -= Zobrist::hand[us][pt];
	}
	else
	{
		Square from = from_sq(m);
		Color them = ~us;
		PieceType ptFrom = type_of(piece_on(from));
		PieceType ptTo = (type_of(m) == PROMOTION ? promoted_type(ptFrom) : ptFrom);
		PieceType captured = type_of(piece_on(to));

		if (captured)
		{
			PieceType ptHand = hand_type(captured);

			k -= Zobrist::psq[them][captured][to];
			h += Zobrist::hand[us][ptHand];
		}

		k += Zobrist::psq[us][ptTo][to] - Zobrist::psq[us][ptFrom][from];
	}

	return k + h;
}

// Position::see_ge (Static Exchange Evaluation Greater or Equal) tests if the
// SEE value of move is greater or equal to the given value. We'll use an
// algorithm similar to alpha-beta pruning with a null window.
bool Position::see_ge(Move m, Value threshold) const {

	assert(is_ok(m));

	bool isDrop = (type_of(m) == DROP);
	Square from = from_sq(m), to = to_sq(m);
	PieceType nextVictim = isDrop ? drop_type(m) : type_of(piece_on(from));
	Color stm = ~color_of(moved_piece(m)); // First consider opponent's move

	// Values of the pieces taken by us minus opponent's one
	Value balance = Eval::CapturePieceValue[piece_on(to)] - threshold;

	if (balance < VALUE_ZERO)
		return false;

	balance -= Eval::CapturePieceValue[nextVictim];

	if (balance >= VALUE_ZERO)
		return true;

	bool opponentToMove = true;
	Bitboard stmAttackers;
	Bitboard occupied = (isDrop ? pieces() ^ to : pieces() ^ from ^ to);

	// Find all attackers to the destination square, with the moving piece removed,
	// but possibly an X-ray attacker added behind it.
	Bitboard attackers = attackers_to(to, occupied) & occupied;

	while (true)
	{
		// The balance is negative only because we assumed we could win
		// the last piece for free. We are truly winning only if we can
		// win the last piece _cheaply enough_. Test if we can actually
		// do this otherwise "give up".
		assert(balance < VALUE_ZERO);

		stmAttackers = attackers & pieces(stm);

		// Don't allow pinned pieces to attack pieces except the king as long all
		// pinners are on their original square.
		if (!(st->pinnersForKing[stm] & ~occupied))
			stmAttackers &= ~st->blockersForKing[stm];

		// If we have no more attackers we must give up
		if (!stmAttackers)
			break;

		// Locate and remove the next least valuable attacker
		nextVictim = min_attacker<PAWN>(*this, stm, to, stmAttackers, occupied, attackers);

		if (nextVictim == KING)
		{
			// Our only attacker is the king. If the opponent still has
			// attackers we must give up. Otherwise we make the move and
			// (having no more attackers) the opponent must give up.
			if (!(attackers & pieces(~stm)))
				opponentToMove = !opponentToMove;
			break;
		}

		// Assume the opponent can win the next piece for free and switch sides
		balance = opponentToMove ? balance + Eval::CapturePieceValue[nextVictim]
								 : balance - Eval::CapturePieceValue[nextVictim];
		opponentToMove = !opponentToMove;

		// If balance is negative after receiving a free piece then give up
		if (balance < VALUE_ZERO)
			break;

		// Complete the process of switching sides. The first line swaps
		// all negative numbers with non-negative numbers. The compiler
		// probably knows that it is just the bitwise negation ~balance.
		balance = -balance - 1;
		stm = ~stm;
	}

	// If the opponent gave up we win, otherwise we lose.
	return opponentToMove;
}

RepetitionState Position::is_repetition(int ply) const {

	assert(st->pliesFromNull >= 0);

	const int repPly = 16;
	int end = std::min(repPly, st->pliesFromNull);

	// 4手かけないと千日手にはならないから、4手前から調べていく。
	if (end < 4)
		return REPETITION_NONE;

	StateInfo* stp = st->previous->previous;

	// 盤上の駒が同一である局面が出現した回数
	int cnt = 0;

	for (int i = 4; i <= end; i += 2)
	{
		stp = stp->previous->previous;

		// board_key : 盤上の駒のみのhash(手駒を除く)
		// 盤上の駒が同じ状態であるかを判定する。
		if (stp->boardKey == st->boardKey)
		{
			// 手駒が一致するなら同一局面である。(2手ずつ遡っているので手番は同じである)
			if (stp->hand == st->hand)
			{
				// root(==ply)より遡る(ply <= i)なら2度出現(cnt == 2)する必要がある。
				// rootより遡らない(ply > i)なら1度目の出現(cnt == 1)で千日手と判定する。
				if (++cnt + (ply > i) == 2)
				{
					// 自分が王手をしている連続王手の千日手なのか？
					if (i <= st->continuousCheck[sideToMove])
						return REPETITION_LOSE;

					// 相手が王手をしている連続王手の千日手なのか？
					if (i <= st->continuousCheck[~sideToMove])
						return REPETITION_WIN;

					return REPETITION_DRAW;
				}
			}
			else 
			{
				// 優等局面か劣等局面であるか。(手番が相手番になっている場合はいま考えない)
				if (hand_is_equal_or_superior(st->hand, stp->hand))
					return REPETITION_SUPERIOR;

				if (hand_is_equal_or_superior(stp->hand, st->hand))
					return REPETITION_INFERIOR;
			}
		}
	}

	// 同じhash keyの局面が見つからなかったので…。
	return REPETITION_NONE;
}

bool Position::pos_is_ok() const {

#if 1
	int pieceCount[PIECE_KIND] = { 0 }; // カウント用の変数

	// 1. 盤上の駒と手駒を合わせて40駒あるか
	int total = 0;
	for (Square s = SQ_11; s <= SQ_99; ++s)
	{
		PieceType pt = type_of(piece_on(s));
		if (pt != NO_PIECE_TYPE)
		{
			++total;
			++pieceCount[hand_type(pt)];
		}
	}
	for (Color c = BLACK; c <= WHITE; ++c)
		for (PieceType pt = PAWN; pt <= GOLD; ++pt)
		{
			int ct = count_of(hand[c], pt);
			total += ct;
			pieceCount[pt] += ct;
		}

	if (total != 40)
		return false;

	// 2. それぞれの駒の枚数は合っているか
	for (PieceType pt = NO_PIECE_TYPE; pt <= GOLD; ++pt)
		if (pieceCount[pt] != get_max_number(pt))
			return false;

#endif

	// 3. 王手している駒
	if (st->checkersBB != (attackers_to(king_square(sideToMove)) & pieces(~sideToMove)))
		return false;

	// 4. 相手玉が取れるということはないか
	if (attackers_to(king_square(~sideToMove)) & pieces(sideToMove))
		return false;

	// 5. occupied bitboardは合っているか
	if ((pieces() != (pieces(BLACK) | pieces(WHITE))) || (pieces(BLACK) & pieces(WHITE)))
		return false;

	return true;
}

bool Position::legal_drop_pawn(Square to) const {

	Color us = sideToMove;
	Color them = ~us;

	// 二歩
	if (pieces(us, PAWN) & FileBB[file_of(to)])
		return false;

	// 二歩でなく、打ち歩王手でないなら良し
	if (pawn_attacks(us, to) != Bitboard(king_square(them)))
		return true;

	// ここから打ち歩詰めのチェック
	// 自駒がtoに利いてなければ詰まない
	if (!(attackers_to(to) & pieces(us)))
		return true;

	// 玉以外の相手駒がtoに利いていて、打ち歩を取ることができるなら詰まない
	if ((attackers_to(to) & (pieces(them) ^ king_square(them)))
		& (~pinned_pieces(them) | FileBB[file_of(to)]))
		return true;

	// 打ち歩を取れなくても避けられれば詰まない
	Square ksq = king_square(them);
	Bitboard escapeBB = (king_attacks(ksq) & ~pieces(them)) ^ to;
	Bitboard occ = pieces() ^ to;
	while (escapeBB)
	{
		Square ksqTo = escapeBB.pop();

		if (!(attackers_to(ksqTo, occ) & pieces(us)))
			return true;
	}

	// 打ち歩詰め
	return false;
}