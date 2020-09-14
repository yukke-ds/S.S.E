#ifndef POSITION_H_INCLUDED
#define POSITION_H_INCLUDED

#include <cassert>
#include <deque>
#include <memory> // For std::unique_ptr
#include <string>

#include "bitboard.h"
#include "evaluate.h"
#include "hand.h"
#include "types.h"

struct MovedPieceInfo {

	ChangedEvalIndex changedPiece[2];
	PieceNumber movedPieceNumber[2];
	int numPieceMoved;
};

struct StateInfo {

	// Copied when making a move
	int pliesFromNull;
	int continuousCheck[COLOR_NB];

	// Not copied when making a move
	Key boardKey;							// 盤面のハッシュキー
	Key handKey;							// 持ち駒のハッシュキー
	Hand hand;								// 手番側の持ち駒
	Bitboard checkersBB;					// 手番側の玉へ王手している駒のビットボード
	PieceType capturedType;					// この局面で捕獲された駒
	Evaluate eval;							// この局面の局面評価クラス
	Value material;							// この局面の駒価値 
	MovedPieceInfo mpi;						// 差分計算のために動かした駒の情報を保存
	StateInfo* previous;					// １つ前の局面に遡るためのポインタ
	Bitboard blockersForKing[COLOR_NB];		// 玉を飛角等から守る駒のビットボード
	Bitboard pinnersForKing[COLOR_NB];		// 他駒の移動によって王手となる駒のビットボード
	Bitboard checkSquares[PIECE_TYPE_NB];	// 敵玉が王手となるマスのbitboard
};

// In a std::deque references to elements are unaffected upon resizing
typedef std::unique_ptr<std::deque<StateInfo>> StateListPtr;

// Position class stores information regarding the board representation as
// pieces, side to move, hash keys, castling info, etc. Important methods are
// do_move() and undo_move(), used by the search to update node info when
// traversing the search tree.
class Thread;

class Position {
public:
	static void init();

	Position() = default;
	Position& operator=(const Position&) = delete;

	// SFEN string input/output
	Position& set(const std::string& sfenStr, StateInfo* si, Thread* th);
	const std::string sfen() const;

	// Position representation
	Bitboard pieces() const;
	Bitboard pieces(PieceType pt) const;
	Bitboard pieces(PieceType pt1, PieceType pt2) const;
	Bitboard pieces(Color c) const;
	Bitboard pieces(Color c, PieceType pt) const;
	Bitboard pieces(Color c, PieceType pt1, PieceType pt2) const;
	Bitboard gold_pieces() const;
	Bitboard gold_pieces(Color c) const;
	Piece piece_on(Square s) const;
	Square king_square(Color c) const;
	bool empty(Square s) const;

	// Checking
	Bitboard checkers() const;
	Bitboard discovered_check_candidates() const;
	Bitboard pinned_pieces(Color c) const;
	Bitboard check_squares(PieceType pt) const;

	// Attacks to/from a given square
	Bitboard attackers_to(Square s) const;
	Bitboard attackers_to(Square s, Bitboard occupied) const;
	Bitboard attacks_from(Piece pc, Square s) const;
	Bitboard slider_blockers(Bitboard sliders, Square s, Bitboard& pinners) const;

	// Properties of moves
	bool legal(Move m) const;
	bool pseudo_legal(const Move m) const;
	bool capture(Move m) const;
	bool capture_or_promotion(Move m) const;
	bool pawn_promotion(Move m) const;
	bool capture_or_pawn_promotion(Move m) const;
	bool gives_check(Move m) const;
	Piece moved_piece(Move m) const;
	PieceType captured_piece() const;

	// Doing and undoing moves
	void do_move(Move m, StateInfo& newSt);
	void do_move(Move m, StateInfo& st, bool givesCheck);
	void undo_move(Move m);
	void do_null_move(StateInfo& st);
	void undo_null_move();

	// Using positional evaluation
	const EvalList* eval_list() const;
	EvalIndex kpp_hand_index_of(Color c, PieceType pt) const;
	PieceNumber piece_no_of(Color c, PieceType pt) const;
	PieceNumber piece_no_of(Square sq) const;

	// Static exchange evaluation
	bool see_ge(Move m, Value threshhold = VALUE_ZERO) const;

	// Accessing hash keys
	Key key() const;
	Key key_after(Move m) const;
	Key board_key() const;
	Key hand_key() const;

	// Repetition check
	RepetitionState is_repetition(int ply) const;

	// Other properties of the position
	Color side_to_move() const;
	Hand hand_of(Color c) const;
	int game_ply() const;
	StateInfo* state() const;
	Thread* this_thread() const;

	// Position consistency check
	bool pos_is_ok() const;
	bool legal_drop_pawn(Square to) const;

private:
	// Initialization helperes (used while setting up a position)
	void set_state(StateInfo* si) const;
	void set_check_info(StateInfo* si) const;

	// Other helpers
	void put_piece(Color c, PieceType pt, Square s);
	void remove_piece(Color c, PieceType pt, Square s);
	void put_hand(Hand& h, PieceType pt, int pn);

	// Data members
	Piece board[SQUARE_NB];
	Bitboard byTypeBB[PIECE_TYPE_NB];
	Bitboard byColorBB[COLOR_NB];
	Hand hand[COLOR_NB];
	Square kingSquare[COLOR_NB];
	int gamePly;
	Color sideToMove;
	Thread* thisThread;
	StateInfo* st;
	EvalList evalList;
};

inline Color Position::side_to_move() const {
	return sideToMove;
}

inline bool Position::empty(Square s) const {
	return board[s] == NO_PIECE;
}

inline Hand Position::hand_of(Color c) const {
	return hand[c];
}

inline Piece Position::piece_on(Square s) const {
	return board[s];
}

inline Piece Position::moved_piece(Move m) const {
	return type_of(m) == DROP ? make_piece(sideToMove, drop_type(m)) : board[from_sq(m)];
}

inline Bitboard Position::pieces() const {
	return byTypeBB[ALL_PIECES];
}

inline Bitboard Position::pieces(PieceType pt) const {
	return byTypeBB[pt];
}

inline Bitboard Position::pieces(PieceType pt1, PieceType pt2) const {
	return byTypeBB[pt1] | byTypeBB[pt2];
}

inline Bitboard Position::pieces(Color c) const {
	return byColorBB[c];
}

inline Bitboard Position::pieces(Color c, PieceType pt) const {
	return byColorBB[c] & byTypeBB[pt];
}

inline Bitboard Position::pieces(Color c, PieceType pt1, PieceType pt2) const {
	return byColorBB[c] & (byTypeBB[pt1] | byTypeBB[pt2]);
}

inline Bitboard Position::gold_pieces() const {
	return (byTypeBB[GOLD] | byTypeBB[PRO_PAWN] | byTypeBB[PRO_LANCE] | byTypeBB[PRO_KNIGHT] | byTypeBB[PRO_SILVER]);
}

inline Bitboard Position::gold_pieces(Color c) const {
	return byColorBB[c] & gold_pieces();
}

inline Square Position::king_square(Color c) const {
	return kingSquare[c];
}

inline Bitboard Position::attacks_from(Piece pc, Square s) const {
	return attacks_bb(pc, s, byTypeBB[ALL_PIECES]);
}

inline Bitboard Position::attackers_to(Square s) const {
	return attackers_to(s, byTypeBB[ALL_PIECES]);
}

inline Bitboard Position::checkers() const {
	return st->checkersBB;
}

inline Bitboard Position::discovered_check_candidates() const {
	return st->blockersForKing[~sideToMove] & pieces(sideToMove);
}

inline Bitboard Position::pinned_pieces(Color c) const {
	return st->blockersForKing[c] & pieces(c);
}

inline Bitboard Position::check_squares(PieceType pt) const {
	return st->checkSquares[pt];
}

inline const EvalList* Position::eval_list() const {
	return &evalList;
}

inline EvalIndex Position::kpp_hand_index_of(Color c, PieceType pt) const {
	int counts = count_of(hand[c], pt);
	assert(counts > 0);
	return (EvalIndex)(KPP_HAND_INDEX[c][pt].fb + counts - 1);
}

inline PieceNumber Position::piece_no_of(Color c, PieceType pt) const {
	return evalList.piece_no_of_hand(kpp_hand_index_of(c, pt));
}

inline PieceNumber Position::piece_no_of(Square sq) const {
	assert(piece_on(sq) != NO_PIECE);
	PieceNumber pn = evalList.piece_no_of_board(sq);
	assert(is_ok(pn));
	return pn;
}

inline Key Position::key() const {
	return st->boardKey + st->handKey;
}

inline Key Position::board_key() const {
	return st->boardKey;
}

inline Key Position::hand_key() const {
	return st->handKey;
}

inline int Position::game_ply() const {
	return gamePly;
}

inline bool Position::capture(Move m) const {

	assert(is_ok(m));
	return (!empty(to_sq(m)) && type_of(m) != DROP);
}

inline bool Position::capture_or_promotion(Move m) const {

	assert(is_ok(m));
	return (m & PROMOTION) || capture(m);
}

inline bool Position::pawn_promotion(Move m) const {

	assert(is_ok(m));
	return (type_of(moved_piece(m)) == PAWN && type_of(m) == PROMOTION);
}

inline bool Position::capture_or_pawn_promotion(Move m) const {

	assert(is_ok(m));
	return pawn_promotion(m) || capture(m);
}

inline PieceType Position::captured_piece() const {
	return st->capturedType;
}

inline StateInfo* Position::state() const {
	return st;
}

inline Thread* Position::this_thread() const {
	return thisThread;
}

// 駒を配置して、内部的に保持しているBitboardも更新する
inline void Position::put_piece(Color c, PieceType pt, Square s) {

	board[s] = make_piece(c, pt);
	byTypeBB[ALL_PIECES] |= s;
	byTypeBB[pt] |= s;
	byColorBB[c] |= s;

	if (pt == KING)
		kingSquare[c] = s;
}

// 駒を盤面から取り除き、内部的に保持しているBitboardも更新する
inline void Position::remove_piece(Color c, PieceType pt, Square s) {

	byTypeBB[ALL_PIECES] ^= s;
	byTypeBB[pt] ^= s;
	byColorBB[c] ^= s;
	board[s] = NO_PIECE;
}

// 持ち駒にpcをc枚加える
inline void Position::put_hand(Hand& h, PieceType pt, int pn) {
	h = (Hand)(h + PieceToHand[pt] * pn);
}

inline void Position::do_move(Move m, StateInfo& newSt) {
	do_move(m, newSt, gives_check(m));
}

#endif // ifndef POSITION_H_INCLUDED