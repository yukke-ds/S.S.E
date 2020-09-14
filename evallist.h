#ifndef EVALLIST_H_INCLUDED
#define EVALLIST_H_INCLUDED

#include "types.h"

enum EvalIndex : int16_t {

	EVAL_INDEX_ZERO = 0,

	F_HAND_PAWN = 1,
	E_HAND_PAWN = F_HAND_PAWN + 19,
	F_HAND_LANCE = E_HAND_PAWN + 19,
	E_HAND_LANCE = F_HAND_LANCE + 5,
	F_HAND_KNIGHT = E_HAND_LANCE + 5,
	E_HAND_KNIGHT = F_HAND_KNIGHT + 5,
	F_HAND_SILVER = E_HAND_KNIGHT + 5,
	E_HAND_SILVER = F_HAND_SILVER + 5,
	F_HAND_GOLD = E_HAND_SILVER + 5,
	E_HAND_GOLD = F_HAND_GOLD + 5,
	F_HAND_BISHOP = E_HAND_GOLD + 5,
	E_HAND_BISHOP = F_HAND_BISHOP + 3,
	F_HAND_ROOK = E_HAND_BISHOP + 3,
	E_HAND_ROOK = F_HAND_ROOK + 3,
	FE_HAND_END = E_HAND_ROOK + 2,

	F_PAWN = FE_HAND_END,
	E_PAWN = F_PAWN + 81,
	F_LANCE = E_PAWN + 81,
	E_LANCE = F_LANCE + 81,
	F_KNIGHT = E_LANCE + 81,
	E_KNIGHT = F_KNIGHT + 81,
	F_SILVER = E_KNIGHT + 81,
	E_SILVER = F_SILVER + 81,
	F_GOLD = E_SILVER + 81,
	E_GOLD = F_GOLD + 81,
	F_BISHOP = E_GOLD + 81,
	E_BISHOP = F_BISHOP + 81,
	F_HORSE = E_BISHOP + 81,
	E_HORSE = F_HORSE + 81,
	F_ROOK = E_HORSE + 81,
	E_ROOK = F_ROOK + 81,
	F_DRAGON = E_ROOK + 81,
	E_DRAGON = F_DRAGON + 81,
	FE_END = E_DRAGON + 81,

	F_KING = FE_END,
	E_KING = F_KING + SQUARE_NB,
};

struct ExtEvalIndex
{
	EvalIndex fb; // from black
	EvalIndex fw; // from white
};

struct ChangedEvalIndex
{
	ExtEvalIndex oldPiece;
	ExtEvalIndex newPiece;
};

// fb = kpp_board_index[pc].fb + sq; // 先手から見たsqにあるpcに対応するEvalIndex
// fw = kpp_board_index[pc].fw + sq; // 後手から見たsqにあるpcに対応するEvalIndex
extern const ExtEvalIndex KPP_BOARD_INDEX[PIECE_NB];
extern const ExtEvalIndex KPP_HAND_INDEX[COLOR_NB][KING];

struct EvalList {

	// 評価関数(FV38型)で用いる駒番号のリスト
	inline EvalIndex* piece_list_fb() const { return const_cast<EvalIndex*>(pieceListFb); }
	inline EvalIndex* piece_list_fw() const { return const_cast<EvalIndex*>(pieceListFw); }
	ExtEvalIndex eval_index(PieceNumber pn) const;

	void put_board_piece(PieceNumber pn, Piece pc, Square sq);
	void put_hand_piece(PieceNumber pn, PieceType pt, Color c, int i);
	PieceNumber piece_no_of_board(Square sq) const;
	PieceNumber piece_no_of_hand(EvalIndex ei) const;

	// pieceListを初期化する。
	// 駒落ちに対応させる時のために、未使用の駒の値はBONA_PIECE_ZEROにしておく。
	// 通常の評価関数を駒落ちの評価関数として流用できる。
	// piece_no_listのほうはデバッグが捗るように-1で初期化。
	void clear()
	{
		for (auto& p : pieceListFb)
			p = EVAL_INDEX_ZERO;
		for (auto& p : pieceListFw)
			p = EVAL_INDEX_ZERO;
		for (auto& h : hand_piece_no_list)
			h = PIECE_NUMBER_NB;
		for (auto& b : board_piece_no_list)
			b = PIECE_NUMBER_NB;
	}

protected:

	// 盤上sqにあるpiece_noの駒のBonaPieceがfb,fwであることを設定する。
	inline void set_piece_on_board(PieceNumber pn, EvalIndex fb, EvalIndex fw, Square sq)
	{
		assert(is_ok(pn));
		pieceListFb[pn] = fb;
		pieceListFw[pn] = fw;
		board_piece_no_list[sq] = pn;
	}

	// 手駒であるpiece_noの駒のBonaPieceがfb,fwであることを設定する。
	inline void set_piece_on_hand(PieceNumber pn, EvalIndex fb, EvalIndex fw)
	{
		assert(is_ok(pn));
		pieceListFb[pn] = fb;
		pieceListFw[pn] = fw;
		hand_piece_no_list[fb] = pn;
	}

	// 駒リスト。駒番号(PieceNumber)いくつの駒がどこにあるのか(EvalIndex)を示す。FV38などで用いる。
	EvalIndex pieceListFb[PIECE_NUMBER_NB];
	EvalIndex pieceListFw[PIECE_NUMBER_NB];

	// あるEvalIndexに対して、その駒番号(PieceNumber)を保持している配列
	PieceNumber hand_piece_no_list[FE_HAND_END];
	PieceNumber board_piece_no_list[SQUARE_NB];
};

inline ExtEvalIndex EvalList::eval_index(PieceNumber pn) const {
	ExtEvalIndex ei;
	ei.fb = pieceListFb[pn];
	ei.fw = pieceListFw[pn];
	return ei;
}

inline void EvalList::put_board_piece(PieceNumber pn, Piece pc, Square sq) {
	set_piece_on_board(pn, EvalIndex(KPP_BOARD_INDEX[pc].fb + sq), EvalIndex(KPP_BOARD_INDEX[pc].fw + inverse(sq)), sq);
}

inline void EvalList::put_hand_piece(PieceNumber pn, PieceType pt, Color c, int i) {
	set_piece_on_hand(pn, EvalIndex(KPP_HAND_INDEX[c][pt].fb + i), EvalIndex(KPP_HAND_INDEX[c][pt].fw + i));
}

inline PieceNumber EvalList::piece_no_of_board(Square sq) const {
	return board_piece_no_list[sq];
}

inline PieceNumber EvalList::piece_no_of_hand(EvalIndex ei) const { 
	return hand_piece_no_list[ei]; 
}

#endif // ifndef EVALLIST_H_INCLUDED