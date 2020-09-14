#ifndef HAND_H_INCLUDED
#define HAND_H_INCLUDED

#include "types.h"

// 手駒のbit位置
constexpr int PieceBits[PIECE_HAND_NB] = {
	0,  0,  8, 12, // 0, 歩, 香, 桂
	16, 20, 24, 28 // 銀, 角, 飛, 金
};

// その持ち駒を表現するのに必要なbit数のmask(例えば3bitなら2の3乗-1で7)
constexpr int PieceBitMask[PIECE_HAND_NB] = {
	0, 31, 7, 7, // 歩:5bit, 香桂:3bit
	7,  3, 3, 7 // 銀:3bit, 角飛:2bit, 金:3bit
};

// 歩:11111, 香:11100000000, ...
constexpr int PieceBitMaskSum[PIECE_HAND_NB] = {
	0, PieceBitMask[PAWN] << PieceBits[PAWN], PieceBitMask[LANCE] << PieceBits[LANCE],
	PieceBitMask[KNIGHT] << PieceBits[KNIGHT], PieceBitMask[SILVER] << PieceBits[SILVER], 
	PieceBitMask[BISHOP] << PieceBits[BISHOP], PieceBitMask[ROOK] << PieceBits[ROOK],
	PieceBitMask[GOLD] << PieceBits[GOLD]
};

// 駒の枚数が格納されているbitが1となっているMASK(駒種を得るときに使う)
constexpr int32_t HandBitMask = PieceBitMaskSum[PAWN] | PieceBitMaskSum[LANCE] | 
								PieceBitMaskSum[KNIGHT] | PieceBitMaskSum[SILVER] |
								PieceBitMaskSum[BISHOP] | PieceBitMaskSum[ROOK] | PieceBitMaskSum[GOLD];

// 余らせてあるbitの集合
constexpr int32_t HandBorrowMask = (HandBitMask << 1) & ~HandBitMask;

// 歩以外の持ち駒のマスク = 0 111 00 11 00 11 0 111 0 111 0 111 000 00000
//							金4    飛2   角2   銀4   桂4   香4      歩18
constexpr uint32_t ExceptHandPawnMask = 
	( (PieceBitMask[LANCE] << PieceBits[LANCE])
	| (PieceBitMask[KNIGHT] << PieceBits[KNIGHT])
	| (PieceBitMask[SILVER] << PieceBits[SILVER])
	| (PieceBitMask[BISHOP] << PieceBits[BISHOP])
	| (PieceBitMask[ROOK] << PieceBits[ROOK])
	| (PieceBitMask[GOLD] << PieceBits[GOLD]));

// Piece(歩, 香, 桂, 銀, 金, 角, 飛)を手駒に変換するテーブル
constexpr Hand PieceToHand[PIECE_HAND_NB] = {
	(Hand)0, (Hand)(1 << PieceBits[PAWN]),
	(Hand)(1 << PieceBits[LANCE]), (Hand)(1 << PieceBits[KNIGHT]),
	(Hand)(1 << PieceBits[SILVER]), (Hand)(1 << PieceBits[BISHOP]),
	(Hand)(1 << PieceBits[ROOK]), (Hand)(1 << PieceBits[GOLD])
};

// 持ち駒ptの枚数を返す
inline int count_of(Hand h, PieceType pt) {
	assert(PAWN <= pt && pt <= GOLD);
	return (h >> PieceBits[pt]) & PieceBitMask[pt];
}

// 特定の駒を持っているか確認
inline uint32_t has(PieceType pt, Hand h) {
	switch (pt) {
	case PAWN: return (h & (PieceBitMask[PAWN] << PieceBits[PAWN]));
	case LANCE: return (h & (PieceBitMask[LANCE] << PieceBits[LANCE]));
	case KNIGHT: return (h & (PieceBitMask[KNIGHT] << PieceBits[KNIGHT]));
	case SILVER: return (h & (PieceBitMask[SILVER] << PieceBits[SILVER]));
	case BISHOP: return (h & (PieceBitMask[BISHOP] << PieceBits[BISHOP]));
	case ROOK: return (h & (PieceBitMask[ROOK] << PieceBits[ROOK]));
	case GOLD: return (h & (PieceBitMask[GOLD] << PieceBits[GOLD]));
	default: assert(false); return 0;
	}
}

inline uint32_t has_except_pawn(Hand h) {
	return (h & ExceptHandPawnMask);
}

// 持ち駒にptをc個増やす（デフォルトはc = 1）
inline void add_hand(Hand& hand, PieceType pt, int c = 1) {
	hand = (Hand)(hand + PieceToHand[pt] * c);
}

// 持ち駒からptをc個減らす（デフォルトはc = 1）
inline void sub_hand(Hand& hand, PieceType pt, int c = 1) {
	hand = (Hand)(hand - PieceToHand[pt] * c);
}

// 手駒h1のほうがh2より優れているか(すべての種類の手駒がh2のそれ以上ある)
inline bool hand_is_equal_or_superior(Hand h1, Hand h2) { 
	return ((h1 - h2) & HandBorrowMask) == 0;
}

#endif // ifndef HAND_H_INCLUDED
