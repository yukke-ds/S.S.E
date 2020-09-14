#ifndef HAND_H_INCLUDED
#define HAND_H_INCLUDED

#include "types.h"

// ����bit�ʒu
constexpr int PieceBits[PIECE_HAND_NB] = {
	0,  0,  8, 12, // 0, ��, ��, �j
	16, 20, 24, 28 // ��, �p, ��, ��
};

// ���̎������\������̂ɕK�v��bit����mask(�Ⴆ��3bit�Ȃ�2��3��-1��7)
constexpr int PieceBitMask[PIECE_HAND_NB] = {
	0, 31, 7, 7, // ��:5bit, ���j:3bit
	7,  3, 3, 7 // ��:3bit, �p��:2bit, ��:3bit
};

// ��:11111, ��:11100000000, ...
constexpr int PieceBitMaskSum[PIECE_HAND_NB] = {
	0, PieceBitMask[PAWN] << PieceBits[PAWN], PieceBitMask[LANCE] << PieceBits[LANCE],
	PieceBitMask[KNIGHT] << PieceBits[KNIGHT], PieceBitMask[SILVER] << PieceBits[SILVER], 
	PieceBitMask[BISHOP] << PieceBits[BISHOP], PieceBitMask[ROOK] << PieceBits[ROOK],
	PieceBitMask[GOLD] << PieceBits[GOLD]
};

// ��̖������i�[����Ă���bit��1�ƂȂ��Ă���MASK(���𓾂�Ƃ��Ɏg��)
constexpr int32_t HandBitMask = PieceBitMaskSum[PAWN] | PieceBitMaskSum[LANCE] | 
								PieceBitMaskSum[KNIGHT] | PieceBitMaskSum[SILVER] |
								PieceBitMaskSum[BISHOP] | PieceBitMaskSum[ROOK] | PieceBitMaskSum[GOLD];

// �]�点�Ă���bit�̏W��
constexpr int32_t HandBorrowMask = (HandBitMask << 1) & ~HandBitMask;

// ���ȊO�̎�����̃}�X�N = 0 111 00 11 00 11 0 111 0 111 0 111 000 00000
//							��4    ��2   �p2   ��4   �j4   ��4      ��18
constexpr uint32_t ExceptHandPawnMask = 
	( (PieceBitMask[LANCE] << PieceBits[LANCE])
	| (PieceBitMask[KNIGHT] << PieceBits[KNIGHT])
	| (PieceBitMask[SILVER] << PieceBits[SILVER])
	| (PieceBitMask[BISHOP] << PieceBits[BISHOP])
	| (PieceBitMask[ROOK] << PieceBits[ROOK])
	| (PieceBitMask[GOLD] << PieceBits[GOLD]));

// Piece(��, ��, �j, ��, ��, �p, ��)�����ɕϊ�����e�[�u��
constexpr Hand PieceToHand[PIECE_HAND_NB] = {
	(Hand)0, (Hand)(1 << PieceBits[PAWN]),
	(Hand)(1 << PieceBits[LANCE]), (Hand)(1 << PieceBits[KNIGHT]),
	(Hand)(1 << PieceBits[SILVER]), (Hand)(1 << PieceBits[BISHOP]),
	(Hand)(1 << PieceBits[ROOK]), (Hand)(1 << PieceBits[GOLD])
};

// ������pt�̖�����Ԃ�
inline int count_of(Hand h, PieceType pt) {
	assert(PAWN <= pt && pt <= GOLD);
	return (h >> PieceBits[pt]) & PieceBitMask[pt];
}

// ����̋�������Ă��邩�m�F
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

// �������pt��c���₷�i�f�t�H���g��c = 1�j
inline void add_hand(Hand& hand, PieceType pt, int c = 1) {
	hand = (Hand)(hand + PieceToHand[pt] * c);
}

// �������pt��c���炷�i�f�t�H���g��c = 1�j
inline void sub_hand(Hand& hand, PieceType pt, int c = 1) {
	hand = (Hand)(hand - PieceToHand[pt] * c);
}

// ���h1�̂ق���h2���D��Ă��邩(���ׂĂ̎�ނ̎�h2�̂���ȏ゠��)
inline bool hand_is_equal_or_superior(Hand h1, Hand h2) { 
	return ((h1 - h2) & HandBorrowMask) == 0;
}

#endif // ifndef HAND_H_INCLUDED
