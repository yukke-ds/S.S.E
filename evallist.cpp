#include "evallist.h"

const ExtEvalIndex KPP_BOARD_INDEX[PIECE_NB] = {
	{ EVAL_INDEX_ZERO, EVAL_INDEX_ZERO },
	{ F_PAWN, E_PAWN },
	{ F_LANCE , E_LANCE },
	{ F_KNIGHT, E_KNIGHT },
	{ F_SILVER, E_SILVER },
	{ F_BISHOP, E_BISHOP },
	{ F_ROOK, E_ROOK },
	{ F_GOLD, E_GOLD },
	{ F_KING, E_KING },
	{ F_GOLD, E_GOLD }, // ����
	{ F_GOLD, E_GOLD }, // ����
	{ F_GOLD, E_GOLD }, // ���j
	{ F_GOLD, E_GOLD }, // ����
	{ F_HORSE, E_HORSE }, // �n
	{ F_DRAGON, E_DRAGON }, // ��
	{ EVAL_INDEX_ZERO, EVAL_INDEX_ZERO }, // ���̐���͂Ȃ�

	// ��肩�猩���ꍇ�Bf��e������ւ��B
	{ EVAL_INDEX_ZERO, EVAL_INDEX_ZERO },
	{ E_PAWN, F_PAWN },
	{ E_LANCE, F_LANCE },
	{ E_KNIGHT, F_KNIGHT },
	{ E_SILVER , F_SILVER },
	{ E_BISHOP, F_BISHOP },
	{ E_ROOK, F_ROOK },
	{ E_GOLD, F_GOLD },
	{ E_KING, F_KING },
	{ E_GOLD, F_GOLD }, // ����
	{ E_GOLD, F_GOLD }, // ����
	{ E_GOLD, F_GOLD }, // ���j
	{ E_GOLD, F_GOLD }, // ����
	{ E_HORSE, F_HORSE }, // �n
	{ E_DRAGON, F_DRAGON }, // ��
	{ EVAL_INDEX_ZERO, EVAL_INDEX_ZERO }, // ���̐���͂Ȃ�
};

const ExtEvalIndex KPP_HAND_INDEX[COLOR_NB][PIECE_HAND_NB] = {
	{
		{ EVAL_INDEX_ZERO, EVAL_INDEX_ZERO },
		{ F_HAND_PAWN, E_HAND_PAWN },
		{ F_HAND_LANCE, E_HAND_LANCE },
		{ F_HAND_KNIGHT, E_HAND_KNIGHT },
		{ F_HAND_SILVER, E_HAND_SILVER },
		{ F_HAND_BISHOP, E_HAND_BISHOP },
		{ F_HAND_ROOK, E_HAND_ROOK },
		{ F_HAND_GOLD, E_HAND_GOLD },
	},
	{
		{ EVAL_INDEX_ZERO, EVAL_INDEX_ZERO },
		{ E_HAND_PAWN, F_HAND_PAWN },
		{ E_HAND_LANCE, F_HAND_LANCE },
		{ E_HAND_KNIGHT, F_HAND_KNIGHT },
		{ E_HAND_SILVER, F_HAND_SILVER },
		{ E_HAND_BISHOP, F_HAND_BISHOP },
		{ E_HAND_ROOK, F_HAND_ROOK },
		{ E_HAND_GOLD, F_HAND_GOLD },
	},
};