#include <algorithm>

#include "position.h"
#include "psq.h"

PsqIndex PsqIndex::hand[COLOR_NB][PIECE_HAND_NB];
PsqIndex PsqIndex::psq[SQUARE_NB][PIECE_NB];
Piece PsqIndex::indexToPiece[PSQ_MAX];
Square PsqIndex::indexToSquare[PSQ_MAX];

PsqPair PsqPair::hand[COLOR_NB][PIECE_HAND_NB][19];
PsqPair PsqPair::psq[SQUARE_NB][PIECE_NB];
PsqPair PsqPair::allPairs[PSQ_MAX];

void PsqIndex::init()
{
	// 1. calculate an index of the hand pieces
	PsqIndex handIndex(0);
	for (Color c = BLACK; c <= WHITE; ++c)
		for (PieceType pt = PAWN; pt <= GOLD; ++pt)
		{
			hand[c][pt] = PsqIndex(handIndex - 1);

			for (int num = 1; num <= get_max_number(pt); ++num) {
				indexToSquare[handIndex] = SQ_NONE;
				indexToPiece[handIndex] = Piece(pt + (16 * c));
				handIndex = PsqIndex(handIndex + 1);
			}
		}

	assert(handIndex == 76);

	// 2. calculate an index of the board pieces
	PsqIndex boardIndex(76);
	for (Color c = BLACK; c <= WHITE; ++c)
		for (PieceType pt = PAWN; pt <= DRAGON; ++pt)
		{
			if (pt == KING)
				continue;

			for (Square s = SQ_11; s <= SQ_99; ++s) 
			{
				if (pt == PAWN || pt == LANCE)
					if (relative_rank(c, rank_of(s)) == RANK_1)
						continue;

				if (pt == KNIGHT && relative_rank(c, rank_of(s)) <= RANK_2)
					continue;

				Piece pc = make_piece(c, pt);
				indexToPiece[boardIndex] = pc;
				indexToSquare[boardIndex] = s;
				psq[s][pc] = boardIndex;
				boardIndex = PsqIndex(boardIndex + 1);
			}
		}

	assert(boardIndex == PSQ_MAX);
}

void PsqPair::init()
{
	PsqIndex::init();

	// 1. hand pieces
	for (Color c = BLACK; c <= WHITE; ++c) {
		for (PieceType pt = PAWN; pt <= GOLD; ++pt) {
			for (int num = 1; num <= get_max_number(pt); ++num)
			{
				PsqIndex indexBlack = PsqIndex::of_hand(c, pt, num);
				PsqIndex indexWhite = PsqIndex::of_hand(~c, pt, num);
				PsqPair psqPair(indexBlack, indexWhite);
				hand[c][pt][num] = psqPair;
				allPairs[psqPair.black()] = psqPair;
			}
		}
	}

	// 2. board pieces
	for (Color c = BLACK; c <= WHITE; ++c)
		for (PieceType pt = PAWN; pt <= DRAGON; ++pt)
		{
			if (pt == KING)
				continue;

			for (Square s = SQ_11; s <= SQ_99; ++s)
			{
				if (pt == PAWN || pt == LANCE)
					if (relative_rank(c, rank_of(s)) == RANK_1)
						continue;

				if (pt == KNIGHT && relative_rank(c, rank_of(s)) <= RANK_2)
					continue;

				Piece pc = make_piece(c, pt);
				Piece invPc = make_piece(~c, pt);
				PsqIndex indexBlack = PsqIndex::of_board(pc, s);
				PsqIndex indexWhite = PsqIndex::of_board(invPc, inverse(s));
				PsqPair psqPair(indexBlack, indexWhite);
				psq[s][pc] = psqPair;
				allPairs[psqPair.black()] = psqPair;
			}
		}
}

PsqList::PsqList(const Position& pos)
	: size(0) {

	// 1. copy the hand pieces
	hand[BLACK] = pos.hand_of(BLACK);
	hand[WHITE] = pos.hand_of(WHITE);

	// 2. add an index of the hand pieces
	for (Color c = BLACK; c <= WHITE; ++c) {
		for (PieceType pt = PAWN; pt <= GOLD; ++pt) {
			for (int i = 1, n = count_of(hand[c], pt); i <= n; ++i)
			{
				list[size] = PsqPair::of_hand(c, pt, i);
				handIndex[c][pt][i] = size;
				++size;
			}
		}
	}

	// 3. add an index of the board pieces
	Bitboard pieces = (pos.pieces() ^ pos.king_square(BLACK)) ^ pos.king_square(WHITE);
	while (pieces)
	{
		Square sq = pieces.pop();
		Piece piece = pos.piece_on(sq);
		list[size] = PsqPair::of_board(piece, sq);
		index[sq] = size;
		++size;
	}

	assert(is_ok());
}

bool PsqList::is_ok() const
{
	// 1. 要素数のチェック
	if (size > maxSize)
		return false;

	// 2. 持ち駒の要素がリストに含まれているかをチェックする
	for (Color c = BLACK; c <= WHITE; ++c) {
		for (PieceType pt = PAWN; pt <= GOLD; ++pt) {
			for (int num = 1; num <= count_of(hand[c], pt); ++num) {
				// 持ち駒の要素がリストに含まれているか探す
				auto iter = std::find_if(begin(), end(), [&](PsqPair item) {
					return item.black() == PsqIndex::of_hand(c, pt, num);
				});

				// a. あるべき要素がリストになかった
				if (iter == end())
					return false;

				// b. hand_index_との整合性がとれていなかった
				if (handIndex[c][pt][num] != (iter - begin()))
					return false;
			}
		}
	}

	// 3. index_が正しい場所を示しているかチェック
	for (size_t i = 0; i < size; ++i)
	{
		PsqPair item = list[i];

		if (item.square() == SQ_NONE)
			continue;

		if (index[item.square()] != static_cast<int>(i))
			return false;
	}

	return true;
}