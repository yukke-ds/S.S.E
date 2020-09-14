#ifndef PSQ_H_INCLUDED
#define PSQ_H_INCLUDED

#include <vector>
#include <cstdlib>

#include "types.h"
using std::vector;

class Position;

enum {
	PSQ_MAX = 2110,
	PSQ_MIN = 0
};

class PsqIndex {
public:
	explicit constexpr PsqIndex(int i = 0) : index(i) {}
	constexpr operator int() const { return index; }

	Piece piece() const;
	Square square() const;

	static PsqIndex of_hand(Color c, PieceType pt, int num);
	static PsqIndex of_board(Piece pc, Square s);

	static void init();

private:
	int index;
	static PsqIndex hand[COLOR_NB][PIECE_HAND_NB];
	static PsqIndex psq[SQUARE_NB][PIECE_NB];
	static Piece indexToPiece[PSQ_MAX];
	static Square indexToSquare[PSQ_MAX];
};

class PsqPair {
public:
	PsqPair() {}

	PsqIndex get(Color c) const;
	PsqIndex black() const;
	PsqIndex white() const;
	Piece piece() const;
	Square square() const;

	static PsqPair of_hand(Color c, PieceType pt, int num);
	static PsqPair of_board(Piece p, Square s);
	PsqPair* all_pairs();

	static void init();

private:

	PsqPair(PsqIndex index_black, PsqIndex index_white) {
		pair[BLACK] = index_black;
		pair[WHITE] = index_white;
	}

	PsqIndex pair[COLOR_NB];

	// 歩は最大18枚持ち駒にできるので、配列は19要素確保する。
	static PsqPair hand[COLOR_NB][PIECE_HAND_NB][19];
	static PsqPair psq[SQUARE_NB][PIECE_NB];
	static PsqPair allPairs[PSQ_MAX];
};

class PsqList {
public:
	PsqList() : size(0) {}
	explicit PsqList(const Position& pos);

	const PsqPair* begin() const;
	const PsqPair* end() const;

	size_t get_size() const;
	PsqPair operator[](size_t i) const;

	void do_move(Move move);
	void undo_move(Move move);
	bool is_ok() const;

	static bool two_lists_have_same_items(const PsqList& list1, const PsqList& list2);

private:
	static constexpr int maxSize = 38;
	size_t size;
	Hand hand[COLOR_NB];
	PsqPair list[maxSize];
	int index[SQUARE_NB];
	int handIndex[COLOR_NB][PIECE_HAND_NB][19];
};

inline Piece PsqIndex::piece() const { 
	return indexToPiece[index]; 
}

inline 	Square PsqIndex::square() const {
	return indexToSquare[index]; 
}

inline PsqIndex PsqIndex::of_hand(Color c, PieceType pt, int num) {
	PsqIndex idx = PsqIndex(hand[c][pt] + num);
	assert(0 <= idx && idx <= 75);
	return idx;
}

inline PsqIndex PsqIndex::of_board(Piece pc, Square s) {
	PsqIndex idx = psq[s][pc];
	assert(76 <= idx && idx < PSQ_MAX);
	return idx;
}

inline PsqIndex PsqPair::get(Color c) const {
	return pair[c];
}

inline PsqIndex PsqPair::black() const {
	return pair[BLACK];
}

inline PsqIndex PsqPair::white() const {
	return pair[WHITE];
}

inline Piece PsqPair::piece() const {
	return pair[BLACK].piece();
}

inline Square PsqPair::square() const {
	return pair[BLACK].square();
}

inline PsqPair PsqPair::of_hand(Color c, PieceType pt, int num) {
	assert(1 <= num && num <= get_max_number(pt));
	return hand[c][pt][num];
}

inline PsqPair PsqPair::of_board(Piece p, Square s) {
	return psq[s][p];
}

inline PsqPair* PsqPair::all_pairs() {
	return allPairs;
}

inline const PsqPair* PsqList::begin() const {
	return &list[0];
}

inline const PsqPair* PsqList::end() const {
	return &list[0] + size;
}

inline size_t PsqList::get_size() const {
	assert(size <= maxSize);
	return size;
}

inline PsqPair PsqList::operator[](size_t i) const {
	assert(i < size && size <= maxSize);
	return list[i];
}

#endif // ifndef PSQ_H_INCLUDED