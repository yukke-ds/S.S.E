#ifndef PROGRESS_H_INCLUDED
#define PROGRESS_H_INCLUDED

#include "position.h"
#include "psq.h"

//#define KIF

class Position;

class Progress {
	
public:
	static constexpr int weightScale = 1 << 16;

	static void read_weights();
	static void learn_parameters();
	static double estimate_progress(const Position& pos, const PsqList& psqList);
	static double estimate_progress(const Position& pos);

	static void init();
	static int relation(Square from, Square to);

	static int32_t weights[SQUARE_NB][PSQ_MAX];
private:

	static int Relation[SQUARE_NB][SQUARE_NB];
};

inline int Progress::relation(Square from, Square to) {
	return Relation[from][to];
}

#endif // ifndef PROGRESS_H_INCLUDED