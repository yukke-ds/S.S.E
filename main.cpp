#include "bitboard.h"
#include "position.h"
#include "progress.h"
#include "psq.h"
#include "search.h"
#include "thread.h"
#include "tt.h"
#include "usi.h"

int main(int argc, char* argv[]) {

	USI::init(Options);
	Bitboards::init();
	Position::init();
	Search::init();
	Evaluate::init();
	Progress::init();
	PsqPair::init();
	TT.resize(Options["Hash"]);
	Threads.set(Options["Threads"]);
	Search::clear(); // After threads are up

	USI::loop(argc, argv);

	Threads.set(0);
	return 0;
}