#include <algorithm> // For std::count
#include <cassert>

#include "ann/ann_evaluator.h"
#include "movegen.h"
#include "search.h"
#include "thread.h"
#include "usi.h"

ThreadPool Threads; // Global object

// Thread constructor launches the thread and waits until it goes to sleep
// in idle_loop(). Note that 'searching' and 'exit' should be alredy set.
Thread::Thread(size_t n) : idx(n), stdThread(&Thread::idle_loop, this) {

	wait_for_search_finished();
}

// Thread destructor wakes up the thread in idle_loop() and waits
// for its termination. Thread should be already waiting.
Thread::~Thread() {

	assert(!searching);

	exit = true;
	start_searching();
	stdThread.join();
}

// Thread::clear() reset histories, usually before a new game
void Thread::clear() {

	counterMoves.fill(MOVE_NONE);
	mainHistory.fill(0);

	for (auto& to : contHistory)
		for (auto& h : to)
			h.get()->fill(0);

	contHistory[NO_PIECE][0].get()->fill(Search::CounterMovePruneThreshold - 1);
}

// Thread::start_searching() wakes up the thread that will start the search
void Thread::start_searching() {

	std::lock_guard<Mutex> lk(mutex);
	searching = true;
	cv.notify_one(); // Wake up the thread in idle_loop()
}

// Thread::wait_for_search_finished() blocks on the condition variable
// until the thread has finished searching.
void Thread::wait_for_search_finished() {
	std::unique_lock<Mutex> lk(mutex);
	cv.wait(lk, [&] { return !searching; });
}

// Thread::idle_loop() is where the thread is parked, blocked on the
// condition variable, when it has no work to do.
void Thread::idle_loop() {

	// If OS already scheduled us on a different group than 0 then don't overwrite
	// the choice, eventually we are one of many one-threaded processes running on
	// some Windows NUMA hardware, for instance in fishtest. To make it simple,
	// just check if running threads are below a threshold, in this case all this
	// NUMA machinery is not needed.
	if (Options["Threads"] >= 8)
		WinProcGroup::bindThisThread(idx);

	while (true)
	{
		std::unique_lock<Mutex> lk(mutex);
		searching = false;
		cv.notify_one(); // Wake up anyone waiting for search finished
		cv.wait(lk, [&] { return searching; });

		if (exit)
			return;

		lk.unlock();

		search();
	}
}

// ThreadPool::set() creates/destroys threads to match the requested number.
// Created and launced threads wil go immediately to sleep in idle_loop.
// Upon resizing, threads are recreated to allow for binding if necessary.
void ThreadPool::set(size_t requested) {

	if (size() > 0) { // destroy any existing thread(s)
		main()->wait_for_search_finished();

		while (size() > 0)
			delete back(), pop_back();
	}

	if (requested > 0) { // create new thread(s)
		push_back(new MainThread(0));

		while (size() < requested)
			push_back(new Thread(size()));
		clear();
	}
}

// ThreadPool::clear() sets threadPool data to initial values.
void ThreadPool::clear() {

	for (Thread* th : *this)
		th->clear();

	main()->callsCnt = 0;
	main()->previousScore = VALUE_INFINITE;
}

// ThreadPool::start_thinking() wakes up main thread waiting in idle_loop() and
// returns immediately. Main thread will wake up other threads and start the search.
void ThreadPool::start_thinking(Position& pos, StateListPtr& states,
								const Search::LimitsType& limits, bool ponderMode) {

	main()->wait_for_search_finished();

	stopOnPonderhit = stop = false;
	ponder = ponderMode;
	Search::Limits = limits;
	Search::RootMoves rootMoves;

	for (const auto& m : MoveList<LEGAL>(pos))
		if (limits.searchmoves.empty()
			|| std::count(limits.searchmoves.begin(), limits.searchmoves.end(), m))
			rootMoves.emplace_back(m);

	// After ownership transfer 'states' becomes empty, so if we stop the search
	// and call 'go' again without setting a new position states.get() == NULL.
	assert(states.get() || setupStates.get());

	if (states.get())
		setupStates = std::move(states); // Ownership transfer, states is now empty

	// We use Position::set() to set root position across threads. But there are
	// some StateInfo fields (previous, pliesFromNull, capturedPiece) that cannot
	// be deduced from a fen string, so set() clears them and to not lose the info
	// we need to backup and later restore setupStates->back(). Note that setupStates
	// is shared by threads but is accessed in read-only mode.
	StateInfo tmp = setupStates->back();

	for (Thread* th : *this)
	{
		th->nodes = th->nmp_ply = th->nmp_odd = 0;
		th->rootDepth = th->completedDepth = DEPTH_ZERO;
		th->rootMoves = rootMoves;
		th->rootPos.set(pos.sfen(), &setupStates->back(), th);
	}

	setupStates->back() = tmp;

	main()->start_searching();
}