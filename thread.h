#ifndef THREAD_H_INCLUDED
#define THREAD_H_INCLUDED

#include <atomic>
#include <bitset>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

#include <omp.h>

#include "movepick.h"
#include "position.h"
#include "search.h"
#include "thread_win32.h"

class ScopedOmpLimiter
{
public:
	ScopedOmpLimiter(int limit) {
		m_originalLimit = omp_get_max_threads();

		if (limit < m_originalLimit)
			omp_set_num_threads(limit);
	}

	~ScopedOmpLimiter() { omp_set_num_threads(m_originalLimit); }

	ScopedOmpLimiter &operator=(const ScopedOmpLimiter &other) = delete;
	ScopedOmpLimiter(const ScopedOmpLimiter &other) = delete;

private:
	int m_originalLimit;
};

class ANNEvaluator;

class Thread {

	Mutex mutex;
	ConditionVariable cv;
	size_t idx;
	bool exit = false, searching = true;
	std::thread stdThread;

public:
	explicit Thread(size_t);
	virtual ~Thread();
	virtual void search();
	void clear();
	void idle_loop();
	void start_searching();
	void wait_for_search_finished();

	size_t PVIdx;
	int selDepth, nmp_ply, nmp_odd;
	std::atomic<uint64_t> nodes;

	Position rootPos;
	Search::RootMoves rootMoves;
	Depth rootDepth, completedDepth;
	CounterMoveHistory counterMoves;
	ButterflyHistory mainHistory;
	ContinuationHistory contHistory;
	ANNEvaluator* annEvaluator;
};

struct MainThread : public Thread {

	using Thread::Thread;

	void search() override;
	void check_time();

	bool easyMovePlayed, failedLow;
	double bestMoveChanges;
	Value previousScore;
	int callsCnt;
};

struct ThreadPool : public std::vector<Thread*> {

	void start_thinking(Position&, StateListPtr&, const Search::LimitsType&, bool = false);
	void clear();
	void set(size_t);

	MainThread* main() { return static_cast<MainThread*>(front()); }
	uint64_t nodes_searched() const { return accumulate(&Thread::nodes); }

	std::atomic_bool stop, ponder, stopOnPonderhit;

private:
	StateListPtr setupStates;

	uint64_t accumulate(std::atomic<uint64_t> Thread::* member) const {

		uint64_t sum = 0;
		for (Thread* th : *this)
			sum += (th->*member).load(std::memory_order_relaxed);
		return sum;
	}
};

extern ThreadPool Threads;

#endif // ifndef THREAD_H_INCLUDED