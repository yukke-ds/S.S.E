#ifndef MISC_H_INCLUDED
#define MISC_H_INCLUDED

#include <cassert>
#include <chrono>
#include <mutex>
#include <ostream>
#include <random>
#include <string>
#include <vector>

#include "types.h"

const std::string engine_info();
void prefetch(void* addr);
void start_logger(bool b);

typedef std::chrono::milliseconds::rep TimePoint; // A value in milliseconds

inline TimePoint now() {
	return std::chrono::duration_cast<std::chrono::milliseconds>
		(std::chrono::steady_clock::now().time_since_epoch()).count();
}

enum SyncCout { IO_LOCK, IO_UNLOCK };
std::ostream& operator<<(std::ostream&, SyncCout);

#define sync_cout std::cout << IO_LOCK
#define sync_endl std::endl << IO_UNLOCK

class PRNG { // Pseudo Random Number Generator

	std::mutex m_mutex;
	std::random_device m_rd;
	uint64_t s;

	uint64_t rand64() {

		s ^= s >> 12, s ^= s << 25, s ^= s >> 27;
		return s * 2685821657736338717LL;
	}

public:
	PRNG() {}
	PRNG(uint64_t seed) : s(seed) { assert(seed); }

	template<typename T> T rand() { return T(rand64()); }

	// Special generator used to fast init magic numbers.
	// Output values only have 1/8th of their bits set on average.
	template<typename T> T sparse_rand()
	{
		return T(rand64() & rand64() & rand64());
	}

	std::random_device::result_type operator()()
	{
		std::lock_guard<std::mutex> l(m_mutex);

		return m_rd();
	}

	std::mt19937 make_mt()
	{
		std::lock_guard<std::mutex> l(m_mutex);

		return std::mt19937(m_rd());
	}
};

extern PRNG gRd;

// Under Windows it is not possible for a process to run on more than one
// logical processor group. This usually means to be limited to use max 64
// cores. To overcome this, some special platform specific API should be
// called to set group affinity for each thread. Original code from Texel by
// Peter Osterlund.
namespace WinProcGroup {
	void bindThisThread(size_t idx);
}

#endif // ifndef MISC_H_INCLUDED