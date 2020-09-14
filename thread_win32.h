#ifndef THREAD_WIN32_H_INCLUDED
#define THREAD_WIN32_H_INCLUDED

#include <condition_variable>
#include <mutex>

#if defined(_WIN32) && !defined(_MSC_VER)

#ifndef NOMINMAX
#  define NOMINMAX // Disable macros min() and max()
#endif

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef WIN32_LEAN_AND_MEAN
#undef NOMINMAX

/// Mutex and ConditionVariable struct are wrappers of the low level locking
/// machinery and are modeled after the corresponding C++11 classes.

struct Mutex {
	Mutex() { InitializeCriticalSection(&cs); }
	~Mutex() { DeleteCriticalSection(&cs); }
	void lock() { EnterCriticalSection(&cs); }
	void unlock() { LeaveCriticalSection(&cs); }

private:
	CRITICAL_SECTION cs;
};

typedef std::condition_variable_any ConditionVariable;

#else // Default case: use STL classes

typedef std::mutex Mutex;
typedef std::condition_variable ConditionVariable;

#endif

#endif // ifndef THREAD_WIN32_H_INCLUDED