#ifdef _WIN32
#if _WIN32_WINNT < 0x0601
#undef  _WIN32_WINNT
#define _WIN32_WINNT 0x0601 // Force to include needed API prototypes
#endif
#include <windows.h>
// The needed Windows API for processor groups could be missed from old Windows
// versions, so instead of calling them directly (forcing the linker to resolve
// the calls at compile time), try to load them at runtime. To do this we need
// first to define the corresponding function pointers.
extern "C" {
	typedef bool(*fun1_t)(LOGICAL_PROCESSOR_RELATIONSHIP,
		PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX, PDWORD);
	typedef bool(*fun2_t)(USHORT, PGROUP_AFFINITY);
	typedef bool(*fun3_t)(HANDLE, CONST GROUP_AFFINITY*, PGROUP_AFFINITY);
}
#endif

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "misc.h"
#include "thread.h"

using namespace std;

PRNG gRd;

namespace {

// Version number. If Version is left empty, then compile date in the format
// DD-MM-YY and show in engine_info.
const string Version = "4.0";

// Our fancy logging facility. The trick here is to replace cin.rdbuf() and
// cout.rdbuf() with two Tie objects that tie cin and cout to a file stream. We
// can toggle the logging of std::cout and std:cin at runtime whilst preserving
// usual I/O functionality, all without changing a single line of code!
// Idea from http://groups.google.com/group/comp.lang.c++/msg/1d941c0f26ea0d81

struct Tie : public streambuf { // MSVC requires split streambuf for cin and cout

	Tie(streambuf* b, streambuf* l) : buf(b), logBuf(l) {}

	int sync() { return logBuf->pubsync(), buf->pubsync(); }
	int overflow(int c) { return log(buf->sputc((char)c), "<< "); }
	int underflow() { return buf->sgetc(); }
	int uflow() { return log(buf->sbumpc(), ">> "); }

	streambuf *buf, *logBuf;

	int log(int c, const char* prefix) {

		static int last = '\n'; // Single log file

		if (last == '\n')
			logBuf->sputn(prefix, 3);

		return last = logBuf->sputc((char)c);
	}
};

class Logger {

	Logger() : in(cin.rdbuf(), file.rdbuf()), out(cout.rdbuf(), file.rdbuf()) {}
	~Logger() { start(false); }

	ofstream file;
	Tie in, out;

public:

	static void start(bool b) {
		static Logger l;

		if (b && !l.file.is_open())
		{
			l.file.open("io_log.txt", ifstream::out);
			cin.rdbuf(&l.in);
			cout.rdbuf(&l.out);
		}
		else if (!b && l.file.is_open())
		{
			cout.rdbuf(l.out.buf);
			cin.rdbuf(l.in.buf);
			l.file.close();
		}
	}
};

} // namespace

// engine_info() returns the full name of the current S.S.E. version. 
const string engine_info() {

	stringstream ss;

	ss << "S.S.E._DNN_ELMO_TIME_PROGRESS_TEST " << Version << setfill('0')
		<< (Is64Bit ? " 64" : "")
		<< (HasPext ? " BMI2" : (HasPopCnt ? " POPCNT" : ""))
		<< "\nid author by Wada Yusuke" << endl;

	return ss.str();
}

// Used to serialize access to std::cout to avoid multiple threads writing at
// the same time.
std::ostream& operator<<(std::ostream& os, SyncCout sc) {

	static Mutex m;

	if (sc == IO_LOCK)
		m.lock();

	if (sc == IO_UNLOCK)
		m.unlock();

	return os;
}

// Trampoline helper to avoid moving Logger to misc.h
void start_logger(bool b) { Logger::start(b); }


// prefetch() preloads the given address in L1/L2 cache. This is a non-blocking
// function that doesn't stall the CPU waiting for data to be loaded from memory,
// which can be quite slow.
#ifdef NO_PREFETCH

void prefetch(void*) {}

#else

void prefetch(void* addr) {

#  if defined(__INTEL_COMPILER)
	// This hack prevents prefetches from being optimized away by
	// Intel compiler. Both MSVC and gcc seem not be affected by this.
	__asm__("");
#  endif

#  if defined(__INTEL_COMPILER) || defined(_MSC_VER)
	_mm_prefetch((char*)addr, _MM_HINT_T0);
#  else
	__builtin_prefetch(addr);
#  endif
}

#endif

namespace WinProcGroup {

#ifndef _WIN32

void bindThisThread(size_t) {}

#else

/// get_group() retrieves logical processor information using Windows specific
/// API and returns the best group id for the thread with index idx. Original
/// code from Texel by Peter Osterlund.

int get_group(size_t idx) {

	int threads = 0;
	int nodes = 0;
	int cores = 0;
	DWORD returnLength = 0;
	DWORD byteOffset = 0;

	// Early exit if the needed API is not available at runtime
	HMODULE k32 = GetModuleHandle(L"Kernel32.dll");
	auto fun1 = (fun1_t)GetProcAddress(k32, "GetLogicalProcessorInformationEx");
	if (!fun1)
		return -1;

	// First call to get returnLength. We expect it to fail due to null buffer
	if (fun1(RelationAll, nullptr, &returnLength))
		return -1;

	// Once we know returnLength, allocate the buffer
	SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX *buffer, *ptr;
	ptr = buffer = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*)malloc(returnLength);

	// Second call, now we expect to succeed
	if (!fun1(RelationAll, buffer, &returnLength))
	{
		free(buffer);
		return -1;
	}

	while (ptr->Size > 0 && byteOffset + ptr->Size <= returnLength)
	{
		if (ptr->Relationship == RelationNumaNode)
			nodes++;

		else if (ptr->Relationship == RelationProcessorCore)
		{
			cores++;
			threads += (ptr->Processor.Flags == LTP_PC_SMT) ? 2 : 1;
		}

		byteOffset += ptr->Size;
		ptr = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*)(((char*)ptr) + ptr->Size);
	}

	free(buffer);

	std::vector<int> groups;

	// Run as many threads as possible on the same node until core limit is
	// reached, then move on filling the next node.
	for (int n = 0; n < nodes; n++)
		for (int i = 0; i < cores / nodes; i++)
			groups.push_back(n);

	// In case a core has more than one logical processor (we assume 2) and we
	// have still threads to allocate, then spread them evenly across available
	// nodes.
	for (int t = 0; t < threads - cores; t++)
		groups.push_back(t % nodes);

	// If we still have more threads than the total number of logical processors
	// then return -1 and let the OS to decide what to do.
	return idx < groups.size() ? groups[idx] : -1;
}


/// bindThisThread() set the group affinity of the current thread

void bindThisThread(size_t idx) {

	// Use only local variables to be thread-safe
	int group = get_group(idx);

	if (group == -1)
		return;

	// Early exit if the needed API are not available at runtime
	HMODULE k32 = GetModuleHandle(L"Kernel32.dll");
	auto fun2 = (fun2_t)GetProcAddress(k32, "GetNumaNodeProcessorMaskEx");
	auto fun3 = (fun3_t)GetProcAddress(k32, "SetThreadGroupAffinity");

	if (!fun2 || !fun3)
		return;

	GROUP_AFFINITY affinity;
	if (fun2(group, &affinity))
		fun3(GetCurrentThread(), &affinity, nullptr);
}

#endif

} // namespace WinProcGroup