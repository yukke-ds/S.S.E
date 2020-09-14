#ifndef TT_H_INCLUDED
#define TT_H_INCLUDED

#include "misc.h"
#include "types.h"

// TTEntry struct is the 10 bytes transposition table entry, defined as below:
//
// key        16 bit
// move       16 bit
// value      16 bit
// eval value 16 bit
// generation  6 bit
// bound type  2 bit
// depth       8 bit

struct TTEntry {

	Move  move()  const { return (Move)move16; }
	Value value() const { return (Value)value16; }
	Value eval()  const { return (Value)eval16; }
	Depth depth() const { return (Depth)depth8; }
	Bound bound() const { return (Bound)(genBound8 & 0x3); }

	void save(Key k, Value v, Bound b, Depth d, Move m, Value ev, uint8_t g) {

		assert(d / ONE_PLY * ONE_PLY == d);

		// Preserve any existing move for the same position
		if (m || (k >> 48) != key16)
			move16 = (uint16_t)m;

		// Don't overwrite more valuable entries
		if (  (k >> 48) != key16
			|| d > depth8 - 4
			/* || g != (genBound8 & 0xFC) // Matching non-zero keys are already refreshed by probe() */
			|| b == BOUND_EXACT)
		{
			key16		= (uint16_t)(k >> 48);
			value16		= (int16_t)v;
			eval16		= (int16_t)ev;
			genBound8	= (uint8_t)(g | b);
			depth8		= (int8_t)d;
		}
	}

private:
	friend class TranspositionTable;

	uint16_t key16;
	uint16_t move16;
	int16_t  value16;
	int16_t  eval16;
	uint8_t  genBound8;
	int8_t   depth8;
};

class TranspositionTable {

	static const int CacheLineSize = 64;
	static const int ClusterSize = 3;

	struct Cluster {
		TTEntry entry[ClusterSize];
		char padding[2]; // Align to a divisor of the cache line size
	};

	static_assert(CacheLineSize % sizeof(Cluster) == 0, "Cluster size incorrect");

public:
	~TranspositionTable() { free(mem); }
	void new_search() { generation8 += 4; } // Lower 2 bits are used by Bound
	uint8_t generation() const { return generation8; }
	TTEntry* probe(const Key key, bool& found) const;
	int hashfull() const;
	void resize(size_t mbSize);
	void clear();

	// The lowest order bits of the key are used to get the index of the cluster
	TTEntry* first_entry(const Key key) const {
		return &table[(size_t)key & (clusterCount - 1)].entry[0];
	}

private:
	size_t clusterCount;
	Cluster* table;
	void* mem;
	uint8_t generation8;
};

extern TranspositionTable TT;

#endif // ifndef TT_H_INCLUDED