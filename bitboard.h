#ifndef BITBOARD_H_INCLUDED
#define BITBOARD_H_INCLUDED

#include <string>

#include "types.h"

namespace Bitboards {

void init();
}

// 各マスのrookが利きを調べる必要があるマスの数
const int RookBlockBits[SQUARE_NB] = {
	14, 13, 13, 13, 13, 13, 13, 13, 14,
	13, 12, 12, 12, 12, 12, 12, 12, 13,
	13, 12, 12, 12, 12, 12, 12, 12, 13,
	13, 12, 12, 12, 12, 12, 12, 12, 13,
	13, 12, 12, 12, 12, 12, 12, 12, 13,
	13, 12, 12, 12, 12, 12, 12, 12, 13,
	13, 12, 12, 12, 12, 12, 12, 12, 13,
	13, 12, 12, 12, 12, 12, 12, 12, 13,
	14, 13, 13, 13, 13, 13, 13, 13, 14
};

// 各マスのbishopが利きを調べる必要があるマスの数
const int BishopBlockBits[SQUARE_NB] = {
	7,  6,  6,  6,  6,  6,  6,  6,  7,
	6,  6,  6,  6,  6,  6,  6,  6,  6,
	6,  6,  8,  8,  8,  8,  8,  6,  6,
	6,  6,  8, 10, 10, 10,  8,  6,  6,
	6,  6,  8, 10, 12, 10,  8,  6,  6,
	6,  6,  8, 10, 10, 10,  8,  6,  6,
	6,  6,  8,  8,  8,  8,  8,  6,  6,
	6,  6,  6,  6,  6,  6,  6,  6,  6,
	7,  6,  6,  6,  6,  6,  6,  6,  7
};

// Magic Bitboard で利きを求める際のシフト量
// PEXT Bitboardを使用する際はシフト量を減らす必要が無い
const int RookShiftBits[SQUARE_NB] = {
	50, 51, 51, 51, 51, 51, 51, 51, 50,
#if defined(USE_PEXT)
	51, 52, 52, 52, 52, 52, 52, 52, 51,
#else
	51, 52, 52, 52, 52, 52, 52, 52, 50, // [17]: 51 -> 50
#endif
	51, 52, 52, 52, 52, 52, 52, 52, 51,
	51, 52, 52, 52, 52, 52, 52, 52, 51,
	51, 52, 52, 52, 52, 52, 52, 52, 51,
#if defined(USE_PEXT)
	51, 52, 52, 52, 52, 52, 52, 52, 51,
#else
	51, 52, 52, 52, 52, 52, 52, 52, 50, // [53]: 51 -> 50
#endif
	51, 52, 52, 52, 52, 52, 52, 52, 51,
	51, 52, 52, 52, 52, 52, 52, 52, 51,
	50, 51, 51, 51, 51, 51, 51, 51, 50
};

// Magic Bitboard で利きを求める際のシフト量
const int BishopShiftBits[SQUARE_NB] = {
	57, 58, 58, 58, 58, 58, 58, 58, 57,
	58, 58, 58, 58, 58, 58, 58, 58, 58,
	58, 58, 56, 56, 56, 56, 56, 58, 58,
	58, 58, 56, 54, 54, 54, 56, 58, 58,
	58, 58, 56, 54, 52, 54, 56, 58, 58,
	58, 58, 56, 54, 54, 54, 56, 58, 58,
	58, 58, 56, 56, 56, 56, 56, 58, 58,
	58, 58, 58, 58, 58, 58, 58, 58, 58,
	57, 58, 58, 58, 58, 58, 58, 58, 57
};

// index を求める為に使用する。
const int Slide[SQUARE_NB] = {
	1,  1,  1,  1,  1,  1,  1,  1,  1,
	10, 10, 10, 10, 10, 10, 10, 10, 10,
	19, 19, 19, 19, 19, 19, 19, 19, 19,
	28, 28, 28, 28, 28, 28, 28, 28, 28,
	37, 37, 37, 37, 37, 37, 37, 37, 37,
	46, 46, 46, 46, 46, 46, 46, 46, 46,
	55, 55, 55, 55, 55, 55, 55, 55, 55,
	1,  1,  1,  1,  1,  1,  1,  1,  1,
	10, 10, 10, 10, 10, 10, 10, 10, 10
};

#if defined(USE_PEXT)
#else
const uint64_t RookMagic[SQUARE_NB] = {
	UINT64_C(0x140000400809300),  UINT64_C(0x1320000902000240), UINT64_C(0x8001910c008180),
	UINT64_C(0x40020004401040),   UINT64_C(0x40010000d01120),   UINT64_C(0x80048020084050),
	UINT64_C(0x40004000080228),   UINT64_C(0x400440000a2a0a),   UINT64_C(0x40003101010102),
	UINT64_C(0x80c4200012108100), UINT64_C(0x4010c00204000c01), UINT64_C(0x220400103250002),
	UINT64_C(0x2600200004001),    UINT64_C(0x40200052400020),   UINT64_C(0xc00100020020008),
	UINT64_C(0x9080201000200004), UINT64_C(0x2200201000080004), UINT64_C(0x80804c0020200191),
	UINT64_C(0x45383000009100),   UINT64_C(0x30002800020040),   UINT64_C(0x40104000988084),
	UINT64_C(0x108001000800415),  UINT64_C(0x14005000400009),   UINT64_C(0xd21001001c00045),
	UINT64_C(0xc0003000200024),   UINT64_C(0x40003000280004),   UINT64_C(0x40021000091102),
	UINT64_C(0x2008a20408000d00), UINT64_C(0x2000100084010040), UINT64_C(0x144080008008001),
	UINT64_C(0x50102400100026a2), UINT64_C(0x1040020008001010), UINT64_C(0x1200200028005010),
	UINT64_C(0x4280030030020898), UINT64_C(0x480081410011004),  UINT64_C(0x34000040800110a),
	UINT64_C(0x101000010c0021),   UINT64_C(0x9210800080082),    UINT64_C(0x6100002000400a7),
	UINT64_C(0xa2240800900800c0), UINT64_C(0x9220082001000801), UINT64_C(0x1040008001140030),
	UINT64_C(0x40002220040008),   UINT64_C(0x28000124008010c),  UINT64_C(0x40008404940002),
	UINT64_C(0x40040800010200),   UINT64_C(0x90000809002100),   UINT64_C(0x2800080001000201),
	UINT64_C(0x1400020001000201), UINT64_C(0x180081014018004),  UINT64_C(0x1100008000400201),
	UINT64_C(0x80004000200201),   UINT64_C(0x420800010000201),  UINT64_C(0x2841c00080200209),
	UINT64_C(0x120002401040001),  UINT64_C(0x14510000101000b),  UINT64_C(0x40080000808001),
	UINT64_C(0x834000188048001),  UINT64_C(0x4001210000800205), UINT64_C(0x4889a8007400201),
	UINT64_C(0x2080044080200062), UINT64_C(0x80004002861002),   UINT64_C(0xc00842049024),
	UINT64_C(0x8040000202020011), UINT64_C(0x400404002c0100),   UINT64_C(0x2080028202000102),
	UINT64_C(0x8100040800590224), UINT64_C(0x2040009004800010), UINT64_C(0x40045000400408),
	UINT64_C(0x2200240020802008), UINT64_C(0x4080042002200204), UINT64_C(0x4000b0000a00a2),
	UINT64_C(0xa600000810100),    UINT64_C(0x1410000d001180),   UINT64_C(0x2200101001080),
	UINT64_C(0x100020014104e120), UINT64_C(0x2407200100004810), UINT64_C(0x80144000a0845050),
	UINT64_C(0x1000200060030c18), UINT64_C(0x4004200020010102), UINT64_C(0x140600021010302)
};

const uint64_t BishopMagic[SQUARE_NB] = {
	UINT64_C(0x20101042c8200428), UINT64_C(0x840240380102),     UINT64_C(0x800800c018108251),
	UINT64_C(0x82428010301000),   UINT64_C(0x481008201000040),  UINT64_C(0x8081020420880800),
	UINT64_C(0x804222110000),     UINT64_C(0xe28301400850),     UINT64_C(0x2010221420800810),
	UINT64_C(0x2600010028801824), UINT64_C(0x8048102102002),    UINT64_C(0x4000248100240402),
	UINT64_C(0x49200200428a2108), UINT64_C(0x460904020844),     UINT64_C(0x2001401020830200),
	UINT64_C(0x1009008120),       UINT64_C(0x4804064008208004), UINT64_C(0x4406000240300ca0),
	UINT64_C(0x222001400803220),  UINT64_C(0x226068400182094),  UINT64_C(0x95208402010d0104),
	UINT64_C(0x4000807500108102), UINT64_C(0xc000200080500500), UINT64_C(0x5211000304038020),
	UINT64_C(0x1108100180400820), UINT64_C(0x10001280a8a21040), UINT64_C(0x100004809408a210),
	UINT64_C(0x202300002041112),  UINT64_C(0x4040a8000460408),  UINT64_C(0x204020021040201),
	UINT64_C(0x8120013180404),    UINT64_C(0xa28400800d020104), UINT64_C(0x200c201000604080),
	UINT64_C(0x1082004000109408), UINT64_C(0x100021c00c410408), UINT64_C(0x880820905004c801),
	UINT64_C(0x1054064080004120), UINT64_C(0x30c0a0224001030),  UINT64_C(0x300060100040821),
	UINT64_C(0x51200801020c006),  UINT64_C(0x2100040042802801), UINT64_C(0x481000820401002),
	UINT64_C(0x40408a0450000801), UINT64_C(0x810104200000a2),   UINT64_C(0x281102102108408),
	UINT64_C(0x804020040280021),  UINT64_C(0x2420401200220040), UINT64_C(0x80010144080c402),
	UINT64_C(0x80104400800002),   UINT64_C(0x1009048080400081), UINT64_C(0x100082000201008c),
	UINT64_C(0x10001008080009),   UINT64_C(0x2a5006b80080004),  UINT64_C(0xc6288018200c2884),
	UINT64_C(0x108100104200a000), UINT64_C(0x141002030814048),  UINT64_C(0x200204080010808),
	UINT64_C(0x200004013922002),  UINT64_C(0x2200000020050815), UINT64_C(0x2011010400040800),
	UINT64_C(0x1020040004220200), UINT64_C(0x944020104840081),  UINT64_C(0x6080a080801c044a),
	UINT64_C(0x2088400811008020), UINT64_C(0xc40aa04208070),    UINT64_C(0x4100800440900220),
	UINT64_C(0x48112050),         UINT64_C(0x818200d062012a10), UINT64_C(0x402008404508302),
	UINT64_C(0x100020101002),     UINT64_C(0x20040420504912),   UINT64_C(0x2004008118814),
	UINT64_C(0x1000810650084024), UINT64_C(0x1002a03002408804), UINT64_C(0x2104294801181420),
	UINT64_C(0x841080240500812),  UINT64_C(0x4406009000004884), UINT64_C(0x80082004012412),
	UINT64_C(0x80090880808183),   UINT64_C(0x300120020400410),  UINT64_C(0x21a090100822002)
};
#endif

struct Bitboard {

	union {
		uint64_t p[2];
		//     p[1]                 p[0]
		//  |  8  7  6  5  4  3  2  1  0
		//ーーーーーーーーーーーーーーーー
		// 0|  9  0 54 45 36 27 18  9  0
		// 1| 10  1 55 46 37 28 19 10  1
		// 2| 11  2 56 47 38 29 20 11  2
		// 3| 12  3 57 48 39 30 21 12  3
		// 4| 13  4 58 49 40 31 22 13  4
		// 5| 14  5 59 50 41 32 23 14  5
		// 6| 15  6 60 51 42 33 24 15  6
		// 7| 16  7 61 52 43 34 25 16  7
		// 8| 17  8 62 53 44 35 26 17  8

		__m128i m;
	};

	Bitboard() {}
	Bitboard(const Bitboard& b) { _mm_store_si128(&this->m, b.m); }
	Bitboard(uint64_t p0, uint64_t p1) { p[0] = p0; p[1] = p1; }
	Bitboard(Square s);
	
	// Bitboard state
	void set(uint64_t p0, uint64_t p1);
	uint64_t merge() const;
	uint64_t cross() const;
	constexpr static int part(const Square s);

	// Extract ones from bitboard
	Square pop();
	Square pop_c() const;
	int pop_count() const;

	// Accessing elements
	template<typename Consumer> void for_each(Consumer) const;

	// Operators
	operator bool() const;
	bool operator == (const Bitboard& b) const;
	bool operator != (const Bitboard& b) const;
	Bitboard& operator = (const Bitboard& b);
	Bitboard& operator |= (const Bitboard& b);
	Bitboard& operator &= (const Bitboard& b);
	Bitboard& operator ^= (const Bitboard& b);
	Bitboard& operator += (const Bitboard& b);
	Bitboard& operator -= (const Bitboard& b);
	Bitboard& operator <<= (int shift);
	Bitboard& operator >>= (int shift);
	Bitboard operator & (const Bitboard& rhs) const { return Bitboard(*this) &= rhs; }
	Bitboard operator | (const Bitboard& rhs) const { return Bitboard(*this) |= rhs; }
	Bitboard operator ^ (const Bitboard& rhs) const { return Bitboard(*this) ^= rhs; }
	Bitboard operator + (const Bitboard& rhs) const { return Bitboard(*this) += rhs; }
	Bitboard operator << (const int i) const { return Bitboard(*this) <<= i; }
	Bitboard operator >> (const int i) const { return Bitboard(*this) >>= i; }
};

const Bitboard ALL1BB = Bitboard(UINT64_C(0x7fffffffffffffff), UINT64_C(0x3ffff));
const Bitboard ALL0BB = Bitboard(0, 0);

const Bitboard File1BB = Bitboard(UINT64_C(0x1ff) << (9 * 0), 0);
const Bitboard File2BB = Bitboard(UINT64_C(0x1ff) << (9 * 1), 0);
const Bitboard File3BB = Bitboard(UINT64_C(0x1ff) << (9 * 2), 0);
const Bitboard File4BB = Bitboard(UINT64_C(0x1ff) << (9 * 3), 0);
const Bitboard File5BB = Bitboard(UINT64_C(0x1ff) << (9 * 4), 0);
const Bitboard File6BB = Bitboard(UINT64_C(0x1ff) << (9 * 5), 0);
const Bitboard File7BB = Bitboard(UINT64_C(0x1ff) << (9 * 6), 0);
const Bitboard File8BB = Bitboard(0, 0x1ff << (9 * 0));
const Bitboard File9BB = Bitboard(0, 0x1ff << (9 * 1));

const Bitboard Rank1BB = Bitboard(UINT64_C(0x40201008040201) << 0, 0x201 << 0);
const Bitboard Rank2BB = Bitboard(UINT64_C(0x40201008040201) << 1, 0x201 << 1);
const Bitboard Rank3BB = Bitboard(UINT64_C(0x40201008040201) << 2, 0x201 << 2);
const Bitboard Rank4BB = Bitboard(UINT64_C(0x40201008040201) << 3, 0x201 << 3);
const Bitboard Rank5BB = Bitboard(UINT64_C(0x40201008040201) << 4, 0x201 << 4);
const Bitboard Rank6BB = Bitboard(UINT64_C(0x40201008040201) << 5, 0x201 << 5);
const Bitboard Rank7BB = Bitboard(UINT64_C(0x40201008040201) << 6, 0x201 << 6);
const Bitboard Rank8BB = Bitboard(UINT64_C(0x40201008040201) << 7, 0x201 << 7);
const Bitboard Rank9BB = Bitboard(UINT64_C(0x40201008040201) << 8, 0x201 << 8);

const Bitboard Rank1_3 = (Rank1BB | Rank2BB | Rank3BB);
const Bitboard Rank1_6 = (Rank1_3 | Rank4BB | Rank5BB | Rank6BB);
const Bitboard Rank1_7 = (Rank1_6 | Rank7BB);
const Bitboard Rank7_9 = (Rank7BB | Rank8BB | Rank9BB);
const Bitboard Rank4_9 = (Rank4BB | Rank5BB | Rank6BB | Rank7_9);
const Bitboard Rank3_9 = (Rank3BB | Rank4_9);

extern Bitboard SquareBB[SQUARE_NB];
extern Bitboard FileBB[FILE_NB];
extern Bitboard RankBB[RANK_NB];

#if defined(USE_PEXT)
extern Bitboard RookAttack[495616];
#else
extern const uint64_t RookMagic[SQUARE_NB];
extern const uint64_t BishopMagic[SQUARE_NB];
extern Bitboard RookAttack[512000];
#endif
extern int RookAttackIndex[SQUARE_NB];
extern Bitboard RookBlockMask[SQUARE_NB];
extern Bitboard BishopAttack[20224];
extern int BishopAttackIndex[SQUARE_NB];
extern Bitboard BishopBlockMask[SQUARE_NB];
extern Bitboard LanceAttack[COLOR_NB][SQUARE_NB][128];
extern Bitboard KingAttack[SQUARE_NB];
extern Bitboard GoldAttack[COLOR_NB][SQUARE_NB];
extern Bitboard SilverAttack[COLOR_NB][SQUARE_NB];
extern Bitboard KnightAttack[COLOR_NB][SQUARE_NB];
extern Bitboard PawnAttack[COLOR_NB][SQUARE_NB];

extern Bitboard InFrontBB[COLOR_NB][RANK_NB];
extern Bitboard LineBB[SQUARE_NB][SQUARE_NB];
extern Bitboard BetweenBB[SQUARE_NB][SQUARE_NB];

extern Bitboard CheckerCandidate[PIECE_NB][SQUARE_NB];
extern Bitboard AdjacentCheckCandidate[PIECE_NB][SQUARE_NB];
extern Direction SquareRelation[SQUARE_NB][SQUARE_NB];

inline Bitboard::Bitboard(Square s) {
	*this = SquareBB[s];
}

inline void Bitboard::set(uint64_t p0, uint64_t p1) {
	p[0] = p0; p[1] = p1;
}

inline uint64_t Bitboard::merge() const {
	return p[0] | p[1]; 
}

inline uint64_t Bitboard::cross() const {
	return p[0] & p[1];
}

inline constexpr int Bitboard::part(const Square s) {
	return static_cast<int>(SQ_79 < s);
}

inline Square Bitboard::pop() {
	return (p[0] != 0) ? Square(pop_lsb(p[0])) : Square(pop_lsb(p[1]) + 63); 
}

inline Square Bitboard::pop_c() const {
	return (p[0] != 0) ? Square(LSB64(p[0])) : Square(LSB64(p[1]) + 63); 
}

inline int Bitboard::pop_count() const {
	return (int)(POPCNT64(p[0]) + POPCNT64(p[1])); 
}

template<typename Consumer>
inline void Bitboard::for_each(Consumer function) const {
	Bitboard clone = *this;
	while (clone) {
		Square s = clone.pop();
		function(s);
	}
}

inline Bitboard::operator bool() const {
	return !(_mm_testz_si128(m, _mm_set1_epi8(static_cast<char>(0xffu)))); 
}

inline bool Bitboard::operator == (const Bitboard& b) const {
	__m128i neq = _mm_xor_si128(this->m, b.m);
	return _mm_test_all_zeros(neq, neq) ? true : false;
}

inline bool Bitboard::operator != (const Bitboard& b) const { 
	return !(*this == b); 
}

inline Bitboard& Bitboard::operator = (const Bitboard& b) {
	_mm_store_si128(&this->m, b.m); 
	return *this; 
}

inline Bitboard& Bitboard::operator |= (const Bitboard& b) {
	this->m = _mm_or_si128(m, b.m);
	return *this;
}

inline Bitboard& Bitboard::operator &= (const Bitboard& b) { 
	this->m = _mm_and_si128(m, b.m); 
	return *this; 
}

inline Bitboard& Bitboard::operator ^= (const Bitboard& b) {
	this->m = _mm_xor_si128(m, b.m); 
	return *this;
}

inline Bitboard& Bitboard::operator += (const Bitboard& b) { 
	this->m = _mm_add_epi64(m, b.m); 
	return *this; 
}

inline Bitboard& Bitboard::operator -= (const Bitboard& b) {
	this->m = _mm_sub_epi64(m, b.m); 
	return *this;
}

inline Bitboard& Bitboard::operator <<= (int shift) {
	m = _mm_slli_epi64(m, shift); 
	return *this; 
}

inline Bitboard& Bitboard::operator >>= (int shift) { 
	m = _mm_srli_epi64(m, shift); 
	return *this;
}

/* Other helpers */
inline Bitboard operator|(const Bitboard& b, Square s) {
	return b | SquareBB[s];
}

inline Bitboard operator&(const Bitboard& b, Square s) {
	return b & SquareBB[s];
}

inline Bitboard operator^(const Bitboard& b, Square s) {
	return b ^ SquareBB[s];
}

inline Bitboard operator~(const Bitboard& b) {
	Bitboard t;
	t.m = _mm_xor_si128(b.m, ALL1BB.m);
	return t;
}

inline std::ostream& operator<<(std::ostream& os, const Bitboard& board) {
	// Bitboardを表示する(デバッグ用)
	for (Rank rank = RANK_1; rank <= RANK_9; ++rank) {
		for (File file = FILE_9; file >= FILE_1; --file)
			os << ((board & (file | rank)) ? " *" : " .");
		os << std::endl;
	}
	os << std::endl;

	return os;
}

inline Bitboard file_fill(Bitboard b) {
	// Kogge-Stone Algorithm
	// (http://chessprogramming.wikispaces.com/Kogge-Stone+Algorithm)
	b.m = _mm_or_si128(b.m, _mm_srli_epi64(b.m, 1));
	b.m = _mm_or_si128(b.m, _mm_srli_epi64(b.m, 2));
	b.m = _mm_or_si128(b.m, _mm_srli_epi64(b.m, 4));
	b.m = _mm_or_si128(b.m, _mm_srli_epi64(b.m, 1));
	b &= Rank1BB;
	b.m = _mm_or_si128(b.m, _mm_slli_epi64(b.m, 1));
	b.m = _mm_or_si128(b.m, _mm_slli_epi64(b.m, 2));
	b.m = _mm_or_si128(b.m, _mm_slli_epi64(b.m, 4));
	b.m = _mm_or_si128(b.m, _mm_slli_epi64(b.m, 1));
	return b & ALL1BB;
}

inline Direction square_relation(const Square s1, const Square s2) {
	return SquareRelation[s1][s2]; 
}

#if defined(USE_PEXT)
inline uint64_t occupied_to_index(const Bitboard& occupied, const Bitboard& mask) { return pext(occupied.merge(), mask.merge()); }
inline Bitboard rookAttack(const Square sq, const Bitboard& occupied) {
	const Bitboard block(occupied & RookBlockMask[sq]);
	return RookAttack[RookAttackIndex[sq] + occupied_to_index(block, RookBlockMask[sq])];
}
inline Bitboard bishopAttack(const Square sq, const Bitboard& occupied) {
	const Bitboard block(occupied & BishopBlockMask[sq]);
	return BishopAttack[BishopAttackIndex[sq] + occupied_to_index(block, BishopBlockMask[sq])];
}
#else
inline uint64_t occupied_to_index(const Bitboard& block, const uint64_t magic, const int shiftBits) {
	return (block.merge() * magic) >> shiftBits;
}

inline Bitboard rook_attacks(const Square s, const Bitboard& occupied) {
	const Bitboard block(occupied & RookBlockMask[s]);
	return RookAttack[RookAttackIndex[s] + occupied_to_index(block, RookMagic[s], RookShiftBits[s])];
}

inline Bitboard bishop_attacks(const Square s, const Bitboard& occupied) {
	const Bitboard block(occupied & BishopBlockMask[s]);
	return BishopAttack[BishopAttackIndex[s] + occupied_to_index(block, BishopMagic[s], BishopShiftBits[s])];
}
#endif

inline Bitboard lance_attacks(const Color c, const Square s, const Bitboard& occupied) {
	const int index = (occupied.p[Bitboard::part(s)] >> Slide[s]) & 127;
	return LanceAttack[c][s][index];
}

inline Bitboard gold_attacks(const Color c, const Square s) { 
	return GoldAttack[c][s]; 
}

inline Bitboard silver_attacks(const Color c, const Square s) { 
	return SilverAttack[c][s]; 
}

inline Bitboard knight_attacks(const Color c, const Square s) {
	return KnightAttack[c][s]; 
}

inline Bitboard pawn_attacks(const Color c, const Square s) {
	return PawnAttack[c][s]; 
}

inline Bitboard king_attacks(const Square s) {
	return KingAttack[s];
}

inline Bitboard dragon_attacks(const Square s, const Bitboard& occupied) {
	return rook_attacks(s, occupied) | king_attacks(s);
}

inline Bitboard horse_attacks(const Square s, const Bitboard& occupied) { 
	return bishop_attacks(s, occupied) | king_attacks(s);
}

inline Bitboard line_bb(const Square s1, const Square s2) {
	return LineBB[s1][s2];
}

inline Bitboard between_bb(const Square s1, const Square s2) {
	return BetweenBB[s1][s2]; 
}

inline Bitboard checker_candidates(const Piece pc, const Square s) {
	return CheckerCandidate[pc][s];
}

inline Bitboard adjacent_check_candidates(const Piece pc, const Square s) {
	return AdjacentCheckCandidate[pc][s];
}

inline Bitboard bishop_step_attacks(const Square s) {
	return silver_attacks(BLACK, s) & silver_attacks(WHITE, s);
}

inline Bitboard attacks_bb(Piece pc, Square s, Bitboard occupied) {

	switch (pc) {
	case B_PAWN: return pawn_attacks(BLACK, s);
	case B_LANCE: return lance_attacks(BLACK, s, occupied);
	case B_KNIGHT: return knight_attacks(BLACK, s);
	case B_SILVER: return silver_attacks(BLACK, s);
	case B_GOLD: case B_PRO_PAWN: case B_PRO_LANCE: case B_PRO_KNIGHT: case B_PRO_SILVER: return gold_attacks(BLACK, s);

	case W_PAWN: return pawn_attacks(WHITE, s);
	case W_LANCE: return lance_attacks(WHITE, s, occupied);
	case W_KNIGHT: return knight_attacks(WHITE, s);
	case W_SILVER: return silver_attacks(WHITE, s);
	case W_GOLD: case W_PRO_PAWN: case W_PRO_LANCE: case W_PRO_KNIGHT: case W_PRO_SILVER: return gold_attacks(WHITE, s);

	case B_BISHOP: case W_BISHOP: return bishop_attacks(s, occupied);
	case B_ROOK:   case W_ROOK:   return rook_attacks(s, occupied);
	case B_HORSE:  case W_HORSE:  return horse_attacks(s, occupied);
	case B_DRAGON: case W_DRAGON: return dragon_attacks(s, occupied);
	case B_KING:   case W_KING:   return king_attacks(s);
	case NO_PIECE: case PIECE_KIND: return ALL0BB;

	default: assert(false); return ALL1BB;
	}
}

inline Bitboard shift_bb(const Bitboard& b, Square s) {
	return s == DELTA_N ? b >> 1 : s == DELTA_S ? b << 1 : ALL0BB;
}

inline Bitboard promotion_area(const Color c) {
	return c == BLACK ? Rank1_3 : Rank7_9;
}

inline bool aligned(Square s1, Square s2, Square s3) {

	Direction d = square_relation(s1, s3);
	return d ? d == square_relation(s2, s3) : false;
}

inline bool more_than_one(const Bitboard& b) {
	assert(!b.cross()); 
	return POPCNT64(b.merge()) > 1;
}

#endif // ifndef BITBOARD_H_INCLUDED