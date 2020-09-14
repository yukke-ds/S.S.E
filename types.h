#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED

#include <cassert>
#include <cctype>
#include <climits>
#include <cstdint>
#include <cstdlib>

#if defined(_MSC_VER)
// Disable some silly and noisy warning from MSVC compiler
#pragma warning(disable: 4127) // Conditional expression is constant
#pragma warning(disable: 4146) // Unary minus operator applied to unsigned type
#pragma warning(disable: 4800) // Forcing value to bool 'true' or 'false'
#endif

// Predefined macros hell:
//
// __GNUC__           Compiler is gcc, Clang or Intel on Linux
// __INTEL_COMPILER   Compiler is Intel
// _MSC_VER           Compiler is MSVC or Intel on Windows
// _WIN32             Building on Windows (any)
// _WIN64             Building on Windows 64 bit

#if defined(_WIN64) && defined(_MSC_VER) // No Makefile used
#  include <intrin.h> // MSVC popcnt and bsfq instrinsics
#  define IS_64BIT
#  define USE_BSFQ
#  define USE_POPCNT // SSE4.2
//#  define USE_PEXT // AVX2
#endif

#if defined(USE_POPCNT) && defined(__INTEL_COMPILER) && defined(_MSC_VER)
#  include <nmmintrin.h> // Intel header for _mm_popcnt_u64() intrinsic
#endif

#if !defined(NO_PREFETCH) && (defined(__INTEL_COMPILER) || defined(_MSC_VER))
#  include <xmmintrin.h> // Intel and Microsoft header for _mm_prefetch()
#endif

#if defined(USE_PEXT)
#  include <immintrin.h> // Header for _pext_u64(), __mm256_...() intrinsic
#  define pext(b, m) _pext_u64(b, m)
#else
#  define pext(b, m) (0)
#endif

#ifdef USE_POPCNT
#  define POPCNT32(a) __popcnt32(a)
#  define POPCNT64(a) __popcnt64(a)
const bool HasPopCnt = true;
#else
const bool HasPopCnt = false;
#endif

#ifdef USE_PEXT
const bool HasPext = true;
#else
const bool HasPext = false;
#endif

#ifdef IS_64BIT
const bool Is64Bit = true;
#else
const bool Is64Bit = false;
#endif

#ifdef IS_64BIT
// 1Ç≈Ç†ÇÈç≈â∫à ÇÃbitÇÃbità íuÇìæÇÈÅB0ÇìnÇµÇƒÇÕÇ»ÇÁÇ»Ç¢ÅB
inline int LSB32(uint32_t v) { assert(v != 0); unsigned long index; _BitScanForward(&index, v); return index; }
inline int LSB64(uint64_t v) { assert(v != 0); unsigned long index; _BitScanForward64(&index, v); return index; }

// 1Ç≈Ç†ÇÈç≈è„à ÇÃbitÇÃbità íuÇìæÇÈÅB0ÇìnÇµÇƒÇÕÇ»ÇÁÇ»Ç¢ÅB
inline int MSB32(uint32_t v) { assert(v != 0); unsigned long index; _BitScanReverse(&index, v); return index; }
inline int MSB64(uint64_t v) { assert(v != 0); unsigned long index; _BitScanReverse64(&index, v); return index; }
#else
// 32bitä¬ã´Ç≈ÇÕ64bitî≈ÇóvãÅÇ≥ÇÍÇΩÇÁ2âÒÇ…ï™ÇØÇƒé¿çsÅB
inline int LSB32(uint32_t v) { assert(v != 0); unsigned long index; _BitScanForward(&index, v); return index; }
inline int LSB64(uint64_t v) { assert(v != 0); return uint32_t(v) ? LSB32(uint32_t(v)) : 32 + LSB32(uint32_t(v >> 32)); }

inline int MSB32(uint32_t v) { assert(v != 0); unsigned long index; _BitScanReverse(&index, v); return index; }
inline int MSB64(uint64_t v) { assert(v != 0); return uint32_t(v >> 32) ? 32 + MSB32(uint32_t(v >> 32)) : MSB32(uint32_t(v)); }
#endif

typedef uint64_t Key;

// ç≈ëÂ593éËÅ®position sfen R8/2K1S1SSk/4B4/9/9/9/9/9/1L1L1L3 b PLNSGBR17p3n3g 1
const int MAX_MOVES = 593;
const int MAX_PLY = 128;
const size_t KB = 1024;
const size_t MB = 1024 * KB;

// 0  0  0000000 0000000
// ê¨ ë≈   from     to
// (ãÓë≈Ç»ÇÁfromÇPieceTypeÇ∆ÇµÇƒÇ¢ÇÈ)
enum Move : uint16_t {
	MOVE_NONE,
	MOVE_NULL = 82
};

enum MoveType {
	NORMAL,
	DROP = 1 << 14,
	PROMOTION = 2 << 14
};

enum Color {
	BLACK, WHITE, COLOR_NB = 2
};

enum Bound {
	BOUND_NONE,
	BOUND_UPPER,
	BOUND_LOWER,
	BOUND_EXACT = BOUND_UPPER | BOUND_LOWER
};

enum Value : int {
	VALUE_ZERO			= 0,
	VALUE_DRAW			= 0,
	VALUE_KNOWN_WIN		= 20000,
	VALUE_MATE			= 32000,
	VALUE_INFINITE		= 32001,
	VALUE_NONE			= 32002,
	VALUE_NOT_EVALUATED	= 32003,

	VALUE_MATE_IN_MAX_PLY = int(VALUE_MATE) - MAX_PLY,
	VALUE_MATED_IN_MAX_PLY = -int(VALUE_MATE) + MAX_PLY,

	PawnValue		= 90,
	LanceValue		= 315,
	KnightValue		= 405,
	SilverValue		= 495,
	GoldValue		= 540,
	BishopValue		= 855,
	RookValue		= 990,
	ProPawnValue	= 540,
	ProLanceValue	= 540,
	ProKnightValue	= 540,
	ProSilverValue	= 540,
	HorseValue		= 945,
	DragonValue		= 1395,
	KingValue		= 15000
};

enum PieceType {
	NO_PIECE_TYPE, PAWN, LANCE, KNIGHT, SILVER, BISHOP, ROOK, GOLD, KING,
	PRO_PAWN, PRO_LANCE, PRO_KNIGHT, PRO_SILVER, HORSE, DRAGON, // = 14
	ALL_PIECES = 0,
	PIECE_TYPE_NB = 15
};

enum Piece {
	NO_PIECE = 0,
	B_PAWN = 1, B_LANCE, B_KNIGHT, B_SILVER, B_BISHOP, B_ROOK, B_GOLD, B_KING,
	B_PRO_PAWN, B_PRO_LANCE, B_PRO_KNIGHT, B_PRO_SILVER, B_HORSE, B_DRAGON,
	W_PAWN = 17, W_LANCE, W_KNIGHT, W_SILVER, W_BISHOP, W_ROOK, W_GOLD, W_KING,
	W_PRO_PAWN, W_PRO_LANCE, W_PRO_KNIGHT, W_PRO_SILVER, W_HORSE, W_DRAGON,
	PIECE_NB = 32,
	
	PIECE_PROMOTION = 8,
	PIECE_HAND_NB = 8,
	PIECE_WHITE = 16,
	PIECE_KIND = 16,
	NUM_PIECES = 52
};

enum PieceNumber {
	PIECE_NUMBER_ZERO = 0,

	PAWN_NUMBER = 0, LANCE_NUMBER = 18, KNIGHT_NUMBER = 22, SILVER_NUMBER = 26,
	GOLD_NUMBER = 30, BISHOP_NUMBER = 34, ROOK_NUMBER = 36, KING_NUMBER = 38,
	BKING_NUMBER = 38, WKING_NUMBER = 39,
	PIECE_NUMBER_NB = 40
};

//  0100  0010  0010  0100  0100  0100  00010010
//   ã‡4    îÚ2   äp2   ã‚4   åj4   çÅ4      ï‡18
enum Hand : uint32_t { HAND_ZERO = 0 };

enum Depth {

	ONE_PLY = 1,

	DEPTH_ZERO			= 0,
	DEPTH_QS_CHECKS		= 0,
	DEPTH_QS_NO_CHECKS	= -1,
	DEPTH_QS_RECAPTURES = -5,

	DEPTH_NONE = -6,
	DEPTH_MAX = MAX_PLY
};

enum Square {
	SQ_11, SQ_12, SQ_13, SQ_14, SQ_15, SQ_16, SQ_17, SQ_18, SQ_19, // 72 63 54 45 36 27 18  9 0
	SQ_21, SQ_22, SQ_23, SQ_24, SQ_25, SQ_26, SQ_27, SQ_28, SQ_29, // 73 64 55 46 37 28 19 10 1
	SQ_31, SQ_32, SQ_33, SQ_34, SQ_35, SQ_36, SQ_37, SQ_38, SQ_39, // 74 65 56 47 38 29 20 11 2
	SQ_41, SQ_42, SQ_43, SQ_44, SQ_45, SQ_46, SQ_47, SQ_48, SQ_49, // 75 66 57 48 39 30 21 12 3
	SQ_51, SQ_52, SQ_53, SQ_54, SQ_55, SQ_56, SQ_57, SQ_58, SQ_59, // 76 67 58 49 40 31 22 13 4
	SQ_61, SQ_62, SQ_63, SQ_64, SQ_65, SQ_66, SQ_67, SQ_68, SQ_69, // 77 68 59 50 41 32 23 14 5
	SQ_71, SQ_72, SQ_73, SQ_74, SQ_75, SQ_76, SQ_77, SQ_78, SQ_79, // 78 69 60 51 42 33 24 15 6
	SQ_81, SQ_82, SQ_83, SQ_84, SQ_85, SQ_86, SQ_87, SQ_88, SQ_89, // 79 70 61 52 43 34 25 16 7
	SQ_91, SQ_92, SQ_93, SQ_94, SQ_95, SQ_96, SQ_97, SQ_98, SQ_99, // 80 71 62 53 44 35 26 17 8
	SQ_NONE,

	SQUARE_NB = 81,
	SQUARE_ZERO = 0,

	DELTA_N = -1,
	DELTA_E = -9,
	DELTA_S = +1,
	DELTA_W = +9,

	DELTA_NN = (int)DELTA_N + (int)DELTA_N,
	DELTA_NE = (int)DELTA_N + (int)DELTA_E,
	DELTA_SE = (int)DELTA_S + (int)DELTA_E,
	DELTA_SS = (int)DELTA_S + (int)DELTA_S,
	DELTA_SW = (int)DELTA_S + (int)DELTA_W,
	DELTA_NW = (int)DELTA_N + (int)DELTA_W
};

enum File {
	FILE_1, FILE_2, FILE_3, FILE_4, FILE_5, FILE_6, FILE_7, FILE_8, FILE_9, FILE_NB
};

enum Rank {
	RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8, RANK_9, RANK_NB
};

const File SquareToFile[SQUARE_NB] = {
	FILE_1, FILE_1, FILE_1, FILE_1, FILE_1, FILE_1, FILE_1, FILE_1, FILE_1,
	FILE_2, FILE_2, FILE_2, FILE_2, FILE_2, FILE_2, FILE_2, FILE_2, FILE_2,
	FILE_3, FILE_3, FILE_3, FILE_3, FILE_3, FILE_3, FILE_3, FILE_3, FILE_3,
	FILE_4, FILE_4, FILE_4, FILE_4, FILE_4, FILE_4, FILE_4, FILE_4, FILE_4,
	FILE_5, FILE_5, FILE_5, FILE_5, FILE_5, FILE_5, FILE_5, FILE_5, FILE_5,
	FILE_6, FILE_6, FILE_6, FILE_6, FILE_6, FILE_6, FILE_6, FILE_6, FILE_6,
	FILE_7, FILE_7, FILE_7, FILE_7, FILE_7, FILE_7, FILE_7, FILE_7, FILE_7,
	FILE_8, FILE_8, FILE_8, FILE_8, FILE_8, FILE_8, FILE_8, FILE_8, FILE_8,
	FILE_9, FILE_9, FILE_9, FILE_9, FILE_9, FILE_9, FILE_9, FILE_9, FILE_9
};

const Rank SquareToRank[SQUARE_NB] = {
	RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8, RANK_9,
	RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8, RANK_9,
	RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8, RANK_9,
	RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8, RANK_9,
	RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8, RANK_9,
	RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8, RANK_9,
	RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8, RANK_9,
	RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8, RANK_9,
	RANK_1, RANK_2, RANK_3, RANK_4, RANK_5, RANK_6, RANK_7, RANK_8, RANK_9
};

enum Direction {
	DIREC_MISC = 0, // ècÅAâ°ÅAéŒÇﬂÇÃà íuÇ…ñ≥Ç¢èÍçá
	DIREC_FILE = 2, // èc
	DIREC_RANK = 3, // â°
	DIREC_DIAG_NESW = 4, // âEè„Ç©ÇÁç∂â∫
	DIREC_DIAG_NWSE = 5, // ç∂è„Ç©ÇÁâEâ∫
	DIREC_CROSS = 2, // ècÅAâ°
	DIREC_DIAG = 4 // éŒÇﬂ
};

enum RepetitionState
{
	REPETITION_NONE,     // êÁì˙éËÇ≈ÇÕÇ»Ç¢
	REPETITION_WIN,      // òAë±â§éËÇÃêÁì˙éËÇ…ÇÊÇÈèüÇø
	REPETITION_LOSE,     // òAë±â§éËÇÃêÁì˙éËÇ…ÇÊÇÈïâÇØ
	REPETITION_DRAW,     // òAë±â§éËÇ≈ÇÕÇ»Ç¢ïÅí ÇÃêÁì˙éË
	REPETITION_SUPERIOR, // óDìôã«ñ (î’è„ÇÃãÓÇ™ìØÇ∂Ç≈éËãÓÇ™ëäéËÇÊÇËóDÇÍÇƒÇ¢ÇÈ)
	REPETITION_INFERIOR, // óÚìôã«ñ (î’è„ÇÃãÓÇ™ìØÇ∂Ç≈éËãÓÇ™ëäéËÇÊÇËóDÇÍÇƒÇ¢ÇÈ)
	REPETITION_NB
};

#define ENABLE_BASE_OPERATORS_ON(T)                             \
inline T operator+(T d1, T d2) { return T(int(d1) + int(d2)); } \
inline T operator-(T d1, T d2) { return T(int(d1) - int(d2)); } \
inline T operator*(int i, T d) { return T(i * int(d)); }        \
inline T operator*(T d, int i) { return T(int(d) * i); }        \
inline T operator-(T d) { return T(-int(d)); }                  \
inline T& operator+=(T& d1, T d2) { return d1 = d1 + d2; }      \
inline T& operator-=(T& d1, T d2) { return d1 = d1 - d2; }      \
inline T& operator*=(T& d, int i) { return d = T(int(d) * i); }

#define ENABLE_FULL_OPERATORS_ON(T)                             \
ENABLE_BASE_OPERATORS_ON(T)                                     \
inline T& operator++(T& d) { return d = T(int(d) + 1); }        \
inline T& operator--(T& d) { return d = T(int(d) - 1); }        \
inline T operator++(T& d, int) { T t = d; d = T(int(d) + 1); return t; }   \
inline T operator--(T& d, int) { T t = d; d = T(int(d) - 1); return t; }   \
inline T operator/(T d, int i) { return T(int(d) / i); }        \
inline int operator/(T d1, T d2) { return int(d1) / int(d2); }  \
inline T& operator/=(T& d, int i) { return d = T(int(d) / i); }

ENABLE_FULL_OPERATORS_ON(Value)
ENABLE_FULL_OPERATORS_ON(PieceType)
ENABLE_FULL_OPERATORS_ON(Piece)
ENABLE_FULL_OPERATORS_ON(PieceNumber)
ENABLE_FULL_OPERATORS_ON(Color)
ENABLE_FULL_OPERATORS_ON(Depth)
ENABLE_FULL_OPERATORS_ON(Square)
ENABLE_FULL_OPERATORS_ON(File)
ENABLE_FULL_OPERATORS_ON(Rank)

#undef ENABLE_FULL_OPERATORS_ON
#undef ENABLE_BASE_OPERATORS_ON

// Additional operators to add integers to a Value
constexpr Value operator+(Value v, int i) { return Value(int(v) + i); }
constexpr Value operator-(Value v, int i) { return Value(int(v) - i); }
inline Value& operator+=(Value& v, int i) { return v = v + i; }
inline Value& operator-=(Value& v, int i) { return v = v - i; }

inline Square operator | (File f, Rank r) {
	return Square(f * 9 + r);
}

constexpr Color operator~(Color c) {
	return Color(c ^ WHITE);
}

template <typename T> // í èÌÇÕ32bit
int pop_lsb(T& b) { 
	int index = LSB32(b);  b = T(b & (b - 1)); return index; 
}

inline int pop_lsb(uint64_t & b) {  // 64bitóp
	int index = LSB64(b);  b &= b - 1; return index;
}

constexpr Value mate_in(int ply) {
	return VALUE_MATE - ply;
}

constexpr Value mated_in(int ply) {
	return -VALUE_MATE + ply;
}

constexpr Piece make_piece(Color c, PieceType pt) {
	return Piece((c << 4) | pt);
}

constexpr PieceType type_of(Piece pc) {
	return PieceType(pc & 15);
}

constexpr MoveType type_of(Move m) {
	return MoveType(m & (3 << 14));
}

inline Color color_of(Piece pc) {
	assert(pc != NO_PIECE);
	return Color(pc >> 4);
}

constexpr bool is_ok(File f) {
	return FILE_1 <= f && f <= FILE_9;
}

constexpr bool is_ok(Rank r) {
	return RANK_1 <= r && r <= RANK_9;
}

constexpr bool is_ok(Square s) {
	return SQ_11 <= s && s <= SQ_99;
}

constexpr bool is_ok(PieceNumber pn) {
	return PIECE_NUMBER_ZERO <= pn && pn < PIECE_NUMBER_NB; 
}

constexpr Rank relative_rank(Color c, Rank r) {
	return c == BLACK ? r : RANK_9 - r;
}

constexpr File to_file(char c) {
	return File(c - '1');
}

constexpr Rank to_rank(char c) {
	return Rank(c - 'a');
}

constexpr File file_of(Square s) {
	return SquareToFile[s];
}

constexpr Rank rank_of(Square s) {
	return SquareToRank[s]; 
}

constexpr Square from_sq(Move m) {
	return Square((m >> 7) & 0x7f);
}

constexpr Square to_sq(Move m) {
	return Square(m & 0x7f);
}

constexpr int from_to(Move m) {
	return (int)(from_sq(m) + (type_of(m) == DROP ? (SQUARE_NB - 1) : 0)) * (int)SQUARE_NB + (int)to_sq(m); 
}

constexpr Square inverse(Square s) {
	return Square((SQUARE_NB - 1) - s);
}

constexpr PieceType drop_type(Move m) {
	return PieceType(from_sq(m));
}

constexpr PieceType hand_type(PieceType pt) {
	return PieceType(pt & 7);
}

constexpr PieceType promoted_type(PieceType pt) {
	return PieceType(pt + PIECE_PROMOTION);
}

inline Move make_move(Square from, Square to) {
	return Move(to + (from << 7));
}

inline Move make_move_promote(Square from, Square to) {
	return Move(to + (from << 7) + PROMOTION);
}

inline Move make_move_drop(PieceType pt, Square to) {
	return Move(to + (pt << 7) + DROP);
}

constexpr bool is_promotable(PieceType pt) {
	return PAWN <= pt && pt <= ROOK;
}

constexpr Piece opponent_piece(Piece pc) {
	assert(pc != NO_PIECE);
	return Piece(pc ^ 0x10);
}

constexpr Piece promoted_piece(Piece pc) {
	assert(pc != NO_PIECE);
	return is_promotable(type_of(pc)) ? Piece(pc | 0x8) : NO_PIECE;
}

constexpr int get_max_number(PieceType pt) {
	const int numPieces[KING] = {
		2, 18, 4, 4, 4, 2, 2, 4 // ã , ï‡, çÅ, åj, ã‚, äp, îÚ, ã‡
	};
	return numPieces[pt];
}

constexpr bool is_ok(Move m) {
	return (m >> 7) != (m & 0x7f);
}

constexpr bool is_promotion_area(Color c, Square s) {
	Rank r = rank_of(s);
	return static_cast<bool>(0x1c00007u & (1u << ((c << 4) + r)));
	// 1 1100 0000 0000 0000 0000 0111
}

constexpr bool can_promote(Piece pc) {
	assert(pc != NO_PIECE);
	return (B_PAWN <= pc && pc <= B_ROOK)
		|| (W_PAWN <= pc && pc <= W_ROOK);
}

#endif // ifndef TYPES_H_INCLUDED