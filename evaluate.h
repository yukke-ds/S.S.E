#ifndef EVALUATE_H_INCLUDED
#define EVALUATE_H_INCLUDED

#include <array>

#include "evallist.h"
#include "types.h"

class Position;

namespace Eval {

// std::array<T,2>に対して += と -= を提供する。
template <typename Tl, typename Tr>
inline std::array<Tl, 2> operator += (std::array<Tl, 2>& lhs, const std::array<Tr, 2>& rhs) {
	lhs[0] += rhs[0];
	lhs[1] += rhs[1];
	return lhs;
}

template <typename Tl, typename Tr>
inline std::array<Tl, 2> operator -= (std::array<Tl, 2>& lhs, const std::array<Tr, 2>& rhs) {
	lhs[0] -= rhs[0];
    lhs[1] -= rhs[1];
    return lhs;
}


void load_eval();
uint64_t calc_check_sum();

Value material(const Position& pos);
Value evaluate(const Position& pos);
Value all_calculate(const Position& pos);
void evaluate_with_no_return(const Position& pos);

const Value Tempo = Value(20); // Must be visible to search
const int FV_SCALE = 32;
const int FULL_SCALE = 10000;

extern std::array<int32_t, 2> kk[SQUARE_NB][SQUARE_NB];
extern std::array<int32_t, 2> kkp[SQUARE_NB][SQUARE_NB][FE_END];
extern std::array<int16_t, 2> kpp[SQUARE_NB][FE_END][FE_END];
extern Value PieceValue[PIECE_NB];
extern Value CapturePieceValue[PIECE_NB];
extern Value PromotionDiff[PIECE_NB];

}

// 手番つきの評価値を足していくときに使うclass
struct alignas(32) Evaluate {

	Evaluate() {}
#if defined(USE_PEXT)
	Evaluate(const Evaluate& es) { _mm256_store_si256(&mm, es.mm); }
#else
	Evaluate::Evaluate(const Evaluate& es) {
		_mm_store_si128(&m[0], es.m[0]);
		_mm_store_si128(&m[1], es.m[1]);
	}
#endif
	static void init();

	int32_t sum(const Color c) const;
	bool calculated() const;

	Evaluate& operator = (const Evaluate& rhs);
	Evaluate& operator += (const Evaluate& rhs);
	Evaluate& operator -= (const Evaluate& rhs);
	Evaluate operator + (const Evaluate& rhs) const;
	Evaluate operator - (const Evaluate& rhs) const;

	union {
		// array<.. , 3>でいいが、この下のstructに合わせてpaddingしておく。
		std::array<std::array<int32_t, 2>, 4> p;

#if defined(USE_PEXT)
		__m256i mm;
		__m128i m[2];
#else
		__m128i m[2];
#endif
	};
};

inline bool Evaluate::calculated() const {
	return p[0][0] != VALUE_NOT_EVALUATED;
}

#if defined(USE_PEXT)
inline Evaluate& Evaluate::operator = (const Evaluate& rhs) {
	_mm256_store_si256(&mm, rhs.mm);
	return *this;
}
#else
inline Evaluate& Evaluate::operator = (const Evaluate& rhs) {
	_mm_store_si128(&m[0], rhs.m[0]);
	_mm_store_si128(&m[1], rhs.m[1]);
	return *this;
}
#endif

inline Evaluate& Evaluate::operator += (const Evaluate& rhs) {
#if defined(USE_PEXT)
	mm = _mm256_add_epi32(mm, rhs.mm);
#else
	m[0] = _mm_add_epi32(m[0], rhs.m[0]);
	m[1] = _mm_add_epi32(m[1], rhs.m[1]);
#endif
	return *this;
}

inline Evaluate& Evaluate::operator -= (const Evaluate& rhs) {
#if defined(USE_PEXT)
	mm = _mm256_sub_epi32(mm, rhs.mm);
#else
	m[0] = _mm_sub_epi32(m[0], rhs.m[0]);
	m[1] = _mm_sub_epi32(m[1], rhs.m[1]);
#endif
	return *this;
}

inline Evaluate Evaluate::operator + (const Evaluate& rhs) const {
	return Evaluate(*this) += rhs;
}

inline Evaluate Evaluate::operator - (const Evaluate& rhs) const { 
	return Evaluate(*this) -= rhs; 
}

#endif // #ifndef EVALUATE_H_INCLUDED