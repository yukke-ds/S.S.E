#ifndef MATH_H_INCLUDED
#define MATH_H_INCLUDED

#include <cmath>

namespace math {

// 符号関数
template <typename T>
inline T sign(T x) {
	return x > T(0) ? T(1) : (x < T(0) ? T(-1) : T(0));
}

// 対数関数
inline float log_f(float x) {
	return log(x);
}

// シグモイド関数
inline double sigmoid(double x) {
	return 1.0 / (1.0 + std::exp(-x));
}

inline float sigmoid_f(float x) {
	return 1.0 / (1.0 + std::exp(-x));
}

// シグモイド関数の導関数
inline double derivative_of_sigmoid(double x) {
	return sigmoid(x) * (1.0 - sigmoid(x));
}

} // namespace math

#endif // ifndef MATH_H_INCLUDED