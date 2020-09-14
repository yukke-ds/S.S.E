#ifndef LEARN_ANN_H_INCLUDED
#define LEARN_ANN_H_INCLUDED

#include "../Eigen/Core"

#include "ann.h"

enum Method {
	TD_LEAF,
	SL
};

namespace LearnAnn
{

template <typename Method T>
EvalNet build_eval_net(int64_t inputDims, int64_t outputDims, bool smallNet);

template <typename Derived1, typename Derived2>
void train_ANN(
	const Eigen::MatrixBase<Derived1> &x,
	const Eigen::MatrixBase<Derived2> &y,
	EvalNet &nn,
	int64_t epochs);

} // namespace LearnAnn

#endif // ifndef LEARN_ANN_H_INCLUDED