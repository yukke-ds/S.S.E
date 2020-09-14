#ifndef ANN_H_INCLUDED
#define ANN_H_INCLUDED

#include <array>
#include <algorithm>
#include <random>
#include <functional>
#include <memory>
#include <exception>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <atomic>
#include <ostream>
#include <istream>

#include <cmath>
#include <cassert>

#include "../common/matrix_ops.h"

#define LINEAR
//#define TANH

#if defined(LINEAR) && defined(TANH)
#error Only select one output activation!
#elif !defined(LINEAR) && !defined(TANH)
#error Must select one output activation!
#endif

enum ActivationFunc
{
	Linear,
	Tanh,
	Relu,
	Softmax,
	Logsig
};

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
class FCANN // Feature Conversion Artificial Neural Network
{
public:
	FCANN() {}

	// initialize with random weights
	FCANN(
		size_t inputs,
		size_t outputs,
		std::vector<size_t> hiddenLayers,
		std::vector<std::vector<Eigen::Triplet<FP> > > &connectionMatrices);

	struct Activations
	{
		// input into each layer
		// 活性化関数に通す前のレイヤの入力×重み行列
		std::vector<NNMatrixRM> act;

		// input into activation functions for each layer
		// 活性化関数に通した後のレイヤの入力×重み行列
		std::vector<NNMatrixRM> actIn;
	};

	struct Gradients
	{
		std::vector<NNVector> biasGradients;
		std::vector<NNMatrix> weightGradients;

		Gradients &operator+=(const Gradients &other)
		{
			assert(biasGradients.size() == other.biasGradients.size());
			assert(weightGradients.size() == other.weightGradients.size());

			for (size_t i = 0; i < biasGradients.size(); ++i)
			{
				biasGradients[i] += other.biasGradients[i];
				weightGradients[i] += other.weightGradients[i];
			}

			return *this;
		}
	};

	void initialize_activations(Activations &act);
	void initialize_gradients(Gradients &grad);

	template <typename Derived>
	NNMatrixRM forward_propagate(const MatrixBase<Derived> &in, Activations &act);

	// same as ForwardPropagate, but doesn't bother with Activations (NOT REENTRANT!!)
	template <typename Derived>
	NNMatrixRM forward_propagate_fast(const MatrixBase<Derived> &in);

	// special case for 1 board and single-valued output - this is used in gameplay (NOT REENTRANT!!)
	template <typename Derived>
	float forward_propagate_single(const MatrixBase<Derived> &vec);

	template <typename Derived>
	void backward_propagate_compute_grad(const MatrixBase<Derived> &err, const Activations &act, Gradients &grad);

	// this is a convenience function that simply runs 1 iteration of GDM
	// GDMはGradient Derivative Matrices?
	template <typename Derived1, typename Derived2>
	float train_GDM(const MatrixBase<Derived1> &x, const MatrixBase<Derived2> &y, float learningRate, float reg);

	void apply_weight_updates(const Gradients &grad, float learningRate, float reg);

	float get_sparsity();

	typedef NNVector BiasType;
	typedef NNMatrix WeightType;
	typedef NNMatrix WeightMaskType;

	// these are used to save and restore nets
	std::vector<BiasType> &get_biases() { return m_params.outputBias; }
	std::vector<WeightType> &get_weights() { m_params.weightsSemiSparseCurrent = false; return m_params.weights; }
	std::vector<WeightMaskType> &get_weight_masks() { return m_params.weightMasks; }

	void notify_weight_masks_changed() { update_weight_masks_regions(); }

	template <typename Derived1, typename Derived2>
	NNMatrixRM error_func(const MatrixBase<Derived1> &pred, const MatrixBase<Derived2> &targets) const;

	template <typename Derived1, typename Derived2, typename Derived3>
	NNMatrixRM error_func_derivative(const MatrixBase<Derived1> &pred, const MatrixBase<Derived2> &targets, const MatrixBase<Derived3> &finalLayerActivations) const;

private:
	template <typename Derived>
	void activate(MatrixBase<Derived> &x, bool last) const;

	template <typename Derived>
	void activate_derivative(MatrixBase<Derived> &x) const;

	// annの実装に無関係だから，ここに実装してある
	void get_thread_block(int64_t numTotal, int64_t &begin, int64_t &num)
	{
		size_t threadId = omp_get_thread_num();
		size_t numThreads = omp_get_num_threads();

		// 入力xをスレッド数で割れば，商が1つ分で余りは分岐で考慮
		size_t rowsPerThread = numTotal / numThreads;
		size_t rem = numTotal % numThreads; // the first "rem" threads get 1 extra row

		if (threadId < rem)
		{
			begin = threadId * (rowsPerThread + 1);
			num = rowsPerThread + 1;
		}
		else
		{
			begin = rem * (rowsPerThread + 1) + (threadId - rem) * rowsPerThread;
			num = rowsPerThread;
		}
	}

	void update_weight_masks_regions();
	void update_weight_semi_sparse();

	// this is used to ensure network stability
	constexpr static FP MAX_WEIGHT = 1000.0f;

	// these are network parameters that should be copied by copy ctor and assignment operator
	struct Params
	{
		// bias, weights, and weightMasks completely define the net
		// weightsMasks：行列内で重みが定義されてれば1.0, そうでなければ0.0のマスク
		std::vector<BiasType> outputBias;
		std::vector<WeightType> weights;
		std::vector<WeightMaskType> weightMasks;

		// optimized form of weight masks (in lists of regions)
		std::vector<std::vector<MatrixRegion> > weightMasksRegions;

		// optimized form of weight matrices (semi-sparse)
		bool weightsSemiSparseCurrent;
		std::vector<SemiSparseMatrix<WeightType>> weightsSemiSparse;

		// these are temporary variables for evaluating the net,
		// so we don't have to keep allocating and de-allocating
		std::vector<NNMatrixRM> evalTmp;
		std::vector<NNVector> evalSingleTmp;

		// the following 2 fields are used by SGD with momentum
		std::vector<NNVector> outputBiasLastUpdate;
		std::vector<NNMatrix> weightsLastUpdate;

		// the following 4 fields are used by ADADELTA
		// RMS(Root Mean Square:二乗平方平方根)
		std::vector<NNVector> outputBiasEg2;
		std::vector<NNMatrix> weightsEg2;
		std::vector<NNVector> outputBiasRMSd2;
		std::vector<NNMatrix> weightsRMSd2;
	} m_params;
};

#ifdef LINEAR
typedef FCANN<Relu, Linear> EvalNet;
#endif
#ifdef TANH
typedef FCANN<Relu, Tanh> EvalNet;
#endif

template <typename T>
void serialize_net(T &net, std::ostream &s);

template <typename T>
void deserialize_net(T &net, std::istream &s);

template <typename T>
void deserialize_net_py(T &net, std::istream &s);

#include "ann_impl.h"

#endif // ANN_H_INCLUDED