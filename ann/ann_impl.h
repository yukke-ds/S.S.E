// ann.cppに実装されてない理由は，VC++のtemplate解釈がinclusion-modelであり，
// 実装しようとするとリンクエラ―になるから．よってann_impl.hとヘッダファイルに実装してある．
// ちなみに，separation-modelならann.cppに実装することが恐らく可能．

#include "ann.h"

#include <iostream>
#include <memory>
#include <algorithm>

#include <cstdint>

#include <omp.h>

// for floating point interrupts
#include <xmmintrin.h>

#include "../misc.h"
#include "../thread.h"

//#define SGDM
#define ADADELTA

#if defined(SGDM) && defined(ADADELTA)
#error Only select one training method!
#elif !defined(SGDM) && !defined(ADADELTA)
#error Must select one training method!
#endif

inline void EnableNanInterrupt()
{
	_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
FCANN<ACTF, ACTFLast>::FCANN(
	size_t inputs,
	size_t outputs,
	std::vector<size_t> hiddenLayers,
	std::vector<std::vector<Eigen::Triplet<FP> > > &connectionMatrices)
{
	assert(connectionMatrices.size() == (hiddenLayers.size() + 1));

	auto mt = gRd.make_mt();

	// then we build the weight and bias vectors, and initialize them
	// ANNを作って，重みとバイアスをランダム値で初期化する．
	for (size_t layer = 0; layer < (hiddenLayers.size() + 1); ++layer)
	{
		size_t inSize = (layer == 0) ? inputs : hiddenLayers[layer - 1];
		size_t outSize = (layer == hiddenLayers.size()) ? outputs : hiddenLayers[layer];

		NNMatrix weightMatrix(inSize, outSize);
		NNVector biasVector(outSize);

		// determine what distribution to use to initialize hidden weights depending on activation func
		// (output layer is always linear)
		std::function<FP()> drawFunc;

		if (ACTF == Linear || layer == hiddenLayers.size())
		{
			std::uniform_real_distribution<FP> dist(-0.01f, 0.01f);
			drawFunc = std::bind(dist, mt);
		}
		else if (ACTF == Tanh)
		{
			// for tanh, we use r = sqrt(6/(fan_in + fan_out)), (-r, r)
			FP r = sqrt(6.0 / (inSize + outSize));
			std::uniform_real_distribution<FP> dist(-r, r);
			drawFunc = std::bind(dist, mt);
		}
		else if (ACTF == Relu)
		{
			// we use the scheme described here - http://arxiv.org/pdf/1502.01852v1.pdf
			std::normal_distribution<FP> dist(0.0f, sqrt(2.0f / outSize));
			drawFunc = std::bind(dist, mt);
		}
		else
			assert(false);

		// actually initialize weights and biases
		for (size_t j = 0; j < outSize; ++j)
			biasVector(j) = 0.0f;

		for (size_t i = 0; i < inSize; ++i)
			for (size_t j = 0; j < outSize; ++j)
				weightMatrix(i, j) = drawFunc();

		m_params.outputBias.push_back(biasVector);
		m_params.weights.push_back(weightMatrix);

		// サイズが0なら全結合, そうでなければスパースレイヤ(重みを一部伝播させないレイヤ)
		// そのためにマスク行列を作って, これを乗ずることで出力を0とする
		if (connectionMatrices[layer].size() != 0)
		{
			// we have a sparse layer
			NNMatrix conn = NNMatrix::Zero(inSize, outSize);

			for (const auto &trip : connectionMatrices[layer])
				conn(trip.row(), trip.col()) = 1.0f;

			m_params.weightMasks.push_back(conn);
		}
		else
			// we have a fully connected layer
			m_params.weightMasks.push_back(NNMatrix::Ones(inSize, outSize));

		// SGD with momentumで利用するベクトル・行列のゼロ初期化
		m_params.outputBiasLastUpdate.push_back(NNVector::Zero(outSize));
		m_params.weightsLastUpdate.push_back(NNMatrix::Zero(inSize, outSize));

		// ADADELTAで利用するベクトル・行列のゼロ初期化
		m_params.outputBiasEg2.push_back(NNVector::Zero(outSize));
		m_params.weightsEg2.push_back(NNMatrix::Zero(inSize, outSize));
		m_params.outputBiasRMSd2.push_back(NNVector::Zero(outSize));
		m_params.weightsRMSd2.push_back(NNMatrix::Zero(inSize, outSize));
	}

	// NNを順伝播させるときに用いる一時変数のリサイズ
	m_params.evalTmp.resize(hiddenLayers.size() + 2);
	m_params.evalSingleTmp.resize(hiddenLayers.size() + 2);

	// 重み行列を最適な形で保持しておく(後々の順伝播における計算が簡単になる?)
	update_weight_masks_regions();
	update_weight_semi_sparse();
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
void FCANN<ACTF, ACTFLast>::initialize_activations(Activations &act)
{
	assert(m_params.weights.size() == m_params.outputBias.size());

	act.act.clear();
	act.actIn.clear();

	// we don't know how many rows these matrices will have yet, since
	// that depends on input batch size
	for (size_t layer = 0; layer < m_params.weights.size(); ++layer)
	{
		// レイヤ毎に入力×重み行列(1×row)が考えうる，全部ゼロ初期化
		act.act.push_back(NNVector::Zero(1, m_params.weights[layer].rows()));
		act.actIn.push_back(NNVector::Zero(1, m_params.weights[layer].rows()));
	}

	// [0]1×row, [1]1×row, [2]1×rowで，最後[2]1×colも考えるみたい
	act.act.push_back(NNVector::Zero(1, m_params.weights[m_params.weights.size() - 1].cols()));
	act.actIn.push_back(NNVector::Zero(1, m_params.weights[m_params.weights.size() - 1].cols()));
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
void FCANN<ACTF, ACTFLast>::initialize_gradients(Gradients &grad)
{
	assert(m_params.weights.size() == m_params.outputBias.size());

	grad.weightGradients.clear();
	grad.biasGradients.clear();

	for (size_t layer = 0; layer < m_params.weights.size(); ++layer)
	{
		// 重み勾配行列はレイヤの大きさ分(row×col)分をゼロ初期化
		grad.weightGradients.push_back(NNMatrix::Zero(m_params.weights[layer].rows(), m_params.weights[layer].cols()));

		// バイアスは各レイヤへの入力につき1つなので(1×col)分をゼロ初期化で良い
		grad.biasGradients.push_back(NNVector::Zero(1, m_params.weights[layer].cols()));
	}
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
template <typename Derived>
NNMatrixRM FCANN<ACTF, ACTFLast>::forward_propagate(const MatrixBase<Derived> &in, Activations &act)
{
	// それぞれANNのレイヤ数分のベクトルがあるはず
	assert(act.act.size() == m_params.weights.size() + 1);
	assert(act.actIn.size() == m_params.weights.size() + 1);

	act.act[0] = in;
	act.actIn[0] = in; // first layer has no activation

	NNMatrixRM x = in;

	// レイヤの層数分の順伝播を行う
	for (size_t layer = 0; layer < m_params.weights.size(); ++layer)
	{
		x *= m_params.weights[layer]; // 入力に重みを掛ける

		x.rowwise() += m_params.outputBias[layer]; // 行ベクトルにバイアスを足す

		act.actIn[layer + 1] = x; // 活性化関数に通す前の[入力×重み(+バイアス)]

		activate(x, layer == (m_params.weights.size() - 1)); // 活性化関数に通す

		act.act[layer + 1] = x; // 活性化関数に通した後の[入力×重み](次の層への入力)
	}

	return x;
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
template <typename Derived>
NNMatrixRM FCANN<ACTF, ACTFLast>::forward_propagate_fast(const MatrixBase<Derived> &in)
{
	//if (!m_params.weightsSemiSparseCurrent)
	//	update_weight_semi_sparse();

	// 簡略化した順伝播(高速?)
	for (size_t layer = 0; layer < m_params.weights.size(); ++layer)
	{
		// 入力×重み+バイアス
		if (layer == 0)
			m_params.evalTmp[layer].noalias() = in * m_params.weights[layer];
		//MatrixMultiplyWithSemiSparse(in, m_params.weightsSemiSparse[layer], m_params.evalTmp[layer]);
		else
			m_params.evalTmp[layer].noalias() = m_params.evalTmp[layer - 1] * m_params.weights[layer];
		//MatrixMultiplyWithSemiSparse(m_params.evalTmp[layer - 1], m_params.weightsSemiSparse[layer], m_params.evalTmp[layer]);

		m_params.evalTmp[layer].rowwise() += m_params.outputBias[layer];

		// 活性化関数
		activate(m_params.evalTmp[layer], layer == (m_params.weights.size() - 1));
	}

	return m_params.evalTmp[m_params.weights.size() - 1];
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
template <typename Derived>
float FCANN<ACTF, ACTFLast>::forward_propagate_single(const MatrixBase<Derived> &vec)
{
	if (!m_params.weightsSemiSparseCurrent)
		update_weight_semi_sparse();

	for (size_t layer = 0; layer < m_params.weights.size(); ++layer)
	{
		if (layer == 0)
			//m_params.evalSingleTmp[layer].noalias() = vec * m_params.weights[layer];
			multiply_with_semi_sparse(vec, m_params.weightsSemiSparse[layer], m_params.evalSingleTmp[layer]);
		else
			//m_params.evalSingleTmp[layer].noalias() = m_params.evalSingleTmp[layer - 1] * m_params.weights[layer];
			multiply_with_semi_sparse(m_params.evalSingleTmp[layer - 1], m_params.weightsSemiSparse[layer], m_params.evalSingleTmp[layer]);

		m_params.evalSingleTmp[layer] += m_params.outputBias[layer];

		activate(m_params.evalSingleTmp[layer], layer == (m_params.weights.size() - 1));
	}

	return m_params.evalSingleTmp[m_params.weights.size() - 1](0, 0);
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
template <typename Derived>
void FCANN<ACTF, ACTFLast>::backward_propagate_compute_grad(const MatrixBase<Derived> &err, const Activations &act, Gradients &grad)
{
	assert(grad.weightGradients.size() == m_params.weights.size());
	assert(grad.biasGradients.size() == m_params.outputBias.size());
	assert(grad.weightGradients.size() == grad.biasGradients.size());

	// currError are the errorTerms of the next layer
	NNMatrixRM errorTerms = err;

	// 逆伝播なので，layerは後ろから前に減らす
	for (int32_t layer = (m_params.weights.size() - 1); layer >= 0; --layer)
	{
		// 行列演算で行と列の関係が崩れていないか確認
		// 重みは入力数×レイヤ数，バイアスは1×レイヤ数
		assert(grad.weightGradients[layer].rows() == m_params.weights[layer].rows());
		assert(grad.weightGradients[layer].cols() == m_params.weights[layer].cols());
		assert(grad.biasGradients[layer].rows() == 1);
		assert(grad.biasGradients[layer].cols() == m_params.outputBias[layer].cols());

		// first we calculate weight gradients for the current layer,
		// which is the transpose of each input to this layer, multiplied by errorTerms
		// NNの最後，act(活性化関数に通した入力×重み(すでに微分した活性化関数に通してある))を
		// 転置させてエラーと乗ずる
		grad.weightGradients[layer].noalias() = act.act[layer].transpose() * errorTerms;

		// bias gradients are just errorTerms
		// バイアス勾配はエラーの列ベクトルの和に同じ
		grad.biasGradients[layer].noalias() = errorTerms.colwise().sum();

		// レイヤの手前(奥)まできたので，次の活性化関数に通す前の入力×重みを微分する
		NNMatrixRM derivatives = act.actIn[layer];
		activate_derivative(derivatives);

		// then we calculate error for the next (previous) layer
		// 次(前)のレイヤで用いるエラーを計算する．重みを転置して掛けている
		errorTerms *= m_params.weights[layer].transpose();
		errorTerms.array() *= derivatives.array();
	}
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
template <typename Derived1, typename Derived2>
float FCANN<ACTF, ACTFLast>::train_GDM(const MatrixBase<Derived1> &x, const MatrixBase<Derived2> &y, float learningRate, float reg)
{
	static std::vector<Gradients> gradLocal; // NN内の重みとバイアス項の勾配を入れる配列
	static std::vector<Activations> actLocal; // NN内で用いる活性化関数f(x)のxを入れる配列
	static bool initialized = false;

	// we limit to 8 threads for the current block size of 256
	// 現在のブロック数256なら8スレッドに制限しておく
	ScopedOmpLimiter tlim(8);

	// 最初の1回のみ初期化を行う
	if (!initialized)
	{
		gradLocal = std::vector<Gradients>(omp_get_max_threads());
		actLocal = std::vector<Activations>(omp_get_max_threads());

		for (int64_t i = 0; i < omp_get_max_threads(); ++i)
		{
			initialize_activations(actLocal[i]);
			initialize_gradients(gradLocal[i]);
		}

		initialized = true;
	}

	// 誤差(エラー)の和
	float errorsMeasureTotal = 0.0f;

	// 順伝搬で求める予測価値predと教師yとの誤差(pred - y)の勾配を逆伝搬で求める
	// ここでは，重みの勾配(更新量)を計算するだけ
#pragma omp parallel
	{
		int64_t begin;
		int64_t numRows;

		// 並列して計算するので，スレッドにどのブロック(x*w)を割り当てるか決める
		get_thread_block(x.rows(), begin, numRows);

		size_t threadId = omp_get_thread_num(); // スレッド番号
		size_t numThreads = omp_get_num_threads(); // スレッド数

												   // スレッドごとに割り当てられたところを順伝搬
		auto pred = forward_propagate(x.block(begin, 0, numRows, x.cols()), actLocal[threadId]);

		// 予測価値predと教師yとの誤差(pred - y)の和を計算する
		// その誤差も活性化関数を通して計算するみたい
		errorsMeasureTotal += error_func(pred, y.block(begin, 0, numRows, y.cols())).sum();

		// 逆伝播の最初で必要な，誤差を通したtanhx(actIn)の微分(エラー)を計算しておく
		NNMatrixRM errorsDerivative = error_func_derivative(pred, y.block(begin, 0, numRows, y.cols()), actLocal[threadId].actIn[actLocal[threadId].actIn.size() - 1]);

		// 逆伝搬（バックプロパゲーション）
		backward_propagate_compute_grad(errorsDerivative, actLocal[threadId], gradLocal[threadId]);

		// reduce all the local gradients into total, using log(n) steps
		// 計算した勾配を集計する(O(log(n))
		for (size_t skip = 2; skip <= numThreads; skip *= 2)
		{
			if ((threadId % skip) == 0 && (threadId + skip / 2) < numThreads)
			{
				gradLocal[threadId] += gradLocal[threadId + skip / 2];
			}
#pragma omp barrier // バリアまでの処理が終わるまで待機させる(明示的同期)
		}
	} // parallel

	  // 求めた更新量を確率的勾配降下法に適用する
	apply_weight_updates(gradLocal[0], learningRate, reg);

	// 単純平均をとって返す
	return errorsMeasureTotal / x.rows();
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
void FCANN<ACTF, ACTFLast>::apply_weight_updates(const Gradients &grad, float /*learningRate*/, float reg)
{
	assert(grad.weightGradients.size() == m_params.weights.size());
	assert(grad.biasGradients.size() == m_params.outputBias.size());
	assert(grad.weightGradients.size() == grad.biasGradients.size());

	// SGD with momentum
	m_params.weightsLastUpdate.resize(m_params.weights.size());
	m_params.outputBiasLastUpdate.resize(m_params.outputBias.size());

	// ADADELTA
	m_params.weightsEg2.resize(m_params.weights.size());
	m_params.outputBiasEg2.resize(m_params.outputBias.size());
	m_params.weightsRMSd2.resize(m_params.weights.size());
	m_params.outputBiasRMSd2.resize(m_params.outputBias.size());

	for (size_t layer = 0; layer < m_params.weights.size(); ++layer)
	{
#pragma omp parallel
		{
			int64_t begin;
			int64_t numCols;

			size_t inSize = m_params.weights[layer].rows();
			size_t outSize = m_params.weights[layer].cols();

			get_thread_block(outSize, begin, numCols);

			if (numCols != 0) // if numCols is less than num threads, some threads won't have anything to do
			{
				// スレッド毎に重み行列のブロック割り当てる
				auto weightsBlock = m_params.weights[layer].block(0, begin, inSize, numCols);
				auto biasBlock = m_params.outputBias[layer].block(0, begin, 1, numCols);

				auto weightsGradientsBlock = grad.weightGradients[layer].block(0, begin, inSize, numCols);
				auto biasGradientsBlock = grad.biasGradients[layer].block(0, begin, 1, numCols);

				auto weightMaskBlock = m_params.weightMasks[layer].block(0, begin, inSize, numCols);

#ifdef ADADELTA
				auto weightsEg2Block = m_params.weightsEg2[layer].block(0, begin, inSize, numCols);
				auto biasEg2Block = m_params.outputBiasEg2[layer].block(0, begin, 1, numCols);
				auto weightsRMSd2Block = m_params.weightsRMSd2[layer].block(0, begin, inSize, numCols);
				auto biasRMSd2Block = m_params.outputBiasRMSd2[layer].block(0, begin, 1, numCols);
#endif

#ifdef SGDM
				auto weightsLastUpdateBlock = m_params.weightsLastUpdate[layer].block(0, begin, inSize, numCols);
				auto outputBiasLastUpdateBlock = m_params.outputBiasLastUpdate[layer].block(0, begin, 1, numCols);
#endif
#define L1_REG // L1正則化
#ifdef L1_REG
				NNMatrix weightReg(weightsBlock.rows(), weightsBlock.cols());

				for (int64_t i = 0; i < (weightReg.rows() * weightReg.cols()); ++i)
				{
					float w = weightsBlock.data()[i];
					float x;

					if (w > 0.0f) {
						if (w > reg) x = -reg;
						else x = -w;
					}
					else {
						if (w < -reg) x = reg;
						else x = -w;
					}

					weightReg.data()[i] = x;
				}
#elif defined(L2_REG) // L2正則化
				NNMatrix weightReg = -reg * weightsBlock;
#else
				NNMatrix weightReg = NNMatrix::Zero(weightsBlock.rows(), weightsBlock.cols());
#endif

#ifdef ADADELTA
				// update Eg2 (Eg2(t) = γEg2(t-1) + (1-γ)g2^2)
				float decay = 0.99f;
				float e = 1e-8f;
				weightsEg2Block.array() *= decay;
				weightsEg2Block.array() += (weightsGradientsBlock.array() * weightsGradientsBlock.array()) * (1.0f - decay);
				biasEg2Block.array() *= decay;
				biasEg2Block.array() += (biasGradientsBlock.array() * biasGradientsBlock.array()) * (1.0f - decay);

				// -g(t) * (RMS[Δw](t-1) / RMS[g](t))
				NNMatrix weightDelta = -weightsGradientsBlock.array() * (weightsRMSd2Block.array() + e).sqrt() / (weightsEg2Block.array() + e).sqrt() + weightReg.array();
				NNVector biasDelta = -biasGradientsBlock.array() * (biasRMSd2Block.array() + e).sqrt() / (biasEg2Block.array() + e).sqrt();
#endif
#ifdef SGDM
				float lr = 0.000001f;
				float momentum = 0.95f;
				NNMatrix weightDelta = -weightsGradientsBlock.array() * lr + momentum * weightsLastUpdateBlock.array()/*+ weightReg.array()*/;
				NNVector biasDelta = -biasGradientsBlock.array() * lr + momentum * outputBiasLastUpdateBlock.array();
#endif

				// 学習率を掛けて重みとバイアスを更新
				weightsBlock += weightDelta/* * learningRate*/;
				weightsBlock.array() *= weightMaskBlock.array();
				biasBlock += biasDelta/* * learningRate*/;

				FP weightMax = std::max(std::max(weightsBlock.maxCoeff(), -weightsBlock.minCoeff()), std::max(biasBlock.maxCoeff(), -biasBlock.minCoeff()));
				if (weightMax > MAX_WEIGHT)
					throw std::runtime_error("Learning rate too high!");

#ifdef ADADELTA
				// RMS[Δw](t)を更新
				weightsRMSd2Block *= decay;
				weightsRMSd2Block.array() += weightDelta.array() * weightDelta.array() * (1.0f - decay);
				biasRMSd2Block *= decay;
				biasRMSd2Block.array() += biasDelta.array() * biasDelta.array() * (1.0f - decay);
#endif
#ifdef SGDM
				weightsLastUpdateBlock = weightDelta;
				outputBiasLastUpdateBlock = biasDelta;
#endif
			}
		} // omp parallel
	}

	m_params.weightsSemiSparseCurrent = false;
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
float FCANN<ACTF, ACTFLast>::get_sparsity() // スパース(希薄さ)推定
{
	uint64_t zCount = 0;
	uint64_t totalCount = 0;

	for (size_t layer = 0; layer < m_params.weights.size(); ++layer)
	{
		totalCount += m_params.weights[layer].rows() * m_params.weights[layer].cols();

		// 重みゼロの部分がどれくらいあるか
		for (int64_t i = 0; i < m_params.weights[layer].size(); ++i)
			if (m_params.weights[layer].data()[i] == 0.0f)
				++zCount;
	}

	return static_cast<float>(zCount) / totalCount;
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
template <typename Derived1, typename Derived2>
NNMatrixRM FCANN<ACTF, ACTFLast>::error_func(
	const MatrixBase<Derived1> &pred,
	const MatrixBase<Derived2> &targets) const
{
	NNMatrixRM ret;

	// for linear and tanh output we use MAE
	if (ACTFLast == Linear)
	{
		ret = (pred - targets).array().abs().matrix();
	}
	else if (ACTFLast == Tanh)
	{
		//ret = (pred - targets).array().abs().matrix();
		ret = ((pred - targets).array() * (pred - targets).array()).matrix();
	}
	// for softmax output we use cross-entropy
	else if (ACTFLast == Softmax)
	{
		// note that for cross-entropy, output is a vector no matter how many classes we have
		ret.resize(pred.rows(), 1);

		for (int64_t i = 0; i < pred.rows(); ++i)
		{
			float e = 0.0f;

			for (int64_t j = 0; j < pred.cols(); ++j)
				if (targets(i, j) == 1.0f)
					e += -log(pred(i, j));

			ret(i, 0) = e;
		}
	}
	else if (ACTFLast == Logsig)
	{
		// cross-entropy
		// - target * log(p) - (1 - target) * log(1-p)

		ret.resize(pred.rows(), pred.cols());

		for (int64_t i = 0; i < pred.rows(); ++i)
			for (int64_t j = 0; j < pred.cols(); ++j)
				ret(i, j) = -targets(i, j) * log(pred(i, j)) - (1.0f - targets(i, j)) * log(1.0f - pred(i, j));
	}
	else assert(false);

	return ret;
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
template <typename Derived1, typename Derived2, typename Derived3>
NNMatrixRM FCANN<ACTF, ACTFLast>::error_func_derivative(
	const MatrixBase<Derived1> &pred,
	const MatrixBase<Derived2> &targets,
	const MatrixBase<Derived3> &finalLayerActivations) const
{
	NNMatrixRM ret;

	if (ACTFLast == Linear)
	{
		NNMatrixRM err = pred - targets;

		ret.resize(err.rows(), err.cols());

		// MAE (Mean Absolute Error, 平均絶対誤差) は数値予測問題における精度評価指標の1つ．
		// https://crowdsolving.jp/node/1029
		for (int64_t i = 0; i < err.rows(); ++i)
			for (int64_t j = 0; j < err.cols(); ++j)
				ret(i, j) = (err(i, j) > 0.0f) ? 1.0f : -1.0f;
	}
	else if (ACTFLast == Tanh)
	{
		ret = pred - targets;

		// now we have to multiply every element by the derivative at that point
		// derivative of tanh is 1-tanh^2(x)
		// 逆伝播のとき，それぞれのレイヤでtanhxの微分(1-tanh^2(x))が必要
		for (int64_t i = 0; i < ret.rows(); ++i)
		{
			for (int64_t j = 0; j < ret.cols(); ++j)
			{
				float tanhx = tanh(finalLayerActivations(i, j));
				ret(i, j) *= 1 - (tanhx * tanhx);
			}
		}
	}
	else if (ACTFLast == Softmax)
	{
		// 未実装
		// cross-entropy
		// curiously,
		// http://www.willamette.edu/~gorr/classes/cs449/classify.html
		ret = pred - targets;
	}
	else if (ACTFLast == Logsig)
	{
		// 未実装
		// with cross-entropy
		// http://cs229.stanford.edu/notes/cs229-notes1.pdf
		ret = pred - targets;
	}
	else assert(false);

	return ret;
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
template <typename Derived>
void FCANN<ACTF, ACTFLast>::activate(MatrixBase<Derived> &x, bool last) const
{
	ActivationFunc actf = last ? ACTFLast : ACTF;

	// these will all be optimized to just be a single case, since
	// ACTF is a template parameter
	if (actf == Linear)
	{
		return; // nothing to do here
	}
	else if (actf == Tanh)
	{
		for (int32_t i = 0; i < x.cols(); ++i)
			for (int32_t j = 0; j < x.rows(); ++j)
				x(j, i) = tanh(x(j, i));
	}
	else if (actf == Relu)
	{
		x = x.array().max(NNMatrix::Zero(x.rows(), x.cols()).array());
	}
	else if (actf == Softmax)
	{
		// the naive implementation is likely to overflow, so we do some shifting first
		// since we are in log space, and dividing in x is subtracting in log(x)
		// dividing all values won't change the distribution

		// we find the max component in each row, and subtract that from each component
		Eigen::Matrix<FP, Eigen::Dynamic, 1> maxElems = x.rowwise().maxCoeff();

		x.colwise() -= maxElems;

		// compute element-wise exp
		x = x.array().exp().matrix();

		// then compute the normalization denominator for each row
		Eigen::Matrix<FP, Eigen::Dynamic, 1> norm = x.rowwise().sum();

		// then normalize all the elements
		x.array().colwise() /= norm.array();
	}
	else if (actf == Logsig)
	{
		// 1 / (exp(-x) + 1)
		x = (1.0f / ((-x).array().exp() + 1)).matrix();
	}
	else assert(false);
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
template <typename Derived>
void FCANN<ACTF, ACTFLast>::activate_derivative(MatrixBase<Derived> &x) const
{
	// these will all be optimized to just be a single case, since
	// ACTF is a template parameter
	if (ACTF == Linear)
	{
		x = NNMatrixRM::Ones(x.rows(), x.cols());
	}
	else if (ACTF == Tanh)
	{
		// derivative of tanh is 1-tanh^2(x)
		for (int64_t i = 0; i < x.cols(); ++i)
		{
			for (int64_t j = 0; j < x.rows(); ++j)
			{
				FP tanhx = tanh(x(j, i));
				x(j, i) = 1 - (tanhx * tanhx);
			}
		}
	}
	else if (ACTF == Relu)
	{
		for (int64_t i = 0; i < x.cols(); ++i)
		{
			for (int64_t j = 0; j < x.rows(); ++j)
				if (x(j, i) > 0)
					x(j, i) = 1.0f;
				else
					x(j, i) = 0.0f;
		}
		// これだとエラーになる
		//x.array() = (x.array() > NNMatrixRM::Zero(x.rows(), x.cols()).array());
	}
	else assert(false);
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
void FCANN<ACTF, ACTFLast>::update_weight_masks_regions()
{
	m_params.weightMasksRegions.resize(m_params.weightMasks.size());

	for (size_t layer = 0; layer < m_params.weightMasks.size(); ++layer)
	{
		WeightMaskType toConvert = m_params.weightMasks[layer];

		m_params.weightMasksRegions[layer] = matrix_to_regions(toConvert);
	}

	m_params.weightsSemiSparseCurrent = false;
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
void FCANN<ACTF, ACTFLast>::update_weight_semi_sparse()
{
	m_params.weightsSemiSparse.resize(m_params.weightMasks.size());

	for (size_t layer = 0; layer < m_params.weightMasks.size(); ++layer)
	{
		WeightType toConvert = m_params.weights[layer];

		m_params.weightsSemiSparse[layer] = to_semi_sparse(toConvert, m_params.weightMasksRegions[layer]);
	}

	m_params.weightsSemiSparseCurrent = true;
}

/* serialization format:
* numLayers
* for each layer:
*		weight matrix
*		weight mask
*		bias
*
* For each matrix:
*	rows
*	cols
*  each field in row major format (rows * cols)
*/

namespace
{

	template <typename Derived>
	void push_matrix(Eigen::MatrixBase<Derived> &m, std::ostream &s)
	{
		s << m.rows() << ' ' << '\n';
		s << m.cols() << ' ' << '\n';

		for (size_t row = 0; row < m.rows(); ++row)
		{
			for (size_t col = 0; col < m.cols(); ++col)
				s << m(row, col) << ' ';
			s << '\n';
		}
	}

	NNMatrix read_matrix(std::istream &s)
	{
		size_t nRows;
		size_t nCols;

		s >> nRows;
		s >> nCols;

		NNMatrix ret(nRows, nCols);

		for (size_t row = 0; row < nRows; ++row)
			for (size_t col = 0; col < nCols; ++col)
				s >> ret(row, col);

		return ret;
	}

	NNMatrix read_weight_py(std::istream &s, size_t& nRows, size_t& nCols)
	{
		std::vector<float> w(nRows * nCols, 0.0);
		s.read((char*)w.data(), sizeof(float) * nRows * nCols);

		//for (size_t i = 0; i < nCols; ++i) {
		//	for (size_t j = 0; j < nRows; ++j)
		//		std::cout << w[i * nRows + j] << " ";
		//	std::cout << std::endl;
		//	std::cout << "---------------------------------------------" << std::endl;
		//}

		NNMatrix ret(nRows, nCols);
		for (size_t col = 0; col < nCols; ++col)
			for (size_t row = 0; row < nRows; ++row)
				ret(row, col) = w[col * nRows + row];

		return ret;
	}

	NNMatrix read_bias_py(std::istream &s, size_t& nRows, size_t& nCols)
	{
		std::vector<float> w(nRows * nCols, 0.0);
		s.read((char*)w.data(), sizeof(float) * nRows * nCols);

		//for (size_t i = 0; i < nRows; ++i) {
		//	for (size_t j = 0; j < nCols; ++j)
		//		std::cout << w[i * nCols + j] << " ";
		//	std::cout << std::endl;
		//	std::cout << "---------------------------------------------" << std::endl;
		//}

		NNMatrix ret(nRows, nCols);
		for (size_t row = 0; row < nRows; ++row)
			for (size_t col = 0; col < nCols; ++col)
				ret(row, col) = w[row * nCols + col];

		return ret;
	}
}

template <typename T>
void serialize_net(T &net, std::ostream &s)
{
	auto weights = net.get_weights();
	auto biases = net.get_biases();
	auto weightMasks = net.get_weight_masks();

	size_t numLayers = weights.size();

	std::vector<size_t> hiddenLayerSizes;

	for (size_t i = 1; i < numLayers; ++i)
		hiddenLayerSizes.push_back(weights[i].rows());

	s << numLayers << '\n';

	for (size_t i = 0; i < numLayers; ++i)
	{
		push_matrix(weights[i], s);
		push_matrix(weightMasks[i], s);
		push_matrix(biases[i], s);
	}
}

template <typename T>
void deserialize_net(T &net, std::istream &s)
{
	std::vector<typename T::WeightType> weights;
	std::vector<typename T::BiasType> biases;
	std::vector<typename T::WeightMaskType> weightMasks;

	size_t numLayers;

	s >> numLayers;

	for (size_t i = 0; i < numLayers; ++i)
	{
		weights.push_back(read_matrix(s));
		weightMasks.push_back(read_matrix(s));
		biases.push_back(read_matrix(s));
	}

	size_t din = weights[0].rows();
	size_t dout = weights[weights.size() - 1].cols();

	std::vector<size_t> hiddenLayerSizes;

	for (size_t i = 1; i < numLayers; ++i)
		hiddenLayerSizes.push_back(weights[i].rows());

	// we just set everything to be fully connected, since we will
	// overwrite the connection matrices anyways
	std::vector<std::vector<Eigen::Triplet<FP> > > connections(hiddenLayerSizes.size() + 1);

	net = T(din, dout, hiddenLayerSizes, connections);

	net.get_weights() = weights;
	net.get_biases() = biases;
	net.get_weight_masks() = weightMasks;

	net.notify_weight_masks_changed();
}

template <typename T>
void deserialize_net_py(T &net, std::istream &s)
{
	std::vector<typename T::WeightType> weights;
	std::vector<typename T::BiasType> biases;
	std::vector<typename T::WeightMaskType> weightMasks;

	const size_t NumLayers = 3;
	const size_t In = 631, Units = 300, Out = 1;
	size_t biasRow = 1;

	for (size_t i = 0; i < NumLayers; ++i)
	{
		size_t n_in = i ? Units : In;
		size_t n_out = i == NumLayers - 1 ? Out : Units;

		weights.push_back(read_weight_py(s, n_in, n_out));
		weightMasks.push_back(NNMatrix::Ones(n_in, n_out));
		biases.push_back(read_bias_py(s, biasRow, n_out));
	}

	size_t din = weights[0].rows();
	size_t dout = weights[weights.size() - 1].cols();

	std::vector<size_t> hiddenLayerSizes;

	for (size_t i = 1; i < NumLayers; ++i)
		hiddenLayerSizes.push_back(weights[i].rows());

	// we just set everything to be fully connected, since we will
	// overwrite the connection matrices anyways
	std::vector<std::vector<Eigen::Triplet<FP> > > connections(hiddenLayerSizes.size() + 1);

	net = T(din, dout, hiddenLayerSizes, connections);

	net.get_weights() = weights;
	net.get_biases() = biases;
	net.get_weight_masks() = weightMasks;

	net.notify_weight_masks_changed();
}