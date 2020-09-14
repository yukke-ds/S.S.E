#include "ann_evaluator.h"

#include <fstream>

#include "../common/math.h"
#include "../usi.h"

ANNEvaluator AnnEvaluator; // Our global ANN evaluator

namespace {
inline std::string sfen_reverse(const std::string sfen) {

	std::istringstream ss(sfen);
	std::string board, stm, hand, gamePly;

	ss >> board >> stm >> hand >> gamePly;

	std::string boardRev;
	for (int64_t j = board.size() - 1; j >= 0; --j) // 反転させるため逆から
	{
		char token = board[j], plus = (j == 0 ? ' ' : board[j - 1]);
		if (std::isdigit(token) || token == '/') // 数字, '/'はそのまま.('+'は下で入れる)
			boardRev.push_back(token);

		else if (isupper(token)) // 大文字は小文字に
		{
			if (plus == '+') { boardRev.push_back(plus); boardRev.push_back(tolower(token)); }
			else boardRev.push_back(tolower(token));
		}
		else if (islower(token)) // 小文字は大文字に
		{
			if (plus == '+') { boardRev.push_back(plus); boardRev.push_back(toupper(token)); }
			else boardRev.push_back(toupper(token));
		}
	}
	board = boardRev;

	std::string handRev;
	for (int64_t j = 0; j < hand.size(); ++j) { // 持ち駒は前からでOK
		char token = hand[j];
		if (std::isdigit(token) || token == '-') handRev.push_back(token);
		else if (isupper(token)) handRev.push_back(tolower(token));
		else if (islower(token)) handRev.push_back(toupper(token));
	}
	hand = handRev;

	std::string newSfen = board + " " + stm + " " + hand + " " + gamePly;
	return newSfen;
}

} // namespace

template <typename Method T>
void ANNEvaluator::build_ANN(int64_t inputDims)
{
	// メインのANN
	m_mainAnn = LearnAnn::build_eval_net<T>(inputDims, 1, false);

	// Lazy Evaluationで用いるANN
	// https://chessprogramming.wikispaces.com/Lazy+Evaluation
	//m_ubAnn = LearnAnn::BuildEvalNet(inputDims, 1, true);
	//m_lbAnn = LearnAnn::BuildEvalNet(inputDims, 1, true);
}

template void ANNEvaluator::build_ANN<TD_LEAF>(int64_t inputDims);
template void ANNEvaluator::build_ANN<SL>(int64_t inputDims);

void ANNEvaluator::serialize(std::ostream &os)
{
	serialize_net(m_mainAnn, os);
	//serialize_net(m_ubAnn, os);
	//serialize_net(m_lbAnn, os);
}

void ANNEvaluator::deserialize(std::istream &is)
{
	deserialize_net(m_mainAnn, is);
	//deserialize_net(m_ubAnn, is);
	//deserialize_net(m_lbAnn, is);
}

void ANNEvaluator::deserialize_py(std::istream &is)
{
	deserialize_net_py(m_mainAnn, is);
}

void ANNEvaluator::train(const std::vector<std::string> &positions, const NNMatrixRM &y, const std::vector<FeaturesConv::FeatureDescription> &featureDescriptions, float learningRate)
{
	auto x = boards_to_feature_representation(positions, featureDescriptions);

	NNMatrixRM predictions;
	EvalNet::Activations act;

	m_mainAnn.initialize_activations(act);

	predictions = m_mainAnn.forward_propagate(x, act);

	NNMatrixRM errorsDerivative = compute_error_derivatives(predictions, y, act.actIn[act.actIn.size() - 1], 1.0f, 1.0f);

	EvalNet::Gradients grad;

	m_mainAnn.initialize_gradients(grad);

	m_mainAnn.backward_propagate_compute_grad(errorsDerivative, act, grad);

	m_mainAnn.apply_weight_updates(grad, learningRate, 0.0f);
}

void ANNEvaluator::train_loop(
	const std::vector<std::string> &positions,
	const NNMatrixRM &y,
	int64_t epochs,
	const std::vector<FeaturesConv::FeatureDescription> &featureDescriptions)
{
	// SFEN数分の局面×特徴数の行列の入力x(float)を作る
	auto x = boards_to_feature_representation(positions, featureDescriptions);

	// 入力x，教師y，ANNの準備ができたので学習を始める
	// (train_ANN<NNMatrixRM, NNMatrixRM>としなくてもよい．関数の引数から自明なため)
	LearnAnn::train_ANN(x, y, m_mainAnn, epochs);
}

float ANNEvaluator::forward_triplet(const std::vector<std::string> &observed, const std::vector<std::string> &parent, const std::vector<std::string> &random, const std::vector<FeaturesConv::FeatureDescription> &featureDescriptions, float learningRate)
{
#ifdef BB_INPUT
	auto xo = boards_to_bitboard_feature(observed);
	auto xp = boards_to_bitboard_feature(parent);
	auto xr = boards_to_bitboard_feature(random);
#endif
#ifdef GIRAFFE_INPUT
	auto xo = boards_to_feature_representation(observed, featureDescriptions);
	auto xp = boards_to_feature_representation(parent, featureDescriptions);
	auto xr = boards_to_feature_representation(random, featureDescriptions);
#endif
	/*std::cout << "# Input features #" << std::endl;
	for (size_t i = 0; i < observed.size(); ++i) {
	std::cout << xo.row(i) << std::endl;
	std::cout << xp.row(i) << std::endl;
	std::cout << xr.row(i) << std::endl;
	}*/

	NNMatrixRM observedPred, parentPred, randomPred;
	EvalNet::Activations act;

	m_mainAnn.initialize_activations(act);

	observedPred = m_mainAnn.forward_propagate(xo, act);
	parentPred = m_mainAnn.forward_propagate(xp, act);
	randomPred = m_mainAnn.forward_propagate(xr, act);

	// loss_a
	NNMatrixRM diff_a = -observedPred + randomPred;
	NNMatrixRM sigDiff_a = diff_a.unaryExpr([&](float x) { 	return 1.0 / (1.0 + std::exp(-x)); });
	NNMatrixRM logDiff_a = sigDiff_a.unaryExpr([&](float x) { return log(x); });
	NNMatrixRM loss_a = logDiff_a;
	float avgLoss_a = logDiff_a.sum() / logDiff_a.size();

	// loss_b
	NNMatrixRM diff_b = observedPred + parentPred;
	NNMatrixRM sigDiff_b = diff_b.unaryExpr([&](float x) { 	return 1.0 / (1.0 + std::exp(-x)); });
	NNMatrixRM logDiff_b = sigDiff_b.unaryExpr([&](float x) { return log(x); });
	NNMatrixRM loss_b = logDiff_b;
	float avgLoss_b = logDiff_b.sum() / logDiff_b.size();

	// loss_c
	NNMatrixRM diff_c = -(observedPred + parentPred);
	NNMatrixRM sigDiff_c = diff_c.unaryExpr([&](float x) { 	return 1.0 / (1.0 + std::exp(-x)); });
	NNMatrixRM logDiff_c = sigDiff_c.unaryExpr([&](float x) { return log(x); });
	NNMatrixRM loss_c = logDiff_c;
	float avgLoss_c = logDiff_c.sum() / logDiff_c.size();

	NNMatrixRM loss = loss_a + loss_b + loss_c;

	float errorsMeasureTotal = avgLoss_a + avgLoss_b + avgLoss_c;

#ifdef LINEAR
	NNMatrixRM errorsDerivative = compute_loss_derivatives<Relu, Linear>(loss, act.actIn[act.actIn.size() - 1], 1.0f, 1.0f);
#endif
#ifdef TANH
	NNMatrixRM errorsDerivative = compute_loss_derivatives<Relu, Tanh>(loss, act.actIn[act.actIn.size() - 1], 1.0f, 1.0f);
#endif
	EvalNet::Gradients grad;

	m_mainAnn.initialize_gradients(grad);

	m_mainAnn.backward_propagate_compute_grad(errorsDerivative, act, grad);

	m_mainAnn.apply_weight_updates(grad, learningRate, 0.0f);

	return errorsMeasureTotal;
}

Value ANNEvaluator::evaluate(const Position& pos)
{
#ifdef BB_INPUT
	std::istringstream ss(pos.sfen());
	std::string board, stm, hand, gamePly;
	ss >> board >> stm >> hand >> gamePly;

	float m_convTmp[FeaturesConv::PositionInputDims] = { 0 };

	FeaturesConv::convert_sfen_to_BB(board, stm, hand, gamePly, m_convTmp);

	// we have to map every time because the vector's buffer could have moved
	Eigen::Map<NNVector> mappedVec(&m_convTmp[0], 1, FeaturesConv::PositionInputDims);
#endif

#ifdef GIRAFFE_INPUT	
#ifdef DEEPPINK_REVERSE
	Position evalPos;
	StateInfo st;
	// ここの手番が先手なら局面は後手の後なので、評価のため盤面を入れ替える必要がある
	std::string evalSfen = (pos.side_to_move() == BLACK ? sfen_reverse(pos.sfen()) : pos.sfen());

	evalPos.set(evalSfen, &st, pos.this_thread());

	std::vector<float> m_convTmp;
	FeaturesConv::convert_pos_to_NN(evalPos, m_convTmp);
#else
	std::vector<float> m_convTmp;
	FeaturesConv::convert_pos_to_NN(pos, m_convTmp);
#endif
	/*for (int i = 0; i < m_convTmp.size(); ++i) {
		std::cout << m_convTmp[i] << " ";
		if (i % 7 == 6)
			std::cout << std::endl;
	}*/

	// we have to map every time because the vector's buffer could have moved
	Eigen::Map<NNVector> mappedVec(&m_convTmp[0], 1, m_convTmp.size());
#endif

	float annOut = m_mainAnn.forward_propagate_single(mappedVec);
	//std::cout << annOut << std::endl;
	
	Value nnRet = Value(int16_t(annOut * EvalFullScale));

	return nnRet;
}

NNMatrixRM ANNEvaluator::boards_to_feature_representation(const std::vector<std::string> &positions, const std::vector<FeaturesConv::FeatureDescription> &featureDescriptions)
{
	// SFEN数分の局面×特徴数の行列
	NNMatrixRM ret(positions.size(), featureDescriptions.size());

	// ompのスレッド数を決定
	ScopedOmpLimiter tlim(8);

#pragma omp parallel
	{
		std::vector<float> features; // each thread reuses a vector to avoid needless allocation/deallocation

									 // 局面から特徴を抽出して入力にする
									 // featuresの型がfloatであることに注意
#pragma omp for
		for (int64_t i = 0; i < positions.size(); ++i)
		{
			Position pos;
			StateInfo st;

			pos.set(positions[i], &st, Threads.main());
			FeaturesConv::convert_pos_to_NN(pos, features);

			assert(features.size() == featureDescriptions.size());

			ret.row(i) = Eigen::Map<NNMatrixRM>(&features[0], 1, static_cast<int64_t>(features.size()));
		}
	}

	return ret;
}

NNMatrixRM ANNEvaluator::boards_to_bitboard_feature(const std::vector<std::string> &positions)
{
	// SFEN数分の局面×特徴数の行列
	NNMatrixRM ret(positions.size(), FeaturesConv::PositionInputDims);

	// ompのスレッド数を決定
	ScopedOmpLimiter tlim(8);

#pragma omp parallel
	{
#pragma omp for
		for (int64_t i = 0; i < positions.size(); ++i)
		{
			float features[FeaturesConv::PositionInputDims] = { 0.0f };

			std::istringstream ss(positions[i]);
			std::string board, stm, hand, gamePly;

			//std::cout << positions[i] << std::endl;
			ss >> board >> stm >> hand >> gamePly;

			if (stm == "w") // 後手なら盤面と持ち駒反転
			{
				std::string boardRev;
				for (int64_t j = board.size() - 1; j >= 0; --j) // 反転させるため逆から
				{
					char token = board[j], plus = (j == 0 ? ' ' : board[j - 1]);
					if (std::isdigit(token) || token == '/') // 数字, '/'はそのまま.('+'は下で入れる)
						boardRev.push_back(token);

					else if (isupper(token)) // 大文字は小文字に
					{
						if (plus == '+') { boardRev.push_back(plus); boardRev.push_back(tolower(token)); }
						else boardRev.push_back(tolower(token));
					}
					else if (islower(token)) // 小文字は大文字に
					{
						if (plus == '+') { boardRev.push_back(plus); boardRev.push_back(toupper(token)); }
						else boardRev.push_back(toupper(token));
					}
				}
				board = boardRev;

				std::string handRev;
				for (int64_t j = 0; j < hand.size(); ++j) { // 持ち駒は前からでOK
					char token = hand[j];
					if (std::isdigit(token) || token == '-') handRev.push_back(token);
					else if (isupper(token)) handRev.push_back(tolower(token));
					else if (islower(token)) handRev.push_back(toupper(token));
				}
				hand = handRev;
			}

			//std::string sfen = board + " " + stm + " " + hand + " " + gamePly;
			//std::cout << sfen << std::endl;

			FeaturesConv::convert_sfen_to_BB(board, stm, hand, gamePly, features);

			ret.row(i) = Eigen::Map<NNMatrixRM>(&features[0], 1, FeaturesConv::PositionInputDims);
		}
	}

	return ret;
}

NNMatrixRM ANNEvaluator::compute_error_derivatives(
	const NNMatrixRM &predictions,
	const NNMatrixRM &targets,
	const NNMatrixRM &finalLayerActivations,
	float positiveWeight,
	float negativeWeight)
{
	int64_t numExamples = predictions.rows();

	NNMatrixRM ret(numExamples, 1);

	// this takes care of everything except the dtanh(act)/dz term, which we can't really vectorize
	ret = (targets - predictions) * -1.0f;

	// derivative of tanh is 1-tanh^2(x)
	for (int64_t i = 0; i < numExamples; ++i)
	{
		float tanhx = tanh(finalLayerActivations(i, 0));
		ret(i, 0) *= 1.0f - tanhx * tanhx;

		if (ret(i, 0) > 0.0f)
			ret(i, 0) *= positiveWeight;
		else
			ret(i, 0) *= negativeWeight;
	}

	return ret;
}

template <ActivationFunc ACTF, ActivationFunc ACTFLast>
NNMatrixRM ANNEvaluator::compute_loss_derivatives(
	const NNMatrixRM &loss,
	const NNMatrixRM &finalLayerActivations,
	float positiveWeight,
	float negativeWeight)
{
	int64_t numExamples = loss.rows();

	NNMatrixRM ret(numExamples, 1);

	// this takes care of everything except the dtanh(act)/dz term, which we can't really vectorize
	ret = loss * -1.0f;

	// derivative of linear is x
	if (ACTFLast == Linear)
	{
		for (int64_t i = 0; i < numExamples; ++i)
		{
			ret(i, 0) *= finalLayerActivations(i, 0); // ここ違う可能性高し

			if (ret(i, 0) > 0.0f)
				ret(i, 0) *= positiveWeight;
			else
				ret(i, 0) *= negativeWeight;
		}
	}

	// derivative of tanh is 1-tanh^2(x)
	else if (ACTFLast == Tanh)
	{
		for (int64_t i = 0; i < numExamples; ++i)
		{
			float tanhx = tanh(finalLayerActivations(i, 0));
			ret(i, 0) *= 1.0f - tanhx * tanhx;

			if (ret(i, 0) > 0.0f)
				ret(i, 0) *= positiveWeight;
			else
				ret(i, 0) *= negativeWeight;
		}
	}

	return ret;
}