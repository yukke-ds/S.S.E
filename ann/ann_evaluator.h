#ifndef ANN_EVALUATOR_H_INCLUDED
#define ANN_EVALUATOR_H_INCLUDED

#include <cmath>
#include <vector>
#include <string>

#include "ann.h"
#include "features_conv.h"
#include "../types.h"

#include "learn_ann.h"

class ANNEvaluator
{
public:

	constexpr static float EvalFullScale = 10000.0f;

	float unscale(Value v) {
		float ret = float(v) / EvalFullScale;

		ret = std::max(ret, -1.0f);
		ret = std::min(ret, 1.0f);

		return ret;
	}

	// ANN related
	template <typename Method T>
	void build_ANN(int64_t inputDims);
	NNMatrixRM boards_to_feature_representation(const std::vector<std::string> &positions, const std::vector<FeaturesConv::FeatureDescription> &featureDescriptions);
	NNMatrixRM boards_to_bitboard_feature(const std::vector<std::string> &positions);

	// for input/output ANN
	void serialize(std::ostream &os);
	void deserialize(std::istream &is);
	void deserialize_py(std::istream &is);

	// for ANN Learning
	void train(const std::vector<std::string> &positions, const NNMatrixRM &y, const std::vector<FeaturesConv::FeatureDescription> &featureDescriptions, float learningRate);
	void train_loop(const std::vector<std::string> &positions, const NNMatrixRM &y, int64_t epochs, const std::vector<FeaturesConv::FeatureDescription> &featureDescriptions);
	float forward_triplet(const std::vector<std::string> &observed, const std::vector<std::string> &parent, const std::vector<std::string> &random, const std::vector<FeaturesConv::FeatureDescription> &featureDescriptions, float learningRate);

	// evaluation
	Value evaluate(const Position& pos);

private:

	NNMatrixRM compute_error_derivatives(
		const NNMatrixRM &predictions,
		const NNMatrixRM &targets,
		const NNMatrixRM &finalLayerActivations,
		float positiveWeight,
		float negativeWeight);

	template <ActivationFunc ACTF, ActivationFunc ACTFLast>
	NNMatrixRM compute_loss_derivatives(
		const NNMatrixRM &loss,
		const NNMatrixRM &finalLayerActivations,
		float positiveWeight,
		float negativeWeight);

	EvalNet m_mainAnn;
	EvalNet m_ubAnn;
	EvalNet m_lbAnn;
};

extern ANNEvaluator AnnEvaluator;

#endif // ifndef ANN_EVALUATOR_H_INCLUDED