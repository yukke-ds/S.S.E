#ifndef FEATURES_CONV_H_INCLUDED
#define FEATURES_CONV_H_INCLUDED

#include <vector>
#include <string>
#include <sstream>
#include <set>

#include "../Eigen/Dense"

#include "ann.h"
#include "../position.h"

//#define BB_INPUT
#define GIRAFFE_INPUT
//#define DEEPPINK_REVERSE

#if defined(BB_INPUT) && defined(GIRAFFE_INPUT)
#error Only select one kind of input feature!
#elif !defined(BB_INPUT) && !defined(GIRAFFE_INPUT)
#error Must select one kind of input feature!
#endif

namespace FeaturesConv
{
	const int64_t BoardInputDims = (14 * 2 * 81);
	const int64_t HandInputDims = (7 * 2);
	const int64_t PositionInputDims = BoardInputDims + HandInputDims/* + 1 Žè”Ô */;

	struct FeatureDescription
	{
		enum FeatureType
		{
			FeatureType_global, // global features are things like side to move, and material counts, and piece lists
			FeatureType_pos // property of a square
		};

		FeatureType featureType;

		// fields for global and pos features
		int32_t group;

		// fields for pos features
		Square sq;
	};

	// convert to NN input format
	// T can either be float (to get actual values) or
	// FeatureDescription (to get feature descriptions)
	template <typename T>
	void convert_pos_to_NN(const Position &pos, std::vector<T> &ret);

	void convert_sfen_to_BB(const std::string& board, const std::string& stm, const std::string& hand, const std::string& gamePly, float ret[PositionInputDims]);

} // namespace FeaturesConv

#endif // ifndef FEATURES_CONV_H_INCLUDED