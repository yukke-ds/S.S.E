#ifndef LEARN_H_INCLUDED
#define LEARN_H_INCLUDED

#include <cstdint>
#include <string>
#include <fstream>
#include <iostream>

#include "ann/ann_evaluator.h"
#include "search.h"
#include "types.h"

#include "Eigen/Core"

//#define DISABLE_TT_PROBE

struct SearchResult
{
	Value score;
	std::vector<Move> pv;
};

namespace Learn
{
	
const static int64_t NumIterations = 100000;
const static float TDLambda = 0.7f; // this is discount due to credit assignment uncertainty
const static float AbsLambda = 0.995f; // this is discount to encourage progress, and account for the snowball effect
const static int64_t HalfMovesToMake = 12;
const static size_t PositionsFirstBatch = 1000000;
const static size_t PositionsPerBatch = 1000;
const static float MaxError = 1.0f;
const static Depth SearchDepth = 3 * ONE_PLY;
const static float LearningRate = 1.0f;
const static float LearningRateSGD = 0.01f;
const static int64_t EvaluatorSerializeInterval = 10;
const static int64_t IterationPrintInterval = 1;
const static int64_t BoundTrainingEpochs = 10;

inline float operator*(Value d1, float d2) { return float(float(d1) * d2); }
inline float normalization(Value score) { return tanh(static_cast<int>(score) / 600.0); }

void make_sfens();
void make_triplet_sfens();

void TD_leaf();
void supervised();

SearchResult search(Position& pos, Depth depth, ANNEvaluator* evaluator);

} // namespace Learn

class Stat
{
public:
	Stat() : m_sum(0.0f), m_count(0) {}

	void reset() 
	{
		m_sum = 0.0f;
		m_count = 0;
	}

	void add_number(float x) 
	{
		m_sum += x;
		++m_count;
	}

	float get_average() 
	{
		if (m_count != 0)
			return m_sum / static_cast<float>(m_count);
		else
			return 0.0f;
	}

private:
	float m_sum;
	uint64_t m_count;
};

#endif // ifndef LEARN_H_INCLUDED