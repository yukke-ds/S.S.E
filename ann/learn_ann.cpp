#include "learn_ann.h"

#include <iostream>
#include <fstream>
#include <random>
#include <thread>
#include <vector>
#include <chrono>
#include <stdexcept>
#include <queue>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <type_traits>
#include <omp.h>

#include <cmath>

#include "ann.h"
#include "features_conv.h"

namespace
{

const size_t MaxBatchSize = 256;
const size_t MaxIterationsPerCheck = 500000 / MaxBatchSize;

typedef std::vector<int32_t> Group;

struct Rows
{
	Rows() {}
	Rows(int64_t begin, int64_t num) : begin(begin), num(num) {}

	int64_t begin;
	int64_t num;
};

template <typename Derived1>
void split_dataset(
	const Eigen::MatrixBase<Derived1> &x,
	Rows &train,
	Rows &val,
	Rows &test)
{
	size_t numExamples = x.rows();

	const float testRatio = 0.2f;
	const size_t MaxTest = 5000;
	const float valRatio = 0.2f;
	const size_t MaxVal = 5000;

	size_t testSize = std::min<size_t>(MaxTest, numExamples * testRatio);
	size_t valSize = std::min<size_t>(MaxVal, numExamples * valRatio);
	size_t trainSize = numExamples - testSize - valSize;

	test = Rows(0, testSize);
	val = Rows(testSize, valSize);
	train = Rows(testSize + valSize, trainSize);
}

template <typename T, typename Derived1, typename Derived2>
void train(
	T &nn,
	int64_t epochs,
	Eigen::MatrixBase<Derived1> &xTrain,
	Eigen::MatrixBase<Derived2> &yTrain,
	Eigen::MatrixBase<Derived1> &xVal,
	Eigen::MatrixBase<Derived2> &yVal,
	Eigen::MatrixBase<Derived1> &/*xTest*/,
	Eigen::MatrixBase<Derived2> &/*yTest*/)
{
	// �w�K�C�e���[�V������
	size_t iteration = 0;

	// ���t�f�[�^�Əo�͒l�Ƃ̌덷�̘a
	float trainingErrorAccum = 0.0f;

	// �w�K���Ԃ��v��
	std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();

	// �����̏��Ȃ�NN���x�X�gNN
	T bestNet = nn;
	FP bestValScore = std::numeric_limits<FP>::max(); // this is updated every time val score improves

	bool done = false;

	// �w�K�f�[�^������o�b�`�T�C�Y�����߂�
	size_t NumBatches = xTrain.rows() / MaxBatchSize;

	if ((xTrain.rows() % MaxBatchSize) != 0)
		++NumBatches;

	// �w�K�̐i�s�x�݂����Ȃ���
	int64_t epoch = 0;

	// we want to check at least once per epoch
	size_t iterationsPerCheck = std::min(MaxIterationsPerCheck, NumBatches);

	size_t examplesSeen = 0;

	// epoch��1�����Ȃ瑱����
	while (!done && epoch < epochs)
	{
		size_t batchNum = iteration % NumBatches;
		size_t begin = batchNum * MaxBatchSize;
		size_t batchSize = std::min(MaxBatchSize, xTrain.rows() - begin);

		examplesSeen += batchSize;

		// ���ꂪ1�ɂȂ�����w�K�I���
		epoch = examplesSeen / xTrain.rows();

		// ���t�Əo�͂̌덷�̘a���W�v���Ă���
		trainingErrorAccum += nn.train_GDM(
			xTrain.block(begin, 0, batchSize, xTrain.cols()),
			yTrain.block(begin, 0, batchSize, yTrain.cols()),
			1.0f,
			0.000001f);

		// epoch��1��͏����o�͂��Ă�������
		if ((iteration % iterationsPerCheck) == 0)
		{
			NNMatrix pred = nn.forward_propagate_fast(xVal);
			NNMatrix errors = nn.error_func(pred, yVal);

			FP valScore = errors.sum() / xVal.rows();
			if (valScore < bestValScore)
			{
				bestValScore = valScore;
				bestNet = nn;
			}

			std::cout << "Iteration: " << iteration << ", ";
			std::cout << "Epoch: " << epoch << ", ";
			std::cout << "Val: " << valScore << ", ";
			std::cout << "Train: " << (trainingErrorAccum / std::min(iteration + 1, iterationsPerCheck)) << ", ";

			std::chrono::seconds t = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - startTime);
			std::cout << "Time: " << (static_cast<float>(t.count()) / 60.0f) << " minutes, ";
			std::cout << "Best Val: " << bestValScore << ", ";
			std::cout << "Sparsity: " << nn.get_sparsity() << std::endl;

			// �����Ō덷�̘a��0�ɖ߂��݂���
			trainingErrorAccum = 0.0f;
		}

		++iteration;
	}

	// ANN�����ʂ̗ǂ��������̂ɒu��������
	nn = bestNet;
}

struct LayerDescription
{
	size_t layerSize;
	std::vector<Eigen::Triplet<float> > connections;

	LayerDescription() : layerSize(0) {}
};

void add_single_nodes_group(
	LayerDescription &layerDescription,
	const Group &groupIn,
	float nodeCountMultiplier)
{
	size_t nodesInGroup = groupIn.size();
	size_t nodesForThisGroup = ceil(nodesInGroup * nodeCountMultiplier);

	for (size_t i = 0; i < nodesForThisGroup; ++i)
	{
		for (auto feature : groupIn)
			layerDescription.connections.push_back(Eigen::Triplet<float>(feature, layerDescription.layerSize, 1.0f));

		++layerDescription.layerSize;
	}
}

void analyze_feature_descriptions(const std::vector<FeaturesConv::FeatureDescription> &featureDescriptions,
	Group &globalGroup, /* global group does not include group 0! */
	Group &squareGroup,
	Group &group0)
{
	// first we make global feature groups
	for (int32_t featureNum = 0; featureNum < featureDescriptions.size(); ++featureNum)
	{
		auto &fd = featureDescriptions[featureNum];

		if (fd.featureType == FeaturesConv::FeatureDescription::FeatureType_global)
		{
			if (fd.group == 0)
				group0.push_back(featureNum);
			else
				globalGroup.push_back(featureNum);
		}
		else if (fd.featureType == FeaturesConv::FeatureDescription::FeatureType_pos)
			squareGroup.push_back(featureNum);
	}

	assert(group0.size() > 5 && group0.size() < 40);
}

} // namespace

namespace LearnAnn
{

template <>
EvalNet build_eval_net<TD_LEAF>(int64_t inputDims, int64_t outputDims, bool smallNet)
{
	// �e���ԑw�̑f�q�����i�[����.�T�C�Y���w���ɂ�����
	std::vector<size_t> layerSizes;

	// 2����float�s��(Triplet=row,col,value)
	// �T�C�Y0�Ȃ�S����,�@�����łȂ���΃X�p�[�X���C��
	std::vector<std::vector<Eigen::Triplet<float> > > connMatrices;

	Group globalGroup;
	Group squareGroup;
	Group group0;

	// get feature descriptions
	std::vector<FeaturesConv::FeatureDescription> featureDescriptions;
	const Position dummyPos;
	FeaturesConv::convert_pos_to_NN(dummyPos, featureDescriptions);

	analyze_feature_descriptions(featureDescriptions, globalGroup, squareGroup, group0);

	// InputLayer{Group0(�i�s�x��), GlobalGroup(������), SquareGroup�܂�}�C
	// SecondLayer�COutputLayer���ꂼ��̓��o�͂̐��C�w�̐����v�Z���z��ɂ܂Ƃ߂�D
	LayerDescription layer0;

	float nodeCountMultiplier = (smallNet ? 0.1f : 0.05f);

	// 1. Input Layer
	// first we add the mixed global group
	add_single_nodes_group(layer0, globalGroup, nodeCountMultiplier);

	// mixed square group
	add_single_nodes_group(layer0, squareGroup, nodeCountMultiplier);

	// pass through group 0 (this contains game phase information)
	add_single_nodes_group(layer0, group0, 1.0f);

	layerSizes.push_back(layer0.layerSize);
	connMatrices.push_back(layer0.connections);

	// 2. Hidden Layer
	// in the second layer, we just fully connect everything
	// this is the output of the second to last layer in the net
	layerSizes.push_back(SQUARE_NB);
	connMatrices.push_back(std::vector<Eigen::Triplet<float> >());

	// 3. Output Layer
	// fully connected output layer
	connMatrices.push_back(std::vector<Eigen::Triplet<float> >());

	std::cout << "group0:" << group0.size()
		<< " GlobalGroup:" << globalGroup.size()
		<< " SquareGroup:" << squareGroup.size() << std::endl;
	std::cout << "layerSizes:" << layerSizes.size()
		<< " -> Hidden1:" << layerSizes[0] 
		<< ", Hidden2:" << layerSizes[1] << std::endl;
	std::cout << "connMatrices:" << connMatrices.size() 
		<< " -> Input-Hidden1:" << connMatrices[0].size() 
		<< ", Hidden1-Hidden2:" << connMatrices[1].size()
		<< ", Hidden2-Output:" << connMatrices[2].size() << std::endl;
		

	return EvalNet(inputDims, outputDims, layerSizes, connMatrices);
}

template <>
EvalNet build_eval_net<SL>(int64_t inputDims, int64_t outputDims, bool smallNet)
{
	const size_t LayerSize = 2;
	const size_t Elements = 300;

	// �e���ԑw�̑f�q�����i�[����.�T�C�Y���w���ɂ�����
	std::vector<size_t> layerSizes;
	for (size_t i = 0; i < LayerSize; ++i)
		layerSizes.push_back(Elements);

	 // 2����float�s��(Triplet=row,col,value)
	 // �T�C�Y0�Ȃ�S����,�@�����łȂ���΃X�p�[�X���C��
	std::vector<std::vector<Eigen::Triplet<float> > > connMatrices;
	for (size_t i = 0; i < LayerSize + 1; ++i)
		connMatrices.push_back(std::vector<Eigen::Triplet<float> >());

	std::cout << "layerSizes:" << layerSizes.size()
		<< ", Hidden1:" << layerSizes[0] 
		<< ", Hidden2:" << layerSizes[1] << std::endl;
	std::cout << "connMatrices:" << connMatrices.size()
		<< " -> Input-Hidden1:" << connMatrices[0].size()
		<< ", Hidden1-Hidden2:" << connMatrices[1].size()
		<< ", Hidden2-Output:" << connMatrices[2].size()
		<< std::endl;
	
	return EvalNet(inputDims, outputDims, layerSizes, connMatrices);
}

template <typename Derived1, typename Derived2>
void train_ANN(
	const Eigen::MatrixBase<Derived1> &x,
	const Eigen::MatrixBase<Derived2> &y,
	EvalNet &nn,
	int64_t epochs)
{
	Rows trainRows, valRows, testRows;

	// �w�K�p�C�]���p�C�e�X�g�p�ɕ�����(���͕]���p=�e�X�g�p)
	split_dataset(x, trainRows, valRows, testRows);

	auto xTrain = x.block(trainRows.begin, 0, trainRows.num, x.cols());
	auto yTrain = y.block(trainRows.begin, 0, trainRows.num, y.cols());
	auto xVal = x.block(valRows.begin, 0, valRows.num, x.cols());
	auto yVal = y.block(valRows.begin, 0, valRows.num, y.cols());
	auto xTest = x.block(testRows.begin, 0, testRows.num, x.cols());
	auto yTest = y.block(testRows.begin, 0, testRows.num, y.cols());

	std::cout << "Train: " << xTrain.rows() << std::endl;
	std::cout << "Val: " << xVal.rows() << std::endl;
	std::cout << "Test: " << xTest.rows() << std::endl;
	std::cout << "Features: " << xTrain.cols() << std::endl;

	std::cout << "Beginning training..." << std::endl;
	train(nn, epochs, xTrain, yTrain, xVal, yVal, xTest, yTest);
}

// here we have to list all instantiations used (except for in this file)
// �e���v���[�g�̖����I���ꉻ
template void train_ANN<NNMatrixRM, NNMatrixRM>(const Eigen::MatrixBase<NNMatrixRM>&, const Eigen::MatrixBase<NNMatrixRM>&, EvalNet &, int64_t);

} // namespace LearnAnn