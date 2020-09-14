#include "learn.h"

#include <fstream>
#include <vector>

#include "ann/features_conv.h"
#include "ann/ann_evaluator.h"
#include "common/matrix_ops.h"
#include "evaluate.h"
#include "gamedb.h"
#include "misc.h"
#include "movepick.h"
#include "search.h"
#include "thread.h"
#include "usi.h"

namespace
{

using namespace Eval;
using namespace Learn;

void init()
{
	std::cout << "# Using " 
		<< omp_get_max_threads() << " OpenMP thread(s), " 
		<< (size_t)Options["Threads"] << " std::thread thread(s), "
		<< "ANNEvaluator " << std::boolalpha << (bool)Options["ANNEvaluator"] 
		<< " #" << std::endl;

	Eigen::initParallel();

	// set Eigen to use 1 thread because we are doing OpenMP here
	Eigen::setNbThreads(1);

	// disable nested parallelism since we don't need it, and disabling it
	// makes managing number of threads easier
	omp_set_nested(0);
}

std::string get_filename(int64_t iter)
{
	std::stringstream filenameSs;

	filenameSs << "TrainingResults/eval" << iter << ".net";

	return filenameSs.str();
}

bool file_exists(const std::string &filename)
{
	std::ifstream is(filename);

	return is.good();
}

}

namespace Learn
{

// Deep Neural NetworkでTDLeaf(λ)
void TD_leaf()
{
	init();

	std::cout << "Starting TDLeaf training..." << std::endl;

	// 棋譜ファイルのストリーム
	const std::string positionsFileName = "TrainingData/sfens.txt";
	std::ifstream positionsFile(positionsFileName);

	assert(positionsFile);

	// these are the root positions for training (they don't change)
	// ルート局面の配列
	std::vector<std::string> rootPositions;

	std::string sfen;

	std::cout << "Reading SFENs..." << std::endl;

	// 棋譜ファイルから1行ずつFENを読み込む
	while (std::getline(positionsFile, sfen))
	{
		// ルート局面の配列に入れる
		rootPositions.push_back(sfen);
		assert(sfen != "");
	}

	std::cout << "Positions read: " << rootPositions.size() << std::endl;

	// these are the leaf positions used in training
	// they are initialized to root positions, but will change in second iteration
	// 学習に用いるリーフ局面の配列.サイズはバッチサイズ分（ただし，1イテレートと2イテレート以降で異なる）
	std::vector<std::string> trainingPositions(PositionsPerBatch);

	// 教師データの配列
	NNMatrixRM trainingTargets(trainingPositions.size(), 1);

	// ANN(Artificial Neural Network)の局面評価クラス
	ANNEvaluator annEvaluator;

	// 特徴を格納する配列．サイズが特徴の数と同じになる
	std::vector<FeaturesConv::FeatureDescription> featureDescriptions;

	// 多分初期化みたいなもの
	const Position dummyPos;
	FeaturesConv::convert_pos_to_NN(dummyPos, featureDescriptions);

	// 学習イテレータ
	int64_t iteration = 0;

	// try to read existing dumps to continue where we left off
	// 学習途中のファイルがあれば，それをANNに適用して学習を開始する
	for (; iteration < NumIterations; iteration += EvaluatorSerializeInterval)
	{
		std::string filename = get_filename(iteration);

		if (!file_exists(filename))
			break;
	}

	if (iteration > 0)
	{
		iteration -= EvaluatorSerializeInterval;

		std::string filename = get_filename(iteration);
		std::ifstream dump(filename);
		annEvaluator.deserialize(dump); // dumpをNNモデルにデシリアライズ
	}

	if (iteration != 0)
		std::cout << "Continuing from iteration " << iteration << std::endl;

	// 学習全体でかかる時間を計測
	TimePoint timeStart = now();

	// TD誤差とかの統計
	Stat errorStat;
	std::ofstream stats("TrainingResults/Stats.txt");

	// 学習イテレーション開始
	for (; iteration < NumIterations; ++iteration)
	{
		// イテレーション開始から時間を計測
		TimePoint iterationStart = now();

		// first we generate new labels
		// if this is the first iteration, we use static material labels
		// まず，学習対象となる局面とその評価値をそれぞれ求めて格納する(学習用，テスト用が含まれる)
		// 最初のイテレーションでは，ルート局面の駒価値を用いる
		if (iteration == 0)
		{
			std::cout << "Labelling using static evaluation..." << std::endl;

			// 最初のバッチサイズ分に学習局面数もサイズ変更．教師データも同様
			trainingPositions.resize(PositionsFirstBatch);
			trainingTargets.resize(trainingPositions.size(), 1);

			// それぞれのルート局面とその駒価値を求めて格納．評価値は-1~1の間にする
			#pragma omp parallel for
			for (int64_t i = 0; i < PositionsFirstBatch; ++i)
			{
				Position pos;
				StateInfo st;
				pos.set(rootPositions[i], &st, Threads.main());
				trainingPositions[i] = rootPositions[i];
				trainingTargets(i, 0) = normalization(material(pos));
			}
		}
		
		// 2回目のイテレーション以降はTD-leaf(λ)を用いる
		else
		{
			// 処理する局面数を記録
			size_t positionsProcessed = 0;

			// 2回目以降のバッチサイズ分に学習局面数，教師データをサイズ変更
			trainingPositions.resize(PositionsPerBatch);
			trainingTargets.resize(trainingPositions.size(), 1);

			// 0~局面数の間で乱数を作るための機能を定義
			auto rng = gRd.make_mt();
			auto positionDist = std::uniform_int_distribution<size_t>(0, rootPositions.size() - 1);
			auto positionDrawFunc = std::bind(positionDist, rng);

			// スレッド並列
			std::vector<std::thread> threads;
			size_t numThreads = (size_t)Options["Threads"];
			Mutex mutex;

			for (size_t i = 0; i < numThreads; ++i)
			{
				threads.emplace_back([&](size_t id) {

				ANNEvaluator thread_annEvaluator = annEvaluator;
				
				size_t sepStart = PositionsPerBatch / numThreads * id + std::min(PositionsPerBatch % numThreads, id);
				size_t sepEnd	= PositionsPerBatch / numThreads * (id+1) + std::min(PositionsPerBatch % numThreads, id+1);

				// バッチサイズ分のループ
				for (size_t j = sepStart; j < sepEnd; ++j)
				{
					// ランダムに選んだFENをrootPosにセットする
					Position rootPos;
					StateInfo stRoot;

					rootPos.set(rootPositions[positionDrawFunc()], &stRoot, Threads[id]);

					if (!rootPos.pos_is_ok())
						continue;

					// make 1 random move
					// it's very important that we make an odd number of moves, so that if the move is something stupid, the
					// opponent can take advantage of it (and we will learn that this position is bad) before we have a chance to fix it
					// 【重要】ランダムな指し手を1手指す．理由については↑の英訳参照．論文P22の真ん中下．
					std::vector<Move> ml;
					for (const auto& m : MoveList<LEGAL>(rootPos))
						ml.push_back(m);

					if (ml.size() == 0)
						continue;

					auto movePickerDist = std::uniform_int_distribution<size_t>(0, ml.size() - 1);

					Move move = ml[movePickerDist(rng)];
					assert(rootPos.legal(move) || rootPos.pseudo_legal(move));

					StateInfo stRand;
					rootPos.do_move(move, stRand, rootPos.gives_check(move));

					// PVを探索
					SearchResult rootResult = Learn::search(rootPos, SearchDepth, &thread_annEvaluator);

					if (rootResult.score == -VALUE_MATE && rootResult.pv.size() == 0)
						continue;

					// PVのリーフ局面を作る
					Position leafPos = rootPos;
					std::vector<StateInfo> stToLeaf(rootResult.pv.size());
					for (size_t k = 0; k < rootResult.pv.size(); ++k)
					{
						Move m = rootResult.pv[k];
						assert(leafPos.legal(m) || leafPos.pseudo_legal(m));
						leafPos.do_move(m, stToLeaf[k], leafPos.gives_check(m));
					}

					// そこのANN評価値を求める
					Value leafScore = thread_annEvaluator.evaluate(leafPos);

					// αβ探索でルートまで上がってきた評価値を先手から見たものにする
					Value rootScoreBlack = rootResult.score * (rootPos.side_to_move() == BLACK ? 1 : -1);

					// リーフ局面のSFENの学習局面として保存
					trainingPositions[j] = leafPos.sfen();

					// -1.0~+1.0の間に正規化
					float leafScoreUnscaled = thread_annEvaluator.unscale(leafScore);

					// PVが0より大きくて，リーフ局面の評価値とルート局面に上がってきた評価値が等しい(詰みでない)なら
					if (rootResult.pv.size() > 0 && (leafScore == rootScoreBlack))
					{
						// PVの1手目を指す
						Move pv = rootResult.pv[0];
						StateInfo stPV;

						rootPos.do_move(pv, stPV, rootPos.gives_check(pv));

						// now we compute the error by making a few moves
						// 数手指していく中で損失を計算する
						float accumulatedError = 0.0f;
						float lastScore = leafScoreUnscaled;
						float tdDiscount = 1.0f;
						float absoluteDiscount = AbsLambda;

						// TD-leaf(λ)のλ収益を深さ12まで計算する
						StateInfo stForLambda[HalfMovesToMake];
						for (size_t m = 0; m < HalfMovesToMake; ++m)
						{
							// 現局面からPVとそのルートに上がってきた評価値を調べる
							SearchResult result = Learn::search(rootPos, SearchDepth, &thread_annEvaluator);

							// 先手からの評価値にしてスケールする
							float scoreBlackUnscaled = thread_annEvaluator.unscale(result.score * (rootPos.side_to_move() == BLACK ? 1 : -1)) * absoluteDiscount; // γVπ(St+1)

							// 割引率を更新
							absoluteDiscount *= AbsLambda;

							// compute error contribution (only if same side)
							// 損失計算(同じ手番のみ)
							if (m % 2 == 1)
							{
								// TD誤差の蓄積 → λj-t(γVπ(St+1) - Vπ(St))
								accumulatedError += tdDiscount * (scoreBlackUnscaled - lastScore);
								lastScore = scoreBlackUnscaled;
								tdDiscount *= TDLambda; // λを更新
							}

							// どちらかが詰みならそこまでで終わり
							if (result.pv.size() == 0)
								break;

							rootPos.do_move(result.pv[0], stForLambda[m], rootPos.gives_check(result.pv[0]));
						}

						// λ収益の絶対値
						float absError = fabs(accumulatedError);

						// 損失の統計
						mutex.lock();
						errorStat.add_number(absError);
						mutex.unlock();

						// -1.0~+1.0に調整
						accumulatedError = std::max(accumulatedError, -MaxError);
						accumulatedError = std::min(accumulatedError, MaxError);

						trainingTargets(j, 0) = leafScoreUnscaled + LearningRate * accumulatedError; // α*(Σj=t,N-1 (∇V*Σ(accumuratedError)))
					}
					else
						// if PV is empty or leaf score is not the same as search score, this is an end position, and we don't need to train it
						trainingTargets(j, 0) = thread_annEvaluator.unscale(leafScore);

					mutex.lock();
					++positionsProcessed;
					mutex.unlock();

				}}, i);
			}

			// スレッドの終了を待機
			for (auto& t : threads)
				t.join();
		}

		// 最初のイテレーション
		if (iteration == 0)
		{
			// 最初はANNを作る
			annEvaluator.build_ANN<TD_LEAF>(featureDescriptions.size());

			// 学習ループへ
			annEvaluator.train_loop(trainingPositions, trainingTargets, 1, featureDescriptions);
		}

		// 2回目以降
		else
			annEvaluator.train(trainingPositions, trainingTargets, featureDescriptions, LearningRateSGD);

		// 色々な出力
		if ((iteration % EvaluatorSerializeInterval) == 0)
		{
			auto mt = gRd.make_mt();
			std::shuffle(rootPositions.begin(), rootPositions.end(), mt);

			std::cout << "Serializing..." << std::endl;

			std::ofstream annOut(get_filename(iteration));

			annEvaluator.serialize(annOut);
		}

		if ((iteration % IterationPrintInterval) == 0)
		{
			std::cout << "Iteration " << iteration << ". ";
			std::cout << "Time: " << (now() - timeStart) / 1000 << " seconds. ";
			std::cout << "Last Iteration took: " << (now() - iterationStart) / 1000 << " seconds. ";

			float tdError = errorStat.get_average();
			std::cout << "TD Error: " << tdError << ". ";
			stats << iteration << "," << tdError << std::endl;

			errorStat.reset();

			std::cout << std::endl;
		}
	}

	std::cout << "Congratulations! Learning is successfully finished!" << std::endl;
}

void supervised()
{
	init();

	std::cout << "Starting Supervised Learning..." << std::endl;

	// 棋譜ファイルのストリーム
	const std::string observedPositionsFileName = "TrainingData/observedSFENs.txt";
	const std::string parentPositionsFileName = "TrainingData/parentSFENs.txt";
	const std::string randomPositionsFileName = "TrainingData/randomSFENs.txt";

	std::ifstream observedPositionsFile(observedPositionsFileName);
	std::ifstream parentPositionsFile(parentPositionsFileName);
	std::ifstream randomPositionsFile(randomPositionsFileName);

	assert(observedPositionsFile);
	assert(parentPositionsFile);
	assert(randomPositionsFile);

	// 教師局面, 親局面, ランダム局面の配列
	std::vector<std::string> observedPositions;
	std::vector<std::string> parentPositions;
	std::vector<std::string> randomPositions;

	std::string sfenObserved, sfenParent, sfenRandom;

	std::cout << "Reading SFENs..." << std::endl;

	// 棋譜ファイルから1行ずつSFENを読み込んで, 配列に入れる
	while (std::getline(observedPositionsFile, sfenObserved)
		&& std::getline(parentPositionsFile, sfenParent)
		&& std::getline(randomPositionsFile, sfenRandom)) 
	{
		observedPositions.push_back(sfenObserved);
		parentPositions.push_back(sfenParent);
		randomPositions.push_back(sfenRandom);
		assert(sfenObserved != "" && sfenParent != "" && sfenRandom != "");
	}

	std::cout << "Observed positions read: " << observedPositions.size() << std::endl;
	std::cout << "Parent positions read: " << parentPositions.size() << std::endl;
	std::cout << "Random positions read: " << randomPositions.size() << std::endl;

	// ANN(Artificial Neural Network)の局面評価クラス
	ANNEvaluator annEvaluator;

	// 特徴を格納する配列．サイズが特徴の数と同じになる
	std::vector<FeaturesConv::FeatureDescription> featureDescriptions;

	// 多分初期化みたいなもの
	const Position dummyPos;
	FeaturesConv::convert_pos_to_NN(dummyPos, featureDescriptions);

	// 学習イテレータ
	int64_t iteration = 0;

	// try to read existing dumps to continue where we left off
	// 学習途中のファイルがあれば，それをANNに適用して学習を開始する
	for (; iteration < NumIterations; iteration += EvaluatorSerializeInterval)
	{
		std::string filename = get_filename(iteration);

		if (!file_exists(filename))
			break;
	}

	if (iteration > 0)
	{
		iteration -= EvaluatorSerializeInterval;

		std::string filename = get_filename(iteration);
		std::ifstream dump(filename);
		annEvaluator.deserialize(dump); // dumpをNNモデルにデシリアライズ
	}

	if (iteration != 0)
		std::cout << "Continuing from iteration " << iteration << std::endl;

	// 学習全体でかかる時間を計測
	TimePoint timeStart = now();

	// 損失とかの統計
	Stat errorStat;
	std::ofstream stats("TrainingResults/supervisedStats.txt");

	// ANN構築
#ifdef BB_INPUT
	annEvaluator.build_ANN<SL>(FeaturesConv::PositionInputDims);
#endif
#ifdef GIRAFFE_INPUT
	annEvaluator.build_ANN<SL>(featureDescriptions.size());
#endif
	// 学習イテレーション開始
	for (; iteration < NumIterations; ++iteration)
	{
		// イテレーション開始から時間を計測
		TimePoint iterationStart = now();

		// 局面をシャッフル(observed-parent-randomの関係は崩さない)
		std::random_device seedGen;
		auto seed = seedGen();
		{
			std::mt19937 mt(seed);
			std::shuffle(observedPositions.begin(), observedPositions.end(), mt);
		}
		{
			std::mt19937 mt(seed);
			std::shuffle(parentPositions.begin(), parentPositions.end(), mt);
		}
		{
			std::mt19937 mt(seed);
			std::shuffle(randomPositions.begin(), randomPositions.end(), mt);
		}

		// 学習に使用する局面の配列
		std::vector<std::string> usingObserved(PositionsPerBatch);
		std::vector<std::string> usingParent(PositionsPerBatch);
		std::vector<std::string> usingRandom(PositionsPerBatch);

		std::copy(observedPositions.begin(), observedPositions.begin() + PositionsPerBatch, usingObserved.begin());
		std::copy(parentPositions.begin(), parentPositions.begin() + PositionsPerBatch, usingParent.begin());
		std::copy(randomPositions.begin(), randomPositions.begin() + PositionsPerBatch, usingRandom.begin());

		/*std::cout << "# Input SFENs #" << std::endl;
		for (size_t i = 0; i < PositionsPerBatch; ++i) {
			std::cout << "[" << i << "]" << std::endl;
			std::cout << usingObserved[i] << std::endl;
			std::cout << usingParent[i] << std::endl;
			std::cout << usingRandom[i] << std::endl;
		}*/

		float accumulatedError = annEvaluator.forward_triplet(usingObserved, usingParent, usingRandom, featureDescriptions, LearningRateSGD);
		
		errorStat.add_number(abs(accumulatedError));

		// 色々な出力
		if ((iteration % EvaluatorSerializeInterval) == 0)
		{
			std::cout << "Serializing..." << std::endl;

			std::ofstream annOut(get_filename(iteration));

			annEvaluator.serialize(annOut);
		}

		if ((iteration % IterationPrintInterval) == 0)
		{
			std::cout << "Iteration " << iteration << ". ";
			std::cout << "Time: " << (now() - timeStart) / 1000 << " seconds. ";
			std::cout << "Last Iteration took: " << (now() - iterationStart) / 1000 << " seconds. ";

			float error = errorStat.get_average();
			std::cout << "Error: " << error << ". ";
			stats << iteration << "," << error << std::endl;

			errorStat.reset();

			std::cout << std::endl;
		}
	}

	std::cout << "Congratulations! Learning is successfully finished!" << std::endl;
}

void make_sfens()
{
	const int numGames = 62510;

	std::ifstream gameDbFile(GameDatabase::defaultDatabaseFile);
	std::ofstream sfensFile(GameDatabase::outputSFENsFile);

	GameDatabase gameDb(gameDbFile);
	std::vector<Game> games, samples;

	std::cout << "start reading games...\n";

	for (int i = 0; i < numGames; ++i)
	{
		Game game;
		if (gameDb.read_one_game(&game)
			&& game.result != Game::DRAW
			&& game.moves.size() < 256)
			games.push_back(game);
	}

	std::cout << "start writing SFENs...\n";

	for (int i = 0; i < games.size(); ++i)
	{
		/*sfensFile << games[i].date << " "
			<< games[i].players[0] << " " << games[i].players[1] << " "
			<< games[i].result << " "
			<< games[i].event 
			<< games[i].opening << std::endl;*/

		StateInfo st;
		Position pos;
		pos.set(StartSFEN, &st, Threads.main());

		std::vector<StateInfo> stEach(games[i].moves.size());
		for (int ply = 0; ply < games[i].moves.size(); ++ply)
		{
			Move move = games[i].moves[ply];

			pos.do_move(move, stEach[ply], pos.gives_check(move));

			if (!pos.pos_is_ok())
				break;

			//sfensFile << (ply ? " " : "") << USI::move(move);
			sfensFile << pos.sfen() << std::endl;
		}
		//sfensFile << std::endl;
	}

	std::cout << "finish writing SFENs." << std::endl;
}

void make_triplet_sfens()
{
	const int numGames = 62510;

	std::ifstream gameDbFile(GameDatabase::defaultDatabaseFile);
	std::ofstream parentSFENsFile(GameDatabase::outputParentSFENsFile);
	std::ofstream observedSFENsFile(GameDatabase::outputObservedSFENsFile);
	std::ofstream randomSFENsFile(GameDatabase::outputRandomSFENsFile);

	GameDatabase gameDb(gameDbFile);
	std::vector<Game> games, samples;

	std::cout << "start reading games...\n";

	for (int i = 0; i < numGames; ++i)
	{
		Game game;
		if (gameDb.read_one_game(&game)
			&& game.result != Game::DRAW
			&& game.moves.size() < 256)
			games.push_back(game);
	}

	std::cout << "start writing SFENs...\n";

	for (int i = 0; i < games.size(); ++i)
	{
		StateInfo st;
		Position parentPos;
		parentPos.set(StartSFEN, &st, Threads.main());

		std::vector<StateInfo> stObserved(games[i].moves.size());
		std::vector<StateInfo> stRandom(games[i].moves.size());
		for (int ply = 0; ply < games[i].moves.size(); ++ply)
		{
			std::vector<Move> parentMoves;
			for (const auto& m : MoveList<LEGAL>(parentPos))
				parentMoves.push_back(m);

			if (parentMoves.size() == 0)
				continue;

			auto rng = gRd.make_mt();
			auto movePickerDist = std::uniform_int_distribution<size_t>(0, parentMoves.size() - 1);

			Position randomPos = parentPos;
			Move observedMove = games[i].moves[ply];
			Move randomMove = parentMoves[movePickerDist(rng)];

			parentSFENsFile << parentPos.sfen() << std::endl;

			parentPos.do_move(observedMove, stObserved[ply], parentPos.gives_check(observedMove));
			randomPos.do_move(randomMove, stRandom[ply], randomPos.gives_check(randomMove));

			observedSFENsFile << parentPos.sfen() << std::endl;
			randomSFENsFile << randomPos.sfen() << std::endl;
		}
	}

	std::cout << "finish writing SFENs." << std::endl;
}

} // namespace Learn