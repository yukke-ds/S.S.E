#define _CRT_SECURE_NO_WARNINGS
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

#include "common/math.h"
#include "gamedb.h"
#include "progress.h"
#include "psq.h"
#include "thread.h"
#include "usi.h"

int32_t Progress::weights[SQUARE_NB][PSQ_MAX];
int Progress::Relation[SQUARE_NB][SQUARE_NB];

void Progress::read_weights() 
{
	FILE* fp = std::fopen("eval/progress.bin", "rb");
	if (fp != NULL)
		std::fread(&weights, sizeof(weights), 1, fp);
	else
		std::printf("info string Failed to open progress.bin.\n");
}


double Progress::estimate_progress(const Position& pos, const PsqList& psqList) 
{
	int64_t sum = 0;
	Square sq_bk = pos.king_square(BLACK);
	Square sq_wk = inverse(pos.king_square(WHITE));
	
	for (const PsqPair& psq : psqList) {
		sum += weights[sq_bk][psq.black()];
		sum += weights[sq_wk][psq.white()];
	}

	return math::sigmoid(double(sum) * double(1.0 / weightScale));
}

double Progress::estimate_progress(const Position& pos) 
{
	const PsqList psqList(pos);
	return estimate_progress(pos, psqList);
}

namespace {

	template <typename T>
	struct WeightsBase {
		void Clear() { std::memset(this, 0, sizeof(*this)); }
		size_t size() const { return sizeof(*this) / sizeof(T); }
		T* begin() { return reinterpret_cast<T*>(this); }
		T* end() { return begin() + size(); }
		const T* begin() const { return reinterpret_cast<const T*>(this); }
		const T* end()   const { return begin() + size(); }
		T operator[](size_t n) const { return *(begin() + n); }
		T& operator[](size_t n) { return *(begin() + n); }
		T psq[PSQ_MAX];
		T relative_kp[PIECE_NB][154];
		T absolute_kp[SQUARE_NB][PSQ_MAX];
	};

	typedef WeightsBase<double> Weights;

	void update_weights(const Position& pos, double inc, Weights* const w) {
		assert(w);
		Square sq_bk = pos.king_square(BLACK);
		Square sq_wk = inverse(pos.king_square(WHITE));
		for (const PsqPair& psq : PsqList(pos)) {
			const PsqIndex i_b = psq.black();
			const PsqIndex i_w = psq.white();
			// 1. 駒の絶対位置
			w->psq[i_b] += inc;
			w->psq[i_w] += inc;
			if (i_b.square() != SQ_NONE) {
				// 2. 玉との相対位置
				w->relative_kp[i_b.piece()][Progress::relation(sq_bk, i_b.square())] += inc;
				w->relative_kp[i_w.piece()][Progress::relation(sq_wk, i_w.square())] += inc;
			}
			// 3. 玉との絶対位置
			w->absolute_kp[sq_bk][i_b] += inc;
			w->absolute_kp[sq_wk][i_w] += inc;
		}
	}

} // namespace

void Progress::learn_parameters()
{
	constexpr int numIterations = 2048;		// パラメタ更新の反復回数
	constexpr int numGames = 62627;			// 学習に利用する対局数
	constexpr int numSamples = 1000;		// 交差検定で用いるサンプル数
	constexpr double L1Penalty = 0.0;		// L1正則化係数
	constexpr double L2Penalty = 256.0;		// L2正則化係数
	constexpr double initialStep = 1e-3;	// 重みベクトルの最初の更新幅
	constexpr double decay = 0.95;			// AdaDelta原論文のρ
	constexpr double momentum = 0.95;		// モーメンタムの強さ
	constexpr double epsilon = 1e-6;		// AdaDelta原論文のε

	const int numThreads = std::max(1U, std::thread::hardware_concurrency());
	omp_set_num_threads(numThreads);
	std::cout << "Set numThreads = " << numThreads << "\n";

	Position startpos;
	StateInfo st;
	startpos.set(StartSFEN, &st, Threads.main());

#ifdef KIF
	std::ifstream gameDbFile(GameDatabase::defaultDatabaseFile);
	GameDatabase gameDb(gameDbFile);
#else
	std::ifstream trainingFile(GameDatabase::defaultTrainingFile);
	std::ifstream testFile(GameDatabase::defaultTestFile);
#endif
	std::vector<Game> games, samples;

	std::cout << "start reading games.\n";

#ifdef KIF
	for (int i = 0; i < numGames; ++i)
	{
		Game game;
		if (gameDb.read_one_game(&game)
			&& game.result != Game::DRAW
			&& game.moves.size() < 256)
			games.push_back(game);
	}
#else
	std::string kifu;
	while (std::getline(trainingFile, kifu))
	{
		std::istringstream is(kifu);
		std::string token;
		Game game;

		while (is >> std::skipws >> token)
		{
			if (token == "startpos" || token == "moves")
				continue;

			Move m = USI::to_move(token);
			game.moves.push_back(m);
		}

		games.push_back(game);
	}
#endif

#ifdef KIF
	for (int i = 0; i < numSamples; ++i)
	{
		Game game;
		if (gameDb.read_one_game(&game)
			&& game.result != Game::DRAW
			&& game.moves.size() < 256)
			samples.push_back(game);
	}
#else
	while (std::getline(testFile, kifu))
	{
		std::istringstream is(kifu);
		std::string token;
		Game game;

		while (is >> std::skipws >> token)
		{
			if (token == "startpos" || token == "moves")
				continue;

			Move m = USI::to_move(token);
			game.moves.push_back(m);
		}

		samples.push_back(game);
	}
#endif

	std::cout << "finish reading games.\n";
	std::cout << "games: " << games.size() << " samples: " << samples.size() << "\n";

	Weights currentWeights, momentums, accumulatedGradients, accumulatedDeltas;
	std::vector<Weights> thread_local_gradients(numThreads);
	momentums.Clear();
	currentWeights.Clear();
	accumulatedGradients.Clear();
	std::fill(accumulatedDeltas.begin(), accumulatedDeltas.end(),
		(1.0 - decay) * std::pow(initialStep, 2.0) / std::pow(1.0 - momentum, 2.0));

	for (int iteration = 1; iteration <= numIterations; ++iteration)
	{
		std::cout << "iteration=" << iteration;

		int numMoves = 0;
		double offset, sumSquareDiff;

		offset = sumSquareDiff = 0;
		for (Weights& g : thread_local_gradients)
			g.Clear();

		// 1. calculate the gradients
#pragma omp parallel for reduction(+:offset, sumSquareDiff, numMoves) schedule(dynamic)
		for (int64_t game_id = 0; game_id < games.size(); ++game_id)
		{
			const Game& game = games[game_id];
			Position pos;
			StateInfo st;

			pos.set(StartSFEN, &st, Threads.main());

			for (int ply = 0, n = game.moves.size(); ply < n; ++ply)
			{
				double teacher = double(ply) / double(n);
				double actual = estimate_progress(pos);
				double diff = actual - teacher;

				offset += diff;
				sumSquareDiff += (diff * diff);
				numMoves += 1;

				int threadId = omp_get_thread_num();
				double delta = diff * math::derivative_of_sigmoid(actual);
				update_weights(pos, delta, &thread_local_gradients[threadId]);

				pos.do_move(game.moves[ply], st, pos.gives_check(game.moves[ply]));
			}
		}

		// 2. aggregate the gradients
		Weights gradients;
		gradients.Clear();
		for (const Weights& g : thread_local_gradients)
			for (size_t i = 0; i < g.size(); ++i)
				gradients[i] += g[i];

		// 3. update the weights
		double lasso = 0, tikhonov = 0;
#pragma omp parallel for reduction(+:lasso, tikhonov) schedule(static)
		for (int64_t i = 0; i < gradients.size(); ++i)
		{
			double& w = currentWeights[i];
			double& r = accumulatedGradients[i];
			double& s = accumulatedDeltas[i];
			double& v = momentums[i];
			double gradient = gradients[i];
			// 損失を計算する
			lasso += L1Penalty * std::abs(w);
			tikhonov += L2Penalty * (w * w);
			// 勾配方向に移動する（AdaDelta + Momentum）
			r = decay * r + (1.0 - decay) * (gradient * gradient);
			double eta = std::sqrt(s) / std::sqrt(r + epsilon);
			v = momentum * v - (1.0 - momentum) * eta * gradient;
			s = decay * s + (1.0 - decay) * (v * v);
			w = w + v;
			// ペナルティをかける（FOBOS）
			// （参考文献）
			//   - John Duchi, Yoram Singer: Efficient Learning with Forward-Backward Splitting,
			//     http://web.stanford.edu/~jduchi/projects/DuchiSi09c_slides.pdf,
			//     pp.12-13, 2009.
			//   - Zachary C. Lipton, Charles Elkan: Efficient Elastic Net Regularization for Sparse Linear Models,
			//     http://zacklipton.com/media/papers/lazy-updates-elastic-net-lipton2015_1.pdf,
			//     p.9, 2015.
			double lambda1 = eta * L1Penalty;
			double lambda2 = eta * L2Penalty;
			w = math::sign(w) * std::max((std::abs(w) - lambda1) / (1.0 + lambda2), 0.0);
		}

		// 4. copy to the weight for the calculation
		for (Square ksq = SQ_11; ksq <= SQ_99; ++ksq) {
			for (PsqIndex i = PsqIndex(PSQ_MIN); i < PSQ_MAX; i = PsqIndex(i + 1)) {
				double w = currentWeights.psq[i];
				if (i.square() != SQ_NONE)
					w += currentWeights.relative_kp[i.piece()][relation(ksq, i.square())];

				w += currentWeights.absolute_kp[ksq][i];
				weights[ksq][i] = int32_t(w * weightScale);
			}
		}

		// 5. carry out cross-validation
		int numSamples = 0;
		double sampleSquareDiff = 0;
#pragma omp parallel for reduction(+:numSamples, sampleSquareDiff) schedule(dynamic)
		for (int64_t sampleId = 0; sampleId < samples.size(); ++sampleId)
		{
			const Game& sample = samples[sampleId];
			Position pos;
			pos.set(StartSFEN, &st, Threads.main());

			for (int ply = 0, n = sample.moves.size(); ply < n; ++ply) {
				double teacher = double(ply) / double(n);
				double diff = teacher - estimate_progress(pos);
				sampleSquareDiff += diff * diff;
				++numSamples;

				pos.do_move(sample.moves[ply], st, pos.gives_check(sample.moves[ply]));
			}
		}

		// 6. information disclosure
		double loss = sumSquareDiff + lasso;
		double stddev = std::sqrt(sumSquareDiff / numMoves);
		double sample_stddev = std::sqrt(sampleSquareDiff / numSamples);
		std::cout << std::fixed << std::setprecision(4)
			<< " Loss=" << loss << " L1=" << lasso << " L2=" << tikhonov
			<< " Stddev=" << stddev << " Validation=" << sample_stddev
			<< " Startpos=" << estimate_progress(startpos)
			<< " moves=" << numMoves << "\n";
	}

	// export to 'progress.bin' file
	FILE* fout = std::fopen("progress.bin", "wb");
	std::fwrite(&weights, sizeof weights, 1, fout);
	std::fclose(fout);

	// print the results
	for (Color c = BLACK; c <= WHITE; ++c)
		for (PieceType pt = PAWN; pt <= DRAGON; ++pt)
		{
			Piece pc = make_piece(c, pt);
			std::cout << "Piece:" << pc << "\n";
			for (Rank r = RANK_1; r <= RANK_9; ++r)
				for (File f = FILE_9; f >= FILE_1; --f) {
					Square ksq = SQ_88;
					PsqIndex i = PsqIndex::of_board(pc, f | r);
					std::cout << " " << std::fixed << std::setprecision(6)
						<< double(weights[ksq][i]) / weightScale;
				}
			std::cout << "\n";
		}
}

void Progress::init() {

	// 進行度パラメータをファイル（"progress.bin"）から読み込む
	read_weights();

	// relation
	for (Square from = SQ_11; from <= SQ_99; ++from) {
		for (Square to = SQ_11; to <= SQ_99; ++to) {
			int x = std::abs(file_of(to) - file_of(from));
			int y = rank_of(to) - rank_of(from) + RANK_9;
			int r = x + 9 * y;
			assert(0 <= r && r <= 152);
			Relation[from][to] = r;
		}
	}
}