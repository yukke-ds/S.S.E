#include <iostream>
#include <fstream>
#include <algorithm>

#include "ann/ann_evaluator.h"
#include "position.h"
#include "evaluate.h"
#include "usi.h"

namespace Eval {

std::array<int32_t, 2> kk[SQUARE_NB][SQUARE_NB];
std::array<int32_t, 2> kkp[SQUARE_NB][SQUARE_NB][FE_END];
std::array<int16_t, 2> kpp[SQUARE_NB][FE_END][FE_END];

Value PieceValue[PIECE_NB] = {
	VALUE_ZERO, PawnValue, LanceValue, KnightValue, SilverValue, BishopValue, RookValue, GoldValue,
	KingValue, ProPawnValue, ProLanceValue, ProKnightValue, ProSilverValue, HorseValue, DragonValue, VALUE_ZERO,

	VALUE_ZERO, -PawnValue, -LanceValue, -KnightValue, -SilverValue, -BishopValue, -RookValue, -GoldValue,
	-KingValue, -ProPawnValue, -ProLanceValue, -ProKnightValue, -ProSilverValue, -HorseValue, -DragonValue, VALUE_ZERO
};

Value CapturePieceValue[PIECE_NB] = {
	VALUE_ZERO, PawnValue * 2, LanceValue * 2, KnightValue * 2, SilverValue * 2,
	BishopValue * 2, RookValue * 2, GoldValue * 2, VALUE_ZERO,
	ProPawnValue + PawnValue, ProLanceValue + LanceValue, ProKnightValue + KnightValue, ProSilverValue + SilverValue,
	HorseValue + BishopValue, DragonValue + RookValue, VALUE_ZERO /* PRO_GOLD */,

	VALUE_ZERO, PawnValue * 2, LanceValue * 2, KnightValue * 2, SilverValue * 2,
	BishopValue * 2, RookValue * 2, GoldValue * 2, VALUE_ZERO,
	ProPawnValue + PawnValue, ProLanceValue + LanceValue, ProKnightValue + KnightValue, ProSilverValue + SilverValue,
	HorseValue + BishopValue, DragonValue + RookValue, VALUE_ZERO /* PRO_GOLD */,
};

Value PromotionDiff[PIECE_NB] = {
	VALUE_ZERO, ProPawnValue - PawnValue, ProLanceValue - LanceValue, ProKnightValue - KnightValue, ProSilverValue - SilverValue, HorseValue - BishopValue, DragonValue - RookValue, VALUE_ZERO,
	VALUE_ZERO, ProPawnValue - PawnValue, ProLanceValue - LanceValue, ProKnightValue - KnightValue, ProSilverValue - SilverValue, HorseValue - BishopValue, DragonValue - RookValue, VALUE_ZERO,
	VALUE_ZERO, ProPawnValue - PawnValue, ProLanceValue - LanceValue, ProKnightValue - KnightValue, ProSilverValue - SilverValue, HorseValue - BishopValue, DragonValue - RookValue, VALUE_ZERO,
	VALUE_ZERO, ProPawnValue - PawnValue, ProLanceValue - LanceValue, ProKnightValue - KnightValue, ProSilverValue - SilverValue, HorseValue - BishopValue, DragonValue - RookValue, VALUE_ZERO,
};

void load_eval()
{
	// ANN
	if (Options["ANNEvaluator"])
	{
#ifdef DEEPPINK_REVERSE
		std::string filename = "RequiredFiles/eval.dat";
		std::ifstream netfIn(filename, ios::in | ios::binary);

		if (!netfIn)
			std::cout << "info string Failed to open " << filename << "." << std::endl;
		else
			AnnEvaluator.deserialize_py(netfIn);
#else
		std::string filename = "RequiredFiles/eval10.net";
		std::ifstream netfIn(filename);

		if (!netfIn)
			std::cout << "info string Failed to open " << filename << "." << std::endl;
		else
			AnnEvaluator.deserialize(netfIn);
#endif
	}

	// elmo
	else
	{
		// KK
		std::ifstream ifsKK("eval/KK_synthesized.bin", std::ios::binary);
		if (ifsKK) ifsKK.read(reinterpret_cast<char*>(kk), sizeof(kk));
		else goto Error;

		// KKP
		std::ifstream ifsKKP("eval/KKP_synthesized.bin", std::ios::binary);
		if (ifsKKP) ifsKKP.read(reinterpret_cast<char*>(kkp), sizeof(kkp));
		else goto Error;

		// KPP
		std::ifstream ifsKPP("eval/KPP_synthesized.bin", std::ios::binary);
		if (ifsKPP) ifsKPP.read(reinterpret_cast<char*>(kpp), sizeof(kpp));
		else goto Error;
	}

	return;

Error:;
	// 評価関数ファイルの読み込みに失敗した場合、思考を開始しないように抑制したほうがいい
	sync_cout << "\ninfo string Error! open evaluation file failed.\n" << sync_endl;
	exit(EXIT_FAILURE);
}

uint64_t calc_check_sum()
{
	uint64_t sum = 0;

	auto add_sum = [&](uint32_t* ptr, size_t t)
	{
		for (size_t i = 0; i < t; ++i)
			sum += ptr[i];
	};

	add_sum(reinterpret_cast<uint32_t*>(kk), sizeof(kk) / sizeof(uint32_t));
	add_sum(reinterpret_cast<uint32_t*>(kkp), sizeof(kkp) / sizeof(uint32_t));
	add_sum(reinterpret_cast<uint32_t*>(kpp), sizeof(kpp) / sizeof(uint32_t));

	return sum;
}

// 先手玉が移動したときに先手側の差分
std::array<int32_t, 2> bk_moved_diff(const Position& pos, const ExtEvalIndex eei)
{
	const Square sq_bk = pos.king_square(BLACK);
	const EvalIndex* list0 = pos.eval_list()->piece_list_fb();

	const auto* pkppb = kpp[sq_bk][eei.fb];
	std::array<int32_t, 2> sum = { { pkppb[list0[0]][0], pkppb[list0[0]][1] } };
	for (int i = 1; i < KING_NUMBER; ++i) {
		sum[0] += pkppb[list0[i]][0];
		sum[1] += pkppb[list0[i]][1];
	}
	return sum;
}

// 後手玉が移動したときの後手側の差分
std::array<int32_t, 2> wk_moved_diff(const Position& pos, const ExtEvalIndex eei)
{
	const Square sq_wk = pos.king_square(WHITE);
	const EvalIndex *list1 = pos.eval_list()->piece_list_fw();

	const auto* pkppw = kpp[inverse(sq_wk)][eei.fw];
	std::array<int32_t, 2> sum = { { pkppw[list1[0]][0], pkppw[list1[0]][1] } };
	for (int i = 1; i < KING_NUMBER; ++i) {
		sum[0] += pkppw[list1[i]][0];
		sum[1] += pkppw[list1[i]][1];
	}
	return sum;
}

// 玉以外の駒が移動したときの差分
Evaluate others_moved_diff(const Position& pos, const ExtEvalIndex eei)
{
	const Square sq_bk = pos.king_square(BLACK);
	const Square sq_wk = pos.king_square(WHITE);
	const EvalIndex* list0 = pos.eval_list()->piece_list_fb();
	const EvalIndex* list1 = pos.eval_list()->piece_list_fw();

	Evaluate eval;
	eval.p[0] = { 0, 0 };
	eval.p[1] = { 0, 0 };
	eval.p[2] = kkp[sq_bk][sq_wk][eei.fb];

	const auto* pkppb = kpp[sq_bk][eei.fb];
	const auto* pkppw = kpp[inverse(sq_wk)][eei.fw];

#if defined (USE_AVX2)

	__m256i zero = _mm256_setzero_si256();
	__m256i sum0 = zero;
	__m256i sum1 = zero;
	int i = 0;
	for (; i + 8 < KING_NUMBER; i += 8) {
		__m256i indexes0 = _mm256_load_si256(reinterpret_cast<const __m256i*>(&list0[i]));
		__m256i indexes1 = _mm256_load_si256(reinterpret_cast<const __m256i*>(&list1[i]));
		__m256i w0 = _mm256_i32gather_epi32(reinterpret_cast<const int*>(pkppb), indexes0, 4);
		__m256i w1 = _mm256_i32gather_epi32(reinterpret_cast<const int*>(pkppw), indexes1, 4);

		__m256i w0lo = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w0, 0));
		sum0 = _mm256_add_epi32(sum0, w0lo);
		__m256i w0hi = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w0, 1));
		sum0 = _mm256_add_epi32(sum0, w0hi);

		__m256i w1lo = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w1, 0));
		sum1 = _mm256_add_epi32(sum1, w1lo);
		__m256i w1hi = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w1, 1));
		sum1 = _mm256_add_epi32(sum1, w1hi);
	}

	for (; i + 4 < KING_NUMBER; i += 4) {
		__m128i indexes0 = _mm_load_si128(reinterpret_cast<const __m128i*>(&list0[i]));
		__m128i indexes1 = _mm_load_si128(reinterpret_cast<const __m128i*>(&list1[i]));
		__m128i w0 = _mm_i32gather_epi32(reinterpret_cast<const int*>(pkppb), indexes0, 4);
		__m128i w1 = _mm_i32gather_epi32(reinterpret_cast<const int*>(pkppw), indexes1, 4);

		__m256i w0lo = _mm256_cvtepi16_epi32(w0);
		sum0 = _mm256_add_epi32(sum0, w0lo);

		__m256i w1lo = _mm256_cvtepi16_epi32(w1);
		sum1 = _mm256_add_epi32(sum1, w1lo);
	}

	for (; i < KING_NUMBER; ++i) {
		sum.p[0] += pkppb[list0[i]];
		sum.p[1] += pkppw[list1[i]];
	}

	// sum0とsum0の上位128ビットと下位128ビットを独立して8バイトシフトしたものを足し合わせる
	sum0 = _mm256_add_epi32(sum0, _mm256_srli_si256(sum0, 8));
	// sum0の上位128ビットと下位128ビットを足しあわせてsum0_128に代入する
	__m128i sum0_128 = _mm_add_epi32(_mm256_extracti128_si256(sum0, 0), _mm256_extracti128_si256(sum0, 1));
	// sum0_128の下位64ビットをdiff.p[1]にストアする
	std::array<int32_t, 2> sum0_array;
	_mm_storel_epi64(reinterpret_cast<__m128i*>(&sum0_array), sum0_128);
	sum.p[0] += sum0_array;

	// sum1とsum1の上位128ビットと下位128ビットを独立して8バイトシフトしたものを足し合わせる
	sum1 = _mm256_add_epi32(sum1, _mm256_srli_si256(sum1, 8));
	// sum1の上位128ビットと下位128ビットを足しあわせてsum1_128に代入する
	__m128i sum1_128 = _mm_add_epi32(_mm256_extracti128_si256(sum1, 0), _mm256_extracti128_si256(sum1, 1));
	// sum1_128の下位64ビットをdiff.p[1]にストアする
	std::array<int32_t, 2> sum1_array;
	_mm_storel_epi64(reinterpret_cast<__m128i*>(&sum1_array), sum1_128);
	sum.p[1] += sum1_array;

#else
	eval.m[0] = _mm_set_epi32(0, 0, *reinterpret_cast<const int32_t*>(&pkppw[list1[0]][0]), *reinterpret_cast<const int32_t*>(&pkppb[list0[0]][0]));
	eval.m[0] = _mm_cvtepi16_epi32(eval.m[0]);
	for (int i = 1; i < KING_NUMBER; ++i) {
		__m128i tmp;
		tmp = _mm_set_epi32(0, 0, *reinterpret_cast<const int32_t*>(&pkppw[list1[i]][0]), *reinterpret_cast<const int32_t*>(&pkppb[list0[i]][0]));
		tmp = _mm_cvtepi16_epi32(tmp);
		eval.m[0] = _mm_add_epi32(eval.m[0], tmp);
	}
#endif

	return eval;
}

Value material(const Position& pos)
{
	Value v = VALUE_ZERO;

	for (Square s = SQ_11; s <= SQ_99; ++s)
		v += PieceValue[pos.piece_on(s)];

	for (Color c = BLACK; c <= WHITE; ++c)
		for (PieceType pt = PAWN; pt < PIECE_HAND_NB; ++pt)
			v += (c == BLACK ? 1 : -1) * Value(count_of(pos.hand_of(c), pt) * PieceValue[pt]);

	return v;
}

Value all_calculate(const Position& pos)
{
	Square sq_bk = pos.king_square(BLACK);
	Square sq_wk = pos.king_square(WHITE);
	const auto* ppkppb = kpp[sq_bk];
	const auto* ppkppw = kpp[inverse(sq_wk)];

	Position& pos_ = *const_cast<Position*>(&pos);

	EvalIndex* list_fb = pos_.eval_list()->piece_list_fb();
	EvalIndex* list_fw = pos_.eval_list()->piece_list_fw();

	int i, j;
	EvalIndex k0, k1, l0, l1;

	// 評価値クラス
	Evaluate eval;

	// eval.p[0](BKPP)とeval.p[1](WKPP)をゼロクリア
	eval.m[0] = _mm_setzero_si128();

	// KK
	eval.p[2] = kk[sq_bk][sq_wk];

	for (i = 0; i < KING_NUMBER; ++i)
	{
		k0 = list_fb[i];
		k1 = list_fw[i];
		const auto* pkppb = ppkppb[k0];
		const auto* pkppw = ppkppw[k1];
		for (j = 0; j < i; ++j)
		{
			l0 = list_fb[j];
			l1 = list_fw[j];

			// pkppw[l1][0],pkppw[l1][1],pkppb[l0][0],pkppb[l0][1]の16bit変数4つを整数拡張で32bit化して足し合わせる
			__m128i tmp;
			tmp = _mm_set_epi32(0, 0, *reinterpret_cast<const int32_t*>(&pkppw[l1][0]), *reinterpret_cast<const int32_t*>(&pkppb[l0][0]));
			// この命令SSE4.1の命令のはず..
			tmp = _mm_cvtepi16_epi32(tmp);
			eval.m[0] = _mm_add_epi32(eval.m[0], tmp);
		}
		eval.p[2] += kkp[sq_bk][sq_wk][k0];
	}

	StateInfo* st = pos.state();
	eval.p[2][0] += st->material * FV_SCALE;

	st->eval = eval;

	return Value(eval.sum(pos.side_to_move()) / FV_SCALE);
}

// evaluate() is the main evaluation function. It returns a static evaluation
// of the position from the point of view of the side to move.
Value evaluate(const Position& pos) {

	// DNN
	if (Options["ANNEvaluator"])
	{
		Value score = pos.this_thread()->annEvaluator->evaluate(pos);

#ifdef DEEPPINK_REVERSE
		return -score;
#else
		return pos.side_to_move() == BLACK ? score : -score;
#endif
	}

	// elmoの評価関数
	else
	{
		StateInfo* st = pos.state();
		Evaluate& eval = st->eval;

		// 現局面は計算済み
		if (eval.calculated())
			return Value(eval.sum(pos.side_to_move()) / FV_SCALE);

		else
		{
			StateInfo* prev = st->previous;

			// 1つ前の局面も未計算なら全計算
			if (!prev->eval.calculated()) 
				return all_calculate(pos);

			// 1つ前は計算済みであるから差分計算
			else
			{
				// 遡るnodeは一つだけ
				// ひとつずつ遡りながらsumKPPがVALUE_NONEでないところまで探してそこからの差分を計算することは出来るが
				// 現状、探索部では毎node、evaluate()を呼び出すから問題ない。

				MovedPieceInfo& mpi = st->mpi;

				// 移動させた駒は最大2つある。その数
				int numPieceMoved = mpi.numPieceMoved;

				EvalIndex* list0 = pos.eval_list()->piece_list_fb();
				EvalIndex* list1 = pos.eval_list()->piece_list_fw();

				PieceNumber movedPieceNumber = mpi.movedPieceNumber[0];

				// 移動させた駒は王か？
				if (movedPieceNumber >= KING_NUMBER)
				{
					// 前のnodeの評価値からの増分を計算していく。
					// (直接この変数に加算していく)
					// この意味においてdiffという名前は少々不適切ではあるが。
					Evaluate diff = prev->eval;

					Square sq_bk = pos.king_square(BLACK);
					Square sq_wk = pos.king_square(WHITE);

					// ΣKKPは最初から全計算するしかないので初期化する。
					diff.p[2] = kk[sq_bk][sq_wk];
					diff.p[2][0] += st->material * FV_SCALE;

					// 後手玉の移動(片側分のKPPを丸ごと求める)
					if (movedPieceNumber == WKING_NUMBER)
					{
						const auto ppkppw = kpp[inverse(sq_wk)];

						// ΣWKPP = 0
						diff.p[1][0] = 0;
						diff.p[1][1] = 0;

#if defined(USE_AVX2)

						__m256i zero = _mm256_setzero_si256();
						__m256i diffp1 = zero;
						for (int i = 0; i < KING_NUMBER; ++i)
						{
							const int k1 = list1[i];
							const auto* pkppw = ppkppw[k1];
							int j = 0;
							for (; j + 8 < i; j += 8)
							{
								// list1[j]から8要素ロードする
								__m256i indexes = _mm256_load_si256(reinterpret_cast<const __m256i*>(&list1[j]));
								// indexesのオフセットに従い、pkppwから8要素ギャザーする
								__m256i w = _mm256_i32gather_epi32(reinterpret_cast<const int*>(pkppw), indexes, 4);
								// 下位128ビットを16ビット整数→32ビット整数に変換する
								__m256i wlo = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w, 0));
								// diffp1に足し合わせる
								diffp1 = _mm256_add_epi32(diffp1, wlo);
								// 上位128ビットを16ビット整数→32ビット整数に変換する
								__m256i whi = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w, 1));
								// diffp1に足し合わせる
								diffp1 = _mm256_add_epi32(diffp1, whi);
							}

							for (; j + 4 < i; j += 4) {
								// list1[j]から4要素ロードする
								__m128i indexes = _mm_load_si128(reinterpret_cast<const __m128i*>(&list1[j]));
								// indexesのオフセットに従い、pkppwから4要素ギャザーする
								__m128i w = _mm_i32gather_epi32(reinterpret_cast<const int*>(pkppw), indexes, 4);
								// 16ビット整数→32ビット整数に変換する
								__m256i wlo = _mm256_cvtepi16_epi32(w);
								// diffp1に足し合わせる
								diffp1 = _mm256_add_epi32(diffp1, wlo);
							}

							for (; j < i; ++j)
							{
								const int l1 = list1[j];
								diff.p[1] += pkppw[l1];
							}

							// KKPのWK分。BKは移動していないから、BK側には影響ない。

							// 後手から見たKKP。後手から見ているのでマイナス
							diff.p[2][0] -= kkp[inverse(sq_wk)][inverse(sq_bk)][k1][0];
							// 後手から見たKKP手番。後手から見るのでマイナスだが、手番は先手から見たスコアを格納するのでさらにマイナスになって、プラス。
							diff.p[2][1] += kkp[inverse(sq_wk)][inverse(sq_bk)][k1][1];
						}

						// diffp1とdiffp1の上位128ビットと下位128ビットを独立して8バイトシフトしたものを足し合わせる
						diffp1 = _mm256_add_epi32(diffp1, _mm256_srli_si256(diffp1, 8));
						// diffp1の上位128ビットと下位128ビットを足しあわせてdiffp1_128に代入する
						__m128i diffp1_128 = _mm_add_epi32(_mm256_extracti128_si256(diffp1, 0), _mm256_extracti128_si256(diffp1, 1));
						// diffp1_128の下位64ビットをdiff.p[1]にストアする
						std::array<int32_t, 2> diffp1_sum;
						_mm_storel_epi64(reinterpret_cast<__m128i*>(&diffp1_sum), diffp1_128);
						diff.p[1] += diffp1_sum;
#else

						for (int i = 0; i < KING_NUMBER; ++i)
						{
							const int k1 = list1[i];
							const auto* pkppw = ppkppw[k1];
							for (int j = 0; j < i; ++j)
							{
								const int l1 = list1[j];
								diff.p[1] += pkppw[l1];
							}

							// KKPのWK分。BKは移動していないから、BK側には影響ない。

							// 後手から見たKKP。後手から見ているのでマイナス
							diff.p[2][0] -= kkp[inverse(sq_wk)][inverse(sq_bk)][k1][0];
							// 後手から見たKKP手番。後手から見るのでマイナスだが、手番は先手から見たスコアを格納するのでさらにマイナスになって、プラス。
							diff.p[2][1] += kkp[inverse(sq_wk)][inverse(sq_bk)][k1][1];
						}
#endif

						// 動かした駒が２つ
						if (numPieceMoved == 2)
						{
							// 瞬間的にeval_listの移動させた駒の番号を変更してしまう。
							// こうすることで前nodeのpiece_listを持たなくて済む。

							const int listIndex_cap = mpi.movedPieceNumber[1];
							diff.p[0] += bk_moved_diff(pos, mpi.changedPiece[1].newPiece);
							list0[listIndex_cap] = mpi.changedPiece[1].oldPiece.fb;
							diff.p[0] -= bk_moved_diff(pos, mpi.changedPiece[1].oldPiece);
							list0[listIndex_cap] = mpi.changedPiece[1].newPiece.fb;
						}
					}
					else {

						// 先手玉の移動
						// さきほどの処理と同様。

						const auto* ppkppb = kpp[sq_bk];
						diff.p[0][0] = 0;
						diff.p[0][1] = 0;

#if defined(USE_AVX2)

						__m256i zero = _mm256_setzero_si256();
						__m256i diffp0 = zero;
						for (int i = 0; i < KING_NUMBER; ++i)
						{
							const int k0 = list0[i];
							const auto* pkppb = ppkppb[k0];
							int j = 0;
							for (; j + 8 < i; j += 8)
							{
								// list0[j]から8要素ロードする
								__m256i indexes = _mm256_load_si256(reinterpret_cast<const __m256i*>(&list0[j]));
								// indexesのオフセットに従い、pkppwから8要素ギャザーする
								__m256i w = _mm256_i32gather_epi32(reinterpret_cast<const int*>(pkppb), indexes, 4);
								// 下位128ビットを16ビット整数→32ビット整数に変換する
								__m256i wlo = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w, 0));
								// diffp0に足し合わせる
								diffp0 = _mm256_add_epi32(diffp0, wlo);
								// 上位128ビットを16ビット整数→32ビット整数に変換する
								__m256i whi = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w, 1));
								// diffp0に足し合わせる
								diffp0 = _mm256_add_epi32(diffp0, whi);
							}

							for (; j + 4 < i; j += 4) {
								// list0[j]から4要素ロードする
								__m128i indexes = _mm_load_si128(reinterpret_cast<const __m128i*>(&list0[j]));
								// indexesのオフセットに従い、pkppwから4要素ギャザーする
								__m128i w = _mm_i32gather_epi32(reinterpret_cast<const int*>(pkppb), indexes, 4);
								// 16ビット整数→32ビット整数に変換する
								__m256i wlo = _mm256_cvtepi16_epi32(w);
								// diffp0に足し合わせる
								diffp0 = _mm256_add_epi32(diffp0, wlo);
							}

							for (; j < i; ++j)
							{
								const int l0 = list0[j];
								diff.p[0] += pkppb[l0];
							}

							diff.p[2] += kkp[sq_bk][sq_wk][k0];
						}

						// diffp0とdiffp0の上位128ビットと下位128ビットを独立して8バイトシフトしたものを足し合わせる
						diffp0 = _mm256_add_epi32(diffp0, _mm256_srli_si256(diffp0, 8));
						// diffp0の上位128ビットと下位128ビットを足しあわせてdiffp0_128に代入する
						__m128i diffp0_128 = _mm_add_epi32(_mm256_extracti128_si256(diffp0, 0), _mm256_extracti128_si256(diffp0, 1));
						// diffp0_128の下位64ビットをdiff.p[1]にストアする
						std::array<int32_t, 2> diffp0_sum;
						_mm_storel_epi64(reinterpret_cast<__m128i*>(&diffp0_sum), diffp0_128);
						diff.p[0] += diffp0_sum;
#else
						for (int i = 0; i < KING_NUMBER; ++i)
						{
							const int k0 = list0[i];
							const auto* pkppb = ppkppb[k0];
							for (int j = 0; j < i; ++j) {
								const int l0 = list0[j];
								diff.p[0] += pkppb[l0];
							}
							diff.p[2] += kkp[sq_bk][sq_wk][k0];
						}
#endif

						if (numPieceMoved == 2) {
							const int listIndex_cap = mpi.movedPieceNumber[1];
							diff.p[1] += wk_moved_diff(pos, mpi.changedPiece[1].newPiece);
							list1[listIndex_cap] = mpi.changedPiece[1].oldPiece.fw;
							diff.p[1] -= wk_moved_diff(pos, mpi.changedPiece[1].oldPiece);
							list1[listIndex_cap] = mpi.changedPiece[1].newPiece.fw;
						}
					}

					// sumの計算が終わったのでpos.state()->sumに反映させておく。(これがこの関数の返し値に相当する。)
					st->eval = diff;

				}
				else {

					// 王以外の駒が移動したケース
					// 今回の差分を計算して、そこに加算する。

					const int listIndex = mpi.movedPieceNumber[0];

					Evaluate diff = others_moved_diff(pos, mpi.changedPiece[0].newPiece);
					if (numPieceMoved == 1) {

						// 動いた駒が1つ。
						list0[listIndex] = mpi.changedPiece[0].oldPiece.fb;
						list1[listIndex] = mpi.changedPiece[0].oldPiece.fw;
						diff -= others_moved_diff(pos, mpi.changedPiece[0].oldPiece);

					}
					else {

						// 動いた駒が2つ。

						auto sq_bk = pos.king_square(BLACK);
						auto sq_wk = pos.king_square(WHITE);

						diff += others_moved_diff(pos, mpi.changedPiece[1].newPiece);
						diff.p[0] -= kpp[sq_bk][mpi.changedPiece[0].newPiece.fb][mpi.changedPiece[1].newPiece.fb];
						diff.p[1] -= kpp[inverse(sq_wk)][mpi.changedPiece[0].newPiece.fw][mpi.changedPiece[1].newPiece.fw];

						const PieceNumber listIndex_cap = mpi.movedPieceNumber[1];
						list0[listIndex_cap] = mpi.changedPiece[1].oldPiece.fb;
						list1[listIndex_cap] = mpi.changedPiece[1].oldPiece.fw;

						list0[listIndex] = mpi.changedPiece[0].oldPiece.fb;
						list1[listIndex] = mpi.changedPiece[0].oldPiece.fw;
						diff -= others_moved_diff(pos, mpi.changedPiece[0].oldPiece);
						diff -= others_moved_diff(pos, mpi.changedPiece[1].oldPiece);

						diff.p[0] += kpp[sq_bk][mpi.changedPiece[0].oldPiece.fb][mpi.changedPiece[1].oldPiece.fb];
						diff.p[1] += kpp[inverse(sq_wk)][mpi.changedPiece[0].oldPiece.fw][mpi.changedPiece[1].oldPiece.fw];
						list0[listIndex_cap] = mpi.changedPiece[1].newPiece.fb;
						list1[listIndex_cap] = mpi.changedPiece[1].newPiece.fw;
					}

					list0[listIndex] = mpi.changedPiece[0].newPiece.fb;
					list1[listIndex] = mpi.changedPiece[0].newPiece.fw;

					// 前nodeからの駒割りの増分を加算。
					diff.p[2][0] += (st->material - prev->material) * FV_SCALE;

					st->eval = diff + prev->eval;
				}
			}

			return Value(st->eval.sum(pos.side_to_move()) / FV_SCALE);
		}
	}
}

void evaluate_with_no_return(const Position& pos)
{
	// まだ評価値が計算されていないなら
	if (!pos.state()->eval.calculated())
		evaluate(pos);
}

} // namespace Eval

// この局面の手番は c側にあるものとする。c側から見た評価値を返す。
int32_t Evaluate::sum(const Color c) const 
{
	// NDF(2014)の手番評価の手法。
	// cf. http://www.computer-shogi.org/wcsc24/appeal/NineDayFever/NDF.txt

	// 手番に依存しない評価値合計
	// p[1][0]はΣWKPPなので符号はマイナス。
	const int32_t scoreBoard = p[0][0] - p[1][0] + p[2][0];
	// 手番に依存する評価値合計
	const int32_t scoreTurn = p[0][1] + p[1][1] + p[2][1];

	// この関数は手番側から見た評価値を返すのでscoreTurnは必ずプラス

	return (c == BLACK ? scoreBoard : -scoreBoard) + scoreTurn;
}

void Evaluate::init() {
}
