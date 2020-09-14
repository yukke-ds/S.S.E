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
	// �]���֐��t�@�C���̓ǂݍ��݂Ɏ��s�����ꍇ�A�v�l���J�n���Ȃ��悤�ɗ}�������ق�������
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

// ���ʂ��ړ������Ƃ��ɐ�葤�̍���
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

// ���ʂ��ړ������Ƃ��̌�葤�̍���
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

// �ʈȊO�̋�ړ������Ƃ��̍���
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

	// sum0��sum0�̏��128�r�b�g�Ɖ���128�r�b�g��Ɨ�����8�o�C�g�V�t�g�������̂𑫂����킹��
	sum0 = _mm256_add_epi32(sum0, _mm256_srli_si256(sum0, 8));
	// sum0�̏��128�r�b�g�Ɖ���128�r�b�g�𑫂����킹��sum0_128�ɑ������
	__m128i sum0_128 = _mm_add_epi32(_mm256_extracti128_si256(sum0, 0), _mm256_extracti128_si256(sum0, 1));
	// sum0_128�̉���64�r�b�g��diff.p[1]�ɃX�g�A����
	std::array<int32_t, 2> sum0_array;
	_mm_storel_epi64(reinterpret_cast<__m128i*>(&sum0_array), sum0_128);
	sum.p[0] += sum0_array;

	// sum1��sum1�̏��128�r�b�g�Ɖ���128�r�b�g��Ɨ�����8�o�C�g�V�t�g�������̂𑫂����킹��
	sum1 = _mm256_add_epi32(sum1, _mm256_srli_si256(sum1, 8));
	// sum1�̏��128�r�b�g�Ɖ���128�r�b�g�𑫂����킹��sum1_128�ɑ������
	__m128i sum1_128 = _mm_add_epi32(_mm256_extracti128_si256(sum1, 0), _mm256_extracti128_si256(sum1, 1));
	// sum1_128�̉���64�r�b�g��diff.p[1]�ɃX�g�A����
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

	// �]���l�N���X
	Evaluate eval;

	// eval.p[0](BKPP)��eval.p[1](WKPP)���[���N���A
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

			// pkppw[l1][0],pkppw[l1][1],pkppb[l0][0],pkppb[l0][1]��16bit�ϐ�4�𐮐��g����32bit�����đ������킹��
			__m128i tmp;
			tmp = _mm_set_epi32(0, 0, *reinterpret_cast<const int32_t*>(&pkppw[l1][0]), *reinterpret_cast<const int32_t*>(&pkppb[l0][0]));
			// ���̖���SSE4.1�̖��߂̂͂�..
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

	// elmo�̕]���֐�
	else
	{
		StateInfo* st = pos.state();
		Evaluate& eval = st->eval;

		// ���ǖʂ͌v�Z�ς�
		if (eval.calculated())
			return Value(eval.sum(pos.side_to_move()) / FV_SCALE);

		else
		{
			StateInfo* prev = st->previous;

			// 1�O�̋ǖʂ����v�Z�Ȃ�S�v�Z
			if (!prev->eval.calculated()) 
				return all_calculate(pos);

			// 1�O�͌v�Z�ς݂ł��邩�獷���v�Z
			else
			{
				// �k��node�͈����
				// �ЂƂ��k��Ȃ���sumKPP��VALUE_NONE�łȂ��Ƃ���܂ŒT���Ă�������̍������v�Z���邱�Ƃ͏o���邪
				// ����A�T�����ł͖�node�Aevaluate()���Ăяo��������Ȃ��B

				MovedPieceInfo& mpi = st->mpi;

				// �ړ���������͍ő�2����B���̐�
				int numPieceMoved = mpi.numPieceMoved;

				EvalIndex* list0 = pos.eval_list()->piece_list_fb();
				EvalIndex* list1 = pos.eval_list()->piece_list_fw();

				PieceNumber movedPieceNumber = mpi.movedPieceNumber[0];

				// �ړ���������͉����H
				if (movedPieceNumber >= KING_NUMBER)
				{
					// �O��node�̕]���l����̑������v�Z���Ă����B
					// (���ڂ��̕ϐ��ɉ��Z���Ă���)
					// ���̈Ӗ��ɂ�����diff�Ƃ������O�͏��X�s�K�؂ł͂��邪�B
					Evaluate diff = prev->eval;

					Square sq_bk = pos.king_square(BLACK);
					Square sq_wk = pos.king_square(WHITE);

					// ��KKP�͍ŏ�����S�v�Z���邵���Ȃ��̂ŏ���������B
					diff.p[2] = kk[sq_bk][sq_wk];
					diff.p[2][0] += st->material * FV_SCALE;

					// ���ʂ̈ړ�(�Б�����KPP���ۂ��Ƌ��߂�)
					if (movedPieceNumber == WKING_NUMBER)
					{
						const auto ppkppw = kpp[inverse(sq_wk)];

						// ��WKPP = 0
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
								// list1[j]����8�v�f���[�h����
								__m256i indexes = _mm256_load_si256(reinterpret_cast<const __m256i*>(&list1[j]));
								// indexes�̃I�t�Z�b�g�ɏ]���Apkppw����8�v�f�M���U�[����
								__m256i w = _mm256_i32gather_epi32(reinterpret_cast<const int*>(pkppw), indexes, 4);
								// ����128�r�b�g��16�r�b�g������32�r�b�g�����ɕϊ�����
								__m256i wlo = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w, 0));
								// diffp1�ɑ������킹��
								diffp1 = _mm256_add_epi32(diffp1, wlo);
								// ���128�r�b�g��16�r�b�g������32�r�b�g�����ɕϊ�����
								__m256i whi = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w, 1));
								// diffp1�ɑ������킹��
								diffp1 = _mm256_add_epi32(diffp1, whi);
							}

							for (; j + 4 < i; j += 4) {
								// list1[j]����4�v�f���[�h����
								__m128i indexes = _mm_load_si128(reinterpret_cast<const __m128i*>(&list1[j]));
								// indexes�̃I�t�Z�b�g�ɏ]���Apkppw����4�v�f�M���U�[����
								__m128i w = _mm_i32gather_epi32(reinterpret_cast<const int*>(pkppw), indexes, 4);
								// 16�r�b�g������32�r�b�g�����ɕϊ�����
								__m256i wlo = _mm256_cvtepi16_epi32(w);
								// diffp1�ɑ������킹��
								diffp1 = _mm256_add_epi32(diffp1, wlo);
							}

							for (; j < i; ++j)
							{
								const int l1 = list1[j];
								diff.p[1] += pkppw[l1];
							}

							// KKP��WK���BBK�͈ړ����Ă��Ȃ�����ABK���ɂ͉e���Ȃ��B

							// ��肩�猩��KKP�B��肩�猩�Ă���̂Ń}�C�i�X
							diff.p[2][0] -= kkp[inverse(sq_wk)][inverse(sq_bk)][k1][0];
							// ��肩�猩��KKP��ԁB��肩�猩��̂Ń}�C�i�X�����A��Ԃ͐�肩�猩���X�R�A���i�[����̂ł���Ƀ}�C�i�X�ɂȂ��āA�v���X�B
							diff.p[2][1] += kkp[inverse(sq_wk)][inverse(sq_bk)][k1][1];
						}

						// diffp1��diffp1�̏��128�r�b�g�Ɖ���128�r�b�g��Ɨ�����8�o�C�g�V�t�g�������̂𑫂����킹��
						diffp1 = _mm256_add_epi32(diffp1, _mm256_srli_si256(diffp1, 8));
						// diffp1�̏��128�r�b�g�Ɖ���128�r�b�g�𑫂����킹��diffp1_128�ɑ������
						__m128i diffp1_128 = _mm_add_epi32(_mm256_extracti128_si256(diffp1, 0), _mm256_extracti128_si256(diffp1, 1));
						// diffp1_128�̉���64�r�b�g��diff.p[1]�ɃX�g�A����
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

							// KKP��WK���BBK�͈ړ����Ă��Ȃ�����ABK���ɂ͉e���Ȃ��B

							// ��肩�猩��KKP�B��肩�猩�Ă���̂Ń}�C�i�X
							diff.p[2][0] -= kkp[inverse(sq_wk)][inverse(sq_bk)][k1][0];
							// ��肩�猩��KKP��ԁB��肩�猩��̂Ń}�C�i�X�����A��Ԃ͐�肩�猩���X�R�A���i�[����̂ł���Ƀ}�C�i�X�ɂȂ��āA�v���X�B
							diff.p[2][1] += kkp[inverse(sq_wk)][inverse(sq_bk)][k1][1];
						}
#endif

						// ����������Q��
						if (numPieceMoved == 2)
						{
							// �u�ԓI��eval_list�̈ړ���������̔ԍ���ύX���Ă��܂��B
							// �������邱�ƂőOnode��piece_list�������Ȃ��čςށB

							const int listIndex_cap = mpi.movedPieceNumber[1];
							diff.p[0] += bk_moved_diff(pos, mpi.changedPiece[1].newPiece);
							list0[listIndex_cap] = mpi.changedPiece[1].oldPiece.fb;
							diff.p[0] -= bk_moved_diff(pos, mpi.changedPiece[1].oldPiece);
							list0[listIndex_cap] = mpi.changedPiece[1].newPiece.fb;
						}
					}
					else {

						// ���ʂ̈ړ�
						// �����قǂ̏����Ɠ��l�B

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
								// list0[j]����8�v�f���[�h����
								__m256i indexes = _mm256_load_si256(reinterpret_cast<const __m256i*>(&list0[j]));
								// indexes�̃I�t�Z�b�g�ɏ]���Apkppw����8�v�f�M���U�[����
								__m256i w = _mm256_i32gather_epi32(reinterpret_cast<const int*>(pkppb), indexes, 4);
								// ����128�r�b�g��16�r�b�g������32�r�b�g�����ɕϊ�����
								__m256i wlo = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w, 0));
								// diffp0�ɑ������킹��
								diffp0 = _mm256_add_epi32(diffp0, wlo);
								// ���128�r�b�g��16�r�b�g������32�r�b�g�����ɕϊ�����
								__m256i whi = _mm256_cvtepi16_epi32(_mm256_extracti128_si256(w, 1));
								// diffp0�ɑ������킹��
								diffp0 = _mm256_add_epi32(diffp0, whi);
							}

							for (; j + 4 < i; j += 4) {
								// list0[j]����4�v�f���[�h����
								__m128i indexes = _mm_load_si128(reinterpret_cast<const __m128i*>(&list0[j]));
								// indexes�̃I�t�Z�b�g�ɏ]���Apkppw����4�v�f�M���U�[����
								__m128i w = _mm_i32gather_epi32(reinterpret_cast<const int*>(pkppb), indexes, 4);
								// 16�r�b�g������32�r�b�g�����ɕϊ�����
								__m256i wlo = _mm256_cvtepi16_epi32(w);
								// diffp0�ɑ������킹��
								diffp0 = _mm256_add_epi32(diffp0, wlo);
							}

							for (; j < i; ++j)
							{
								const int l0 = list0[j];
								diff.p[0] += pkppb[l0];
							}

							diff.p[2] += kkp[sq_bk][sq_wk][k0];
						}

						// diffp0��diffp0�̏��128�r�b�g�Ɖ���128�r�b�g��Ɨ�����8�o�C�g�V�t�g�������̂𑫂����킹��
						diffp0 = _mm256_add_epi32(diffp0, _mm256_srli_si256(diffp0, 8));
						// diffp0�̏��128�r�b�g�Ɖ���128�r�b�g�𑫂����킹��diffp0_128�ɑ������
						__m128i diffp0_128 = _mm_add_epi32(_mm256_extracti128_si256(diffp0, 0), _mm256_extracti128_si256(diffp0, 1));
						// diffp0_128�̉���64�r�b�g��diff.p[1]�ɃX�g�A����
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

					// sum�̌v�Z���I������̂�pos.state()->sum�ɔ��f�����Ă����B(���ꂪ���̊֐��̕Ԃ��l�ɑ�������B)
					st->eval = diff;

				}
				else {

					// ���ȊO�̋�ړ������P�[�X
					// ����̍������v�Z���āA�����ɉ��Z����B

					const int listIndex = mpi.movedPieceNumber[0];

					Evaluate diff = others_moved_diff(pos, mpi.changedPiece[0].newPiece);
					if (numPieceMoved == 1) {

						// �������1�B
						list0[listIndex] = mpi.changedPiece[0].oldPiece.fb;
						list1[listIndex] = mpi.changedPiece[0].oldPiece.fw;
						diff -= others_moved_diff(pos, mpi.changedPiece[0].oldPiece);

					}
					else {

						// �������2�B

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

					// �Onode����̋��̑��������Z�B
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
	// �܂��]���l���v�Z����Ă��Ȃ��Ȃ�
	if (!pos.state()->eval.calculated())
		evaluate(pos);
}

} // namespace Eval

// ���̋ǖʂ̎�Ԃ� c���ɂ�����̂Ƃ���Bc�����猩���]���l��Ԃ��B
int32_t Evaluate::sum(const Color c) const 
{
	// NDF(2014)�̎�ԕ]���̎�@�B
	// cf. http://www.computer-shogi.org/wcsc24/appeal/NineDayFever/NDF.txt

	// ��ԂɈˑ����Ȃ��]���l���v
	// p[1][0]�̓�WKPP�Ȃ̂ŕ����̓}�C�i�X�B
	const int32_t scoreBoard = p[0][0] - p[1][0] + p[2][0];
	// ��ԂɈˑ�����]���l���v
	const int32_t scoreTurn = p[0][1] + p[1][1] + p[2][1];

	// ���̊֐��͎�ԑ����猩���]���l��Ԃ��̂�scoreTurn�͕K���v���X

	return (c == BLACK ? scoreBoard : -scoreBoard) + scoreTurn;
}

void Evaluate::init() {
}
