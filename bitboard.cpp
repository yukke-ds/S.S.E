#include <algorithm>

#include "bitboard.h"
#include "misc.h"

#if defined(USE_PEXT)
Bitboard RookAttack[495616];
#else
Bitboard RookAttack[512000];
#endif
Bitboard RookBlockMask[SQUARE_NB];
int RookAttackIndex[SQUARE_NB];

Bitboard BishopAttack[20224];
Bitboard BishopBlockMask[SQUARE_NB];
int BishopAttackIndex[SQUARE_NB];

Bitboard LanceAttack[COLOR_NB][SQUARE_NB][128];
Bitboard KingAttack[SQUARE_NB];
Bitboard GoldAttack[COLOR_NB][SQUARE_NB];
Bitboard SilverAttack[COLOR_NB][SQUARE_NB];
Bitboard KnightAttack[COLOR_NB][SQUARE_NB];
Bitboard PawnAttack[COLOR_NB][SQUARE_NB];

Bitboard SquareBB[SQUARE_NB];
Bitboard FileBB[FILE_NB];
Bitboard RankBB[RANK_NB];
Bitboard LineBB[SQUARE_NB][SQUARE_NB];			// 2地点の直線上が1のビットボード
Bitboard BetweenBB[SQUARE_NB][SQUARE_NB];		// 2地点の線分が1のビットボード
Bitboard InFrontBB[COLOR_NB][RANK_NB];			// RANK_NBより前が1のビットボード

Bitboard CheckerCandidate[PIECE_NB][SQUARE_NB];
Bitboard AdjacentCheckCandidate[PIECE_NB][SQUARE_NB];
Direction SquareRelation[SQUARE_NB][SQUARE_NB];	// 2地点の位置関係を入れるテーブル

// 関数プロトタイプ宣言
Bitboard index_to_occupied(const int, const int, const Bitboard&);
Bitboard calc_attack_mask(Square, PieceType);
Bitboard lance_block_mask(const Square);
Bitboard attack_calc(const Square, const Bitboard&, const PieceType);
Bitboard lance_attack_calc(const Color, const Square, const Bitboard&);
Bitboard compute_checker_candidates(const Piece, const Square, const bool);

// ビットボード関連の初期化
void Bitboards::init() {

	// SquareBB
	for (Square s = SQ_11; s <= SQ_99; ++s)
	{
		Rank r = rank_of(s);
		File f = file_of(s);
		SquareBB[s].p[0] = (f <= FILE_7) ? ((uint64_t)1 << (f * 9 + r)) : 0;
		SquareBB[s].p[1] = (f >= FILE_8) ? ((uint64_t)1 << ((f - FILE_8) * 9 + r)) : 0;
	}

	// FileBB
	for (File f = FILE_1; f <= FILE_9; ++f)
		FileBB[f] = f > FILE_1 ? (f <= FILE_7 ? FileBB[f - 1] << 9 : File8BB << ((f - FILE_8) * 9)) : File1BB;

	// RankBB
	for (Rank r = RANK_1; r <= RANK_9; ++r)
		RankBB[r] = r > RANK_1 ? RankBB[r - 1] << 1 : Rank1BB;

	// InFrontBB
	for (Rank r = RANK_1; r < RANK_9; ++r)
		InFrontBB[WHITE][r] = ~(InFrontBB[BLACK][r + 1] = InFrontBB[BLACK][r] | RankBB[r]);

	// SquareRelation
	for (Square s1 = SQ_11; s1 <= SQ_99; ++s1)
	{
		File file1 = file_of(s1);
		Rank rank1 = rank_of(s1);
		for (Square s2 = SQ_11; s2 <= SQ_99; ++s2) {
			File file2 = file_of(s2);
			Rank rank2 = rank_of(s2);
			SquareRelation[s1][s2] = DIREC_MISC;
			if (s1 == s2) continue;

			if (file1 == file2)
				SquareRelation[s1][s2] = DIREC_FILE;
			else if (rank1 == rank2)
				SquareRelation[s1][s2] = DIREC_RANK;
			else if (static_cast<int>(rank1 - rank2) == static_cast<int>(file1 - file2))
				SquareRelation[s1][s2] = DIREC_DIAG_NESW;
			else if (static_cast<int>(rank1 - rank2) == static_cast<int>(file2 - file1))
				SquareRelation[s1][s2] = DIREC_DIAG_NWSE;
		}
	}

	// 飛角
	for (PieceType pt : { BISHOP, ROOK })
	{
		auto* attacks = (pt == BISHOP) ? BishopAttack : RookAttack;
		auto* attackIndex = (pt == BISHOP) ? BishopAttackIndex : RookAttackIndex;
		auto* blockMask = (pt == BISHOP) ? BishopBlockMask : RookBlockMask;
		auto* shift = (pt == BISHOP) ? BishopShiftBits : RookShiftBits;
#if defined(USE_PEXT)
#else
		auto* magic = (pt == BISHOP) ? BishopMagic : RookMagic;
#endif

		int index = 0;
		for (Square s = SQ_11; s <= SQ_99; ++s) {
			attackIndex[s] = index;
			blockMask[s] = calc_attack_mask(s, pt);

			assert(!(blockMask[s].p[0] & blockMask[s].p[1]));

			const int bits = blockMask[s].pop_count();
			const int num = 1 << bits;

			for (int i = 0; i < num; ++i) {
				Bitboard occupied = index_to_occupied(i, bits, blockMask[s]);
#if defined(USE_PEXT)
				attacks[index + occupied_to_index(occupied & blockMask[s], blockMask[s])] = attackCalc(s, occupied, pt);
#else
				attacks[index + occupied_to_index(occupied, magic[s], shift[s])] = attack_calc(s, occupied, pt);
#endif
			}

			index += 1 << (64 - shift[s]);
		}
	}

	// 玉
	for (Square s = SQ_11; s <= SQ_99; ++s)
		KingAttack[s] = rook_attacks(s, ALL1BB) | bishop_attacks(s, ALL1BB);

	// 金
	for (Color c = BLACK; c <= WHITE; ++c)
		for (Square s = SQ_11; s <= SQ_99; ++s)
			GoldAttack[c][s] = (king_attacks(s) & InFrontBB[c][rank_of(s)]) | rook_attacks(s, ALL1BB);

	// 銀
	for (Color c = BLACK; c <= WHITE; ++c)
		for (Square s = SQ_11; s <= SQ_99; ++s)
			SilverAttack[c][s] = (king_attacks(s) & InFrontBB[c][rank_of(s)]) | bishop_attacks(s, ALL1BB);

	// 歩
	for (Color c = BLACK; c <= WHITE; ++c)
		for (Square s = SQ_11; s <= SQ_99; ++s)
			PawnAttack[c][s] = silver_attacks(c, s) ^ bishop_attacks(s, ALL1BB);

	// 桂
	for (Color c = BLACK; c <= WHITE; ++c)
		for (Square s = SQ_11; s <= SQ_99; ++s) 
		{
			KnightAttack[c][s] = ALL0BB;
			const Bitboard b = pawn_attacks(c, s);
			if (b)
				KnightAttack[c][s] = bishop_step_attacks(b.pop_c()) & InFrontBB[c][rank_of(s)];
		}

	// 香
	for (Color c = BLACK; c <= WHITE; ++c)
		for (Square s = SQ_11; s <= SQ_99; ++s)
		{
			const Bitboard blockMask = lance_block_mask(s);
			const int bits = blockMask.pop_count();
			const int num = 1 << bits;

			for (int i = 0; i < num; ++i) 
			{
				Bitboard occupied = index_to_occupied(i, bits, blockMask);
				LanceAttack[c][s][i] = lance_attack_calc(c, s, occupied);
			}
		}

	// BetweenBBとLineBB
	for (Square s1 = SQ_11; s1 <= SQ_99; ++s1)
		for (Square s2 = SQ_11; s2 <= SQ_99; ++s2)
		{
			BetweenBB[s1][s2] = ALL0BB;
			if (s1 == s2) continue;
			const Direction direc = square_relation(s1, s2);
			if (direc & DIREC_CROSS)
			{
				BetweenBB[s1][s2] = rook_attacks(s1, SquareBB[s2]) & rook_attacks(s2, SquareBB[s1]);
				LineBB[s1][s2] = (rook_attacks(s1, ALL0BB) & rook_attacks(s2, ALL0BB)) | s1 | s2;
			}
			else if (direc & DIREC_DIAG)
			{
				BetweenBB[s1][s2] = bishop_attacks(s1, SquareBB[s2]) & bishop_attacks(s2, SquareBB[s1]);
				LineBB[s1][s2] = (bishop_attacks(s1, ALL0BB) & bishop_attacks(s2, ALL0BB)) | s1 | s2;
			}
		}

	// CheckerCandidateとAdjacentCheckCandidate
	for (Color c = BLACK; c <= WHITE; ++c)
		for (PieceType pt = PAWN; pt <= DRAGON; ++pt)
		{
			Piece pc = make_piece(c, pt);
			for (Square s = SQ_11; s <= SQ_99; ++s)
			{
				CheckerCandidate[pc][s] = compute_checker_candidates(pc, s, false);
				AdjacentCheckCandidate[pc][s] = compute_checker_candidates(pc, s, true);
			}
		}
}

Bitboard index_to_occupied(const int index, const int bits, const Bitboard& mask_) {

	auto mask = mask_;
	auto result = ALL0BB;

	for (int i = 0; i < bits; ++i)
	{
		const Square s = mask.pop();
		if (index & (1 << i)) result ^= s;
	}

	return result;
}

// square のマスにおける、障害物を調べる必要がある場所を調べて Bitboard で返す。
Bitboard calc_attack_mask(Square s, PieceType pt) {

	Bitboard result;
	if (pt == BISHOP) {

		result = ALL0BB;

		for (Rank r = RANK_2; r <= RANK_8; ++r) // 外周は除外
			for (File f = FILE_2; f <= FILE_8; ++f)
				if (abs(rank_of(s) - r) == abs(file_of(s) - f))
					result ^= (f | r);
	}
	else {
		assert(pt == ROOK);

		result = RankBB[rank_of(s)] ^ FileBB[file_of(s)];

		// 外周に居ない限り、その外周升は利きの計算には関係ない。
		if (file_of(s) != FILE_1) { result &= ~File1BB; }
		if (file_of(s) != FILE_9) { result &= ~File9BB; }
		if (rank_of(s) != RANK_1) { result &= ~Rank1BB; }
		if (rank_of(s) != RANK_9) { result &= ~Rank9BB; }
	}

	// sqの地点は関係ないのでクリアしておく
	result &= ~Bitboard(s);

	return result;
}

Bitboard lance_block_mask(const Square s) {

	return FileBB[file_of(s)] & ~(Rank1BB | Rank9BB);
}

Bitboard attack_calc(const Square s, const Bitboard& occupied, const PieceType pt) {

	Bitboard result = ALL0BB;
	const Square DeltaArray[2][4] = { { DELTA_NE, DELTA_SE, DELTA_SW, DELTA_NW }, { DELTA_N, DELTA_S, DELTA_E, DELTA_W } };
	for (Square delta : DeltaArray[(pt == BISHOP) ? 0 : 1])
	{
		for (Square sq = s + delta; is_ok(sq) && abs(rank_of(sq - delta) - rank_of(sq)) <= 1; sq += delta)
		{
			result ^= sq;
			if (occupied & SquareBB[sq]) break;
		}
	}
	return result;
}

Bitboard lance_attack_calc(const Color c, const Square s, const Bitboard& occupied) {
	return rook_attacks(s, occupied) & InFrontBB[c][rank_of(s)];
}

Bitboard compute_checker_candidates(const Piece pc, const Square s, const bool isAdjacentCheckOnly) {

	if (pc == B_KING || pc == W_KING)
		return ALL0BB;

	Bitboard neighborhood12 = (king_attacks(s) | knight_attacks(BLACK, s) | knight_attacks(WHITE, s));
	Bitboard target = isAdjacentCheckOnly ? neighborhood12 : ALL1BB;
	Bitboard result = ALL0BB;

	// 1. 不成の王手
	Piece opponentPiece = opponent_piece(pc);
	Bitboard nonPromotionTarget = attacks_bb(opponentPiece, s, ALL0BB) & target;
	nonPromotionTarget.for_each([&](Square to) {
		result |= attacks_bb(opponentPiece, to, ALL0BB);
	});

	// 2. 成る王手
	if (can_promote(pc))
	{
		Color c = color_of(pc);
		Piece promotedPiece = promoted_piece(opponentPiece);
		Bitboard promotionTarget = attacks_bb(promotedPiece, s, ALL0BB) & target;
		promotionTarget.for_each([&](Square to) {
			if (is_promotion_area(c, to))
				result |= attacks_bb(opponentPiece, to, ALL0BB);
			else
				result |= attacks_bb(opponentPiece, to, ALL0BB) & promotion_area(c);
		});
	}

	// 3. すでに玉に効きをつけている駒は、候補から除く
	result &= ~attacks_bb(opponentPiece, s, ALL1BB);

	// 4. 受け方の玉があるマスは、取り除いておく
	result &= ~Bitboard(s);

	return result;
}