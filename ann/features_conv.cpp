#include "features_conv.h"

#include "../progress.h"

namespace
{

using namespace FeaturesConv;

struct AttackMaps
{
	PieceType blackLeastValuableAttackers[SQUARE_NB];
	PieceType whiteLeastValuableAttackers[SQUARE_NB];

	uint8_t blackNumAttackers[SQUARE_NB];
	uint8_t whiteNumAttackers[SQUARE_NB];
	
	// normalized maximum value of pieces white and black can put on each square
	float blackCtrl[SQUARE_NB];
	float whiteCtrl[SQUARE_NB];
	
	// is it safe to move a piece of type pt to sq?
	bool is_safe(Color c, PieceType pt, Square sq)
	{
		// here we are doing a very simple form of SEE
		// if the opponent has no attacker, the piece is safe
		// if the opponent has an attacker and it's lower valued, we are not safe
		// if the opponent has an attacker and it's equal or higher valued, we are
		// safe as long as we also have an attacker (that's not ourselves)
		
		// we don't have to worry about winning captures, because qsearch will take care of that
		// here we are only looking at moving to empty squares
		PieceType opponentAttacker = (c == BLACK) ? whiteLeastValuableAttackers[sq] : blackLeastValuableAttackers[sq];
		PieceType friendlyAttacker = (c == BLACK) ? blackLeastValuableAttackers[sq] : whiteLeastValuableAttackers[sq];
		uint8_t numFriendlyAttackers = (c == BLACK) ? blackNumAttackers[sq] : whiteNumAttackers[sq];
		
		if (opponentAttacker == NO_PIECE_TYPE)
			return true;
		else if (Eval::PieceValue[opponentAttacker] < Eval::PieceValue[pt])
			return false;
		else
			return friendlyAttacker != NO_PIECE_TYPE && (numFriendlyAttackers > 1);
	}
};

float normalize_coord(int x)
{
	// map x from 1 - 9 to 0 - 1
	return 0.111 * x;
}

float normalize_count(int x, int typicalMaxCount)
{
	return static_cast<float>(x) / static_cast<float>(typicalMaxCount);
}

template <typename T> void push_global_bool(std::vector<T> &ret, bool x, int32_t group);
template<> void push_global_bool<float>(std::vector<float> &ret, bool x, int32_t /*group*/)
{
	if (x)
		ret.push_back(1.0f);
	else
		ret.push_back(0.0f);
}

template<> void push_global_bool<FeatureDescription>(std::vector<FeatureDescription> &ret, bool /*x*/, int32_t group)
{
	FeatureDescription fd;
	fd.featureType = FeatureDescription::FeatureType_global;
	fd.group = group;
	ret.push_back(fd);
}

template <typename T> void push_global_float(std::vector<T> &ret, float x, int32_t group);
template<> void push_global_float<float>(std::vector<float> &ret, float x, int32_t /*group*/)
{
	ret.push_back(x);
}

template<> void push_global_float<FeatureDescription>(std::vector<FeatureDescription> &ret, float /*x*/, int32_t group)
{
	FeatureDescription fd;
	fd.featureType = FeatureDescription::FeatureType_global;
	fd.group = group;
	ret.push_back(fd);
}

template <typename T>
void push_global_coords(std::vector<T> &ret, bool exists, Square sq, int32_t group, bool mustExist = false)
{
	if (!mustExist)
		push_global_bool(ret, exists, group);

	File f = file_of(sq);
	Rank r = rank_of(sq);

	push_global_float(ret, exists ? normalize_coord(f + 1) : 0.0f, group);
	push_global_float(ret, exists ? normalize_coord(r + 1) : 0.0f, group);
}

template <typename T> void push_mobility(
	std::vector<T> &ret,
	float mob,
	int32_t group);
template<> void push_mobility<float>(
	std::vector<float> &ret,
	float mob,
	int32_t /*group*/)
{
	ret.push_back(mob);
}

template<> void push_mobility<FeatureDescription>(
	std::vector<FeatureDescription> &ret,
	float /*mob*/,
	int32_t group)
{
	FeatureDescription fd;

	fd.featureType = FeatureDescription::FeatureType_global;
	fd.group = group;

	ret.push_back(fd);
}

template <typename T> void push_pos_float(std::vector<T> &ret, Square sq, float x, int32_t group);
template<> void push_pos_float<float>(std::vector<float> &ret, Square /*sq*/, float x, int32_t /*group*/)
{
	ret.push_back(x);
}

template<> void push_pos_float<FeatureDescription>(std::vector<FeatureDescription> &ret, Square sq, float /*x*/, int32_t group)
{
	FeatureDescription fd;
	fd.featureType = FeatureDescription::FeatureType_pos;
	fd.sq = sq;
	fd.group = group;
	ret.push_back(fd);
}

template <Color Us, typename T>
void push_attacks(std::vector<T> &ret, Square sq, PieceType pt, bool exists, const Position &pos, AttackMaps &atkMaps, int32_t group)
{
	int32_t safeMovesCount = 0;

	File fStart = file_of(sq);
	Rank rStart = rank_of(sq);
	
	// 香車
	if (pt == LANCE)
	{
		// figure out how far we can go in each direction
		Rank up = (Us == BLACK ? -RANK_2 : RANK_2);
		int32_t count = 0;
		File f = fStart;
		Rank r = rStart + up;

		while (is_ok(f) && is_ok(r) && exists)
		{
			++count;

			if (atkMaps.is_safe(Us, pt, f | r))
				++safeMovesCount;
		
			if (!pos.empty(f | r))
				break;

			r += up;
		}
	
		push_mobility(ret, normalize_count(count, 8), group);
	}

	// 角と馬
	else if (pt == BISHOP || pt == HORSE)
	{
		// figure out how far we can go in each direction
		// 左下，右下，左上，右上
		const File DirFileOffsets[4] = { FILE_2, -FILE_2, FILE_2, -FILE_2 };
		const Rank DirRankOffsets[4] = { RANK_2, RANK_2, -RANK_2, -RANK_2 };
		
		for (int32_t i = 0; i < 4; ++i)
		{
			// for each direction, keep going until we hit either the board edge or another piece
			// sqから斜めそれぞれ，どこまで(いくつ)利きがあるかを調べる
			int32_t count = 0;
			File f = fStart + DirFileOffsets[i];
			Rank r = rStart + DirRankOffsets[i];
			
			while (is_ok(f) && is_ok(r) && exists)
			{				
				++count;
				
				if (atkMaps.is_safe(Us, pt, f | r))
					++safeMovesCount;
				
				if (!pos.empty(f | r))
					break;

				f += DirFileOffsets[i];
				r += DirRankOffsets[i];
			}

			push_mobility(ret, normalize_count(count, 8), group);
		}
	}

	// 飛車と龍
	else if (pt == ROOK || pt == DRAGON)
	{
		// figure out how far we can go in each direction
		// 左，右，下，上
		const File DirFileOffsets[4] = { FILE_2, -FILE_2, FILE_1, FILE_1 };
		const Rank DirRankOffsets[4] = { RANK_1, RANK_1, RANK_2, -RANK_2 };
		
		for (int32_t i = 0; i < 4; ++i)
		{
			// for each direction, keep going until we hit either the board edge or another piece
			// sqから斜めそれぞれ，どこまで(いくつ)利きがあるかを調べる
			int32_t count = 0;
			File f = fStart + DirFileOffsets[i];
			Rank r = rStart + DirRankOffsets[i];
			
			while (is_ok(f) && is_ok(r) && exists)
			{				
				++count;
				
				if (atkMaps.is_safe(Us, pt, f | r))
					++safeMovesCount;
			
				if (!pos.empty(f | r))
					break;

				f += DirFileOffsets[i];
				r += DirRankOffsets[i];
			}
			
			push_mobility(ret, normalize_count(count, 8), group);
		}
	}

	// 歩
	//else if (pt == PAWN) {}

	// 桂馬
	else if (pt == KNIGHT)
	{
		// 先手:左上，右上　後手:左下，右下
		const File DirFileOffsets[2] = { FILE_2, -FILE_2 };
		const Rank DirRankOffsets[2] = { -RANK_3, -RANK_3 };
		int32_t count = 0;
		
		for (int i = 0; i < 2; ++i)
		{
			File f = fStart + DirFileOffsets[i];
			Rank r = rStart + (Us == BLACK ? DirRankOffsets[i] : -DirRankOffsets[i]);
			
			if (is_ok(f) && is_ok(r) && exists)
			{				
				count++;
				
				if (atkMaps.is_safe(Us, pt, f | r))
					++safeMovesCount;

				if (!pos.empty(f | r))
					continue;
			}
		}
		
		push_mobility(ret, normalize_count(count, 2), group);
	}

	// 銀
	else if (pt == SILVER)
	{
		// 先手:左上，上，右上，左下，右下　後手:左上，下，右上，左下，右下
		const File DirFileOffsets[5] = { FILE_2, FILE_1, -FILE_2, FILE_2, -FILE_2 };
		const Rank DirRankOffsets[5] = { -RANK_2, -RANK_2, -RANK_2, RANK_2, RANK_2 };
		int32_t count = 0;
		
		for (int i = 0; i < 5; ++i)
		{
			File f = fStart + DirFileOffsets[i];
			Rank r = rStart + ((Us == WHITE && i == 1) ? -DirRankOffsets[i] : DirRankOffsets[i]);
			
			if (is_ok(f) && is_ok(r) && exists)
			{				
				count++;
				
				if (atkMaps.is_safe(Us, pt, f | r))
					++safeMovesCount;

				if (!pos.empty(f | r))
					continue;
			}
		}
		
		push_mobility(ret, normalize_count(count, 5), group);
	}
	
	// 金
	else if (pt == GOLD)
	{
		// 先手:左上，上，右上，左，右，下　後手:左下，上，右下，左，右，下
		const File DirFileOffsets[6] = { FILE_2, FILE_1, -FILE_2, FILE_2, -FILE_2, FILE_1 };
		const Rank DirRankOffsets[6] = { -RANK_2, -RANK_2, -RANK_2, RANK_1, RANK_1, RANK_2 };
		int32_t count = 0;
		
		for (int i = 0; i < 6; ++i)
		{
			File f = fStart + DirFileOffsets[i];
			Rank r = rStart + ((Us == WHITE && (i == 0 || i == 2)) ? -DirRankOffsets[i] : DirRankOffsets[i]);
			
			if (is_ok(f) && is_ok(r) && exists)
			{				
				count++;
			
				if (atkMaps.is_safe(Us, pt, f | r))
					++safeMovesCount;

				if (!pos.empty(f | r))
					continue;
			}
		}

		push_mobility(ret, normalize_count(count, 6), group);
	}
	
	// 16 is the "reasonably maximum", though queens in the centre of an empty board can have up to 27. That's fine.
	// 最大利き:香(8)，桂(2)，銀(5)，金(6)，飛角(16)，馬龍(20) Total:57
	push_mobility(ret, normalize_count(exists ? safeMovesCount : 0, 28), group);
}

template <typename T>
void push_square_features(std::vector<T> &ret, const Position &/*pos*/, AttackMaps &atkMaps, int32_t &group)
{
	for (Square s = SQ_11; s <= SQ_99; ++s)
	{
		push_pos_float(ret, s, atkMaps.blackCtrl[s], group);
		push_pos_float(ret, s, atkMaps.whiteCtrl[s], group + 1);
	}

	group += 2;
}

template <Color Us, typename T>
void push_pawns(std::vector<T> &ret, Bitboard pawns, AttackMaps &atkMaps, int32_t &group)
{
	std::tuple<bool, Square> assignments[FILE_NB];

	for (File f = FILE_1; f <= FILE_9; ++f)
		std::get<0>(assignments[f]) = false;
	
	// in the first pass, we assign each pawn to the corresponding file if possible,
	// and keep a list (in a bitboard) of pawns that still need to be assigned
	while (pawns)
	{
		Square sq = pawns.pop();

		File f = file_of(sq);
		
		if (std::get<0>(assignments[f]) == false)
		{
			std::get<0>(assignments[f]) = true;
			std::get<1>(assignments[f]) = sq;
		}
		else
			assert(false); // 二歩局面はあってはならない
	}

	for (File f = FILE_1; f <= FILE_9; ++f)
	{
		bool exists = std::get<0>(assignments[f]);
		Square sq = std::get<1>(assignments[f]);

		push_global_coords(ret, exists, sq, group);
		push_threat(ret, sq, exists, atkMaps, group);
	}
}

template <typename T>
void push_threat(
	std::vector<T> &ret,
	Square sq,
	bool exists,
	AttackMaps &atkMaps,
	int32_t group)
{
	if (exists)
	{
		// we push both black and white control because one would be defending the piece,
		// and one attacking
		push_global_float(ret, atkMaps.blackCtrl[sq], group);
		push_global_float(ret, atkMaps.whiteCtrl[sq], group);
	}
	else
	{
		push_global_float(ret, 0.0f, group);
		push_global_float(ret, 0.0f, group);
	}
}

template <Color Us, typename T>
void push_pieces(
	std::vector<T> &ret,
	Bitboard pieces,
	PieceType pt,
	const Position &pos,
	int32_t group,
	AttackMaps &atkMaps)
{
	// 歩~龍までの存在可能枚数
	const size_t MaxPieces[PIECE_TYPE_NB] = {
		0, 18, 4, 4, 4, 2, 2, 4, 1, 18, 4, 4, 4, 2, 2
	};

	// 香(，桂，銀，金)，角，飛(，成駒)，馬，龍の座標と利きを作る
	std::vector<std::tuple<bool, Square>> assignments(MaxPieces[pt]);

	for (size_t i = 0; i < MaxPieces[pt]; ++i)
		std::get<0>(assignments[i]) = false;

	for (size_t i = 0; pieces; ++i)
	{
		Square sq = pieces.pop();

		std::get<0>(assignments[i]) = true;
		std::get<1>(assignments[i]) = sq;
	}

	for (size_t i = 0; i < MaxPieces[pt]; ++i)
	{
		bool exists = std::get<0>(assignments[i]);
		Square sq = std::get<1>(assignments[i]);

		push_global_coords(ret, exists, sq, group);
		push_attacks<Us>(ret, sq, pt, exists, pos, atkMaps, group);
		push_threat(ret, sq, exists, atkMaps, group);
	}
}

void compute_least_valuable_attackers(const Position& pos, PieceType attackers[64], uint8_t numAttackers[64], Color stm)
{
	// initialize everything to empty
	for (Square s = SQ_11; s <= SQ_99; ++s)
	{
		attackers[s] = NO_PIECE_TYPE;
		numAttackers[s] = 0;
	}

	// now we start from the most valuable and go to least, and just keep overwriting
	// 盤上の利きを駒価値の高いものから低いものへ上書きしていく
	// 1. 玉
	Bitboard toBB = king_attacks(pos.king_square(stm));
	toBB.for_each([&](Square to) {
		attackers[to] = KING;
		++numAttackers[to];
	});

	// 2. 龍
	pos.pieces(stm, DRAGON).for_each([&](Square from) {
		toBB = dragon_attacks(from, pos.pieces());
		toBB.for_each([&](Square to) {
			attackers[to] = DRAGON;
			++numAttackers[to];
		});
	});

	// 3. 馬
	pos.pieces(stm, HORSE).for_each([&](Square from) {
		toBB = horse_attacks(from, pos.pieces());
		toBB.for_each([&](Square to) {
			attackers[to] = HORSE;
			++numAttackers[to];
		});
	});

	// 4. 飛車
	pos.pieces(stm, ROOK).for_each([&](Square from) {
		toBB = rook_attacks(from, pos.pieces());
		toBB.for_each([&](Square to) {
			attackers[to] = ROOK;
			++numAttackers[to];
		});
	});

	// 5. 角
	pos.pieces(stm, BISHOP).for_each([&](Square from) {
		toBB = bishop_attacks(from, pos.pieces());
		toBB.for_each([&](Square to) {
			attackers[to] = BISHOP;
			++numAttackers[to];
		});
	});

	// 6. と金
	pos.pieces(stm, PRO_PAWN).for_each([&](Square from) {
		toBB = gold_attacks(stm, from);
		toBB.for_each([&](Square to) {
			attackers[to] = PRO_PAWN;
			++numAttackers[to];
		});
	});

	// 7. 成桂
	pos.pieces(stm, PRO_KNIGHT).for_each([&](Square from) {
		toBB = gold_attacks(stm, from);
		toBB.for_each([&](Square to) {
			attackers[to] = PRO_KNIGHT;
			++numAttackers[to];
		});
	});

	// 8. 成香
	pos.pieces(stm, PRO_LANCE).for_each([&](Square from) {
		toBB = gold_attacks(stm, from);
		toBB.for_each([&](Square to) {
			attackers[to] = PRO_LANCE;
			++numAttackers[to];
		});
	});

	// 9. 成銀
	pos.pieces(stm, PRO_SILVER).for_each([&](Square from) {
		toBB = gold_attacks(stm, from);
		toBB.for_each([&](Square to) {
			attackers[to] = PRO_SILVER;
			++numAttackers[to];
		});
	});

	// 10. 金
	pos.pieces(stm, GOLD).for_each([&](Square from) {
		toBB = gold_attacks(stm, from);
		toBB.for_each([&](Square to) {
			attackers[to] = GOLD;
			++numAttackers[to];
		});
	});

	// 11. 銀
	pos.pieces(stm, SILVER).for_each([&](Square from) {
		toBB = silver_attacks(stm, from);
		toBB.for_each([&](Square to) {
			attackers[to] = SILVER;
			++numAttackers[to];
		});
	});

	// 12. 桂馬
	pos.pieces(stm, KNIGHT).for_each([&](Square from) {
		toBB = knight_attacks(stm, from);
		toBB.for_each([&](Square to) {
			attackers[to] = KNIGHT;
			++numAttackers[to];
		});
	});

	// 13. 香車
	pos.pieces(stm, LANCE).for_each([&](Square from) {
		toBB = lance_attacks(stm, from, pos.pieces());
		toBB.for_each([&](Square to) {
			attackers[to] = LANCE;
			++numAttackers[to];
		});
	});

	// 14. 歩
	Square up = (stm == BLACK ? DELTA_N : DELTA_S);
	toBB = shift_bb(pos.pieces(stm, PAWN), up);
	toBB.for_each([&](Square to) {
		attackers[to] = PAWN;
		++numAttackers[to];
	});
}

AttackMaps compute_attack_maps(const Position& pos)
{
	AttackMaps ret;

	compute_least_valuable_attackers(pos, ret.blackLeastValuableAttackers, ret.blackNumAttackers, BLACK);
	compute_least_valuable_attackers(pos, ret.whiteLeastValuableAttackers, ret.whiteNumAttackers, WHITE);

	// convert them to control values
	for (Square s = SQ_11; s <= SQ_99; ++s)
	{
		PieceType blackPt = ret.blackLeastValuableAttackers[s];
		PieceType whitePt = ret.whiteLeastValuableAttackers[s];

		// if a side doesn't attack the square, control is 0
		// if a side attacks with a piece, control is higher the lower valued the piece is
		ret.blackCtrl[s] = (blackPt == NO_PIECE_TYPE) ? 0.0f : normalize_count(Eval::PieceValue[KING] + Eval::PieceValue[KING] / 2 - Eval::PieceValue[blackPt], Eval::PieceValue[KING] * 2);
		ret.whiteCtrl[s] = (whitePt == NO_PIECE_TYPE) ? 0.0f : normalize_count(Eval::PieceValue[KING] + Eval::PieceValue[KING] / 2 - Eval::PieceValue[whitePt], Eval::PieceValue[KING] * 2);
	}

	return ret;
}

} // namespace

namespace FeaturesConv
{

template <typename T>
void convert_pos_to_NN(const Position &pos, std::vector<T> &ret)
{
	ret.clear(); // this shouldn't actually deallocate memory

	// we start by computing values that will be used later
	int BHPCount = count_of(pos.hand_of(BLACK), PAWN);
	int BHLCount = count_of(pos.hand_of(BLACK), LANCE);
	int BHNCount = count_of(pos.hand_of(BLACK), KNIGHT);
	int BHSCount = count_of(pos.hand_of(BLACK), SILVER);
	int BHGCount = count_of(pos.hand_of(BLACK), GOLD);
	int BHBCount = count_of(pos.hand_of(BLACK), BISHOP);
	int BHRCount = count_of(pos.hand_of(BLACK), ROOK);

	int WHPCount = count_of(pos.hand_of(WHITE), PAWN);
	int WHLCount = count_of(pos.hand_of(WHITE), LANCE);
	int WHNCount = count_of(pos.hand_of(WHITE), KNIGHT);
	int WHSCount = count_of(pos.hand_of(WHITE), SILVER);
	int WHGCount = count_of(pos.hand_of(WHITE), GOLD);
	int WHBCount = count_of(pos.hand_of(WHITE), BISHOP);
	int WHRCount = count_of(pos.hand_of(WHITE), ROOK);

	// Attack and Defend Maps
	AttackMaps atkMaps = compute_attack_maps(pos);

	// now we can start actually forming the groups
	int32_t group = 0;

	// 1. Global Features (6 inputs)
	// チェスと違い，将棋は駒の数で進行度を把握できない．
	//push_global_float(ret, Progress::estimate_progress(pos), group);

	// which side to move
	push_global_bool(ret, pos.side_to_move() == BLACK, group);

	// king positions
	push_global_coords(ret, true, pos.king_square(BLACK), group, true);
	push_global_coords(ret, true, pos.king_square(WHITE), group, true);

	// 2. Piece-Centric Features (472 inputs)
	// 利き(Pieces Mobility)は基本的に香角飛馬龍のみ．歩桂銀金はあとで実験．
	// pawns (all pawns are in the same group) (90 inputs)
	++group;
	push_pawns<BLACK>(ret, pos.pieces(BLACK, PAWN), atkMaps, group);
	push_global_float(ret, normalize_count(BHPCount, 6), group);
	push_pawns<WHITE>(ret, pos.pieces(WHITE, PAWN), atkMaps, group);
	push_global_float(ret, normalize_count(WHPCount, 6), group);

	// lances (58 inputs)
	++group;
	push_pieces<BLACK>(ret, pos.pieces(BLACK, LANCE), LANCE, pos, group, atkMaps);
	push_global_float(ret, normalize_count(BHLCount, 4), group);
	push_pieces<WHITE>(ret, pos.pieces(WHITE, LANCE), LANCE, pos, group, atkMaps);
	push_global_float(ret, normalize_count(WHLCount, 4), group);

	// knights (58 inputs)
	++group;
	push_pieces<BLACK>(ret, pos.pieces(BLACK, KNIGHT), KNIGHT, pos, group, atkMaps);
	push_global_float(ret, normalize_count(BHNCount, 4), group);
	push_pieces<WHITE>(ret, pos.pieces(WHITE, KNIGHT), KNIGHT, pos, group, atkMaps);
	push_global_float(ret, normalize_count(WHNCount, 4), group);

	// silvers (58 inputs)
	++group;
	push_pieces<BLACK>(ret, pos.pieces(BLACK, SILVER), SILVER, pos, group, atkMaps);
	push_global_float(ret, normalize_count(BHSCount, 4), group);
	push_pieces<WHITE>(ret, pos.pieces(WHITE, SILVER), SILVER, pos, group, atkMaps);
	push_global_float(ret, normalize_count(WHSCount, 4), group);

	// golds (58 inputs) TODO:成駒どうするか
	++group;
	push_pieces<BLACK>(ret, pos.pieces(BLACK, GOLD), GOLD, pos, group, atkMaps);
	push_global_float(ret, normalize_count(BHGCount, 4), group);
	push_pieces<WHITE>(ret, pos.pieces(WHITE, GOLD), GOLD, pos, group, atkMaps);
	push_global_float(ret, normalize_count(WHGCount, 4), group);

	// bishops (42 inputs)
	++group;
	push_pieces<BLACK>(ret, pos.pieces(BLACK, BISHOP), BISHOP, pos, group, atkMaps);
	push_global_float(ret, normalize_count(BHBCount, 2), group);
	push_pieces<WHITE>(ret, pos.pieces(WHITE, BISHOP), BISHOP, pos, group, atkMaps);
	push_global_float(ret, normalize_count(WHBCount, 2), group);

	// rooks (42 inputs)
	++group;
	push_pieces<BLACK>(ret, pos.pieces(BLACK, ROOK), ROOK, pos, group, atkMaps);
	push_global_float(ret, normalize_count(BHRCount, 2), group);
	push_pieces<WHITE>(ret, pos.pieces(WHITE, ROOK), ROOK, pos, group, atkMaps);
	push_global_float(ret, normalize_count(WHRCount, 2), group);

	// horses (40 inputs)
	++group;
	push_pieces<BLACK>(ret, pos.pieces(BLACK, HORSE), HORSE, pos, group, atkMaps);
	push_pieces<WHITE>(ret, pos.pieces(WHITE, HORSE), HORSE, pos, group, atkMaps);

	// dragons (40 inputs)
	++group;
	push_pieces<BLACK>(ret, pos.pieces(BLACK, DRAGON), DRAGON, pos, group, atkMaps);
	push_pieces<WHITE>(ret, pos.pieces(WHITE, DRAGON), DRAGON, pos, group, atkMaps);

	// 3. Square-Centric Features (162 input)
	push_square_features(ret, pos, atkMaps, group);

	// total 656 inputs
}

void convert_sfen_to_BB(const std::string& board, const std::string& stm, const std::string& hand, const std::string& /*gamePly*/, float ret[PositionInputDims])
{
	unsigned char token;
	bool promotion = false;
	const std::string PieceToChar = "PLNSGBRK      plnsgbrk";
	size_t idx;
	File f = FILE_1;
	Rank r = RANK_1;

	if (PositionInputDims == 2282)
	{
		// 1. 駒の配置
		for (size_t i = 0; i < board.size(); ++i)
		{
			token = board[i];

			if (isdigit(token))
				f += File(token - '0');
			
			else if (token == '/')
			{
				f = FILE_1;
				++r;
			}

			else if (token == '+')
				promotion = true;
			
			else if ((idx = PieceToChar.find(token)) != std::string::npos)
			{
				// 角飛の成りは金の成りがない関係で-1して詰める
				bool HDCheck = (idx == 5 || idx == 6 || idx == 19 || idx == 20);
				size_t place = (idx + (promotion ? (HDCheck ? PIECE_PROMOTION - 1 : PIECE_PROMOTION) : 0)) * 81;
				ret[place + (r * 9 + f)] = 1.0f;
				
				promotion = false;
				++f;
			}
		}
		
		// 2. 持ち駒
		int pn = 0;
		for (size_t i = 0; i < hand.size(); ++i)
		{
			token = hand[i];
			
			if (token == '-')
				break;
			
			if (isdigit(token))
				pn = (token - '0') + pn * 10;
			
			else if ((idx = PieceToChar.find(token)) != std::string::npos)
			{
				pn = std::max(pn, 1);
				ret[BoardInputDims + (idx < 7 ? idx : idx - 7)] = float(pn);
				
				pn = 0;
			}
		}

		// 3. 手番
		//ret[PositionInputDims - 1] = (stm == "b" ? 0.0f : 1.0f);
	}
	else if (PositionInputDims == 518)
	{
		// 1. 駒の配置
		for (size_t i = 0; i < board.size(); ++i)
		{
			token = board[i];

			if (isdigit(token))
				f += File(token - '0');

			else if (token == '/')
			{
				f = FILE_1;
				++r;
			}

			else if (token == '+')
				promotion = true;

			else if ((idx = PieceToChar.find(token)) != std::string::npos)
			{
				// 角飛の成りは金の成りがない関係で-1して詰める
				bool HDCheck = (idx == 5 || idx == 6 || idx == 19 || idx == 20);
				size_t place = (idx + (promotion ? (HDCheck ? PIECE_PROMOTION - 1 : PIECE_PROMOTION) : 0)) * 18;
				ret[place + f] = 1.0f;
				ret[place + 9 + r] = 1.0f;

				promotion = false;
				++f;
			}
		}

		// 2. 持ち駒
		int pn = 0;
		for (size_t i = 0; i < hand.size(); ++i)
		{
			token = hand[i];
			
			if (token == '-')
				break;
			
			if (isdigit(token))
				pn = (token - '0') + pn * 10;
			
			else if ((idx = PieceToChar.find(token)) != std::string::npos)
			{
				pn = std::max(pn, 1);
				ret[BoardInputDims + (idx < 7 ? idx : idx - 7)] = float(pn);
				
				pn = 0;
			}
		}

		// 3. 手番
		//ret[PositionInputDims - 1] = (stm == "b" ? 0.0f : 1.0f);
	}
}

// 明示的特殊化
template void convert_pos_to_NN<float>(const Position &pos, std::vector<float> &ret);
template void convert_pos_to_NN<FeatureDescription>(const Position &pos, std::vector<FeatureDescription> &ret);

} // namespace FeaturesConv