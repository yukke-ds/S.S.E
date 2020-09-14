#include "movegen.h"
#include "position.h"

namespace {

ExtMove* generate_moves(const Position& pos, ExtMove* moveList, Color us, Bitboard target) {

	// 1. 歩
	Square up = (us == BLACK ? DELTA_N : DELTA_S);
	Bitboard toBB = shift_bb(pos.pieces(us, PAWN), up) & target;
	toBB.for_each([&](Square to) {
		if (is_promotion_area(us, to))
			*moveList++ = make_move_promote(to - up, to);
		else
			*moveList++ = make_move(to - up, to);
	});

	// 2. 香車
	pos.pieces(us, LANCE).for_each([&](Square from)
	{
		toBB = lance_attacks(us, from, pos.pieces()) & target;
		toBB.for_each([&](Square to) {
			if (is_promotion_area(us, to))
				*moveList++ = make_move_promote(from, to);
			if (relative_rank(us, rank_of(to)) >= RANK_3)
				*moveList++ = make_move(from, to);
		});
	});

	// 3. 桂馬
	pos.pieces(us, KNIGHT).for_each([&](Square from) {
		toBB = knight_attacks(us, from) & target;
		toBB.for_each([&](Square to) {
			if (is_promotion_area(us, to))
				*moveList++ = make_move_promote(from, to);
			if (relative_rank(us, rank_of(to)) >= RANK_3)
				*moveList++ = make_move(from, to);
		});
	});

	// 4. 銀
	pos.pieces(us, SILVER).for_each([&](Square from) {
		toBB = silver_attacks(us, from) & target;
		toBB.for_each([&](Square to) {
			if (is_promotion_area(us, to) || is_promotion_area(us, from))
				*moveList++ = make_move_promote(from, to);
			*moveList++ = make_move(from, to);
		});
	});

	// 5. 金
	(pos.gold_pieces() & pos.pieces(us)).for_each([&](Square from) {
		toBB = gold_attacks(us, from) & target;
		toBB.for_each([&](Square to) {
			*moveList++ = make_move(from, to);
		});
	});

	// 6. 角
	pos.pieces(us, BISHOP).for_each([&](Square from) {
		toBB = bishop_attacks(from, pos.pieces()) & target;
		toBB.for_each([&](Square to) {
			if (is_promotion_area(us, to) || is_promotion_area(us, from))
				*moveList++ = make_move_promote(from, to);
			else
				*moveList++ = make_move(from, to);
		});
	});

	// 7. 飛車
	pos.pieces(us, ROOK).for_each([&](Square from) {
		toBB = rook_attacks(from, pos.pieces()) & target;
		toBB.for_each([&](Square to) {
			if (is_promotion_area(us, to) || is_promotion_area(us, from))
				*moveList++ = make_move_promote(from, to);
			else
				*moveList++ = make_move(from, to);
		});
	});

	// 8. 馬
	pos.pieces(us, HORSE).for_each([&](Square from) {
		toBB = horse_attacks(from, pos.pieces()) & target;
		toBB.for_each([&](Square to) {
			*moveList++ = make_move(from, to);
		});
	});

	// 9. 龍
	pos.pieces(us, DRAGON).for_each([&](Square from) {
		toBB = dragon_attacks(from, pos.pieces()) & target;
		toBB.for_each([&](Square to) {
			*moveList++ = make_move(from, to);
		});
	});

	return moveList;
}

ExtMove* generate_drops(const Position& pos, ExtMove* moveList, Color us, Bitboard target) {

	if (!target)
		return moveList;

	Hand hand = pos.hand_of(us);

	// 歩
	if (has(PAWN, hand))
	{
		Bitboard toBB = target;
		Bitboard pawnCheck = pawn_attacks(~us, pos.king_square(~us));

		// 1段目を除く
		toBB &= ~(us == BLACK ? Rank1BB : Rank9BB);

		// 二歩になる列を除く
		toBB &= ~file_fill(pos.pieces(us, PAWN));

		// 打ち歩詰めを除く
		if (pawnCheck & toBB) {
			Square to = pawnCheck.pop_c();

			if (!pos.legal_drop_pawn(to))
				toBB ^= pawnCheck;
		}

		toBB.for_each([&](Square to) {
			*moveList++ = make_move_drop(PAWN, to);
		});
	}

	// 歩以外
	if (has_except_pawn(hand))
	{
		PieceType pieceHand[6];
		int pieceNum = 0;

		if (has(KNIGHT, hand)) pieceHand[pieceNum++] = KNIGHT;
		int noKnightIdx = pieceNum;

		if (has(LANCE, hand)) pieceHand[pieceNum++] = LANCE;
		int noLanceIdx = pieceNum;

		if (has(SILVER, hand)) pieceHand[pieceNum++] = SILVER;
		if (has(GOLD, hand)) pieceHand[pieceNum++] = GOLD;
		if (has(BISHOP, hand)) pieceHand[pieceNum++] = BISHOP;
		if (has(ROOK, hand)) pieceHand[pieceNum++] = ROOK;

		Bitboard dropRank1 = target & (us == BLACK ? Rank1BB : Rank9BB);
		Bitboard dropRank2 = target & (us == BLACK ? Rank2BB : Rank8BB);
		Bitboard dropOther = target & ~(dropRank1 | dropRank2);

		// 1段目
		dropRank1.for_each([&](Square to) {
			for (int i = noLanceIdx; i < pieceNum; ++i)
				*moveList++ = make_move_drop(pieceHand[i], to);
		});

		// 2段目
		dropRank2.for_each([&](Square to) {
			for (int i = noKnightIdx; i < pieceNum; ++i)
				*moveList++ = make_move_drop(pieceHand[i], to);
		});

		// 3~9段目
		dropOther.for_each([&](Square to) {
			for (int i = 0; i < pieceNum; ++i)
				*moveList++ = make_move_drop(pieceHand[i], to);
		});
	}

	return moveList;
}

ExtMove* discovered_checks(const Position& pos, ExtMove* moveList, Color us, Bitboard target) {

	Square ksq = pos.king_square(~us);
	Bitboard dc = pos.discovered_check_candidates();
	Bitboard promotionArea = promotion_area(us);
	Bitboard rank3_9 = (us == BLACK ? Rank3_9 : Rank1_7);
	Bitboard rank4_9 = (us == BLACK ? Rank4_9 : Rank1_6);

	while (dc)
	{
		Square from = dc.pop();
		Piece pc = pos.piece_on(from);
		PieceType pt = type_of(pc);
		Bitboard attacks = pos.attacks_from(pc, from);
		Bitboard discoveredChecks = attacks & target & ~line_bb(ksq, from);

		// 不成 (香の2段目、直接の王手は除く)
		Piece opponentPiece = opponent_piece(pc);
		Bitboard directChecks = pos.attacks_from(opponentPiece, ksq);
		Bitboard nonPromotions = discoveredChecks & ~directChecks;
		if (pt == BISHOP || pt == ROOK)
		{
			if (!is_promotion_area(us, from)) {
				nonPromotions &= rank4_9;
				nonPromotions.for_each([&](Square to) {
					*moveList++ = make_move(from, to);
				});
			}
		}
		else 
		{
			if (pt == PAWN)
				nonPromotions &= rank4_9;
			else if (pt == LANCE || pt == KNIGHT)
				nonPromotions &= rank3_9;

			nonPromotions.for_each([&](Square to) {
				*moveList++ = make_move(from, to);
			});
		}	

		// 成り (直接の王手は除く)
		if (can_promote(pc)) 
		{
			Piece promotedPiece = promoted_piece(opponentPiece);
			Bitboard directChecks = pos.attacks_from(promotedPiece, ksq);
			Bitboard promotions = discoveredChecks & ~directChecks;
			if (pt == SILVER || pt == BISHOP || pt == ROOK) 
			{
				if (!is_promotion_area(us, from))
					promotions &= promotionArea;
			}
			else
			{
				assert(pt == PAWN || pt == LANCE || pt == KNIGHT);
				promotions &= promotionArea;
			}
			promotions.for_each([&](Square to) {
				*moveList++ = make_move_promote(from, to);
			});
		} 
	}

	return moveList;
}

ExtMove* direct_checks(const Position& pos, ExtMove* moveList, Color us, Bitboard target) {

	Square ksq = pos.king_square(~us);
	Bitboard ownPieces = pos.pieces(us);
	Bitboard promotionArea = promotion_area(us);
	Bitboard rank3_9 = (us == BLACK ? Rank3_9 : Rank1_7);
	Bitboard rank4_9 = (us == BLACK ? Rank4_9 : Rank1_6);

	auto minor_piece_generator = [&](PieceType pt) -> void
	{
		Piece pc = make_piece(us, pt);
		Piece opponentPiece = opponent_piece(pc);
		Piece promotedPiece = promoted_piece(opponentPiece);

		Bitboard pieces = (pt == GOLD ? pos.gold_pieces(us) : pos.pieces(us, pt));
		Bitboard candidates = checker_candidates(pc, ksq);

		if (!candidates)
			return;

		Bitboard promotionTarget = pos.attacks_from(promotedPiece, ksq) & ~ownPieces;
		Bitboard nonPromotionTarget = pos.attacks_from(opponentPiece, ksq) & ~ownPieces;

		if (pt == PAWN)
			nonPromotionTarget &= rank4_9;
		else if (pt == LANCE || pt == KNIGHT)
			nonPromotionTarget &= rank3_9;

		if (pt != SILVER)
			promotionTarget &= promotionArea;

		Bitboard checker_candidates = pieces & candidates;
		checker_candidates.for_each([&](Square from)
		{
			Bitboard attacks = pos.attacks_from(pc, from) & target;
			Bitboard nonPromotions = attacks & nonPromotionTarget;
			nonPromotions.for_each([&](Square to) {
				*moveList++ = make_move(from, to);
			});

			if (can_promote(pc))
			{
				Bitboard promotions = attacks & promotionTarget;
				if (pt == SILVER && relative_rank(us, rank_of(from)) >= RANK_4)
					promotions &= promotionArea;
				promotions.for_each([&](Square to) {
					*moveList++ = make_move_promote(from, to);
				});
			}
		});
	};

	auto major_piece_generator = [&](PieceType pt) -> void
	{
		if (pt == HORSE || pt == DRAGON) 
		{
			Piece pc = make_piece(us, pt);
			Bitboard nonPromotionTarget = pos.attacks_from(pc, ksq) & target;
			pos.pieces(us, pt).for_each([&](Square from) {
				if (attacks_bb(pc, from, ALL0BB) & nonPromotionTarget) {
					Bitboard attacks = pos.attacks_from(pc, from);
					Bitboard nonPromotions = attacks & nonPromotionTarget;
					nonPromotions.for_each([&](Square to) {
						*moveList++ = make_move(from, to);
					});
				}
			});
		}
		else
		{
			Piece pc = make_piece(us, pt);
			Bitboard nonPromotionTarget = pos.attacks_from(pc, ksq) & target;
			Bitboard promotionTarget = pos.attacks_from(promoted_piece(pc), ksq) & target;
			Bitboard pieces = pos.pieces(us, pt);
			(pieces & promotionArea).for_each([&](Square from) {
				if (attacks_bb(pc, from, ALL0BB) & promotionTarget) {
					Bitboard attacks = pos.attacks_from(pc, from);
					Bitboard promotions = attacks & promotionTarget;
					promotions.for_each([&](Square to) {
						*moveList++ = make_move_promote(from, to);
					});
				}
			});
			(pieces & ~promotionArea).for_each([&](Square from) {
				if (attacks_bb(pc, from, ALL0BB) & promotionTarget) {
					Bitboard attacks = pos.attacks_from(pc, from);
					Bitboard promotions = attacks & promotionTarget & promotionArea;
					Bitboard nonPromotions = attacks & nonPromotionTarget & ~promotionArea;
					promotions.for_each([&](Square to) {
						*moveList++ = make_move_promote(from, to);
					});
					nonPromotions.for_each([&](Square to) {
						*moveList++ = make_move(from, to);
					});
				}
			});
		}
	};

	// 小駒
	minor_piece_generator(PAWN);
	minor_piece_generator(LANCE);
	minor_piece_generator(KNIGHT);
	minor_piece_generator(SILVER);
	minor_piece_generator(GOLD); // 成駒も含む

	// 大駒
	major_piece_generator(BISHOP);
	major_piece_generator(HORSE);
	major_piece_generator(ROOK);
	major_piece_generator(DRAGON);

	return moveList;
}

ExtMove* drop_checks(const Position& pos, ExtMove* moveList, Color us, Bitboard target) {

	Hand hand = pos.hand_of(us);
	Square ksq = pos.king_square(~us);

	auto generator = [&](PieceType pt) -> void
	{
		if (!has(pt, hand))
			return;

		Piece opponentPiece = make_piece(~us, pt);
		Bitboard dropTarget = pos.attacks_from(opponentPiece, ksq) & target;

		if (!dropTarget)
			return;

		if (pt == PAWN && !pos.legal_drop_pawn(dropTarget.pop_c()))
			return;

		dropTarget.for_each([&](Square to) {
			*moveList++ = make_move_drop(pt, to);
		});
	};

	generator(PAWN);
	generator(LANCE);
	generator(KNIGHT);
	generator(SILVER);
	generator(GOLD);
	generator(BISHOP);
	generator(ROOK);

	return moveList;
}

template<Color Us, GenType Type>
ExtMove* generate_all(const Position& pos, ExtMove* moveList, Bitboard target) {

	moveList = generate_moves(pos, moveList, Us, target);
	moveList = generate_drops(pos, moveList, Us, target & ~pos.pieces());

	if (Type != EVASIONS)
	{
		Square ksq = pos.king_square(Us);
		Bitboard b = king_attacks(ksq) & target;
		b.for_each([&](Square to) {
			*moveList++ = make_move(ksq, to);
		});
	}

	return moveList;
}

} // namespace

// generate<CAPTURES> generates all pseudo-legal captures and queen
// promotions. Returns a pointer to the end of the move list.
//
// generate<QUIETS> generates all pseudo-legal non-captures.
// Returns a pointer to the end of the move list.
//
// generate<NON_EVASIONS> generates all pseudo-legal captures and
// non-captures. Returns a pointer to the end of the move list.
template<GenType Type>
ExtMove* generate(const Position& pos, ExtMove* moveList) {

	assert(Type == CAPTURES || Type == QUIETS || Type == NON_EVASIONS);
	assert(!pos.checkers());

	Color us = pos.side_to_move();

	Bitboard target = Type == CAPTURES	   ? pos.pieces(~us)			// 相手駒の位置BB
					: Type == QUIETS	   ? ~pos.pieces()				// 駒のない場所BB
					: Type == NON_EVASIONS ? ~pos.pieces(us) : ALL0BB;	// 自駒のない場所BB

	return us == BLACK	? generate_all<BLACK, Type>(pos, moveList, target)
						: generate_all<WHITE, Type>(pos, moveList, target);
}

// Explicit template instantiations
template ExtMove* generate<CAPTURES>(const Position&, ExtMove*);
template ExtMove* generate<QUIETS>(const Position&, ExtMove*);
template ExtMove* generate<NON_EVASIONS>(const Position&, ExtMove*);

// generate_checks generates all pseudo-legal checks.
// Returns a pointer to the end of the move list.
template<GenType Type>
ExtMove* generate_checks(const Position& pos, ExtMove* moveList) {

	assert(Type == CHECKS || Type == QUIET_CHECKS);
	assert(!pos.checkers());

	Color us = pos.side_to_move();
	Bitboard target = Type == CHECKS ? ~pos.pieces(us)
					: Type == QUIET_CHECKS ? ~pos.pieces() : ALL0BB;

	moveList = discovered_checks(pos, moveList, us, target);
	moveList = direct_checks(pos, moveList, us, target);
	moveList = drop_checks(pos, moveList, us, ~pos.pieces());

	return moveList;
}

// Explicit template instantiations
template ExtMove* generate_checks<CHECKS>(const Position&, ExtMove*);
template ExtMove* generate_checks<QUIET_CHECKS>(const Position&, ExtMove*);

// generate<EVASIONS> generates all pseudo-legal check evasions when the side
// to move is in check. Returns a pointer to the end of the move list.
template<>
ExtMove* generate<EVASIONS>(const Position& pos, ExtMove* moveList) {

	assert(pos.checkers());

	Color us = pos.side_to_move();
	Square ksq = pos.king_square(us);
	Bitboard sliderAttacks = ALL0BB;
	Bitboard checkers = pos.checkers();
	Bitboard occupied = pos.pieces() ^ Bitboard(ksq);

	// Find all the squares attacked by slider checkers. We will remove them from
	// the king evasions in order to skip known illegal moves, which avoids any
	// useless legality checks later on.
	do {
		Square checksq = checkers.pop();
		sliderAttacks |= attacks_bb(pos.piece_on(checksq), checksq, occupied);
	} while (checkers);

	// Generate evasions for king, capture and non capture moves
	Bitboard b = king_attacks(ksq) & ~pos.pieces(us) & ~sliderAttacks;
	while (b)
		*moveList++ = make_move(ksq, b.pop());

	// Double check, only a king move can save the day
	if (more_than_one(pos.checkers()))
		return moveList;

	// Generate blocking evasions or captures of the checking piece
	Square checksq = pos.checkers().pop_c();
	Bitboard target = between_bb(checksq, ksq) | checksq;

	return us == BLACK ? generate_all<BLACK, EVASIONS>(pos, moveList, target)
					   : generate_all<WHITE, EVASIONS>(pos, moveList, target);
}

// generate<LEGAL> generates all the legal moves in the given position
template<>
ExtMove* generate<LEGAL>(const Position& pos, ExtMove* moveList) {

	Bitboard pinned = pos.pinned_pieces(pos.side_to_move());
	Square ksq = pos.king_square(pos.side_to_move());
	ExtMove* cur = moveList;

	moveList = pos.checkers() ? generate<EVASIONS>(pos, moveList)
							  : generate<NON_EVASIONS>(pos, moveList);

	while (cur != moveList)
		if ((pinned || from_sq(*cur) == ksq) && !pos.legal(*cur))
			*cur = (--moveList)->move;
		else
			++cur;

	return moveList;
}