#include "mate1ply.h"

namespace {

inline Bitboard new_pinned(const Position& pos, Color c, Square exclusion) {

	Square ksq = pos.king_square(c);
	Bitboard excl = ~Bitboard(exclusion);
	Bitboard result = ALL0BB;

	// Squareにある駒をsnipers(pinners)とblockersに含まないようにする
	// これで王手駒の移動によるpinnersとblockersの変化に対応する
	Bitboard snipers = ((pos.pieces(~c, ROOK, DRAGON) & rook_attacks(ksq, ALL0BB))
		| (pos.pieces(~c, BISHOP, HORSE) & bishop_attacks(ksq, ALL0BB))
		| (pos.pieces(~c, LANCE) & lance_attacks(c, ksq, ALL0BB))) & excl;
	Bitboard pieces = pos.pieces() & excl;

	while (snipers) 
	{
		Square sniperSq = snipers.pop();
		Bitboard blocker = between_bb(ksq, sniperSq) & pieces;
		if (!more_than_one(blocker))
			result |= blocker;
	}

	return result;
}

inline bool can_evade_king(const Position& pos, Color us, Piece checkPc, Square from, Square to, Bitboard occupied) {

	Square ksq = pos.king_square(us);

	occupied |= Bitboard(to);  // fromからtoへ王手駒が移動したことにする
	occupied ^= Bitboard(ksq); // 玉がないと考えて利きを考える必要がある

	Bitboard attacks = attacks_bb(checkPc, to, occupied);
	Bitboard evasions = king_attacks(ksq) & ~(attacks | to | pos.pieces(us));

	// 駒打ちのとき
	if (from == SQ_NONE) {
		while (evasions) {
			Square escape = evasions.pop();
			if (!(pos.attackers_to(escape, occupied) & pos.pieces(~us)))
				return true;
		}
	} 
	// 駒移動のとき
	else {
		while (evasions) {
			Square escape = evasions.pop();
			// ~Bitboard(from)がないと、attackers_toでfromにあるcheckPcの利きが残ってしまう
			if (!(pos.attackers_to(escape, occupied) & pos.pieces(~us) & ~Bitboard(from)))
				return true;
		}
	}

	return false;
}

inline bool can_catch_checker(const Position& pos, Color us, Square to, const Bitboard& pinned, Bitboard defenders) {

	Square ksq = pos.king_square(us);
	defenders &= ~Bitboard(ksq);

	while (defenders) {
		Square catcher = defenders.pop();
		if (!(pinned & catcher) || aligned(catcher, to, ksq))
			return true;
	}

	return false;
}

template<Color Us, PieceType Pt>
inline bool is_checkmate_drop(const Position& pos, Square to) {

	Color them = ~Us;

	// 玉の位置
	Square ksq = pos.king_square(them);

	// 盤面の状態
	Bitboard occupied = pos.pieces();

	// pinされた駒
	Bitboard pinned = pos.pinned_pieces(them);

	// 詰ます側と守る側
	Bitboard target = pos.attackers_to(to, occupied);
	Bitboard attackers = target & pos.pieces(Us);
	Bitboard defenders = target & pos.pieces(them);

	// 1. 攻撃駒以外の味方駒がtoに利いていなければ詰まない
	if (Pt != KNIGHT)
		if (!attackers)
			return false;

	// 2. 玉以外の防御駒が攻撃駒を取れれば詰まない
	if (can_catch_checker(pos, them, to, pinned, defenders)) 
		return false;

	// 3. 玉が逃げられれば詰まない
	if (can_evade_king(pos, them, make_piece(Us, Pt), SQ_NONE, to, occupied))
		return false;

	return true;
}

template<Color Us, PieceType Pt>
inline bool find_mate_drop(const Position& pos, Bitboard target, Move& mateMove) {

	Hand hand = pos.hand_of(Us);

	if (!has(Pt, hand))
		return false;

	if (Pt != KNIGHT && !target)
		return false;

	if (Pt == LANCE && has(ROOK, hand))
		return false;

	if (Pt == SILVER && has(BISHOP, hand) && has(GOLD, hand))
		return false;

	Color them = ~Us;
	Piece pc = make_piece(Us, Pt);
	Piece opponentPiece = opponent_piece(pc);
	Square ksq = pos.king_square(them);
	Bitboard occupied = pos.pieces() & ~Bitboard(ksq);
	Bitboard dropTarget = Pt == KNIGHT ? knight_attacks(them, ksq) & ~occupied
									   : attacks_bb(opponentPiece, ksq, ALL1BB) & target;

	// 相手玉の尻に打つなら金 < 飛
	if (Pt == GOLD && has(ROOK, hand))
		dropTarget &= ~pawn_attacks(Us, ksq);

	else if (Pt == SILVER)
	{
		// 相手玉の頭3方向に打つなら銀 < 金
		if (has(GOLD, hand))
			dropTarget &= ~(silver_attacks(them, ksq) & gold_attacks(them, ksq));

		// 相手玉の尻斜め2方向に打つなら銀 < 角
		else if (has(BISHOP, hand))
			dropTarget &= ~(silver_attacks(Us, ksq) & gold_attacks(Us, ksq));
	}

	while (dropTarget) {
		Square to = dropTarget.pop();
		if (is_checkmate_drop<Us, Pt>(pos, to)) {
			mateMove = make_move_drop(Pt, to);
			assert(pos.legal(mateMove) && pos.pseudo_legal(mateMove));
			return true;
		}
	}

	return false;
}

template<Color Us, PieceType Pt, bool IsPromotion>
inline bool is_checkmate_move(const Position& pos, Square from, Square to) {

	// 王手する駒
	Piece checkPc = IsPromotion ? promoted_piece(make_piece(Us, Pt))
								: make_piece(Us, Pt);
	Color them = ~Us;

	// 玉の位置
	Square ourKsq = pos.king_square(Us);
	Square oppKsq = pos.king_square(them);

	// 盤面の状態
	Bitboard occupied = pos.pieces() ^ from;
	Bitboard dcCandidates = pos.discovered_check_candidates();

	// pinされた駒
	Bitboard ourPinned = pos.pinned_pieces(Us);
	Bitboard oppPinned = pos.pinned_pieces(them);
	Bitboard newPinned = new_pinned(pos, them, from);

	// 詰ます側と守る側
	Bitboard target = pos.attackers_to(to, occupied);
	Bitboard attackers = target & pos.pieces(Us);
	Bitboard defenders = target & pos.pieces(them);

	// 1. 攻撃駒以外の味方駒がtoに利いていなければ詰まない
	// 2. 玉以外の防御駒が攻撃駒を取れれば詰まない(両王手のときは飛ばす)
	// 3. 玉が逃げられれば詰まない
	// 4. 詰ます手が自殺手であるなら指せない

	switch (Pt) {
	case PAWN: // 1, 2(oppPinned), 3, 4
		if (!(attackers ^ from)) return false;
		if (can_catch_checker(pos, them, to, oppPinned, defenders)) return false;
		if (can_evade_king(pos, them, checkPc, from, to, occupied)) return false;
		if ((ourPinned & from) && !aligned(from, to, ourKsq)) return false;
		return true;

	case LANCE: // 1, 2(dc, oppPinned), 3, 4
		if (!(attackers ^ from)) return false;
		if (!(dcCandidates & from))
			if (can_catch_checker(pos, them, to, oppPinned, defenders)) return false;
		if (can_evade_king(pos, them, checkPc, from, to, occupied)) return false;
		if ((ourPinned & from) && !aligned(from, to, ourKsq)) return false;
		return true;

	case KNIGHT: // 1(promotion only), 2(dc, newPinned), 3, 4
		if (IsPromotion)
			if (!(attackers ^ from)) return false;
		if (!(dcCandidates & from))
			if (can_catch_checker(pos, them, to, newPinned, defenders)) return false;
		if (can_evade_king(pos, them, checkPc, from, to, occupied)) return false;
		if ((ourPinned & from) && !aligned(from, to, ourKsq)) return false;
		return true;

	case SILVER: case GOLD: case HORSE: // 1, 2(dcとalign, newPinned), 3, 4
		if (!(attackers ^ from)) return false;
		if (!(dcCandidates & from) || aligned(from, to, oppKsq))
			if (can_catch_checker(pos, them, to, newPinned, defenders)) return false;
		if (can_evade_king(pos, them, checkPc, from, to, occupied)) return false;
		if ((ourPinned & from) && !aligned(from, to, ourKsq)) return false;
		return true;

	case BISHOP: case ROOK: // 1, 2(dc, newPinned), 3, 4
		if (!(attackers ^ from)) return false;
		if (!(dcCandidates & from))
			if (can_catch_checker(pos, them, to, newPinned, defenders)) return false;
		if (can_evade_king(pos, them, checkPc, from, to, occupied)) return false;
		if ((ourPinned & from) && !aligned(from, to, ourKsq)) return false;
		return true;
		
	case DRAGON: // 1, 2(newPinned), 3, 4
		if (!(attackers ^ from)) return false;
		if (can_catch_checker(pos, them, to, newPinned, defenders)) return false;
		if (can_evade_king(pos, them, checkPc, from, to, occupied)) return false;
		if ((ourPinned & from) && !aligned(from, to, ourKsq)) return false;
		return true;
	}

	assert(false);
	return false;
}

template<Color Us, PieceType Pt>
inline bool find_mate_move(const Position& pos, Bitboard target, Move& mateMove) {

	Piece pc = make_piece(Us, Pt);
	Piece opponentPiece = opponent_piece(pc);
	Piece promotedPiece = promoted_piece(opponentPiece);
	Square ksq = pos.king_square(~Us);

	Bitboard pieces = (Pt == GOLD ? pos.gold_pieces(Us) : pos.pieces(Us, Pt));
	Bitboard candidates = adjacent_check_candidates(pc, ksq);

	if (!candidates || !target)
		return false;

	Bitboard promotionTarget = attacks_bb(promotedPiece, ksq, ALL1BB) & target;
	Bitboard nonPromotionTarget = attacks_bb(opponentPiece, ksq, ALL1BB);
	nonPromotionTarget &= (Pt == KNIGHT ? ~pos.pieces(Us) : target); // 桂成らずの王手は8近傍でない

	Bitboard checkerCandidates = pieces & candidates;
	while (checkerCandidates)
	{
		Square from = checkerCandidates.pop();
		Bitboard attacks = pos.attacks_from(pc, from);

		// 成り
		if (can_promote(pc))
		{
			Bitboard promotionChecks = attacks & promotionTarget;
			if ((Pt == PAWN || Pt == LANCE || Pt == KNIGHT) || !is_promotion_area(Us, from))
				promotionChecks &= promotion_area(Us);

			while (promotionChecks) {
				Square to = promotionChecks.pop();
				if (is_checkmate_move<Us, Pt, true>(pos, from, to)) {
					mateMove = make_move_promote(from, to);
					assert(pos.legal(mateMove) && pos.pseudo_legal(mateMove));
					return true;
				}
			}
		}

		// 成らず
		Bitboard nonPromotionChecks = attacks & nonPromotionTarget;
		if (Pt == PAWN || Pt == BISHOP || Pt == ROOK)
			nonPromotionChecks &= (Us == BLACK ? Rank4_9 : Rank1_6);
		else if (Pt == LANCE || Pt == KNIGHT)
			nonPromotionChecks &= (Us == BLACK ? Rank3_9 : Rank1_7);

		while (nonPromotionChecks) {
			Square to = nonPromotionChecks.pop();
			if (is_checkmate_move<Us, Pt, false>(pos, from, to)) {
				mateMove = make_move(from, to);
				assert(pos.legal(mateMove) && pos.pseudo_legal(mateMove));
				return true;
			}
		}
	}

	return false;
}

template<Color Us>
bool is_mate_in_1ply(const Position& pos, Move& mateMove) {

	assert(!pos.checkers());

	Square ksq = pos.king_square(~Us);
	Bitboard target = king_attacks(ksq) & ~pos.pieces(Us);

	// 駒の移動手
	if (find_mate_move<Us, PAWN>(pos, target, mateMove)) return true;
	if (find_mate_move<Us, LANCE>(pos, target, mateMove)) return true;
	if (find_mate_move<Us, KNIGHT>(pos, target, mateMove)) return true;
	if (find_mate_move<Us, SILVER>(pos, target, mateMove)) return true;
	if (find_mate_move<Us, GOLD>(pos, target, mateMove)) return true;
	if (find_mate_move<Us, BISHOP>(pos, target, mateMove)) return true;
	if (find_mate_move<Us, ROOK>(pos, target, mateMove)) return true;
	if (find_mate_move<Us, HORSE>(pos, target, mateMove)) return true;
	if (find_mate_move<Us, DRAGON>(pos, target, mateMove)) return true;

	// 駒打ち
	if (has_except_pawn(pos.hand_of(Us)))
	{
		target &= ~pos.pieces(~Us);

		if (find_mate_drop<Us, LANCE>(pos, target, mateMove)) return true;
		if (find_mate_drop<Us, KNIGHT>(pos, target, mateMove)) return true;
		if (find_mate_drop<Us, SILVER>(pos, target, mateMove)) return true;
		if (find_mate_drop<Us, GOLD>(pos, target, mateMove)) return true;
		if (find_mate_drop<Us, BISHOP>(pos, target, mateMove)) return true;
		if (find_mate_drop<Us, ROOK>(pos, target, mateMove)) return true;
	}

	return false;
}

} // namespace

bool is_mate_in_1ply(const Position& pos, Move& mateMove)
{
	Color us = pos.side_to_move();
	return us == BLACK ? is_mate_in_1ply<BLACK>(pos, mateMove)
					   : is_mate_in_1ply<WHITE>(pos, mateMove);
}