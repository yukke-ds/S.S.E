#include <map>

#include "notations.h"

namespace {

const std::map<std::string, PieceType> PieceTypeFromCsa = {
	{ "FU", PAWN },
	{ "KY", LANCE },
	{ "KE", KNIGHT },
	{ "GI", SILVER },
	{ "KI", GOLD },
	{ "KA", BISHOP },
	{ "HI", ROOK },
	{ "OU", KING },
	{ "TO", PRO_PAWN },
	{ "NY", PRO_LANCE },
	{ "NK", PRO_KNIGHT },
	{ "NG", PRO_SILVER },
	{ "UM", HORSE },
	{ "RY", DRAGON },
};

Square square_from_csa(const std::string& csa) {
	assert(csa.size() == 2);
	File file = static_cast<File>(std::stoi(csa.substr(0, 1)) - 1 + FILE_1);
	Rank rank = static_cast<Rank>(std::stoi(csa.substr(1, 1)) - 1 + RANK_1);
	assert(FILE_1 <= file && file <= FILE_9);
	assert(RANK_1 <= rank && rank <= RANK_9);

	return file | rank;
}

} // namespace

Move Csa::parse_move(Position& pos, std::string& csa)
{
	assert(csa.size() == 6 || csa.size() == 7);

	// 1. divide the csa string to each part
	size_t offset = csa.size() - 6;
	std::string fromStr = csa.substr(offset + 0, 2);
	std::string toStr = csa.substr(offset + 2, 2);
	std::string pieceStr = csa.substr(offset + 4, 2);

	// 2. convert the string to 'Square' and 'PieceType'
	bool isDrop = (fromStr == "00");
	Square to = square_from_csa(toStr);
	PieceType pt = PieceTypeFromCsa.at(pieceStr);

	// 3. make move(drop, promotion or normal)
	if (isDrop)
		return make_move_drop(pt, to);

	else {
		Square from = square_from_csa(fromStr);
		PieceType ptFrom = type_of(pos.piece_on(from));
		bool isPromotion = ptFrom != pt;

		if (isPromotion)
			return make_move_promote(from, to);
		else
			return make_move(from, to);
	}
}