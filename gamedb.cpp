#include <iomanip>
#include <sstream>
#include <unordered_set>

#include "gamedb.h"
#include "notations.h"
#include "thread.h"
#include "usi.h"

namespace {

const std::unordered_set<std::string> TitleMatches = {
	"—³‰¤í", "–¼lí", "‰¤ˆÊí", "‰¤Àí", "Šû‰¤í", "‰¤«í", "Šû¹í", "‡ˆÊí"
};

} // namespace

bool GameDatabase::read_one_game(Game* game)
{
	assert(game != nullptr);

START:
	std::string line;

	// Step 1. reads the header
	if (!std::getline(inputStream, line))
		return false;

	// Step 2. analyze this header
	std::istringstream headerInput(line);
	int length, id, result;
	headerInput >> id >> game->date
		>> game->players[BLACK] >> game->players[WHITE]
		>> result >> length >> game->event >> game->opening;
	game->result = static_cast<Game::Result>(result);

	// Step 3. read the moves
	if (!std::getline(inputStream, line))
		return false;

	// if specified, skipping than bout
	if (titleMatchesOnly && TitleMatches.count(game->event) == 0)
		goto START;

	// Step 4. parse the moves
	std::istringstream movesInput(line);
	StateInfo st;
	Position pos;
	pos.set(StartSFEN, &st, Threads.main());

	std::vector<StateInfo> stEach(length);
	game->moves.clear();

	for (int ply = 0; ply < length; ++ply)
	{
		std::string moveStr;
		movesInput >> std::setw(6) >> moveStr;

		if (!movesInput)
			break;

		Move move = Csa::parse_move(pos, moveStr);

		if (!pos.legal(move) || !pos.pseudo_legal(move))
			break;

		game->moves.push_back(move);
		pos.do_move(move, stEach[ply], pos.gives_check(move));
	}

	return true;
}