#include <iostream>
#include <sstream>
#include <string>

#include "usi.h"

#include "mate1ply.h"
#include "ann/ann_evaluator.h";
#include "learn.h"
#include "movegen.h"
#include "position.h"
#include "progress.h"
#include "search.h"
#include "thread.h"
#include "timeman.h"

using namespace std;

// SFEN string of the initial position, hirate shogi
const char* StartSFEN = "lnsgkgsnl/1r5b1/ppppppppp/9/9/9/PPPPPPPPP/1B5R1/LNSGKGSNL b - 1";

namespace {

void isready() {

	static bool first = true;

	if (first)
	{
		Eval::load_eval();

		first = false;
	}

	sync_cout << "readyok" << sync_endl;
}

// position() is called when engine receives the "position" UCI command.
// The function sets up the position described in the given FEN string ("fen")
// or the starting position ("startpos") and then makes the moves given in the
// following move list ("moves").
void position(Position& pos, istringstream& is, StateListPtr& states) {

	Move m;
	string token, sfen;

	is >> token;

	if (token == "startpos")
	{
		sfen = StartSFEN;
		is >> token; // Consume "moves" token if any
	}
	else if (token == "sfen")
		while (is >> token && token != "moves")
			sfen += token + " ";
	else
		return;

	states = StateListPtr(new std::deque<StateInfo>(1));
	pos.set(sfen, &states->back(), Threads.main());

	// Parse move list (if any)
	while (is >> token && (m = USI::to_move(token)) != MOVE_NONE)
	{
		states->emplace_back();
		pos.do_move(m, states->back());
	}
}

void setoption(istringstream& is) {

	string token, name, value;

	is >> token; // Consume "name" token
					 
	// Read option name (can contain spaces)
	while (is >> token && token != "value")
		name += string(" ", name.empty() ? 0 : 1) + token;

	// Read option value (can contain spaces)
	while (is >> token)
		value += string(" ", value.empty() ? 0 : 1) + token;

	if (Options.count(name))
		Options[name] = value;
	else
		sync_cout << "No such option: " << name << sync_endl;
}

// go() is called when engine receives the "go" UCI command. The function sets
// the thinking time and other parameters from the input string, then starts
// the search.
void go(Position& pos, istringstream& is, StateListPtr& states) {

	Search::LimitsType limits;
	string token;
	bool ponderMode = false;

	Time.start(); // As early as possible!

	while (is >> token)
		if (token == "searchmoves")
			while (is >> token)
				limits.searchmoves.push_back(USI::to_move(token));

		else if (token == "btime")     is >> limits.time[BLACK];
		else if (token == "wtime")     is >> limits.time[WHITE];
		else if (token == "winc")      is >> limits.inc[WHITE];
		else if (token == "binc")      is >> limits.inc[BLACK];
		else if (token == "byoyomi")   is >> limits.byoyomi;
		else if (token == "movestogo") is >> limits.movestogo;
		else if (token == "depth")     is >> limits.depth;
		else if (token == "nodes")     is >> limits.nodes;
		else if (token == "movetime")  is >> limits.movetime;
		else if (token == "mate")      is >> limits.mate;
		else if (token == "infinite")  limits.infinite = 1;
		else if (token == "ponder")    ponderMode = true;

		Threads.start_thinking(pos, states, limits, ponderMode);
}

} // namespace

// USI::loop() waits for a command from stdin, parses it and calls the appropriate
// function. Also intercepts EOF from stdin to ensure gracefully exiting if the
// GUI dies unexpectedly. When called with some command line arguments, e.g. to
// run 'bench', once the command is executed the function returns immediately.
// In addition to the USI ones, also some additional debug commands are supported.
void USI::loop(int argc, char* argv[]) {

	Position pos;
	string token, cmd;
	StateListPtr states(new std::deque<StateInfo>(1));
	auto uiThread = std::make_shared<Thread>(0);

	pos.set(StartSFEN, &states->back(), uiThread.get());

	for (int i = 1; i < argc; ++i)
		cmd += std::string(argv[i]) + " ";
	
	do {
		if (argc == 1 && !getline(cin, cmd)) // Block here waiting for input or EOF
			cmd = "quit";

		istringstream is(cmd);

		token.clear(); // Avoid a stale if getline() returns empty or blank line
		is >> skipws >> token;

		// The GUI sends 'ponderhit' to tell us the user has played the expected move.
		// So 'ponderhit' will be sent if we were told to ponder on the same move the
		// user has played. We should continue searching but switch from pondering to
		// normal search. In case Threads.stopOnPonderhit is set we are waiting for
		// 'ponderhit' to stop the search, for instance if max search depth is reached.
		if (	token == "quit"
			||  token == "stop"
			|| (token == "ponderhit" && Threads.stopOnPonderhit)
			||	token == "gameover")
			Threads.stop = true;

		else if (token == "ponderhit")
		{
			Time.start_ponderhit();
			Threads.ponder = false; // Switch to normal search
		}
		else if (token == "usi")
			sync_cout << "id name " << engine_info()
					  << "\n"		<< Options
					  << "\nusiok"	<< sync_endl;

		else if (token == "setoption")	setoption(is);
		else if (token == "go")			go(pos, is, states);
		else if (token == "position")	position(pos, is, states);
		else if (token == "usinewgame") Search::clear();
		else if (token == "isready")	isready();

		// Additional custom non-USI commands
		// for debugging
		else if (token == "flip");
		else if (token == "bench");
		else if (token == "d");
		else if (token == "eval")	sync_cout << "Evaluation: " << Eval::all_calculate(pos) << sync_endl;
		else if (token == "perft");

		// for making SFENs
		else if (token == "make-sfens")		Learn::make_sfens();
		else if (token == "make-sfens-tlp") Learn::make_triplet_sfens();

		// for learning
		else if (token == "learn-progress")	Progress::learn_parameters();
		else if (token == "learn-tdl")		Learn::TD_leaf();
		else if (token == "learn-sl")		Learn::supervised();
		else
			sync_cout << "Unknown command: " << cmd << sync_endl;

	} while (token != "quit" && argc == 1); // Command line args are one-shot
}

// USI::value() converts a Value to a string suitable for use with the USI
// protocol specification:
//
// cp <x>    The score from the engine's point of view in centipawns.
// mate <y>  Mate in y moves, not plies. If the engine is getting mated
//           use negative values for y.
string USI::value(Value v) {

	assert(-VALUE_INFINITE < v && v < VALUE_INFINITE);

	stringstream ss;

	if (abs(v) < VALUE_MATE_IN_MAX_PLY)
		ss << "cp " << v * 100 / PawnValue;
	else
		ss << "mate " << (v > 0 ? VALUE_MATE - v : -VALUE_MATE - v);

	return ss.str();
}

string USI::square(Square s) {
	return string{ char('1' + file_of(s)), char('a' + rank_of(s)) };
}

// USI::move() converts a Move to a string in coordinate notation (7g7f, S*5b).
string USI::move(Move m) {

	Square from = from_sq(m);
	Square to = to_sq(m);
	string move;

	if (m == MOVE_NONE)
		return "resign";

	if (m == MOVE_NULL)
		return "0000";

	if (type_of(m) == DROP) {
		move = " PLNSBRG"[drop_type(m)];
		move += '*';
		move += USI::square(to);
		return move;
	}

	move = USI::square(from) + USI::square(to);

	if (type_of(m) == PROMOTION)
		move += "+";

	return move;
}

// USI::to_move() converts a string representing a move in coordinate notation
// (7g7f, S*5b) to the corresponding legal Move, if any.
Move USI::to_move(string& str) {

	//for (const auto& m : MoveList<LEGAL>(pos))
	//	if (str == USI::move(m))
	//		return m;

	// Žw‚µŽè¶¬‚Å‚Í¬‚ç‚¸‚ð¶¬‚µ‚È‚¢‚½‚ß•Ê‚ÉŽÀ‘•‚µ‚Ä‚¨‚­
	if (str == "resign")
		return MOVE_NONE;

	File toFile = to_file(str[2]);
	Rank toRank = to_rank(str[3]);

	assert((toFile >= FILE_1 && toFile <= FILE_9)
		&& (toRank >= RANK_1 && toRank <= RANK_9));

	Square to = toFile | toRank;
	Move move = MOVE_NONE;
	bool promotion = str[4] == '+';
	bool drop = str[1] == '*';

	if (!drop)
	{
		File fromFile = to_file(str[0]);
		Rank fromRank = to_rank(str[1]);
		assert((fromFile >= FILE_1 && fromFile <= FILE_9) && (fromRank >= RANK_1 && fromRank <= RANK_9));
		Square from = fromFile | fromRank;
		move = promotion ? make_move_promote(from, to) : make_move(from, to);
	}
	else {
		for (int i = 1; i <= 7; ++i)
			if (" PLNSBRG"[i] == str[0]) {
				move = make_move_drop((PieceType)i, to);
				break;
			}
	}

	return move;
}