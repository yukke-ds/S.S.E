#ifndef NOTATIONS_H_INCLUDED
#define NOTATIONS_H_INCLUDED

#include "position.h"

class Csa {
public:
	static Move parse_move(Position& pos, std::string& csa);
};

#endif // ifndef NOTATIONS_H_INCLUDED