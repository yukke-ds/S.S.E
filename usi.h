#ifndef USI_H_INCLUDED
#define USI_H_INCLUDED

#include <map>
#include <string>

#include "types.h"

class Position;

namespace USI {

class Option;

// Custom comparator because USI options should be case insensitive
struct CaseInsensitiveLess {
	bool operator() (const std::string&, const std::string&) const;
};

// Our options container is actually a std::map
typedef std::map<std::string, Option, CaseInsensitiveLess> OptionsMap;

class Option {

	typedef void(*OnChange)(const Option&);

public:
	Option(OnChange = nullptr);
	Option(bool v, OnChange = nullptr);
	Option(const char* v, OnChange = nullptr);
	Option(int v, int min, int max, OnChange = nullptr);

	Option& operator=(const std::string&);
	void operator<<(const Option&);
	operator int() const;
	operator std::string() const;

private:
	friend std::ostream& operator<<(std::ostream&, const OptionsMap&);

	std::string defaultValue, currentValue, type;
	int min, max;
	size_t idx;
	OnChange on_change;
};

void init(OptionsMap&);
void loop(int argc, char* argv[]);
std::string value(Value v);
std::string square(Square s);
std::string move(Move m);
std::string pv(const Position& pos, Depth depth, Value alpha, Value beta);
Move to_move(std::string& str);

} // namespace USI

extern USI::OptionsMap Options;
extern const char* StartSFEN;

#endif // ifndef USI_H_INCLUDED