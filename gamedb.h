#ifndef GAMEDB_H_INCLUDED
#define GAMEDB_H_INCLUDED

#include <istream>
#include <string>
#include <vector>
#include "types.h"

struct Game {

	enum Result {
		BLACK_WIN = 1,
		WHITE_WIN,
		DRAW
	};

	std::string players[COLOR_NB];			// �΋ǎ�
	Result result;							// �΋ǌ���
	std::string date;						// �΋Ǔ�
	std::string event;						// ����̖��O
	std::string opening;					// ��^�̖��O
	std::vector<Move> moves;				// �����̖��O
};

class GameDatabase {
public:
	static constexpr const char* defaultDatabaseFile = "TrainingData/kifu_nonComputers.txt";
	static constexpr const char* defaultTrainingFile = "TrainingData/kifu/for_training.sfen";
	static constexpr const char* defaultTestFile = "TrainingData/kifu/for_test.sfen";
	static constexpr const char* outputSFENsFile = "TrainingData/SFENs.txt";
	static constexpr const char* outputParentSFENsFile = "TrainingData/parentSFENs.txt";
	static constexpr const char* outputObservedSFENsFile = "TrainingData/observedSFENs.txt";
	static constexpr const char* outputRandomSFENsFile = "TrainingData/randomSFENs.txt";

	GameDatabase(std::istream& is) : inputStream(is) {}

	bool read_one_game(Game* game);
	void set_title_matches_only(bool titleMatchesOnly);

private:
	std::istream& inputStream;
	bool titleMatchesOnly = false;
};

inline void GameDatabase::set_title_matches_only(bool b) {
	titleMatchesOnly = b;
}

#endif