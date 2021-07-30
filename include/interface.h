
#ifndef INTERFACE
#define INTERFACE
#include "molecule.h"
#include "CalcSim.h"
#include "StrucMatch.h"
#include "SubSMatch.h"
#include "socketConn.h"
#include <map>

extern std::string molEnd;
extern std::string SDFdelim;

//get hit molecules from molPool, by send in one string formate mol, return 1 (hit) or 0 (no hit)
	//matchType: 0 similarity match, 1 structure match, 2 sub-structure match
	//st: similarity threshold
	//molType: 0 MDL mol formate, 1 ctable
int molMatch(std::string mol, std::vector<std::string>& molIDs,
		int matchType = 0, float st = 0.9, int molType = 0);
int molMatch(std::string& mol, std::vector<molecule>& molPool, std::vector<simResult>& resultSet,
		int matchType = 0, int min_sim = 90, int molType = 0);

bool molMatch(std::string mol, std::vector<molecule> &oneMolfromPool, int &resultValue,
	   	int matchType = 0, int min_sim = 90, int molType = 0);

//get hit molecules from molPool, by send a file, return 1 (hit) or 0 (no hit)
	//matchType: 0 similarity match, 1 structure match, 2 sub-structure match
	//st: similarity threshold
	//molType: 0 MDL mol formate, 1 ctable
int molsMatch(std::string fileName, std::vector<std::string>& molIDs,
			  int matchType = 0, float st = 0.9, int molType = 0);

//get unique molecule set from two molecule set
bool mergeMols(std::vector<molecule>& vA, std::vector<molecule>& vB, std::vector<molecule>& retv);

//load molecule from file, create molPool, never check duplication
int loadMol(std::string fileName);

//load molecule from file, create molPool
//	it put molecule name at map->first
//	if molecule name meet duplecation, it push_back this molcule in map->second and remove duplicate molecule
int loadMol(std::string fileName, std::map<std::string, std::vector<molecule> >  &molPool);

//load molecule from file, create molPool
//	created molecule pool contains unique eliment(molecule)
//	merge duplicated molecule's molName
int loadMol(std::string fileName, std::vector<molecule> &molPool);


//make a molecule pool server by socket
int makeSocketMolServer(int lintenPort, std::vector<molecule> &molPoll, std::vector<simResult> &rs);
#endif
