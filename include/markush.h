#ifndef MARKUSH
#define MARKUSH

#include "atom.h"
#include "bond.h"
#include <iostream>

extern std::string mdlRGroupMark;
extern std::string CtabRLOGMark;
extern std::string CRLF;

class rGroupMDL
{
private:
	std::string mdlRgroup;	//mdl RGroup line
public:
	int firstNO;	//first number in mdl rgroup line
	int serNum;		//atom serial number in molecule, second number in mdl rgroup line
	int groupNO;	//group number, third number in mdl rgroup line
	rGroupMDL();
	rGroupMDL(std::string oneLine);
	void setRGroup();
	void showRGP();
	std::string toString();
};

class rGroupCtab
{
private:
	std::string ctabLine;	//accelrys ctable one R LOG rgroup line
public:
	int firstNO;	//first number
	int secondNO;	//second number

	rGroupCtab(std::string& oneLine);
	void setRGroup();
	void showRGP();
};

class markush
{
public:
	atom pointAtom;		//R atom
	int groupNO;		//rgroup NO. this atom belong to
	std::vector<atom> atoms;	//atoms around pointAtom
	std::vector<bond> bonds;	//bonds bonded to pointAtom
	markush();
	markush(atom &pA, int &gNO, std::vector<atom>& as);
	bool operator == (const markush& m) const;
	bool getBond(int &serNum, bond& retBond);	//get bond around pointAtom by serNum
	void getBonds(std::vector<bond>& retBonds);	//get bonds around pointAtom
	void getAtoms(std::vector<atom>& retAtoms);	//get this->atoms;
	bool chkASN(int serNum);	//if serNum exist in pointAtom and arounded atoms, true
	void eraseAB();				//erase each atom's bond which bond to pointAtom in atoms
	void show();
};

#endif