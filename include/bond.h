#ifndef BOND_H
#define BOND_H

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
class bond
{
private:
	//std::string bondInfor; //one bond line
public:
	char aACS;			   //atom A chemical symbol
	char aACSE;			   //atom A chemical symbol extention
	char aBCS;			   //atom B chemical symbol
	char aBCSE;			   //atom B chemical symbol extention
	int bType;             //bond type, 1 single, 2 double, 3 triple bond, 4 aromic bond, 5 aromic triple bond
	int bA;                //bonded atom order A
	int bB;                //bonded atom order B
	int stereo;            //0 = not stereo,1 = Up, 4 = Either,6 = Down, Double bonds: 0 = Use x-, y-, z-coords from atom block to determine cis or trans,3 = Cis or trans (either) double bond
	bool used;
	bool bAused;
	bool bBused;
	bond(std::string oneLine);
	bond(int bt, int ba, int bb);
	bond();
	void setBond(std::string& bondInfor);
	void standardize();		//make sure bA less than bB
	bool operator == (const bond& b) const;
	void showBond();
	std::string toString();
};

#endif
