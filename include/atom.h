#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "bond.h"

class atom
{
public:
	//std::string atomInfor;	//one atom line
	char chemSym;			//chemical symbol
	char chemSymEx;			//chemical symbol extension
	int atomNum;		//atomic number
	float x;
	float y;
	float z;
	int massD;				//mass difference
	int charge;				//0 if uncharged; 1, 2, 3 if positive charges where 1 = +3, 2 = +2, 3 = +1; 4 if doublet radical; 5, 6, 7 if negative charges where 5 = -1, 6 = -2, 7 = -3
	int stereo;				//atom stereo parity
	bool star;				//if *, true
	bool atomH;				//if this atom is H, equal true
	bool nonC;				//if NOT sp3 C, true
	bool hAtom;				//if heteroatom, true
	bool aroma;				//if aromatic, true
	bool tAtom;				//if terminal atom, true, border atom, bCount = bCount2 = 1
	bool RGPAtom;			//if rgroup point atom, true
	int rgpNO;				//rgroup NO.
	bool used;				//if this atom was used, true
	bool newAdded;			//if this atom was new added, true;
	int serNum;				//atom serial number in molecule
	int bCount;				//count of bonds bonded to atom in molecule
	int bCount2;			//count of bonds bonded to atom in frag. New added atom should set 1, cause there must be one another atom bond to it in fragment. New added atom should not try to calculate how many another atoms it bond to, see fragment::setFrag()
	std::vector<bond> bonds;	//bonds bonded to atom
	std::vector<int> bondTs;	//bond types around atom in a fragment or molecule
	std::vector<int> dists;		//distance between start and end atom, bcount2 equal 1 is start or end atom
	atom(std::string oneLine, int serialNum, bool H=false);
	atom();
	void setAtom(std::string& atomInfor);
	void setBonds(std::vector<bond> bonds);		//get all bonds bonded to this atom
	bool operator == (const atom& a) const;
	bool operator > (const atom& a) const;
	bool chkAround(atom &a, std::vector<int> &retASNs);	//check bonded atoms, if this atom contains all atoms, true
	bool eraseBond(int serNum);	//delete bond by atom serial number
	bool reviseBond(int serNum, int newSerNum, bool mark=false);	//revise bond by atom serial number
	bool getBond(int serNum, bond& retBond);	//get bond by serial number
	void showAtom();
	void resetUsed();
	std::string toString();
	void clear();
};

inline bool equalV(const std::vector<int>& VA, const std::vector<int>& VB)
{
	int i=0;
	for(i=0; i<VA.size(); i++)
	{
		if(VA[i] != VB[i])
		{
			return false;
		}
	}
	return true;
}

inline bool largerAtom(const atom& a, const atom& b)
{
	if(a.atomNum > b.atomNum)
		return true;
	else if(a.atomNum == b.atomNum)
	{
		if(a.aroma && !b.aroma)
			return true;
		else if(!a.aroma && b.aroma)
			return false;
		if(a.aroma == a.aroma)
		{
			if(a.bCount2 < b.bCount2)   //!!!! bCount or bCount2
				return true;
			else if(a.bCount2 == b.bCount2)
			{
				if(a.serNum < b.serNum)
					return true;
				else 
					return false;
			}
			else
				return false;
		}
		else
			return false;

	}
	else
		return false;

}
#endif
