#ifndef FRAGMENT_H
#define FRAGMENT_H

#include "atom.h"
#include "bond.h"
#include "markush.h"
#include <vector>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <algorithm>

class fragment
{
public:
	std::vector<atom> fAtoms;
	std::vector<bond> fBonds;
	std::vector<atom> hAtoms;	//hetero atoms
	std::vector<markush> markushs;	//markush information, see class markush
	std::string fASN;	//fragment atom serial number string
	int atomCount;
	int bondCount;
	int maxASN;		//max atom serial number in fragment
	bool linear;	//line
	bool cyclic;	//ring
	bool branched;	//with branch
	int rCount;		//ring count
	int charge;		//chemical charge which is different from atom charge, see atom charge
	bool ioIon;		//if inorganic salt, true
	bool oIon;		//if organic salt, true
	bool baseRGP;	//if base group, true
	bool aroma;		//if aromatic, true
	bool hetero;	//if hetero fragment, true
	bool sameHA;	//if same hetero aomts, true
	int hAtomC;		//hetero atom count
	int aAtomC;		//aromatic atom count
	int maxRGNO;	//max rgroup NO.
	bool used;		//if this is used, true

	fragment();
	fragment(std::vector<atom> Atoms, std::vector<bond> Bonds);
	void setFrag();			//check line, branch and ring, set bount2, add bond if missed, get bonds and atoms ready
	void setHetero();		//set hetero information
	bool addAtom(atom &a);	//add atom by serial number, if success, return true
	bool addAtom2(atom &a, bool markUsed);	//all most same with addAtom(), but for sub structure match
	void setMarkush();		//get markush ready
	void setSalt();			//set charge, check organic or inorganic
	void setFASN();			//set all atoms' serial number into a string fragment::fASN
	bool chkASN(int asn);	//check all atoms' serial number, if already exist return true
	bool chkBAB(bond b);	//check bond b's bA and bB to determine whether same bond exists int this->fBonds, if exist return true
	bool chkAtomSym(atom &a);	//if this fragment contain atom chemical symble same as atom a,true
	bool chkAroma();		//check whether this fragment is aromatic, for molecule::setAroma
	bool chkMarkush(markush &m);		//if m exist in markushs, true
	void setMaxASN();		//find out the max atom serial number
	void setABTs();			//set atoms' bonds type, see atom::bondTs
	void setADs();			//set atoms's distance, see atom::dist
	void findTA();			//find out terminal atoms
	void findAA();			//find out aromatic atoms
	void getASN(std::vector<int>& asns) const;	//get atom serial number array
	bool getAtom(atom &a, std::vector<atom>& returnAtoms, bool oneTime=false);	//if this fragment has one atom as same as a, return this atom, one atom can be got one time cause "it be used"
	bool getAtom(int serNum, atom& returnAtom, bool oneTime=true); //get one atom by serNum, one atom can be got one time
	bool nextAtom(atom &a, std::vector<atom>& retAtoms, int &count,bool oneTime=true);	//get atoms bonded to atom a
	void getNMAtoms(std::vector<atom>& retAtoms);		//get non-markush atoms
	void getNMBonds(std::vector<bond>& retBonds);		//get non-markush bonds
	bool operator == (fragment& frag);	//setFrag(), setABTs(), setADs() should run before
	void clear();			//reset all member variable
	bool eraseAtom(atom &a);	//delete one atom in the fragment
	bool eraseBond(int serNum);	//delete bonds by serial number
	bool reviseAtom(atom &targetAtom);		//revise atom into targetAtom
	void reviseASN(int asn);//revise all atoms' atom serial number, start at asn
	void formateASN();		//revise all atoms' atom serial number by fAtoms's order
	void clearABTs();		//clear atom bond types
	void clearADs();		//clear atom distances
	void resetUsed();		//set all atoms' used = false
	void resetUsed2();		//set all bonds' used = false;
	void resetNewAdded();	//set all atoms' newAdded = false;
	void sortAtoms();		//sort atoms by atom serial number
	void sortAtoms2();		//sort atoms by atom type, bonds count
	void showFrag();		//show this fragment
	void showMol(std::string& name, std::string& molInfor, int baseCount);
	void showHetero();		//show hetero information
};

inline bool chkABTs(std::vector<int>& large, std::vector<int>& sub)
{
	if(large.size() < sub.size())
		return false;
	else
	{
		int i=0;
		for(i=0; i<sub.size(); i++)
		{
			if(sub[i] != large[i])
				return false;
		}
		return true;
	}
}

#endif
