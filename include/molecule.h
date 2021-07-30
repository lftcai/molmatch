#ifndef MOLECULE_H
#define MOLECULE_H

#include <string>
#include <vector>
#include <iomanip>
#include "atom.h"
#include "bond.h"
#include "fragment.h"
#include "markush.h"
#include "utility.h"

bool chkMol(std::string inputTmp);
bool chkCtab(std::string inputTmp);
void CtabToMol(CTable& ctab, std::string& retMol);

class molecule
{
public:
	std::string molName;    //first line molecule name
	std::string molBuilder; //the second line
	std::string molComment; //the third line
	std::string molInfor;   //the forth line, the counts line
	std::string molFileVer; //file version
	std::string molBody;	//other lines beside 1,2,3wline
	std::string molMap;          //molecule coordinate and connection
	std::string molExtInfor;     //informiation between connection map and mol end mark
	std::string molProps;        //molecule property 
	std::vector<std::string> propNames;
	std::vector<std::string> propValues;
	int propCount;
	int atomCount;
	int bondCount;
	long fragCount;	//fragment count
	int sIonCount;	//salt ion count, include non-ion
	int arCount;	//aromatic rings count
	bool noMolStru;
	bool molExt;     //if has molExtInfor return true
	bool rgroup;	//if rgroup molecule, true;
	bool ctable;	//if got from ctable file, true
	bool salt;       //if exist as salt, true
	bool ion;	//if all salt ion are iorganic, true	
	bool aroma;      //if exist aromatic ring, true
	bool polymer;	//if polymer (contain * atom), true
	//int saltFrag;    //how many fragments (salt ion, water and etc.) consist of this molecule
	int charge;      //molecule charge value
	long molID;
	static long molCount;	//molecule count
	static long saltCount;	//salt ion count, include non-salt fragment
	static long errMolCount;
	std::vector<atom> atoms;
	std::vector<bond> bonds;
	std::vector<fragment> frags;	//setFrag() generated fragments saved here
	std::vector<fragment> salts;	//setSalt() found salt ions (parts) saved here
	std::vector<atom> osAtoms;		//oganic salt ions (parts)' atom
	std::vector<fragment> aRings;	//setAroma() found aromatic rings saved here, large rings
	std::vector<rGroupMDL> MDLRGP;	//MDL rgroup information
	std::vector<rGroupCtab> CtabRGP;//CTable markush R LOG
	std::vector<fragment> RGPFras;	//rgroup represented all frags;
	
	molecule(std::string tmp);
	molecule();
	void createMol(std::string tmp);
	bool init();    //call the next three function, which is basicly get information for file
	void setABC();  //set atom count, bond count and other information
	void setCCP();  //set molecule strutcure information and molecule property
	void setProp(); //set propertyies
	void prepare();	//call several functions to get ready for CalcSim, StrucMatch and SubSMatch
	bool chkStarAtom();	//check if contain * atom
	void addRAtom();	//add a r atom, use in ctable
	void setRGroup();	//decode rgroup information, save represented frags in RGPfrags
	void setSalt(); //save salt ion as fragment
	bool setOSAtoms();	//if if exist as salt, get all organic atoms
	bool chkIon();	//check if all salt ion are iorganic
	void setAroma();	//save aromatic ring as fragment, set bond type as 4
	void setSaltAroma();//set fragment::aroma of each salt ion
	void delH();    //delete H
	void sortAB();	//sort atoms, sort bonds
	void setAB();   //set atom and its bonds
	void setBC();	//set bonds' chemical symbol
	void setStereo();	//set atoms stereo parity according to bond connection map
	void setFragEx(int num = 0);	//set frag's detail information 
	void getProp(int index, std::string &propN, std::string &propV);//index >=1
	void getProp(std::vector<std::string> &propN, std::vector<std::string> &propV);
	void showSalt();
	void showAroma();
	void showMol();
	void outMol(std::fstream& outFile);
	void saveSDF(bool salt, bool removeEmpty = true, bool reviseBond = true);	//if salt, save each salt ion, else save deleted h molecule
	static void getCount(long &total, long &errMol)//get molecule count
	{
		total = molCount;
		errMol = errMolCount;
	}
};

#endif
