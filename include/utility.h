#ifndef UTILITY_H
#define UTILITY_H
#pragma warning (disable: 4786)

#include <string>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include "atom.h"
#include "bond.h"
#include "fragment.h"
#include "markush.h"
class fileReader
{
private:
	std::string filePath;
	std::string fileName;
	long fileSize;
	long readedPos;
public:
	int fileType;	//0 unknow, 1 sdf, 2 ctable
	long readedCount;
	bool readToEnd;
	std::vector<std::string> readedObjects;
	fileReader(std::string path);
	void read(char delim='$', long limit=-1);	//limit = -1 means no limit
	void Read();
};

class RGP
{
public:
	std::string all;
	std::string rgpNO;
	std::vector<std::string> frags;

	RGP(std::string &str);
	void setRGP();
	void show();
};

class CTable
{
public:
	std::string all;
	std::string head;	//HDR to END HDR
	std::string baseFrag;	//the first Ctable
	std::vector<RGP> rgps;
	
	CTable(std::string &str);
	void setCTable();
	void createFrag();
	void show();
};

struct simResult
{
	std::string idA;
	std::string idB;
	float sim;
	bool operator < (const simResult& sF) const
	{
		if(this->sim < sF.sim)
			return true;
		else
			return false;
	}
};
inline bool compSR(const simResult &srA, const simResult &srB)
{
	if(srA.sim > srB.sim)
		return true;
	else
		return false;
}

inline bool more(const int& a, const int& b)
{
	if(a > b)
		return true;
	else
		return false;
}


inline bool lessASN(const atom& a, const atom& b)
{
	if(a.serNum < b.serNum)
		return true;
	else 
		return false;
}

inline bool lessBASN(const bond& a, const bond& b)
{
	if(a.bA < b.bA)
		return true;
	else 
		return false;
}


inline bool largerFragAC(const fragment& fA, const fragment& fB)
{
	if(fA.atomCount > fB.atomCount)
		return true;
	else
		return false;
}

inline bool equalFragASN(const fragment& fA, const fragment& fB)
{
	if(fA.atomCount != fB.atomCount)
		return false;

	std::vector<int> asnA;
	std::vector<int> asnB;
	fA.getASN(asnA);
	fB.getASN(asnB);
	std::sort(asnA.begin(), asnA.end(), more);
	std::sort(asnB.begin(), asnB.end(), more);

	if(equalV(asnA, asnB))
		return true;
	else
		return false;
}

class connector
{
public:
	void connectAtom(atom &atomA, atom &atomB,
							bond &retBond, int bondType, bool revise = true);
	bool connectRGP(fragment &fragA, fragment fragB,
					std::vector<atom>& retAtoms, std::vector<bond>& retBonds);
	bool connectRGP(fragment &fragA, fragment fragB, fragment &regFrag);
};

void showVector(std::vector<int> v, int width=10000, int delim=0);
void showVector(std::vector<double> v, int width=10000, int delim=0);
void showVector(std::vector<long> v, int width=10000, int delim=0);
void showVector(std::vector<std::string> v, int width=10000, int delim=0);
void showVector(std::vector<simResult> &v);
long timeNow();
std::string timeElapse(long time1,long time2);
bool getMolFiles(std::string dir, std::vector<std::string>& files); //get mol-file names from dir

//check whether new linear fragment is duplicate
inline bool chkLine(fragment &tmpFrag, std::vector<fragment> &frags, long startPos, long endPos)
{
	if(tmpFrag.linear)
	{
		for(long l=startPos;l<endPos;l++)
		{
			if(tmpFrag.fASN == frags[l].fASN)
				return true;
			else if(l == endPos)
				return false;
		}
	}
	return false;
}

//check whether new ring fragment is duplicate
inline bool chkRing(fragment &tmpFrag, std::vector<fragment> &frags, long startPos, long endPos)
{
	if(tmpFrag.cyclic)
	{
		bool oneAtomExist;
		for(long i=startPos;i<endPos;i++)
		{
			if(frags[i].cyclic && tmpFrag.fASN.size() == frags[i].fASN.size())
			{
				for(int j=0;j<tmpFrag.atomCount;j++)
				{
					oneAtomExist = false;
					for(int k=0;k<frags[i].atomCount;k++)
					{
						if(tmpFrag.fAtoms[j].serNum == frags[i].fAtoms[k].serNum
							&& tmpFrag.fAtoms[j].bCount2 == frags[i].fAtoms[k].bCount2)
						{
							oneAtomExist = true;
							break;
						}
					}
					if(!oneAtomExist)
						break;
					else if(j == tmpFrag.atomCount-1)
						return true;
				}
			}
		}
	}
	return false;
}

inline std::string& trim(std::string& str)
{ 
	std::string buff(str); 
	char space = ' '; 
	str.assign(buff.begin() + buff.find_first_not_of(space), buff.begin() + buff.find_last_not_of(space) + 1); 
	return str;
} 

inline void showStr(std::string str, std::string delim="\r\n")
{
	std::cout<<str<<delim;
}

//form fragment by using atoms, if they can connect to each other
//if each atom is single, return false
inline bool formFrag(std::vector<atom> &atoms, fragment &frag)
{
	int i=0,j=0,k=0,l=0;
	for(i=0;i<atoms.size();i++)
		if(frag.addAtom(atoms[i]))
			for(j=0;j<frag.atomCount;j++)
				for(k=0;k<frag.fAtoms[j].bCount;k++)
					for(l=0;l<atoms.size();l++)
						if(frag.fAtoms[j].bonds[k].bB == atoms[l].serNum)
						{
							frag.addAtom(atoms[l]);
							break;
						}
	if(frag.atomCount > 1)
		return true;
	else
		return false;
}
inline void formFrag(std::vector<atom> &atoms, std::vector<bond> &bonds, fragment &retFrag, bool clear=false)
{
	int i=0, j=0;
	retFrag.fAtoms = atoms;
	retFrag.fBonds = bonds;
	retFrag.atomCount = atoms.size();
	retFrag.bondCount = bonds.size();
	for(i=0;i<atoms.size();i++)
	{
		if(atoms[i].chemSym != 'C' || atoms[i].chemSymEx !=' ')
		{
			atoms[i].nonC = true;
			atoms[i].hAtom = true;
		}
		if(atoms[i].chemSymEx == ' ' && atoms[i].chemSym == 'C')
			for(j=0;j<atoms[i].bCount;j++)
				if(atoms[i].bonds[j].bType != 1)
				{
					atoms[i].nonC = true;
					break;
				}
		if(atoms[i].bCount == 1)
			atoms[i].tAtom = true;
	}
	retFrag.setFrag();
	if(clear)
	{
		atoms.clear();
		bonds.clear();
	}
}

//set A not set B by atom serial number
inline std::vector<atom> AnotB(std::vector<atom> &atomsA, std::vector<atom> &atomsB)
{
	int i=0,j=0;
	std::vector<atom>::iterator it;
	for(i=0;i<atomsA.size();i++)
		for(j=0;j<atomsB.size();j++)
			if(atomsA[i].serNum == atomsB[j].serNum)
			{
				it = atomsA.begin();
				//atomsA.erase(&atomsA[i]);
				atomsA.erase(it+i);
				i = -1;
				break;
			}
	return atomsA;
}
inline void AnotB(std::vector<atom> &atomsA, std::vector<atom> &atomsB, std::vector<atom> &retAtoms)
{
	int i=0, j=0, k=0;
	int sizeA = atomsA.size();
	int sizeB = atomsB.size();
	bool find = false;
	for(i=0;i<sizeA;i++)
	{
		//if(atomsA[i].serNum < atomsB[0].serNum) continue;
		find = false;
		for(j=0;j<sizeB;j++)
		{
			if(atomsA[i].serNum == atomsB[j].serNum)
			{
				find = true;
				break;
			}
		}
		if(!find)
		{
			if(!atomsA[i].RGPAtom)
				retAtoms.push_back(atomsA[i]);
			else
			{
				find = false;
				for(j=0;j<atomsA[i].bCount;j++)
				{
					for(k=0;k<sizeB;k++)
					{
						if(atomsA[i].bonds[j].bB == atomsB[k].serNum)
						{
							find = true;
							break;
						}
					}
					if(find) break;
				}
				if(!find)
					retAtoms.push_back(atomsA[i]);
			}
		}
	}
}

//get rid of none-bonded and one-bonded atoms
inline void shake(std::vector<atom> &atoms)
{
	int i=0,j=0,k=0;
	int count=0;
	std::vector<atom>::iterator it;
	for(i=0;i<atoms.size();i++)
	{
		count = 0;
		for(j=0;j<atoms[i].bCount;j++)
		{
			for(k=0;k<atoms.size();k++)
			{
				if(atoms[i].bonds[j].bB == atoms[k].serNum)
				{
					count++;
					break;
				}
			}
		}
		if(count <= 1)
		{		
			//atoms.erase(&atoms[i]);
			it = atoms.begin();
			atoms.erase(it+i);
			i = -1;
		}
	}
}

//split fragment into two small fragment by cut off one bond
//if success split into two small fragment, return true
inline bool splitFrag(fragment &frag, bond b, fragment &smallF1, fragment &smallF2)
{
	int i=0,j=0,k=0;
	int tmpbB = 0;
	int tmpASN = 0;
	fragment tmpFrag = frag;
	bool recover = false;
	for(i=0;i<frag.atomCount;i++)
	{
		for(j=0;j<frag.fAtoms[i].bCount;j++)
			if((frag.fAtoms[i].bonds[j].bA == b.bA && frag.fAtoms[i].bonds[j].bB == b.bB)
				|| (frag.fAtoms[i].bonds[j].bA == b.bB && frag.fAtoms[i].bonds[j].bB == b.bA))
			{
				tmpbB = frag.fAtoms[i].bonds[j].bB;
				tmpASN = frag.fAtoms[i].serNum;
				frag.fAtoms[i].bonds[j].bB = 0;
				break;
			}
		if(tmpbB != 0)
			break;
	}

	if(tmpbB != 0)
	{
		formFrag(frag.fAtoms,smallF1);
		frag.fAtoms = AnotB(frag.fAtoms,smallF1.fAtoms);
		formFrag(frag.fAtoms,smallF2);
		frag = tmpFrag;
	}

	if(smallF1.atomCount > 0 && smallF2.atomCount > 0)
	{
		bool revise = false;
		for(i=0;i<smallF1.atomCount;i++)
		{
			if(smallF1.fAtoms[i].serNum == tmpASN)
			{
				for(j=0;j<smallF1.fAtoms[i].bCount;j++)
				{
					if(smallF1.fAtoms[i].bonds[j].bB == 0)
					{
						smallF1.fAtoms[i].bonds[j].bB = tmpbB;
						revise = true;
						break;
					}
				}
				if(revise) break;
			}
		}
		if(!revise)
		{
			for(i=0;i<smallF2.atomCount;i++)
			{
				if(smallF2.fAtoms[i].serNum == tmpASN)
				{
					for(j=0;j<smallF2.fAtoms[i].bCount;j++)
					{
						if(smallF2.fAtoms[i].bonds[j].bB == 0)
						{
							smallF2.fAtoms[i].bonds[j].bB = tmpbB;
							revise = true;
							break;
						}
					}
					if(revise) break;
				}
			}
		}

		return true;
	}
	else
		return false;
}

//split fragment into pieces by branch-point bonds, saved in retrunFrags
//if split into no less than two pieces, return true
inline bool shake(fragment &frag, std::vector<fragment> &retrunFrags)
{
	int i=0,j=0;
	int size = 0;
	bool cut = false;
	int eraseTarget = 0;
	fragment frag1;
	fragment frag2;
	fragment tmpFrag;
	std::vector<fragment>::iterator it = retrunFrags.begin();
	for(i=0;i<frag.atomCount;i++)
	{
		if(frag.fAtoms[i].bCount >= 3)
		{
			for(j=0;j<frag.fAtoms[i].bCount;j++)
			{
				if(frag.chkBAB(frag.fAtoms[i].bonds[j]))
				{
					frag1.clear();
					frag2.clear();
					if(splitFrag(frag,frag.fAtoms[i].bonds[j],frag1,frag2))
					{
						cut = true;
						shake(frag1.fAtoms);
						if(formFrag(frag1.fAtoms,tmpFrag))
							retrunFrags.push_back(tmpFrag);
						tmpFrag.clear();
						shake(frag2.fAtoms);
						if(formFrag(frag2.fAtoms,tmpFrag))
							retrunFrags.push_back(tmpFrag);
						tmpFrag.clear();
						break;
					}
				}
			}
		}
		if(cut)
			break;
	}

	if(cut)
	{	std::vector<fragment>::iterator itF;
		for(i=0;i<retrunFrags.size();i++)
		{	
			itF = retrunFrags.begin();
			tmpFrag = retrunFrags[i];
			retrunFrags.erase(itF + i);
			//retrunFrags.erase(&retrunFrags[i]);
			shake(tmpFrag,retrunFrags);
		}
	}
	else
		retrunFrags.insert(it,frag);
	
	return cut;
}

//show information of how to use 
inline void showHtU(int mode=1)
{
	extern std::string CRLF;

	std::cout<<CRLF<<"Molengine 2012-12, east linden co,. ltd.";
	std::cout<<CRLF<<"For chemical structure retrieve";
	std::cout<<CRLF<<"Three main utilities:";
	std::cout<<CRLF<<"\tcalcSim: similarity calculation";
	std::cout<<CRLF<<"\tsubSMatch: sub-structure match";
	std::cout<<CRLF<<"\tSMatch: structure match";
	std::cout<<CRLF<<CRLF<<"Use: ";
	if(mode == 1)
	{
		std::cout<<"molengine mode [threshold] file1 file2 result";
		std::cout<<CRLF<<"\tmolengine calcSim 0.9 A.sdf B.sdf output.txt";
		std::cout<<CRLF<<"\tmolengine subSMatch A.sdf sub.sdf output.txt";
		std::cout<<CRLF<<"\tmolengine SMatch A.sdf B.sdf output.txt\r\n";
	}
	if(mode == 2)
	{
		std::cout<<CRLF<<"Follow next guide to creat your molecule database";
		std::cout<<CRLF<<"All three utilities will base on that";
	}
}

//check input arguments, if legal true
inline bool chkInput(int argc, char* argv[])
{
	bool OK = false;
	if(argc < 5 || argc > 6)
	{
		showHtU();
	}
	else if(argc == 5)
	{
		if(strcasecmp(argv[1],"subSMatch") != 0 && strcasecmp(argv[1],"SMatch") != 0)
		{
			showHtU();
		}
		else
			OK = true;
	}
	else if(argc == 6)
	{
		if(strcasecmp(argv[1],"calcSim") != 0)
		{
			showHtU();
		}
		else
			OK = true;
	}

	return OK;
}

//check two int vector, if contains at least one same element, true
inline bool chkMix(std::vector<int>& v1, std::vector<int>& v2)
{
	int i=0, j=0;
	bool result = false;
	for(i=0;i<v1.size();i++)
	{
		for(j=0;j<v2.size();j++)
		{
			if(v1[i] == v2[j])
			{
				result = true;
				break;
			}
		}
		if(result)
			break;
	}
	return result;
}

//get max atom serial number in atom array
inline int getMaxAsn(std::vector<atom>& as)
{
	int result = 0;
	if(as.size() == 0) return result;

	int i=0;
	result = as[0].serNum;
	for(i=0;i<as.size();i++)
		result = as[i].serNum > result ? as[i].serNum : result;

	return result;
}

//put large fragment in the first
inline void sortFrag(std::vector<fragment>& frags)
{
/*
	int i=0, j=0;
	fragment tmpFrag;
	for(i=0;i<frags.size();i++)
	{
		for(j=i+1;j<frags.size();j++)
		{
			if(frags[i].atomCount < frags[j].atomCount)
			{
				tmpFrag.clear();
				tmpFrag = frags[i];
				frags[i] = frags[j];
				frags[j] = tmpFrag;
			}
		}
	}
*/
	std::sort(frags.begin(), frags.end(), largerFragAC);

}

//check two fragments' asn, if asn of fragA contain all asn of fragB, true
inline bool chkFragASN(fragment& fragLarge, fragment& fragSmall)
{
	int i=0;
	for(i=0;i<fragSmall.fAtoms.size();i++)
	{
		if(!fragLarge.chkASN(fragSmall.fAtoms[i].serNum))
			return false;
	}
	return true;
}

//change bond bA and bB, base on one number
inline void changeB(std::vector<bond>& bonds, int offset)
{
	int i=0;
	for(i=0;i<bonds.size();i++)
	{
		bonds[i].bA = bonds[i].bA + offset;
		bonds[i].bB = bonds[i].bB + offset;
	}
}

//convert to string
inline std::string convertToStr(int &number)
{
	char tmp[999999];
	sprintf(tmp,"%d",number);
	return tmp;
}

//save in file
inline void saveFile(const char *file, std::vector<simResult>& vsr)
{
	long i=0;
	std::fstream savefile;
	extern std::string CRLF;
	savefile.open(file, std::ios::out | std::ios::trunc);

	for(i=0;i<vsr.size();i++)
	{
		savefile<<vsr[i].idA<<"\t"<<vsr[i].idB<<"\t"<<vsr[i].sim<<CRLF;
	}

	savefile.close();
}

inline void splitStr(std::string inputStr, std::vector<std::string>& retStrs, std::string delim)
{
	if(inputStr.find(delim) != std::string::npos)
	{
		std::string tmpStr = "";
		do
		{
			tmpStr = inputStr.substr(0, inputStr.find(delim));
			retStrs.push_back(tmpStr);
			inputStr = inputStr.substr(inputStr.find(delim) + delim.length());
		}while(inputStr.find(delim) != std::string::npos);
	}
	retStrs.push_back(inputStr);
}

//remove duplication in string
inline void modifyStr(std::string &inputStr, std::string &outputStr, std::string delim)
{
	std::vector<std::string> tmpStrV;
	std::vector<std::string>::iterator its;
	std::vector<std::string>::iterator its2;
	bool same = false;
	
	splitStr(inputStr, tmpStrV, delim);
	
	for(its = tmpStrV.begin(); its != tmpStrV.end(); its++)
	{
		its2 = its + 1;
		same = false;
		if(its2 != tmpStrV.end())
		{
			for(; its2 != tmpStrV.end(); its2++)
			{
				if(*its == *its2)
				{
					same = true;
					break;
				}
			}
		}
		if(!same)
		{
			if(its != tmpStrV.begin())
				outputStr.append(delim);
			outputStr.append(*its);
		}
	}
}

#endif
