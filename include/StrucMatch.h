#ifndef STRUCMATCH
#define STRUCMATCH

#include "fragment.h"
#include "molecule.h"
#include "utility.h"
extern std::string CRLF;

inline void getReadyMatch(molecule &mol)
{
	mol.init();
	mol.delH();
	mol.setAroma();
	mol.setSalt();
	mol.setProp();
	if(mol.salt)
		for(int i=0;i<mol.sIonCount;i++)
		{
			mol.salts[i].setABTs();
			mol.salts[i].setHetero();
			mol.salts[i].findAA();
			mol.salts[i].setSalt();
		}
	else
	{
		mol.salts[0].setABTs();
		mol.salts[0].setHetero();
		mol.salts[0].findAA();
		mol.salts[0].setSalt();
	}
}

//check some basic property before calulate, if false, the next calculation is needless
inline bool filter(molecule &molA, molecule &molB)
{
	if( molA.rgroup || molB.rgroup)
		return true;
	else if(  
		molA.atomCount != molB.atomCount
		||molA.bondCount != molB.bondCount
		||molA.salt != molB.salt
		||molA.sIonCount != molB.sIonCount
		||molA.aroma != molB.aroma
		||molA.arCount != molB.arCount
		||molA.polymer != molB.polymer
		)
		return false;

	return true;
}
//check some basic property before calulate, if false, the next calculation is needless
inline bool filter(fragment &fa, fragment &fb)
{
	if(    fa.atomCount != fb.atomCount
		|| fa.bondCount != fb.bondCount
		|| fa.linear != fb.linear
		|| fa.cyclic != fb.cyclic
		|| fa.branched != fb.branched
	//	|| fa.rCount != fb.rCount
	//	|| fa.charge != fb.charge
		|| fa.oIon != fb.oIon
		|| fa.ioIon != fb.ioIon
		|| fa.aAtomC != fb.aAtomC
		|| fa.hAtomC != fb.hAtomC
		)
		return false;
	return true;
}


inline bool Compare(fragment &fragA, fragment &fragB)
{
	if(fragA == fragB)
	{
		fragment tmpFragA = fragA;
		fragment tmpFragB = fragB;
		tmpFragA.setADs();
		tmpFragB.setADs();
		if(tmpFragA == tmpFragB)
		{
			return true;
		}
		else
			return false;
	}
	else
		return false;
}

inline bool CompareRGP(std::vector<fragment>& fragsA, std::vector<fragment>& fragsB)
{
	int i=0, j=0;
	for(i=0;i<fragsA.size();i++)
	{
		for(j=0;j<fragsB.size();j++)
		{
			if(!filter(fragsA[i],fragsB[j])) continue;
			else
			{
				if(Compare(fragsA[i],fragsB[j]))
					return true;
			}
		}
	}
	return false;
}

inline bool Compare(std::vector<fragment>& fragsA, std::vector<fragment>& fragsB)
{
	int i=0, j=0;
	for(i=0;i<fragsA.size();i++)
	{
		for(j=i;j<fragsB.size();j++)
		{
			if(!filter(fragsA[i],fragsB[j])) return false;
			else
			{
				if(Compare(fragsA[i],fragsB[j]))
					break;
				else
					return false;
			}
		}
	}
	return true;
}

inline void StrucMatch(std::vector<molecule> &molAs, std::vector<molecule> &molBs, std::vector<simResult> &results)
{
	//std::cout<<"Processing...";

	long i=0, j=0;
	simResult tmp;
	for(i=0;i<molAs.size();i++)
	{
		for(j=0;j<molBs.size();j++)
		{
			if(!molAs[i].rgroup && !molBs[j].rgroup)
			{
				if(!filter(molAs[i],molBs[j]))
					continue;
				if(Compare(molAs[i].salts, molBs[j].salts))
				{
					if(molAs[i].propValues.size() > 0 && molBs[j].propValues.size() >0)
					{
						tmp.idA = molAs[i].propValues[0];
						tmp.idB = molBs[j].propValues[0];
					}
					else
					{
						tmp.idA = molAs[i].molName;
						tmp.idB = molBs[j].molName;
					}
					tmp.sim = 1;
					results.push_back(tmp);
				}
			}
			else
			{
				if(CompareRGP(molAs[i].salts, molBs[j].salts))
				{
					if(molAs[i].propValues.size() > 0 && molBs[j].propValues.size() >0)
					{
						tmp.idA = molAs[i].propValues[0];
						tmp.idB = molBs[j].propValues[0];
					}
					else
					{
						tmp.idA = molAs[i].molName;
						tmp.idB = molBs[j].molName;
					}
					tmp.sim = 1;
					results.push_back(tmp);
				}
			}
		}
	}

	//std::cout<<"done"<<CRLF;
}

//check if two molecule set have same molecule
inline bool chkMolExist(std::vector<molecule> &molAs, std::vector<molecule> &molBs, bool reviseName = false)
{
	long i=0, j=0;
	for(i=0;i<molAs.size();i++)
	{
		for(j=0;j<molBs.size();j++)
		{
			if(!molAs[i].rgroup && !molBs[j].rgroup)
			{
				if(!filter(molAs[i],molBs[j]))
					continue;
				if(Compare(molAs[i].salts, molBs[j].salts))
				{
					if(reviseName)
					{
						std::string names = "";
						names = molAs[i].molName + "+" + molBs[j].molName;
						molAs[i].molName = names;
						molBs[j].molName = names;
					}
					return true;
				}
			}
			else if(molAs[i].rgroup && molBs[j].rgroup && (molAs[i].sIonCount == molBs[j].sIonCount))
			{
				if(Compare(molAs[i].salts, molBs[j].salts))
				{
					if(reviseName)
					{
						std::string names = "";
						names = molAs[i].molName + "+" + molBs[j].molName;
						molAs[i].molName = names;
						molBs[j].molName = names;
					}
					return true;
				}
			}
		}
	}
	return false;
}
#endif
