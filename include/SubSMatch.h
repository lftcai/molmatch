#ifndef SUBSMATCH
#define SUBSMATCH

#include "molecule.h"
#include "fragment.h"
#include "utility.h"
extern std::string CRLF;
extern int beginPos;

inline void getReadySubMatch(molecule &mol)
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
inline bool filterSub(molecule &Structure, molecule &SubStructure)
{
	if(Structure.rgroup || SubStructure.rgroup)
		return true;
	else if(Structure.atomCount < SubStructure.atomCount
		|| Structure.bondCount < SubStructure.bondCount)
		return false;
	else if(SubStructure.aroma && !Structure.aroma)
		return false;
	else if(SubStructure.polymer && !Structure.polymer)
		return false;
	else if(SubStructure.arCount > Structure.arCount)
		return false;
	else if(SubStructure.atomCount > 20 || Structure.atomCount > 30)
		return false;
	//else if(SubStructure.salt && !Structure.salt)
		//return false;
	return true;
}

//check some basic property before calulate, if false, the next calculation is needless
inline bool filterSub(fragment &frag, fragment &subfrag)
{
	if(frag.atomCount < subfrag.atomCount
		|| frag.bondCount < subfrag.bondCount)
		return false;
	else if(subfrag.cyclic && !frag.cyclic)
		return false;
	else if(subfrag.branched && !frag.branched)
		return false;
	else if(subfrag.oIon && !frag.oIon)
		return false;
	else if(subfrag.ioIon && !frag.ioIon)
		return false;
	else if(subfrag.aroma && !frag.aroma)
		return false;
	else if(subfrag.aAtomC > frag.aAtomC)
		return false;
	else if(subfrag.hetero && !frag.hetero)
		return false;
	else if(subfrag.hAtomC > frag.hAtomC)
		return false;

	return true;
}


//get start of sub structure match calculation
//get subfrag's first atom
//get frag's same atom like subfrag's first atom, construct all this atoms into frags
inline bool begin(fragment &frag, fragment &subfrag,
				  std::vector<fragment> &frags, std::vector<atom> &atoms)
{
//std::cout<<"\nin begin\n";
	int i=0;
	std::vector<atom> tmpAtoms;
	fragment tmpFrag;
	atoms.push_back(subfrag.fAtoms[0]);
	subfrag.fAtoms[0].used = true;


//std::c(ut<<"subFrag:\n";
//subfrag.fAtoms[0].showAtom();
//std::cout<<"\nLfrag:\n";
	if(frag.getAtom(subfrag.fAtoms[0],tmpAtoms,false))
	{
		for(i=0;i<tmpAtoms.size();i++)
		{
			tmpFrag.addAtom(tmpAtoms[i]);
			frags.push_back(tmpFrag);
//tmpFrag.showFrag();
			tmpFrag.clear();
		}
		return true;
	}
	else
		return false;
}

//get subfrag's next atoms
inline bool move(fragment &subfrag, std::vector<atom> &atoms, int &count)
{
//std::cout<<"\nin move single\n";
	count=0;
	subfrag.nextAtom(atoms[0], atoms, count, true);
//std::cout<<"subFrag:\n";
//for(int i=0;i<atoms.size();i++) atoms[i].showAtom();
	if(count > 0)
		return true;
	else
		return false;
}

//small frag grow up if atom a can be added into
inline bool move(fragment frag, atom &a, fragment &Lfrag,
				 std::vector<fragment> &frags, bool &wait, bool &markUsed)
{
//std::cout<<"Lfrag:\n";
	int i=0, j=0;
	bool added = false;
	int pos = 0;
	atom tmpAtom;
	fragment tmpFrag;
	std::vector<int> ASN;
	for(i=0;i<frag.atomCount;i++)
	{
		if(wait && frag.fAtoms[i].newAdded) continue;
		if(!frag.fAtoms[i].used && frag.fAtoms[i].chkAround(a, ASN))
		{
			for(j=0;j<ASN.size();j++)
			{
				if(Lfrag.getAtom(ASN[j], tmpAtom, false))
				{
					tmpFrag.clear();
					tmpFrag = frag;
					if(tmpFrag.addAtom2(tmpAtom, markUsed))
					{
						added = true;
						frags.push_back(tmpFrag);
//tmpFrag.showFrag();
					}
					tmpAtom.clear();
				}
			}
			ASN.clear();
		}
	}
	return added;
}

//remove duplicate frag
inline void remDup(std::vector<fragment>& frags)
{
	if(frags.size() <= 1) return;
	
	int i=0, j=0;
	std::vector<fragment>::iterator itF;
	for(i=0; i<frags.size(); i++)
	{
		for(j=i+1; j<frags.size(); j++)
		{
			if(equalFragASN(frags[i], frags[j]))
			{
				itF = frags.begin();
				//frags.erase(&frags[j]);
				frags.erase(itF + j);
				j -= 1;
			}
		}
	}
}

//move a small frag same with subfrag
inline bool moveBoth(fragment &subfrag, std::vector<atom> &atoms,
					 fragment &Lfrag, std::vector<fragment> &frags)
{
	int i=0, j=0;
	int count = 0;
	int growSize=0;
	bool result = false;
	bool wait = false;
	bool markUsed = false;
	std::vector<atom>::iterator itA;
	std::vector<fragment>::iterator itF;
	while(atoms.size() > 0)
	{
		if(move(subfrag,atoms,count))
		{
			for(i=0;i<count;i++)
			{
				result = false;
				growSize = frags.size();
				if(i == 0)
					wait = false;
				else
					wait = true;
				if(count > 1)
				{
					if(i==count-1)
						markUsed = true;
					else
						markUsed = false;
				}
				else
					markUsed = true;
				for(j=0;j<growSize;j++)
				{
					move(frags[0],atoms[atoms.size()-1-i],Lfrag,frags,wait,markUsed);
					//frags.erase(&frags[0]);
					itF = frags.begin();
					frags.erase(itF);
				}
			}
			remDup(frags);
			if(frags.size() == 0)
			{
				result = false;
				break;
			}
			else
				result = true;
		}
		//atoms.erase(&atoms[0]);
		itA = atoms.begin();
		atoms.erase(itA);
	}
	return result; 
}

//move all small frags same with subfrag
inline bool moveAll(fragment &subfrag, std::vector<atom> &atoms,
					 fragment &Lfrag, std::vector<fragment> &frags)
{
	if(frags.size() == 0) return false;

	std::vector<fragment> tmpFrags;
	std::vector<atom> tmpAtoms;
	fragment tmpSubFrag;
	fragment tmpFrag;

	int i=0, j=0;
	bool result = false;

	//tmpSubFrag.setADs();
	for(i=0;i<frags.size();i++)
	{
		tmpSubFrag = subfrag;
		tmpAtoms = atoms;
		tmpFrags.clear();
		tmpFrags.push_back(frags[i]);
		if(moveBoth(tmpSubFrag,tmpAtoms,Lfrag,tmpFrags))
		{
			for(j=0;j<tmpFrags.size();j++)
			{
				tmpFrags[j].setFrag();
				tmpFrags[j].resetUsed();
				tmpFrags[j].sortAtoms2();
				tmpFrags[j].setABTs();
				if(tmpSubFrag == tmpFrags[j])
				{
					tmpFrag = tmpSubFrag;
					tmpFrag.setADs();
					tmpFrag.resetUsed();
					tmpFrags[j].setADs();
					if(tmpFrag == tmpFrags[j])
					{
						result = true;
						break;
					}
				}
			}
		}
		if(result)
			break;
	}

	return result;
}

inline bool compare(fragment &Lfrag, fragment &subFrag)
{
	if(Lfrag.atomCount == subFrag.atomCount)
		if(Lfrag == subFrag)
			return true;

	if(subFrag.atomCount == 1)
		if(Lfrag.chkAtomSym(subFrag.fAtoms[0]))
			return true;
		else
			return false;

	std::vector<fragment> frags;
	std::vector<atom> atoms;
	bool result = false;
	if(begin(Lfrag,subFrag,frags,atoms))
	{
		if(moveAll(subFrag,atoms,Lfrag,frags))
		{
			result = true;
		}
		else
			result = false;
	}
	else
		result = false;
	subFrag.fAtoms[0].used = false;
	return result;
}
inline bool compare(std::vector<fragment>& Lfrags, std::vector<fragment>& subFrags)
{
	int i=0, j=0;
	for(i=0;i<Lfrags.size();i++)
	{
		for(j=0;j<subFrags.size();j++)
		{
			if(!filterSub(Lfrags[i],subFrags[j])) continue;
			else
			{
				if(compare(Lfrags[i],subFrags[j]))
					return true;
			}
		}
	}
	return false;
}

inline void SubSMatch(std::vector<molecule>& Lmols, std::vector<molecule>& subMols, std::vector<simResult>& results)
{
	//std::cout<<"Processing...";
	long i=0, j=0, k=0, l=0;
	simResult tmp;
	int sameCount = 0;

	for(i=0;i<Lmols.size();i++)
	{
		for(j=0;j<subMols.size();j++)
		{
			if(!Lmols[i].rgroup && !subMols[j].rgroup)
			{
				if(!filterSub(Lmols[i],subMols[j]))
					continue;
				if(!subMols[j].salt)
				{
					for(k=0;k<Lmols[i].sIonCount;k++)
					{
						if(!filterSub(Lmols[i].salts[k],subMols[j].salts[0]))
							continue;
						if(compare(Lmols[i].salts[k],subMols[j].salts[0]))
						{
							if(Lmols[i].propValues.size() > 0 && subMols[j].propValues.size() > 0)
							{
								tmp.idA = Lmols[i].propValues[0];
								tmp.idB = subMols[j].propValues[0];
							}
							else
							{
								tmp.idA = Lmols[i].molName;
								tmp.idB = subMols[j].molName;
							}
							tmp.sim = -1;
							results.push_back(tmp);
							break;
						}
					}
				}
				else if(subMols[j].salt)
				{
					sameCount = 0;
					for(k=0;k<Lmols[i].sIonCount;k++)
					{
						for(l=0;l<subMols[j].sIonCount;l++)
						{
							if(!filterSub(Lmols[i].salts[k],subMols[j].salts[l]))
								continue;
							if(compare(Lmols[i].salts[k],subMols[j].salts[l]))
								sameCount++;
						}
					}
					if(sameCount >= subMols[j].sIonCount)
					{
						if(Lmols[i].propValues.size() > 0 && subMols[j].propValues.size() > 0)
						{
							tmp.idA = Lmols[i].propValues[0];
							tmp.idB = subMols[j].propValues[0];
						}
						else
						{
							tmp.idA = Lmols[i].molName;
							tmp.idB = subMols[j].molName;
						}
						tmp.sim = -1;
						results.push_back(tmp);
					}
				}
			}
			else
			{
				if(compare(Lmols[i].salts,subMols[j].salts))
				{
					if(Lmols[i].propValues.size() > 0 && subMols[j].propValues.size() > 0)
					{
						tmp.idA = Lmols[i].propValues[0];
						tmp.idB = subMols[j].propValues[0];
					}
					else
					{
						tmp.idA = Lmols[i].molName;
						tmp.idB = subMols[j].molName;
					}
					tmp.sim = -1;
					results.push_back(tmp);
				}
			}
		}
	}

	//std::cout<<"done"<<CRLF;
}

#endif
