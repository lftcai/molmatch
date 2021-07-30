#ifndef CALCSIM
#define CALCSIM

#include "atom.h"
#include "fragment.h"
#include "molecule.h"
#include "utility.h"

extern std::string CRLF;
inline void getReadySim(molecule &mol)
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
			mol.setOSAtoms();
		}
	else
		mol.salts[0].setABTs();
}


//calculate the similarity between molecule A and B
//by check their atom array
inline float compare(std::vector<atom> &atomsA, std::vector<atom> &atomsB, float threshold)
{
	int i=0, j=0;
	int count=0;
	int total = atomsA.size() + atomsB.size();
	int min = atomsA.size() < atomsB.size() ? atomsA.size() : atomsB.size();
	float sim = 0;
	bool needless = false;

	if((1.0 * min / (total - min)) < threshold)
		return sim;

	for(i=0;i<atomsA.size();i++)
	{
		for(j=0;j<atomsB.size();j++)
		{
			if(atomsB[j].used) continue;
			if(atomsA[i] == atomsB[j])
			{
				atomsB[j].used = true;
				count++;
				break;
			}
		}
		if(1.0*(count+atomsA.size()-i-1)/(total-count-atomsA.size()+i+1) < threshold)
		{
			needless = true;
			break;
		}
	}
	for(i=0;i<atomsB.size();i++)
		atomsB[i].used = false;

	if(needless)
		sim = 0;
	else
		sim = (1.0 * count / (total - count));

	return sim;
}

//calculate the similarity between molecule A and B
//by check their largest fragments
//only two fragments comparation
inline float compare(fragment &fragA, fragment &fragB, float threshold)
{
	int i=0, j=0;
	int count=0;
	int total = fragA.atomCount + fragB.atomCount;
	int min = fragA.atomCount < fragB.atomCount ? fragA.atomCount : fragB.atomCount;
	float sim = 0;
	bool needless = false;

	if((1.0 * min / (total - min)) < threshold)
		return sim;

	for(i=0;i<fragA.atomCount;i++)
	{
		//if(this->fAtoms[i].used) continue;
		for(j=0;j<fragB.atomCount;j++)
		{
			if(fragB.fAtoms[j].used) continue;
			if(fragA.fAtoms[i] > fragB.fAtoms[j]) break;
			if(fragA.fAtoms[i] == fragB.fAtoms[j])
			{
				//this->fAtoms[i].used = true;
				fragB.fAtoms[j].used = true;
				count++;
				break;
			}
		}
		
		if(1.0*(count+fragA.atomCount-i-1)/(total-count-fragA.atomCount+i+1) < threshold)
		{
			needless = true;
			break;
		}
		
	}

	for(i=0;i<fragB.atomCount;i++)
		fragB.fAtoms[i].used = false;

	if(needless)
		sim = 0;
	else
		sim = (1.0 * count / (total - count));

	return sim;
}
inline float compare(std::vector<fragment>& fragsA, std::vector<fragment>& fragsB, float threshold)
{
	float s = 0;
	float tmp = 0;
	int i=0, j=0;
	for(i=0;i<fragsA.size();i++)
	{
		for(j=0;j<fragsB.size();j++)
		{
			tmp = compare(fragsA[i],fragsB[j],threshold);
			s = tmp > s ? tmp : s;
		}
	}
	return s;
}

//calculate similarity between two molecule array
//saved result in a array
inline void CalcSim(std::vector<molecule> &molAs, std::vector<molecule> &molBs, std::vector<simResult> &simResults, float threshold)
{
	//std::cout<<"Processing...";
	simResult tmp;
	long i=0, j=0;
	fragment fragA, fragB;

	for(i=0;i<molAs.size();i++)
	{
		for(j=0;j<molBs.size();j++)
		{
			if(!molAs[i].rgroup && !molBs[j].rgroup)
			{
				if(!molAs[i].salt && !molBs[j].salt)
				{
					tmp.sim = compare(molAs[i].salts[0], molBs[j].salts[0], threshold);
					if(tmp.sim == 1)
					{
						fragA.clear();
						fragB.clear();
						fragA = molAs[i].salts[0];
						fragB = molBs[j].salts[0];
						fragA.setADs();
						fragB.setADs();
						tmp.sim = compare(fragA, fragB, threshold);
					}
					if(tmp.sim >= threshold)
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
						simResults.push_back(tmp);
					}
				}
				else if(molAs[i].salt && !molBs[j].salt)
				{
					tmp.sim = compare(molAs[i].osAtoms, molBs[j].salts[0].fAtoms, threshold);
					if(tmp.sim == 1)
					{
						tmp.sim = (float)0.99999;
					}
					if(tmp.sim >= threshold)
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
						simResults.push_back(tmp);
					}
				}
				else if(!molAs[i].salt && molBs[j].salt)
				{
					tmp.sim = compare(molAs[i].salts[0].fAtoms, molBs[j].osAtoms, threshold);
					if(tmp.sim == 1)
					{
						tmp.sim = (float)0.99999;
					}
					if(tmp.sim >= threshold)
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
						simResults.push_back(tmp);
					}
				}
				else if(molAs[i].salt && molBs[j].salt)
				{
					if(molAs[i].ion && molBs[j].ion)
					{
						tmp.sim = compare(molAs[i].atoms, molBs[j].atoms, threshold);
					}
					else
					{
						tmp.sim = compare(molAs[i].osAtoms, molBs[j].osAtoms, threshold);
					}
					if(tmp.sim == 1)
					{
						tmp.sim = (float)0.99999;
					}
					if(tmp.sim >= threshold)
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
						simResults.push_back(tmp);
					}
				}
			}
			else
			{
				tmp.sim = compare(molAs[i].salts, molBs[j].salts, threshold);
				if(tmp.sim >= threshold)
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
					simResults.push_back(tmp);
				}
			}
		}
	}
	//std::cout<<"done"<<CRLF;
}

inline void sort(std::vector<simResult> &simResults)
{
	std::sort(simResults.begin(), simResults.end(), compSR);
}

#endif
