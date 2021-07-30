#include "fragment.h"

fragment::fragment(std::vector<atom> Atoms, std::vector<bond> Bonds)
{
	this->branched = false;
	this->cyclic = false;
	this->linear = false;
	this->atomCount = 0;
	this->bondCount = 0;
	this->fASN = "";
	this->maxASN = 0;
	this->rCount = 0;
	this->charge = 0;
	this->ioIon = false;
	this->oIon = false;
	this->baseRGP = false;
	this->aroma = false;
	this->hetero = false;
	this->sameHA = true;
	this->hAtomC = 0;
	this->aAtomC = 0;
	this->maxRGNO = 0;
	this->used = false;

	this->fAtoms = Atoms;
	this->fBonds = Bonds;
	this->atomCount = fragment::fAtoms.size();
	this->bondCount = fragment::fBonds.size();
	this->setMaxASN();
	this->setFASN();
	//this->showFrag();
}

fragment::fragment()
{
	this->branched = false;
	this->cyclic = false;
	this->linear = false;
	this->atomCount = 0;
	this->bondCount = 0;
	this->fASN = "";
	this->maxASN = 0;
	this->rCount = 0;
	this->charge = 0;
	this->ioIon = false;
	this->oIon = false;
	this->baseRGP = false;
	this->aroma = false;
	this->hetero = false;
	this->sameHA = true;
	this->hAtomC = 0;
	this->aAtomC = 0;
	this->maxRGNO = 0;
	this->used = false;
}

bool fragment::addAtom(atom &a)
{
	if(this->chkASN(a.serNum))
		return false;

	int tmp = a.bCount2;
	int i=0,j=0;
	std::vector<atom>::iterator it = this->fAtoms.begin();
	bool added = false;
	for(i=0;i<this->atomCount;i++)
	{
		for(j=0;j<this->fAtoms[i].bCount;j++)
		{
			if(this->fAtoms[i].bonds[j].bB == a.serNum /*&& !this->chkASN(a.serNum)*/)
			{
				//if(!added)
				//{
					this->fAtoms[i].bCount2 = this->fAtoms[i].bCount2 + 1;
					a.bCount2 = 1;
					this->fAtoms.insert(it+i+1,a);
					added = true;
					this->fBonds.push_back(this->fAtoms[i].bonds[j]);
					break;
				//}
			}
		}
		if(added)
			break;
	}
	if(!added && this->atomCount==0)
	{
		a.bCount2 = 0;
		this->fAtoms.push_back(a);
		added = true;
	}
	this->atomCount = this->fAtoms.size();
	this->bondCount = this->fBonds.size();

	//this->showFrag();system("pause");
	a.bCount2 = tmp;
	return added;
}

bool fragment::addAtom2(atom &a, bool markUsed)
{
	if(this->chkASN(a.serNum))
		return false;

	int tmp = a.bCount2;
	int i=0,j=0;
	std::vector<atom>::iterator it = this->fAtoms.begin();
	bool added = false;
	for(i=0;i<this->atomCount;i++)
	{
		if(this->fAtoms[i].used) continue;
		for(j=0;j<this->fAtoms[i].bCount;j++)
		{
			if(this->fAtoms[i].bonds[j].bB == a.serNum /*&& !this->chkASN(a.serNum)*/)
			{
				//if(!added)
				//{
					this->fAtoms[i].bCount2 = this->fAtoms[i].bCount2 + 1;
					this->fAtoms[i].newAdded = false;
					if(markUsed)
					{
						if(this->fAtoms[i].bCount2 == this->fAtoms[i].bCount)
							this->fAtoms[i].used = true;
						else
							this->fAtoms[i].used = false;
					}
					a.bCount2 = 1;
					if(!a.tAtom)
						a.newAdded = true;
					else
					{
						a.newAdded = false;
						a.used = true;
					}
					this->fAtoms.insert(it+i+1,a);
					added = true;
					this->fBonds.push_back(this->fAtoms[i].bonds[j]);
					break;
				//}
			}
		}
		if(added)
			break;
	}
	if(!added && this->atomCount==0)
	{
		a.bCount2 = 0;
		this->fAtoms.push_back(a);
		added = true;
	}
	this->atomCount = this->fAtoms.size();
	this->bondCount = this->fBonds.size();

	//this->showFrag();system("pause");
	a.bCount2 = tmp;
	return added;
}

void fragment::setMarkush()
{
	this->markushs.clear();
	int i=0;
	int count;
	std::vector<atom> tmpAtoms;
	for(i=0;i<this->fAtoms.size();i++)
	{
		if(this->fAtoms[i].RGPAtom)
		{
			this->maxRGNO = this->fAtoms[i].rgpNO > this->maxRGNO ? this->fAtoms[i].rgpNO : this->maxRGNO;
			tmpAtoms.clear();
			this->nextAtom(this->fAtoms[i],tmpAtoms,count,false);
			this->markushs.push_back(markush(this->fAtoms[i],this->fAtoms[i].rgpNO,tmpAtoms));
		}
	}
}

void fragment::setFrag()
{
	int i=0,j=0;
	/*
	std::vector<atom>::iterator itA = this->fAtoms.begin();
	std::vector<bond>::iterator itB = this->fBonds.begin();
	for(i=0;i<this->atomCount;i++)
	{
		if(this->chkASN(this->fAtoms[i].serNum))
		{
			this->fAtoms.erase(itA+i,&this->fAtoms[i]);
			this->atomCount = this->fAtoms.size();
		}
	}
	for(i=0;i<this->bondCount;i++)
	{
		if(this->chkBAB(this->fBonds[i]))
		{
			this->fBonds.erase(itB+i,&this->fBonds[i]);
			this->bondCount = this->fBonds.size();
		}
	}
	*/
	if(this->atomCount == 1 && this->bondCount == 0)
	{
		this->branched = false;
		this->cyclic = false;
		this->linear = false;
	}
	else if(this->atomCount == 2 && this->bondCount == 1)
	{
		this->branched = false;
		this->cyclic = false;
		this->linear = true;
	}
	else if(this->atomCount >= 3 && this->bondCount >= 2)
	{
		int bondC=0;
		for(i=0;i<this->atomCount;i++)
		{
			if(this->fAtoms[i].bCount - this->fAtoms[i].bCount2 >0)
			{
				bondC = 0;
				for(j=0;j<this->fAtoms[i].bCount;j++)
				{
					if(this->chkASN(this->fAtoms[i].bonds[j].bB))
						bondC++;
					if(!this->chkBAB(this->fAtoms[i].bonds[j]) && this->chkASN(this->fAtoms[i].bonds[j].bB))
						this->fBonds.push_back(this->fAtoms[i].bonds[j]);
				}
				if(bondC > this->fAtoms[i].bCount2)
				{
					this->cyclic = true;
					this->linear = false;
					this->rCount += bondC - this->fAtoms[i].bCount2;
					this->fAtoms[i].bCount2 = bondC;
				}
			}
			if(this->fAtoms[i].bCount2 >= 3)
				this->branched = true;
			if(this->fAtoms[i].charge != 0)
				this->charge += (4 - this->fAtoms[i].charge);
		}
		this->linear = this->cyclic ? false : true;
		this->bondCount = this->fBonds.size();

	}

 	for(i=0;i<this->bondCount;i++)
		this->fBonds[i].standardize();
	//std::cout<<this->branched<<this->cyclic<<this->linear<<"\r\n";
}

void fragment::setHetero()
{
	this->hAtoms.clear();
	this->hAtomC = 0;

	int i=0, j=0;
	for(i=0;i<this->atomCount;i++)
		if(this->fAtoms[i].hAtom)
			this->hAtoms.push_back(this->fAtoms[i]);

	if(hAtoms.size() > 0)
	{
		this->hAtomC = hAtoms.size();
		this->hetero = true;
	}

	if(hAtoms.size() > 1)
	{
		for(i=1;i<hAtoms.size();i++)
		{
			if(hAtoms[0].chemSym != hAtoms[i].chemSym || hAtoms[0].chemSymEx != hAtoms[i].chemSymEx)
			{
				this->sameHA = false;
				break;
			}
		}
	}
}

void fragment::setSalt()
{
	int i=0;
	this->charge = 0;
	for(i=0;i<this->atomCount;i++)
	{
		if(this->fAtoms[i].charge != 0)
			this->charge += (4 - this->fAtoms[i].charge);
	}
	for(i=0;i<this->atomCount;i++)
	{
		if(this->fAtoms[i].chemSym == 'C' && this->fAtoms[i].chemSymEx == ' ')
		{
			this->oIon = true;
			this->ioIon = false;
			break;
		}
		else
			this->ioIon = true;
	}
}

void fragment::setFASN()
{
	for(int i=0;i<this->atomCount;i++)
	{
		char tmp[3];
		//itoa(this->fAtoms[i].serNum,tmp,10);
		sprintf(tmp,"%d",this->fAtoms[i].serNum);
		this->fASN += tmp;
	}
}

bool fragment::chkASN(int asn)
{
	for(int i=0;i<this->fAtoms.size();i++)
	{
		if(asn == this->fAtoms[i].serNum)
			return true;
	}
	return false;
}

bool fragment::chkBAB(bond b)
{

	for(int i=0;i<this->fBonds.size();i++)
	{
		if(
			(this->fBonds[i].bA == b.bA && this->fBonds[i].bB == b.bB)
			|| (this->fBonds[i].bA == b.bB && this->fBonds[i].bB == b.bA)
			)
			return true;
	}
	return false;
}

bool fragment::chkAtomSym(atom &a)
{
	for(int i=0;i<this->fAtoms.size();i++)
	{
		if(a.chemSym == this->fAtoms[i].chemSym && a.chemSymEx == this->fAtoms[i].chemSymEx)
			return true;
	}
	return false;
}

bool fragment::chkAroma()
{
	int count=0;
	int minElectron = 0;
	int maxElectron = 0;
	int i=0;

	for(i=0;i<this->atomCount;i++)
		if(this->fAtoms[i].chemSym != 'C' && this->fAtoms[i].chemSym != 'R'  && this->fAtoms[i].chemSymEx ==' ')
			count++;

	minElectron = this->atomCount;
	maxElectron = minElectron + count;

	for(i=minElectron;i<=maxElectron;i++)
		if((i-2)%4 == 0)  //4n+2 rule
		{
			this->aroma = true;
			break;
		}

	if(this->aroma)
	{
		for(i=0;i<this->atomCount;i++)
			this->fAtoms[i].aroma = true;
		for(i=0;i<this->bondCount;i++)
		{
			if(this->fBonds[i].bType == 1 || this->fBonds[i].bType == 2)
				this->fBonds[i].bType = 4;
			else
				this->fBonds[i].bType = 5;
		}
	}
	return this->aroma;
}

bool fragment::chkMarkush(markush &m)
{
	if(this->markushs.size() == 0) return true;

	int i=0;
	bool exist = false;
	for(i=0;i<this->markushs.size();i++)
	{
		if(this->markushs[i] == m)
		{
			exist = true;
			break;
		}
	}
	return exist;
}

void fragment::setMaxASN()
{
	this->maxASN = 0;
	for(int i=0;i<this->atomCount;i++)
	{
		this->maxASN = (int)this->fAtoms[i].serNum > this->maxASN ? this->fAtoms[i].serNum : this->maxASN;
	}
	//std::cout<<"maxASN: "<<this->maxASN<<"\r\n";
}

void fragment::setABTs()
{
	int i=0, j=0;
	for(i=0;i<this->atomCount;i++)
	{
		this->fAtoms[i].bondTs.clear();
		for(j=0;j<this->bondCount;j++)
			if(this->fAtoms[i].serNum == this->fBonds[j].bA || this->fAtoms[i].serNum == this->fBonds[j].bB)
				this->fAtoms[i].bondTs.push_back(this->fBonds[j].bType);
	}
	for(i=0; i<this->atomCount; i++)
	{
		std::sort(this->fAtoms[i].bondTs.begin(), this->fAtoms[i].bondTs.end());
	}
}

//check whether atom already exist in a atom vector
inline bool chkAtom(std::vector<atom> &atoms, atom &a)
{
	int i=0;
	for(i=0;i<atoms.size();i++)
		if(atoms[i].serNum == a.serNum)
			return true;
	return false;
}

void fragment::setADs()
{
	int i=0, j=0, k=0, l=0, m=0;
	std::vector<atom> atoms;
	for(i=0;i<this->atomCount;i++)
		this->fAtoms[i].dists.clear();

	for(i=0;i<this->atomCount;i++)
	{
		if(this->fAtoms[i].bCount2 == 1)
		{
			this->fAtoms[i].dists.push_back(0);
			atoms.push_back(this->fAtoms[i]);
			for(j=0;j<atoms.size();j++)
			{
				for(k=0;k<atoms[j].bCount;k++)
				{
					for(l=0;l<this->atomCount;l++)
					{
						if(atoms[j].bonds[k].bB == this->fAtoms[l].serNum && !chkAtom(atoms,this->fAtoms[l]))
						{
							this->fAtoms[l].dists.push_back(atoms[j].dists[atoms[j].dists.size()-1] + 1);
							atoms.push_back(this->fAtoms[l]);
							break;
						}
					}
				}
			}
			atoms.clear();
		}
	}
	for(i=0; i<this->atomCount; i++)
	{
		std::sort(this->fAtoms[i].dists.begin(), this->fAtoms[i].dists.end());
	}
}

void fragment::findTA()
{
	int i=0;
	for(i=0;i<this->atomCount;i++)
	{
		if(this->fAtoms[i].bCount2 == 1)
			this->fAtoms[i].tAtom = true;
		else
			this->fAtoms[i].tAtom = false;
	}
}

void fragment::findAA()
{
	int i=0;
	this->aAtomC = 0;
	for(i=0;i<this->atomCount;i++)
		if(this->fAtoms[i].aroma)
			this->aAtomC++;
}

void fragment::getASN(std::vector<int>& asns) const
{
	int i=0;
	for(i=0; i<this->atomCount; i++)
	{
		asns.push_back(this->fAtoms[i].serNum);
	}
}


bool fragment::operator == (fragment& frag)
{
	bool result;
	if(this->atomCount == frag.atomCount)
	{
		int i=0, j=0;
		bool find = false;
		for(i=0;i<this->atomCount;i++)
		{
			find = false;
			for(j=0;j<frag.atomCount;j++)
			{
				if(frag.fAtoms[j].used) continue;
				if(this->fAtoms[i] > frag.fAtoms[j]) break; //fragment atoms should have been sorted
				else if(this->fAtoms[i] == frag.fAtoms[j])
				{
					find = true;
					frag.fAtoms[j].used = true;
					break;
				}
			}
			if(!find)
				break;
		}
		if(find)
			result = true;
		else
			result = false;

		for(i=0;i<frag.atomCount;i++)
			frag.fAtoms[i].used = false;
	}
	else
		result = false;

	return result;
}


/*
bool fragment::operator == (fragment& frag) const
{
	if(this->atomCount == frag.atomCount)
	{
		int i=0;
		for(i=0; i<this->atomCount; i++)
		{
			if(!(this->fAtoms[i] ==  frag.fAtoms[i]))
			{
				return false;
			}
		}
		return true;
	}
	else
		return false;
}
*/
bool fragment::getAtom(atom &a, std::vector<atom>& returnAtoms, bool oneTime)
{
	int i=0;
	int size = returnAtoms.size();
	for(i=0;i<this->atomCount;i++)
	{
		if(oneTime && this->fAtoms[i].used) continue;
		//if(this->fAtoms[i] == a)
		if(this->fAtoms[i].chemSym == a.chemSym && this->fAtoms[i].chemSymEx == a.chemSymEx && chkABTs(this->fAtoms[i].bondTs, a.bondTs))
		{
			returnAtoms.push_back(this->fAtoms[i]);
			if(oneTime)
				this->fAtoms[i].used = true;
		}
	}

	if(returnAtoms.size() > size)
		return true;
	else
		return false;
}

bool fragment::getAtom(int serNum, atom& returnAtom, bool oneTime)
{
	int i=0;
	for(i=0;i<this->atomCount;i++)
	{
		if(oneTime && this->fAtoms[i].used) continue;
		if(this->fAtoms[i].serNum == serNum)
		{
			returnAtom = this->fAtoms[i];
			if(oneTime)
				this->fAtoms[i].used = true;
			return true;
		}
	}
	return false;
}

bool fragment::nextAtom(atom &a, std::vector<atom>& retAtoms, int &count, bool oneTime)
{
	int i=0;
	int size = retAtoms.size();
	count = 0;
	std::vector<atom> atoms;
	atom tmpAtom;
	for(i=0;i<a.bCount;i++)
	{
		if(this->getAtom(a.bonds[i].bB, tmpAtom, oneTime))
		{
			atoms.push_back(tmpAtom);
			tmpAtom.clear();
			count++;
		}
	}
	
	for(i=0;i<atoms.size();i++) retAtoms.push_back(atoms[i]);

	if(retAtoms.size() > size)
		return true;
	else
		return false;
}

void fragment::getNMAtoms(std::vector<atom>& retAtoms)
{
	int i=0, j=0;
	bool exist = false;
	for(i=0;i<this->fAtoms.size();i++)
	{
		exist = false;
		for(j=0;j<this->markushs.size();j++)
		{
			if(this->markushs[j].chkASN(this->fAtoms[i].serNum))
			{
				exist = true;
				break;
			}
		}
		if(!exist)
			retAtoms.push_back(this->fAtoms[i]);
	}
}

void fragment::getNMBonds(std::vector<bond>& retBonds)
{
	int i=0, j=0;
	bool exist = false;
	bond tmpB;
	for(i=0;i<this->bondCount;i++)
	{
		exist = false;
		for(j=0;j<this->markushs.size();j++)
		{
			if(this->markushs[j].pointAtom.serNum == this->fBonds[i].bA
				|| this->markushs[j].pointAtom.serNum == this->fBonds[i].bB)
			{
				exist = true;
				break;
			}
		}
		if(!exist)
			retBonds.push_back(this->fBonds[i]);
	}
	for(i=0;i<this->markushs.size();i++)
	{
		for(j=0;j<this->markushs[i].atoms.size();j++)
		{
			if(this->markushs[i].atoms[j].getBond(this->markushs[i].pointAtom.serNum, tmpB))
				retBonds.push_back(tmpB);
		}
	}
}

void fragment::sortAtoms()
{
	int i=0;
	int start = 0;
	atom tmpAtom;
	std::vector<atom> atoms;

	this->setMaxASN();
	start = this->maxASN - this->atomCount + 1;
	for(i=0;i<this->atomCount;i++)
	{
		this->getAtom(start++,tmpAtom,false);
		atoms.push_back(tmpAtom);
	}
	this->fAtoms = atoms;
}

void fragment::sortAtoms2()
{
/*
	int i=0, j=0;
	atom tmpAtom;
	for(i=0;i<this->atomCount;i++)
	{
		for(j=i+1;j<this->atomCount;j++)
		{
			if(this->fAtoms[j] > this->fAtoms[i])
			{
				tmpAtom = this->fAtoms[i];
				this->fAtoms[i] = this->fAtoms[j];
				this->fAtoms[j] = tmpAtom;
			}
		}
	}
*/
	std::sort(this->fAtoms.begin(), this->fAtoms.end(), largerAtom);

}

void fragment::showFrag()
{
	int i=0;
	extern std::string CRLF;
	for(i=0;i<this->atomCount;i++)
		std::cout<<this->fAtoms[i].serNum<<"("<<this->fAtoms[i].bCount2<<")";
	std::cout<<" LCB:"<<this->linear<<this->cyclic<<this->branched<<CRLF;
	for(i=0;i<this->atomCount;i++)
		this->fAtoms[i].showAtom();
	
	for(i=0;i<this->bondCount;i++)
	{
		this->fBonds[i].showBond();
	}
	//system("pause");
}

void fragment::showMol(std::string& name, std::string& molInfor, int baseCount)
{
	this->sortAtoms();

	extern std::string CRLF;
	int i=0;
	std::cout<<name<<CRLF;
	std::cout<<"\tMolengineFragmentSave"<<CRLF;
	std::cout<<CRLF;
	std::cout<<std::setw(3)<<this->atomCount;
	std::cout<<std::setw(3)<<this->bondCount;
	std::cout<<molInfor.substr(6)<<CRLF;
	for(i=0;i<this->atomCount;i++)
	{
		std::cout<<std::setiosflags(std::ios::fixed);
		std::cout<<std::setw(10)<<std::setprecision(4)<<this->fAtoms[i].x;
		std::cout<<std::setw(10)<<std::setprecision(4)<<this->fAtoms[i].y;
		std::cout<<std::setw(10)<<std::setprecision(4)<<this->fAtoms[i].z;
		std::cout<<' ';
		std::cout<<this->fAtoms[i].chemSym<<this->fAtoms[i].chemSymEx<<' ';
		std::cout<<std::setw(2)<<this->fAtoms[i].massD;
		std::cout<<std::setw(3)<<this->fAtoms[i].charge;
		//std::cout<<std::setw(3)<<this->fAtoms[i].stereo;
		std::cout<<CRLF;
	}
	for(i=0;i<this->bondCount;i++)
	{
		std::cout<<std::setw(3)<<this->fBonds[i].bA - baseCount;
		std::cout<<std::setw(3)<<this->fBonds[i].bB - baseCount;
		std::cout<<std::setw(3)<<this->fBonds[i].bType;
		std::cout<<std::setw(3)<<this->fBonds[i].stereo;
		std::cout<<CRLF;
	}
}

void fragment::showHetero()
{
	extern std::string CRLF;
	int i=0;
	if(this->sameHA)
		std::cout<<"Hetero atom count: "<<this->hAtomC<<", all same"<<CRLF;
	else
		std::cout<<"Hetero atom count: "<<this->hAtomC<<", different"<<CRLF;
	for(i=0;i<this->hAtomC;i++)
		this->hAtoms[i].showAtom();
}

bool fragment::eraseAtom(atom &a)
{
	int i=0;
	bool deleted = false;
	std::vector<atom>::iterator it;
	for(i=0;i<this->atomCount;i++)
	{
		if(a.serNum == this->fAtoms[i].serNum)
		{
			//this->fAtoms.erase(&this->fAtoms[i]);
			it = this->fAtoms.begin();
			this->fAtoms.erase(it + i);
			this->eraseBond(a.serNum);
			deleted = true;
			break;
		}
	}

	if(deleted)
		this->atomCount = this->fAtoms.size();

	return deleted;
}

bool fragment::eraseBond(int serNum)
{
	int i=0;
	bool deleted = false;
	std::vector<bond>::iterator it;
	for(i=0;i<this->fBonds.size();i++)
	{
		if(this->fBonds[i].bA == serNum || this->fBonds[i].bB == serNum)
		{
			//this->fBonds.erase(&this->fBonds[i]);
			it = this->fBonds.begin();
			this->fBonds.erase(it + i);
			i -= 1;
			deleted = true;
		}
	}

	if(deleted)
	{
		this->bondCount = this->fBonds.size();
		for(i=0;i<this->fAtoms.size();i++)
			this->fAtoms[i].eraseBond(serNum);
	}

	return deleted;
}

bool fragment::reviseAtom(atom &targetAtom)
{
	int i=0;
	bool revised = false;
	for(i=0;i<this->atomCount;i++)
	{
		if(this->fAtoms[i].serNum == targetAtom.serNum)
		{
			this->fAtoms[i] = targetAtom;
			break;
			revised = true;
		}
	}
	return revised;
}

void fragment::reviseASN(int asn)
{
	if(this->atomCount == 0) return;

	int i=0, j=0;
	int m = 0;
	int tmp = this->fAtoms[0].serNum;
	for(i=0;i<this->fAtoms.size();i++)
		tmp = this->fAtoms[i].serNum < tmp ? this->fAtoms[i].serNum : tmp;

	m = asn - tmp;
	if(m != 0)
	{
		for(i=0;i<this->fAtoms.size();i++)
		{
			tmp = this->fAtoms[i].serNum + m;
			this->fAtoms[i].serNum = tmp;
			for(j=0;j<this->fAtoms[i].bCount;j++)
			{
				this->fAtoms[i].bonds[j].bA = tmp;
				this->fAtoms[i].bonds[j].bB = this->fAtoms[i].bonds[j].bB + m;
			}
		}
		for(i=0;i<this->fBonds.size();i++)
		{
			this->fBonds[i].bA = this->fBonds[i].bA + m;
			this->fBonds[i].bB = this->fBonds[i].bB + m;
		}
	}
}

void fragment::formateASN()
{
	int i=0, j=0, k=0;
	int tmp=0;
	for(i=0;i<this->atomCount;i++)
	{
		tmp = this->fAtoms[i].serNum;
		this->fAtoms[i].serNum = i+1;
		for(j=0;j<this->fAtoms[i].bCount;j++)
		{
			this->fAtoms[i].bonds[j].bA = i+1;
		}
		for(j=0;j<this->atomCount;j++)
		{
			if(j == i) continue;
			for(k=0;k<this->fAtoms[j].bCount;k++)
			{
				if(this->fAtoms[j].bonds[k].bB == tmp && !this->fAtoms[j].bonds[k].used)
				{
					this->fAtoms[j].bonds[k].bB = i+1;
					this->fAtoms[j].bonds[k].used = true;
				}
			}
		}
		for(j=0;j<this->bondCount;j++)
		{
			if(this->fBonds[j].bA == tmp && !this->fBonds[j].bAused)
			{
				this->fBonds[j].bA = i+1;
				this->fBonds[j].bAused = true;
			}
			if(this->fBonds[j].bB == tmp && !this->fBonds[j].bBused)
			{
				this->fBonds[j].bB = i+1;
				this->fBonds[j].bBused = true;
			}
		}
	}
	this->resetUsed2();
}

void fragment::clear()
{
	this->fAtoms.clear();
	this->fBonds.clear();
	this->hAtoms.clear();
	this->markushs.clear();

	this->branched = false;
	this->cyclic = false;
	this->linear = false;
	this->atomCount = 0;
	this->bondCount = 0;
	this->fASN = "";
	this->maxASN = 0;
	this->rCount = 0;
	this->charge = 0;
	this->ioIon = false;
	this->oIon = false;
	this->baseRGP = false;
	this->aroma = false;
	this->hetero = false;
	this->sameHA = true;
	this->hAtomC = 0;
	this->aAtomC = 0;
	this->maxRGNO = 0;
	this->used = false;
}

void fragment::clearABTs()
{
	int i=0;
	for(i=0;i<this->atomCount;i++)
		this->fAtoms[i].bondTs.clear();
}

void fragment::clearADs()
{
	int i=0;
	for(i=0;i<this->atomCount;i++)
		this->fAtoms[i].dists.clear();
}

void fragment::resetUsed()
{
	int i=0;
	for(i=0;i<this->atomCount;i++)
		this->fAtoms[i].used = false;
}

void fragment::resetUsed2()
{
	int i=0, j=0;
	for(i=0;i<this->atomCount;i++)
	{
		for(j=0;j<this->fAtoms[i].bCount;j++)
		{
			this->fAtoms[i].bonds[j].bAused = false;
			this->fAtoms[i].bonds[j].bBused = false;
			this->fAtoms[i].bonds[j].used = false;
		}
	}
	for(i=0;i<this->bondCount;i++)
	{
		this->fBonds[i].used = false;
		this->fBonds[i].bAused = false;
		this->fBonds[i].bBused = false;
	}
}

void fragment::resetNewAdded()
{
	int i=0;
	for(i=0;i<this->atomCount;i++)
		this->fAtoms[i].newAdded = false;
}
