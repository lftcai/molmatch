#include "atom.h"

atom::atom(std::string oneLine, int serialNum, bool H)
{
	this->bCount = 0;
	this->bCount2 = 0;
	//this->atomInfor = oneLine;
	this->serNum = serialNum;
	this->star = true;
	this->atomH = H;
	this->nonC = false;
	this->hAtom = false;
	this->aroma = false;
	this->tAtom = false;
	this->RGPAtom = false;
	this->rgpNO = 0;
	this->used = false;
	this->newAdded = false;
	this->massD = 0;
	this->charge = 0;
	this->stereo = 0;
	this->chemSym = '\0';
	this->chemSymEx = '\0';
	this->atomNum = 0;
	this->setAtom(oneLine);
	//this->showAtom();
}
atom::atom()
{
	this->bCount = 0;
	this->bCount2 = 0;
	this->star = true;
	this->atomH = false;
	this->nonC = false;
	this->aroma = false;
	this->hAtom = false;
	this->tAtom = false;
	this->RGPAtom = false;
	this->rgpNO = 0;
	this->used = false;
	this->newAdded = false;
	//this->atomInfor = "";
	this->x = 0;
	this->y = 0;
	this->z = 0;
	this->massD = 0;
	this->charge = 0;
	this->stereo = 0;
	this->chemSym = '\0';
	this->chemSymEx = '\0';
	this->atomNum = 0;
}

void atom::setAtom(std::string& atomInfor)
{
	this->x = atof(atomInfor.substr(0,10).c_str());
	this->y = atof(atomInfor.substr(10,10).c_str());
	this->z = atof(atomInfor.substr(20,10).c_str());
	this->chemSym = atomInfor[31];
	this->chemSymEx = atomInfor[32];
	this->massD = atoi(atomInfor.substr(34,2).c_str());
	this->charge = atoi(atomInfor.substr(36,3).c_str());
	//this->stereo = atoi(atomInfor.substr(39,3).c_str());

	if(chemSym == '*')
	{
		this->star = true;
		this->atomNum = 999;
	}
	else if(chemSym == 'C' && chemSymEx == ' ')
		this->atomNum = 6;
	else if(chemSym == 'O' && chemSymEx == ' ')
		this->atomNum = 8;
	else if(chemSym == 'N' && chemSymEx == ' ')
		this->atomNum = 7;
	else if(chemSym == 'C' && chemSymEx == 'l')
		this->atomNum = 17;
	else if(chemSym == 'P' && chemSymEx == ' ')
		this->atomNum = 15;
	else if(chemSym == 'S' && chemSymEx == ' ')
		this->atomNum = 16;
	else
		this->atomNum = 102;

	if(this->chemSymEx != ' ' || this->chemSym != 'C')
	{
		this->nonC = true;
		this->hAtom = true;
	}

	if(this->chemSym == 'R' && this->chemSymEx == '#')
		this->RGPAtom = true;
	else if(this->chemSym == 'R' && this->chemSymEx == ' ')
	{
		this->RGPAtom = true;
		this->rgpNO = 1;
	}
}

void atom::setBonds(std::vector<bond> b)
{
	this->bCount = b.size();
	this->bonds = b;

	if(this->chemSymEx == ' ' && this->chemSym == 'C')
		for(int i=0;i<this->bCount;i++)
			if(this->bonds[i].bType != 1)
			{
				this->nonC = true;
				break;
			}
	if(this->bCount == 1)
		this->tAtom = true;
}

bool atom::operator == (const atom& a) const
{
	if(	this->atomNum != a.atomNum
		||this->nonC != a.nonC
		||this->hAtom != a.hAtom
		||this->aroma != a.aroma
		||this->stereo != a.stereo
		)
		return false;

	bool result = false;
	if(this->chemSym == a.chemSym
		&& this->chemSymEx == a.chemSymEx
		&& this->massD == a.massD)
	{
		if(this->bondTs.size() == a.bondTs.size()
			&& this->bondTs.size() !=0
			&& a.bondTs.size() != 0)
		{
			if(equalV(this->bondTs, a.bondTs))
				result = true;
			else
				result = false;
		}
		else if(this->bondTs.size() == 0 || a.bondTs.size() == 0)
			result = true;
		else
			result = false;
	}
	else
		result = false;

	if(result == true && this->dists.size() > 0 && a.dists.size() > 0)
	{
		if(this->dists.size() == a.dists.size())
		{
			if(equalV(this->dists, a.dists))
				result = true;
			else
				result = false;
		}
		else
			result = false;
	}

	return result;
}

bool atom::operator > (const atom& a) const
{
	//others > Cl > S > P > O > N > C

	if(this->atomNum > a.atomNum)
		return true;
	else if(this->atomNum == a.atomNum)
	{
		if(this->aroma && !a.aroma)
			return true;
		else if(!this->aroma && a.aroma)
			return false;

		if(this->aroma == a.aroma)
		{
			if(this->bCount2 < a.bCount2)   //!!!! bCount or bCount2
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

bool atom::chkAround(atom &a, std::vector<int> &retASNs)
{
	int i=0, j=0;
	int size = retASNs.size();
	for(i=0;i<this->bCount;i++)
	{
		for(j=0;j<a.bCount;j++)
		{
			//if(a.bonds[j].used) continue;
			if(this->bonds[i] == a.bonds[j])
			{
				//a.bonds[j].used = true;
				retASNs.push_back(this->bonds[i].bB);
				break;
			}
		}
	}
	a.resetUsed();
	if(retASNs.size() > size)
		return true;
	else
		return false;
}

bool atom::eraseBond(int serNum)
{
	int i=0;
	bool deleted = false;
	std::vector<bond>::iterator it;
	for(i=0;i<this->bonds.size();i++)
	{
		if(this->bonds[i].bB == serNum)
		{
			//this->bonds.erase(&this->bonds[i]);
			it = this->bonds.begin();
			this->bonds.erase(it + i);
			deleted = true;
			break;
		}
	}

	if(deleted)
	{
		this->bCount = this->bonds.size();
		this->bCount2 = (this->bCount2 -1 >= 0 ? this->bCount2 -1 : 0);
		if(this->bCount == 1) this->tAtom = true;
	}

	return deleted;
}

bool atom::reviseBond(int serNum, int newSerNum, bool mark)
{
	int i=0;
	bool revised = false;
	for(i=0;i<this->bonds.size();i++)
	{
		if(this->bonds[i].bB == serNum && !this->bonds[i].used)
		{
			this->bonds[i].bB = newSerNum;
			if(mark) this->bonds[i].used = true;
			revised = true;
			break;
		}
	}
	return revised;
}

bool atom::getBond(int serNum, bond& retBond)
{
	int i=0;
	bool exist = false;
	for(i=0;i<this->bonds.size();i++)
	{
		if(this->bonds[i].bB == serNum)
		{
			retBond = this->bonds[i];
			exist = true;
			break;
		}
	}
	return exist;
}

void atom::showAtom()
{
	int i=0;
	extern std::string CRLF;
	std::cout
		<<this->serNum<<" "
		<<this->chemSym<<this->chemSymEx<<" "
		<<this->x<<" "<<this->y<<" "<<this->z
		//<<" mD:"<<this->massD
		<<" c:"<<this->charge
		<<" s:"<<this->stereo
		<<" a:"<<this->aroma
		<<" bC:"<<bCount<<" bC2:"<<bCount2
		<<" tA:"<<this->tAtom
		<<" rNO:"<<this->rgpNO
		<<" bT:";
	for(i=0;i<this->bondTs.size();i++)
	{
		std::cout<<this->bondTs[i];
		if(i != this->bondTs.size()-1)
			std::cout<<"-";
	}
	std::cout<<" ds:";
	for(i=0;i<this->dists.size();i++)
	{
		std::cout<<this->dists[i];
		if(i != this->dists.size()-1)
			std::cout<<"-";
	}
	std::cout<<" u:";std::cout<<this->used;
	std::cout<<" na:";std::cout<<this->newAdded;
	//std::cout<<" h:";std::cout<<this->hAtom;
	//std::cout<<" nC:";std::cout<<this->nonC;
	std::cout<<CRLF;
	//for(i=0;i<this->bCount;i++) this->bonds[i].showBond();
}

void atom::clear()
{
	this->bCount = 0;
	this->bCount2 = 0;
	this->star = false;	
	this->atomH = false;
	this->nonC = false;
	this->aroma = false;
	this->hAtom = false;
	this->tAtom = false;
	this->RGPAtom = false;
	this->rgpNO = 0;
	this->used = false;
	this->newAdded = false;
	//this->atomInfor = "";
	this->x = 0;
	this->y = 0;
	this->z = 0;
	this->massD = 0;
	this->charge = 0;
	this->stereo = 0;
	this->chemSym = '\0';
	this->chemSymEx = '\0';
	this->atomNum = 0;

	this->bonds.clear();
	this->bondTs.clear();
	this->dists.clear();
}

void atom::resetUsed()
{
	int i=0;
	for(i=0;i<this->bCount;i++)
		this->bonds[i].used = false;
}

std::string atom::toString()
{
	std::string retStr;
	extern std::string CRLF;
	char tmp[45];

	sprintf(tmp,"% 10.4f",this->x);
	retStr.append(tmp);
	sprintf(tmp,"% 10.4f",this->y);
	retStr.append(tmp);
	sprintf(tmp,"% 10.4f",this->z);
	retStr.append(tmp);
	retStr.append(" ");
	tmp[0]  = this->chemSym;
	tmp[1]  = this->chemSymEx;
	tmp[2]  = '\0';
	retStr.append(tmp);
	sprintf(tmp,"% 2d",this->massD);
	retStr.append(tmp);
	sprintf(tmp,"% 3d",this->charge);
	retStr.append(tmp);
	sprintf(tmp,"% 3d",this->stereo);
	//retStr.append(tmp);
	retStr.append(CRLF);
	return retStr;
}
