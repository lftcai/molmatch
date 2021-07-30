#include "markush.h"

rGroupMDL::rGroupMDL()
{
	this->mdlRgroup = "";
	this->firstNO = 0;
	this->serNum = 0;
	this->groupNO = 0;
}

rGroupMDL::rGroupMDL(std::string oneLine)
{
	this->mdlRgroup = oneLine;
	this->setRGroup();
}

void rGroupMDL::setRGroup()
{
	int l = mdlRGroupMark.length();
	this->firstNO = atoi(this->mdlRgroup.substr(l,3).c_str());
	this->serNum = atoi(this->mdlRgroup.substr(l+3,4).c_str());
	this->groupNO = atoi(this->mdlRgroup.substr(l+7,4).c_str());
}

void rGroupMDL::showRGP()
{
	std::cout<<this->mdlRgroup<<CRLF;
	std::cout<<"f: "<<this->firstNO;
	std::cout<<" s: "<<this->serNum;
	std::cout<<" g: "<<this->groupNO<<CRLF;
}

std::string rGroupMDL::toString()
{
	std::string retStr;
	extern std::string CRLF;
	char tmp[20];

	retStr.append(mdlRGroupMark);
	sprintf(tmp,"% 3d",this->firstNO);
	retStr.append(tmp);
	sprintf(tmp,"% 3d",this->serNum);
	retStr.append(tmp);
	sprintf(tmp,"% 3d",this->groupNO);
	retStr.append(tmp);
	retStr.append(CRLF);
	return retStr;
}

rGroupCtab::rGroupCtab(std::string& oneLine)
{
	this->ctabLine = oneLine;
	this->setRGroup();
}

void rGroupCtab::setRGroup()
{
	int l = CtabRLOGMark.length();
	this->firstNO = atoi(this->ctabLine.substr(l,3).c_str());
	this->secondNO = atoi(this->ctabLine.substr(l+3,4).c_str());
}

void rGroupCtab::showRGP()
{
	std::cout<<this->ctabLine<<CRLF;
	std::cout<<"f: "<<this->firstNO;
	std::cout<<" s: "<<this->secondNO<<CRLF;
}

markush::markush()
{
	this->groupNO = 0;
	this->atoms.clear();
}

markush::markush(atom &pA, int &gNO, std::vector<atom>& as)
{
	this->pointAtom = pA;
	this->groupNO = gNO;
	this->atoms = as;
	int i=0;
	for(i=0;i<this->pointAtom.bCount;i++)
		this->bonds.push_back(this->pointAtom.bonds[i]);
	for(i=0;i<this->bonds.size();i++)
		this->bonds[i].standardize();
}

bool markush::operator == (const markush& m) const
{
	if(this->groupNO == m.groupNO)
		return true;
	else
		return false;
}

bool markush::getBond(int &serNum, bond& retBond)
{
	int i=0;
	bool exist = false;
	for(i=0;i<this->pointAtom.bonds.size();i++)
	{
		if(this->pointAtom.bonds[i].bB == serNum)
		{
			retBond = this->pointAtom.bonds[i];
			exist = true;
			break;
		}
	}
	return exist;
}

void markush::getBonds(std::vector<bond>& retBonds)
{
	int i=0;
	for(i=0;i<this->bonds.size();i++)
		retBonds.push_back(this->bonds[i]);
}

void markush::getAtoms(std::vector<atom>& retAtoms)
{
	int i=0;
	for(i=0;i<this->atoms.size();i++)
		retAtoms.push_back(this->atoms[i]);
}

bool markush::chkASN(int serNum)
{
	bool exist = false;
	int i=0;
	for(i=0;i<this->atoms.size();i++)
	{
		if(this->atoms[i].serNum == serNum)
		{
			exist = true;
			break;
		}
	}
	if(this->pointAtom.serNum != serNum && !exist)
		return false;
	else 
		return true;
}

void markush::eraseAB()
{
	int i=0;
	for(i=0;i<this->atoms.size();i++)
		this->atoms[i].eraseBond(this->pointAtom.serNum);
}

void markush::show()
{
	std::cout<<"Group: "<<this->groupNO<<" PointAtom:"<<CRLF;
	this->pointAtom.showAtom();
	std::cout<<"Atoms arounded:"<<CRLF;
	int i=0;
	for(i=0;i<this->atoms.size();i++)
		this->atoms[i].showAtom();
}
