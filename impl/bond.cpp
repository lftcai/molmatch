#include "bond.h"

bond::bond(std::string oneLine)
{
	//this->bondInfor = oneLine;
	this->aACS = '\0';
	this->aACSE = '\0';
	this->aBCS = '\0';
	this->aBCSE = '\0';
	this->used = false;
	this->bAused = false;
	this->bBused = false;
	bond::setBond(oneLine);
	//this->showBond();
}

bond::bond(int bt, int ba, int bb)
{
	this->bType = bt;
	this->bA = ba;
	this->bB = bb;
}

bond::bond()
{
	//this->bondInfor = "";
	this->bType = 0;
	this->bA = 0;
	this->bB = 0;
	this->aACS = '\0';
	this->aACSE = '\0';
	this->aBCS = '\0';
	this->aBCSE = '\0';
	this->stereo = 0;
	this->used = false;
	this->bAused = false;
	this->bBused = false;
}

void bond::setBond(std::string& bondInfor)
{
	this->bA = atoi(bondInfor.substr(0,3).c_str());
	this->bB = atoi(bondInfor.substr(3,3).c_str());
	this->bType = atoi(bondInfor.substr(6,3).c_str());
	this->stereo = atoi(bondInfor.substr(9,3).c_str());
}

void bond::standardize()
{
	if(this->bB < this->bA)
	{
		char tmpCS, tmpCSE;
		this->bA = this->bA + this->bB;
		this->bB = this->bA - this->bB;
		this->bA = this->bA -this->bB;
		tmpCS = this->aACS;
		tmpCSE = this->aACSE;
		this->aACS = this->aBCS;
		this->aACSE = this->aBCSE;
		this->aBCS = tmpCS;
		this->aBCSE = tmpCSE;
	}
}

bool bond::operator == (const bond& b) const
{
	if(this->bType == b.bType
		&& this->stereo == b.stereo)
	{
		if(this->aACS == b.aACS && this->aACSE == b.aACSE
			&& this->aBCS == b.aBCS && this->aBCSE == b.aBCSE)
			return true;
		else if(this->aACS == b.aBCS && this->aACSE == b.aBCSE
			&& this->aBCS == b.aACS && this->aBCSE == b.aACSE)
			return true;
		else
			return false;
	}
	else
		return false;
}

void bond::showBond()
{
	extern std::string CRLF;
	std::cout
		<<this->bA<<"("<<this->aACS<<this->aACSE<<") "
		<<this->bB<<"("<<this->aBCS<<this->aBCSE<<") "
		<<this->bType
		<<" s:"<<this->stereo
		<<CRLF;
}

std::string bond::toString()
{
	std::string retStr;
	extern std::string CRLF;
	char tmp[15];

	sprintf(tmp,"% 3d",this->bA);
	retStr.append(tmp);
	sprintf(tmp,"% 3d",this->bB);
	retStr.append(tmp);
	sprintf(tmp,"% 3d",this->bType);
	retStr.append(tmp);
	sprintf(tmp,"% 3d",this->stereo);
	retStr.append(tmp);
	retStr.append(CRLF);
	return retStr;
}
