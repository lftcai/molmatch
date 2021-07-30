#include "molecule.h"

std::string molEnd = "M  END";
std::string CTABmark = "CTAB";
std::string SDFdelim = "$$$$";
std::string CRLF = "\n";//in string, c++ is \n, java is \r\n
std::string molPropMark = ">  <";
std::string molPropMarkB = ">";
std::string mdlRGroupMark = "M  RGP";
std::string CtabRLOGMark = "M  LOG";

bool chkMol(std::string inputTmp)
{
	if(inputTmp.size() > 10)
	{
		if(inputTmp.find(molEnd) != std::string::npos
			&& inputTmp.find("$"+CTABmark) == std::string::npos)
			return true;
		else
			return false;
	}
	else
		return false;
}

bool chkCtab(std::string inputTmp)
{
	if(inputTmp.size() > 10)
	{
		if(inputTmp.find("$"+CTABmark) != std::string::npos
			&& inputTmp.find(molEnd) != std::string::npos)
			return true;
		else
			return false;
	}
	else
		return false;
}

void CtabToMol(CTable& ctab, std::string& retMol)
{
	int i=0, j=0;
	int start = 0;
	char arr[4];
	rGroupMDL MDLrgp;
	std::vector<molecule> tmpMols;
	std::vector<atom> atoms;
	std::vector<bond> bonds;
	std::vector<rGroupMDL> MDLrgps;

	tmpMols.push_back(molecule(ctab.baseFrag));
	tmpMols[0].setABC();
	tmpMols[0].setCCP();

	if(tmpMols[0].CtabRGP.size() >= 0 && tmpMols[0].MDLRGP.size() == 1)
	{
		tmpMols[0].MDLRGP.clear();
		for(i=0;i<tmpMols[0].CtabRGP.size();i++)
		{
			for(j=start;j<tmpMols[0].atomCount;j++)
			{
				if(tmpMols[0].atoms[j].RGPAtom)
				{
					MDLrgp.firstNO = tmpMols[0].CtabRGP[i].firstNO;
					MDLrgp.serNum = tmpMols[0].atoms[j].serNum;
					MDLrgp.groupNO = tmpMols[0].CtabRGP[i].secondNO;
					tmpMols[0].MDLRGP.push_back(MDLrgp);
					start = j + 1;
					break;
				}
			}
		}
	}

	for(i=0;i<ctab.rgps.size();i++)
	{
		for(j=0;j<ctab.rgps[i].frags.size();j++)
		{
			tmpMols.push_back(molecule(ctab.rgps[i].frags[j]));
			tmpMols[tmpMols.size()-1].rgroup = true;
			tmpMols[tmpMols.size()-1].setABC();
			tmpMols[tmpMols.size()-1].setCCP();
			tmpMols[tmpMols.size()-1].addRAtom();
		}
	}

	int order = 0;
	order = tmpMols[0].atomCount;
	for(i=1;i<tmpMols.size();i++)
	{
		changeB(tmpMols[i].bonds, order);
		tmpMols[i].MDLRGP[0].serNum = order + 1;
		order += tmpMols[i].atomCount;
	}

	for(i=0;i<tmpMols.size();i++)
	{
		for(j=0;j<tmpMols[i].atomCount;j++)
			atoms.push_back(tmpMols[i].atoms[j]);
		for(j=0;j<tmpMols[i].bondCount;j++)
			bonds.push_back(tmpMols[i].bonds[j]);
		for(j=0;j<tmpMols[i].MDLRGP.size();j++)
			MDLrgps.push_back(tmpMols[i].MDLRGP[j]);
	}
	
	retMol.append(ctab.head);
	sprintf(arr,"% 3d",atoms.size());
	retMol.append(arr);
	sprintf(arr,"% 3d",bonds.size());
	retMol.append(arr);
	retMol.append(tmpMols[0].molInfor.substr(6));
	retMol.append(CRLF);
	for(i=0;i<atoms.size();i++)
		retMol.append(atoms[i].toString());
	for(i=0;i<bonds.size();i++)
		retMol.append(bonds[i].toString());
	for(i=0;i<MDLrgps.size();i++)
		retMol.append(MDLrgps[i].toString());
	retMol.append(molEnd);
	retMol.append(CRLF);
}

long molecule::molCount = 0;
long molecule::saltCount = 0;
long molecule::errMolCount = 0;

molecule::molecule()
{
	this->molFileVer = "";
	this->molName = "";
	this->molBuilder = "";
	this->molComment = "";
	this->molInfor = "";
	this->molFileVer = "";
	this->molBody = "";
	this->molMap = "";
	this->molExtInfor = "";
	this->molProps = "";
	this->propCount = 0;
	this->atomCount = 0;
	this->bondCount = 0;
	this->fragCount = 0;
	this->sIonCount = 0;
	this->arCount = 0;
	this->noMolStru = false;
	this->molExt = false;
	this->rgroup = false;
	this->ctable = false;
	this->ion = false;
	this->salt = false;
	this->aroma = false;
	this->polymer = false;
	//this->saltFrag = 0;
	this->charge = 0;
	this->noMolStru = false;
	this->molID = 0;
}

molecule::molecule(std::string tmp)
{
	//std::cout<<"constructing molecule...";

	this->molFileVer = "";
	this->molName = "";
	this->molBuilder = "";
	this->molComment = "";
	this->molInfor = "";
	this->molFileVer = "";
	this->molBody = "";
	this->molMap = "";
	this->molExtInfor = "";
	this->molProps = "";
	this->propCount = 0;
	this->atomCount = 0;
	this->bondCount = 0;
	this->fragCount = 0;
	this->sIonCount = 0;
	this->arCount = 0;
	this->noMolStru = false;
	this->molExt = false;
	this->rgroup = false;
	this->ctable = false;
	this->ion = false;
	this->salt = false;
	this->aroma = false;
	this->polymer = false;
	//this->saltFrag = 0;
	this->charge = 0;
	this->noMolStru = false;
	this->molID = 0;

	this->createMol(tmp);
}

void molecule::createMol(std::string tmp)
{
	std::vector<std::string> oneLine;
	std::vector<std::string>::iterator it;
	int count = -1;
	while(++count<=3)
	{
		oneLine.push_back(tmp.substr(0,tmp.find(CRLF)));
		tmp = tmp.substr(tmp.find(CRLF)+CRLF.length());
	}
	if(oneLine[3].find("V2000") != std::string::npos)
		this->molFileVer = "V2000";
/*	else if(oneLine[3].find("V3000") != std::string::npos)
		this->molFileVer = "V3000";   */
	else
	{
		oneLine.push_back(tmp.substr(0,tmp.find(CRLF)));
		tmp = tmp.substr(tmp.find(CRLF)+CRLF.length());
		if(oneLine[4].find("V2000") != std::string::npos)
			this->molFileVer = "V2000";
/*		else if(oneLine[4].find("V3000") != std::string::npos)
			this->molFileVer = "V3000";  */
		if(this->molFileVer != "")
		{
			it = oneLine.begin();
			//oneLine.erase(&oneLine[0]);
			oneLine.erase(it);
		}
	}

	if(this->molFileVer == "V2000"/* || this->molFileVer == "V3000"*/)
	{
		this->molCount++;
		this->molID = molCount;
		this->molName = oneLine[0];
		this->molBuilder = oneLine[1];
		this->molComment = oneLine[2];
		this->molInfor = oneLine[3];
		std::string head = "";
		this->molMap = tmp.substr(0,tmp.find(molEnd));
		this->molBody = this->molInfor + CRLF + tmp;
		//std::cout<<"\r\n"<<this->molMap;system("pause");
	}
	else
	{
		this->molFileVer = "Wrong MOL formate";
		this->noMolStru = true;
		this->errMolCount++;
		//this->molCount++;
	}
	//std::cout<<"done\r\n";
}

void molecule::setABC()
{
	//std::cout<<"seting "<<this->molName<<"'s atom bond count...";
	//if(this->noMolStru) {std::cout<<"Wrong MOL formate\r\n";return;};

	this->atomCount = atoi(this->molInfor.substr(0,3).c_str());
	this->bondCount = atoi(this->molInfor.substr(3,3).c_str());
	if(this->atomCount == 0) this->noMolStru = true;
	//std::cout<<this->atomCount<<std::endl;system("pause");

	//std::cout<<"done\r\n";
}

void molecule::setCCP()
{
	//std::cout<<"seting "<<this->molName<<"'s coordinate connection...";
	//if(this->noMolStru) {std::cout<<"Wrong MOL formate\r\n";return;};
	long pos = this->molBody.find(molEnd) + molEnd.length() + CRLF.length();
	if(this->molBody.length() > pos)
		this->molProps = this->molBody.substr(pos);
	//std::cout<<this->molProps<<"$"<<std::endl;
	if(this->molFileVer!="V2000" || this->molMap.length()<10)
	{
		std::cout<<"Wrong file formate, cann't set Atom";
		exit(0);
	}
	else
	{
		int tmpCount = 0;
		std::string tmp = this->molMap;
		this->molMap = "";
		std::vector<std::string> oneLine;
		while(tmp.find(CRLF) != std::string::npos)
		{
			oneLine.push_back(tmp.substr(0,tmp.find(CRLF)));
			tmp = tmp.substr(tmp.find(CRLF)+CRLF.length());
			if(tmpCount++ < this->atomCount + this->bondCount)
			{
				this->molMap.append(oneLine[oneLine.size()-1]);
				this->molMap.append(CRLF);
			}
		}
		if(oneLine.size() >= (this->atomCount + this->bondCount))
		{
			int i=0;
			atom tmpAtom;
			for(i=0;i<this->atomCount;i++)
			{	
				tmp = oneLine[i].substr(30,3);	
				if(trim(tmp) != "H")
					atoms.push_back(atom(oneLine[i], i+1, false));
				else
					atoms.push_back(atom(oneLine[i], i+1, true));

				if(!this->rgroup)
				{
					if(atoms[atoms.size()-1].RGPAtom)
						this->rgroup = true;
				}
			}
			for(i=this->atomCount;i<(this->atomCount + this->bondCount);i++)
				bonds.push_back(bond(oneLine[i]));
			for(i=(this->atomCount + this->bondCount);i<oneLine.size();i++)
			{
				this->molExtInfor.append(oneLine[i]);
				if(oneLine[i].find(mdlRGroupMark) != std::string::npos)
					this->MDLRGP.push_back(rGroupMDL(oneLine[i]));
				else if(oneLine[i].find(CtabRLOGMark) != std::string::npos)
					this->CtabRGP.push_back(rGroupCtab(oneLine[i]));
			}
			if(this->MDLRGP.size() > 0 && this->CtabRGP.size() == 0)
			{
				for(i=0;i<this->MDLRGP.size();i++)
					this->atoms[this->MDLRGP[i].serNum-1].rgpNO
					= this->MDLRGP[i].groupNO;
			}

			if(this->molExtInfor != "")
				this->molExt = true;
		}
		else
			this->noMolStru = true;
	}

	//std::cout<<"done\r\n";
}

void molecule::setProp()
{
	//std::cout<<"seting "<<this->molName<<"'s property...";
	//if(this->noMolStru) {std::cout<<"Wrong MOL formate\r\n";return;};

	std::string tmp = this->molProps;
	std::string tmpP = "";
	std::vector<std::string> oneLine;

	while(tmp.find(CRLF) != std::string::npos)
	{
		oneLine.push_back(tmp.substr(0,tmp.find(CRLF)));
		tmp = tmp.substr(tmp.find(CRLF) + CRLF.length());
	}
	
	for(int i=0;i<oneLine.size();i++)
	{
		std::string tmpL = oneLine[i];
		if(tmpL.find(molPropMark) != std::string::npos)
		{
			tmpL = tmpL.substr(molPropMark.length());
			this->propNames.push_back(tmpL.substr(0,tmpL.find(">")));
			if(tmpP != "")
			{
				this->propValues.push_back(tmpP);
				tmpP = "";
			}
		}
		else if(tmpL.length() > 0)
		{
			tmpP += oneLine[i];
		}
	}

	if(this->propNames.size() - this->propValues.size() == 1)
		this->propValues.push_back(tmpP);

	if(this->propNames.size() == this->propValues.size())
		this->propCount = propNames.size();
	else
		this->propCount = 0;
	//std::cout<<this->propCount<<"$"<<std::endl;
	//showVector(this->propNames);
	//showVector(this->propValues);

	//std::cout<<"done\r\n";
}

bool molecule::init()
{
	//std::cout<<this->molID<<CRLF;

	//if(this->noMolStru) return;		//wrong format will be no structure
	molecule::setABC();
	if(this->noMolStru) return false;		//atoms and bonds 0

	molecule::setCCP();	
	molecule::setProp();

	return true;
}

void molecule::setSalt()
{
	//this->setAB();//re do this after aroma set

	int i=0, j=0, k=0, l=0;
	this->charge = 0;
	fragment tmpFrag;
	std::vector<atom> tmpAtoms;
	bool exist = false;
	for(i=0;i<this->atomCount;i++)
	{
		exist = false;
		for(l=0;l<tmpAtoms.size();l++)
		{
			if(tmpAtoms[l].serNum == this->atoms[i].serNum)
			{	
				exist = true;
				break;
			}
		}
		if(exist)
			continue;

		if(tmpFrag.addAtom(this->atoms[i]))
		{
			tmpAtoms.push_back(this->atoms[i]);
			for(j=0;j<tmpFrag.atomCount;j++)
			{
				for(k=0;k<tmpFrag.fAtoms[j].bCount;k++)
				{
					tmpFrag.addAtom(this->atoms[tmpFrag.fAtoms[j].bonds[k].bB-1]);
					tmpAtoms.push_back(this->atoms[tmpFrag.fAtoms[j].bonds[k].bB-1]);
				}
			}
			tmpFrag.setFrag();
			tmpFrag.setSalt();
			this->charge += tmpFrag.charge;
			this->salts.push_back(tmpFrag);
			tmpFrag.clear();
		}
	}

	this->sIonCount = this->salts.size();
	if(this->sIonCount > 1) this->salt = true;
	saltCount += this->sIonCount;
	if(this->salt && !this->ctable)
		sortFrag(this->salts);
}

bool molecule::setOSAtoms()
{
	if(this->salt)
	{
		int i=0, j=0;
		this->osAtoms.clear();
		for(i=0;i<this->sIonCount;i++)
		{
			if(this->salts[i].oIon)
			{
				for(j=0;j<this->salts[i].atomCount;j++)
					this->osAtoms.push_back(this->salts[i].fAtoms[j]);
			}
		}
		return true;
	}
	else
		return false;
}

bool molecule::chkIon()
{
	if(!this->salt)
	{
		this->ion = false;
		return false;
	}
	else
	{
		int i=0;
		for(i=0; i<this->sIonCount; i++)
		{
			if(this->salts[i].oIon)
			{
				this->ion = false;
				return false;
			}
		}
		this->ion = true;
		return true;
	}
}

void molecule::addRAtom()
{
	int i=0;
	rGroupMDL tmpMDLrgp;
	atom Ratom = this->atoms[0];
	bond Rbond = this->bonds[0];
	std::vector<atom>::iterator it = this->atoms.begin();
	std::vector<bond>::iterator it2 = this->bonds.begin();

	Ratom.chemSym = 'R';
	Ratom.chemSymEx = '#';
	Rbond.bA = 1;
	Rbond.bB = 2;
	Rbond.bType = 1;
	for(i=0;i<this->bondCount;i++)
	{
		this->bonds[i].bA++;
		this->bonds[i].bB++;
	}
	tmpMDLrgp.firstNO = 1;
	tmpMDLrgp.serNum = 1;
	tmpMDLrgp.groupNO =atoi(this->molName.c_str());
	this->MDLRGP.push_back(tmpMDLrgp);

	this->atoms.insert(it, Ratom);
	this->bonds.insert(it2, Rbond);
	this->atomCount = this->atomCount + 1;
	this->bondCount = this->bondCount + 1;
}

void molecule::setRGroup()
{
	if(this->rgroup && this->MDLRGP.size() ==0)
		this->rgroup = false;
	if(this->rgroup && this->salt)
	{
		int i=0, j=0, k=0, l=0;
		bool find = false;
		bool added = false;
		fragment tmpFrag;
		std::vector<fragment> tmpFrags = this->salts;
		std::vector<rGroupMDL>::iterator itRGmdl;
		std::vector<markush>::iterator itMarku;
		std::vector<fragment>::iterator itF;
		connector conn;

		//set r-group information
		for(i=0;i<this->sIonCount;i++) 	this->salts[i].setMarkush();

		//locate base r-group
		for(i=1;i<this->sIonCount;i++)
		{
			if(this->salts[0].markushs.size() < this->salts[i].markushs.size())
			{
				tmpFrag.clear();
				tmpFrag = this->salts[0];
				this->salts[0] = this->salts[i];
				this->salts[i] = tmpFrag;
			}
		}
		this->salts[0].baseRGP = true;
		
		//delete empty r-group
		for(i=0;i<this->sIonCount;i++)
		{
			for(j=0;j<this->salts[i].markushs.size();j++)
			{
				find = false; 
				for(k=0;k<this->sIonCount;k++)
				{
					if(k == i) continue;
					if(this->salts[k].chkMarkush(this->salts[i].markushs[j]))
					{
						find = true;
						break;
					}
				}
				if(!find)
				{
					this->salts[i].eraseAtom(this->salts[i].markushs[j].pointAtom);
					for(k=0;k<this->MDLRGP.size();k++)
					{
						if(this->MDLRGP[k].groupNO == this->salts[i].markushs[j].groupNO)
						{
							//this->MDLRGP.erase(&this->MDLRGP[k]);
							itRGmdl = this->MDLRGP.begin();
							this->MDLRGP.erase(itRGmdl + k);
							k -= 1;
						}
					}
					//this->salts[i].markushs.erase(&this->salts[i].markushs[j]);
					itMarku = this->salts[i].markushs.begin();
					this->salts[i].markushs.erase(itMarku + j);
					j -= 1;
				}
			}
		}
		
		//change more-point (not base) r-group into one-point r-group
		for(i=1;i<this->sIonCount;i++)
		{
			if(this->salts[i].markushs.size() >= 2)
			{
				added = false;
				for(j=0;j<this->salts[i].markushs.size();j++)
				{
					for(k=1;k<this->sIonCount;k++)
					{
						if(this->salts[i].used && this->salts[i].maxRGNO > this->salts[k].markushs[0].groupNO) continue;
						if(k == i || this->salts[k].markushs.size() > 1 || this->salts[k].used) continue;
						else
						{
							if(this->salts[i].markushs[j] == this->salts[k].markushs[0])
							{
								conn.connectRGP(this->salts[i],this->salts[k],tmpFrag);
								tmpFrag.setMarkush();
								tmpFrag.used = true;
								this->salts.push_back(tmpFrag);
								added = true;
								tmpFrag.clear();
							}
						}
					}
				}
				if(added)
				{
					//this->salts.erase(&this->salts[i]);
					itF = this->salts.begin();
					this->salts.erase(itF + i);
					this->sIonCount = this->salts.size();
					i -= 1;
				}
			}
		}

		for(i=0;i<this->sIonCount;i++) this->salts[i].used = false;

		//generate represented fragments
		for(i=0;i<this->salts[0].markushs.size();i++)
		{
			for(j=1;j<this->sIonCount;j++)
			{
				if(this->salts[0].markushs[i] == this->salts[j].markushs[0])
				{
					conn.connectRGP(this->salts[0],this->salts[j],tmpFrag);
					tmpFrag.setMarkush();
					this->RGPFras.push_back(tmpFrag);
					tmpFrag.clear();
				}
			}
			break;
		}
		for(i=0;i<this->RGPFras.size();i++)
		{
			for(j=0;j<this->RGPFras[i].markushs.size();j++)
			{
				for(k=1;k<this->sIonCount;k++)
				{
					if(this->RGPFras[i].markushs[j] == this->salts[k].markushs[0])
					{
						conn.connectRGP(this->RGPFras[i],this->salts[k],tmpFrag);
						tmpFrag.setMarkush();
						this->RGPFras.push_back(tmpFrag);
						tmpFrag.clear();
					}
				}
				break;
			}
			if(this->RGPFras[i].markushs.size() > 0)
			{
				//this->RGPFras.erase(&this->RGPFras[i]);
				itF = this->RGPFras.begin();
				this->RGPFras.erase(itF + i);
				i -= 1;
			}
		}

		//modify and save
		for(i=0;i<this->RGPFras.size();i++) this->RGPFras[i].formateASN();
		this->salts = this->RGPFras;
		sortFrag(this->salts);
		this->sIonCount = this->salts.size();
		this->RGPFras = tmpFrags;
		saltCount -= this->RGPFras.size();
		saltCount += this->sIonCount;
	}
}

void molecule::setAroma()
{
	//this->setAB();

	int i=0,j=0;
	std::vector<atom> tmpAtoms;
	std::vector<fragment> tmpFrags;
	std::vector<fragment> tmpFrags2;
	fragment tmpFrag;

	for(i=0;i<this->atomCount;i++)
		if(this->atoms[i].nonC)					//get non-sp3 C atom, include other atom like O, N, S etc.
			tmpAtoms.push_back(this->atoms[i]);

	shake(tmpAtoms);

	while(tmpAtoms.size()>0)
	{
		formFrag(tmpAtoms,tmpFrag);
		tmpFrags.push_back(tmpFrag);
		//std::sort(tmpFrag.fAtoms.begin(), tmpFrag.fAtoms.end(), lessASN);
		tmpAtoms = AnotB(tmpAtoms,tmpFrag.fAtoms);
		tmpFrag.clear();
	}

	for(i=0;i<tmpFrags.size();i++)
	{
		shake(tmpFrags[i],tmpFrags2);
		if(tmpFrags2.size() == 0)
			this->aRings.push_back(tmpFrags[i]);
		else
		{
			for(j=0;j<tmpFrags2.size();j++)
				this->aRings.push_back(tmpFrags2[j]);
			tmpFrags2.clear();
		}
	}
/*
	for(i=0;i<this->aRings.size();i++)
	{
		this->aRings[i].setFrag();
		this->aRings[i].chkAroma();
		if(this->aRings[i].aroma)
		{
			this->arCount++;
			for(j=0;j<this->bondCount;j++)
			{
				if(this->aRings[i].chkBAB(this->bonds[j]))
				{
					if(this->bonds[j].bType == 1 || this->bonds[j].bType == 2)
						this->bonds[j].bType = 4;
					else
						this->bonds[j].bType = 5;
				}
			}
			for(j=0;j<this->atomCount;j++)
			{
				if(this->aRings[i].chkASN(this->atoms[j].serNum))
					this->atoms[j].aroma = true;
			}
		}
	}
*/

	for(i=0;i<this->aRings.size();i++)
	{
		this->aRings[i].setFrag();
		this->aRings[i].chkAroma();
		if(this->aRings[i].aroma)
		{
			this->arCount++;
			std::sort(this->aRings[i].fAtoms.begin(), this->aRings[i].fAtoms.end(), lessASN);
			std::sort(this->aRings[i].fBonds.begin(), this->aRings[i].fBonds.end(), lessBASN);

			for(j=0;j<this->bondCount;j++)
			{
				if(this->bonds[j].bA < this->aRings[i].fBonds[0].bA) continue;
				else if(this->aRings[i].chkBAB(this->bonds[j]))
				{
					if(this->bonds[j].bType == 1 || this->bonds[j].bType == 2)
						this->bonds[j].bType = 4;
					else
						this->bonds[j].bType = 5;
				}
			}
			for(j=0;j<this->atomCount;j++)
			{
				if(this->atoms[j].serNum < this->aRings[i].fAtoms[0].serNum) continue;
				else if(this->aRings[i].chkASN(this->atoms[j].serNum))
					this->atoms[j].aroma = true;
			}
		}
	}
	for(i=0; i<this->atomCount; i++)
	{
		if(this->atoms[i].aroma)
		{
			for(j=0; j<this->atoms[i].bCount; j++)
			{
				if(this->atoms[this->atoms[i].bonds[j].bB-1].aroma)
				{
					if(this->atoms[i].bonds[j].bType == 1 || this->atoms[i].bonds[j].bType == 2)
						this->atoms[i].bonds[j].bType = 4;
					else
						this->atoms[i].bonds[j].bType = 5;
				}
			}
		}
	}



	if(this->arCount > 0)
		this->aroma = true;
}

void molecule::setSaltAroma()
{
	if(this->aRings.size() > 0 && this->salts.size() > 0)
	{
		int i=0, j=0;
		for(i=0;i<this->salts.size();i++)
		{
			for(j=0;j<this->aRings.size();j++)
			{
				if(chkFragASN(this->salts[i],this->aRings[j]))
					this->salts[i].aroma = true;
			}
		}
	}
}

void molecule::delH()
{
	//if(this->noMolStru) {std::cout<<"Wrong MOL formate\r\n";return;};

	std::vector<atom> tmpAtoms;
	std::vector<bond>::iterator it;
	atom tmpAtom;
	bond tmpBond;
	int i=0;
	int j=0;
	int Hcount=0;

	for(i=0;i<this->atomCount;i++)
	{
		if(!this->atoms[i].atomH)
		{
			tmpAtom = this->atoms[i];
			tmpAtom.serNum = tmpAtom.serNum - Hcount;
			tmpAtoms.push_back(tmpAtom);
			for(j=0;j<this->bondCount;j++)
			{
				if(this->bonds[j].bA == i+1)
					this->bonds[j].bA = tmpAtom.serNum;
				else if(this->bonds[j].bB == i+1)
					this->bonds[j].bB = tmpAtom.serNum;
			}
		}
		else
		{
			Hcount++;
			for(it=this->bonds.begin();it<this->bonds.end();it++)
			{
				tmpBond = *it;
				if(tmpBond.bA == i+1 || tmpBond.bB == i+1)
				{
					this->bonds.erase(it);
					this->bondCount = this->bonds.size();
				}
			}
		}
	}

	if(Hcount>0)
	{
		this->atoms = tmpAtoms;
		this->atomCount = this->atoms.size();

		char atomC[3];
		char bondC[3];
		sprintf(atomC,"% 3d",this->atomCount);
		sprintf(bondC,"% 3d",this->bondCount);
		std::string a = atomC;
		std::string b = bondC;
		std::string head="";
		this->molComment = "MolEngine delH";
		this->molInfor = a.append(b).append(this->molInfor.substr(6));
	}
}

void molecule::setAB()
{
	if(this->noMolStru) {std::cout<<"Wrong MOL formate\r\n";return;};
	this->setBC();
	this->sortAB();
	int i=0, j=0;
	bond tmpBond;
	std::vector<bond> tmpB;

	for(i=0;i<this->atomCount;i++)
	{
		tmpB.clear();
		for(j=0;j<this->bondCount;j++)
		{
			if(this->atoms[i].serNum < this->bonds[j].bA) break;
			else if(this->atoms[i].serNum == bonds[j].bA)
			{
				tmpB.push_back(bonds[j]);
			}
			else if(this->atoms[i].serNum == bonds[j].bB)
			{
				tmpBond = bonds[j];
				tmpBond.bA = this->atoms[i].serNum;
				tmpBond.bB = bonds[j].bA;
				tmpBond.aACS = this->atoms[i].chemSym;
				tmpBond.aACSE = this->atoms[i].chemSymEx;
				tmpBond.aBCS = bonds[j].aACS;
				tmpBond.aBCSE = bonds[j].aACSE;
				tmpB.push_back(tmpBond);
			}
		}
		this->atoms[i].setBonds(tmpB);
	}
}

void molecule::setBC()
{
	int i=0, j=0;
	for(i=0;i<this->atomCount;i++)
	{
		for(j=0;j<this->bondCount;j++)
		{
			if(this->atoms[i].serNum == this->bonds[j].bA)
			{
				this->bonds[j].aACS = this->atoms[i].chemSym;
				this->bonds[j].aACSE = this->atoms[i].chemSymEx;
			}
			else if(this->atoms[i].serNum == this->bonds[j].bB)
			{
				this->bonds[j].aBCS = this->atoms[i].chemSym;
				this->bonds[j].aBCSE = this->atoms[i].chemSymEx;
			}
		}
	}
}

void molecule::setStereo()
{
	int i=0;
	for(i=0; i<this->bondCount; i++)
	{
		if(this->bonds[i].stereo !=0 && this->bonds[i].bB - 1 < this->atomCount)
		{
			if(this->atoms[this->bonds[i].bB - 1].stereo == 0)
				this->atoms[this->bonds[i].bB - 1].stereo = this->bonds[i].stereo;
			else
				this->atoms[this->bonds[i].bB - 1].stereo = 0;
		}
	}
}

void molecule::sortAB()
{
	if(this->atomCount >= 2)
		std::sort(this->atoms.begin(), this->atoms.end(), lessASN);
	if(this->bondCount >=2)
	{
		int i=0;
		for(i=0; i<this->bondCount; i++)
			this->bonds[i].standardize();	
		std::sort(this->bonds.begin(), this->bonds.end(), lessBASN);
	}
}

void molecule::setFragEx(int num)
{
	int i=0;
	switch(num)
	{
	case 0:
		for(i=0;i<this->fragCount;i++)
		{
			this->frags[i].setABTs();
			this->frags[i].setADs();
			this->frags[i].findAA();
			this->frags[i].setHetero();
		}
		break;
	case 1:
		for(i=0;i<this->fragCount;i++)
		{
			this->frags[i].setABTs();
			//this->frags[i].setADs();
			//this->frags[i].findAA();
			//this->frags[i].setHetero();
		}
		break;
	case 2:
		for(i=0;i<this->fragCount;i++)
		{
			//this->frags[i].setABTs();
			this->frags[i].setADs();
			//this->frags[i].findAA();
			//this->frags[i].setHetero();
		}
		break;
	case 3:
		for(i=0;i<this->fragCount;i++)
		{
			//this->frags[i].setABTs();
			//this->frags[i].setADs();
			this->frags[i].findAA();
			this->frags[i].setHetero();
		}
		break;
	default:
		break;
	}
}

void molecule::getProp(int index, std::string &propN, std::string &propV)
{
	//if(this->noMolStru) {std::cout<<"Wrong MOL formate\r\n";return;};
	propN = this->propNames[index-1];
	propV = this->propValues[index-1];
}
void molecule::getProp(std::vector<std::string> &propN, std::vector<std::string> &propV)
{
	//if(this->noMolStru) {std::cout<<"Wrong MOL formate\r\n";return;};
	propN = this->propNames;
	propV = this->propValues;
}

void molecule::prepare()
{
	if(!this->init()) return;
	
	int i=0;
	if(this->atomCount > 1)
		this->delH();
	this->setAB();
	//this->setStereo();
	this->setAroma();
	this->setSalt();
	this->chkIon();
	this->chkStarAtom();
	this->setRGroup();
	this->setSaltAroma();
	this->setProp();
	if(this->salt)
	{
		this->setOSAtoms();
		for(i=0;i<this->sIonCount;i++)
		{
			this->salts[i].setABTs();
			this->salts[i].setHetero();
			this->salts[i].findAA();
			this->salts[i].setSalt();
			this->salts[i].sortAtoms2();
		}
	}
	else
	{
		this->salts[0].setABTs();
		this->salts[0].setHetero();
		this->salts[0].findAA();
		this->salts[0].setSalt();
		this->salts[0].sortAtoms2();
	}
}

bool molecule::chkStarAtom()
{
	int i=0;
	for(i=0; i<this->atomCount; i++)
	{
		if(this->atoms[i].star)
		{
			return true;
		}
	}
	return false;
}

void molecule::showSalt()
{
	if(this->sIonCount <= 1 && this->charge==0)
		std::cout<<this->molName<<" is not salt"<<CRLF;
	else
	{
		std::cout<<this->molName<<" is salt consist of "<<this->sIonCount<<" ions (parts)"<<CRLF;
		for(int i=0;i<this->sIonCount;i++)
		{
			this->salts[i].showFrag();
			if(this->salts[i].ioIon)
				std::cout<<"Inorganic, ";
			if(this->salts[i].oIon)
				std::cout<<"Organic, ";
			std::cout<<"Charge value: ";
			std::cout<<this->salts[i].charge;
			std::cout<<CRLF;
		}
	//system("pause");
	} 
}

void molecule::showAroma()
{
	if(this->aroma)
	{
		std::cout<<this->molName<<" is aromatic, has "<<this->arCount<<" aromatic fragment(s)"<<CRLF;
		for(int i=0;i<this->aRings.size();i++)
			this->aRings[i].showFrag();
	}
	else
	{
		std::cout<<this->molName<<" is not aromatic"<<CRLF;
	}
}

void molecule::showMol()
{
	//if(this->noMolStru) {std::cout<<"Wrong MOL formate\r\n";return;};
	std::cout<<this->molName<<CRLF;
	std::cout<<this->molBuilder<<CRLF;
	std::cout<<this->molComment<<CRLF;
	std::cout<<std::setw(3)<<this->atomCount;
	std::cout<<std::setw(3)<<this->bondCount;
	std::cout<<this->molInfor.substr(6)<<CRLF;
	int i=0;
	for(i=0;i<this->atomCount;i++)
	{
		std::cout<<std::setiosflags(std::ios::fixed);
		std::cout<<std::setw(10)<<std::setprecision(4)<<this->atoms[i].x;
		std::cout<<std::setw(10)<<std::setprecision(4)<<this->atoms[i].y;
		std::cout<<std::setw(10)<<std::setprecision(4)<<this->atoms[i].z;
		std::cout<<' ';
		std::cout<<this->atoms[i].chemSym<<this->atoms[i].chemSymEx<<' ';
		std::cout<<std::setw(2)<<this->atoms[i].massD;
		std::cout<<std::setw(3)<<this->atoms[i].charge;
		//std::cout<<std::setw(3)<<this->atoms[i].stereo;
		std::cout<<CRLF;
	}
	for(i=0;i<this->bondCount;i++)
	{
		std::cout<<std::setw(3)<<this->bonds[i].bA;
		std::cout<<std::setw(3)<<this->bonds[i].bB;
		std::cout<<std::setw(3)<<this->bonds[i].bType;
		std::cout<<std::setw(3)<<this->bonds[i].stereo;
		std::cout<<CRLF;
	}
	if(this->rgroup)
	{
		for(i=0;i<this->MDLRGP.size();i++)
		{
			std::cout<<mdlRGroupMark<<std::setw(3)<<this->MDLRGP[i].firstNO
				<<std::setw(4)<<this->MDLRGP[i].serNum
				<<std::setw(4)<<this->MDLRGP[i].groupNO<<CRLF;
		}
	}
	std::cout<<molEnd<<CRLF;
}

void molecule::outMol(std::fstream& outFile)
{
	//if(this->noMolStru) {std::cout<<"Wrong MOL formate\r\n";return;};
	outFile<<this->molName<<CRLF;
	outFile<<this->molBuilder<<CRLF;
	outFile<<this->molComment<<CRLF;
	outFile<<std::setw(3)<<this->atomCount;
	outFile<<std::setw(3)<<this->bondCount;
	outFile<<this->molInfor.substr(6)<<CRLF;
	int i=0;
	for(i=0;i<this->atomCount;i++)
	{
		outFile<<std::setiosflags(std::ios::fixed);
		outFile<<std::setw(10)<<std::setprecision(4)<<this->atoms[i].x;
		outFile<<std::setw(10)<<std::setprecision(4)<<this->atoms[i].y;
		outFile<<std::setw(10)<<std::setprecision(4)<<this->atoms[i].z;
		outFile<<' ';
		outFile<<this->atoms[i].chemSym<<this->atoms[i].chemSymEx<<' ';
		outFile<<std::setw(2)<<this->atoms[i].massD;
		outFile<<std::setw(3)<<this->atoms[i].charge;
		//outFile<<std::setw(3)<<this->atoms[i].stereo;
		outFile<<CRLF;
	}
	for(i=0;i<this->bondCount;i++)
	{
		outFile<<std::setw(3)<<this->bonds[i].bA;
		outFile<<std::setw(3)<<this->bonds[i].bB;
		outFile<<std::setw(3)<<this->bonds[i].bType;
		outFile<<std::setw(3)<<this->bonds[i].stereo;
		outFile<<CRLF;
	}
	if(this->rgroup)
	{
		for(i=0;i<this->MDLRGP.size();i++)
		{
			outFile<<mdlRGroupMark<<std::setw(3)<<this->MDLRGP[i].firstNO
				<<std::setw(4)<<this->MDLRGP[i].serNum
				<<std::setw(4)<<this->MDLRGP[i].groupNO<<CRLF;
		}
	}
	outFile<<molEnd<<CRLF;
}

void molecule::saveSDF(bool salt, bool removeEmptyMol, bool reviseBond)
{
	if(!salt)
	{
		if(removeEmptyMol)
		{
			if(!this->noMolStru)
			{
				this->showMol();
				std::cout<<SDFdelim<<CRLF;
			}
		}
		else
		{
			this->showMol();
			std::cout<<SDFdelim<<CRLF;
		}
	}
	else
	{
		if(!this->noMolStru)
		{
			int i=0;
			int count = 0;
			char tmp[100];
			for(i=0;i<this->sIonCount;i++)
			{
				if(this->salts[i].oIon)
				{
					//itoa(i+1,tmp,10);
					sprintf(tmp,"%d",i+1);
					std::string tmpStr = this->molName + "-" + tmp;
					this->salts[i].showMol(tmpStr,this->molInfor, count);
					if(reviseBond)
						count = this->salts[i].atomCount;
					std::cout<<molEnd<<CRLF;
					std::cout<<SDFdelim<<CRLF;
				}
			}
		}
	}
}

