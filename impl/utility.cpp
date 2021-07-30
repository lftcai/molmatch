#include "utility.h"

fileReader::fileReader(std::string path)
{
	extern std::string CRLF;
	this->filePath = path;

	std::fstream openFile;
	openFile.open(path.c_str(),std::ios::in);
	if(!openFile)
	{
		std::cout<<"Open \""<<path<<"\" failed"<<CRLF;
		exit(0);
	}
	openFile.seekg(0,std::ios::end);
	this->fileSize = openFile.tellg(); 

	int tmp = path.find_last_of("\\");
	if( tmp > 0)
		this->fileName = path.substr(tmp+1);
	else
		this->fileName = path;

	if(this->fileName.find(".sdf") != std::string::npos)
		this->fileType = 1;
	else if(this->fileName.find(".ctb") != std::string::npos)
		this->fileType = 2;
	else
		this->fileType = 0;

	this->readedCount = 0;
	this->readedPos = 0;
	this->readToEnd = false;
	openFile.close();
}

void fileReader::read(char delim, long limit)
{
	extern std::string CRLF;
	int len = CRLF.length();
	std::string tmp;
	std::fstream openFile;
	openFile.open(this->filePath.c_str(),std::ios::in);
	if(!openFile)
	{
		std::cout<<"Open file "<<this->fileName<<" failed"<<CRLF;
		exit(0);
	}
	openFile.seekg(this->readedPos);
	if(limit == -1)
	{
		while(!openFile.eof())
		{
			std::getline(openFile,tmp,delim);
			if(tmp.size() <= len) continue;
			this->readedCount++;
			this->readedObjects.push_back(tmp);
		}
	}
	else
	{
		while(!openFile.eof() && this->readedCount < limit)
		{
			std::getline(openFile,tmp,delim);
			this->readedCount++;
			this->readedObjects.push_back(tmp);
		}
	}
	this->readToEnd = openFile.eof();
	this->readedPos = openFile.tellg();
	openFile.close();
}

void fileReader::Read()
{
	if(this->fileType == 1)
		this->read('$', -1);
	else if(this->fileType == 2)
		this->read('%', -1);
	else
		this->read('~', -1);
}

RGP::RGP(std::string &str)
{
	this->all = "";
	this->all.append(str);
	this->setRGP();
}

void RGP::setRGP()
{
	int i=0;
	long start = 0;
	long end = 0;
	extern std::string CRLF;

	end = this->all.find(CRLF);
	this->rgpNO = this->all.substr(0, end).c_str();

	this->all = this->all.substr(end + CRLF.length());
	
	while(this->all.find("$CTAB") != std::string::npos)
	{
		start = this->all.find("$CTAB") + strlen("$CTAB") + CRLF.length();
		end = this->all.find("$END CTAB");
		this->frags.push_back(this->all.substr(start,end-start));
		this->all = this->all.substr(end + strlen("$END CTAB") + CRLF.length());
	}
}

void RGP::show()
{
	extern std::string CRLF;
	int i=0;
	std::cout<<"Group: "<<this->rgpNO<<CRLF;
	for(i=0;i<this->frags.size();i++)
		std::cout<<frags[i];
}

CTable::CTable(std::string &str)
{
	this->all = "";
	this->all.append(str);
}

void CTable::setCTable()
{
	int i=0;
	long start = 0;
	long end = 0;
	extern std::string CRLF;
	std::vector<std::string> RGPs;

	start = this->all.find("$HDR") + strlen("$HDR") + CRLF.length();
	end = this->all.find("$END HDR");
	this->head = this->all.substr(start,end-start);
	this->head = this->head.substr(this->head.find(CRLF) + CRLF.length());
	this->head = "Molengine-CtableToMol" + CRLF + this->head;

	start = this->all.find("$CTAB") + strlen("$CTAB") + CRLF.length();
	end = this->all.find("$END CTAB");
	this->baseFrag = this->all.substr(start,end-start);

	this->all = this->all.substr(this->all.find("$RGP"));
	while(this->all.find("$RGP") != std::string::npos)
	{
		start = this->all.find("$RGP") + strlen("$RGP") + CRLF.length();
		end = this->all.find("$END RGP");
		RGPs.push_back(this->all.substr(start,end-start));
		this->all = this->all.substr(end + strlen("$END RGP") + CRLF.length());
	}

	for(i=0;i<RGPs.size();i++)
		this->rgps.push_back(RGP(RGPs[i]));
}

void CTable::createFrag()
{
	int i=0, j=0;
	extern std::string CRLF;
	std::string tmp = "";
	std::string tmp2 = "";
	std::string Head = this->head.substr(this->head.find(CRLF) + CRLF.length());

	tmp.append("base");
	tmp.append(CRLF);
	tmp.append(Head);
	tmp.append(this->baseFrag);
	this->baseFrag = tmp;

	for(i=0;i<this->rgps.size();i++)
	{
		tmp = "";
		tmp.append(this->rgps[i].rgpNO);
		tmp.append(CRLF);
		tmp.append(Head);
		for(j=0;j<this->rgps[i].frags.size();j++)
		{
			tmp2 = tmp;
			tmp2.append(rgps[i].frags[j]);
			rgps[i].frags[j] = tmp2;
		}
	}
}

void CTable::show()
{
	int i=0;
	extern std::string CRLF;
	std::cout<<this->head;
	std::cout<<"Base Fragment: "<<CRLF;
	std::cout<<this->baseFrag;
	for(i=0;i<this->rgps.size();i++)
	{
		this->rgps[i].show();
	}
}

void showVector(std::vector<int> v, int width, int delim)
{
	extern std::string CRLF;
	int end = 0;
	end = delim == 0 ? v.size() : delim;
	end = end >= v.size() ? v.size() : end;
	int control=1;
	for(int i=0;i<end;i++)
	{
		std::cout<<v[i]<<" ";
		if(control++ == width)
		{
			control = 1;
			std::cout<<CRLF;
		}
	}
	std::cout<<CRLF;
}

void showVector(std::vector<double> v, int width, int delim)
{
	extern std::string CRLF;
	int end = 0;
	end = delim == 0 ? v.size() : delim;
	end = end >= v.size() ? v.size() : end;

	int control=1;
	for(int i=0;i<end;i++)
	{
		std::cout<<v[i]<<" ";
		if(control++ == width)
		{
			control = 1;
			std::cout<<CRLF;
		}
	}
	std::cout<<CRLF;
}

void showVector(std::vector<long> v, int width, int delim)
{
	extern std::string CRLF;
	long end = 0;
	end = delim == 0 ? v.size() : delim;
	end = end >= v.size() ? v.size() : end;

	int control=1;
	for(long i=0;i<end;i++)
	{
		std::cout<<v[i]<<" ";
		if(control++ == width)
		{
			control = 1;
			std::cout<<CRLF;
		}
	}
	std::cout<<CRLF;
}

void showVector(std::vector<std::string> v, int width, int delim)
{
	extern std::string CRLF;
	int end = 0;
	end = delim == 0 ? v.size() : delim;
	end = end >= v.size() ? v.size() : end;

	int control=1;
	for(int i=0;i<end;i++)
	{
		std::cout<<v[i]<<" ";
		if(control++ == width)
		{
			control = 1;
			std::cout<<CRLF;
		}
	}
	std::cout<<CRLF;
}

long timeNow()
{
	time_t now;
	time(&now);
	return now;
}

void showVector(std::vector<simResult> &v)
{
	extern std::string CRLF;
	for(int i=0;i<v.size();i++)
	{
		std::cout<<v[i].idA<<"\t"<<v[i].idB<<"\t"<<v[i].sim<<CRLF;
	}
}

std::string timeElapse(long time1,long time2)
{
	int sec=0;
	int min=0;
	int hou=0;
	char second[5];
	char minite[5];
	char hour[5];
	long time_elapsed=0;
	std::string return_str;
	time_elapsed=time2-time1;
	hou = time_elapsed/3600 > 0 ? time_elapsed/3600 : 0;
	min = hou >0 ? ((time_elapsed%3600)/60 >0 ? (time_elapsed%3600)/60 : time_elapsed%3600) : (time_elapsed/60 > 0 ? (time_elapsed/60) : 0);
	sec = time_elapsed - min*60 - hou*3600;
	sec = sec ==0 ? 1 : sec;
	//itoa(sec,second,10);
	//itoa(min,minite,10);
	//itoa(hou,hour,10);
	sprintf(second,"%d",sec);
	sprintf(minite,"%d",min);
	sprintf(hour,"%d",hou);
	return_str.append(hour).append(":").append(minite).append(":").append(second);
	return return_str;
}

void connector::connectAtom(atom &atomA, atom &atomB,
							bond &retBond, int bondType, bool revise)
{
	int i=0, j=0;
	bond tmpBond;
	tmpBond.bType = bondType;

	tmpBond.aACS = atomA.chemSym;
	tmpBond.aACSE = atomA.chemSymEx;
	tmpBond.bA = atomA.serNum;
	tmpBond.aBCS = atomB.chemSym;
	tmpBond.aBCSE = atomB.chemSymEx;
	tmpBond.bB = atomB.serNum;
	atomA.bonds.push_back(tmpBond);
	if(revise)
	{
		atomA.bCount = atomA.bCount + 1;
		atomA.bCount2 = atomA.bCount2 + 1;
	}

	tmpBond.aACS = atomB.chemSym;
	tmpBond.aACSE = atomB.chemSymEx;
	tmpBond.bA = atomB.serNum;
	tmpBond.aBCS = atomA.chemSym;
	tmpBond.aBCSE = atomA.chemSymEx;
	tmpBond.bB = atomA.serNum;
	atomB.bonds.push_back(tmpBond);
	if(revise)
	{
		atomB.bCount = atomB.bCount + 1;
		atomB.bCount2 = atomB.bCount2 + 1;
	}

	tmpBond.standardize();
	retBond = tmpBond;
}

bool connector::connectRGP(fragment &fragA, fragment fragB,
					std::vector<atom>& retAtoms, std::vector<bond>& retBonds)
{
	fragA.setMarkush();

	int i=0, j=0, k=0, l=0;
	std::vector<markush> fragAMarkush = fragA.markushs;
	std::vector<markush> fragBMarkush = fragB.markushs;
	bond newBond;
	bond gotBond;
	int start = 0;
	int tmp = 0;
	bool added = false;
	for(i=0;i<fragA.markushs.size();i++)
	{
		for(j=0;j<fragB.markushs.size();j++)
		{
			if(fragA.markushs[i] == fragB.markushs[j]
				&& (fragA.markushs[i].atoms.size() == 1
					|| fragB.markushs[j].atoms.size() == 1))
			{
				added = false;
				fragA.setMaxASN();
				tmp = getMaxAsn(retAtoms);
				start = fragA.maxASN + 1 > tmp ? fragA.maxASN + 1 : tmp;
				fragB.reviseASN(start);
				fragB.setMarkush();
				for(k=0;k<fragA.markushs[i].atoms.size();k++)
				{
					fragA.markushs[i].getBond(fragA.markushs[i].atoms[k].serNum,gotBond);
					for(l=0;l<fragB.markushs[j].atoms.size();l++)
					{
						connectAtom(fragA.markushs[i].atoms[k],
									fragB.markushs[j].atoms[l],
									newBond,gotBond.bType,true);
						retBonds.push_back(newBond);
						added = true;
					}
				}
				if(added)
				{
					fragA.markushs[i].eraseAB();
					fragB.markushs[j].eraseAB();
					fragA.markushs[i].getAtoms(retAtoms);
					fragB.markushs[j].getAtoms(retAtoms);
					AnotB(fragB.fAtoms,retAtoms,retAtoms);
					fragB.getNMBonds(retBonds);
				}
			}
		}
	}

	//for(k=0;k<retAtoms.size();k++) retAtoms[k].showAtom();
	//for(k=0;k<retBonds.size();k++) retBonds[k].showBond();
	return added;
}

bool connector::connectRGP(fragment &fragA, fragment fragB, fragment &regFrag)
{
	fragA.setMarkush();

	int i=0, j=0, k=0, l=0;
	int start = 0;
	int tmp = 0;
	bool added = false;
	bond newBond;
	bond gotBond;
	std::vector<atom> retAtoms;
	for(i=0;i<fragA.markushs.size();i++)
	{
		for(j=0;j<fragB.markushs.size();j++)
		{
			if(fragA.markushs[i] == fragB.markushs[j]
				&& (fragA.markushs[i].atoms.size() == 1
					|| fragB.markushs[j].atoms.size() == 1))
			{
				added = false;
				fragA.setMaxASN();
				regFrag.setMaxASN();
				tmp = regFrag.maxASN;
				start = fragA.maxASN + 1 > tmp ? fragA.maxASN + 1 : tmp;
				fragB.reviseASN(start);
				fragB.setMarkush();
				for(k=0;k<fragA.markushs[i].atoms.size();k++)
				{
					fragA.markushs[i].getBond(fragA.markushs[i].atoms[k].serNum,gotBond);
					for(l=0;l<fragB.markushs[j].atoms.size();l++)
					{
						connectAtom(fragA.markushs[i].atoms[k],
									fragB.markushs[j].atoms[l],
									newBond,gotBond.bType,true);
						added = true;
					}
				}
				if(added)
				{
					fragA.markushs[i].eraseAB();
					fragB.markushs[j].eraseAB();
					fragA.markushs[i].getAtoms(retAtoms);
					fragB.markushs[j].getAtoms(retAtoms);
					AnotB(fragB.fAtoms,retAtoms,retAtoms);
				}
			}
		}
	}
	AnotB(fragA.fAtoms,retAtoms,retAtoms);
	formFrag(retAtoms,regFrag);
	regFrag.setFrag();
	retAtoms.clear();

	return added;
}
