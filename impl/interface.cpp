#include "interface.h"

std::vector<molecule> molPool;

bool getStrMols(std::string& mols, std::vector<std::string>& strMols)
{
	if(mols.find(SDFdelim) != std::string::npos)
	{
		std::string tmpStr = "";
		long pos = 0;
		do
		{
			tmpStr = mols.substr(0, mols.find(SDFdelim));
			if(tmpStr.find(molEnd) != std::string::npos)
				strMols.push_back(tmpStr);
			mols = mols.substr(mols.find(SDFdelim) + SDFdelim.length());
		}while(mols.find(SDFdelim) != std::string::npos);

		if(mols.find(molEnd) != std::string::npos)
			strMols.push_back(mols);
		if(!strMols.empty())
			return true;
		else
			return false;
	}
	else
		return false;
}

int molMatch(std::string mol, std::vector<std::string>& molIDs, int matchType, float st, int molType)
{
	std::vector<std::string> strMols;
	std::vector<molecule> mols;
	std::vector<simResult> rs;
	molecule newMol;
	long i=0;

	switch(molType)
	{
	case 0:
		if(getStrMols(mol, strMols))
		{
			for(i=0;i<strMols.size();i++)
			{
				molecule m(strMols[i]);
				m.prepare();
				mols.push_back(m);
			}
		}
		else
		{
			newMol.createMol(mol);
			newMol.prepare();
			mols.push_back(newMol);
		}
		break;
	case 1:
		CTable ct(mol);	
		ct.setCTable();
		ct.createFrag();
		std::string tmp = "";
		CtabToMol(ct,tmp);
		newMol.createMol(tmp);
		newMol.prepare();
		mols.push_back(newMol);
		break;
	}
	
	switch(matchType)
	{
	case 0:
		CalcSim(molPool,mols,rs,st);
		sort(rs);
		break;
	case 1:
		StrucMatch(molPool,mols,rs);
		break;
	case 2:
		SubSMatch(molPool,mols,rs);
		break;
	}

	//showVector(rs);
	if(rs.size() == 0)
	{
		return 0;
	}

	for(i=0;i<rs.size();i++)
	{
		molIDs.push_back(rs[i].idA);
	}

	return 1;
}

int molMatch(std::string& mol, std::vector<molecule>& molPool, std::vector<simResult>& rs, int matchType, int min_sim, int molType)
{
	if(mol.find(molEnd) == std::string::npos)
		return 0;

	std::vector<std::string> strMols;
	std::vector<molecule> mols;
	molecule newMol;
	long i=0;
	
	float st = (float)min_sim/100;

	switch(molType)
	{
	case 0:
		if(getStrMols(mol, strMols))
		{
			for(i=0;i<strMols.size();i++)
			{
				molecule m(strMols[i]);
				m.prepare();
				mols.push_back(m);
			}
		}
		else
		{
			newMol.createMol(mol);
			newMol.prepare();
			mols.push_back(newMol);
		}
		break;
	case 1:
		CTable ct(mol);	
		ct.setCTable();
		ct.createFrag();
		std::string tmp = "";
		CtabToMol(ct,tmp);
		newMol.createMol(tmp);
		newMol.prepare();
		mols.push_back(newMol);
		break;
	}
	
	switch(matchType)
	{
	case 0:
		CalcSim(molPool,mols,rs,st);
		sort(rs);
		break;
	case 1:
		StrucMatch(molPool,mols,rs);
		break;
	case 2:
		SubSMatch(molPool,mols,rs);
		break;
	}

	//showVector(rs);

	return rs.size();
}

bool molMatch(std::string mol, std::vector<molecule> &oneMolfromPool, int &resultValue, int matchType, int min_sim, int molType)
{
	std::vector<std::string> strMols;
	std::vector<molecule> mols;
	std::vector<molecule> tmpMolPool;
	std::vector<simResult> rs;
	molecule newMol;
	bool matchResult;

	tmpMolPool = oneMolfromPool;
	float st = (float)min_sim/100;

	switch(molType)
	{
	case 0:
		if(getStrMols(mol, strMols))
		{
			for(int i=0;i<strMols.size();i++)
			{
				molecule m(strMols[i]);
				m.prepare();
				mols.push_back(m);
			}
		}
		else
		{
			newMol.createMol(mol);
			newMol.prepare();
			mols.push_back(newMol);
		}
		break;
	case 1:
		CTable ct(mol);	
		ct.setCTable();
		ct.createFrag();
		std::string tmp = "";
		CtabToMol(ct,tmp);
		newMol.createMol(tmp);
		newMol.prepare();
		mols.push_back(newMol);
		break;
	}
	
	switch(matchType)
	{
	case 0:
		CalcSim(tmpMolPool,mols,rs,st);
		sort(rs);
		break;
	case 1:
		StrucMatch(tmpMolPool,mols,rs);
		break;
	case 2:
		SubSMatch(tmpMolPool,mols,rs);
		break;
	}

	if(!rs.empty() && (rs[0].sim >= st || rs[0].sim == -1))
	{
		resultValue = (int)rs[0].sim * 100;
		matchResult = true;
	}
	else 
		matchResult = false;

	strMols.clear();
	tmpMolPool.clear();
	mols.clear();
	rs.clear();

	return matchResult;
}

int molsMatch(std::string fileName, std::vector<std::string>& molIDs, int matchType, float st, int molType)
{
	std::vector<molecule> mols;
	std::vector<simResult> rs;
	long i=0;

	fileReader fr(fileName);
	fr.Read();
	for(i=0;i<fr.readedCount;i++)
	{
		//std::cout<<i<<CRLF;
		if(fr.fileType == 1)
		{
			molecule a(fr.readedObjects[i]);
			a.prepare();
			mols.push_back(a);
		}
		else if(fr.fileType == 2)
		{
			CTable ct(fr.readedObjects[i]);	
			ct.setCTable();
			ct.createFrag();
			std::string tmp = "";
			CtabToMol(ct,tmp);
			molecule a(tmp);
			a.prepare();
			mols.push_back(a);
		}
	}

	switch(matchType)
	{
	case 0:
		CalcSim(molPool,mols,rs,st);
		sort(rs);
		break;
	case 1:
		StrucMatch(molPool,mols,rs);
		break;
	case 2:
		SubSMatch(molPool,mols,rs);
		break;
	}

	//showVector(rs);
	if(rs.size() == 0)
	{
		return 0;
	}

	for(i=0;i<rs.size();i++)
	{
		molIDs.push_back(rs[i].idA);
	}

	return 1;
}

int loadMol(std::string fileName)
{
	long i=0;
	fileReader fr(fileName);
	fr.Read();
	for(i=0;i<fr.readedCount;i++)
	{
		//std::cout<<i<<CRLF;
		if(fr.fileType == 1)
		{
			molecule a(fr.readedObjects[i]);
			a.prepare();
			molPool.push_back(a);
		}
		else if(fr.fileType == 2)
		{
			CTable ct(fr.readedObjects[0]);	
			ct.setCTable();
			ct.createFrag();
			std::string tmp = "";
			CtabToMol(ct,tmp);
			molecule a(tmp);
			a.prepare();
			molPool.push_back(a);
		}
	}

	if(molPool.size() > 0)
		return 1;
	else
		return 0;
}

//int loadMol(std::string fileName, std::map<doc_id_t, molecule>& molPool)
//int loadMol(std::string fileName, std::map<std::string, molecule>& molPool)
int loadMol(std::string fileName, std::map<std::string, std::vector<molecule> > &molPool)
{
	long i=0;
	std::vector<molecule> tmpMol;
	std::map<std::string, std::vector<molecule> >::iterator pIt;
	fileReader fr(fileName);
	fr.Read();
	for(i=0;i<fr.readedCount;i++)
	{
		//std::cout<<i<<CRLF;
		if(fr.fileType == 1)
		{
			molecule a(fr.readedObjects[i]);
			a.prepare();
			if(a.molName.length() > 2)
			{
				pIt = molPool.find(a.molName);
				if(pIt != molPool.end())
				{
					tmpMol.push_back(a);
					mergeMols(tmpMol, pIt->second, pIt->second);
					tmpMol.clear();
				}
				else
				{
					tmpMol.push_back(a);
					molPool.insert(std::pair<std::string, std::vector<molecule> >(a.molName, tmpMol));
					tmpMol.clear();
				}

			}
		}
		else if(fr.fileType == 2)
		{
			CTable ct(fr.readedObjects[0]);	
			ct.setCTable();
			ct.createFrag();
			std::string tmp = "";
			CtabToMol(ct,tmp);
			molecule a(tmp);
			a.prepare();
			if(a.molName.length() > 2)
			{
				pIt = molPool.find(a.molName);
				if(pIt != molPool.end())
				{
					tmpMol.push_back(a);
					mergeMols(tmpMol, pIt->second, pIt->second);
					tmpMol.clear();
				}
				else
				{
					tmpMol.push_back(a);
					molPool.insert(std::pair<std::string, std::vector<molecule> >(a.molName, tmpMol));
					tmpMol.clear();
				}

			}
		}
	}
	if(molPool.size() > 0)
		return 1;
	else
		return 0;
}

bool mergeMols(std::vector<molecule>& vA, std::vector<molecule>& vB, std::vector<molecule>& retv)
{
	std::vector<molecule> tmpV;
	std::vector<molecule> tmpV2;
	std::vector<molecule> tmpV3;
	std::vector<molecule>::iterator itA;
	std::vector<molecule>::iterator itB;
	
	tmpV.insert(tmpV.begin(), vA.begin(), vA.end());
	tmpV.insert(tmpV.end(), vB.begin(), vB.end());

	retv.clear();

	for(itA=tmpV.begin(); itA!=tmpV.end(); itA++)
	{
		tmpV2.push_back(*itA);
		itB = itA + 1;
		if(itB != tmpV.end())
		{
			tmpV3.insert(tmpV3.begin(), itB, tmpV.end());
			if(!chkMolExist(tmpV2, tmpV3))
				retv.push_back(*itA);
		}
		else
			retv.push_back(*itA);
		tmpV2.clear();
		tmpV3.clear();
	}
	
	if(!retv.empty())
		return true;
	else
		return false;
}

int loadMol(std::string fileName, std::vector<molecule> &molPool)
{
	long i=0;
	std::vector<molecule> tmpMol;
	std::vector<molecule>::iterator it;
	std::string tmpStr = "";

	fileReader fr(fileName);
	fr.Read();
	for(i=0;i<fr.readedCount;i++)
	{
		if(fr.fileType == 1)
		{
			molecule a(fr.readedObjects[i]);
			a.prepare();
			if(a.molName.length() > 2)
			{
				tmpMol.push_back(a);
				if(!chkMolExist(tmpMol, molPool, true))
					molPool.push_back(a);
				tmpMol.clear();
			}
		}
		else if(fr.fileType == 2)
		{
			CTable ct(fr.readedObjects[0]);	
			ct.setCTable();
			ct.createFrag();
			std::string tmp = "";
			CtabToMol(ct,tmp);

			molecule a(tmp);
			a.prepare();
			if(a.molName.length() > 2)
			{
				tmpMol.push_back(a);
				if(!chkMolExist(tmpMol, molPool, true))
					molPool.push_back(a);
				tmpMol.clear();
			}
		}
	}
	
	for(it = molPool.begin(); it != molPool.end(); it++)
	{
		modifyStr(it->molName, tmpStr, "+");
		it->molName = tmpStr;
		tmpStr = "";
	}
	return molPool.size();
}

int makeSocketMolServer(int listenPort, std::vector<molecule> &molPool, std::vector<simResult> &rs)
{
	socketConn sc;
	int r = sc.serverListen(listenPort);
	if(r)
	{
		do
		{
			r = sc.fetchConnection();
			if(r)
			{
				char buff[1024];
				bzero(tmpbuf, 1024);
				r = sc.receiveData(buff, sizeof(buff));
				if(r > 0)
				{
					std::cout<<"server rec: "<<buff<<std::endl;
					std::string tmpstr = buff;
					int n = 0;
					n = molMatch(tmpstr, molPool, rs, 0);
					char tmpbuf[1024];
					bzero(tmpbuf, 1024);
					sprintf(tmpbuf, "%d", n+'\0');
					sc.writeOnConnection(tmpbuf, r);
					sc.closeConnection();
				}
			}
		}while(1);
	}
	sc.closeConnection();
	sc.closeSocket();
	return -1;
}

