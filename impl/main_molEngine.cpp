#include "interface.h"
int main(int argc, char *argv[])
{
//test interface

	if(argc < 2)
	{
		std::cout<<"Need file input\n";
		exit(0);
	}

	fileReader fr(argv[1]);// simulate molecule input
	fr.Read();
	std::vector<molecule> mpool;

	std::cout<<"Initializing...";
	//if(loadMol("all.sdf", mpool))// create mol pool
	if(loadMol("wtmmol.sdf", mpool))// create mol pool
		std::cout<<"done"<<CRLF;
	else
	{
		std::cout<<"failed"<<CRLF;
		exit(0);
	}
	
	std::cout<<"Total: "<<mpool.size()<<" unique molecules"<<std::endl;	
	
	std::vector<simResult> rs;
	std::cout<<"Hit: "<<molMatch(fr.readedObjects[0], mpool, rs, 0)<<std::endl;
/*	
	extern std::string molPropMark;
	extern std::string molPropMarkB;
	std::string tmpStr;
	std::vector<molecule>::iterator it;
	for(it = mpool.begin(); it != mpool.end(); it++)
	{
		tmpStr = it->molName;
		it->molName = "";
		it->showMol();
		std::cout<<molPropMark<<"PNS"<<molPropMarkB<<CRLF;
		std::cout<<tmpStr<<CRLF<<CRLF<<SDFdelim<<CRLF;
	}	
*/
	return 0;
}

