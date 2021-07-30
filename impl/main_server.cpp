#include "interface.h"
int main()
{
	std::vector<molecule> mpool;
	std::cout<<"Initializing...";
	if(loadMol("/Users/lftcai/mol/impl/wtmmol.sdf", mpool))// create mol pool
		std::cout<<"done"<<CRLF;
	else
	{
		std::cout<<"failed"<<CRLF;
		exit(0);
	}
	
	std::cout<<"Total: "<<mpool.size()<<" unique molecules"<<std::endl;	
	
	std::vector<simResult> rs;
	int a = makeSocketMolServer(9988, mpool, rs);
	return 0;
}

