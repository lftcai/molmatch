#include "interface.h"
int main()
{
	fileReader fr("/Users/lftcai/mol/impl/demo.sdf");// simulate molecule input
	fr.Read();
	socketConn sc;
	sc.clientConnect("127.0.0.1", 9988);
	int l = sc.sendData(fr.readedObjects[0].c_str());
	std::cout<<"sent len: "<<l<<",data len: "<<fr.readedObjects[0].length()<<"\n";
	char tmpbuf[1024];
	bzero(tmpbuf,1024);
	std::cout<<"server response:"<<std::endl;
	sc.readFromSocket(tmpbuf,1024);
	std::cout<<tmpbuf<<std::endl;
	sc.closeConnection();
	sc.closeSocket();
	return 0;
}

