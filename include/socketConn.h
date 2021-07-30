/*************************************************************************
    > File Name: include/socketConn.h
    > Author: Frank Lu
    > Mail: lufengfrank@didichuxing.com
    > Created Time: Tue 15 Sep 2015 07:39:34 PM CST
    > Utility: connet by socket
 ************************************************************************/

#ifndef SOCKETCONN
#define SOCKETCONN
#include <unistd.h> //close
#include <arpa/inet.h>//socket
#include <signal.h> //signal(SIGPIPE)
#include <string.h>


class socketConn
{
private:
	int m_socketfd;
	int m_connection;
	int m_max_conn;
public:
	socketConn();
	bool serverListen(int port);
	bool clientBind();
	bool clientConnect(char *serverIP, int serverPort);
	bool fetchConnection();
	int receiveData(char *receivedData, size_t buffSize);
	int sendData(const char *sendData);
	int writeOnConnection(const void *buf, size_t size);
	int readOnConnection(void *buf, size_t size);
	int readFromSocket(void *buff, size_t size);
	void closeSocket();
	void closeConnection();

};

#endif
