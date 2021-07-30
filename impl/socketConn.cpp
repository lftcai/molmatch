/*************************************************************************
    > File Name: impl/socketConn.cpp
    > Author: Frank Lu
    > Mail: lufengfrank@didichuxing.com
    > Created Time: Tue 15 Sep 2015 08:09:32 PM CST
    > Utility: connet by socket
 ************************************************************************/

#include "socketConn.h"

socketConn::socketConn()
{
	m_socketfd = socket(AF_INET,SOCK_STREAM, 0);
	m_connection = 0;
	m_max_conn = 100;
	signal(SIGPIPE, SIG_IGN);
}

bool socketConn::serverListen(int port)
{
	int err = 0;
	struct sockaddr_in server_sockaddr;
	server_sockaddr.sin_family = AF_INET;
	server_sockaddr.sin_port = htons(port);
	server_sockaddr.sin_addr.s_addr = htonl(INADDR_ANY);

	err = bind(m_socketfd, (struct sockaddr *)&server_sockaddr,sizeof(server_sockaddr));

	if(err != 0)
		return false;
	
	err = listen(m_socketfd, m_max_conn);

	if(err == -1)
		return false;

	return true;
}


bool socketConn::clientBind()
{
	int err = 0;
	struct sockaddr_in client_sockaddr;
	client_sockaddr.sin_family = AF_INET;
	client_sockaddr.sin_port = htons(0);
	client_sockaddr.sin_addr.s_addr = htonl(INADDR_ANY);

	err = bind(m_socketfd, (struct sockaddr *)&client_sockaddr,sizeof(client_sockaddr));

	if(err != 0)
		return false;
	
	return true;
}

bool socketConn::clientConnect(char *serverIP, int serverPort)
{
	int err = 0;
	struct sockaddr_in server_addr;
	memset(&server_addr, 0, sizeof(server_addr));
	server_addr.sin_family = AF_INET;
	server_addr.sin_port = htons(serverPort);
	server_addr.sin_addr.s_addr = inet_addr(serverIP);
	err = connect(m_socketfd, (struct sockaddr *)&server_addr, sizeof(server_addr));
	if(err != 0)
		return false;

	return true;
}

bool socketConn::fetchConnection()
{
	struct sockaddr_in client_addr;
	memset(&client_addr, 0, sizeof(client_addr));
	socklen_t length = sizeof(client_addr);
	int conn = accept(m_socketfd, (struct sockaddr*)&client_addr, &length);
	if(conn < 0)
		return false;
	else
		m_connection = conn;
}

int socketConn::receiveData(char *receivedData, size_t buffSize)
{
	memset(receivedData,0,sizeof(receivedData));
	if(m_connection > 0)
	{
		int receivedLength = recv(m_connection, receivedData, buffSize, 0);
		return receivedLength;
	}
	else
		return -1;
}

int socketConn::sendData(const char *sendData)
{
	return send(m_socketfd, sendData, strlen(sendData), 0);
}

int socketConn::writeOnConnection(const void *buf, size_t size)
{
	return write(m_connection, buf, size);
}

int socketConn::readOnConnection(void *buff, size_t size)
{
	return read(m_connection, buff, size);
}

int socketConn::readFromSocket(void *buff, size_t size)
{
	return read(m_socketfd, buff, size);
}

void socketConn::closeSocket()
{
	close(m_socketfd);
}

void socketConn::closeConnection()
{
	close(m_connection);
}
