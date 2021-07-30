
ifeq ($(shell getconf LONG_BIT), 32)
        WORDSIZE = 32
else
        WORDSIZE = 64
endif

CC = g++

CFLAGS = -I../include/ \

OBJS = atom.o bond.o fragment.o \
     molecule.o markush.o utility.o interface.o socketConn.o
OBJS_ENGINE = main_molEngine.o
OBJS_SERVER = main_server.o
OBJS_CLIENT = main_client.o

BINNAME = molEngine
LIBNAME = libmol.a
SERVER = mols
CLIENT = molc

RM = rm -f

engine: $(BINNAME) $(LIBNAME)
	cp $(LIBNAME) ../lib$(WORDSIZE)/
socket: $(SERVER) $(CLIENT)

$(LIBNAME): $(OBJS)
	ar rcv $@ $^

$(BINNAME): $(OBJS) $(OBJS_ENGINE)
	$(CC) -g $(CFLAGS) -o $(BINNAME) $(OBJS) $(OBJS_ENGINE)
$(SERVER): $(OBJS) $(OBJS_SERVER)
	$(CC) -g $(CFLAGS) -o $(SERVER) $(OBJS) $(OBJS_SERVER)
$(CLIENT): $(OBJS) $(OBJS_CLIENT)
	$(CC) -g $(CFLAGS) -o $(CLIENT) $(OBJS) $(OBJS_CLIENT)

$(OBJS):%.o:%.cpp
	$(CC) -g -c $(CFLAGS) $^ -o $@
$(OBJS_ENGINE):%.o:%.cpp
	$(CC) -g -c $(CFLAGS) $^ -o $@
$(OBJS_SERVER):%.o:%.cpp
	$(CC) -g -c $(CFLAGS) $^ -o $@
$(OBJS_CLIENT):%.o:%.cpp
	$(CC) -g -c $(CFLAGS) $^ -o $@

clean:
	$(RM) $(BINNAME) $(SERVER) $(CLIENT) *.o *.a ../lib$(WORDSIZE)/*.a
