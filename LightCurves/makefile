CC=g++
CFLAGS= -Wall -g -fPIC -o
INCLUDES= -I. -I../LensCore -I../Laguerre -I../Images -I../CCC -I../Utils

#g++ -shared -W1,-soname,ccc.so -o ccc.so ccc.o
#g++ -shared -o ccc.so ccc.o

#lcbase.o:
#	$(CC) -c $(INCLUDES) $(CFLAGS) lcbase.o lcbase.cc

lcirs.o: lcbase.o
	$(CC) -c $(INCLUDES) $(CFLAGS) lcirs.o lcirs.cc

-include $(INCLUDES)

clean:
	rm *.o

