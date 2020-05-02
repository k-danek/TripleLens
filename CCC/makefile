CC=g++
CFLAGS= -Wall -g -fPIC -o
INCLUDES= -I. -I../LensCore -I../Laguerre

#g++ -shared -W1,-soname,ccc.so -o ccc.so ccc.o
#g++ -shared -o ccc.so ccc.o


#ccc.o: ccc.cc ccc.h
#	$(CC) -c $(INCLUDES) $(CFLAGS) ccc.o ccc.cc

ccc.so: ccc.o
	$(CC) -shared -Wl,--no-undefined -o ccc.so ccc.o

-include $(INCLUDES)

clean:
	rm *.o

