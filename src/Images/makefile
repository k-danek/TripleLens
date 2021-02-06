CC=g++
CFLAGS= -Wall -g -fPIC -o
INCLUDES= -I. -I../LensCore -I../Laguerre

#g++ -shared -W1,-soname,ccc.so -o ccc.so ccc.o
#g++ -shared -o ccc.so ccc.o


imgpoint.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) imgpoint.o imgpoint.cc


-include $(INCLUDES)

clean:
	rm *.o

