CC=g++
CFLAGS= -Wall -g -o
INCLUDES= -I. -I../LensCore -I../Laguerre -I../Images -I../CCC -I../Utils

#lcbase.o:
#	$(CC) -c $(INCLUDES) $(CFLAGS) lcbase.o lcbase.cc

lcirs.o: lcbase.o
	$(CC) -c $(INCLUDES) $(CFLAGS) lcirs.o lcirs.cc

-include $(INCLUDES)

clean:
	rm *.o

