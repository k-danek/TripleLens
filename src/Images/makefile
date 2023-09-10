CC=g++
#CFLAGS= -Wall -g -fPIC -o
INCLUDES= -I. -I../LensCore -I../Laguerre
CFLAGS= -Wl,--no-undefined -fPIC -Ofast -g -std=c++17
#CFLAGS= -Wl,--no-undefined -g
#g++ -shared -W1,-soname,ccc.so -o ccc.so ccc.o
#g++ -shared -o ccc.so ccc.o

SOURCE_DIR=../
INC_LENS=$(SOURCE_DIR)/LensCore
INC_LAGUERRE=$(SOURCE_DIR)/Laguerre

kuangtest: imgpoint.o laguerre.o lens.o kuangtest.o kuangcoef.o
	$(CC) $(INCLUDES) $(CFLAGS) -o kuangtest ./kuangtest.o ./kuangcoef.o ./imgpoint.o ./laguerre.o ./lens.o

kuangcoef.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) -o kuangcoef.o imgpointcoeffkuang.cc

kuangtest.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) -o kuangtest.o test_kuang.cc

imgpoint.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) -o imgpoint.o imgpoint.cc

laguerre.o: 
	$(CC) -c $(INCLUDES) $(CFLAGS) -o ./laguerre.o $(INC_LAGUERRE)/laguerre.cc

lens.o: 
	$(CC) -c $(INCLUDES) $(CFLAGS) -o ./lens.o $(INC_LENS)/lens.cc

-include $(INCLUDES)

clean:
	rm *.o kuangtest

