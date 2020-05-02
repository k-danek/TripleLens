CC=g++
# Note flag -Ofast as it made execution cca 10 times faster
# Flag -fPIC makes o files usable in the shared library
CFLAGS= -Wl,--no-undefined -fPIC -Ofast -g -o
CFLAGS_SHARED= -shared
INC_CCC=./CCC
INC_LENS=./LensCore
INC_LAGUERRE=./Laguerre
BUILD_TARGET=./bin
INCLUDES= -I. -I$(INC_CCC) -I$(INC_LENS) -I$(INC_LAGUERRE) -I$(BUILD_TARGET)

ccc_test: main.o lens.o laguerre.o ccc.o
	$(CC) $(CFLAGS) $(BUILD_TARGET)/ccc_test $(BUILD_TARGET)/main.o $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/laguerre.o $(BUILD_TARGET)/lens.o -lm -ltiff

ccc.o: $(INC_CCC)/ccc.cc $(INC_CCC)/ccc.h lens.o laguerre.o
	$(CC) -c $(INCLUDES) $(CFLAGS) $(CFLAGS) $(BUILD_TARGET)/ccc.o $(INC_CCC)/ccc.cc

ccc.so: $(BUILD_TARGET)/lens.o laguerre.o ccc.o
	$(CC) $(CFLAGS_SHARED) $(CFLAGS) $(BUILD_TARGET)/ccc.so $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/laguerre.o $(BUILD_TARGET)/lens.o 

laguerre.o: $(INC_LAGUERRE)/laguerre.cc $(INC_LAGUERRE)/laguerre.h
	$(CC) -c $(INCLUDES) $(CFLAGS) $(BUILD_TARGET)/laguerre.o $(INC_LAGUERRE)/laguerre.cc

lens.o: $(INC_LENS)/lens.cc $(INC_LENS)/lens.h
	$(CC) -c $(INCLUDES) $(CFLAGS) $(BUILD_TARGET)/lens.o $(INC_LENS)/lens.cc

main.o: main.cc ccc.o lens.o laguerre.o ccc.so
	$(CC) -c $(INCLUDES) $(CFLAGS) $(BUILD_TARGET)/main.o main.cc

-include $(INCLUDES)

clean:
	rm $(BUILD_TARGET)/*.o $(BUILD_TARGET)/*.so $(BUILD_TARGET)/ccc_test

