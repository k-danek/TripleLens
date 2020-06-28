CC=g++
# Note flag -Ofast as it made execution cca 10 times faster
# Flag -fPIC makes o files usable in the shared library
CFLAGS= -Wl,--no-undefined -fPIC -Ofast -g -std=c++11
CFLAGS_SHARED= -shared
INC_CCC=./CCC
INC_LENS=./LensCore
INC_LAGUERRE=./Laguerre
INC_LC=./LightCurves
INC_IMG=./Images
INC_UTILS=./Utils
BUILD_TARGET=./bin
INCLUDES= -I. -I$(INC_CCC) -I$(INC_LENS) -I$(INC_LAGUERRE) -I$(INC_IMG) -I$(INC_UTILS) -I$(INC_LC) -I$(BUILD_TARGET)

ccc_test: main.o ccc.so imgpoint.so lcirs.so
	$(CC) $(CFLAGS) -o $(BUILD_TARGET)/ccc_test $(BUILD_TARGET)/main.o $(BUILD_TARGET)/amoeba.o $(BUILD_TARGET)/lcbase.o $(BUILD_TARGET)/lcirs.o $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/liblaguerre.a $(BUILD_TARGET)/lens.o

lcbase.o: lens.o imgpoint.o
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/lcbase.o $(BUILD_TARGET)/imgpoint.o $(INC_LC)/lcbase.cc 

lcirs.o: lens.o ccc.o amoeba.o imgpoint.o
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/lcirs.o $(BUILD_TARGET)/amoeba.o $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/lens.o $(INC_LC)/lcirs.cc 

lcirs.so: lens.o ccc.o amoeba.o imgpoint.o lcirs.o
	$(CC) -c $(CFLAGS_SHARED) $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/lcirs.so $(BUILD_TARGET)/amoeba.o $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/lens.o

amoeba.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/amoeba.o $(INC_UTILS)/amoeba.cc 

imgpoint.o: $(INC_IMG)/imgpoint.cc $(INC_IMG)/imgpoint.h lens.o laguerre.o
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/liblaguerre.a $(INC_IMG)/imgpoint.cc 

imgpoint.so: imgpoint.o liblaguerre.a imgpoint.o 
	$(CC) -c $(CFLAGS_SHARED) $(CFLAGS) -o $(BUILD_TARGET)/imgpoint.so $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/liblaguerre.a $(BUILD_TARGET)/lens.o

ccc.o: $(INC_CCC)/ccc.cc $(INC_CCC)/ccc.h lens.o liblaguerre.a
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/liblaguerre.a $(INC_CCC)/ccc.cc 

ccc.so: lens.o liblaguerre.a ccc.o
	$(CC) $(CFLAGS_SHARED) $(CFLAGS) -o $(BUILD_TARGET)/ccc.so $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/liblaguerre.a $(BUILD_TARGET)/lens.o 

liblaguerre.a: laguerre.o
	ar -cvq $(BUILD_TARGET)/liblaguerre.a $(BUILD_TARGET)/laguerre.o

laguerre.o: $(INC_LAGUERRE)/laguerre.cc $(INC_LAGUERRE)/laguerre.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/laguerre.o $(INC_LAGUERRE)/laguerre.cc

lens.o: $(INC_LENS)/lens.cc $(INC_LENS)/lens.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/lens.o $(INC_LENS)/lens.cc

main.o: main.cc ccc.o lens.o liblaguerre.a imgpoint.o lcbase.o lcirs.o
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/main.o main.cc

-include $(INCLUDES)

clean:
	rm $(BUILD_TARGET)/*.o $(BUILD_TARGET)/*.so $(BUILD_TARGET)/*.a $(BUILD_TARGET)/ccc_test

