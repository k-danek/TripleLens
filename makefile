CC=g++
# Note flag -Ofast as it made execution cca 10 times faster than -O2 on CCC 
# Flag -fPIC makes o files usable in the shared library
CFLAGS= -Wl,--no-undefined -fPIC -Ofast -g -std=c++17
CFLAGS_SHARED= -shared
INC_CCC=./CCC
INC_LENS=./LensCore
INC_LAGUERRE=./Laguerre
INC_LC=./LightCurves
INC_IMG=./Images
INC_UTILS=./Utils
BUILD_TARGET=./bin
INCLUDES= -I. -I$(INC_CCC) -I$(INC_LENS) -I$(INC_LAGUERRE) -I$(INC_IMG) -I$(INC_UTILS) -I$(INC_LC) -I$(BUILD_TARGET)

all: imgpoint.so ccc.so lcirs.so ccc_test

# Executables
ccc_test: main.o ccc.so imgpoint.so lcirs.so  liblaguerre.a 
	$(CC) $(CFLAGS) -o $(BUILD_TARGET)/ccc_test $(BUILD_TARGET)/main.o $(BUILD_TARGET)/amoeba.o $(BUILD_TARGET)/lcbase.o $(BUILD_TARGET)/lcirs.o $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/lens.o $(BUILD_TARGET)/liblaguerre.a

# Shared libraries
lcirs.so: amoeba.o lens.o ccc.o imgpoint.o lcbase.o lcirs.o liblaguerre.a 
	$(CC) $(CFLAGS_SHARED) $(CFLAGS) -o $(BUILD_TARGET)/lcirs.so $(BUILD_TARGET)/amoeba.o $(BUILD_TARGET)/lens.o $(BUILD_TARGET)/ccc.o  $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/lcbase.o $(BUILD_TARGET)/lcirs.o $(BUILD_TARGET)/liblaguerre.a 

imgpoint.so: lens.o imgpoint.o liblaguerre.a  
	$(CC) $(CFLAGS_SHARED) $(CFLAGS) -o $(BUILD_TARGET)/imgpoint.so $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/lens.o $(BUILD_TARGET)/liblaguerre.a

ccc.so: lens.o ccc.o liblaguerre.a 
	$(CC) $(CFLAGS_SHARED) $(CFLAGS) -o $(BUILD_TARGET)/ccc.so $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/lens.o $(BUILD_TARGET)/liblaguerre.a 

# Static libraries
liblaguerre.a: laguerre.o
	ar -cvq $(BUILD_TARGET)/liblaguerre.a $(BUILD_TARGET)/laguerre.o

# Object files
main.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/main.o main.cc

lcirs.o: $(INC_LC)/lcirs.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/lcirs.o $(INC_LC)/lcirs.cc 

amoeba.o: $(INC_UTILS)/amoeba.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/amoeba.o $(INC_UTILS)/amoeba.cc 

lcbase.o: $(INC_LC)/lcbase.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/lcbase.o $(INC_LC)/lcbase.cc 

imgpoint.o: $(INC_IMG)/imgpoint.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/imgpoint.o $(INC_IMG)/imgpoint.cc 

ccc.o: $(INC_CCC)/ccc.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/ccc.o $(INC_CCC)/ccc.cc 

laguerre.o: $(INC_LAGUERRE)/laguerre.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/laguerre.o $(INC_LAGUERRE)/laguerre.cc

lens.o: $(INC_LENS)/lens.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/lens.o $(INC_LENS)/lens.cc

-include $(INCLUDES)

clean:
	rm $(BUILD_TARGET)/*.o $(BUILD_TARGET)/*.so $(BUILD_TARGET)/*.a $(BUILD_TARGET)/ccc_test

