CC=g++
# Note flag -Ofast as it made execution cca 10 times faster than -O2 on CCC 
# Flag -fPIC makes o files usable in the shared library
CFLAGS= -Wl,--no-undefined -fPIC -Ofast -g -std=c++17
CPROFFLAGS= -pg -no-pie -fno-builtin
CFLAGS_SHARED= -shared

CFCOMPATIBLE=

# CUDA root dir:
CUDA_ROOT_DIR=/usr/local/cuda
# CUDA library directory:
CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
# CUDA include directory:
CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include

NVCC=nvcc
# Note flag -Ofast as it made execution cca 10 times faster than -O2 on CCC 
# Flag -fPIC makes o files usable in the shared library
CUDAFLAGS= --relocatable-device-code=true --compiler-options '-fPIC'


SOURCE_DIR=./src
INC_CCC=$(SOURCE_DIR)/CCC
INC_LENS=$(SOURCE_DIR)/LensCore
INC_LAGUERRE=$(SOURCE_DIR)/Laguerre
INC_LC=$(SOURCE_DIR)/LightCurves
INC_IMG=$(SOURCE_DIR)/Images
INC_CUDA=$(SOURCE_DIR)/Cuda
INC_UTILS=$(SOURCE_DIR)/Utils
BUILD_TARGET=./bin
INCLUDES= -I. -I$(INC_CCC) -I$(INC_LENS) -I$(INC_LAGUERRE) -I$(INC_IMG) -I$(INC_UTILS) -I$(INC_CUDA) -I$(INC_LC) $(CUDA_LIB_DIR) $(CUDA_INC_DIR) -I$(BUILD_TARGET)

all: libimgpoint.so libccc.so liblcirs.so ccc_test

profile: libimgpoint.so libccc.so liblcirs.so ccc_profile

# Executables
ccc_test: main.o libccc.so libimgpoint.so liblcirs.so laguerre.o cudapointcollector.o cudalink.o
	$(CC) $(CFLAGS) -o $(BUILD_TARGET)/ccc_test $(BUILD_TARGET)/main.o $(BUILD_TARGET)/amoeba.o $(BUILD_TARGET)/lcbase.o $(BUILD_TARGET)/lcirs.o $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/lens.o $(BUILD_TARGET)/laguerre.o $(BUILD_TARGET)/cudapointcollector.o $(BUILD_TARGET)/cudalink.o $(BUILD_TARGET)/cudairs.o -lcudadevrt -lcudart

ccc_profile: main.o ccc.so libimgpoint.so liblcirs.so laguerre.o cudapointcollector.o cudalink.o
	$(CC) $(CFLAGS) $(CPROFFLAGS) -o $(BUILD_TARGET)/ccc_profile $(BUILD_TARGET)/main.o $(BUILD_TARGET)/amoeba.o $(BUILD_TARGET)/lcbase.o $(BUILD_TARGET)/lcirs.o $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/lens.o $(BUILD_TARGET)/laguerre.o $(BUILD_TARGET)/cudapointcollector.o $(BUILD_TARGET)/cudalink.o $(BUILD_TARGET)/cudairs.o -lcudadevrt -lcudart

# Shared libraries
liblcirs.so: amoeba.o lens.o ccc.o imgpoint.o lcbase.o lcirs.o laguerre.o cudapointcollector.o
	$(CC) $(CFLAGS_SHARED) $(CFLAGS) -o $(BUILD_TARGET)/liblcirs.so $(BUILD_TARGET)/amoeba.o $(BUILD_TARGET)/lens.o $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/lcbase.o $(BUILD_TARGET)/lcirs.o $(BUILD_TARGET)/laguerre.o $(BUILD_TARGET)/cudapointcollector.o $(BUILD_TARGET)/cudalink.o $(BUILD_TARGET)/cudairs.o -lcudadevrt -lcudart

libimgpoint.so: lens.o imgpoint.o laguerre.o  
	$(CC) $(CFLAGS_SHARED) $(CFLAGS) -o $(BUILD_TARGET)/libimgpoint.so $(BUILD_TARGET)/imgpoint.o $(BUILD_TARGET)/lens.o $(BUILD_TARGET)/laguerre.o

libccc.so: lens.o ccc.o laguerre.o 
	$(CC) $(CFLAGS_SHARED) $(CFLAGS) -o $(BUILD_TARGET)/libccc.so $(BUILD_TARGET)/ccc.o $(BUILD_TARGET)/lens.o $(BUILD_TARGET)/laguerre.o 

# Object files
main.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/main.o $(SOURCE_DIR)/main.cc

lcirs.o: $(INC_LC)/lcirs.h cudapointcollector.o
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/lcirs.o $(INC_LC)/lcirs.cc

amoeba.o: $(INC_UTILS)/amoeba.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/amoeba.o $(INC_UTILS)/amoeba.cc 

cudapointcollector.o: $(INC_CUDA)/cudapointcollector.h cudalink.o
	$(CC) -c $(INCLUDES) $(CFLAGS) -o $(BUILD_TARGET)/cudapointcollector.o $(INC_CUDA)/cudapointcollector.cc

cudalink.o: cudairs.o
	$(NVCC) -dlink --compiler-options '-fPIC' -o $(BUILD_TARGET)/cudalink.o $(BUILD_TARGET)/cudairs.o -lcudadevrt -lcudart

cudairs.o: $(INC_CUDA)/cudairs.cu
	$(NVCC) -c $(CUDAFLAGS) -dc -o $(BUILD_TARGET)/cudairs.o $(INC_CUDA)/cudairs.cu

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
	rm $(BUILD_TARGET)/*.o $(BUILD_TARGET)/*.so $(BUILD_TARGET)/*.a $(BUILD_TARGET)/ccc_test $(BUILD_TARGET)/ccc_profile


