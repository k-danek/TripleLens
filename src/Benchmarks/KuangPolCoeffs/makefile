CC=g++
#https://stackoverflow.com/questions/7420665/what-does-gccs-ffast-math-actually-do
#https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html
#-Ofast flags and their impact:
# 0) all seems to be under -ffast-math
# 1) -funsafe-math-optimizations makes the performance of my code on the same level as Kuang's
# 2) -ffinite-math-only has no impact at all
# 3) -fno-rounding-math no impact
# 4) -fno-signaling-nans no impact
# 5) -fcx-limited-range some impact but not too much
# 6) -fexcess-precision=fast no measurable impact, probably already in place.
# 7) -freciprocal-math seems to mess up accuracy at 15-16th figure.
# 9) -fassociative-math -fno-signed-zeros -fno-trapping-math combination leads
#     to correct result while leveling the performance.
# Also see https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html


#OPTFLAGS= -O3 -fassociative-math -fno-signed-zeros -fno-trapping-math
#OPTFLAGS= -O3
OPTFLAGS= -Ofast
CFLAGS= -Wl,--no-undefined -fPIC $(OPTFLAGS) -g -std=c++17

SOURCE_DIR=../..
INC_LENS=$(SOURCE_DIR)/LensCore
INC_LAGUERRE=$(SOURCE_DIR)/Laguerre
INC_IMG=$(SOURCE_DIR)/Images

INCLUDES= -I. -I$(INC_LENS) -I$(INC_LENS) -I$(INC_LAGUERRE) -I$(INC_IMG)

kuangtest:  laguerre.o lens.o imgpoint.o kuangtest.o kuangcoef.o
	$(CC) $(INCLUDES) $(CFLAGS) -o kuangtest ./kuangtest.o ./kuangcoef.o ./imgpoint.o ./laguerre.o ./lens.o

kuangcoef.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) -o kuangcoef.o imgpointcoeffkuang.cc

kuangtest.o:
	$(CC) -c $(INCLUDES) $(CFLAGS) -o kuangtest.o test_kuang.cc

imgpoint.o: $(INC_IMG)/imgpoint.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o ./imgpoint.o $(INC_IMG)/imgpoint.cc

laguerre.o: $(INC_LAGUERRE)/laguerre.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o ./laguerre.o $(INC_LAGUERRE)/laguerre.cc

lens.o: $(INC_LENS)/lens.h
	$(CC) -c $(INCLUDES) $(CFLAGS) -o ./lens.o $(INC_LENS)/lens.cc

-include $(INCLUDES)

clean:
	rm *.o kuangtest

