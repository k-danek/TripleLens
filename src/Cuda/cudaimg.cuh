#include <iostream>
#include <vector>
#include <complex>
#include <math.h>
#include <unordered_map>
#include <limits>
#include <list>
#include <typeinfo>
#include <functional>

#include <thrust/complex.h>
#include <thrust/device_reference.h>

#include "lens.h"
#include "laguerre.cuh"

#ifndef GRIDLINE
#define GRIDLINE
/* xright - leftmost x coordinate occupied by the lensed image
   xleft - rightmost x coordinate occupied by the lensed image
*/
struct GridLine
{
  thrust::complex<double> start;
  thrust::complex<double> end;
  unsigned int steps;
};
#endif


#ifndef POINTIMAGE_CUH
#define POINTIMAGE_CUH

class ImgPointCUDA
{
  public:

    ImgPointCUDA(double a,
                 double b,
                 double th,
                 double m2,
                 double m3,
                 double sourceSize
                );

    // Binary only version
    ImgPointCUDA(double a,
                 double th,
                 double m,
                 double sourceSize
                );


    void freeAll();

    double syncAndReturn(int lcStep);

    void trigger(double      sourcePosX,
                 double      sourcePosY);

    // Image position caulculation
    void getRoots(bool forceNewRoots,
                  bool isBinaryLens); 

    // update and return images in one functional call
    void getRootsPrecalculated(bool forceNewRoots);

    void allocateHost(int size);
    void allocateCuda();
    void setConstantPars();

  private:
    const int _numOfBlocks = 128;

    double *_tempParams;
    double *_ampsHost, *_ampsDeviceA, *_ampsDeviceB, *_ampsDeviceC;
    GridLine *_trajectoryHost, *_trajectoryDeviceA, *_trajectoryDeviceB, *_trajectoryDeviceC;
    double _a, _b, _th, _m2, _m3, _sourceSize;
    void _setConstantPar();
    void _invokeKernelDouble(double* amps, std::vector<GridLine> trajectories); 
    void _invokeKernelTriple(double* amps, std::vector<GridLine> trajectories);
};


__global__
void getAmps(double* amps, GridLine* trajectories);

__device__
void getImageMask(thrust::complex<double>* coeffs);

__device__
void getCoeffs(thrust::complex<double>* coeffs, thrust::complex<double> zeta);

__device__
void getCoeffsBinOpt(thrust::complex<double>* coeffs, thrust::complex<double> zeta);


#endif
