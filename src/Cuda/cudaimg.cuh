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

//#include "lens.h"
#include "laguerre.cuh"

#ifndef GRIDLINE
#define GRIDLINE
/* xright - leftmost x coordinate occupied by the lensed image
   xleft - rightmost x coordinate occupied by the lensed image
*/
struct GridLine
{
  // Actual first point on the trajectory
  thrust::complex<double> start;
  // Actual end on the rajectory
  thrust::complex<double> end;
  //unsigned int steps; // commenting out as steps need to be _threadsPerBlock
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

    void trigger(std::vector<GridLine> trajectories);
    
    // Version that uses locally stored trajectories
    void trigger();

    std::vector<std::vector<float>> syncAndReturn();

    // Image position caulculation
    void getRoots(bool forceNewRoots,
                  bool isBinaryLens); 

    // update and return images in one functional call
    void getRootsPrecalculated(bool forceNewRoots);

    void allocateHost();
    void allocateCuda();
    void setConstantPars();

    std::vector<GridLine> getPolarTrajectories(std::vector<double> impactParameters,
                                               std::vector<double> angles);

  private:
    // correspond to number to trajectories to run in one call
    const int _numOfBlocks = 16;

    // correspond to number of steps per trajectory 
    const int _threadsPerBlock = 128;

    double *_tempParams;
    float  *_ampsHost, *_ampsDeviceA, *_ampsDeviceB, *_ampsDeviceC;
    GridLine *_trajectoryHost, *_trajectoryDeviceA, *_trajectoryDeviceB, *_trajectoryDeviceC;
    std::vector<GridLine> _storedTrajectories;
    double _a, _b, _th, _m2, _m3, _sourceSize;
    void _setConstantPar();
    void _invokeKernelDouble(double* amps, std::vector<GridLine> trajectories); 
    void _invokeKernelTriple();
};
#endif

__global__
void getAmps(double* amps, GridLine* trajectories);

__device__
void rootsToAmps(double* amps,
                 thrust::complex<double>* roots,
                 const thrust::complex<double> &sourcePos,
                 const thrust::complex<double> &z2,
                 const thrust::complex<double> &z3);

__device__
void getImageMask(thrust::complex<double>* coeffs);

__device__
void getCoeffs(thrust::complex<double>* coeffs, thrust::complex<double> zeta);

__device__
void getCoeffsBinOpt(thrust::complex<double>* coeffs, thrust::complex<double> zeta);
