#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <cuda_profiler_api.h>

#include "cudaimg.cuh"

__constant__ double params[8];


ImgPointCUDA::ImgPointCUDA(double a,
                           double b,
                           double th,
                           double m2,
                           double m3,
                           double sourceSize
                          ): _a(a), _b(b), _th(th), _m2(m2), _m3(m3), _sourceSize(sourceSize) 
{
    _setConstantPar();
};

// Version for double lens
ImgPointCUDA::ImgPointCUDA(double a,
                           double th,
                           double m,
                           double sourceSize
                          ):_a(a), _m2(m), _sourceSize(sourceSize) 
{
    _setConstantPar();
};    


void ImgPointCUDA::_setConstantPar()
{
  double _tempParams[8];

  _tempParams[0] = double(1.0-_m2-_m3);  // m1
  _tempParams[1] = double(_m2);         // m2
  _tempParams[2] = double(_m3);         // m3
  _tempParams[3] = double(_a);           // z2x
  _tempParams[4] = double(0.0);          // z2y
  _tempParams[5] = double(_b*cos(_th));  // z3x
  _tempParams[6] = double(_b*sin(_th));  // z3y
  _tempParams[7] = double(_sourceSize); // sourceSize

  cudaMemcpyToSymbol(params, _tempParams, sizeof(double)*8);
}

void ImgPointCUDA::allocateCuda()
{
  // Allocate bufferts for CUDA. The size is fixed and set up in constants.

  clock_t beginTime, endTime; 

  //beginTime = clock();
  // Allocation of device buffers
  cudaMalloc((void**)&_trajectoryDeviceA, _numOfBlocks*sizeof(GridLine));

  // Allocatin device buffers
  cudaMalloc((void**)&_ampsDeviceA, _numOfBlocks*sizeof(double));
  //endTime = clock();
  //_gpuMallocTime += double(endTime-beginTime);
};


void ImgPointCUDA::freeAll()
{
  //  double *_ampsHost, *_ampsDeviceA, *_ampsDeviceB, *_ampsDeviceC;
  //  GridLine *_trajectoryHost, *_trajectoryDeviceA, *_trajectoryDeviceB, *_trajectoryDeviceC;

  if(_ampsHost != nullptr) {
    cudaFreeHost(_ampsHost);
    _ampsHost = nullptr;
  }

  if(_ampsDeviceA != nullptr) {
    cudaFreeHost(_ampsDeviceA);
    _ampsDeviceA = nullptr;
  }

  if(_trajectoryHost != nullptr) {
    cudaFreeHost(_trajectoryHost);
    _trajectoryHost = nullptr;
  }

  if(_trajectoryDeviceA != nullptr) {
    cudaFreeHost(_trajectoryDeviceA);
    _trajectoryDeviceA = nullptr;
  }

  if(_tempParams != nullptr) {
    free(_tempParams);
    _tempParams = nullptr;
  }
}

void ImgPointCUDA::allocateHost(int size)
{
  // size should be something like numberOfNodesExtended
  //_numberOfNodesBufferSize = size;

  clock_t beginTime, endTime; 

  //beginTime = clock();
  cudaHostAlloc((void**)& _trajectoryHost,
                sizeof(GridLine)*size,
                cudaHostAllocDefault);
  //endTime = clock();
  //_gpuMallocTime += double(endTime-beginTime);

  // Allocating in pinned memory
  // Always check whether these need to be initialized correctly. 
  cudaHostAlloc((void**)&_ampsHost,sizeof(double)*size, cudaHostAllocDefault);

  //beginTime = clock();
  _tempParams = (double*)malloc(sizeof(double)*8);
  //endTime = clock();
  //_gpuConstMemTime += double(endTime-beginTime);
};

