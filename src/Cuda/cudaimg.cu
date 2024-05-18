#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <cuda_profiler_api.h>

#include "cudaimg.cuh"

__constant__ double paramsImg[8];


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


std::vector<GridLine> ImgPointCUDA::getPolarTrajectories(
                      std::vector<double> impactParameters,
                      std::vector<double> angles)
{
  // for the treajectory x = cos(alpha)*t-q*sin(alpha)
  //                     y = sin(alpha)*t+q*sin(alpha)
  //  t goes from -1 to 1 to cover 2 Einstein Radii. That is about to get changed 
  
  int trajectoryCounter = 0;
  std::vector<GridLine> outputTrajectories;
  const double iniTime = -1.0;
  const double finTime =  1.0;
  
  for(auto q: impactParameters)
  {
    for(auto alpha: angles)
    {
      complex<double> ini(cos(alpha)*iniTime-q*sin(alpha), sin(alpha)*iniTime+q*cos(alpha));
      complex<double> fin(cos(alpha)*finTime-q*sin(alpha), sin(alpha)*finTime+q*cos(alpha));
      GridLine trajectory = {ini, fin};
      outputTrajectories.push_back(trajectory);
      trajectoryCounter++;
      if(trajectoryCounter == _numOfBlocks)
      {
        break;
      }
    }
  }

  return outputTrajectories;
}




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
  _tempParams[7] = double(_sourceSize*_sourceSize); // sourceSize squared

  cudaMemcpyToSymbol(paramsImg, _tempParams, sizeof(double)*8);
}

void ImgPointCUDA::allocateCuda()
{
  // Allocate bufferts for CUDA. The size is fixed and set up in constants.

  clock_t beginTime, endTime; 

  //beginTime = clock();
  // Allocation of device buffers
  cudaMalloc((void**)&_trajectoryDeviceA, _numOfBlocks*sizeof(GridLine));

  // Allocatin device buffers
  cudaMalloc((void**)&_ampsDeviceA, sizeof(float)*_numOfBlocks*_threadsPerBlock);
  //endTime = clock();
  //_gpuMallocTime += double(endTime-beginTime);
};

void ImgPointCUDA::allocateHost()
{
  // size should be something like numberOfNodesExtended
  clock_t beginTime, endTime; 
  
  // Allocate pinned memory for input trajectories
  cudaHostAlloc((void**)& _trajectoryHost,
                sizeof(GridLine)*_numOfBlocks,
                cudaHostAllocDefault);

  // Allocating pinned memory for output amplifications
  // Always check whether these need to be initialized correctly. 
  cudaHostAlloc((void**)&_ampsHost,
                sizeof(float)*_numOfBlocks*_threadsPerBlock,
                cudaHostAllocDefault);

  _tempParams = (double*)malloc(sizeof(double)*8);
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

// Variable threads per points
__global__
void trajectoriesToAmps(GridLine* gridLine,
                        float*    amps)
{
  const thrust::complex<double> z2 = thrust::complex<double>(paramsImg[8],paramsImg[9]);
  const thrust::complex<double> z3 = thrust::complex<double>(paramsImg[10],paramsImg[11]);

  // Taking block index as index for the of one trajectory
  GridLine locLine = gridLine[blockIdx.x];
  const thrust::complex<double> start = locLine.start;
  const thrust::complex<double> stop = locLine.end;
  
  // Steps now serve as threads-per-block as block index is per one trajectory.
  const double stepRatio = __int2double_rn(threadIdx.x)/__int2double_rn(blockDim.x- 1);
  // actual index of a thread
  //int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

  // Make sure these are already shifted to the bottom left corner of a pixel! 
  //
  //threadIdx.x % subgridSize; that is subgrid x
  //threadIdx.x / subgridSize; that is subgrid y
  //double xShift = paramsImg[13] + __int2double_rn(threadIdx.x % sgSize[0])*paramsImg[12];
  //double yShift = paramsImg[14] + __ll2double_rn(gridY)*paramsImg[7] + __int2double_rn(threadIdx.x / sgSize[0])*paramsImg[12];
  
  const thrust::complex<double> sourcePos = start*(1.0-stepRatio)+stop*stepRatio;

  thrust::complex<double> coeffs[11];

  getCoeffs(coeffs, sourcePos);

  thrust::complex<double> roots[10];

  // roots, coeffs, order, max iterations
  solveRootsCUDA(roots, coeffs, 10, 30);

  double imgAmps[10];

  // amps, roots
  rootsToAmps(imgAmps, roots, sourcePos, z2, z3);

  double tempAmps = 0.0;


  // In the first stage stop after coeff generation
  for(unsigned int i = 0; i <= 10; i++)
  {
    tempAmps += imgAmps[i]; 
  }

  // For sm_30 there is no atomicAdd that would accept doubles.
  // Change this if you evet lay your hands on sm_60.
  //atomicAdd(&amps[blockIdx.x], __double2float_rn(tempAmps));
  atomicAdd(&amps[blockIdx.x * blockDim.x + threadIdx.x], __double2float_rn(tempAmps));
  //atomicAdd(&amps[blockIdx.x * blockDim.x + threadIdx.x], __double2float_rn(coeffs[0].real()));
  //atomicAdd(&amps[blockIdx.x * blockDim.x + threadIdx.x], 1.0*threadIdx.x);
  //atomicAdd(&amps[blockIdx.x * blockDim.x + threadIdx.x], sourcePos.real());
};


void ImgPointCUDA::_invokeKernelTriple()
{
    // Num of blocks ~ num of trajectories
    // threadsPerBlock ~ steps per a trajectory
    trajectoriesToAmps<<<_numOfBlocks, _threadsPerBlock, 0>>>(_trajectoryDeviceA,
                                                              _ampsDeviceA);
};

std::vector<std::vector<float>> ImgPointCUDA::syncAndReturn()
{
  std::cout << "syncAndReturn visited \n";
  std::vector<std::vector<float>> amps(_numOfBlocks, std::vector<float>(_threadsPerBlock, 0.0));

  std::cout << "syncAndReturn sizes; amps: " 
            << amps.size() << ", each amp: "
            << amps[0].size() << ", _ampHost: "
            << sizeof(_ampsHost)/sizeof(float) << "\n"
            << "some ampHost values:" << " " << _ampsHost[0] << " " << _ampsHost[100] << "\n";

  // Transfer the data from the 1D array to the 2D vector
  for (int i = 0; i < _numOfBlocks; ++i)
  {
    std::copy(_ampsHost+i*_threadsPerBlock, _ampsHost+(i+1)*_threadsPerBlock, amps[i].begin());
  }

  freeAll();

  std::cout << "syncAndReturn all moved \n";
  return amps;
}

void ImgPointCUDA::trigger(std::vector<GridLine> trajectories)
{
  _storedTrajectories = trajectories;
  trigger();
  return;
}

void ImgPointCUDA::trigger()
{
  // I might easily run out of available blocks per grid.
  // Supposed size of the number of blocks is 65535.
  // https://en.wikipedia.org/wiki/Thread_block_(CUDA_programming)#Dimensions
  // Please note that Device query claims following:
  // Max dimension size of a thread block (x,y,z): (1024, 1024, 64)
  // Max dimension size of a grid size    (x,y,z): (2147483647, 65535, 65535)
  // It equals number of blocks used for a single GPU execution. 
  //const int numOfBlocks = 128;
  
  //// Segment size is a number of nodes put into buffer for calculation with the streams.
  //const int segmentSize = 3 * _numOfBlocks;

  //const int leftOverNodes = _numOfNodes % segmentSize;
  //// Adding number of dummy nodes in order to keep 
  //const int numOfNodesExtended = _numOfNodes + (segmentSize - leftOverNodes);


  // creating the streams
  //cudaStreamCreateWithFlags(&_streamA,cudaStreamNonBlocking);
  //cudaStreamCreateWithFlags(&_streamB,cudaStreamNonBlocking);
  //cudaStreamCreateWithFlags(&_streamC,cudaStreamNonBlocking);

  //std::cout << "Streams created \n";

  allocateCuda();
  allocateHost();

  if(_storedTrajectories.size() != _numOfBlocks)
  {
    std::cout << "ERROR:: wrong size of trajectories vector; _numOfBlocks " << _numOfBlocks
              << " ; trajectories.size() " << _storedTrajectories.size() << "\n";
  }

  // copy the vector to pinned host array 
  std::copy(_storedTrajectories.begin(), _storedTrajectories.end(), _trajectoryHost);

  // copy from host to device
  cudaMemcpy(_trajectoryDeviceA,_trajectoryHost, sizeof(GridLine)*_numOfBlocks,cudaMemcpyHostToDevice);

  // initialize outputs
  cudaMemset(_ampsDeviceA,0,sizeof(float)*_numOfBlocks*_threadsPerBlock);

  // invoke kernel
  _invokeKernelTriple();

  // copy from device to host
  cudaMemcpy(_ampsHost,_ampsDeviceA, sizeof(float)*_numOfBlocks*_threadsPerBlock,cudaMemcpyDeviceToHost);

  //freeAll();

};

__device__
void rootsToAmps(double* amps, 
                 thrust::complex<double>* roots,
                 const thrust::complex<double> &sourcePos,
                 const thrust::complex<double> &z2,
                 const thrust::complex<double> &z3)
{
  for(unsigned int i = 0; i < 10; i++)
  {
    thrust::complex<double> testSourcePos = conj(roots[i])
                                            -conj(sourcePos)
                                            -paramsImg[3]/(roots[i])
                                            -paramsImg[4]/(roots[i]-z2)
                                            -paramsImg[5]/(roots[i]-z3);
  
    double detJac = thrust::norm(paramsImg[3]/(roots[i])/(roots[i]) +
                                 paramsImg[4]/(roots[i]-z2)/(roots[i]-z2)+ 
                                 paramsImg[5]/(roots[i]-z3)/(roots[i]-z3));

    amps[i] = double(thrust::norm(testSourcePos) < paramsImg[7])/sqrt(detJac); 
    //amps[i] = 1.0/sqrt(detJac); 
  }

  return;
}


// In order to improve readability of the .cc files the calculation of the coeffs is stored in a separate file.
__device__
void getCoeffs(thrust::complex<double>* coeffs, thrust::complex<double> zeta)
{
  thrust::complex<double> z2(paramsImg[3], paramsImg[4]);   
  thrust::complex<double> z2c(paramsImg[3], -paramsImg[4]);   
  thrust::complex<double> z3(paramsImg[5], paramsImg[6]);   
  thrust::complex<double> z3c(paramsImg[5], -paramsImg[6]);   
  double m2 = paramsImg[1]; 
  double m3 = paramsImg[2]; 
  thrust::complex<double> zetac= thrust::conj(zeta);

  coeffs[10] = (-zetac * z3c + zetac * zetac) * z2 + zetac * zetac * z3c - pow(zetac, 0.3e1);
  coeffs[9] = (-0.3e1 * zetac * zetac + 0.3e1 * zetac * z3c) * z2 * z2 + ((-0.3e1 * zetac * zetac + 0.3e1 * zetac * z3c) * z3 + (-m3 - 0.3e1 * zetac * zetac - m2 + zeta * zetac) * z3c + zetac * (0.3e1 * zetac * zetac - zeta * zetac + 0.1e1 + m2)) * z2 + (-0.3e1 * zetac * zetac * z3c + 0.3e1 * pow(zetac, 0.3e1)) * z3 + zetac * (-zeta * zetac + 0.1e1 + m3) * z3c + zetac * zetac * (zeta * zetac - 0.2e1);
  coeffs[8] = (-0.3e1 * zetac * z3c + 0.3e1 * zetac * zetac) * pow(z2, 0.3e1) + ((-0.9e1 * zetac * z3c + 0.9e1 * zetac * zetac) * z3 + (0.3e1 * zetac * zetac - 0.3e1 * zeta * zetac + 0.3e1 * m3 + 0.2e1 * m2) * z3c - zetac * (0.3e1 * zetac * zetac - 0.3e1 * zeta * zetac + 0.3e1 + m2)) * z2 * z2 + ((-0.3e1 * zetac * z3c + 0.3e1 * zetac * zetac) * z3 * z3 + ((0.9e1 * zetac * zetac + 0.3e1 * m2 - 0.3e1 * zeta * zetac + 0.2e1 * m3) * z3c - zetac * (0.9e1 * zetac * zetac - 0.3e1 * zeta * zetac + 0.3e1 - 0.2e1 * m3 + 0.3e1 * m2)) * z3 + (0.3e1 * zeta * zetac * zetac + 0.2e1 * zetac * m2 - 0.3e1 * zetac - 0.3e1 * zetac * m3 + zeta) * z3c + m2 + 0.6e1 * zetac * zetac - 0.3e1 * zetac * zetac * m2 - 0.3e1 * zeta * pow(zetac, 0.3e1) - 0.2e1 * zeta * zetac) * z2 + (0.3e1 * zetac * zetac * z3c - 0.3e1 * pow(zetac, 0.3e1)) * z3 * z3 + (-zetac * (-0.3e1 * zeta * zetac + 0.3e1 + m3) * z3c - 0.3e1 * zetac * zetac * (zeta * zetac - 0.2e1 + m3)) * z3 + (m3 - 0.2e1 * zeta * zetac) * z3c + zetac * (0.3e1 * zeta * zetac - 0.1e1);
  coeffs[7] = (-zetac * zetac + zetac * z3c) * pow(z2, 0.4e1) + ((-0.9e1 * zetac * zetac + 0.9e1 * zetac * z3c) * z3 + (0.3e1 * zeta * zetac - 0.3e1 * m3 - m2 - zetac * zetac) * z3c - zetac * (-zetac * zetac + 0.3e1 * zeta * zetac - 0.3e1 + m2)) * pow(z2, 0.3e1) + ((-0.9e1 * zetac * zetac + 0.9e1 * zetac * z3c) * z3 * z3 + ((0.9e1 * zeta * zetac - 0.6e1 * m3 - 0.9e1 * zetac * zetac - 0.6e1 * m2) * z3c + 0.3e1 * zetac * (0.3e1 * zetac * zetac - 0.3e1 * zeta * zetac + m2 - 0.2e1 * m3 + 0.3e1)) * z3 + (0.3e1 * zetac * m3 + zeta * m2 - 0.3e1 * zeta - 0.4e1 * zetac * m2 + 0.3e1 * zetac - 0.3e1 * zeta * zetac * zetac) * z3c + 0.6e1 * zeta * zetac - 0.2e1 * m2 - 0.2e1 * zeta * zetac * m2 + 0.6e1 * zetac * zetac * m2 + m2 * m2 - 0.6e1 * zetac * zetac + 0.3e1 * zeta * pow(zetac, 0.3e1)) * z2 * z2 + ((-zetac * zetac + zetac * z3c) * pow(z3, 0.3e1) + ((0.3e1 * zeta * zetac - 0.9e1 * zetac * zetac - m3 - 0.3e1 * m2) * z3c + zetac * (0.9e1 * zetac * zetac - 0.3e1 * zeta * zetac + 0.3e1 - 0.4e1 * m3 + 0.3e1 * m2)) * z3 * z3 + ((-0.9e1 * zeta * zetac * zetac - 0.3e1 * zeta + 0.3e1 * zetac * m3 + 0.9e1 * zetac - 0.6e1 * zetac * m2 + zeta * m3) * z3c + m3 + 0.9e1 * zetac * zetac * m3 - 0.3e1 * m2 - 0.18e2 * zetac * zetac + m3 * m2 - 0.2e1 * zeta * zetac * m3 + 0.9e1 * zeta * pow(zetac, 0.3e1) + 0.6e1 * zeta * zetac + 0.9e1 * zetac * zetac * m2) * z3 + (0.6e1 * zeta * zetac + m3 * m2 - 0.3e1 * m3 - 0.2e1 * zeta * zetac * m2 + m2) * z3c - 0.9e1 * zeta * zetac * zetac - 0.4e1 * zetac * m2 - zeta + 0.3e1 * zetac + 0.3e1 * zeta * zetac * zetac * m2) * z2 + (-zetac * zetac * z3c + pow(zetac, 0.3e1)) * pow(z3, 0.3e1) + (-zetac * (0.3e1 * zeta * zetac - 0.3e1 + m3) * z3c + 0.3e1 * zetac * zetac * (zeta * zetac - 0.2e1 + 0.2e1 * m3)) * z3 * z3 + ((0.6e1 * zeta * zetac - 0.2e1 * zeta * zetac * m3 + m3 * m3 - 0.2e1 * m3) * z3c + zetac * (0.3e1 * zeta * zetac * m3 - 0.9e1 * zeta * zetac - 0.4e1 * m3 + 0.3e1)) * z3 + 0.3e1 * zeta * zetac - zeta * z3c;
  coeffs[6] = ((-0.3e1 * zetac * z3c + 0.3e1 * zetac * zetac) * z3 + (-zeta * zetac + m3) * z3c + zetac * (zeta * zetac + m2 - 0.1e1)) * pow(z2, 0.4e1) + ((-0.9e1 * zetac * z3c + 0.9e1 * zetac * zetac) * z3 * z3 + ((0.6e1 * m3 - 0.9e1 * zeta * zetac + 0.3e1 * zetac * zetac + 0.3e1 * m2) * z3c + 0.3e1 * zetac * (-zetac * zetac + 0.3e1 * zeta * zetac + 0.2e1 * m3 - 0.3e1 + m2)) * z3 + (0.3e1 * zeta + 0.2e1 * zetac * m2 - zetac * m3 + zeta * zetac * zetac - 0.2e1 * zeta * m2 - zetac) * z3c - m2 * m2 - 0.6e1 * zeta * zetac - zeta * pow(zetac, 0.3e1) + 0.4e1 * zeta * zetac * m2 - 0.3e1 * zetac * zetac * m2 + m2 + 0.2e1 * zetac * zetac) * pow(z2, 0.3e1) + ((-0.3e1 * zetac * z3c + 0.3e1 * zetac * zetac) * pow(z3, 0.3e1) + ((0.3e1 * m3 - 0.9e1 * zeta * zetac + 0.6e1 * m2 + 0.9e1 * zetac * zetac) * z3c - 0.3e1 * zetac * (0.3e1 * zetac * zetac - 0.3e1 * zeta * zetac + m2 + 0.3e1 - 0.4e1 * m3)) * z3 * z3 + ((0.9e1 * zeta * zetac * zetac - 0.9e1 * zetac - 0.3e1 * zeta * m2 + 0.9e1 * zeta - 0.3e1 * zetac * m3 + 0.12e2 * zetac * m2 - 0.3e1 * zeta * m3) * z3c - 0.3e1 * m3 - 0.9e1 * zeta * pow(zetac, 0.3e1) - m3 * m2 + 0.6e1 * m2 + 0.18e2 * zetac * zetac + 0.6e1 * zeta * zetac * m2 + 0.6e1 * zeta * zetac * m3 - 0.9e1 * zetac * zetac * m3 - 0.3e1 * m2 * m2 - 0.18e2 * zetac * zetac * m2 - 0.18e2 * zeta * zetac) * z3 + (-0.6e1 * zeta * zetac - 0.2e1 * m2 + 0.3e1 * m3 - 0.2e1 * m3 * m2 + m2 * m2 + 0.4e1 * zeta * zetac * m2) * z3c + 0.8e1 * zetac * m2 - 0.3e1 * zetac + 0.3e1 * zeta + 0.9e1 * zeta * zetac * zetac - 0.6e1 * zeta * zetac * zetac * m2 - 0.2e1 * zeta * m2 - 0.3e1 * zetac * m2 * m2) * z2 * z2 + (((-zeta * zetac + 0.3e1 * zetac * zetac + m2) * z3c - zetac * (0.3e1 * zetac * zetac - zeta * zetac + 0.1e1 - 0.2e1 * m3 + m2)) * pow(z3, 0.3e1) + ((0.9e1 * zeta * zetac * zetac + 0.3e1 * zetac * m3 + 0.3e1 * zeta + 0.6e1 * zetac * m2 - 0.9e1 * zetac - 0.2e1 * zeta * m3) * z3c - 0.9e1 * zetac * zetac * m2 + 0.18e2 * zetac * zetac - 0.18e2 * zetac * zetac * m3 - 0.2e1 * m3 - 0.6e1 * zeta * zetac + m3 * m3 + 0.3e1 * m2 + 0.4e1 * zeta * zetac * m3 - 0.2e1 * m3 * m2 - 0.9e1 * zeta * pow(zetac, 0.3e1)) * z3 * z3 + ((0.6e1 * zeta * zetac * m2 + 0.6e1 * m3 - 0.18e2 * zeta * zetac - 0.3e1 * m2 - 0.3e1 * m3 * m3 + 0.6e1 * zeta * zetac * m3 - m3 * m2) * z3c + 0.27e2 * zeta * zetac * zetac - 0.2e1 * zeta * m3 - 0.9e1 * zeta * zetac * zetac * m3 - 0.6e1 * zetac * m3 * m2 - 0.9e1 * zeta * zetac * zetac * m2 + 0.3e1 * zeta + 0.12e2 * zetac * m2 + 0.12e2 * zetac * m3 - 0.9e1 * zetac) * z3 - zeta * (-0.3e1 + 0.2e1 * m2) * z3c + 0.6e1 * zeta * zetac * m2 - m2 - 0.9e1 * zeta * zetac) * z2 + (zetac * (zeta * zetac - 0.1e1 + m3) * z3c - zetac * zetac * (zeta * zetac + 0.3e1 * m3 - 0.2e1)) * pow(z3, 0.3e1) + ((-0.6e1 * zeta * zetac + 0.4e1 * zeta * zetac * m3 + m3 - m3 * m3) * z3c - zetac * (-0.9e1 * zeta * zetac + 0.6e1 * zeta * zetac * m3 + 0.3e1 * m3 * m3 + 0.3e1 - 0.8e1 * m3)) * z3 * z3 + (-zeta * (-0.3e1 + 0.2e1 * m3) * z3c - 0.9e1 * zeta * zetac + 0.6e1 * zeta * zetac * m3 - m3) * z3 + zeta;
  coeffs[5] = ((-0.3e1 * zetac * zetac + 0.3e1 * zetac * z3c) * z3 * z3 + ((0.3e1 * zeta * zetac - 0.2e1 * m3) * z3c - zetac * (0.3e1 * zeta * zetac + 0.2e1 * m3 + 0.3e1 * m2 - 0.3e1)) * z3 + zeta * (m2 - 0.1e1) * z3c - 0.2e1 * zeta * zetac * (m2 - 0.1e1)) * pow(z2, 0.4e1) + ((-0.3e1 * zetac * zetac + 0.3e1 * zetac * z3c) * pow(z3, 0.3e1) + ((0.9e1 * zeta * zetac - 0.3e1 * zetac * zetac - 0.3e1 * m3 - 0.3e1 * m2) * z3c - 0.3e1 * zetac * (-zetac * zetac + 0.3e1 * zeta * zetac - 0.3e1 + 0.4e1 * m3 + m2)) * z3 * z3 + ((-0.3e1 * zeta * zetac * zetac + 0.3e1 * zetac + 0.6e1 * zeta * m2 + zetac * m3 - 0.6e1 * zetac * m2 + 0.3e1 * zeta * m3 - 0.9e1 * zeta) * z3c + 0.3e1 * zeta * pow(zetac, 0.3e1) - 0.12e2 * zeta * zetac * m2 - 0.3e1 * m2 - 0.6e1 * zetac * zetac - m3 * m2 + 0.18e2 * zeta * zetac + 0.3e1 * zetac * zetac * m3 + 0.3e1 * m2 * m2 + 0.3e1 * m3 + 0.9e1 * zetac * zetac * m2 - 0.6e1 * zeta * zetac * m3) * z3 - (m2 - 0.1e1) * (m2 - m3 + 0.2e1 * zeta * zetac) * z3c - (m2 - 0.1e1) * (zeta * m2 - 0.3e1 * zetac * m2 - 0.3e1 * zeta * zetac * zetac + zetac - 0.3e1 * zeta)) * pow(z2, 0.3e1) + (((0.3e1 * zeta * zetac - 0.2e1 * m2 - 0.3e1 * zetac * zetac) * z3c + zetac * (0.3e1 * zetac * zetac - 0.3e1 * zeta * zetac + 0.3e1 - 0.6e1 * m3 + m2)) * pow(z3, 0.3e1) + ((-0.12e2 * zetac * m2 - 0.9e1 * zeta + 0.9e1 * zetac - 0.3e1 * zetac * m3 + 0.3e1 * zeta * m2 - 0.9e1 * zeta * zetac * zetac + 0.6e1 * zeta * m3) * z3c + 0.6e1 * m3 + 0.9e1 * zeta * pow(zetac, 0.3e1) - 0.6e1 * zeta * zetac * m2 - 0.6e1 * m2 + 0.18e2 * zeta * zetac + 0.2e1 * m3 * m2 - 0.3e1 * m3 * m3 + 0.18e2 * zetac * zetac * m3 - 0.12e2 * zeta * zetac * m3 - 0.18e2 * zetac * zetac + 0.18e2 * zetac * zetac * m2 + 0.3e1 * m2 * m2) * z3 * z3 + ((0.6e1 * m2 - 0.3e1 * m2 * m2 + 0.18e2 * zeta * zetac + 0.3e1 * m3 * m3 - 0.6e1 * zeta * zetac * m3 - 0.12e2 * zeta * zetac * m2 - 0.6e1 * m3 + 0.2e1 * m3 * m2) * z3c + 0.6e1 * zeta * m2 + 0.9e1 * zetac + 0.12e2 * zetac * m3 * m2 - 0.27e2 * zeta * zetac * zetac - 0.24e2 * zetac * m2 - 0.2e1 * zeta * m2 * m3 - 0.12e2 * zetac * m3 + 0.6e1 * zeta * m3 - 0.9e1 * zeta + 0.9e1 * zeta * zetac * zetac * m3 + 0.18e2 * zeta * zetac * zetac * m2 + 0.9e1 * zetac * m2 * m2) * z3 - zeta * (m2 - 0.1e1) * (m2 - 0.3e1) * z3c + (m2 - 0.1e1) * (0.3e1 * zeta * zetac * m2 - 0.2e1 * m2 - 0.9e1 * zeta * zetac)) * z2 * z2 + (((-zeta - 0.2e1 * zetac * m2 + 0.3e1 * zetac - 0.3e1 * zetac * m3 - 0.3e1 * zeta * zetac * zetac + zeta * m3) * z3c + 0.9e1 * zetac * zetac * m3 - m2 + 0.3e1 * zeta * pow(zetac, 0.3e1) + m3 + 0.2e1 * zeta * zetac - 0.2e1 * zeta * zetac * m3 + m3 * m2 + 0.3e1 * zetac * zetac * m2 - m3 * m3 - 0.6e1 * zetac * zetac) * pow(z3, 0.3e1) + ((0.3e1 * m3 * m3 + 0.3e1 * m2 + 0.18e2 * zeta * zetac - 0.6e1 * zeta * zetac * m2 - 0.3e1 * m3 - 0.12e2 * zeta * zetac * m3 - m3 * m2) * z3c + 0.9e1 * zetac * m3 * m3 - 0.27e2 * zeta * zetac * zetac + 0.18e2 * zeta * zetac * zetac * m3 + 0.12e2 * zetac * m3 * m2 - zeta * m3 * m3 - 0.3e1 * zeta - 0.12e2 * zetac * m2 + 0.9e1 * zeta * zetac * zetac * m2 + 0.9e1 * zetac - 0.24e2 * zetac * m3 + 0.4e1 * zeta * m3) * z3 * z3 + (-zeta * (0.2e1 * m3 * m2 + 0.9e1 - 0.6e1 * m3 - 0.6e1 * m2) * z3c - 0.18e2 * zeta * zetac * m2 + 0.3e1 * m2 + 0.3e1 * m3 + 0.6e1 * zeta * zetac * m2 * m3 - 0.18e2 * zeta * zetac * m3 - 0.4e1 * m3 * m2 + 0.27e2 * zeta * zetac) * z3 + 0.3e1 * zeta * (m2 - 0.1e1)) * z2 + (-0.2e1 * zeta * zetac * (-0.1e1 + m3) * z3c + (-0.1e1 + m3) * zetac * (0.3e1 * m3 + 0.3e1 * zeta * zetac - 0.1e1)) * pow(z3, 0.3e1) + (-zeta * (-0.1e1 + m3) * (-0.3e1 + m3) * z3c + (-0.1e1 + m3) * (0.3e1 * zeta * zetac * m3 - 0.2e1 * m3 - 0.9e1 * zeta * zetac)) * z3 * z3 + 0.3e1 * zeta * (-0.1e1 + m3) * z3;
  coeffs[4] = ((-zetac * z3c + zetac * zetac) * pow(z3, 0.3e1) + ((-0.3e1 * zeta * zetac + m3) * z3c + zetac * (0.3e1 * zeta * zetac - 0.3e1 + 0.3e1 * m2 + 0.4e1 * m3)) * z3 * z3 + (-zeta * (m3 - 0.3e1 + 0.3e1 * m2) * z3c + 0.6e1 * zeta * zetac * m2 + 0.2e1 * zeta * zetac * m3 + m3 * m2 - m3 - 0.6e1 * zeta * zetac) * z3 + zeta * pow(m2 - 0.1e1, 0.2e1)) * pow(z2, 0.4e1) + (((-0.3e1 * zeta * zetac + zetac * zetac + m2) * z3c + zetac * (-zetac * zetac + 0.3e1 * zeta * zetac - 0.3e1 + 0.6e1 * m3 + m2)) * pow(z3, 0.3e1) + ((0.9e1 * zeta - 0.6e1 * zeta * m2 + zetac * m3 + 0.3e1 * zeta * zetac * zetac + 0.6e1 * zetac * m2 - 0.6e1 * zeta * m3 - 0.3e1 * zetac) * z3c - 0.6e1 * m3 + 0.12e2 * zeta * zetac * m2 + 0.3e1 * m2 + 0.2e1 * m3 * m2 - 0.3e1 * zeta * pow(zetac, 0.3e1) + 0.3e1 * m3 * m3 - 0.6e1 * zetac * zetac * m3 - 0.18e2 * zeta * zetac + 0.12e2 * zeta * zetac * m3 + 0.6e1 * zetac * zetac - 0.9e1 * zetac * zetac * m2 - 0.3e1 * m2 * m2) * z3 * z3 + ((0.3e1 * m2 * m2 - 0.6e1 * zeta * zetac - 0.3e1 * m2 + 0.2e1 * zeta * zetac * m3 - m3 * m3 + 0.2e1 * m3 - m3 * m2 + 0.6e1 * zeta * zetac * m2) * z3c - 0.12e2 * zeta * m2 - 0.3e1 * zetac - 0.6e1 * zetac * m3 * m2 + 0.9e1 * zeta * zetac * zetac + 0.9e1 * zeta + 0.3e1 * zeta * m2 * m2 + 0.4e1 * zetac * m3 - 0.6e1 * zeta * m3 + 0.12e2 * zetac * m2 + 0.4e1 * zeta * m2 * m3 - 0.3e1 * zeta * zetac * zetac * m3 - 0.9e1 * zeta * zetac * zetac * m2 - 0.9e1 * zetac * m2 * m2) * z3 + zeta * pow(m2 - 0.1e1, 0.2e1) * z3c - pow(m2 - 0.1e1, 0.2e1) * (m2 + 0.3e1 * zeta * zetac)) * pow(z2, 0.3e1) + (((-0.3e1 * zetac - 0.3e1 * zeta * m3 - zeta * m2 + 0.4e1 * zetac * m2 + 0.3e1 * zetac * m3 + 0.3e1 * zeta * zetac * zetac + 0.3e1 * zeta) * z3c - 0.3e1 * m3 + 0.3e1 * m3 * m3 - 0.3e1 * zeta * pow(zetac, 0.3e1) - 0.6e1 * zeta * zetac + 0.2e1 * m2 + 0.2e1 * zeta * zetac * m2 - 0.6e1 * zetac * zetac * m2 + 0.6e1 * zetac * zetac - 0.9e1 * zetac * zetac * m3 - m2 * m2 - m3 * m2 + 0.6e1 * zeta * zetac * m3) * pow(z3, 0.3e1) + ((0.12e2 * zeta * zetac * m3 + 0.12e2 * zeta * zetac * m2 - 0.18e2 * zeta * zetac + 0.3e1 * m3 + 0.3e1 * m2 * m2 + 0.2e1 * m3 * m2 - 0.3e1 * m3 * m3 - 0.6e1 * m2) * z3c + 0.9e1 * zeta - 0.24e2 * zetac * m3 * m2 - 0.6e1 * zeta * m2 + 0.24e2 * zetac * m3 - 0.9e1 * zetac * m3 * m3 + 0.27e2 * zeta * zetac * zetac + 0.3e1 * zeta * m3 * m3 + 0.4e1 * zeta * m2 * m3 - 0.12e2 * zeta * m3 - 0.9e1 * zetac - 0.18e2 * zeta * zetac * zetac * m3 + 0.24e2 * zetac * m2 - 0.9e1 * zetac * m2 * m2 - 0.18e2 * zeta * zetac * zetac * m2) * z3 * z3 + (zeta * (-0.6e1 * m3 - 0.12e2 * m2 + 0.4e1 * m3 * m2 + 0.9e1 + 0.3e1 * m2 * m2) * z3c - 0.3e1 * m3 + 0.18e2 * zeta * zetac * m3 - 0.12e2 * zeta * zetac * m2 * m3 + 0.6e1 * m2 * m2 - 0.3e1 * m2 * m2 * m3 + 0.36e2 * zeta * zetac * m2 - 0.27e2 * zeta * zetac + 0.8e1 * m3 * m2 - 0.9e1 * zeta * zetac * m2 * m2 - 0.6e1 * m2) * z3 + 0.3e1 * zeta * pow(m2 - 0.1e1, 0.2e1)) * z2 * z2 + (((m3 * m2 - 0.6e1 * zeta * zetac + 0.6e1 * zeta * zetac * m3 + 0.2e1 * zeta * zetac * m2 - m2) * z3c - 0.3e1 * zeta * zetac * zetac * m2 + 0.4e1 * zetac * m2 - 0.9e1 * zetac * m3 * m3 - 0.9e1 * zeta * zetac * zetac * m3 - 0.6e1 * zetac * m3 * m2 + zeta * m3 * m3 + 0.9e1 * zeta * zetac * zetac - 0.3e1 * zetac + 0.12e2 * zetac * m3 + zeta - 0.2e1 * zeta * m3) * pow(z3, 0.3e1) + (zeta * (-0.12e2 * m3 + 0.3e1 * m3 * m3 + 0.9e1 + 0.4e1 * m3 * m2 - 0.6e1 * m2) * z3c - 0.9e1 * zeta * zetac * m3 * m3 + 0.6e1 * m3 * m3 - 0.27e2 * zeta * zetac - 0.12e2 * zeta * zetac * m2 * m3 - 0.3e1 * m2 + 0.18e2 * zeta * zetac * m2 + 0.36e2 * zeta * zetac * m3 - 0.6e1 * m3 - 0.3e1 * m3 * m3 * m2 + 0.8e1 * m3 * m2) * z3 * z3 + 0.3e1 * zeta * (0.3e1 - 0.3e1 * m2 + 0.2e1 * m3 * m2 - 0.3e1 * m3) * z3) * z2 + (zeta * pow(-0.1e1 + m3, 0.2e1) * z3c - pow(-0.1e1 + m3, 0.2e1) * (m3 + 0.3e1 * zeta * zetac)) * pow(z3, 0.3e1) + 0.3e1 * zeta * pow(-0.1e1 + m3, 0.2e1) * z3 * z3;
  coeffs[3] = ((zeta * zetac * z3c - zetac * (zeta * zetac + m2 + 0.2e1 * m3 - 0.1e1)) * pow(z3, 0.3e1) + (zeta * (0.2e1 * m3 - 0.3e1 + 0.3e1 * m2) * z3c - m3 * m3 + 0.6e1 * zeta * zetac - 0.4e1 * zeta * zetac * m3 - 0.6e1 * zeta * zetac * m2 + 0.2e1 * m3 - 0.2e1 * m3 * m2) * z3 * z3 - zeta * (m2 - 0.1e1) * (0.2e1 * m3 - 0.3e1 + 0.3e1 * m2) * z3) * pow(z2, 0.4e1) + (((0.2e1 * zeta * m2 + zetac - zeta * zetac * zetac - 0.2e1 * zetac * m2 - zetac * m3 - 0.3e1 * zeta + 0.3e1 * zeta * m3) * z3c - m3 * m2 + 0.3e1 * m3 - 0.4e1 * zeta * zetac * m2 - m2 - 0.3e1 * m3 * m3 + 0.6e1 * zeta * zetac - 0.6e1 * zeta * zetac * m3 - 0.2e1 * zetac * zetac + zeta * pow(zetac, 0.3e1) + 0.3e1 * zetac * zetac * m2 + m2 * m2 + 0.3e1 * zetac * zetac * m3) * pow(z3, 0.3e1) + ((0.3e1 * m2 + 0.6e1 * zeta * zetac - 0.6e1 * zeta * zetac * m2 - 0.3e1 * m2 * m2 - m3 + m3 * m3 - 0.4e1 * zeta * zetac * m3 - m3 * m2) * z3c + 0.3e1 * zetac * m3 * m3 - 0.8e1 * zetac * m3 - 0.8e1 * zeta * m2 * m3 - 0.9e1 * zeta + 0.6e1 * zeta * zetac * zetac * m3 + 0.9e1 * zetac * m2 * m2 + 0.12e2 * zeta * m2 + 0.3e1 * zetac - 0.12e2 * zetac * m2 - 0.9e1 * zeta * zetac * zetac - 0.3e1 * zeta * m2 * m2 - 0.3e1 * zeta * m3 * m3 + 0.12e2 * zetac * m3 * m2 + 0.12e2 * zeta * m3 + 0.9e1 * zeta * zetac * zetac * m2) * z3 * z3 + (-zeta * (m2 - 0.1e1) * (0.2e1 * m3 - 0.3e1 + 0.3e1 * m2) * z3c + (m2 - 0.1e1) * (0.3e1 * m2 * m2 - 0.3e1 * m2 + 0.9e1 * zeta * zetac * m2 + 0.3e1 * m3 * m2 - 0.9e1 * zeta * zetac + 0.6e1 * zeta * zetac * m3 - m3)) * z3 + zeta * pow(m2 - 0.1e1, 0.3e1)) * pow(z2, 0.3e1) + (((-m2 * m2 - 0.6e1 * zeta * zetac * m3 + 0.6e1 * zeta * zetac - 0.4e1 * zeta * zetac * m2 - 0.2e1 * m3 * m2 + 0.2e1 * m2) * z3c - 0.2e1 * zeta * m2 * m3 + 0.6e1 * zeta * zetac * zetac * m2 + 0.3e1 * zetac - 0.3e1 * zeta - 0.12e2 * zetac * m3 + 0.9e1 * zetac * m3 * m3 + 0.3e1 * zetac * m2 * m2 + 0.6e1 * zeta * m3 + 0.2e1 * zeta * m2 - 0.9e1 * zeta * zetac * zetac + 0.9e1 * zeta * zetac * zetac * m3 + 0.12e2 * zetac * m3 * m2 - 0.8e1 * zetac * m2 - 0.3e1 * zeta * m3 * m3) * pow(z3, 0.3e1) + (-zeta * (0.8e1 * m3 * m2 + 0.3e1 * m3 * m3 + 0.9e1 - 0.12e2 * m2 - 0.12e2 * m3 + 0.3e1 * m2 * m2) * z3c - 0.36e2 * zeta * zetac * m3 + 0.24e2 * zeta * zetac * m2 * m3 + 0.6e1 * m2 + 0.6e1 * m2 * m2 * m3 - 0.36e2 * zeta * zetac * m2 - 0.6e1 * m3 * m3 - 0.6e1 * m2 * m2 - 0.16e2 * m3 * m2 + 0.6e1 * m3 + 0.6e1 * m3 * m3 * m2 + 0.9e1 * zeta * zetac * m3 * m3 + 0.9e1 * zeta * zetac * m2 * m2 + 0.27e2 * zeta * zetac) * z3 * z3 + 0.3e1 * zeta * (m2 - 0.1e1) * (-0.3e1 * m2 + m3 * m2 - 0.3e1 * m3 + 0.3e1) * z3) * z2 * z2 + ((-zeta * (-0.1e1 + m3) * (0.3e1 * m3 - 0.3e1 + 0.2e1 * m2) * z3c + (-0.1e1 + m3) * (0.3e1 * m3 * m3 - 0.3e1 * m3 + 0.3e1 * m3 * m2 + 0.9e1 * zeta * zetac * m3 + 0.6e1 * zeta * zetac * m2 - m2 - 0.9e1 * zeta * zetac)) * pow(z3, 0.3e1) + 0.3e1 * zeta * (-0.1e1 + m3) * (-0.3e1 * m2 + m3 * m2 - 0.3e1 * m3 + 0.3e1) * z3 * z3) * z2 + zeta * pow(-0.1e1 + m3, 0.3e1) * pow(z3, 0.3e1);
  coeffs[2] = ((-(-0.1e1 + m2 + m3) * zeta * z3c + (-0.1e1 + m2 + m3) * (0.2e1 * zeta * zetac + m3)) * pow(z3, 0.3e1) + (-0.1e1 + m2 + m3) * zeta * (m3 - 0.3e1 + 0.3e1 * m2) * z3 * z3) * pow(z2, 0.4e1) + (((-0.1e1 + m2 + m3) * (0.2e1 * zeta * zetac + m2) * z3c + (-0.1e1 + m2 + m3) * (0.3e1 * zeta * m3 - 0.3e1 * zeta + zetac - 0.3e1 * zetac * m3 - 0.3e1 * zeta * zetac * zetac + zeta * m2 - 0.3e1 * zetac * m2)) * pow(z3, 0.3e1) + ((-0.1e1 + m2 + m3) * zeta * (m3 - 0.3e1 + 0.3e1 * m2) * z3c - (-0.1e1 + m2 + m3) * (0.9e1 * zeta * zetac * m2 + 0.3e1 * m3 * m2 + 0.3e1 * zeta * zetac * m3 + 0.3e1 * m2 * m2 - 0.9e1 * zeta * zetac - 0.2e1 * m3 - 0.3e1 * m2)) * z3 * z3 - 0.3e1 * (-0.1e1 + m2 + m3) * zeta * pow(m2 - 0.1e1, 0.2e1) * z3) * pow(z2, 0.3e1) + (((-0.1e1 + m2 + m3) * zeta * (-0.3e1 + 0.3e1 * m3 + m2) * z3c - (-0.1e1 + m2 + m3) * (0.9e1 * zeta * zetac * m3 + 0.3e1 * m3 * m3 - 0.9e1 * zeta * zetac + 0.3e1 * zeta * zetac * m2 - 0.3e1 * m3 + 0.3e1 * m3 * m2 - 0.2e1 * m2)) * pow(z3, 0.3e1) - 0.3e1 * (-0.1e1 + m2 + m3) * zeta * (0.3e1 - 0.3e1 * m2 + 0.2e1 * m3 * m2 - 0.3e1 * m3) * z3 * z3) * z2 * z2 - 0.3e1 * (-0.1e1 + m2 + m3) * zeta * pow(-0.1e1 + m3, 0.2e1) * pow(z3, 0.3e1) * z2;
  coeffs[1] = -pow(z3, 0.3e1) * pow(-0.1e1 + m2 + m3, 0.2e1) * zeta * pow(z2, 0.4e1) + ((-pow(-0.1e1 + m2 + m3, 0.2e1) * zeta * z3c + pow(-0.1e1 + m2 + m3, 0.2e1) * (m2 + m3 + 0.3e1 * zeta * zetac)) * pow(z3, 0.3e1) + 0.3e1 * pow(-0.1e1 + m2 + m3, 0.2e1) * zeta * (m2 - 0.1e1) * z3 * z3) * pow(z2, 0.3e1) + 0.3e1 * pow(-0.1e1 + m2 + m3, 0.2e1) * zeta * (-0.1e1 + m3) * pow(z3, 0.3e1) * z2 * z2;
  coeffs[0] = -zeta * pow(z2, 0.3e1) * pow(z3, 0.3e1) * pow(-0.1e1 + m2 + m3, 0.3e1);
};

__device__
void getCoeffsBinOpt(thrust::complex<double>* coeffs, thrust::complex<double> zeta)
{
  // m2 paramsImg[1]
  // d  paramsImg[3]

  thrust::complex<double> z2(paramsImg[3], paramsImg[4]);   
  thrust::complex<double> z2c(paramsImg[3], -paramsImg[4]);   
  thrust::complex<double> z3(paramsImg[5], paramsImg[6]);   
  thrust::complex<double> z3c(paramsImg[5], -paramsImg[6]);
  thrust::complex<double> zetac = conj(zeta);
  double m2 = paramsImg[1]; 
  double m3 = paramsImg[2]; 

  thrust::complex<double> t1 = zetac*zetac;
  thrust::complex<double> a5 = -t1+zetac*paramsImg[3];
  thrust::complex<double> t6 = 2.0*t1*paramsImg[3];
  thrust::complex<double> t7 = paramsImg[3]*paramsImg[3];
  thrust::complex<double> t11 = paramsImg[1]*paramsImg[3];
  thrust::complex<double> a4 = -a5*zeta+t6+(-1.0-2.0*t7)*zetac+paramsImg[3]-t11;
  thrust::complex<double> t17 = t1*t7;
  thrust::complex<double> t18 = 2.0*t11;
  thrust::complex<double> t19 = t7*paramsImg[3];
  thrust::complex<double> t22 = paramsImg[1]*t7;
  thrust::complex<double> a3 = (-t6+(2.0*t7+2.0)*zetac-paramsImg[3])*zeta-t17+(t18+t19)*zetac-t7+t22;
  thrust::complex<double> a2 = (t17+(-2.0*paramsImg[3]-t19-t18)*zetac+t22+1.0+t7)*zeta+(-2.0*t22+t7)*zetac-paramsImg[3]+t11;
  thrust::complex<double> a1 = (2.0*t7*zetac*paramsImg[1]+t11*(-2.0-t7))*zeta+t11*(paramsImg[3]-t11);
  thrust::complex<double> a0 = zeta*paramsImg[1]*paramsImg[1]*t7;

  coeffs[0] = a0;
  coeffs[1] = a1;
  coeffs[2] = a2;
  coeffs[3] = a3;
  coeffs[4] = a4;
  coeffs[5] = a5;
}


