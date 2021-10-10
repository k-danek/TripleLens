#include <fstream>
#include <iomanip>
#include <cuda_profiler_api.h>

#include "cudairs.cuh"

// Check whether source position is hit 
// In context of the point source, that source radius is just an error term.
//z1 = {0.0,0.0};
//z2 = {a,0.0};
//z3 = {b*cos(th), b*sin(th)};
//a =  params[0];
//b =  params[1];
//th = params[2];
//m1 = params[3];
//m2 = params[4];
//m3 = params[5];
//sourceSize = params[6] 
//sourcePosX = params[7] 
//sourcePosY = params[8] 
//imgPixSize = params[9]
// Kernel
float getAmpKernel(float* collectedPoints,
                   double a,
                   double b,
                   double th,
                   double m2,
                   double m3,
                   double sourceSize,
                   double sourcePosX,
                   double sourcePosY,
                   double imgPixSize)
{

  int size = 384; 
  float* amps;
  float* params;
  float* collectedPointsShared;

  // Allocate Unified Memory – accessible from CPU or GPU
  cudaMallocManaged(&collectedPointsShared, 2*size*sizeof(float));
  cudaMallocManaged(&params, 10*sizeof(float));
  cudaMallocManaged(&amps, size*sizeof(float));

  params[0]  = float(a);          // a
  params[1]  = float(b);          // b
  params[2]  = float(th);         // th
  params[3]  = float(1.0-m2-m3);  // m1
  params[4]  = float(m2);         // m2
  params[5]  = float(m3);         // m3
  params[6]  = float(sourceSize); // img origin X
  params[7]  = float(sourcePosX); // img origin Y
  params[8]  = float(sourcePosY); // img range
  params[9]  = float(imgPixSize); // size of pixel

  // initialize x and y arrays on the host
  for (int i = 0; i < 2*size; i++)
  {
    collectedPointsShared[i] = collectedPoints[i];
  }


  // Run kernel on 1M elements on the GPU
  int threadsPerBlock = 1<<5;
  int numBlocks = (size + threadsPerBlock - 1) / threadsPerBlock;

  cudaProfilerStart(); 

  arrangeShooting<<<numBlocks, threadsPerBlock>>>(collectedPointsShared, params, amps);

  cudaProfilerStop();

  cudaDeviceSynchronize();

  double totalAmp = 0.0;
  for(int i = 0; i < size; i++)
  {
    totalAmp += amps[i];
  }

  // Free memory
  cudaFree(amps);
  cudaFree(params);
  cudaFree(collectedPointsShared);

  return totalAmp;
};


__global__
void arrangeShooting(float* collectedPoints,
                     float* params,
                     float* amps)
{
  // Each pixel will be subdivided into finer grid.
  // subgridSize determines how fine the subgrid should be.
  const int subgridSize = 8;
  
  const thrust::complex<float> z1 = thrust::complex<float>(0.0,0.0);
  const thrust::complex<float> z2 = thrust::complex<float>(params[0],0.0);
  const thrust::complex<float> z3 = thrust::complex<float>(params[1]*cos(params[2]),
                                                           params[1]*sin(params[2]));
 
  // increment in image plane iteration
  const float inc = params[9]/__int2float_rn(subgridSize);

  //const int maxMapInd = subgridSize*subgridSize - 1;

  // actual index of a thread
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  //int stride = blockDim.x * gridDim.x;

  thrust::complex<float> sourcePos = thrust::complex<float>(params[7], params[8]);
  thrust::complex<float> imgPos = thrust::complex<float>(
                                                     collectedPoints[2*index]-params[9]/2,
                                                     collectedPoints[2*index+1]-params[9]/2);
  
  for (int i = index; i < subgridSize; i++)
  { 
    imgPos.real(imgPos.real() + inc);
    //float imgY = __int2float_rn(i) * inc;    
    for (int j = 0; j < subgridSize; j++)
    {
      imgPos.imag(imgPos.imag() + inc);
      
      amps[index] = irs(params, z2, z3, imgPos, sourcePos); 

    }
  }

};


__device__
float irs(const float*                  params,
          //const thrust::complex<float>& z1,
          const thrust::complex<float>& z2,
          const thrust::complex<float>& z3,
          const thrust::complex<float>& img,
          const thrust::complex<float>& sourcePos)
{

    thrust::complex<float> impact = img-params[3]/conj(img)
                                    -params[4]/conj(img-z2)
                                    -params[5]/conj(img-z3);

    float r = thrust::abs(impact-sourcePos)/params[6];

    if(r <= 1.0)
    {
        return 0.6+0.4*sqrt(1.0-r*r); 
    }
    else
    {
        return 0.0;
    }
};

