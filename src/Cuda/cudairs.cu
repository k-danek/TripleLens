#include <fstream>
#include <iomanip>
#include <vector>
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
float getAmpKernel(const std::vector<float>& collectedPoints,
                   double                    a,
                   double                    b,
                   double                    th,
                   double                    m2,
                   double                    m3,
                   double                    sourceSize,
                   double                    sourcePosX,
                   double                    sourcePosY,
                   double                    imgPixSize)
{

  int size = collectedPoints.size(); 
  int numOfPoints;
  if( size % 2 != 0)
  {
    std::cout << "getAmpKernel: Odd size of points in collectedPoints.\n";
    numOfPoints = 0;
    abort();
  }
  else
  {
    numOfPoints = size / 2;
  }
  
  float* amps;
  float* params;
  float* collectedPointsShared;

  // Allocate Unified Memory – accessible from CPU or GPU
  cudaMallocManaged(&collectedPointsShared, size*sizeof(float));
  cudaMallocManaged(&params, 10*sizeof(float));
  cudaMallocManaged(&amps, numOfPoints*sizeof(float));

  params[0]  = float(a);          // a
  params[1]  = float(b);          // b
  params[2]  = float(th);         // th
  params[3]  = float(1.0-m2-m3);  // m1
  params[4]  = float(m2);         // m2
  params[5]  = float(m3);         // m3
  params[6]  = float(sourceSize); // sourceSize
  params[7]  = float(sourcePosX); // source position X
  params[8]  = float(sourcePosY); // source position Y
  params[9]  = float(imgPixSize); // size of pixel

  // initialize x and y arrays on the host
  for (int i = 0; i < size; i++)
  {
    collectedPointsShared[i] = collectedPoints[i];
  }


  // Run kernel on 1M elements on the GPU
  int threadsPerBlock = 1<<5;
  int numBlocks = (numOfPoints + threadsPerBlock - 1) / threadsPerBlock;

  cudaProfilerStart(); 

  arrangeShooting<<<numBlocks, threadsPerBlock>>>(collectedPointsShared, params, amps);

  cudaProfilerStop();

  cudaDeviceSynchronize();

  double totalAmp = 0.0;
  for(int i = 0; i < numOfPoints; i++)
  {
    totalAmp += amps[i];
  }

  //// Some debugging stuff
  //const std::complex<float> z2 = std::complex<float>(params[0],0.0);
  //const std::complex<float> z3 = std::complex<float>(params[1]*cos(params[2]),
  //                                                         params[1]*sin(params[2]));
 
  //// increment in image plane iteration
  ////const float inc = params[9]/8.0;

  //std::complex<float> sourcePos = std::complex<float>(params[7], params[8]);
  //std::complex<float> imgPos = std::complex<float>(collectedPoints[0],
  //                                                 collectedPoints[1]);

  //float irsAmp = irsCPU(params,z2,z3,imgPos,sourcePos);

  //std::cout << "Testing one IRS: " << irsAmp << "\n";
  //std::cout << "Used points were: " << collectedPoints[0] << " " << collectedPoints[1] << "\n";
  //std::cout << "Source pos: " << sourcePos.real() << " " << sourcePos.imag() << "\n";
  //std::cout << "Source size and img plane par: " << params[6] << " " << params[9] << "\n";

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
  const int subgridSize = 32;
  
  //const thrust::complex<float> z1 = thrust::complex<float>(0.0,0.0);
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
  
  for (int i = 0; i < subgridSize; i++)
  { 
    imgPos.real(imgPos.real() + inc);
    //float imgY = __int2float_rn(i) * inc;    
    for (int j = 0; j < subgridSize; j++)
    {
      imgPos.imag(imgPos.imag() + inc);
      
      amps[index] += irs(params, z2, z3, imgPos, sourcePos); 

      //amps[index] += 1.0; 

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

    //return 1.0/(1.0+r*r);
};

float irsCPU(const float*                  params,
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

  //  return 1.0/(1.0+r*r);
};