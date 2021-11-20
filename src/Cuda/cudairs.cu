#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <cuda_profiler_api.h>

#include "cudairs.cuh"


__constant__ float params[15];


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
// This actually triggers a kernel run
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

  // Each pixel will be subdivided into finer grid.
  // subgridSize determines how fine the subgrid should be.
  const int subgridSize = 16;
  
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
  //float* collectedPointsShared;
  float* collectedPointsDevice;

  cudaMalloc( (void**)&collectedPointsDevice, size*sizeof(float));

  // Allocate Unified Memory – accessible from CPU or GPU
  //cudaMallocManaged(&collectedPointsShared, size*sizeof(float));
  //cudaMallocManaged(&params, 10*sizeof(float));
  cudaMallocManaged(&amps, numOfPoints*sizeof(float));

  float* tempParams = (float*)malloc(sizeof(float)*15);

  tempParams[0]  = float(a);          // a
  tempParams[1]  = float(b);          // b
  tempParams[2]  = float(th);         // th
  tempParams[3]  = float(1.0-m2-m3);  // m1
  tempParams[4]  = float(m2);         // m2
  tempParams[5]  = float(m3);         // m3
  tempParams[6]  = float(sourceSize); // sourceSize
  tempParams[7]  = float(sourcePosX); // source position X
  tempParams[8]  = float(sourcePosY); // source position Y
  tempParams[9]  = float(imgPixSize); // size of pixel
  tempParams[10] = float(a);          // z2x
  tempParams[11] = float(0.0);        // z2y
  tempParams[12] = float(b*cos(th));  // z3x
  tempParams[13] = float(b*sin(th));  // z3y
  tempParams[14] = float(imgPixSize/float(subgridSize)); // subgrid increment

  cudaMemcpyToSymbol(params, tempParams, sizeof(float)*15);
  free(tempParams);

  // initialize x and y arrays on the host
  //for (int i = 0; i < size; i++)
  //{
  //  collectedPointsShared[i] = collectedPoints[i];
  //}

  cudaMemcpy(collectedPointsDevice,
             &collectedPoints[0],
             sizeof(float)*size,
             cudaMemcpyHostToDevice);

  // Run kernel on 1M elements on the GPU
  int threadsPerBlock = 1<<5;

  // I might easily run out of available blocks per grid.
  // Supposed size of the number of blocks is 65535.
  // https://en.wikipedia.org/wiki/Thread_block_(CUDA_programming)#Dimensions
  // Please note that Device query claims following:
  // Max dimension size of a grid size    (x,y,z): (2147483647, 65535, 65535)  
  int numBlocks = std::min(65535,(numOfPoints + threadsPerBlock - 1) / threadsPerBlock);
 

  cudaProfilerStart(); 

  arrangeShooting<<<numBlocks, threadsPerBlock>>>(collectedPointsDevice,
                                                  amps,
                                                  subgridSize,
                                                  numOfPoints);

  cudaProfilerStop();

  cudaDeviceSynchronize();

  double totalAmp = 0.0;
  for(int i = 0; i < numOfPoints; i++)
  {
    totalAmp += amps[i];
  }

  // Free memory
  cudaFree(amps);
  //cudaFree(params);
  //cudaFree(collectedPointsShared);
  cudaFree(collectedPointsDevice);

  return totalAmp/float(subgridSize*subgridSize);
};


__global__
void arrangeShooting(float*    collectedPoints,
                     float*    amps,
                     const int subgridSize,
                     const int numOfPoints)
{
  //const thrust::complex<float> z1 = thrust::complex<float>(0.0,0.0);
  const thrust::complex<float> z2 = thrust::complex<float>(params[10],params[11]);
  const thrust::complex<float> z3 = thrust::complex<float>(params[12],params[13]);
 
  // increment in image plane iteration
  //const float inc = params[9]/__int2float_rn(subgridSize);

  // actual index of a thread
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  while(index < numOfPoints) // in the case too many threads are launched
  {
    thrust::complex<float> sourcePos = thrust::complex<float>(params[7], params[8]);
    thrust::complex<float> imgPos = thrust::complex<float>(
                                                       collectedPoints[2*index]-params[9]/2,
                                                       collectedPoints[2*index+1]-params[9]/2);
  
    for (int i = 0; i < subgridSize; i++)
    { 
      imgPos.real(imgPos.real() + params[14]);
      for (int j = 0; j < subgridSize; j++)
      {
        imgPos.imag(imgPos.imag() + params[14]);
        
        amps[index] += irs(z2, z3, imgPos, sourcePos); 
      }
    }

    index += blockDim.x * gridDim.x;
  }

};


__device__
float irs(const thrust::complex<float>& z2,
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

    //return 1.0/(1.0+r*r);
};