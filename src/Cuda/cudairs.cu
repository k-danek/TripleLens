#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <cuda_profiler_api.h>

#include "cudairs.cuh"


__constant__ double params[17];
//texture<float,cudaTextureType1D,cudaReadModeElementType> collectedPointsTexture;

//const int numberOfCores = 384; 

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
float getAmpKernel(amoebae_t&                amoeba,
                   double                    a,
                   double                    b,
                   double                    th,
                   double                    m2,
                   double                    m3,
                   double                    sourceSize,
                   double                    sourcePosX,
                   double                    sourcePosY,
                   double                    imgPixSize,
                   std::complex<double>      imgPlaneOrigin)
{

  int size = amoeba.size(); 
  //int numOfPoints;

  std::vector<Node> nodes;
  int numOfNodes = 0;

  for(auto lines: amoeba)
  {
    for(auto segment: lines.second)
    {
      nodes.push_back(Node(lines.first, segment.xleft, segment.xright));
      numOfNodes++;
    }
  }

  // Each pixel will be subdivided into finer grid.
  // subgridSize determines how fine the subgrid should be.
  const int subgridSize = 8;

  Node* nodesDevice;

  cudaMalloc( (void**)&nodesDevice, numOfNodes*sizeof(Node));
  cudaMemcpy(nodesDevice,
             &nodes[0],
             sizeof(Node)*numOfNodes,
             cudaMemcpyHostToDevice);
 
  float* amps;
  float* ampsOut;
  ampsOut = (float*) malloc(sizeof(float)*numOfNodes);

  cudaMalloc((void**)&amps, numOfNodes*sizeof(float));
  
  // initialize amps
  for(int i = 0; i < numOfNodes; i++)
  {
    ampsOut[i] = 0.0;
  }
  cudaMemcpy(amps,
             ampsOut,
             sizeof(Node)*numOfNodes,
             cudaMemcpyHostToDevice);

  double* tempParams = (double*)malloc(sizeof(double)*17);

  tempParams[0]  = double(a);          // a
  tempParams[1]  = double(b);          // b
  tempParams[2]  = double(th);         // th
  tempParams[3]  = double(1.0-m2-m3);  // m1
  tempParams[4]  = double(m2);         // m2
  tempParams[5]  = double(m3);         // m3
  tempParams[6]  = double(sourceSize); // sourceSize
  tempParams[7]  = double(sourcePosX); // source position X
  tempParams[8]  = double(sourcePosY); // source position Y
  tempParams[9]  = double(imgPixSize); // size of pixel
  tempParams[10] = double(a);          // z2x
  tempParams[11] = double(0.0);        // z2y
  tempParams[12] = double(b*cos(th));  // z3x
  tempParams[13] = double(b*sin(th));  // z3y
  tempParams[14] = double(imgPixSize/float(subgridSize)); // subgrid increment
  tempParams[15] = double(imgPlaneOrigin.real()-imgPixSize/2.0); // x-origin of coordinates in image plane 
  tempParams[16] = double(imgPlaneOrigin.imag()-imgPixSize/2.0); // y-origin of coordinates in image plane

  cudaMemcpyToSymbol(params, tempParams, sizeof(double)*17);
  free(tempParams);

  //cudaMemcpy(collectedPointsDevice,
  //           &collectedPoints[0],
  //           sizeof(float)*size,
  //           cudaMemcpyHostToDevice);

  // Run kernel on 1M elements on the GPU
  int threadsPerBlock = 1<<6;

  // I might easily run out of available blocks per grid.
  // Supposed size of the number of blocks is 65535.
  // https://en.wikipedia.org/wiki/Thread_block_(CUDA_programming)#Dimensions
  // Please note that Device query claims following:

  //int numBlocks = std::min(65535, 
  //                        (numOfPoints * subgridSize * subgridSize + threadsPerBlock)
  //                        / threadsPerBlock);

  //int numBlocks = std::min(65535, 
  //                         (numOfNodes + threadsPerBlock)
  //                         / threadsPerBlock);

  // Number of blocks correspons to number of nodes.
  int numBlocks = std::min(65535, numOfNodes);

  cudaProfilerStart(); 

  arrangeShootingAmoeba<<<numBlocks, threadsPerBlock>>>(nodesDevice,
                                                        amps,
                                                        subgridSize,
                                                        numOfNodes);

  cudaProfilerStop();

  cudaDeviceSynchronize();

  cudaMemcpy(ampsOut,amps,numOfNodes*sizeof(float),cudaMemcpyDeviceToHost);

  float totalAmp = 0.0;
  for(int i = 0; i < numOfNodes; i++)
  {
    totalAmp += ampsOut[i];
  }

  // Free memory
  cudaFree(amps);
  free(ampsOut);
  //cudaFree(params);
  //cudaFree(collectedPointsShared);
  //cudaFree(collectedPointsDevice);
  cudaFree(nodesDevice);

  return totalAmp/float(subgridSize*subgridSize);
};

// Variable threads per points
__global__
void arrangeShootingAmoeba(Node*     nodes,
                           float*    amps,
                           const int subgridSize,
                           const int numOfNodes)
{
  //const thrust::complex<float> z1 = thrust::complex<float>(0.0,0.0);
  const thrust::complex<double> z2 = thrust::complex<double>(params[10],params[11]);
  const thrust::complex<double> z3 = thrust::complex<double>(params[12],params[13]);

  // use blockIdx.x as a node index
  Node locNode = nodes[blockIdx.x];
  long int gridY  = locNode.y;
  long int gridXl = locNode.xL;
  long int gridXr = locNode.xR; 

  // actual index of a thread
  //int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

  // Make sure these are already shifted to the bottom left corner of a pixel! 
  //
  //threadIdx.x % subgridSize; - subgrid x
  //threadIdx.x / subgridSize; - subgrid y
  double xShift = params[15] + __int2double_rn(threadIdx.x % subgridSize)*params[9];
  double yShift = params[16] + __ll2double_rn(gridY)*params[9] + __int2double_rn(threadIdx.x / subgridSize)*params[9];

  thrust::complex<double> sourcePos = thrust::complex<double>(params[7], params[8]);

  double tempAmp = 0.0;

  for(long int gridX = gridXl; gridX <= gridXr; gridX++)
  {

    thrust::complex<double> imgPos = thrust::complex<double>(
      // origin + position of the pixel + position of subpixel
      xShift + __ll2double_rn(gridX)*params[9],
      yShift
    );

    tempAmp += irs(z2, z3, imgPos, sourcePos);
  }

  atomicAdd(&amps[blockIdx.x], __double2float_rn(tempAmp));
};

__device__
double irs(const thrust::complex<double>& z2,
           const thrust::complex<double>& z3,
           const thrust::complex<double>& img,
           const thrust::complex<double>& sourcePos)
{

    thrust::complex<double> impact = img-params[3]/conj(img)
                                    -params[4]/conj(img-z2)
                                    -params[5]/conj(img-z3);

    double r = thrust::abs(impact-sourcePos)/params[6];

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