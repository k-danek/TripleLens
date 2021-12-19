#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <cuda_profiler_api.h>

#include "cudairs.cuh"


__constant__ float params[15];
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
  
  float* amps;
  //float* collectedPointsShared;
  //float* collectedPointsDevice;

  Node* nodesDevice;

  cudaMalloc( (void**)&nodesDevice, numOfNodes*sizeof(Node));

  //collectedPointsDevice = &collectedPoints[0];

  // Allocate Unified Memory – accessible from CPU or GPU
  //cudaMallocManaged(&collectedPointsShared, size*sizeof(float));
  //cudaMallocManaged(&params, 10*sizeof(float));
  cudaMallocManaged(&amps, numOfNodes*sizeof(float));

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

  //cudaMemcpy(collectedPointsDevice,
  //           &collectedPoints[0],
  //           sizeof(float)*size,
  //           cudaMemcpyHostToDevice);

  cudaMemcpy(nodesDevice,
             &nodes[0],
             sizeof(Node)*numOfNodes,
             cudaMemcpyHostToDevice);

  // Shifted by half-pixel to safe myself some GPU time
  ImgPlane imgPlane(imgPixSize,
                    real(imgPlaneOrigin)-imgPixSize/2.0,
                    imag(imgPlaneOrigin-imgPixSize/2.0)
                   );

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
                                                        imgPlane,
                                                        amps,
                                                        subgridSize,
                                                        numOfNodes);

  cudaProfilerStop();

  cudaDeviceSynchronize();

  double totalAmp = 0.0;
  for(int i = 0; i < numOfNodes; i++)
  {
    totalAmp += amps[i];
  }

  //cudaUnbindTexture(collectedPointsTexture);

  // Free memory
  cudaFree(amps);
  //cudaFree(params);
  //cudaFree(collectedPointsShared);
  //cudaFree(collectedPointsDevice);
  cudaFree(nodesDevice);

  return totalAmp/float(subgridSize*subgridSize);
};

//// One thread per subgrid point
//__global__
//void arrangeShootingThreadPerSubgrid(float*    amps,
//                                     const int subgridSize,
//                                     const int numOfPoints)
//{
//  //const thrust::complex<float> z1 = thrust::complex<float>(0.0,0.0);
//  const thrust::complex<float> z2 = thrust::complex<float>(params[10],params[11]);
//  const thrust::complex<float> z3 = thrust::complex<float>(params[12],params[13]);
// 
//  // actual index of a thread
//  int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
//
//  // sppt are subgrid points per thread
//  // const int sppt = 1;
//
//  // threadsPerPoint
//  const int threadsPerPoint = subgridSize * subgridSize;
//
//  // The index of a collected point
//  int pointIndex = (threadIndex + threadsPerPoint -1) / threadsPerPoint;
//  
//  // The index of the subbox of mutliple threads within the subgrid
//  // Now the subbox does have just one point. Can this be replaced with threadIndex?
//  const int subBoxIndex = (threadIndex + threadsPerPoint -1) % threadsPerPoint; 
//
//  // does this always work?
//  const int gridPointShift = (blockDim.x * gridDim.x) / threadsPerPoint;
//
//  while (pointIndex < numOfPoints)
//  {
//    thrust::complex<float> sourcePos = thrust::complex<float>(params[7], params[8]);
//    thrust::complex<float> imgPos = thrust::complex<float>(
//                                 tex1Dfetch(collectedPointsTexture,2*pointIndex)-params[9]/2,
//                                 tex1Dfetch(collectedPointsTexture,2*pointIndex+1)-params[9]/2
//                                                          );
//
//    thrust::complex<float> tempImgPos; 
//
//    tempImgPos.real(imgPos.real() + (subBoxIndex % subgridSize) * params[14]);
//    tempImgPos.imag(imgPos.imag() + (subBoxIndex / subgridSize) * params[14]);
//
//    atomicAdd(&amps[pointIndex], irs(z2, z3, tempImgPos, sourcePos));
//
//    pointIndex += gridPointShift;
//  }
//
//};

// Variable threads per points
__global__
void arrangeShootingAmoeba(Node*     nodes,
                           ImgPlane& imgPlane,
                           float*    amps,
                           const int subgridSize,
                           const int numOfNodes)
{
  //const thrust::complex<float> z1 = thrust::complex<float>(0.0,0.0);
  const thrust::complex<float> z2 = thrust::complex<float>(params[10],params[11]);
  const thrust::complex<float> z3 = thrust::complex<float>(params[12],params[13]);

  // use blockIdx.x as a node index
  Node locNode = nodes[blockIdx.x];
  long int gridY  = locNode.y;
  long int gridXl = locNode.xL;
  long int gridXr = locNode.xR; 

  // actual index of a thread
  //int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

  double subgridPixSize = imgPlane.step / __int2double_rn(subgridSize);
  
  // Make sure these are already shifted to the bottom left corner of a pixel! 
  //
  //threadIdx.x % subgridSize; - subgrid x
  //threadIdx.x / subgridSize; - subgrid y
  double xShift = imgPlane.originX + __int2float_rn(threadIdx.x % subgridSize)*subgridPixSize;
  double yShift = imgPlane.originY + __ll2float_rn(gridY)*params[9] + __int2float_rn(threadIdx.x / subgridSize)*subgridPixSize;

  thrust::complex<float> sourcePos = thrust::complex<float>(params[7], params[8]);

  float tempAmp = 0.0;

  for(long int gridX = gridXl; gridX <= gridXr; gridX++)
  {

    thrust::complex<float> imgPos = thrust::complex<float>(
      // origin + position of the pixel + position of subpixel
      xShift + __ll2float_rn(gridX)*params[9],
      yShift
    );

    tempAmp += irs(z2, z3, imgPos, sourcePos);
  }

  atomicAdd(&amps[blockIdx.x], tempAmp);
};


// Variable threads per points
//__global__
//void arrangeShooting(float*    amps,
//                     const int subgridSize,
//                     const int numOfPoints,
//                     const int threadsPerPoint)
//{
//  //const thrust::complex<float> z1 = thrust::complex<float>(0.0,0.0);
//  const thrust::complex<float> z2 = thrust::complex<float>(params[10],params[11]);
//  const thrust::complex<float> z3 = thrust::complex<float>(params[12],params[13]);
// 
//  // actual index of a thread
//  int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;
//
//  // sppt are subgrid points per thread
//  const int sppt = subgridSize*subgridSize/threadsPerPoint;
//
//  int pointIndex = (threadIndex + threadsPerPoint -1) / threadsPerPoint;
//  const int subBoxIndex = (threadIndex + threadsPerPoint -1) % threadsPerPoint; 
//
//  // does this always work?
//  const int gridPointShift = (blockDim.x * gridDim.x) / threadsPerPoint;
//
//  while (pointIndex < numOfPoints)
//  {
//    thrust::complex<float> sourcePos = thrust::complex<float>(params[7], params[8]);
//    thrust::complex<float> imgPos = thrust::complex<float>(
//                                 tex1Dfetch(collectedPointsTexture,2*pointIndex)-params[9]/2,
//                                 tex1Dfetch(collectedPointsTexture,2*pointIndex+1)-params[9]/2
//                                                          );
//
//    float tempAmp = 0;
//    thrust::complex<float> tempImgPos; 
//
//    // sgi for subgrid index
//    for(int sgi = subBoxIndex * sppt; sgi < (subBoxIndex + 1) * sppt; sgi++)
//    {
//      tempImgPos.real(imgPos.real() + (sgi % subgridSize) * params[14]);
//      tempImgPos.imag(imgPos.imag() + (sgi / subgridSize) * params[14]);
//      tempAmp += irs(z2, z3, tempImgPos, sourcePos); 
//    }
//
//    atomicAdd(&amps[pointIndex], tempAmp);
//
//    pointIndex += gridPointShift;
//  }
//
//};

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