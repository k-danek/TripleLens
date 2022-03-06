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
  std::vector<Node> nodes;
  int numOfNodes = 0;
  long int numOfAmoebaPoints = 0;

  for(auto lines: amoeba)
  {
    for(auto segment: lines.second)
    {
      nodes.push_back(Node(lines.first, segment.xleft, segment.xright));
      numOfAmoebaPoints += 1 + (segment.xright - segment.xleft);     
      numOfNodes++;
    }
  }

  std::cout << "Copied amoeba with " << numOfAmoebaPoints << " points.\n";

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
  double* ampsCPU;
  float* ampsOut;
  ampsOut = (float*) malloc(sizeof(float)*numOfNodes);

  cudaMalloc((void**)&amps, numOfNodes*sizeof(float));
  
  // initialize amps
  for(int i = 0; i < numOfNodes; i++)
  {
    ampsOut[i] = 0.0;
    //ampsCPU[i] = 0.0;
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
  tempParams[14] = double(imgPixSize/double(subgridSize)); // subgrid increment
  tempParams[15] = double(imgPlaneOrigin.real()-imgPixSize*(0.5-0.5/double(subgridSize))); // x-origin of coordinates in image plane 
  tempParams[16] = double(imgPlaneOrigin.imag()-imgPixSize*(0.5-0.5/double(subgridSize))); // y-origin of coordinates in image plane

  cudaMemcpyToSymbol(params, tempParams, sizeof(double)*17);

  //cudaMemcpy(collectedPointsDevice,
  //           &collectedPoints[0],
  //           sizeof(float)*size,
  //           cudaMemcpyHostToDevice);

  // Run kernel on 1M elements on the GPU
  int threadsPerBlock = subgridSize * subgridSize;

  // I might easily run out of available blocks per grid.
  // Supposed size of the number of blocks is 65535.
  // https://en.wikipedia.org/wiki/Thread_block_(CUDA_programming)#Dimensions
  // Please note that Device query claims following:
  // Max dimension size of a thread block (x,y,z): (1024, 1024, 64)
  // Max dimension size of a grid size    (x,y,z): (2147483647, 65535, 65535)


  // Number of blocks correspons to number of nodes.
  int numBlocks = std::min(1<<31-1, numOfNodes);

  cudaProfilerStart(); 

  arrangeShootingAmoeba<<<numBlocks, threadsPerBlock>>>(nodesDevice,
                                                        amps,
                                                        subgridSize,
                                                        numOfNodes);

  cudaProfilerStop();

  cudaDeviceSynchronize();

  cudaMemcpy(ampsOut,amps,numOfNodes*sizeof(float),cudaMemcpyDeviceToHost);

  double totalAmpCUDA = 0.0;
  for(int i = 0; i < numOfNodes; i++)
  {
    totalAmpCUDA += ampsOut[i];
  }

  std::cout << "total cuda amp =" << totalAmpCUDA/float(subgridSize*subgridSize)*0.000146912 << "\n";

  // Free memory
  cudaFree(amps);
  free(ampsOut);
  //free(ampsCPU);
  free(tempParams);
  cudaFree(nodesDevice);

  return totalAmpCUDA/float(subgridSize*subgridSize);
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
  const long int gridY  = locNode.y;
  const long int gridXl = locNode.xL;
  const long int gridXr = locNode.xR; 

  // actual index of a thread
  //int threadIndex = blockIdx.x * blockDim.x + threadIdx.x;

  // Make sure these are already shifted to the bottom left corner of a pixel! 
  //
  //threadIdx.x % subgridSize; that is subgrid x
  //threadIdx.x / subgridSize; that is subgrid y
  double xShift = params[15] + __int2double_rn(threadIdx.x % subgridSize)*params[14];
  double yShift = params[16] + __ll2double_rn(gridY)*params[9] + __int2double_rn(threadIdx.x / subgridSize)*params[14];

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

double irsCPU(const double*                  params,
             //const thrust::complex<float>& z1,
             const std::complex<double>& z2,
             const std::complex<double>& z3,
             const std::complex<double>& img,
             const std::complex<double>& sourcePos)
{

    std::complex<double> impact = img-params[3]/conj(img)
                                    -params[4]/conj(img-z2)
                                    -params[5]/conj(img-z3);

    double r = std::abs(impact-sourcePos)/params[6];

    if(r <= 1.0)
    {
        return 0.6+0.4*sqrt(1.0-r*r); 
    }
    else
    {
        return 0.0;
    }
};



void arrangeShootingCPU(std::vector<Node>     nodes,
                        double*               amps,
                        double*               params,
                        const int             subgridSize,
                        const int             numOfNodes)
{
  //const thrust::complex<float> z1 = thrust::complex<float>(0.0,0.0);
  const std::complex<double> z2 = std::complex<double>(params[10],params[11]);
  const std::complex<double> z3 = std::complex<double>(params[12],params[13]);

  const int numOfThreads = subgridSize *subgridSize;

  for(int blockInd = 0; blockInd < numOfNodes; blockInd++)
  {

    // use blockIdx.x as a node index
    Node locNode = nodes[blockInd];
    long int gridY  = locNode.y;
    long int gridXl = locNode.xL;
    long int gridXr = locNode.xR; 

    double tempAmp = 0.0;

    for(int threadInd = 0; threadInd < numOfThreads; threadInd++)
    {
      double xShift = params[15] + double(threadInd % subgridSize)*params[14];
      double yShift = params[16] + double(gridY)*params[9] + double(threadInd / subgridSize)*params[14];

      std::complex<double> sourcePos = std::complex<double>(params[7], params[8]);

      for(long int gridX = gridXl; gridX <= gridXr; gridX++)
      {
        std::complex<double> imgPos = std::complex<double>(
          // origin + position of the pixel + position of subpixel
          xShift + double(gridX)*params[9],
          yShift
        );

        tempAmp += double(irsCPU(params, z2, z3, imgPos, sourcePos));
      }
    }

    amps[blockInd] += tempAmp;

  }
};