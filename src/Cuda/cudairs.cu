#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <cuda_profiler_api.h>

#include "cudairs.cuh"


__constant__ cudaFloat params[17];

SyncerCUDA::SyncerCUDA(double                     a,
                       double                     b,
                       double                     th,
                       double                     m2,
                       double                     m3,
                       double                     sourceSize,
                       double                     imgPixSize,
                       std::complex<double>       imgPlaneOrigin)
{
  _a = a;
  _b = b;
  _th = th;
  _m2 = m2;
  _m3 = m3;
  _sourceSize = sourceSize;
  _imgPixSize = imgPixSize;
  _imgPlaneOrigin = imgPlaneOrigin;
  _numOfNodes = 0;
};

double SyncerCUDA::syncAndReturn(int lcStep)
{
  if(lcStep < 0)
  {
    return 0.0;
  }

  cudaStreamSynchronize(_streamA);
  cudaStreamSynchronize(_streamB);
  cudaStreamSynchronize(_streamC);

  cudaFloat totalAmpCUDA = 0.0;
  for(int i = 0; i < _numOfNodes; i++)
  {
    totalAmpCUDA += _ampsHost[i];
  }

  std::cout << "total cuda amp =" << totalAmpCUDA/cudaFloat(subgridSize*subgridSize)*0.000146912 << "\n";

  // Free memory
  cudaFree(_ampsDeviceA);
  cudaFree(_ampsDeviceB);
  cudaFree(_ampsDeviceC);
  cudaFreeHost(_ampsHost);
  cudaFreeHost(_nodesHost);
  free(_tempParams);
  cudaFree(_nodesDeviceA);
  cudaFree(_nodesDeviceB);
  cudaFree(_nodesDeviceC);
  cudaStreamDestroy(_streamA);
  cudaStreamDestroy(_streamB);
  cudaStreamDestroy(_streamC);

  return totalAmpCUDA/cudaFloat(subgridSize*subgridSize);

}

void SyncerCUDA::trigger(amoebae_t& amoeba,
                         double     sourcePosX,
                         double     sourcePosY)
{
  cudaFloat* tempParams = (cudaFloat*)malloc(sizeof(cudaFloat)*17);

  tempParams[0]  = cudaFloat(_a);          // a
  tempParams[1]  = cudaFloat(_b);          // b
  tempParams[2]  = cudaFloat(_th);         // th
  tempParams[3]  = cudaFloat(1.0-_m2-_m3);  // m1
  tempParams[4]  = cudaFloat(_m2);         // m2
  tempParams[5]  = cudaFloat(_m3);         // m3
  tempParams[6]  = cudaFloat(_sourceSize); // sourceSize
  tempParams[7]  = cudaFloat(sourcePosX); // source position X
  tempParams[8]  = cudaFloat(sourcePosY); // source position Y
  tempParams[9]  = cudaFloat(_imgPixSize); // size of pixel
  tempParams[10] = cudaFloat(_a);          // z2x
  tempParams[11] = cudaFloat(0.0);        // z2y
  tempParams[12] = cudaFloat(_b*cos(_th));  // z3x
  tempParams[13] = cudaFloat(_b*sin(_th));  // z3y
  tempParams[14] = cudaFloat(_imgPixSize/double(subgridSize)); // subgrid increment
  tempParams[15] = cudaFloat(_imgPlaneOrigin.real()-_imgPixSize*(0.5-0.5/cudaFloat(subgridSize))); // x-origin of coordinates in image plane 
  tempParams[16] = cudaFloat(_imgPlaneOrigin.imag()-_imgPixSize*(0.5-0.5/cudaFloat(subgridSize))); // y-origin of coordinates in image plane

  cudaMemcpyToSymbol(params, tempParams, sizeof(cudaFloat)*17); 

  std::vector<Node> nodes;
  long int numOfAmoebaPoints = 0;

  _numOfNodes = 0;

  for(auto lines: amoeba)
  {
    for(auto segment: lines.second)
    {
      nodes.push_back(Node(lines.first, segment.xleft, segment.xright));
      numOfAmoebaPoints += 1 + (segment.xright - segment.xleft);     
      _numOfNodes++;
    }
  }

  // I might easily run out of available blocks per grid.
  // Supposed size of the number of blocks is 65535.
  // https://en.wikipedia.org/wiki/Thread_block_(CUDA_programming)#Dimensions
  // Please note that Device query claims following:
  // Max dimension size of a thread block (x,y,z): (1024, 1024, 64)
  // Max dimension size of a grid size    (x,y,z): (2147483647, 65535, 65535)
  // It equals number of blocks used for a single GPU execution. 
  const int numOfBlocks = 128;
  
  // Segment size is a number of nodes put into buffer for calculation with the streams.
  const int segmentSize = 3 * numOfBlocks;

  const int leftOverNodes = _numOfNodes % segmentSize;
  // Adding number of dummy nodes in order to keep 
  const int numOfNodesExtended = _numOfNodes + (segmentSize - leftOverNodes);

  std::cout << "Copied amoeba with " << numOfAmoebaPoints << " points. Segment size is " 
            << segmentSize << " .\n"
            << "Left over after segmentaion is " << numOfNodesExtended % segmentSize
            << ", number of segments is " << _numOfNodes / segmentSize << "\n";

  // creating the streams
  cudaStreamCreateWithFlags(&_streamA,cudaStreamNonBlocking);
  cudaStreamCreateWithFlags(&_streamB,cudaStreamNonBlocking);
  cudaStreamCreateWithFlags(&_streamC,cudaStreamNonBlocking);

  cudaHostAlloc((void**)& _nodesHost,
                sizeof(Node)*numOfNodesExtended,
                cudaHostAllocDefault);


  // TODO: replace this with actual one-structure directly from amoeba
  for(int i = 0; i < _numOfNodes; i++)
  {
    _nodesHost[i] = nodes[i];
  }
  for(int i = _numOfNodes; i < numOfNodesExtended; i++)
  {
    _nodesHost[i] = Node(0,0,0);
  }

  // Allocation of device buffers
  cudaMalloc((void**)&_nodesDeviceA, numOfBlocks*sizeof(Node));
  cudaMalloc((void**)&_nodesDeviceB, numOfBlocks*sizeof(Node));
  cudaMalloc((void**)&_nodesDeviceC, numOfBlocks*sizeof(Node));

  // Allocating in pinned memory
  // Always check whether these need to be initialized correctly. 
  cudaHostAlloc((void**)&_ampsHost,sizeof(float)*numOfNodesExtended, cudaHostAllocDefault);

  // Allocatin device buffers
  cudaMalloc((void**)&_ampsDeviceA, numOfBlocks*sizeof(float));
  cudaMalloc((void**)&_ampsDeviceB, numOfBlocks*sizeof(float));
  cudaMalloc((void**)&_ampsDeviceC, numOfBlocks*sizeof(float));

  // Run kernel on 1M elements on the GPU
  const int threadsPerBlock = subgridSize * subgridSize;

  // Filling the stream queues
  for(int i = 0; i < numOfNodesExtended; i += segmentSize)
  {
    // copy host -> device
    cudaMemcpyAsync(_nodesDeviceA,_nodesHost+i,              sizeof(Node)*numOfBlocks,cudaMemcpyHostToDevice,_streamA);
    cudaMemcpyAsync(_nodesDeviceB,_nodesHost+i+numOfBlocks,  sizeof(Node)*numOfBlocks,cudaMemcpyHostToDevice,_streamB);
    cudaMemcpyAsync(_nodesDeviceC,_nodesHost+i+2*numOfBlocks,sizeof(Node)*numOfBlocks,cudaMemcpyHostToDevice,_streamC);

    // initializing outputs
    cudaMemsetAsync(_ampsDeviceA,0,sizeof(float)*numOfBlocks,_streamA);
    cudaMemsetAsync(_ampsDeviceB,0,sizeof(float)*numOfBlocks,_streamB);
    cudaMemsetAsync(_ampsDeviceC,0,sizeof(float)*numOfBlocks,_streamC);

    // invoking the kernel
    arrangeShootingAmoeba<<<numOfBlocks, threadsPerBlock, 0, _streamA>>>(_nodesDeviceA,
                                                                         _ampsDeviceA,
                                                                         subgridSize,
                                                                         numOfBlocks);

    arrangeShootingAmoeba<<<numOfBlocks, threadsPerBlock, 0, _streamB>>>(_nodesDeviceB,
                                                                         _ampsDeviceB,
                                                                         subgridSize,
                                                                         numOfBlocks);

    arrangeShootingAmoeba<<<numOfBlocks, threadsPerBlock, 0, _streamC>>>(_nodesDeviceC,
                                                                         _ampsDeviceC,
                                                                         subgridSize,
                                                                         numOfBlocks);

    // copy device -> host
    cudaMemcpyAsync(_ampsHost+i              ,_ampsDeviceA,numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost,_streamA);
    cudaMemcpyAsync(_ampsHost+i+  numOfBlocks,_ampsDeviceB,numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost,_streamB);
    cudaMemcpyAsync(_ampsHost+i+2*numOfBlocks,_ampsDeviceC,numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost,_streamC);
  }
}


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

  // Each pixel will be subdivided into finer grid.
  // subgridSize determines how fine the subgrid should be.
  const int subgridSize = 8;

  cudaFloat* tempParams = (cudaFloat*)malloc(sizeof(cudaFloat)*17);

  tempParams[0]  = cudaFloat(a);          // a
  tempParams[1]  = cudaFloat(b);          // b
  tempParams[2]  = cudaFloat(th);         // th
  tempParams[3]  = cudaFloat(1.0-m2-m3);  // m1
  tempParams[4]  = cudaFloat(m2);         // m2
  tempParams[5]  = cudaFloat(m3);         // m3
  tempParams[6]  = cudaFloat(sourceSize); // sourceSize
  tempParams[7]  = cudaFloat(sourcePosX); // source position X
  tempParams[8]  = cudaFloat(sourcePosY); // source position Y
  tempParams[9]  = cudaFloat(imgPixSize); // size of pixel
  tempParams[10] = cudaFloat(a);          // z2x
  tempParams[11] = cudaFloat(0.0);        // z2y
  tempParams[12] = cudaFloat(b*cos(th));  // z3x
  tempParams[13] = cudaFloat(b*sin(th));  // z3y
  tempParams[14] = cudaFloat(imgPixSize/double(subgridSize)); // subgrid increment
  tempParams[15] = cudaFloat(imgPlaneOrigin.real()-imgPixSize*(0.5-0.5/cudaFloat(subgridSize))); // x-origin of coordinates in image plane 
  tempParams[16] = cudaFloat(imgPlaneOrigin.imag()-imgPixSize*(0.5-0.5/cudaFloat(subgridSize))); // y-origin of coordinates in image plane

  cudaMemcpyToSymbol(params, tempParams, sizeof(cudaFloat)*17);

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

  // I might easily run out of available blocks per grid.
  // Supposed size of the number of blocks is 65535.
  // https://en.wikipedia.org/wiki/Thread_block_(CUDA_programming)#Dimensions
  // Please note that Device query claims following:
  // Max dimension size of a thread block (x,y,z): (1024, 1024, 64)
  // Max dimension size of a grid size    (x,y,z): (2147483647, 65535, 65535)
  // It equals number of blocks used for a single GPU execution. 
  const int numOfBlocks = 128;
  
  // Segment size is a number of nodes put into buffer for calculation with the streams.
  const int segmentSize = 3 * numOfBlocks;

  const int leftOverNodes = numOfNodes % segmentSize;
  // Adding number of dummy nodes in order to keep 
  const int numOfNodesExtended = numOfNodes + (segmentSize - leftOverNodes);

  std::cout << "Copied amoeba with " << numOfAmoebaPoints << " points. Segment size is " 
            << segmentSize << " .\n"
            << "Left over after segmentaion is " << numOfNodesExtended % segmentSize
            << ", number of segments is " << numOfNodes / segmentSize << "\n";

  cudaStream_t streamA, streamB, streamC;
  cudaStreamCreateWithFlags(&streamA,cudaStreamNonBlocking);
  cudaStreamCreateWithFlags(&streamB,cudaStreamNonBlocking);
  cudaStreamCreateWithFlags(&streamC,cudaStreamNonBlocking);

  Node* nodesDeviceA; 
  Node* nodesDeviceB; 
  Node* nodesDeviceC; 

  Node* nodesHost;

  float *ampsDeviceA, *ampsDeviceB, *ampsDeviceC;
  float *ampsHost;

  cudaHostAlloc((void**)& nodesHost,
                sizeof(Node)*numOfNodesExtended,
                cudaHostAllocDefault);

  // TODO: replace this with actual one-structure directly from amoeba
  for(int i = 0; i < numOfNodes; i++)
  {
    nodesHost[i] = nodes[i];
  }
  for(int i = numOfNodes; i < numOfNodesExtended; i++)
  {
    nodesHost[i] = Node(0,0,0);
  }

  // Allocation of device buffers
  cudaMalloc((void**)&nodesDeviceA, numOfBlocks*sizeof(Node));
  cudaMalloc((void**)&nodesDeviceB, numOfBlocks*sizeof(Node));
  cudaMalloc((void**)&nodesDeviceC, numOfBlocks*sizeof(Node));

  // Allocating in pinned memory
  // Always check whether these need to be initialized correctly. 
  cudaHostAlloc((void**)&ampsHost,sizeof(float)*numOfNodesExtended, cudaHostAllocDefault);

  // Allocatin device buffers
  cudaMalloc((void**)&ampsDeviceA, numOfBlocks*sizeof(float));
  cudaMalloc((void**)&ampsDeviceB, numOfBlocks*sizeof(float));
  cudaMalloc((void**)&ampsDeviceC, numOfBlocks*sizeof(float));

  // Run kernel on 1M elements on the GPU
  const int threadsPerBlock = subgridSize * subgridSize;

  // Filling the stream queues
  for(int i = 0; i < numOfNodesExtended; i += segmentSize)
  {
    // copy host -> device
    cudaMemcpyAsync(nodesDeviceA,nodesHost+i,              sizeof(Node)*numOfBlocks,cudaMemcpyHostToDevice,streamA);
    cudaMemcpyAsync(nodesDeviceB,nodesHost+i+numOfBlocks,  sizeof(Node)*numOfBlocks,cudaMemcpyHostToDevice,streamB);
    cudaMemcpyAsync(nodesDeviceC,nodesHost+i+2*numOfBlocks,sizeof(Node)*numOfBlocks,cudaMemcpyHostToDevice,streamC);

    // initializing outputs
    cudaMemsetAsync(ampsDeviceA,0,sizeof(float)*numOfBlocks,streamA);
    cudaMemsetAsync(ampsDeviceB,0,sizeof(float)*numOfBlocks,streamB);
    cudaMemsetAsync(ampsDeviceC,0,sizeof(float)*numOfBlocks,streamC);

    // invoking the kernel
    arrangeShootingAmoeba<<<numOfBlocks, threadsPerBlock, 0, streamA>>>(nodesDeviceA,
                                                                          ampsDeviceA,
                                                                          subgridSize,
                                                                          numOfBlocks);

    arrangeShootingAmoeba<<<numOfBlocks, threadsPerBlock, 0, streamB>>>(nodesDeviceB,
                                                                          ampsDeviceB,
                                                                          subgridSize,
                                                                          numOfBlocks);

    arrangeShootingAmoeba<<<numOfBlocks, threadsPerBlock, 0, streamC>>>(nodesDeviceC,
                                                                          ampsDeviceC,
                                                                          subgridSize,
                                                                          numOfBlocks);

    // copy device -> host
    cudaMemcpyAsync(ampsHost+i              ,ampsDeviceA,numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost,streamA);
    cudaMemcpyAsync(ampsHost+i+  numOfBlocks,ampsDeviceB,numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost,streamB);
    cudaMemcpyAsync(ampsHost+i+2*numOfBlocks,ampsDeviceC,numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost,streamC);
}

  cudaStreamSynchronize(streamA);
  cudaStreamSynchronize(streamB);
  cudaStreamSynchronize(streamC);

  cudaFloat totalAmpCUDA = 0.0;
  for(int i = 0; i < numOfNodes; i++)
  {
    totalAmpCUDA += ampsHost[i];
  }

  std::cout << "total cuda amp =" << totalAmpCUDA/cudaFloat(subgridSize*subgridSize)*0.000146912 << "\n";

  // Free memory
  cudaFree(ampsDeviceA);
  cudaFree(ampsDeviceB);
  cudaFree(ampsDeviceC);
  cudaFreeHost(ampsHost);
  cudaFreeHost(nodesHost);
  free(tempParams);
  cudaFree(nodesDeviceA);
  cudaFree(nodesDeviceB);
  cudaFree(nodesDeviceC);
  cudaStreamDestroy(streamA);
  cudaStreamDestroy(streamB);
  cudaStreamDestroy(streamC);

  return totalAmpCUDA/cudaFloat(subgridSize*subgridSize);
};

// Variable threads per points
__global__
void arrangeShootingAmoeba(Node*     nodes,
                           float*    amps,
                           const int subgridSize,
                           const int numOfNodes)
{
  const thrust::complex<cudaFloat> z2 = thrust::complex<cudaFloat>(params[10],params[11]);
  const thrust::complex<cudaFloat> z3 = thrust::complex<cudaFloat>(params[12],params[13]);

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

  thrust::complex<cudaFloat> sourcePos = thrust::complex<cudaFloat>(params[7], params[8]);

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

  // For sm_30 there is no atomicAdd that would accept doubles.
  // Change this if you evet lay your hands on sm_60.
  atomicAdd(&amps[blockIdx.x], __double2float_rn(tempAmp));
};

__device__
cudaFloat irs(const thrust::complex<cudaFloat>& z2,
              const thrust::complex<cudaFloat>& z3,
              const thrust::complex<cudaFloat>& img,
              const thrust::complex<cudaFloat>& sourcePos)
{

    thrust::complex<cudaFloat> impact = img-params[3]/conj(img)
                                    -params[4]/conj(img-z2)
                                    -params[5]/conj(img-z3);

    cudaFloat r = thrust::abs(impact-sourcePos)/params[6];

    if(r <= 1.0)
    {
        return 0.6+0.4*sqrt(1.0-r*r); 
    }
    else
    {
        return 0.0;
    }
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