#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <cuda_profiler_api.h>

#include "cudairs.cuh"


__constant__ cudaFloat params[18];
__constant__ cudaFloat sourcePosParams[2];
// Subgrid size
__constant__ int       sgSize[1];

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
  allocateHost(_numberOfNodesBufferSize);
  allocateCuda();
  setConstantPars();
  _linLimbDarkeningA = 0.6;
  _linLimbDarkeningB = 0.8;
  if(m3 == 0.0)
  {
    _isBinary = true;
    _invokeKernel = &SyncerCUDA::_invokeKernelDouble;
  } else {
    _isBinary = false;
    _invokeKernel = &SyncerCUDA::_invokeKernelTriple;
  }
};

SyncerCUDA::SyncerCUDA(double                     a,
                       double                     b,
                       double                     th,
                       double                     m2,
                       double                     m3,
                       double                     sourceSize,
                       double                     imgPixSize,
                       std::complex<double>       imgPlaneOrigin,
                       double                     limbDarkeningA,
                       double                     limbDarkeningB)
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
  _linLimbDarkeningA = limbDarkeningA;
  _linLimbDarkeningB = limbDarkeningB;
  if(m3 == 0.0)
  {
    _isBinary = true;
    _invokeKernel = &SyncerCUDA::_invokeKernelDouble;
  } else {
    _isBinary = false;
    _invokeKernel = &SyncerCUDA::_invokeKernelTriple;
  }
  allocateHost(_numberOfNodesBufferSize);
  allocateCuda();
  setConstantPars();
};

void SyncerCUDA::setConstantPars()
{
  clock_t beginTime, endTime; 
  beginTime = clock();

  _tempParams[0]  = cudaFloat(_a);          // a
  _tempParams[1]  = cudaFloat(_b);          // b
  _tempParams[2]  = cudaFloat(_th);         // th
  _tempParams[3]  = cudaFloat(1.0-_m2-_m3);  // m1
  _tempParams[4]  = cudaFloat(_m2);         // m2
  _tempParams[5]  = cudaFloat(_m3);         // m3
  _tempParams[6]  = cudaFloat(_sourceSize); // sourceSize
  _tempParams[7]  = cudaFloat(_imgPixSize);  // size of pixel
  _tempParams[8]  = cudaFloat(_a);           // z2x
  _tempParams[9]  = cudaFloat(0.0);          // z2y
  _tempParams[10] = cudaFloat(_b*cos(_th));  // z3x
  _tempParams[11] = cudaFloat(_b*sin(_th));  // z3y
  _tempParams[12] = cudaFloat(_imgPixSize/double(subgridSize)); // subgrid increment
  _tempParams[13] = cudaFloat(_imgPlaneOrigin.real()-_imgPixSize*(0.5-0.5/cudaFloat(subgridSize))); // x-origin of coordinates in image plane 
  _tempParams[14] = cudaFloat(_imgPlaneOrigin.imag()-_imgPixSize*(0.5-0.5/cudaFloat(subgridSize))); // y-origin of coordinates in image plane
  _tempParams[15] = cudaFloat(_linLimbDarkeningA);
  _tempParams[16] = cudaFloat(_linLimbDarkeningB/_sourceSize);
  _tempParams[17] = cudaFloat(_sourceSize*_sourceSize);

  cudaMemcpyToSymbol(params, _tempParams, sizeof(cudaFloat)*18);
  // putting subgridsize to constant memory
  cudaMemcpyToSymbol(sgSize, (const void*)&subgridSize, sizeof(int)); 
  endTime = clock();
  _gpuConstMemTime += double(endTime-beginTime);
};

void SyncerCUDA::freeAll()
{
  if(_ampsHost != nullptr) {
    cudaFreeHost(_ampsHost);
    _ampsHost = nullptr;
  }

  if(_nodesHost != nullptr) {
    cudaFreeHost(_nodesHost);
    _nodesHost = nullptr;
  }

  if(_ampsDeviceA != nullptr) {
    cudaFree(_ampsDeviceA);
    _ampsDeviceA = nullptr;
  }

  if(_ampsDeviceB != nullptr) {
    cudaFree(_ampsDeviceB);
    _ampsDeviceB = nullptr;
  }

  if(_ampsDeviceC != nullptr) {
    cudaFree(_ampsDeviceC);
    _ampsDeviceC = nullptr;
  }

  if(_nodesDeviceA != nullptr) {
    cudaFree(_nodesDeviceA);
    _nodesDeviceA = nullptr;
  }
  
  if(_nodesDeviceB != nullptr) {
    cudaFree(_nodesDeviceB);
    _nodesDeviceB = nullptr;
  }
  
  if(_nodesDeviceC != nullptr) {
    cudaFree(_nodesDeviceC);
    _nodesDeviceC = nullptr;
  }

  if(_tempParams != nullptr) {
    free(_tempParams);
    _tempParams = nullptr;
  }
}

void SyncerCUDA::printOutTimes()
{
  //std::cout << "From CUDA part: \n"
  //          << "GPU malloc time:" << _gpuMallocTime / CLOCKS_PER_SEC  << "\n"
  //          << "GPU malloc host time:" << _gpuMallocHostTime / CLOCKS_PER_SEC  << "\n"
  //          << "GPU amoeba time:" << _gpuAmoebaTime / CLOCKS_PER_SEC  << "\n"
  //          << "GPU const mem time:" << _gpuConstMemTime / CLOCKS_PER_SEC  << "\n"
  //          << "GPU sync time:" << _gpuSyncTime / CLOCKS_PER_SEC  << "\n"
  //          << "GPU free time:" << _gpuFreeTime / CLOCKS_PER_SEC  << "\n"
  //          << "GPU query time:" << _gpuQueryTime / CLOCKS_PER_SEC  << "\n"
  //          << "GPU copy up time:" << _gpuCopyUpTime / CLOCKS_PER_SEC  << "\n"
  //          << "GPU copy down time:" << _gpuCopyDownTime / CLOCKS_PER_SEC  << "\n";
}

void SyncerCUDA::allocateHost(int size)
{
  // size should be something like numberOfNodesExtended
  _numberOfNodesBufferSize = size;

  clock_t beginTime, endTime; 

  beginTime = clock();
  cudaHostAlloc((void**)& _nodesHost,
                sizeof(Node)*size,
                cudaHostAllocDefault);
  endTime = clock();
  _gpuMallocTime += double(endTime-beginTime);

  // Allocating in pinned memory
  // Always check whether these need to be initialized correctly. 
  cudaHostAlloc((void**)&_ampsHost,sizeof(float)*size, cudaHostAllocDefault);

  beginTime = clock();
  _tempParams = (cudaFloat*)malloc(sizeof(cudaFloat)*15);
  endTime = clock();
  _gpuConstMemTime += double(endTime-beginTime);
};

void SyncerCUDA::allocateCuda()
{
  // Allocate bufferts for CUDA. The size is fixed and set up in constants.

  clock_t beginTime, endTime; 

  beginTime = clock();
  // Allocation of device buffers
  cudaMalloc((void**)&_nodesDeviceA, _numOfBlocks*sizeof(Node));
  cudaMalloc((void**)&_nodesDeviceB, _numOfBlocks*sizeof(Node));
  cudaMalloc((void**)&_nodesDeviceC, _numOfBlocks*sizeof(Node));

  // Allocatin device buffers
  cudaMalloc((void**)&_ampsDeviceA, _numOfBlocks*sizeof(float));
  cudaMalloc((void**)&_ampsDeviceB, _numOfBlocks*sizeof(float));
  cudaMalloc((void**)&_ampsDeviceC, _numOfBlocks*sizeof(float));
  endTime = clock();
  _gpuMallocTime += double(endTime-beginTime);
};

double SyncerCUDA::syncAndReturn(int lcStep)
{
  //std::cout << "Sync and return\n";

  clock_t beginTime, endTime; 

  if(lcStep < 0)
  {
    return 0.0;
  }

  beginTime = clock();
  cudaStreamSynchronize(_streamA);
  cudaStreamSynchronize(_streamB);
  cudaStreamSynchronize(_streamC);
  endTime = clock();
  _gpuSyncTime += double(endTime-beginTime);

  //std::cout << "Streams synchronised\n";

  cudaFloat totalAmpCUDA = 0.0;
  for(int i = 0; i < _numOfNodes; i++)
  {
    totalAmpCUDA += _ampsHost[i];
  }

  //std::cout << "total cuda amp =" << totalAmpCUDA/cudaFloat(subgridSize*subgridSize)*0.000146912 << "\n";

  beginTime = clock();
  cudaStreamDestroy(_streamA);
  cudaStreamDestroy(_streamB);
  cudaStreamDestroy(_streamC);
  endTime = clock();
  _gpuFreeTime += double(endTime-beginTime);
  
  return totalAmpCUDA/cudaFloat(subgridSize*subgridSize);
}

void SyncerCUDA::trigger(amoebae_t& amoeba,
                         double     sourcePosX,
                         double     sourcePosY)
{

  clock_t beginTime, endTime; 

  beginTime = clock();

  cudaFloat tempSourcePos[2];
  tempSourcePos[0] = cudaFloat(sourcePosX);
  tempSourcePos[1] = cudaFloat(sourcePosY);

  cudaMemcpyToSymbol(sourcePosParams, tempSourcePos, sizeof(cudaFloat)*2); 
  endTime = clock();
  _gpuConstMemTime += double(endTime-beginTime);

  //std::cout << "Constant memory copied over \n";

  std::vector<Node> nodes;
  long int numOfAmoebaPoints = 0;

  _numOfNodes = 0;

  beginTime = clock();
  for(auto lines: amoeba)
  {
    for(auto segment: lines.second)
    {
      nodes.push_back(Node(lines.first, segment.xleft, segment.xright));
      numOfAmoebaPoints += 1 + (segment.xright - segment.xleft);     
      _numOfNodes++;
    }
  }
  endTime = clock();
  _gpuAmoebaTime += double(endTime-beginTime);


  //std::cout << "Nodes assigned \n";

  // I might easily run out of available blocks per grid.
  // Supposed size of the number of blocks is 65535.
  // https://en.wikipedia.org/wiki/Thread_block_(CUDA_programming)#Dimensions
  // Please note that Device query claims following:
  // Max dimension size of a thread block (x,y,z): (1024, 1024, 64)
  // Max dimension size of a grid size    (x,y,z): (2147483647, 65535, 65535)
  // It equals number of blocks used for a single GPU execution. 
  //const int numOfBlocks = 128;
  
  // Segment size is a number of nodes put into buffer for calculation with the streams.
  const int segmentSize = 3 * _numOfBlocks;

  const int leftOverNodes = _numOfNodes % segmentSize;
  // Adding number of dummy nodes in order to keep 
  const int numOfNodesExtended = _numOfNodes + (segmentSize - leftOverNodes);

  //std::cout << "Copied amoeba with " << numOfAmoebaPoints << " points. Segment size is " 
  //          << segmentSize << " .\n"
  //          << "Left over after segmentaion is " << numOfNodesExtended % segmentSize
  //          << ", number of segments is " << _numOfNodes / segmentSize << "\n";

  // Pinned memory is allocated as a buffer. If the needed size is higher, realocate the host memory!
  if(numOfNodesExtended > _numberOfNodesBufferSize)
  {
    allocateHost(numOfNodesExtended);
  }

  // creating the streams
  cudaStreamCreateWithFlags(&_streamA,cudaStreamNonBlocking);
  cudaStreamCreateWithFlags(&_streamB,cudaStreamNonBlocking);
  cudaStreamCreateWithFlags(&_streamC,cudaStreamNonBlocking);

  //std::cout << "Streams created \n";

  beginTime = clock();
  // TODO: replace this with actual one-structure directly from amoeba
  for(int i = 0; i < _numOfNodes; i++)
  {
    _nodesHost[i] = nodes[i];
  }
  for(int i = _numOfNodes; i < numOfNodesExtended; i++)
  {
    _nodesHost[i] = Node(0,0,0);
  }
  endTime = clock();
  _gpuAmoebaTime += double(endTime-beginTime);
  
  //std::cout << "Everything allocated\n";


  // Run kernel on 1M elements on the GPU
  const int threadsPerBlock = subgridSize * subgridSize;

  beginTime = clock();
  // Filling the stream queues
  for(int i = 0; i < numOfNodesExtended; i += segmentSize)
  {
    // copy host -> device
    cudaMemcpyAsync(_nodesDeviceA,_nodesHost+i,               sizeof(Node)*_numOfBlocks,cudaMemcpyHostToDevice,_streamA);
    cudaMemcpyAsync(_nodesDeviceB,_nodesHost+i+_numOfBlocks,  sizeof(Node)*_numOfBlocks,cudaMemcpyHostToDevice,_streamB);
    cudaMemcpyAsync(_nodesDeviceC,_nodesHost+i+2*_numOfBlocks,sizeof(Node)*_numOfBlocks,cudaMemcpyHostToDevice,_streamC);

    // initializing outputs
    cudaMemsetAsync(_ampsDeviceA,0,sizeof(float)*_numOfBlocks,_streamA);
    cudaMemsetAsync(_ampsDeviceB,0,sizeof(float)*_numOfBlocks,_streamB);
    cudaMemsetAsync(_ampsDeviceC,0,sizeof(float)*_numOfBlocks,_streamC);

    // invoking the kernel
    (this->*_invokeKernel)(threadsPerBlock);

    // copy device -> host
    cudaMemcpyAsync(_ampsHost+i               ,_ampsDeviceA,_numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost,_streamA);
    cudaMemcpyAsync(_ampsHost+i+  _numOfBlocks,_ampsDeviceB,_numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost,_streamB);
    cudaMemcpyAsync(_ampsHost+i+2*_numOfBlocks,_ampsDeviceC,_numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost,_streamC);
  }
  endTime = clock();
  _gpuQueryTime += double(endTime-beginTime);
  //std::cout << "Queues filled\n";

}

void SyncerCUDA::_invokeKernelDouble(int threadsPerBlock)
{
  arrangeShootingAmoebaBinary<<<_numOfBlocks, threadsPerBlock, 0, _streamA>>>(_nodesDeviceA,
                                                                              _ampsDeviceA);

  arrangeShootingAmoebaBinary<<<_numOfBlocks, threadsPerBlock, 0, _streamB>>>(_nodesDeviceB,
                                                                              _ampsDeviceB);

  arrangeShootingAmoebaBinary<<<_numOfBlocks, threadsPerBlock, 0, _streamC>>>(_nodesDeviceC,
                                                                              _ampsDeviceC);
}

void SyncerCUDA::_invokeKernelTriple(int threadsPerBlock)
{
  arrangeShootingAmoeba<<<_numOfBlocks, threadsPerBlock, 0, _streamA>>>(_nodesDeviceA,
                                                                        _ampsDeviceA);

  arrangeShootingAmoeba<<<_numOfBlocks, threadsPerBlock, 0, _streamB>>>(_nodesDeviceB,
                                                                        _ampsDeviceB);

  arrangeShootingAmoeba<<<_numOfBlocks, threadsPerBlock, 0, _streamC>>>(_nodesDeviceC,
                                                                        _ampsDeviceC);
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
float SyncerCUDA::getAmpSync(
  amoebae_t&                amoeba,
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

  clock_t beginTime, endTime; 

  // Each pixel will be subdivided into finer grid.
  // subgridSize determines how fine the subgrid should be.
  const int subgridSize = 8;

  beginTime = clock();

  cudaFloat* tempParams = (cudaFloat*)malloc(sizeof(cudaFloat)*18);

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
  tempParams[17] = cudaFloat(sourceSize*sourceSize); // y-origin of coordinates in image plane

  cudaMemcpyToSymbol(params, tempParams, sizeof(cudaFloat)*18);

  endTime = clock();
  _gpuConstMemTime += double(endTime-beginTime);


  beginTime = clock();

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
  endTime = clock();
  _gpuAmoebaTime += double(endTime-beginTime);

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

  //cudaStream_t streamA, streamB, streamC;
  //cudaStreamCreateWithFlags(&streamA,cudaStreamNonBlocking);
  //cudaStreamCreateWithFlags(&streamB,cudaStreamNonBlocking);
  //cudaStreamCreateWithFlags(&streamC,cudaStreamNonBlocking);

  Node* nodesDeviceA; 
  Node* nodesDeviceB; 
  Node* nodesDeviceC; 

  Node* nodesHost;

  float *ampsDeviceA, *ampsDeviceB, *ampsDeviceC;
  float *ampsHost;

  beginTime = clock();
  cudaHostAlloc((void**)& nodesHost,
                sizeof(Node)*numOfNodesExtended,
                cudaHostAllocDefault);

  endTime = clock();
  _gpuMallocHostTime += double(endTime-beginTime);

  // TODO: replace this with actual one-structure directly from amoeba
  for(int i = 0; i < numOfNodes; i++)
  {
    nodesHost[i] = nodes[i];
  }
  for(int i = numOfNodes; i < numOfNodesExtended; i++)
  {
    nodesHost[i] = Node(0,0,0);
  }

  beginTime = clock();
  // Allocation of device buffers
  cudaMalloc((void**)&nodesDeviceA, numOfBlocks*sizeof(Node));
  cudaMalloc((void**)&nodesDeviceB, numOfBlocks*sizeof(Node));
  cudaMalloc((void**)&nodesDeviceC, numOfBlocks*sizeof(Node));
  endTime = clock();
  _gpuMallocTime += double(endTime-beginTime);


  beginTime = clock();
  // Allocating in pinned memory
  // Always check whether these need to be initialized correctly. 
  cudaHostAlloc((void**)&ampsHost,sizeof(float)*numOfNodesExtended, cudaHostAllocDefault);
  endTime = clock();
  _gpuMallocHostTime += double(endTime-beginTime);

  // Allocatin device buffers
  cudaMalloc((void**)&ampsDeviceA, numOfBlocks*sizeof(float));
  cudaMalloc((void**)&ampsDeviceB, numOfBlocks*sizeof(float));
  cudaMalloc((void**)&ampsDeviceC, numOfBlocks*sizeof(float));
  endTime = clock();
  _gpuMallocTime += double(endTime-beginTime);

  // Run kernel on 1M elements on the GPU
  const int threadsPerBlock = subgridSize * subgridSize;

  // Filling the stream queues
  for(int i = 0; i < numOfNodesExtended; i += segmentSize)
  {
 
    beginTime = clock();
    // copy host -> device
    cudaMemcpy(nodesDeviceA,nodesHost+i,              sizeof(Node)*numOfBlocks,cudaMemcpyHostToDevice);
    cudaMemcpy(nodesDeviceB,nodesHost+i+numOfBlocks,  sizeof(Node)*numOfBlocks,cudaMemcpyHostToDevice);
    cudaMemcpy(nodesDeviceC,nodesHost+i+2*numOfBlocks,sizeof(Node)*numOfBlocks,cudaMemcpyHostToDevice);
    endTime = clock();
    _gpuCopyUpTime += double(endTime-beginTime);

    beginTime = clock();
    // initializing outputs
    cudaMemset(ampsDeviceA,0,sizeof(float)*numOfBlocks);
    cudaMemset(ampsDeviceB,0,sizeof(float)*numOfBlocks);
    cudaMemset(ampsDeviceC,0,sizeof(float)*numOfBlocks);
    endTime = clock();
    _gpuMemsetTime += double(endTime-beginTime);
    
    beginTime = clock();
    // invoking the kernel
    arrangeShootingAmoeba<<<numOfBlocks, threadsPerBlock, 0>>>(nodesDeviceA,
                                                               ampsDeviceA);

    arrangeShootingAmoeba<<<numOfBlocks, threadsPerBlock, 0>>>(nodesDeviceB,
                                                               ampsDeviceB);

    arrangeShootingAmoeba<<<numOfBlocks, threadsPerBlock, 0>>>(nodesDeviceC,
                                                               ampsDeviceC);

    cudaDeviceSynchronize(); 
    endTime = clock();
    _gpuSyncTime += double(endTime-beginTime);

    beginTime = clock();
    // copy device -> host
    cudaMemcpy(ampsHost+i              ,ampsDeviceA,numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(ampsHost+i+  numOfBlocks,ampsDeviceB,numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost);
    cudaMemcpy(ampsHost+i+2*numOfBlocks,ampsDeviceC,numOfBlocks*sizeof(float),cudaMemcpyDeviceToHost);
    endTime = clock();
    _gpuCopyDownTime += double(endTime-beginTime);

  }

  beginTime = clock();
  cudaFloat totalAmpCUDA = 0.0;
  for(int i = 0; i < numOfNodes; i++)
  {
    totalAmpCUDA += ampsHost[i];
  }
  endTime = clock();
  _gpuCopyDownTime += double(endTime-beginTime);
  
  std::cout << "total cuda amp =" << totalAmpCUDA/cudaFloat(subgridSize*subgridSize)*0.000146912 << "\n";


  beginTime = clock();
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
  endTime = clock();
  _gpuFreeTime += double(endTime-beginTime);

  return totalAmpCUDA/cudaFloat(subgridSize*subgridSize);
};

// Variable threads per points
__global__
void arrangeShootingAmoeba(Node*     nodes,
                           float*    amps)
{
  const thrust::complex<cudaFloat> z2 = thrust::complex<cudaFloat>(params[8],params[9]);
  const thrust::complex<cudaFloat> z3 = thrust::complex<cudaFloat>(params[10],params[11]);

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
  double xShift = params[13] + __int2double_rn(threadIdx.x % sgSize[0])*params[12];
  double yShift = params[14] + __ll2double_rn(gridY)*params[7] + __int2double_rn(threadIdx.x / sgSize[0])*params[12];

  thrust::complex<cudaFloat> sourcePos = thrust::complex<cudaFloat>(sourcePosParams[0], sourcePosParams[1]);

  double tempAmp = 0.0;

  for(long int gridX = gridXl; gridX <= gridXr; gridX++)
  {
    thrust::complex<double> imgPos = thrust::complex<double>(
      // origin + position of the pixel + position of subpixel
      xShift + __ll2double_rn(gridX)*params[7],
      yShift
    );

    tempAmp += irs(z2, z3, imgPos, sourcePos);
  }

  // For sm_30 there is no atomicAdd that would accept doubles.
  // Change this if you evet lay your hands on sm_60.
  atomicAdd(&amps[blockIdx.x], __double2float_rn(tempAmp));
  //atomicAdd(&amps[blockIdx.x], __double2float_rn(sourcePosParams[0]));
};

// Dobule lens version
// TODO: Remove code duplication without cost to performance
__global__
void arrangeShootingAmoebaBinary(Node*     nodes,
                                 float*    amps)
{
  const thrust::complex<cudaFloat> z2 = thrust::complex<cudaFloat>(params[8],params[9]);
  const thrust::complex<cudaFloat> z3 = thrust::complex<cudaFloat>(params[10],params[11]);

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
  double xShift = params[13] + __int2double_rn(threadIdx.x % sgSize[0])*params[12];
  double yShift = params[14] + __ll2double_rn(gridY)*params[7] + __int2double_rn(threadIdx.x / sgSize[0])*params[12];

  thrust::complex<cudaFloat> sourcePos = thrust::complex<cudaFloat>(sourcePosParams[0], sourcePosParams[1]);

  double tempAmp = 0.0;

  for(long int gridX = gridXl; gridX <= gridXr; gridX++)
  {
    thrust::complex<double> imgPos = thrust::complex<double>(
      // origin + position of the pixel + position of subpixel
      xShift + __ll2double_rn(gridX)*params[7],
      yShift
    );

    tempAmp += irsBinary(z2, imgPos, sourcePos);
  }

  atomicAdd(&amps[blockIdx.x], __double2float_rn(tempAmp));
};


__device__
cudaFloat irs(const thrust::complex<cudaFloat>& z2,
              const thrust::complex<cudaFloat>& z3,
              const thrust::complex<cudaFloat>& img,
              const thrust::complex<cudaFloat>& sourcePos)
{
    thrust::complex<cudaFloat> rC = img-sourcePos-params[3]/conj(img)
                                    -params[4]/conj(img-z2)
                                    -params[5]/conj(img-z3);

    cudaFloat rSq = (rC.real()*rC.real()+rC.imag()*rC.imag());

    cudaFloat step = cudaFloat(rSq<=params[17]);
    return (params[15]+params[16]*sqrt(params[17]-rSq*step))*step;
};

__device__
cudaFloat irsBinary(const thrust::complex<cudaFloat>& z,
                    const thrust::complex<cudaFloat>& img,
                    const thrust::complex<cudaFloat>& sourcePos)
{
    thrust::complex<cudaFloat> rC = img-sourcePos-params[3]/conj(img)
                                    -params[4]/conj(img-z);

    cudaFloat rSq = (rC.real()*rC.real()+rC.imag()*rC.imag());

    cudaFloat step = cudaFloat(rSq<=params[17]);
    return (params[15]+params[16]*sqrt(params[17]-rSq*step))*step;
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