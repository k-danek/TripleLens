#include <iostream>
#include <vector>
#include <complex>
#include <math.h>
#include <unordered_map>
#include <list>

#include <thrust/complex.h>
#include <thrust/device_reference.h>

#ifndef XRANGE
#define XRANGE
/* xright - leftmost x coordinate occupied by the lensed image
   xleft - rightmost x coordinate occupied by the lensed image
*/
struct XRange
{
  long int xleft;
  long int xright;
  
  XRange(long int xLeftMost=0, long int xRightMost=0)
  {
    xleft = xLeftMost;
    xright= xRightMost;
  }

  // useful for sorting non-overlapping ranges
  bool operator<(const XRange& xRange) const
  {
      return (xleft < xRange.xleft);
  }

};
#endif

#ifndef NODE
#define NODE
struct Node 
{
    long int y;
    long int xL;
    long int xR;

    Node(long int inY, long int inL, long int inR)
    {
      y  = inY;
      xL = inL;    
      xR = inR;
    }
};
#endif

typedef std::unordered_map<long int, std::list<XRange>> amoebae_t;

// This is to be able between single and double precision
typedef float cudaFloat;

class SyncerCUDA
{
  public:
    SyncerCUDA(double                     a,
               double                     b,
               double                     th,
               double                     m2,
               double                     m3,
               double                     sourceSize,
               double                     imgPixSize,
               std::complex<double>       imgPlaneOrigin);

    SyncerCUDA(double                     a,
               double                     b,
               double                     th,
               double                     m2,
               double                     m3,
               double                     sourceSize,
               double                     imgPixSize,
               std::complex<double>       imgPlaneOrigin,
               double                     limbDarkeningA,
               double                     limbDarkeningB);

    void freeAll();

    double syncAndReturn(int lcStep);

    void trigger(amoebae_t&  amoebae,
                 double      sourcePosX,
                 double      sourcePosY);

    // This is a simple sync version to make time estimates
    float getAmpSync(amoebae_t&                 amoebae,
                     double                     a,
                     double                     b,
                     double                     th,
                     double                     m2,
                     double                     m3,
                     double                     sourceSize,
                     double                     sourcePosX,
                     double                     sourcePosY,
                     double                     imgPixSize,
                     std::complex<double>       imgPlaneOrigin);

    void allocateHost(int size);
    void allocateCuda();
    void setConstantPars();

    void printOutTimes();

    // holds index of lightcurve that is being calculated
    int        currentStep;

    // Each pixel will be subdivided into finer grid.
    // subgridSize determines how fine the subgrid should be.
    const int subgridSize = 8;

  private:
    double     _sourcePosX;
    double     _sourcePosY;
    int        _numOfNodes;
    cudaStream_t _streamA, _streamB, _streamC;
    Node *_nodesHost, *_nodesDeviceA, *_nodesDeviceB, *_nodesDeviceC;
    float *_ampsHost, *_ampsDeviceA, *_ampsDeviceB, *_ampsDeviceC;
    cudaFloat *_tempParams;

    double               _a;
    double               _b;
    double               _th;
    double               _m2;
    double               _m3;
    double               _sourceSize;
    double               _imgPixSize;
    std::complex<double> _imgPlaneOrigin;

    // Limb Darkening
    double               _linLimbDarkeningA;
    double               _linLimbDarkeningB;

    // Constant size
    int       _numberOfNodesBufferSize = 768;
    const int _numOfBlocks = 128;


    // time analysis variables
    double               _gpuMallocTime = 0.0; 
    double               _gpuMallocHostTime = 0.0; 
    double               _gpuMemsetTime = 0.0; 
    double               _gpuCopyUpTime = 0.0; 
    double               _gpuCopyDownTime = 0.0; 
    double               _gpuAmoebaTime = 0.0; 
    double               _gpuSyncTime = 0.0; 
    double               _gpuQueryTime = 0.0; 
    double               _gpuFreeTime = 0.0; 
    double               _gpuConstMemTime = 0.0; 
};

__global__
void arrangeShootingAmoeba(Node*     nodes,
                           float*    amps);

__device__
cudaFloat irs(const thrust::complex<cudaFloat>& z2,
              const thrust::complex<cudaFloat>& z3,
              const thrust::complex<cudaFloat>& img,
              const thrust::complex<cudaFloat>& sourcePos);

__device__
cudaFloat irsBinary(const thrust::complex<cudaFloat>& z2,
                    const thrust::complex<cudaFloat>& z3,
                    const thrust::complex<cudaFloat>& img,
                    const thrust::complex<cudaFloat>& sourcePos);

float irsCPU(const float*                  params,
             const thrust::complex<float>& z2,
             const thrust::complex<float>& z3,
             const thrust::complex<float>& img,
             const thrust::complex<float>& sourcePos);

void arrangeShootingCPU(std::vector<Node>     nodes,
                        double*               amps,
                        double*               params,
                        const int             subGridSize,
                        const int             numOfNodes);
