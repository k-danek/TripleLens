#include <iostream>
#include <vector>
#include <math.h>
//#include <nvfunctional>
#include <unordered_map>
#include <list>

#include <thrust/complex.h>
#include <thrust/device_reference.h>
//#include <unordered_map>

//#include <amoeba.h>

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

typedef std::unordered_map<long int, std::list<XRange>> amoebae_t;


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

// Kernel
float getAmpKernel(amoebae_t&                 amoebae,
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

//__global__
//void arrangeShooting(float*    amps,
//                     const int subGridSize,
//                     const int numOfPoint,
//                     const int pointsPerThread);

__global__
void arrangeShootingAmoeba(Node*     nodes,
                           float*    amps,
                           const int subGridSize,
                           const int numOfNodes);

//__global__
//void arrangeShootingThreadPerSubgrid(
//                     float*    amps,
//                     const int subGridSize,
//                     const int numOfPoint);

__device__
double irs(//const thrust::complex<float>& z1,
          const thrust::complex<double>& z2,
          const thrust::complex<double>& z3,
          const thrust::complex<double>& img,
          const thrust::complex<double>& sourcePos);

float irsCPU(const float*                  params,
             //const thrust::complex<float>& z1,
             const thrust::complex<float>& z2,
             const thrust::complex<float>& z3,
             const thrust::complex<float>& img,
             const thrust::complex<float>& sourcePos);

void arrangeShootingCPU(std::vector<Node>     nodes,
                        double*    amps,
                        double*   params,
                        const int subGridSize,
                        const int numOfNodes);


