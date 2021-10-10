#include <iostream>
#include <math.h>
//#include <nvfunctional>

#include <thrust/complex.h>
#include <thrust/device_reference.h>

// Kernel
float getAmpKernel(float* collectedPoints,
                   double a,
                   double b,
                   double th,
                   double m2,
                   double m3,
                   double sourceSize,
                   double sourcePosX,
                   double sourcePosY);

__global__
void arrangeShooting(float* collectionPoints,
                      float* params,
                      float* amps);

__device__
float irs(const float*                  params,
          //const thrust::complex<float>& z1,
          const thrust::complex<float>& z2,
          const thrust::complex<float>& z3,
          const thrust::complex<float>& img,
          const thrust::complex<float>& sourcePos);


