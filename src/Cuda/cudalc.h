/*****************************************************************************
*
* Class for obtaining and printing out of Light Curve of
* Gravitational Microlens using in Inverse-Ray Shooting method
* for extended sources.
*
*
* Note: Includes a wrapper readable by Python ctypes module.
*
*****************************************************************************/


#ifndef LIGHTCURVECUDA_H
#define LIGHTCURVECUDA_H

#include<vector>

// cmath is here just for M_PI
#include<cmath>

#include<ccc.h>
#include<lcbase.h>
#include<lcirs.h>
#include<imgpoint.h>
#include<lens.h>
#include<limbdarkeningmodel.h>

#include "cudairs.cuh"

class LightCurveCUDA: public LightCurveIRS
{
  
  public:
    LightCurveCUDA(double          a,
                   double          b,
                   double          th,
                   double          m2,
                   double          m3,
                   double          sourceSize,
                   unsigned int    lcLength,
                   long int        imgPlaneSize
                  );

    ~LightCurveCUDA();

    void getLCCUDA(complex<double> startPos,
                   complex<double> endPos 
                  );

    std::vector<complex<double>> getSeeds(complex<double> pos);

    // All inverse-ray-shooting functions bool whether ray hit the source
    bool irsCheck(double imgX,
                  double imgY,
                  complex<double> sPos
                 );

    // Basic check without optimization
    bool irsCheckBase(double imgX,
                      double imgY,
                      complex<double> sPos
                     );

    // IRS optimized for Intel CPUs with AVX 256 bit registers 
    bool irsCheckSIMD(double imgX,
                      double imgY,
                      complex<double> sPos
                     );

    // IRS optimized for Intel CPUs with AVX 256 bit registers 
    bool irsCheckBinary(double imgX,
                        double imgY,
                        complex<double> sPos
                       );


    void lineFloodFillCUDA(long int        nx,
                           long int        ny,
                           complex<double> sPos,
                           bool            checked = false);

    void lineFloodFillIRS(long int        nx,
                          long int        ny,
                          complex<double> sPos,
                          bool            checked = false);

    CriticalCurveCaustic ccc;
    Amoeba amoebae;

  private:

    // Whole point of using an envelope is to make sure that 
    // thin parts of the image are included in the amoeba that
    // is given to cuda filling.
    double _envelopeRadiusSq;

    double _cpuSeeds = 0.0;
    double _cpuFloodFill = 0.0;
    double _gpuTrigger = 0.0;
    double _gpuSync = 0.0;
    double _gpuInit = 0.0;

    // Define function pointer to store pointer
    typedef bool (LightCurveCUDA::*_irsChecker)(double, double, complex<double>);
    _irsChecker _irsCheck;
};

#endif                

