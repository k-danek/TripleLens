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
//#include<amoeba.h>
#include<cudapointcollector.h>


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
    
    void getLCCUDA(complex<double> startPos,
                      complex<double> endPos 
                     );

    bool irsCheck(double imgX,
                  double imgY,
                  complex<double> sPos
                 );

    void lineFloodFillCUDA(long int        nx,
                           long int        ny,
                           complex<double> sPos,
                           bool            checked = false);

    CriticalCurveCaustic ccc;
    Amoeba amoebae;

  private:
    // scaling that give amplification per ray shoot
    // it takes into the account surface brightness of the source.
    // Amplification is multiplied by this factor before it is added
    // to the light curve.
    double _ampScale = 1.0;

    CudaPointCollector _cudaPointCollector;
};

#endif                

