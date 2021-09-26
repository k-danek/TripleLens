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


#ifndef LIGHTCURVEIRS_H
#define LIGHTCURVEIRS_H

#include<vector>

// cmath is here just for M_PI
#include<cmath>

#include<ccc.h>
#include<lcbase.h>
#include<imgpoint.h>
#include<lens.h>
#include<amoeba.h>
#include<cudapointcollector.h>

template<class T>
using pair = std::pair<T, T>;


class LightCurveIRS: public LightCurveBase
{
  
  public:
    LightCurveIRS(double          a,
                  double          b,
                  double          th,
                  double          m2,
                  double          m3,
                  double          sourceSize,
                  unsigned int    lcLength,
                  long int        imgPlaneSize
                 );
    
    void getLCIRS(complex<double> startPos,
                  complex<double> endPos 
                 );

    double irs(double imgX,
               double imgY,
               complex<double> sPos
              );

    // grid point to a position    
    double nxToX(long int nx);
    double nyToY(long int ny);

    // position to grid point
    long int xToNx(double x);
    long int yToNy(double y);

    double sourceBrightness(double r);

    void lineFloodFill(long int        nx,
                       long int        ny,
                       complex<double> sPos,
                       bool            checked = false);

    void lineFloodFillCUDA(long int        nx,
                           long int        ny,
                           complex<double> sPos,
                           bool            checked = false);

    bool intersectionCheck(complex<double>  pos,
                           complex<double>  pointA,
                           complex<double>  pointB,
                           complex<double>& trialPoint
                          );

    CriticalCurveCaustic ccc;
    Amoeba amoebae;

  private:
    unsigned int _lcLength;
    double _sourceRadius;
    double _amplification;
    unsigned long int _irsCount = 0;

    // parameters defining connected with image plane
    // points per radius of source - determines grid density
    int _pointsPerRadius;
    // total number of grid points per edge of the image-plane grid
    long int _imgPlaneSize;
    // lenght of the edge of the image-plane (storing for convenience)
    double _imgPlaneSizeDouble;
    // corners of the image plane
    complex<double> _bottomLeftCornerImg = {-1.0, -1.0};
    complex<double> _topRightCornerImg = {1.0, 1.0};
    // scaling that give amplification per ray shoot
    // it takes into the account surface brightness of the source.
    // Amplification is multiplied by this factor before it is added
    // to the light curve.
    double _ampScale = 1.0;

    CudaPointCollector _cudaPointCollector;

    vector<pair<complex<double>>> caBoxes;

    // constants to control surface brightness
    const double _vFactor = 0.4;
    const double _OneMvFactor = 1.0 - _vFactor;

    void _getCaBoxes();
 
    // sets Image plane to be square enveloping all the Critical Curve
    void _getImgPlanePars(); 
};

#endif                

