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

#include<ccc.h>
#include<lcbase.h>
#include<imgpoint.h>
#include<lens.h>
#include<amoeba.h>

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
                       complex<double> sPos);

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
    long int _imgPlaneSize;
    complex<double> _bottomLeftCornerImg = {-1.0, -1.0};
    complex<double> _topRightCornerImg = {1.0, 1.0};
    double _imgPlaneSizeDouble;
    vector<pair<complex<double>, complex<double>>> caBoxes;

    void _getCaBoxes();
 
    // sets Image plane to be square enveloping all the Critical Curve
    void _getImgPlanePars(); 

};

#endif                

