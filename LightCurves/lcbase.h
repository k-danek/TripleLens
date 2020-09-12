/*****************************************************************************
*
* Main class for obtaining and printing out of Light Curve of
* Gravitational Microlens.
*
* This is a base class with point-source amplification.
* Other classes with extended-source amplification inherit
* from this class. 
*
*
*****************************************************************************/


#ifndef LIGHTCURVEBASE_H
#define LIGHTCURVEBASE_H

#include<imgpoint.h>
#include<lens.h>

class LightCurveBase: public Lens
{
  
  public:
    LightCurveBase(double          a,
                   double          b,
                   double          th,
                   double          m2,
                   double          m3,
                   unsigned int    lcLength
                  );

    vector<double> lcVec;
    
    void getLC(complex<double> startPos,
               complex<double> endPos 
              );
    
    double getPointAmp(complex<double> sourcePos);

  protected:
    ImgPoint _pointImages;
    bool _hasLightCurve = false;
    unsigned int _lcLength;
};

#endif                

