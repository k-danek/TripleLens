#include <limbdarkeningmodel.h>
#include <math.h>

// intensityCoeff is a chosen so that an integral over surface of the source with
// limb darkening is equal to the same integral over constant surface brightness equal to 1 

// Default is linear limb darkening model with a=0.6 and b=0.4
LimbDarkeningModel::LimbDarkeningModel()
{
  // contribution from one pixel is multiplied by A+B*sqrt(1-r^2)
  _linearLimbDarkeningA = 0.6;
  _linearLimbDarkeningB = 0.4;
  intensityCoeff = 1.0 / (_linearLimbDarkeningA + 2.0/3.0*_linearLimbDarkeningB);
  model = linear;
}

// Normalized linear model just takes one parameters
LimbDarkeningModel::LimbDarkeningModel(double v)
{
  // LB formula I0*(1.0-v*(1-sqrt(1.0-r*r))) 
  // Integral of the formula: PI/3*I0*(3-v)*r^2 == PI*r^2
  double i0 = 1.0 / ((1.0-v/3.0));
  _linearLimbDarkeningA = i0*(1.0-v);
  _linearLimbDarkeningB = i0*v;
  // Inverse ratio of integral of brightness profile
  intensityCoeff = 1.0;
  model = linearNormalized;
}

// Linear model just takes one parameters
LimbDarkeningModel::LimbDarkeningModel(double i0,double v)
{
  _linearLimbDarkeningA = i0*(1.0-v);
  _linearLimbDarkeningB = i0*v;
  intensityCoeff = 1.0 / (i0 * (1.0-v/3.0));
  model = linear;
}

double LimbDarkeningModel::sourceBrightness(double r)
{
  //return _OneMvFactor+_vFactor*sqrt(1.0-r*r);
  return _linearLimbDarkeningA+_linearLimbDarkeningB*sqrt(1.0-r*r);
}

