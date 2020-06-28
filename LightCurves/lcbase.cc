#include<lcbase.h>

LightCurveBase::LightCurveBase(
                            double       a,
                            double       b,
                            double       th,
                            double       m2,
                            double       m3,
                            unsigned int lcLength = 50
                        ): Lens(a, b, th, m2, m3),
                           _pointImages(lensPar) 
{
  _length = lcLength;
  lcVec.resize(_length);
};

double LightCurveBase::getPointAmp(complex<double> sourcePos)
{
  vector<complex<double>> imgPos = _pointImages.getImages(sourcePos);
  double amp = 0.0;

  for(auto img: imgPos)
  {
    amp +=1.0/abs(1.0-norm(m1/pow(img-z1,2)+m2/pow(img-z2,2)+m3/pow(img-z3,2)));  
  }  

  return amp;
};

void LightCurveBase::getLC(complex<double> startPoint,
                           complex<double> endPoint
                          )
{

  complex<double> pos;
  for(unsigned int i = 0; i < _length; i++)
  {
    pos = startPoint + (endPoint-startPoint)*(i/(_length-1.0));
    lcVec[i] = getPointAmp(pos);
  }

  _hasLightCurve = true;
};
