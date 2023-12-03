#include<lcbase.h>

LightCurveBase::LightCurveBase(
                            double       a,
                            double       b,
                            double       th,
                            double       m2,
                            double       m3,
                            unsigned int lcLength = 100
                        ): Lens(a, b, th, m2, m3),
                           _pointImages(lensPar) 
{
  _lcLength = lcLength;
  lcVec.resize(_lcLength);
};

double LightCurveBase::getPointAmp(complex<double> sourcePos)
{
  vector<complex<double>> imgPos = _pointImages.getImages(sourcePos);
  double amp = 0.0;
  double ampTemp = 0.0;

  for(auto img: imgPos)
  {
    // beware of using abs, it tends to use math.h version returning ints for some reason. 
    ampTemp = 1.0/std::abs(1.0-norm(m1/pow(img-z1,2.0)+m2/pow(img-z2,2.0)+m3/pow(img-z3,2.0)));  

    // NaNs are not equal to themselves
    if(ampTemp != ampTemp)
    {
      ampTemp = 0.0;
      cout << "nan amp, printing img: [" << img.real() << "," << img.imag() << "]\n"; 
    }
    amp += ampTemp;
  }  

  return amp;
};

void LightCurveBase::getLC(complex<double> startPoint,
                           complex<double> endPoint
                          )
{
  complex<double> pos;
  for(unsigned int i = 0; i < _lcLength; i++)
  {
    pos = startPoint + (endPoint-startPoint)*(i/(_lcLength-1.0));
    lcVec[i] = getPointAmp(pos);
  }

  _hasLightCurve = true;
};
