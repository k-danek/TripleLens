#ifndef POINTCOLLECTOR_H
#define POINTCOLLECTOR_H

#include<lens.h>
#include "cudairs.cuh"

class CudaPointCollector: public LensPar 
{
  public:
    CudaPointCollector(double                a,
                       double                b,
                       double                th,
                       double                m2,
                       double                m3,
                       double                sourceSize,
                       double                pixelSize
                      );

    // Explicitly defining destructor
    ~CudaPointCollector();

    // updates source position
    void setSourcePos(double sourceX,
                      double sourceY
                     );

    void addPoint(double seedX,
                  double seedY 
                 );
    
    int getNumberOfPoints();

    void updateOrigin(std::complex<double> origin);

    void reset();

    // calls cuda, synchronises the Cuda, deletes the collected points;
    //double getAmp();
    double getAmp(amoebae_t& amoebae);

  private:

    // vector to hold the points
    std::vector<float> _collectedPoints;

    float     _sSize;
    float     _sX;
    float     _sY;

    // number of points already collected
    int       _numberOfPoints;

    // size of an pixel in image plane
    const double _pixelSize;

    std::complex<double> _imgPlaneOrigin;
};

#endif
