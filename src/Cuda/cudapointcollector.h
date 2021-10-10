#ifndef POINTCOLLECTOR_H
#define POINTCOLLECTOR_H

#include<lens.h>
#include "cudairs.cuh"

class CudaPointCollector: public LensPar 
{
  public:
    CudaPointCollector(double          a,
                       double          b,
                       double          th,
                       double          m2,
                       double          m3,
                       double          sourceSize,
                       double          pixelSize,
                       int             collectionSize
                      );

    // Explicitly defining destructor
    ~CudaPointCollector();

    // updates source position
    void setSourcePos(double sourceX,
                      double sourceY
                     );

    // adds a point to _collectedPoints and checks if the array is not full
    bool addPoint(double seedX,
                  double seedY 
                 );

    void reset();

    // calls cuda, synchronises the Cuda, deletes the collected points;
    double getAmp();

  private:

    // array to hold the points
    float*    _collectedPoints;

    float     _sSize;
    float     _sX;
    float     _sY;

    // number of points already collected
    int       _numberOfPoints;

    // total number of points to be collected
    const int _collectionSize;
};

#endif
