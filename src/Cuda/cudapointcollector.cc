#include<cudapointcollector.h>

CudaPointCollector::CudaPointCollector(double          a,
                                       double          b,
                                       double          th,
                                       double          m2,
                                       double          m3,
                                       double          sourceSize,
                                       int             collectionSize
                                      ): LensPar(a, b, th, m2, m3),
                                         _collectionSize(collectionSize)
{

  _sSize = float(sourceSize);
  
  // Source positions are initialised to 0.0 as they need to be updated before floodfilling
  _sX    = float(0.0);
  _sY    = float(0.0);

  // Allocate Unified Memory – accessible from CPU or GPU
  //cudaMallocManaged(&_collectedPoints, 2*collectionSize*sizeof(float));
  //cudaMallocManaged(&params, 8*sizeof(float));

  _collectedPoints = (float *) malloc(2*_collectionSize*sizeof(float));

  // sets _collectedPoints and _numberOfPoints to zero
  reset();

}

CudaPointCollector::~CudaPointCollector()
{
  delete(_collectedPoints);
}

void CudaPointCollector::setSourcePos(double sourceX,
                                      double sourceY)
{
  _sX = float(sourceX);
  _sY = float(sourceY);
}

bool CudaPointCollector::addPoint(double seedX,
                                  double seedY)
{
  if (_numberOfPoints < _collectionSize)
  {
    int ind = 2*_numberOfPoints;
    _collectedPoints[ind] = seedX;
    _collectedPoints[ind] = seedY;
    _numberOfPoints++;
    return true;
  }
  else
  {
    return false; 
  } 
}

void CudaPointCollector::reset()
{
  _numberOfPoints = 0;

  for(int i = 0; i < 2*_collectionSize; i++)
  {
    _collectedPoints[i] = 0.0;
  
  }
}

// this is called then I collect enough points
double CudaPointCollector::getAmp()
{
  //std::cout << "getAmp called with " << _numberOfPoints << "points"; 
  return 0.0;
}




