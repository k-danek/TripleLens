#include<cudapointcollector.h>

CudaPointCollector::CudaPointCollector(double          a,
                                       double          b,
                                       double          th,
                                       double          m2,
                                       double          m3,
                                       double          sourceSize,
                                       double          pixelSize
                                      ): LensPar(a, b, th, m2, m3),
                                         _pixelSize(pixelSize)
{

  _sSize = float(sourceSize);
  
  // Source positions are initialised to 0.0 as they need to be updated before floodfilling
  _sX    = float(0.0);
  _sY    = float(0.0);

  // initializing collectedPoints
  //_collectedPoints = std::vector<float>{};

  // Allocate Unified Memory – accessible from CPU or GPU
  //cudaMallocManaged(&_collectedPoints, 2*collectionSize*sizeof(float));
  //cudaMallocManaged(&params, 8*sizeof(float));

  // sets _collectedPoints and _numberOfPoints to zero
  reset();

}

CudaPointCollector::~CudaPointCollector()
{
  _collectedPoints.clear();
}

void CudaPointCollector::setSourcePos(double sourceX,
                                      double sourceY)
{
  _sX = float(sourceX);
  _sY = float(sourceY);
}

void CudaPointCollector::addPoint(double seedX,
                                  double seedY)
{
  _collectedPoints.push_back(float(seedX));
  _collectedPoints.push_back(float(seedY));
  _numberOfPoints++;
}

int CudaPointCollector::getNumberOfPoints()
{
  return _numberOfPoints;
}

void CudaPointCollector::reset()
{
  _numberOfPoints = 0;
  _collectedPoints.clear();
}

// this is called then I collect enough points
double CudaPointCollector::getAmp()
{

  //std::cout << "getAmp called with " << _numberOfPoints << "points"; 
  
  double amp = getAmpKernel(_collectedPoints, a, b, th, m2, m3, _sSize, _sX, _sY, _pixelSize);
  reset();
  return amp;
}




