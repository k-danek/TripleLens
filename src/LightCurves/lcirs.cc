#include<lcirs.h>

LightCurveIRS::LightCurveIRS(
                             double       a,
                             double       b,
                             double       th,
                             double       m2,
                             double       m3,
                             double       sourceSize,
                             unsigned int lcLength = 100,
                             long int     pointsPerRadius = 300
                            ): LightCurveBase(a, b, th, m2, m3, lcLength),
                               ccc(lensPar, 500),
                               amoebae(pointsPerRadius),
                               _cudaPointCollector(a,
                                                   b,
                                                   th,
                                                   m2,
                                                   m3,
                                                   sourceSize,
                                                   sourceSize/double(pointsPerRadius),
                                                   384)
{
  _lcLength = lcLength;
  _sourceRadius = sourceSize;
  _pointsPerRadius = pointsPerRadius;
  ccc.getCa();
  _getCaBoxes();
  _getImgPlanePars();
};

void LightCurveIRS::_getCaBoxes()
{
  complex<double> min = {0.0,0.0}, max = {0.0,0.0};

  for(auto rootVec: ccc.caVec)
  {
    minmax<vector<complex<double>>>(rootVec, min, max);
    caBoxes.push_back(std::make_pair(min, max));
  }
};


void LightCurveIRS::_getImgPlanePars()
{
  const double resizingFactor = 1.2;
  complex<double> min = {0.0,0.0}, max = {0.0,0.0};
  complex<double> centre;
  double halfEdge;

  for(auto rootVec: ccc.ccVec)
  {
    minmax<vector<complex<double>>>(rootVec, min, max);
  }

  halfEdge = max.real()-min.real();
  if(max.imag()-min.imag() > halfEdge)
    halfEdge = max.imag() - min.imag(); 

  halfEdge *= resizingFactor/2.0;

  centre.real((max.real() + min.real())/2.0);
  centre.imag((max.imag() + min.imag())/2.0);
  
  _topRightCornerImg.real(centre.real() + halfEdge);
  _topRightCornerImg.imag(centre.imag() + halfEdge);
  _bottomLeftCornerImg.real(centre.real() - halfEdge);
  _bottomLeftCornerImg.imag(centre.imag() - halfEdge);


  _imgPlaneSize = static_cast<long int>(halfEdge/_sourceRadius*_pointsPerRadius);
  amoebae.resize(_imgPlaneSize);

  _imgPlaneSizeDouble = 2*halfEdge;

  // Set the ampScale as number of grid point of the image plane that would fit
  // in the source radius. Also divide by a factor from source brightness
  // integration (for source of radius 1).
  _ampScale = 1.0/M_PI/pow(_pointsPerRadius,2.0)/(1.0-_vFactor/3.0);

  //cout << "Image plane parameters:\n";
  //cout << "_imgPlaneSize:" << _imgPlaneSize << "\n";
  //cout << "_imgPlaneSizeDouble:" << _imgPlaneSizeDouble << "\n";
  //cout << "topRightCornerImg:(" << _topRightCornerImg.real() << ","
  //     << _topRightCornerImg.imag() << ")\n"; 
  //cout << "bottomLeftCornerImg:(" << _bottomLeftCornerImg.real() << ","
  //     << _bottomLeftCornerImg.imag() << ")\n"; 

};


void LightCurveIRS::getLCIRS(complex<double> startPoint,
                             complex<double> endPoint
                            )
{
  
  cout << "IRS called with imgPlaneSize:" << _imgPlaneSize << "\n";
  complex<double> pos = startPoint;

  // Looping over source positions
  for(unsigned int i = 0; i <= _lcLength; i++)
  {
    //pos = (endPoint-startPoint)*(i/(_lcLength-1.0));
    pos = startPoint + (endPoint-startPoint)*(double(i)/double(_lcLength));
    cout << "started pos:" << i << ", (" << pos.real()
         << "," << pos.imag() <<")\n";
    vector<complex<double>> imgPos = _pointImages.getImages(pos);
    complex<double> trialPoint;
    bool pointTaken = false;

    // erase amoeba by initialising it again.
    amoebae = Amoeba(_imgPlaneSize);

    // this whole bit can be in an additional method
    for(unsigned int rootInd = 0; rootInd < 6; rootInd++)
    {
      // check if the point is in the caustic box
      if(pos.real()+_sourceRadius >= caBoxes[rootInd].first.real()  &&
         pos.real()-_sourceRadius <= caBoxes[rootInd].second.real() &&
         pos.imag()+_sourceRadius >= caBoxes[rootInd].first.imag()  &&
         pos.imag()-_sourceRadius <= caBoxes[rootInd].second.imag() 
        )
      {
        pointTaken = false;
        
        // looping over points on ccc to get a closest points 
        for(unsigned int solInd = 0; solInd < ccc.caVec[rootInd].size()-1; solInd++)
        {
          // If a point is taken the next positive intersection will not add other.
          // Next point then will be added ony if there the previous was not an 
          // intersection. 
          if(intersectionCheck(pos,
                               ccc.caVec[rootInd][solInd],
                               ccc.caVec[rootInd][solInd+1],
                               trialPoint))
          {
            if(!pointTaken)
            {
              for(auto trialImg: _pointImages.getImages(trialPoint))
                imgPos.push_back(trialImg);
              
              pointTaken = true;
            }
          }
          else
            pointTaken = false;
        }
      }  
    }

    //cout << "number of trial seeds: " << imgPos.size() << "\n"; 

    // Add up amplification from all the images/seeds
    _amplification = 0.0;
    _irsCount = 0;

    // Updating point collector for CUDA with new source position
    _cudaPointCollector.reset();
    _cudaPointCollector.setSourcePos(pos.real(), pos.imag());

    for(auto imgSeed: imgPos)
    {
      lineFloodFillCUDA(xToNx(imgSeed.real()), yToNy(imgSeed.imag()), pos);
    }
    cout << "amplification: " << _amplification*_ampScale << " and the count " << _irsCount << "\n";

    // As the size of the lcVec is determined at the initialisation of LightCurveIRS class
    // we use looping over the indices rather than push_back.
    lcVec[i] = _amplification*_ampScale;

  }

  _hasLightCurve = true;
};

// For two points and source ring check whether they do have intersection. If they do returns a point whose position should serve as trial for supplementary floods
// from http://doswa.com/2009/07/13/circle-segment-intersectioncollision.html
bool LightCurveIRS::intersectionCheck(complex<double>  sourcePos,
                                      complex<double>  pointA,
                                      complex<double>  pointB,
                                      complex<double>& trialPoint)
{
    double sourceRadiusSq=pow(_sourceRadius, 2.0);

    complex<double> bs = pointB - sourcePos; 
    double bsDistance = norm(bs);

    complex<double> as = pointA - sourcePos; 
    double asDistance=norm(as);

    bool aOut=(asDistance>sourceRadiusSq);
    bool bOut=(bsDistance>sourceRadiusSq);

    if(aOut != bOut)
    {
      trialPoint = bOut ? pointA : pointB;

      return true;
    }
    else if ( aOut && bOut )
    {
      complex<double> ab = pointA - pointB;   
      double abDistance = norm(ab);
      // distance on AB of projected SA times length of AB
      double projTimesDistance = as.real()*ab.real()+as.imag()*ab.imag();

      if(projTimesDistance<0)
        return 0;

      if(projTimesDistance>abDistance)
        return 0;

      if(projTimesDistance<_sourceRadius)
      {
       
        // A-ab_vector*projected_length*len(AB)/len(AB)^2
        // TODO: figure more elegant type deduction here
        trialPoint = pointA-pointB*static_cast<double>(std::abs(projTimesDistance)
                                                           /abDistance);

        complex<double> stPoint = trialPoint - sourcePos;
        double closestDist = norm(stPoint);
        
        return (closestDist < sourceRadiusSq);

      }
      else
       return 0;
    }
    else
      return 0;

  std::cout << "Error: Intersection check went all way through without arriving to"
            << " conclusion \n";
  
  return 0;

}

// Check whether source position is hit 
// In context of the point source, that source radius is just an error term.
double LightCurveIRS::irs(double imgX,
                          double imgY,
                          complex<double> sourcePos)
{
   // Computationally heavy part, optimise as much as possible!
   complex<double> img(imgX,imgY);  
   complex<double> testSourcePos=img-m1/conj(img-z1)-m2/conj(img-z2)-m3/conj(img-z3);
   double R = std::abs(testSourcePos-sourcePos);

   if( R < _sourceRadius)
   {
     double r = R/_sourceRadius; 
     return sourceBrightness(r);
   }  
   else
     return 0.0;
};


double LightCurveIRS::nxToX(long int nx)
{
  double xFrac = nx/double(_imgPlaneSize-1.0);
  double pos = _bottomLeftCornerImg.real();
  pos += xFrac*_imgPlaneSizeDouble;
  return pos;
}


double LightCurveIRS::nyToY(long int ny)
{
  double yFrac = ny/double(_imgPlaneSize-1.0);
  double pos = _bottomLeftCornerImg.imag();
  pos += yFrac*_imgPlaneSizeDouble;
  return pos;
}


long int LightCurveIRS::xToNx(double x)
{
  return static_cast<long int>(
      (x - _bottomLeftCornerImg.real())*_imgPlaneSize/_imgPlaneSizeDouble+0.5);
}


long int LightCurveIRS::yToNy(double y)
{
  return static_cast<long int>(
      (y - _bottomLeftCornerImg.imag())*_imgPlaneSize/_imgPlaneSizeDouble+0.5);
}


double LightCurveIRS::sourceBrightness(double r)
{
  return _OneMvFactor+_vFactor*sqrt(1.0-r*r);
}


void LightCurveIRS::lineFloodFill(long int nx,
                                  long int ny,
                                  complex<double> sPos,
                                  bool checked)
{
    //cout << "line floodfill run with nx " << nx << " ny " << ny << "\n";

    if(!checked)
    {
      if (ny <= 0 || ny >= _imgPlaneSize)
      {
        //cout << "Row outside image plane: " << ny << " \n";
        return;
      }

      if (!amoebae.checkLine(ny, nx))
      {
        //cout << "Amoebae check failed: " << nx << " , " << ny << " \n";
        return;
      }
    }
    // need to get the y only once as the fill stays withing a line
    double y = nyToY(ny), amp = irs(nxToX(nx), y, sPos); 

    if (amp <= 0.0) return;
    else
    {
      _amplification += amp;
      _irsCount++;
    }

    long int nL, nR, nn;

    // scan right
    for (nR = nx+1; nR < _imgPlaneSize; nR++)
    {
      amp = irs(nxToX(nR), y, sPos);
      
      if (amp <= 0.0)
      {
        nR--;
        break;
      }
      else
      {
        _amplification += amp;
        _irsCount++;
      }
    }

    // scan left
    for (nL = nx-1; nL > 0; nL--)
    {
      amp = irs(nxToX(nL), y, sPos);
      
      if (amp <= 0.0)
      {
        nL++;
        break;
      }
      else
      {
        _amplification += amp;
        _irsCount++;
      }
    }

    amoebae.addNode(nL, nR, ny);
    //cout << "got out of addNode \n";

    // trying a good position to move one row up/down
 
    nn = nL;
    // upper line
    while (nn <= nR) {
      if (amoebae.checkLineShift(ny+1, nn))
      {  
        lineFloodFill(nn, ny+1, sPos, true);
        nn++;
      }
    }
    nn = nL;
    // lower line
    while (nn <= nR) {
      if (amoebae.checkLineShift(ny-1, nn))
      {
        lineFloodFill(nn, ny-1, sPos, true);
        nn++;
      }
    }

    //cout << "finished one fill \n";

    return;
}

void LightCurveIRS::lineFloodFillCUDA(long int nx,
                                      long int ny,
                                      complex<double> sPos,
                                      bool checked)
{
    //cout << "line floodfill run with nx " << nx << " ny " << ny << "\n";

    if(!checked)
    {
      if (ny <= 0 || ny >= _imgPlaneSize)
      {
        //cout << "Row outside image plane: " << ny << " \n";
        return;
      }

      if (!amoebae.checkLine(ny, nx))
      {
        //cout << "Amoebae check failed: " << nx << " , " << ny << " \n";
        return;
      }
    }
    // need to get the y only once as the fill stays withing a line
    double x = nxToX(nx);
    double y = nyToY(ny);
    double amp = irs(x, y, sPos); 

    if (!_cudaPointCollector.addPoint(x,y))
    {
      //std::cout << "finished the point collection";
      _amplification += _cudaPointCollector.getAmp();
    };

    if (amp <= 0.0) return;
    else
    {
      //_amplification += amp;
      _irsCount++;
    }

    long int nL, nR, nn;

    // scan right
    for (nR = nx+1; nR < _imgPlaneSize; nR++)
    {
      x = nxToX(nR); 
      amp = irs(x, y, sPos);
      
      if (amp <= 0.0)
      {
        nR--;
        break;
      }
      else
      {
        //_amplification += amp;
        _irsCount++;
        if (!_cudaPointCollector.addPoint(x,y))
        {
          
          //std::cout << "finished the point collection";
          _amplification = _cudaPointCollector.getAmp();
        };
      }
    }

    // scan left
    for (nL = nx-1; nL > 0; nL--)
    {
      x = nxToX(nL); 
      amp = irs(x, y, sPos);
      
      if (amp <= 0.0)
      {
        nL++;
        break;
      }
      else
      {
        //_amplification += amp;
        _irsCount++;
        if (!_cudaPointCollector.addPoint(x,y))
        {
          //std::cout << "finished the point collection";
          _amplification += _cudaPointCollector.getAmp();
        };
      }
    }

    amoebae.addNode(nL, nR, ny);
    //cout << "got out of addNode \n";

    // trying a good position to move one row up/down
 
    nn = nL;
    // upper line
    while (nn <= nR) {
      if (amoebae.checkLineShift(ny+1, nn))
      {  
        lineFloodFillCUDA(nn, ny+1, sPos, true);
        nn++;
      }
    }
    nn = nL;
    // lower line
    while (nn <= nR) {
      if (amoebae.checkLineShift(ny-1, nn))
      {
        lineFloodFillCUDA(nn, ny-1, sPos, true);
        nn++;
      }
    }

    //cout << "finished one fill \n";

    return;
}



// Python Wrapper for ctypes module
extern "C"
{
  LightCurveIRS* lcirs_new(double       a,
                           double       b,
                           double       th,
                           double       m2,
                           double       m3,
                           double       sourceSize,
                           unsigned int lcLength,
                           long int     pointsPerRadius 
                          )
  {
    return new LightCurveIRS(a,
                             b,
                             th,
                             m2,
                             m3,
                             sourceSize,
                             lcLength,
                             pointsPerRadius
                            );
  }

  // Point source light curve
  void get_lc(LightCurveIRS* lcirs,
              double         iniX,
              double         iniY,
              double         finX,
              double         finY
             )
  {
    lcirs->getLC(complex<double>{iniX,iniY},
                 complex<double>{finX,finY}
                );
  }

  // Extended source light curve using Inverse Ray Shooting
  void get_lc_irs(LightCurveIRS* lcirs,
                  double         iniX,
                  double         iniY,
                  double         finX,
                  double         finY
                 )
  {
    lcirs->getLCIRS(complex<double>{iniX,iniY},
                    complex<double>{finX,finY}
                   );
  }

  // In order to access the data in python, 
  // we copy them to array of complex<double>
  void copy_lc(LightCurveIRS* lc,
               double*        lcArray        
              )
  {
    unsigned int length = lc->lcVec.size();

    for(unsigned int i = 0; i < length; i++)
    {
      lcArray[i] = lc->lcVec[i];
    }
  }
}

