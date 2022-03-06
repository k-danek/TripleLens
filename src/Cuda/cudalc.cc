#include<cudalc.h>

LightCurveCUDA::LightCurveCUDA(
                             double       a,
                             double       b,
                             double       th,
                             double       m2,
                             double       m3,
                             double       sourceSize,
                             unsigned int lcLength = 100,
                             long int     pointsPerRadius = 300
                              ): LightCurveIRS(a, b, th, m2, m3, sourceSize, lcLength, pointsPerRadius),
                                 ccc(lensPar, 500),
                                 amoebae(pointsPerRadius)
{
  _lcLength = lcLength;
  _sourceRadius = sourceSize;
  _pointsPerRadius = pointsPerRadius;
  _envelopeRadius = sourceSize * 1.1;
  ccc.getCa();
  _getCaBoxes();
  _getImgPlanePars();
};


void LightCurveCUDA::getLCCUDA(complex<double> startPoint,
                               complex<double> endPoint
                              )
{
  
  cout << "IRS called with imgPlaneSize:" << _imgPlaneSize << "\n";
  cout << "IRS called with sourceRadius:" << _sourceRadius << "\n";
  cout << "IRS called with _imgPlaneSizeDouble:" << _imgPlaneSizeDouble << "\n";
  cout << "IRS called with _ampScale:" << _ampScale << "\n";
  
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

    // Add up amplification from all the images/seeds
    _amplification = 0.0;
    _irsCount = 0;

    for(auto imgSeed: imgPos)
    {
      lineFloodFillCUDA(xToNx(imgSeed.real()), yToNy(imgSeed.imag()), pos);
    }

    _amplification = getAmpKernel(amoebae.amoebae,
                                   a,
                                   b,
                                   th,
                                   m2,
                                   m3,
                                   _sourceRadius,
                                   pos.real(),
                                   pos.imag(),
                                   _imgPlaneSizeDouble/double(_imgPlaneSize-1.0),
                                   _bottomLeftCornerImg);

    cout << "cuda amplification: " << _amplification*_ampScale << " and the count "
         << _irsCount << " and scale " << _ampScale << "\n";
    // As the size of the lcVec is determined at the initialisation of LightCurveIRS class
    // we use looping over the indices rather than push_back.
    lcVec[i] = _amplification*_ampScale;
  }

  _hasLightCurve = true;
};

// Check whether source position is hit 
// In context of the point source, that source radius is just an error term.
bool LightCurveCUDA::irsCheck(double imgX,
                              double imgY,
                              complex<double> sourcePos)
{
   // Computationally heavy part, optimise as much as possible!
   complex<double> img(imgX,imgY);  
   complex<double> testSourcePos=img-m1/conj(img-z1)-m2/conj(img-z2)-m3/conj(img-z3);
   double r = std::abs(testSourcePos-sourcePos);

   // Beware.Evelope radius is used instead of source radius.
   // Resulting amoeba correspons to images of the envelope size.   
   if( r <= _envelopeRadius)
   {
     return true;
   }  

   return false;
};

void LightCurveCUDA::lineFloodFillIRS(long int nx,
                                  long int ny,
                                  complex<double> sPos,
                                  bool checked)
{
    if(!checked)
    {
      if (ny <= 0 || ny >= _imgPlaneSize)
      {
        return;
      }

      if (!amoebae.checkLine(ny, nx))
      {
        return;
      }
    }
    // need to get the y only once as the fill stays withing a line
    double y = nyToY(ny), amp = irs(nxToX(nx), y, sPos); 

    if (amp <= 0.0)
    {
      return;
    }
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

    // trying a good position to move one row up/down
    nn = nL;
    // upper line
    while (nn <= nR) {
      if (amoebae.checkLineShift(ny+1, nn))
      {  
        lineFloodFillIRS(nn, ny+1, sPos, true);
        nn++;
      }
    }
    nn = nL;
    // lower line
    while (nn <= nR) {
      if (amoebae.checkLineShift(ny-1, nn))
      {
        lineFloodFillIRS(nn, ny-1, sPos, true);
        nn++;
      }
    }

    return;
}

void LightCurveCUDA::lineFloodFillCUDA(long int nx,
                                       long int ny,
                                       complex<double> sPos,
                                       bool checked)
{
    const int lenghtOfStep = 5;

    if(!checked)
    {
      if (ny <= 0 || ny >= _imgPlaneSize)
      {
        cout << "Row outside image plane: " << ny << " \n";
        return;
      }

      if (!amoebae.checkLine(ny, nx))
      {
        return;
      }
    }
    // need to get the y only once as the fill stays withing a line
    double x = nxToX(nx);
    double y = nyToY(ny);

    if (!irsCheck(x,y,sPos)) {
      return;
    }

    long int nL, nR, nn;

    // scan right
    for (nR = nx+lenghtOfStep; nR < _imgPlaneSize; nR+=lenghtOfStep)
    {
      x = nxToX(nR); 
      
      if (!irsCheck(x,y,sPos))
      {
        nR--;
        break;
      }
    }

    // scan left
    for (nL = nx-lenghtOfStep; nL > 0; nL-=lenghtOfStep)
    {
      x = nxToX(nL); 
      
      if (!irsCheck(x,y,sPos))
      {
        nL++;
        break;
      }

    }

    amoebae.addNode(nL, nR, ny);
    _irsCount += nR - nL + 1;

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

    return;
}


// Python Wrapper for ctypes module
extern "C"
{
  LightCurveCUDA* lccuda_new(double       a,
                             double       b,
                             double       th,
                             double       m2,
                             double       m3,
                             double       sourceSize,
                             unsigned int lcLength,
                             long int     pointsPerRadius 
                            )
  {
    return new LightCurveCUDA(a,
                              b,
                              th,
                              m2,
                              m3,
                              sourceSize,
                              lcLength,
                              pointsPerRadius
                             );
  }

  // Extedned source calculation with cuda
  void get_lc_cuda(LightCurveCUDA* lccuda,
                  double          iniX,
                  double          iniY,
                  double          finX,
                  double          finY
                  )
  {
    lccuda->getLCCUDA(complex<double>{iniX,iniY},
                      complex<double>{finX,finY}
                     );
  }

  // In order to access the data in python, 
  // we copy them to array of complex<double>
  void copycuda_lc(LightCurveCUDA* lc,
                   double*         lcArray        
                  )
  {
    unsigned int length = lc->lcVec.size();

    for(unsigned int i = 0; i < length; i++)
    {
      lcArray[i] = lc->lcVec[i];
    }
  }
}
