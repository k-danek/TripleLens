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
                                 amoebae(pointsPerRadius),
                                 _cudaPointCollector(a,
                                                     b,
                                                     th,
                                                     m2,
                                                     m3,
                                                     sourceSize,
                                                     sourceSize/double(pointsPerRadius))
{
  _lcLength = lcLength;
  _sourceRadius = sourceSize;
  _pointsPerRadius = pointsPerRadius;
  ccc.getCa();
  _getCaBoxes();
  _getImgPlanePars();
  _cudaPointCollector.updateOrigin(_bottomLeftCornerImg);
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

    //cout << "number of trial seeds: " << imgPos.size() << "\n"; 

    // Add up amplification from all the images/seeds
    _amplification = 0.0;
    _irsCount = 0;

    // Updating point collector for CUDA with new source position
    _cudaPointCollector.reset();
    _cudaPointCollector.setSourcePos(pos.real(), pos.imag());

    //std::cout << "Just before the floodfill for " << imgPos.size() << " seeds\n";

    for(auto imgSeed: imgPos)
    {
      //std::cout << "Started filling for a seed with "
      //          << _cudaPointCollector.getNumberOfPoints()
      //          << " points in the collector and "
      //          << _irsCount << " in count \n";
      lineFloodFillCUDA(xToNx(imgSeed.real()), yToNy(imgSeed.imag()), pos);
    }

    //std::cout << "Finished filling for the seeds\n";

    _amplification += _cudaPointCollector.getAmp(amoebae.amoebae);

    std::cout << "Got back to lcIRS with amplification " << _amplification << "\n";

    cout << "cuda amplification: " << _amplification*_ampScale << " and the count "
         << _irsCount << " and scale " << _ampScale << "\n";
    // As the size of the lcVec is determined at the initialisation of LightCurveIRS class
    // we use looping over the indices rather than push_back.
    lcVec[i] = _amplification*_ampScale;

    std::cout << "Assigned the amp to lc amoeba size was " << amoebae.amoebae.size() << "\n";

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
   double R = std::abs(testSourcePos-sourcePos);

   if( R < _sourceRadius)
   {
     return true;
   }  

   return false;
};

void LightCurveCUDA::lineFloodFillCUDA(long int nx,
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

    if (!irsCheck(x,y,sPos)) {
      return;
    }
    else
    {
      //_amplification += amp;
      _irsCount++;
      _cudaPointCollector.addPoint(x,y);
    }

    long int nL, nR, nn;

    // scan right
    for (nR = nx+1; nR < _imgPlaneSize; nR++)
    {
      x = nxToX(nR); 
      
      if (!irsCheck(x,y,sPos))
      {
        nR--;
        break;
      }
      else
      {
        //_amplification += amp;
        _irsCount++;
        _cudaPointCollector.addPoint(x,y);
      }
    }

    // scan left
    for (nL = nx-1; nL > 0; nL--)
    {
      x = nxToX(nL); 
      
      if (!irsCheck(x,y,sPos))
      {
        nL++;
        break;
      }
      else
      {
        //_amplification += amp;
        _irsCount++;
        _cudaPointCollector.addPoint(x,y);
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
