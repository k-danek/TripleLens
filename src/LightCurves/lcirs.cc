#include<lcirs.h>
#include <string.h>

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
                               amoebae(pointsPerRadius)
{
  _lcLength = lcLength;
  _sourceRadius = sourceSize;
  _sourceRadiusSq = sourceSize*sourceSize;
  _pointsPerRadius = pointsPerRadius;
  ccc.getCa();
  _getCaBoxes();
  _getImgPlanePars();
  _limbDarkeningModel = LimbDarkeningModel();
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
  const double resizingFactor = 3.0;
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


  _imgPlaneSize = static_cast<long int>(2*halfEdge/_sourceRadius*_pointsPerRadius);
  amoebae.resize(_imgPlaneSize);

  _imgPlaneSizeDouble = 2*halfEdge;

  // Set the ampScale as number of grid point of the image plane that would fit
  // in the source radius. Also divide by a factor from source brightness
  // integration (for source of radius 1).
  //_ampScale = 1.0/M_PI/pow(_pointsPerRadius,2.0)/(1.0-_vFactor/3.0);
  _ampScale = _limbDarkeningModel.intensityCoeff/M_PI/pow(_pointsPerRadius,2.0);

  //cout << "Image plane parameters:\n";
  //cout << "_imgPlaneSize:" << _imgPlaneSize << "\n";
  //cout << "_imgPlaneSizeDouble:" << _imgPlaneSizeDouble << "\n";
  //cout << "topRightCornerImg:(" << _topRightCornerImg.real() << ","
  //     << _topRightCornerImg.imag() << ")\n"; 
  //cout << "bottomLeftCornerImg:(" << _bottomLeftCornerImg.real() << ","
  //     << _bottomLeftCornerImg.imag() << ")\n"; 

};

void LightCurveIRS::setLimbDarkeningModel(LimbDarkeningModel ldm)
{
  _limbDarkeningModel = ldm;
  // Needs to re-calculate scalling amp parameter according to limb-darkening model
  _ampScale = _limbDarkeningModel.intensityCoeff/M_PI/pow(_pointsPerRadius,2.0);

  // I should also update IRS function.
};

void LightCurveIRS::getLCIRS(complex<double> startPoint,
                             complex<double> endPoint
                            )
{

  if(_amoebaPrintOut) 
  {
    printOutLCParameters();
  }

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

    // Add up amplification from all the images/seeds
    _amplification = 0.0;
    _irsCount = 0;

    for(auto imgSeed: imgPos)
    {
      lineFloodFill(xToNx(imgSeed.real()), yToNy(imgSeed.imag()), pos);
    }

    cout << "irs amplification: " << _amplification*_ampScale << " and the count " << _irsCount
         << " and amp scale " << _ampScale << "\n";

    // As the size of the lcVec is determined at the initialisation of LightCurveIRS class
    // we use looping over the indices rather than push_back.
    lcVec[i] = _amplification*_ampScale;

    if(_amoebaPrintOut)
    {
      amoebae.printOut(_amoebaFilename+std::to_string(i));
    }

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
   double r = std::abs(testSourcePos-sourcePos);

   if( r < _sourceRadius)
   {
     double rn = r/_sourceRadius; 
     return _limbDarkeningModel.sourceBrightness(rn);
   }  
   else
   {
     return 0.0;
   }
};

double LightCurveIRS::irsOptimized(double imgX,
                                   double imgY,
                                   complex<double> sourcePos)
{
   std::complex<double> z(imgX, imgY);
   // Computationally heavy part, optimise as much as possible!
   std::complex<double> z1s = z1 - z;
   std::complex<double> z2s = z2 - z;
   std::complex<double> z3s = z3 - z;
   sourcePos -= z; 
   double rsq = std::norm(sourcePos-m1*z1s/std::norm(z1s)-m2*z2s/std::norm(z2s)-m3*z3s/std::norm(z3s));

   if( rsq < _sourceRadiusSq)
   {
     rsq = rsq/_sourceRadiusSq; 
     return _limbDarkeningModel.sourceBrightnessRSq(rsq);
   }  
   else
   {
     return 0.0;
   }
}


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
    double y = nyToY(ny), amp = irsOptimized(nxToX(nx), y, sPos); 

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
      amp = irsOptimized(nxToX(nR), y, sPos);
      
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
      amp = irsOptimized(nxToX(nL), y, sPos);
      
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

void LightCurveIRS::setAmoebaPrintOut(
                                  bool printOutAmoebae,
                                  std::string amoebaFilename,
                                  std::string parFilename)
{
  _amoebaPrintOut = printOutAmoebae;
  _amoebaFilename = amoebaFilename;
  _parFilename =    parFilename;
}

void LightCurveIRS::printOutLCParameters()
{
  std::ofstream parOutputFile;
  parOutputFile.open(_parFilename);
  std::string indent = "\t"; 

  parOutputFile << "{" << '\n';
  parOutputFile << indent << "\"sourceRadius\": " << _sourceRadius << ",\n";
  parOutputFile << indent <<"\"pointsPerRadius\": " << _pointsPerRadius << ",\n";
  parOutputFile << indent <<"\"imgPlaneSize\": " << _imgPlaneSize << ",\n";
  parOutputFile << indent <<"\"imgPlaneSizeDouble\": " << _imgPlaneSizeDouble << ",\n";
  parOutputFile << indent <<"\"bottomLeftCornerImgX\": " << _bottomLeftCornerImg.real() << ",\n";
  parOutputFile << indent <<"\"bottomLeftCornerImgY\": " << _bottomLeftCornerImg.imag()<< ",\n";
  parOutputFile << indent <<"\"topRightCornerImgX\": " << _topRightCornerImg.real() << ",\n";
  parOutputFile << indent <<"\"topRightCornerImgY\": " << _topRightCornerImg.imag()<< ",\n";

  parOutputFile << indent <<"\"a\": " << a << ",\n";
  parOutputFile << indent <<"\"b\": " << b << ",\n";
  parOutputFile << indent <<"\"th\": " << th << ",\n";
  parOutputFile << indent <<"\"m1\": " << m1 << ",\n";
  parOutputFile << indent <<"\"m2\": " << m2 << ",\n";
  parOutputFile << indent <<"\"m3\": " << m2 << ",\n";

  parOutputFile << indent <<"\"z1x\": " << z1.real() << ",\n";
  parOutputFile << indent <<"\"z1y\": " << z1.imag() << ",\n";
  parOutputFile << indent <<"\"z2x\": " << z2.real() << ",\n";
  parOutputFile << indent <<"\"z2y\": " << z2.imag() << ",\n";
  parOutputFile << indent <<"\"z3x\": " << z3.real() << ",\n";
  parOutputFile << indent <<"\"z3y\": " << z3.imag() << "\n";

  parOutputFile << "}";

  parOutputFile.close();
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

  // Api to set limb darkening model.
  void set_limb_darkening(LightCurveIRS* lc,
                          const char*    model,
                          double         v        
                         )
  {
    LimbDarkeningModel ldm;

    if(strcmp(model, "linear") == 0) 
    {
      std::cout<< "model is linear\n";
      ldm = LimbDarkeningModel(v);
    }
    else 
    {
      std::cout<< "model is default\n";
      ldm = LimbDarkeningModel();
    }

    lc->setLimbDarkeningModel(ldm);
  }

  // Api to set limb darkening model.
  void set_amoeba_printout(LightCurveIRS* lc,
                           const char*    amoebaFilename,
                           const char*    parFilename        
                          )
  {
    lc->setAmoebaPrintOut(true, (std::string)amoebaFilename, (std::string)parFilename);
  }

}
