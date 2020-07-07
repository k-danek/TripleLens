#include<lcirs.h>

LightCurveIRS::LightCurveIRS(
                            double       a,
                            double       b,
                            double       th,
                            double       m2,
                            double       m3,
                            double       sourceSize,
                            unsigned int lcLength = 100,
                            long int     imgPlaneSize = 100000
                            ): LightCurveBase(a, b, th, m2, m3, lcLength),
                               ccc(lensPar, 500),
                               amoebae(imgPlaneSize)
{
  _lcLength = lcLength;
  lcVec.resize(lcLength);
  _sourceRadius = sourceSize;
  ccc.getCa();
  _getCaBoxes();
  _getImgPlanePars();
  _imgPlaneSizeDouble = _topRightCornerImg.real() - _bottomLeftCornerImg.real(); 
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

  _imgPlaneSize = static_cast<long int>(halfEdge/_sourceRadius*800);
  amoebae.amoebae.resize(_imgPlaneSize);
};


void LightCurveIRS::getLCIRS(complex<double> startPoint,
                             complex<double> endPoint
                            )
{
  complex<double> pos = startPoint;
  for(unsigned int i = 0; i < _lcLength; i++)
  {

    pos = (endPoint-startPoint)*(i/(_lcLength-1.0));
    vector<complex<double>> imgPos = _pointImages.getImages(pos);
    complex<double> trialPoint;
    bool pointTaken = false;

    // This check should be put in once the issue of missing images is resolved
    //if(imgPos.size()<4 || imgPos.size() % 2 == 1 )
    //  cout << "Incorrect number of images: " << imgPos.size() << "\n";

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
        
        // looping over points on ccc for particular to get a closest points 
        for(unsigned int solInd = 0; solInd < ccc.caVec[rootInd].size()-1; solInd++)
        {
          // if a point is taken the next positive intersection will not add other.
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

    _amplification = 0.0;
    for(auto imgSeed: imgPos)
    {
      lineFloodFill(xToNx(imgSeed.real()), yToNy(imgSeed.imag()), pos);
    }
    
    lcVec.push_back(_amplification);
    
    if(_amplification > 0.0)
      cout << "non-zero amplification:" << _amplification << "\n";

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

    //double BCx=Bx-SourceRe;
    //double BCy=By-SourceIm;
    complex<double> bs = pointB - sourcePos; 
    double bsDistance = norm(bs);

    //double ACx=Ax-SourceRe;
    //double ACy=Ay-SourceIm;
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
      //   double BAx=Bx-Ax;
      //   double BAy=By-Ay;
      //   double ABdistance=BAx*BAx+BAy*BAy;
      //   double proj_times_distance=(SourceRe-Ax)*BAx+(SourceIm-Ay)*BAy;
      double abDistance = norm(ab);
      // distance on AB of projected SA times length of AB
      double projTimesDistance = as.real()*ab.real()+as.imag()*ab.imag();

      if(projTimesDistance<0)
        return 0;

      if(projTimesDistance>abDistance)
        return 0;

      if(projTimesDistance<_sourceRadius)
      {
        //TrialX=BAx*abs(proj_times_distance)/ABdistance+Ax;
        //TrialY=BAy*abs(proj_times_distance)/ABdistance+Ay;

        //double CTx=TrialX-SourceRe;
        //double CTy=TrialY-SourceIm;
        //double closest_dist=CTx*CTx+CTy*CTy;
        //if(closest_dist<SourceRadiusSq)
        // return 1;
        //else
        // return 0;
        
        // A-ab_vector*projected_length*len(AB)/len(AB)^2
        // TODO: figure more elegant type deduction here
        trialPoint = pointA-pointB*static_cast<double>(abs(projTimesDistance)
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
   double R = abs(testSourcePos-sourcePos);

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
  double posDiff = _topRightCornerImg.real() - _bottomLeftCornerImg.real();
  pos +=xFrac*posDiff;
  return pos;
}


double LightCurveIRS::nyToY(long int ny)
{
  double yFrac = ny/double(_imgPlaneSize-1.0);
  double pos = _bottomLeftCornerImg.imag();
  double posDiff = _topRightCornerImg.imag() - _bottomLeftCornerImg.imag();
  pos += yFrac*posDiff;
  return pos;
}


long int LightCurveIRS::xToNx(double x)
{
  return static_cast<long int>(
      (x - _bottomLeftCornerImg.real())/_imgPlaneSizeDouble)*_imgPlaneSize;
}


long int LightCurveIRS::yToNy(double y)
{
  return static_cast<long int>(
      (y - _bottomLeftCornerImg.imag())/_imgPlaneSizeDouble)*_imgPlaneSize;
}


double LightCurveIRS::sourceBrightness(double r)
{
  // v_factor = 0.4
  return 0.6+0.4*sqrt(1.0-r*r);
};


void LightCurveIRS::lineFloodFill(long int nx, long int ny, complex<double> sPos) {
    
    if (ny <= 0 || ny >= _imgPlaneSize) return;

    long int nL, nR, nn;

    if (!amoebae.checkLine(ny, nx)) return;

    double y = nyToY(ny), amp = irs(nxToX(nx), y, sPos); 

    if( amp <= 0.0) return;
    else _amplification += amp;

    // scan right
    for (nR = nx+1; nR < _imgPlaneSize; nR++) {
        amp = irs(nxToX(nR), y, sPos);
        
        if (amp <= 0.0) {
          nR--;
          break;
        }
        else
          _amplification += amp;
    }

    // scan left
    for (nL = nx-1; nL > 0; nL--) {
        amp = irs(nxToX(nL), y, sPos);
       
        if (amp <= 0.0) {
          nL++;
          break;
        }
        else
          _amplification += amp;
    
    }

    amoebae.addNode(nL, nR, ny);

    nn = nL;

    // trying a good position to move one row up/down
    while (nn <= nR) {
      if (!amoebae.checkLine(ny+1, nn))
        lineFloodFill(nn, ny+1, sPos);
    
      if (!amoebae.checkLine(ny-1, nn))
        lineFloodFill(nn, ny-1, sPos);
      nn++;
    }
    return;
}





