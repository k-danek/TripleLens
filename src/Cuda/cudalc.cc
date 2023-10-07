#include<cudalc.h>
#include <immintrin.h>

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
                                 ccc(lensPar, 1500),
                                 amoebae(pointsPerRadius)
{
  std::cout << "Initializing LCCuda";

  _lcLength = lcLength;
  _sourceRadius = sourceSize;
  _sourceRadiusSq = sourceSize*sourceSize;
  _pointsPerRadius = pointsPerRadius;
  _envelopeRadiusSq = _sourceRadiusSq * 1.4;
  ccc.getCa();
  _getCaBoxes();
  _getImgPlanePars();
  _limbDarkeningModel = LimbDarkeningModel();

  // Make spurious images check based on actual source size
  _pointImages.setSourceSize(sourceSize);
  std::cout << " Initializing LCCuda - finished\n";
};

LightCurveCUDA::~LightCurveCUDA()
{
  std::cout << "LightCurve CUDA destroyed\n";
};


std::vector<complex<double>> LightCurveCUDA::getSeeds(complex<double> pos)
{
    // Get images for the image position
    vector<complex<double>> imgPos = _pointImages.getImages(pos);

    bool pointTaken = false;
    complex<double> trialPoint;

    // check for instersections with caustic.
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
            // Take just 1 point. This fails for 1-root self-intersection when we would have
            // two intersections from just 1 root.
            if(!pointTaken)
            {
              for(auto trialImg: _pointImages.getImages(trialPoint))
                imgPos.push_back(trialImg);
              
              pointTaken = true;
            }

            // Add points on Critical curve that might have images appearing on them.
            // If there is an intersection add the surrounding points.
            imgPos.push_back(ccc.ccVec[rootInd][solInd]);
            imgPos.push_back(ccc.ccVec[rootInd][solInd+1]);
          
          }
          else
            pointTaken = false;
        }
      }  
    }

    return imgPos;
}


void LightCurveCUDA::getLCCUDA(complex<double> startPoint,
                               complex<double> endPoint
                              )
{
  clock_t beginTime, endTime;  

  cout << "IRS called with imgPlaneSize:" << _imgPlaneSize << "\n";
  cout << "IRS called with sourceRadius:" << _sourceRadius << "\n";
  cout << "IRS called with _imgPlaneSizeDouble:" << _imgPlaneSizeDouble << "\n";
  cout << "IRS called with _ampScale:" << _ampScale << "\n";
  
  complex<double> pos = startPoint;
  
  beginTime = clock(); 
  SyncerCUDA cudaSyncer(a,
                        b,
                        th,
                        m2,
                        m3,
                        _sourceRadius,
                        _imgPlaneSizeDouble/double(_imgPlaneSize-1.0),
                        _bottomLeftCornerImg,
                        _limbDarkeningModel.getParA(),
                        _limbDarkeningModel.getParB());
  endTime = clock(); 
  _gpuInit += double(endTime - beginTime);

  cout << "Syncer run" << "\n";

  // Looping over source positions
  for(unsigned int i = 0; i <= _lcLength; i++)
  {
    pos = startPoint + (endPoint-startPoint)*(double(i)/double(_lcLength));
    //cout << "started pos:" << i << ", (" << pos.real()
    //     << "," << pos.imag() <<")\n";
    
    beginTime = clock(); 
    // This includes both images of source center and source-caustic intersections
    vector<complex<double>> imgSeeds = getSeeds(pos);
    endTime = clock(); 
    _cpuSeeds += double(endTime - beginTime);

    // erase amoeba by initialising it again.
    amoebae = Amoeba(_imgPlaneSize);

    // Add up amplification from all the images/seeds
    _amplification = 0.0;
    _irsCount = 0;


    beginTime = clock(); 
    for(auto imgSeed: imgSeeds)
    {
      lineFloodFillCUDA(xToNx(imgSeed.real()), yToNy(imgSeed.imag()), pos);
    }
    endTime = clock(); 
    _cpuFloodFill += double(endTime - beginTime);

    if(i > 0) 
    {
      beginTime = clock(); 
      // sync and copy the calculation from the previous step
      _amplification = cudaSyncer.syncAndReturn(i-1);
      endTime = clock(); 
      _gpuSync += double(endTime - beginTime);
      
      //cout << "cuda amplification: " << _amplification*_ampScale << " and the count "
      //   << _irsCount << " and scale " << _ampScale << "\n";
      // As the size of the lcVec is determined at the initialisation of LightCurveIRS class
      // we use looping over the indices rather than push_back.
      lcVec[i-1] = _amplification*_ampScale;    
    };

    beginTime = clock(); 
    // Fill the GPU queue
    cudaSyncer.trigger(amoebae.amoebae,pos.real(),pos.imag());
    endTime = clock();
    _gpuTrigger += double(endTime - beginTime);
  }

  cout << "LC completed" << "\n";

  //  // Kernel
  //  float getAmpSync(amoebae_t&                 amoebae,
  //                   double                     a,
  //                   double                     b,
  //                   double                     th,
  //                   double                     m2,
  //                   double                     m3,
  //                   double                     sourceSize,
  //                   double                     sourcePosX,
  //                   double                     sourcePosY,
  //                   double                     imgPixSize,
  //                   std::complex<double>       imgPlaneOrigin);


  beginTime = clock(); 
  _amplification = cudaSyncer.syncAndReturn(_lcLength);


  cout << "LC synced" << "\n";

  endTime = clock();
  _gpuSync += double(endTime - beginTime);

  cudaSyncer.freeAll();

  cout << "Syncer freed" << "\n";

  lcVec[_lcLength] = _amplification*_ampScale;    
  //cout << "last cuda amplification: " << lcVec[_lcLength] << " and the count "
  //     << _irsCount << " and scale " << _ampScale << "\n";

  cout << "CPU Seed time: " << _cpuSeeds / CLOCKS_PER_SEC << "\n"
       << "CPU Flood Fill time: " << _cpuFloodFill / CLOCKS_PER_SEC << "\n"
       << "GPU Trigger time " << _gpuTrigger / CLOCKS_PER_SEC << "\n"
       << "GPU Sync time: " << _gpuSync / CLOCKS_PER_SEC << "\n"
       << "GPU Init time: " << _gpuInit / CLOCKS_PER_SEC << "\n"
       << "GPU of the overal time: " << (_gpuTrigger+_gpuSync)/(_cpuSeeds+_cpuFloodFill+_gpuTrigger+_gpuSync) << "\n";


  cout << "Syncer Cuda set with:" << _limbDarkeningModel.getParA() << " " << _limbDarkeningModel.getParB() << "\n";

  cudaSyncer.printOutTimes();

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
   double rSq = std::norm(testSourcePos-sourcePos);

   // Beware.Evelope radius is used instead of source radius.
   // Resulting amoeba correspons to images of the envelope size.   
   return rSq <= _envelopeRadiusSq;
};

bool LightCurveCUDA::irsCheckSIMD(double imgX,
                                  double imgY,
                                  complex<double> sourcePos)
{
  __m256d m0123 = _mm256_set_pd(1.0, m1, m2, m3);

  __m256d zetaz123R = _mm256_set_pd(sourcePos.real(), z1.real(), z2.real(), z3.real());
  __m256d zetaz123I = _mm256_set_pd(sourcePos.imag(), z1.imag(), z2.imag(), z3.imag());
  __m256d imgR      = _mm256_set1_pd(imgX);
  __m256d imgI      = _mm256_set1_pd(imgY);

  // Work with subtracteed values  
  zetaz123R = _mm256_sub_pd(zetaz123R, imgR);
  zetaz123I = _mm256_sub_pd(zetaz123I, imgI);

  // Try to get normalization

  // Multiplication, 3/4 effective
  __m256d zetaz123Rsq = _mm256_mul_pd(zetaz123R, zetaz123R);
  __m256d zetaz123Isq = _mm256_mul_pd(zetaz123I, zetaz123I);

  // Now, I will need to sum re-im pairs. 3/4 effective
  __m256d realConstants = _mm256_add_pd(zetaz123Rsq, zetaz123Isq);
  // Normalize first element after multiplication. Please note the negative signs to fix the formula
  realConstants[3] = -1.0;
  // Divide masses with obtained sum-vector. 3/4 effective
  realConstants = _mm256_div_pd(m0123, realConstants);

  // Multiply complex numbers with real-number contants.
  zetaz123R = _mm256_mul_pd(zetaz123R, realConstants);
  zetaz123I = _mm256_mul_pd(zetaz123I, realConstants);

  // https://stackoverflow.com/questions/13422747/reverse-a-avx-register-containing-doubles-using-a-single-avx-intrinsic
  // https://stackoverflow.com/questions/6996764/fastest-way-to-do-horizontal-sse-vector-sum-or-other-reduction

  __m128d zetaz123Rlo = _mm256_extractf128_pd(zetaz123R, 0); // Latenct/throughput 3/1
  __m128d zetaz123Rhi = _mm256_extractf128_pd(zetaz123R, 1); // 3/1
  __m128d zetaz123Ilo = _mm256_extractf128_pd(zetaz123I, 0); // 3/1
  __m128d zetaz123Ihi = _mm256_extractf128_pd(zetaz123I, 1); // 3/1

  double tempR[2];
  double tempI[2];

  __m128d rSumR = _mm_add_pd(zetaz123Rlo, zetaz123Rhi); // 3/1 
  __m128d rSumI = _mm_add_pd(zetaz123Ilo, zetaz123Ihi); // 3/1

  _mm_storeu_pd(tempR, rSumR); // 1/1
  _mm_storeu_pd(tempI, rSumI); // 1/1

  // distance from source
  return (tempR[0]+tempR[1])*(tempR[0]+tempR[1])+(tempI[0]+tempI[1])*(tempI[0]+tempI[1]) <= _envelopeRadiusSq;
}

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

    if (!irsCheckSIMD(x,y,sPos)) {
      return;
    }

    long int nL, nR, nn;

    // scan right
    for (nR = nx+lenghtOfStep; nR < _imgPlaneSize; nR+=lenghtOfStep)
    {
      x = nxToX(nR); 
      
      if (!irsCheckSIMD(x,y,sPos))
      {
        nR--;
        break;
      }
    }

    // scan left
    for (nL = nx-lenghtOfStep; nL > 0; nL-=lenghtOfStep)
    {
      x = nxToX(nL); 
      
      if (!irsCheckSIMD(x,y,sPos))
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

  // Api to set limb darkening model.
  void set_limb_darkening_cuda(LightCurveCUDA* lc,
                               const char*    model,
                               double         v        
                              )
  {
    LimbDarkeningModel ldm;

    if(strcmp(model, "linear") == 0) 
    {
      ldm = LimbDarkeningModel(v);
    }
    else 
    {
      ldm = LimbDarkeningModel();
    }

    lc->setLimbDarkeningModel(ldm);
  }

}
