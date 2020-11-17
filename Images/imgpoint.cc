#include <imgpoint.h>
#include <imgpointcoeff.cc>

// Constructor uses constructor of Lens class
ImgPoint::ImgPoint(
                   double       a,
                   double       b,
                   double       th,
                   double       m2,
                   double       m3,
                   double       posX = 0.0,
                   double       posY = 0.0
                  ): Lens(a, b, th, m2, m3)
{
  setPos(posX,posY);
};

ImgPoint::ImgPoint(const LensPar &lensParam): Lens(lensParam.a,
                                                   lensParam.b,
                                                   lensParam.th,
                                                   lensParam.m2,
                                                   lensParam.m3
                                                  )
{
  setPos(0.0, 0.0);
};


// Normal check if an image falls within a source radius. 
// In context of the point source, that source radius is just an error term.
bool ImgPoint::imgCheck(complex<double> img, double sourceRadius)
{
   complex<double> testSourcePos=img-m1/conj(img-z1)-m2/conj(img-z2)-m3/conj(img-z3);

   if(abs(testSourcePos-_sourcePos)<(sourceRadius))
   {
     return 1;
   } 
   else
   {
     return 0;
   }  
};

// Sets the image position
void ImgPoint::setPos(double posX,
                      double posY
                     )
{
  _sourcePos      = {posX, posY};  
  // Invalidates the previous calculation!
  _rootsAvailable = false;
  _imgsAvailable  = false;
};

// Sets the image position - overload
void ImgPoint::setPos(complex<double> pos)
{
  _sourcePos      = pos;  
  // Invalidates the previous calculation!
  _rootsAvailable = false;
  _imgsAvailable  = false;
};

// Main method to initialise vector with critical curve
void ImgPoint::getRoots(bool forceNewRoots = false)
{
  vector<complex<double>> imgCoef = getCoeffs();

  //Beware! This is complex conjugate of the Critical Curve
  Laguerre laguerre(imgCoef);

  //This part decides whether to polish previous or calculate new roots.
  //Doing just polishing should speed up the process cca twice.
  //Also, it leaves the roots in the same order as in the previous step.
  //<10 condition is there to check whether tempRoots holds all 10 roots    
  if(_tempRoots.size() < 10 || forceNewRoots)
  {
    _tempRoots = laguerre.solveRoots();
    //cout << "Roots recalculated \n";
  }
  else
  {
    _tempRoots = laguerre.solveRoots(_tempRoots);  

    if(!laguerre.checkRoots(_tempRoots))
    {
      cout << "Roots off for img " << "\n";
      _tempRoots = laguerre.solveRoots();
    }        
  };  

  // no cleaning of roots as we just assing the roots
  roots = _tempRoots;

  // Calculation finished, the critical curve is now available.
  _rootsAvailable = true;
};


// Given all the image candidates, we establish which one is a real root.
void ImgPoint::getImages()
{
  // clears all the elements of the vector
  imgs.clear();
  isImg.clear();
  double err = 1.0e-6;

  if(!_rootsAvailable)
  {
    getRoots();
  }

  for(auto root: roots)
  {
    if(imgCheck(root, err))
    {
      imgs.push_back(root);
      isImg.push_back(true);
    }  
    else
      isImg.push_back(false);
  }  

  // This is an important check for correct number of images
  // Unfortunatelly, this fails a way too often!
  // Come back to this after full testing of other functionality. 
  if(imgs.size() % 2 != 0 || imgs.size() < 4)
  { 
    
    getRoots(true);
    // clears all the elements of the vector
    imgs.clear();
    isImg.clear();
    
    for(auto root: roots)
    {
      if(imgCheck(root, err))
      {
        imgs.push_back(root);
        isImg.push_back(true);
      }  
      else
        isImg.push_back(false);
    }    
    
    if(imgs.size() % 2 != 0 || imgs.size() < 4)
    { 
      vector<complex<double>>  coeffs = getCoeffsOpt();
      cout << "Wrong number of images after the check: " << imgs.size() << "\n";

      for(auto coeff: coeffs)
      {
        cout << std::scientific << std::setw(4) << "[" << coeff.real() 
             << "," << coeff.imag() << "]";
      }
      cout << "\n"; 
    }
  }
  _imgsAvailable = true;

};

vector<complex<double>> ImgPoint::getImages(complex<double> pos)
{
  setPos(pos);
  getImages();
  return imgs;
};  

// Python Wrapper for ctypes module
extern "C"
{
  ImgPoint* img_new(double a,
                    double b,
                    double th,
                    double m2,
                    double m3,
                    double posX,
                    double posY
                   )
  {
    return new ImgPoint(a,
                        b,
                        th,
                        m2,
                        m3,
                        posX,
                        posY
                       );
  };

  void get_roots(ImgPoint* img)
  {
    img->getRoots();
  };

  void get_images(ImgPoint* img)
  {
    img->getImages();
  };

  void set_pos(ImgPoint* img,
               double    posX,
               double    posY)
  {
    img->setPos(posX, posY);
  }

  // In order to access the data in python, 
  // we copy them to array of complex<double>
  void copy_images(ImgPoint*         img,
                   complex<double>*  roots,
                   bool*             isImg)
  {
    unsigned int length = img->roots.size();

    for(unsigned int i = 0; i < length; i++)
    {
      roots[i] = img->roots[i];
      isImg[i] = img->isImg[i];
    }
  }

}


