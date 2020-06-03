#include <ccc.h>

// Constructor uses constructor of Lens class
CriticalCurveCaustic::CriticalCurveCaustic(
                            double       a,
                            double       b,
                            double       th,
                            double       m2,
                            double       m3,
                            unsigned int cccLength = 500
                        ): Lens(
                             a,
                             b,
                             th,
                             m2,
                             m3
                           )
{
  _length = cccLength;
  z2c = conj(z2);
  z3c = conj(z3);
};


// Main method to initialise vector with critical curve
void CriticalCurveCaustic::getCC()
{

  complex<double> cc[7];
  complex<double> eiphi;

  // Memory allocation as we know how many roots we will have
  ccVec.resize(6);
  for(auto vec: ccVec) vec.reserve(_length);

  cc[1]=2.0*z2c*z3c*(z2c+z3c)*(m2+m3-1.0);
  cc[0]=-z2c*z2c*z3c*z3c*(m2+m3-1.0);
  
  for (unsigned int j=0; j< _length ;j++)
  {
    double phi=2.0*3.14159265359*j/_length;
    
    eiphi = complex<double>(cos(phi),sin(phi));

    cc[6]=-eiphi;
    cc[5]=2.0*eiphi*(z2c+z3c);
    cc[4]=1.0-eiphi*(z2c*z2c+z3c*z3c+4.0*z2c*z3c);
    cc[3]=2.0*eiphi*z2c*z3c*(z2c+z3c)+2.0*z2c*(m2-1.0)+2.0*z3c*(m3-1.0);
    cc[2]=-eiphi*z2c*z2c*z3c*z3c+4.0*(1.0-m3-m2)*z2c*z3c-z2c*z2c*(m2-1.0)-z3c*z3c*(m3-1.0);

    vector<complex<double>> ccCoef (cc, cc + sizeof(cc) / sizeof(complex<double>) );

    //Beware! This is complex conjugate of the Critical Curve
    Laguerre laguerre(ccCoef);

    //This part decides whether to polish previous or calculate new roots.
    //Doing just polishing should speed up the process cca twice.
    //Also, it leaves the roots in the same order as in the previous step.
    //<6 condition is there to check whether    
    if(_tempRoots.size() < 6)
    {
      _tempRoots = laguerre.solveRoots(); 
    }
    else
    {
      _tempRoots = laguerre.polishRoots(_tempRoots);   
      if(!laguerre.checkRoots(_tempRoots))
      {
        cout << "Roots off for j = " << j << "\n";
        _tempRoots = laguerre.solveRoots(); 
      }        
    };   

    for(unsigned int k=0; k<6; k++)
    {
      ccVec[k].push_back(_tempRoots[k]);
    } 
  }

  // Calculation finished, the critical curve is now available.
  _ccAvailable = true;
};

void CriticalCurveCaustic::getCa()
{
  complex<double> ccSol; 

  // Caustic is obrained from Critical Curve, so make sure you have one.
  if(!_ccAvailable)
  {
    getCC();
  }

  // Memory allocation as we know how big our vectors are goint to be.
  caVec.resize(6);
  for(auto vec: caVec) vec.reserve(_length);

  // looping over 6 ccVec roots
  for(unsigned int i = 0; i < 6; i++)
  {
    // looping over length solutions for particular vec
    for(unsigned int j = 0; j < _length; j++)
    {
      ccSol = ccVec[i][j];
      caVec[i].push_back(conj(ccSol-(1.0-m2-m3)/conj(ccSol)
                                   -m2/conj((ccSol-z2c))
                                   -m3/conj((ccSol-z3c))));
    }
  }
  // Caustic is now available
  _caAvailable = true;
};

void CriticalCurveCaustic::printCCC(std::string fileName)
{
  std::ofstream outFile;
  outFile.open(fileName);

  outFile.precision(6);
  // Print Header
  outFile << "#Critical Curve and Caustic\n";
  outFile << "# (a,b,theta,m2,m3) = (" << a << ", " << b << ", " << th << ", "
          << m2 << ", " << m3 << ")\n";

  // Print Critical Curve
  for(unsigned int i = 0; i < _length; i++)
  {
    for(unsigned int j = 0; j < 6; j++)
    {
      outFile << ccVec[j][i].real() << " ";
      outFile << ccVec[j][i].imag() << " ";
    }  
    outFile << "\n";
  }  

  // Print Caustic
  outFile << "\n###########\n";
  for(unsigned int i = 0; i < _length; i++)
  {
    for(unsigned int j = 0; j < 6; j++)
    {
      outFile << caVec[j][i].real() << " ";
      outFile << caVec[j][i].imag() << " ";
    }  
    outFile << "\n";
  }  
  outFile.close();
};



// Python Wrapper for ctypes module
extern "C"
{
   //Foo* Foo_new(int n) {return new Foo(n);}
   //void Foo_bar(Foo* foo) {foo->bar();}
   //int Foo_foobar(Foo* foo, int n) {return foo->foobar(n);}
  CriticalCurveCaustic* ccc_new(double a,
                                double b,
                                double th,
                                double m2,
                                double m3,
                                int    cccLength = 500
                               )
  {
    return new CriticalCurveCaustic(a,
                                    b,
                                    th,
                                    m2,
                                    m3,
                                    cccLength
                                   );
  };

  void get_cc(CriticalCurveCaustic* ccc)
  {
    ccc->getCC();
  };

  void get_ca(CriticalCurveCaustic* ccc)
  {
    ccc->getCa();
  };

  void print_ccc(CriticalCurveCaustic* ccc, char* fileName)
  {
    ccc->printCCC(fileName);
  };

  // In order to access the data in python, 
  // we copy them to array of complex<double>
  void copy_cc_ca(CriticalCurveCaustic* ccc,
                  complex<double>*      cc,
                  complex<double>*      ca)
  {
    unsigned int length = ccc->ccVec[0].size();
    unsigned int idx = 0;

    for(unsigned int root = 0; root < 6; root++)
      for(unsigned int step = 0; step < length; step++)
      {
        idx = root * length + step;
        cc[idx] = ccc->ccVec[root][step];
        ca[idx] = ccc->caVec[root][step];
      }
  }

  void copy_lenses(CriticalCurveCaustics* ccc,
                   complex<double>*       z1,
                   complex<double>*       z2,
                   complex<double>*       z3
                  )
  {
    z1 = ccc->z1;
    z2 = ccc->z2;
    z3 = ccc->z3;
  }


  void getBoundingBox(CriticalCurveCaustic* ccc,
                      complex<double>*      ccMin,
                      complex<double>*      ccMax,
                      complex<double>*      caMin,
                      complex<double>*      caMax
                     )
  {
    for(auto ccRoot: ccc->ccVec)
      minmax<vector<complex<double>>>(ccRoot, *ccMin, *ccMax);

    for(auto caRoot: ccc->caVec)
      minmax<vector<complex<double>>>(caRoot, *caMin, *caMax);
  }


}


