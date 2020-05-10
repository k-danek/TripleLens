#include <imgpoint.h>

// Constructor uses constructor of Lens class
ImgPoint::ImgPoint(
                   double       a,
                   double       b,
                   double       th,
                   double       m2,
                   double       m3,
                   double       posX,
                   double       posY
                  ): Lens(
                       a,
                       b,
                       th,
                       m2,
                       m3
                     )
{
  setPos(posX,posY);
};


// Normal check if an image falls within a source radius. 
// In context of the point source, that source radius is just an error term.
bool ImgPoint::imgCheck(complex<double> img, double sourceRadius)
{
   complex<double> testSourcePos=_sourcePos-(1.0-m2-m3)/conj(_sourcePos-z1)-m2/conj(_sourcePos-z2)-m3/conj(_sourcePos-z3);

   if(abs(testSourcePos-_sourcePos)<(sourceRadius))
    return 1;
   else
    return 0;
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
void ImgPoint::getRoots()
{
  vector<complex<double>> imgCoef = getCoeffs();

  //Beware! This is complex conjugate of the Critical Curve
  Laguerre laguerre(imgCoef);

  //This part decides whether to polish previous or calculate new roots.
  //Doing just polishing should speed up the process cca twice.
  //Also, it leaves the roots in the same order as in the previous step.
  //<10 condition is there to check whether    
  if(_tempRoots.size() < 10)
  {
    _tempRoots = laguerre.solveRoots();
  }
  else
  {
    _tempRoots = laguerre.polishRoots(_tempRoots);  
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
  double err = 1.0e-8;

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

};

vector<complex<double>> ImgPoint::getCoeffs()
{
  //complex<double>         z2c = z2;
  complex<double>         z3c = conj(z3);
  complex<double>         zeta = _sourcePos;
  complex<double>         zetac= conj(zeta);
  vector<complex<double>> coeffs(11);

  coeffs[10] = (-zetac * z3c + zetac * zetac) * z2 + zetac * zetac * z3c - pow(zetac, 0.3e1);
  coeffs[9] = (-0.3e1 * zetac * zetac + 0.3e1 * zetac * z3c) * z2 * z2 + ((-0.3e1 * zetac * zetac + 0.3e1 * zetac * z3c) * z3 + (-m3 - 0.3e1 * zetac * zetac - m2 + zeta * zetac) * z3c + zetac * (0.3e1 * zetac * zetac - zeta * zetac + 0.1e1 + m2)) * z2 + (-0.3e1 * zetac * zetac * z3c + 0.3e1 * pow(zetac, 0.3e1)) * z3 + zetac * (-zeta * zetac + 0.1e1 + m3) * z3c + zetac * zetac * (zeta * zetac - 0.2e1);
  coeffs[8] = (-0.3e1 * zetac * z3c + 0.3e1 * zetac * zetac) * pow(z2, 0.3e1) + ((-0.9e1 * zetac * z3c + 0.9e1 * zetac * zetac) * z3 + (0.3e1 * zetac * zetac - 0.3e1 * zeta * zetac + 0.3e1 * m3 + 0.2e1 * m2) * z3c - zetac * (0.3e1 * zetac * zetac - 0.3e1 * zeta * zetac + 0.3e1 + m2)) * z2 * z2 + ((-0.3e1 * zetac * z3c + 0.3e1 * zetac * zetac) * z3 * z3 + ((0.9e1 * zetac * zetac + 0.3e1 * m2 - 0.3e1 * zeta * zetac + 0.2e1 * m3) * z3c - zetac * (0.9e1 * zetac * zetac - 0.3e1 * zeta * zetac + 0.3e1 - 0.2e1 * m3 + 0.3e1 * m2)) * z3 + (0.3e1 * zeta * zetac * zetac + 0.2e1 * zetac * m2 - 0.3e1 * zetac - 0.3e1 * zetac * m3 + zeta) * z3c + m2 + 0.6e1 * zetac * zetac - 0.3e1 * zetac * zetac * m2 - 0.3e1 * zeta * pow(zetac, 0.3e1) - 0.2e1 * zeta * zetac) * z2 + (0.3e1 * zetac * zetac * z3c - 0.3e1 * pow(zetac, 0.3e1)) * z3 * z3 + (-zetac * (-0.3e1 * zeta * zetac + 0.3e1 + m3) * z3c - 0.3e1 * zetac * zetac * (zeta * zetac - 0.2e1 + m3)) * z3 + (m3 - 0.2e1 * zeta * zetac) * z3c + zetac * (0.3e1 * zeta * zetac - 0.1e1);
  coeffs[7] = (-zetac * zetac + zetac * z3c) * pow(z2, 0.4e1) + ((-0.9e1 * zetac * zetac + 0.9e1 * zetac * z3c) * z3 + (0.3e1 * zeta * zetac - 0.3e1 * m3 - m2 - zetac * zetac) * z3c - zetac * (-zetac * zetac + 0.3e1 * zeta * zetac - 0.3e1 + m2)) * pow(z2, 0.3e1) + ((-0.9e1 * zetac * zetac + 0.9e1 * zetac * z3c) * z3 * z3 + ((0.9e1 * zeta * zetac - 0.6e1 * m3 - 0.9e1 * zetac * zetac - 0.6e1 * m2) * z3c + 0.3e1 * zetac * (0.3e1 * zetac * zetac - 0.3e1 * zeta * zetac + m2 - 0.2e1 * m3 + 0.3e1)) * z3 + (0.3e1 * zetac * m3 + zeta * m2 - 0.3e1 * zeta - 0.4e1 * zetac * m2 + 0.3e1 * zetac - 0.3e1 * zeta * zetac * zetac) * z3c + 0.6e1 * zeta * zetac - 0.2e1 * m2 - 0.2e1 * zeta * zetac * m2 + 0.6e1 * zetac * zetac * m2 + m2 * m2 - 0.6e1 * zetac * zetac + 0.3e1 * zeta * pow(zetac, 0.3e1)) * z2 * z2 + ((-zetac * zetac + zetac * z3c) * pow(z3, 0.3e1) + ((0.3e1 * zeta * zetac - 0.9e1 * zetac * zetac - m3 - 0.3e1 * m2) * z3c + zetac * (0.9e1 * zetac * zetac - 0.3e1 * zeta * zetac + 0.3e1 - 0.4e1 * m3 + 0.3e1 * m2)) * z3 * z3 + ((-0.9e1 * zeta * zetac * zetac - 0.3e1 * zeta + 0.3e1 * zetac * m3 + 0.9e1 * zetac - 0.6e1 * zetac * m2 + zeta * m3) * z3c + m3 + 0.9e1 * zetac * zetac * m3 - 0.3e1 * m2 - 0.18e2 * zetac * zetac + m3 * m2 - 0.2e1 * zeta * zetac * m3 + 0.9e1 * zeta * pow(zetac, 0.3e1) + 0.6e1 * zeta * zetac + 0.9e1 * zetac * zetac * m2) * z3 + (0.6e1 * zeta * zetac + m3 * m2 - 0.3e1 * m3 - 0.2e1 * zeta * zetac * m2 + m2) * z3c - 0.9e1 * zeta * zetac * zetac - 0.4e1 * zetac * m2 - zeta + 0.3e1 * zetac + 0.3e1 * zeta * zetac * zetac * m2) * z2 + (-zetac * zetac * z3c + pow(zetac, 0.3e1)) * pow(z3, 0.3e1) + (-zetac * (0.3e1 * zeta * zetac - 0.3e1 + m3) * z3c + 0.3e1 * zetac * zetac * (zeta * zetac - 0.2e1 + 0.2e1 * m3)) * z3 * z3 + ((0.6e1 * zeta * zetac - 0.2e1 * zeta * zetac * m3 + m3 * m3 - 0.2e1 * m3) * z3c + zetac * (0.3e1 * zeta * zetac * m3 - 0.9e1 * zeta * zetac - 0.4e1 * m3 + 0.3e1)) * z3 + 0.3e1 * zeta * zetac - zeta * z3c;
  coeffs[6] = ((-0.3e1 * zetac * z3c + 0.3e1 * zetac * zetac) * z3 + (-zeta * zetac + m3) * z3c + zetac * (zeta * zetac + m2 - 0.1e1)) * pow(z2, 0.4e1) + ((-0.9e1 * zetac * z3c + 0.9e1 * zetac * zetac) * z3 * z3 + ((0.6e1 * m3 - 0.9e1 * zeta * zetac + 0.3e1 * zetac * zetac + 0.3e1 * m2) * z3c + 0.3e1 * zetac * (-zetac * zetac + 0.3e1 * zeta * zetac + 0.2e1 * m3 - 0.3e1 + m2)) * z3 + (0.3e1 * zeta + 0.2e1 * zetac * m2 - zetac * m3 + zeta * zetac * zetac - 0.2e1 * zeta * m2 - zetac) * z3c - m2 * m2 - 0.6e1 * zeta * zetac - zeta * pow(zetac, 0.3e1) + 0.4e1 * zeta * zetac * m2 - 0.3e1 * zetac * zetac * m2 + m2 + 0.2e1 * zetac * zetac) * pow(z2, 0.3e1) + ((-0.3e1 * zetac * z3c + 0.3e1 * zetac * zetac) * pow(z3, 0.3e1) + ((0.3e1 * m3 - 0.9e1 * zeta * zetac + 0.6e1 * m2 + 0.9e1 * zetac * zetac) * z3c - 0.3e1 * zetac * (0.3e1 * zetac * zetac - 0.3e1 * zeta * zetac + m2 + 0.3e1 - 0.4e1 * m3)) * z3 * z3 + ((0.9e1 * zeta * zetac * zetac - 0.9e1 * zetac - 0.3e1 * zeta * m2 + 0.9e1 * zeta - 0.3e1 * zetac * m3 + 0.12e2 * zetac * m2 - 0.3e1 * zeta * m3) * z3c - 0.3e1 * m3 - 0.9e1 * zeta * pow(zetac, 0.3e1) - m3 * m2 + 0.6e1 * m2 + 0.18e2 * zetac * zetac + 0.6e1 * zeta * zetac * m2 + 0.6e1 * zeta * zetac * m3 - 0.9e1 * zetac * zetac * m3 - 0.3e1 * m2 * m2 - 0.18e2 * zetac * zetac * m2 - 0.18e2 * zeta * zetac) * z3 + (-0.6e1 * zeta * zetac - 0.2e1 * m2 + 0.3e1 * m3 - 0.2e1 * m3 * m2 + m2 * m2 + 0.4e1 * zeta * zetac * m2) * z3c + 0.8e1 * zetac * m2 - 0.3e1 * zetac + 0.3e1 * zeta + 0.9e1 * zeta * zetac * zetac - 0.6e1 * zeta * zetac * zetac * m2 - 0.2e1 * zeta * m2 - 0.3e1 * zetac * m2 * m2) * z2 * z2 + (((-zeta * zetac + 0.3e1 * zetac * zetac + m2) * z3c - zetac * (0.3e1 * zetac * zetac - zeta * zetac + 0.1e1 - 0.2e1 * m3 + m2)) * pow(z3, 0.3e1) + ((0.9e1 * zeta * zetac * zetac + 0.3e1 * zetac * m3 + 0.3e1 * zeta + 0.6e1 * zetac * m2 - 0.9e1 * zetac - 0.2e1 * zeta * m3) * z3c - 0.9e1 * zetac * zetac * m2 + 0.18e2 * zetac * zetac - 0.18e2 * zetac * zetac * m3 - 0.2e1 * m3 - 0.6e1 * zeta * zetac + m3 * m3 + 0.3e1 * m2 + 0.4e1 * zeta * zetac * m3 - 0.2e1 * m3 * m2 - 0.9e1 * zeta * pow(zetac, 0.3e1)) * z3 * z3 + ((0.6e1 * zeta * zetac * m2 + 0.6e1 * m3 - 0.18e2 * zeta * zetac - 0.3e1 * m2 - 0.3e1 * m3 * m3 + 0.6e1 * zeta * zetac * m3 - m3 * m2) * z3c + 0.27e2 * zeta * zetac * zetac - 0.2e1 * zeta * m3 - 0.9e1 * zeta * zetac * zetac * m3 - 0.6e1 * zetac * m3 * m2 - 0.9e1 * zeta * zetac * zetac * m2 + 0.3e1 * zeta + 0.12e2 * zetac * m2 + 0.12e2 * zetac * m3 - 0.9e1 * zetac) * z3 - zeta * (-0.3e1 + 0.2e1 * m2) * z3c + 0.6e1 * zeta * zetac * m2 - m2 - 0.9e1 * zeta * zetac) * z2 + (zetac * (zeta * zetac - 0.1e1 + m3) * z3c - zetac * zetac * (zeta * zetac + 0.3e1 * m3 - 0.2e1)) * pow(z3, 0.3e1) + ((-0.6e1 * zeta * zetac + 0.4e1 * zeta * zetac * m3 + m3 - m3 * m3) * z3c - zetac * (-0.9e1 * zeta * zetac + 0.6e1 * zeta * zetac * m3 + 0.3e1 * m3 * m3 + 0.3e1 - 0.8e1 * m3)) * z3 * z3 + (-zeta * (-0.3e1 + 0.2e1 * m3) * z3c - 0.9e1 * zeta * zetac + 0.6e1 * zeta * zetac * m3 - m3) * z3 + zeta;
  coeffs[5] = ((-0.3e1 * zetac * zetac + 0.3e1 * zetac * z3c) * z3 * z3 + ((0.3e1 * zeta * zetac - 0.2e1 * m3) * z3c - zetac * (0.3e1 * zeta * zetac + 0.2e1 * m3 + 0.3e1 * m2 - 0.3e1)) * z3 + zeta * (m2 - 0.1e1) * z3c - 0.2e1 * zeta * zetac * (m2 - 0.1e1)) * pow(z2, 0.4e1) + ((-0.3e1 * zetac * zetac + 0.3e1 * zetac * z3c) * pow(z3, 0.3e1) + ((0.9e1 * zeta * zetac - 0.3e1 * zetac * zetac - 0.3e1 * m3 - 0.3e1 * m2) * z3c - 0.3e1 * zetac * (-zetac * zetac + 0.3e1 * zeta * zetac - 0.3e1 + 0.4e1 * m3 + m2)) * z3 * z3 + ((-0.3e1 * zeta * zetac * zetac + 0.3e1 * zetac + 0.6e1 * zeta * m2 + zetac * m3 - 0.6e1 * zetac * m2 + 0.3e1 * zeta * m3 - 0.9e1 * zeta) * z3c + 0.3e1 * zeta * pow(zetac, 0.3e1) - 0.12e2 * zeta * zetac * m2 - 0.3e1 * m2 - 0.6e1 * zetac * zetac - m3 * m2 + 0.18e2 * zeta * zetac + 0.3e1 * zetac * zetac * m3 + 0.3e1 * m2 * m2 + 0.3e1 * m3 + 0.9e1 * zetac * zetac * m2 - 0.6e1 * zeta * zetac * m3) * z3 - (m2 - 0.1e1) * (m2 - m3 + 0.2e1 * zeta * zetac) * z3c - (m2 - 0.1e1) * (zeta * m2 - 0.3e1 * zetac * m2 - 0.3e1 * zeta * zetac * zetac + zetac - 0.3e1 * zeta)) * pow(z2, 0.3e1) + (((0.3e1 * zeta * zetac - 0.2e1 * m2 - 0.3e1 * zetac * zetac) * z3c + zetac * (0.3e1 * zetac * zetac - 0.3e1 * zeta * zetac + 0.3e1 - 0.6e1 * m3 + m2)) * pow(z3, 0.3e1) + ((-0.12e2 * zetac * m2 - 0.9e1 * zeta + 0.9e1 * zetac - 0.3e1 * zetac * m3 + 0.3e1 * zeta * m2 - 0.9e1 * zeta * zetac * zetac + 0.6e1 * zeta * m3) * z3c + 0.6e1 * m3 + 0.9e1 * zeta * pow(zetac, 0.3e1) - 0.6e1 * zeta * zetac * m2 - 0.6e1 * m2 + 0.18e2 * zeta * zetac + 0.2e1 * m3 * m2 - 0.3e1 * m3 * m3 + 0.18e2 * zetac * zetac * m3 - 0.12e2 * zeta * zetac * m3 - 0.18e2 * zetac * zetac + 0.18e2 * zetac * zetac * m2 + 0.3e1 * m2 * m2) * z3 * z3 + ((0.6e1 * m2 - 0.3e1 * m2 * m2 + 0.18e2 * zeta * zetac + 0.3e1 * m3 * m3 - 0.6e1 * zeta * zetac * m3 - 0.12e2 * zeta * zetac * m2 - 0.6e1 * m3 + 0.2e1 * m3 * m2) * z3c + 0.6e1 * zeta * m2 + 0.9e1 * zetac + 0.12e2 * zetac * m3 * m2 - 0.27e2 * zeta * zetac * zetac - 0.24e2 * zetac * m2 - 0.2e1 * zeta * m2 * m3 - 0.12e2 * zetac * m3 + 0.6e1 * zeta * m3 - 0.9e1 * zeta + 0.9e1 * zeta * zetac * zetac * m3 + 0.18e2 * zeta * zetac * zetac * m2 + 0.9e1 * zetac * m2 * m2) * z3 - zeta * (m2 - 0.1e1) * (m2 - 0.3e1) * z3c + (m2 - 0.1e1) * (0.3e1 * zeta * zetac * m2 - 0.2e1 * m2 - 0.9e1 * zeta * zetac)) * z2 * z2 + (((-zeta - 0.2e1 * zetac * m2 + 0.3e1 * zetac - 0.3e1 * zetac * m3 - 0.3e1 * zeta * zetac * zetac + zeta * m3) * z3c + 0.9e1 * zetac * zetac * m3 - m2 + 0.3e1 * zeta * pow(zetac, 0.3e1) + m3 + 0.2e1 * zeta * zetac - 0.2e1 * zeta * zetac * m3 + m3 * m2 + 0.3e1 * zetac * zetac * m2 - m3 * m3 - 0.6e1 * zetac * zetac) * pow(z3, 0.3e1) + ((0.3e1 * m3 * m3 + 0.3e1 * m2 + 0.18e2 * zeta * zetac - 0.6e1 * zeta * zetac * m2 - 0.3e1 * m3 - 0.12e2 * zeta * zetac * m3 - m3 * m2) * z3c + 0.9e1 * zetac * m3 * m3 - 0.27e2 * zeta * zetac * zetac + 0.18e2 * zeta * zetac * zetac * m3 + 0.12e2 * zetac * m3 * m2 - zeta * m3 * m3 - 0.3e1 * zeta - 0.12e2 * zetac * m2 + 0.9e1 * zeta * zetac * zetac * m2 + 0.9e1 * zetac - 0.24e2 * zetac * m3 + 0.4e1 * zeta * m3) * z3 * z3 + (-zeta * (0.2e1 * m3 * m2 + 0.9e1 - 0.6e1 * m3 - 0.6e1 * m2) * z3c - 0.18e2 * zeta * zetac * m2 + 0.3e1 * m2 + 0.3e1 * m3 + 0.6e1 * zeta * zetac * m2 * m3 - 0.18e2 * zeta * zetac * m3 - 0.4e1 * m3 * m2 + 0.27e2 * zeta * zetac) * z3 + 0.3e1 * zeta * (m2 - 0.1e1)) * z2 + (-0.2e1 * zeta * zetac * (-0.1e1 + m3) * z3c + (-0.1e1 + m3) * zetac * (0.3e1 * m3 + 0.3e1 * zeta * zetac - 0.1e1)) * pow(z3, 0.3e1) + (-zeta * (-0.1e1 + m3) * (-0.3e1 + m3) * z3c + (-0.1e1 + m3) * (0.3e1 * zeta * zetac * m3 - 0.2e1 * m3 - 0.9e1 * zeta * zetac)) * z3 * z3 + 0.3e1 * zeta * (-0.1e1 + m3) * z3;
  coeffs[4] = ((-zetac * z3c + zetac * zetac) * pow(z3, 0.3e1) + ((-0.3e1 * zeta * zetac + m3) * z3c + zetac * (0.3e1 * zeta * zetac - 0.3e1 + 0.3e1 * m2 + 0.4e1 * m3)) * z3 * z3 + (-zeta * (m3 - 0.3e1 + 0.3e1 * m2) * z3c + 0.6e1 * zeta * zetac * m2 + 0.2e1 * zeta * zetac * m3 + m3 * m2 - m3 - 0.6e1 * zeta * zetac) * z3 + zeta * pow(m2 - 0.1e1, 0.2e1)) * pow(z2, 0.4e1) + (((-0.3e1 * zeta * zetac + zetac * zetac + m2) * z3c + zetac * (-zetac * zetac + 0.3e1 * zeta * zetac - 0.3e1 + 0.6e1 * m3 + m2)) * pow(z3, 0.3e1) + ((0.9e1 * zeta - 0.6e1 * zeta * m2 + zetac * m3 + 0.3e1 * zeta * zetac * zetac + 0.6e1 * zetac * m2 - 0.6e1 * zeta * m3 - 0.3e1 * zetac) * z3c - 0.6e1 * m3 + 0.12e2 * zeta * zetac * m2 + 0.3e1 * m2 + 0.2e1 * m3 * m2 - 0.3e1 * zeta * pow(zetac, 0.3e1) + 0.3e1 * m3 * m3 - 0.6e1 * zetac * zetac * m3 - 0.18e2 * zeta * zetac + 0.12e2 * zeta * zetac * m3 + 0.6e1 * zetac * zetac - 0.9e1 * zetac * zetac * m2 - 0.3e1 * m2 * m2) * z3 * z3 + ((0.3e1 * m2 * m2 - 0.6e1 * zeta * zetac - 0.3e1 * m2 + 0.2e1 * zeta * zetac * m3 - m3 * m3 + 0.2e1 * m3 - m3 * m2 + 0.6e1 * zeta * zetac * m2) * z3c - 0.12e2 * zeta * m2 - 0.3e1 * zetac - 0.6e1 * zetac * m3 * m2 + 0.9e1 * zeta * zetac * zetac + 0.9e1 * zeta + 0.3e1 * zeta * m2 * m2 + 0.4e1 * zetac * m3 - 0.6e1 * zeta * m3 + 0.12e2 * zetac * m2 + 0.4e1 * zeta * m2 * m3 - 0.3e1 * zeta * zetac * zetac * m3 - 0.9e1 * zeta * zetac * zetac * m2 - 0.9e1 * zetac * m2 * m2) * z3 + zeta * pow(m2 - 0.1e1, 0.2e1) * z3c - pow(m2 - 0.1e1, 0.2e1) * (m2 + 0.3e1 * zeta * zetac)) * pow(z2, 0.3e1) + (((-0.3e1 * zetac - 0.3e1 * zeta * m3 - zeta * m2 + 0.4e1 * zetac * m2 + 0.3e1 * zetac * m3 + 0.3e1 * zeta * zetac * zetac + 0.3e1 * zeta) * z3c - 0.3e1 * m3 + 0.3e1 * m3 * m3 - 0.3e1 * zeta * pow(zetac, 0.3e1) - 0.6e1 * zeta * zetac + 0.2e1 * m2 + 0.2e1 * zeta * zetac * m2 - 0.6e1 * zetac * zetac * m2 + 0.6e1 * zetac * zetac - 0.9e1 * zetac * zetac * m3 - m2 * m2 - m3 * m2 + 0.6e1 * zeta * zetac * m3) * pow(z3, 0.3e1) + ((0.12e2 * zeta * zetac * m3 + 0.12e2 * zeta * zetac * m2 - 0.18e2 * zeta * zetac + 0.3e1 * m3 + 0.3e1 * m2 * m2 + 0.2e1 * m3 * m2 - 0.3e1 * m3 * m3 - 0.6e1 * m2) * z3c + 0.9e1 * zeta - 0.24e2 * zetac * m3 * m2 - 0.6e1 * zeta * m2 + 0.24e2 * zetac * m3 - 0.9e1 * zetac * m3 * m3 + 0.27e2 * zeta * zetac * zetac + 0.3e1 * zeta * m3 * m3 + 0.4e1 * zeta * m2 * m3 - 0.12e2 * zeta * m3 - 0.9e1 * zetac - 0.18e2 * zeta * zetac * zetac * m3 + 0.24e2 * zetac * m2 - 0.9e1 * zetac * m2 * m2 - 0.18e2 * zeta * zetac * zetac * m2) * z3 * z3 + (zeta * (-0.6e1 * m3 - 0.12e2 * m2 + 0.4e1 * m3 * m2 + 0.9e1 + 0.3e1 * m2 * m2) * z3c - 0.3e1 * m3 + 0.18e2 * zeta * zetac * m3 - 0.12e2 * zeta * zetac * m2 * m3 + 0.6e1 * m2 * m2 - 0.3e1 * m2 * m2 * m3 + 0.36e2 * zeta * zetac * m2 - 0.27e2 * zeta * zetac + 0.8e1 * m3 * m2 - 0.9e1 * zeta * zetac * m2 * m2 - 0.6e1 * m2) * z3 + 0.3e1 * zeta * pow(m2 - 0.1e1, 0.2e1)) * z2 * z2 + (((m3 * m2 - 0.6e1 * zeta * zetac + 0.6e1 * zeta * zetac * m3 + 0.2e1 * zeta * zetac * m2 - m2) * z3c - 0.3e1 * zeta * zetac * zetac * m2 + 0.4e1 * zetac * m2 - 0.9e1 * zetac * m3 * m3 - 0.9e1 * zeta * zetac * zetac * m3 - 0.6e1 * zetac * m3 * m2 + zeta * m3 * m3 + 0.9e1 * zeta * zetac * zetac - 0.3e1 * zetac + 0.12e2 * zetac * m3 + zeta - 0.2e1 * zeta * m3) * pow(z3, 0.3e1) + (zeta * (-0.12e2 * m3 + 0.3e1 * m3 * m3 + 0.9e1 + 0.4e1 * m3 * m2 - 0.6e1 * m2) * z3c - 0.9e1 * zeta * zetac * m3 * m3 + 0.6e1 * m3 * m3 - 0.27e2 * zeta * zetac - 0.12e2 * zeta * zetac * m2 * m3 - 0.3e1 * m2 + 0.18e2 * zeta * zetac * m2 + 0.36e2 * zeta * zetac * m3 - 0.6e1 * m3 - 0.3e1 * m3 * m3 * m2 + 0.8e1 * m3 * m2) * z3 * z3 + 0.3e1 * zeta * (0.3e1 - 0.3e1 * m2 + 0.2e1 * m3 * m2 - 0.3e1 * m3) * z3) * z2 + (zeta * pow(-0.1e1 + m3, 0.2e1) * z3c - pow(-0.1e1 + m3, 0.2e1) * (m3 + 0.3e1 * zeta * zetac)) * pow(z3, 0.3e1) + 0.3e1 * zeta * pow(-0.1e1 + m3, 0.2e1) * z3 * z3;
  coeffs[3] = ((zeta * zetac * z3c - zetac * (zeta * zetac + m2 + 0.2e1 * m3 - 0.1e1)) * pow(z3, 0.3e1) + (zeta * (0.2e1 * m3 - 0.3e1 + 0.3e1 * m2) * z3c - m3 * m3 + 0.6e1 * zeta * zetac - 0.4e1 * zeta * zetac * m3 - 0.6e1 * zeta * zetac * m2 + 0.2e1 * m3 - 0.2e1 * m3 * m2) * z3 * z3 - zeta * (m2 - 0.1e1) * (0.2e1 * m3 - 0.3e1 + 0.3e1 * m2) * z3) * pow(z2, 0.4e1) + (((0.2e1 * zeta * m2 + zetac - zeta * zetac * zetac - 0.2e1 * zetac * m2 - zetac * m3 - 0.3e1 * zeta + 0.3e1 * zeta * m3) * z3c - m3 * m2 + 0.3e1 * m3 - 0.4e1 * zeta * zetac * m2 - m2 - 0.3e1 * m3 * m3 + 0.6e1 * zeta * zetac - 0.6e1 * zeta * zetac * m3 - 0.2e1 * zetac * zetac + zeta * pow(zetac, 0.3e1) + 0.3e1 * zetac * zetac * m2 + m2 * m2 + 0.3e1 * zetac * zetac * m3) * pow(z3, 0.3e1) + ((0.3e1 * m2 + 0.6e1 * zeta * zetac - 0.6e1 * zeta * zetac * m2 - 0.3e1 * m2 * m2 - m3 + m3 * m3 - 0.4e1 * zeta * zetac * m3 - m3 * m2) * z3c + 0.3e1 * zetac * m3 * m3 - 0.8e1 * zetac * m3 - 0.8e1 * zeta * m2 * m3 - 0.9e1 * zeta + 0.6e1 * zeta * zetac * zetac * m3 + 0.9e1 * zetac * m2 * m2 + 0.12e2 * zeta * m2 + 0.3e1 * zetac - 0.12e2 * zetac * m2 - 0.9e1 * zeta * zetac * zetac - 0.3e1 * zeta * m2 * m2 - 0.3e1 * zeta * m3 * m3 + 0.12e2 * zetac * m3 * m2 + 0.12e2 * zeta * m3 + 0.9e1 * zeta * zetac * zetac * m2) * z3 * z3 + (-zeta * (m2 - 0.1e1) * (0.2e1 * m3 - 0.3e1 + 0.3e1 * m2) * z3c + (m2 - 0.1e1) * (0.3e1 * m2 * m2 - 0.3e1 * m2 + 0.9e1 * zeta * zetac * m2 + 0.3e1 * m3 * m2 - 0.9e1 * zeta * zetac + 0.6e1 * zeta * zetac * m3 - m3)) * z3 + zeta * pow(m2 - 0.1e1, 0.3e1)) * pow(z2, 0.3e1) + (((-m2 * m2 - 0.6e1 * zeta * zetac * m3 + 0.6e1 * zeta * zetac - 0.4e1 * zeta * zetac * m2 - 0.2e1 * m3 * m2 + 0.2e1 * m2) * z3c - 0.2e1 * zeta * m2 * m3 + 0.6e1 * zeta * zetac * zetac * m2 + 0.3e1 * zetac - 0.3e1 * zeta - 0.12e2 * zetac * m3 + 0.9e1 * zetac * m3 * m3 + 0.3e1 * zetac * m2 * m2 + 0.6e1 * zeta * m3 + 0.2e1 * zeta * m2 - 0.9e1 * zeta * zetac * zetac + 0.9e1 * zeta * zetac * zetac * m3 + 0.12e2 * zetac * m3 * m2 - 0.8e1 * zetac * m2 - 0.3e1 * zeta * m3 * m3) * pow(z3, 0.3e1) + (-zeta * (0.8e1 * m3 * m2 + 0.3e1 * m3 * m3 + 0.9e1 - 0.12e2 * m2 - 0.12e2 * m3 + 0.3e1 * m2 * m2) * z3c - 0.36e2 * zeta * zetac * m3 + 0.24e2 * zeta * zetac * m2 * m3 + 0.6e1 * m2 + 0.6e1 * m2 * m2 * m3 - 0.36e2 * zeta * zetac * m2 - 0.6e1 * m3 * m3 - 0.6e1 * m2 * m2 - 0.16e2 * m3 * m2 + 0.6e1 * m3 + 0.6e1 * m3 * m3 * m2 + 0.9e1 * zeta * zetac * m3 * m3 + 0.9e1 * zeta * zetac * m2 * m2 + 0.27e2 * zeta * zetac) * z3 * z3 + 0.3e1 * zeta * (m2 - 0.1e1) * (-0.3e1 * m2 + m3 * m2 - 0.3e1 * m3 + 0.3e1) * z3) * z2 * z2 + ((-zeta * (-0.1e1 + m3) * (0.3e1 * m3 - 0.3e1 + 0.2e1 * m2) * z3c + (-0.1e1 + m3) * (0.3e1 * m3 * m3 - 0.3e1 * m3 + 0.3e1 * m3 * m2 + 0.9e1 * zeta * zetac * m3 + 0.6e1 * zeta * zetac * m2 - m2 - 0.9e1 * zeta * zetac)) * pow(z3, 0.3e1) + 0.3e1 * zeta * (-0.1e1 + m3) * (-0.3e1 * m2 + m3 * m2 - 0.3e1 * m3 + 0.3e1) * z3 * z3) * z2 + zeta * pow(-0.1e1 + m3, 0.3e1) * pow(z3, 0.3e1);
  coeffs[2] = ((-(-0.1e1 + m2 + m3) * zeta * z3c + (-0.1e1 + m2 + m3) * (0.2e1 * zeta * zetac + m3)) * pow(z3, 0.3e1) + (-0.1e1 + m2 + m3) * zeta * (m3 - 0.3e1 + 0.3e1 * m2) * z3 * z3) * pow(z2, 0.4e1) + (((-0.1e1 + m2 + m3) * (0.2e1 * zeta * zetac + m2) * z3c + (-0.1e1 + m2 + m3) * (0.3e1 * zeta * m3 - 0.3e1 * zeta + zetac - 0.3e1 * zetac * m3 - 0.3e1 * zeta * zetac * zetac + zeta * m2 - 0.3e1 * zetac * m2)) * pow(z3, 0.3e1) + ((-0.1e1 + m2 + m3) * zeta * (m3 - 0.3e1 + 0.3e1 * m2) * z3c - (-0.1e1 + m2 + m3) * (0.9e1 * zeta * zetac * m2 + 0.3e1 * m3 * m2 + 0.3e1 * zeta * zetac * m3 + 0.3e1 * m2 * m2 - 0.9e1 * zeta * zetac - 0.2e1 * m3 - 0.3e1 * m2)) * z3 * z3 - 0.3e1 * (-0.1e1 + m2 + m3) * zeta * pow(m2 - 0.1e1, 0.2e1) * z3) * pow(z2, 0.3e1) + (((-0.1e1 + m2 + m3) * zeta * (-0.3e1 + 0.3e1 * m3 + m2) * z3c - (-0.1e1 + m2 + m3) * (0.9e1 * zeta * zetac * m3 + 0.3e1 * m3 * m3 - 0.9e1 * zeta * zetac + 0.3e1 * zeta * zetac * m2 - 0.3e1 * m3 + 0.3e1 * m3 * m2 - 0.2e1 * m2)) * pow(z3, 0.3e1) - 0.3e1 * (-0.1e1 + m2 + m3) * zeta * (0.3e1 - 0.3e1 * m2 + 0.2e1 * m3 * m2 - 0.3e1 * m3) * z3 * z3) * z2 * z2 - 0.3e1 * (-0.1e1 + m2 + m3) * zeta * pow(-0.1e1 + m3, 0.2e1) * pow(z3, 0.3e1) * z2;
  coeffs[1] = -pow(z3, 0.3e1) * pow(-0.1e1 + m2 + m3, 0.2e1) * zeta * pow(z2, 0.4e1) + ((-pow(-0.1e1 + m2 + m3, 0.2e1) * zeta * z3c + pow(-0.1e1 + m2 + m3, 0.2e1) * (m2 + m3 + 0.3e1 * zeta * zetac)) * pow(z3, 0.3e1) + 0.3e1 * pow(-0.1e1 + m2 + m3, 0.2e1) * zeta * (m2 - 0.1e1) * z3 * z3) * pow(z2, 0.3e1) + 0.3e1 * pow(-0.1e1 + m2 + m3, 0.2e1) * zeta * (-0.1e1 + m3) * pow(z3, 0.3e1) * z2 * z2;
  coeffs[0] = -zeta * pow(z2, 0.3e1) * pow(z3, 0.3e1) * pow(-0.1e1 + m2 + m3, 0.3e1);

  return coeffs;

}

vector<complex<double>> ImgPoint::getCoeffsOpt()
{
  //complex<double>         z2c = z2;
  complex<double>         z3c = conj(z3);
  complex<double>         zeta = _sourcePos;
  complex<double>         zetac= conj(zeta);
  vector<complex<double>> coeffs(11);

	complex<double> t2 = zetac * zetac;
	complex<double>	t3 = -zetac * z3c + t2;
	complex<double>	t5 = t2 * z3c;
	complex<double>	t6 = t2 * zetac;
	complex<double>	a10 = t3 * z2 + t5 - t6;
	complex<double>	t7 = -t3;
	complex<double>	t8 = 0.3e1 * t7;
	complex<double>	t9 = z2 * z2;
	complex<double>	t12 = zeta * zetac;
	complex<double>	t13 = 0.3e1 * t2;
	complex<double>	t20 = -t5 + t6;
	complex<double>	t21 = 0.3e1 * t20;
	complex<double>	a9 = t8 * t9 + (t8 * z3 + (-m2 + t12 - t13 - m3) * z3c + zetac * (t13 - t12 + m2 + 0.1e1)) * z2 + t21 * z3 + zetac * (-t12 + 0.1e1 + m3) * z3c + t2 * (t12 - 0.2e1);
	complex<double>	t28 = -t8;
	complex<double>	t29 = t9 * z2;
	complex<double>	t31 = 0.9e1 * t3;
	complex<double>	t33 = 0.3e1 * m3;
	complex<double>	t34 = 0.2e1 * m2;
	complex<double>	t35 = 0.3e1 * t12;
	complex<double>	t42 = z3 * z3;
	complex<double>	t44 = 0.2e1 * m3;
	complex<double>	t45 = 0.9e1 * t2;
	complex<double>	t46 = 0.3e1 * m2;
	complex<double>	t53 = zetac * m3;
	complex<double>	t54 = 0.3e1 * t53;
	complex<double>	t55 = zetac * m2;
	complex<double>	t56 = 0.2e1 * t55;
	complex<double>	t57 = zeta * t2;
	complex<double>	t58 = 0.3e1 * t57;
	complex<double>	t59 = 0.3e1 * zetac;
	complex<double>	t62 = 0.2e1 * t12;
	complex<double>	t63 = zeta * t6;
	complex<double>	t64 = 0.3e1 * t63;
	complex<double>	t65 = t2 * m2;
	complex<double>	t66 = 0.3e1 * t65;
	complex<double>	t67 = 0.6e1 * t2;
	complex<double>	a8 = t28 * t29 + (t31 * z3 + (t13 + t33 + t34 - t35) * z3c - zetac * (t13 - t35 + 0.3e1 + m2)) * t9 + (t28 * t42 + ((-t35 + t44 + t45 + t46) * z3c - zetac * (t45 - t35 - t44 + t46 + 0.3e1)) * z3 + (-t54 + t56 + t58 + zeta - t59) * z3c + m2 - t62 - t64 - t66 + t67) * z2 - t21 * t42 + (-zetac * (-t35 + m3 + 0.3e1) * z3c - 0.3e1 * t2 * (t12 - 0.2e1 + m3)) * z3 + (m3 - t62) * z3c + zetac * (t35 - 0.1e1);
	complex<double>	t84 = t9 * t9;
	complex<double>	t86 = -t31;
	complex<double>	t95 = 0.6e1 * m2;
	complex<double>	t96 = 0.6e1 * m3;
	complex<double>	t97 = 0.9e1 * t12;
	complex<double>	t105 = 0.3e1 * zeta;
	complex<double>	t106 = zeta * m2;
	complex<double>	t107 = 0.4e1 * t55;
	complex<double>	t110 = 0.6e1 * t65;
	complex<double>	t111 = t12 * m2;
	complex<double>	t112 = 0.2e1 * t111;
	complex<double>	t113 = m2 * m2;
	complex<double>	t114 = 0.6e1 * t12;
	complex<double>	t117 = t42 * z3;
	complex<double>	t121 = 0.4e1 * m3;
	complex<double>	t126 = 0.9e1 * t57;
	complex<double>	t127 = zeta * m3;
	complex<double>	t128 = 0.6e1 * t55;
	complex<double>	t129 = 0.9e1 * zetac;
	complex<double>	t132 = 0.9e1 * t63;
	complex<double>	t133 = m2 * m3;
	complex<double>	t134 = t2 * m3;
	complex<double>	t135 = 0.9e1 * t134;
	complex<double>	t136 = t12 * m3;
	complex<double>	t137 = 0.2e1 * t136;
	complex<double>	t138 = 0.9e1 * t65;
	complex<double>	t139 = 0.18e2 * t2;
	complex<double>	t144 = t57 * m2;
	complex<double>	t145 = 0.3e1 * t144;
	complex<double>	t157 = m3 * m3;
	complex<double>	t160 = 0.3e1 * t136;
	complex<double>	a7 = t7 * t84 + (t86 * z3 + (-t2 + t35 - t33 - m2) * z3c - zetac * (-t2 + t35 + m2 - 0.3e1)) * t29 + (t86 * t42 + ((-t95 - t96 - t45 + t97) * z3c + 0.3e1 * zetac * (t13 - t35 + 0.3e1 + m2 - t44)) * z3 + (-t58 - t105 + t54 + t59 + t106 - t107) * z3c + t110 - t34 - t112 + t64 + t113 + t114 - t67) * t9 + (t7 * t117 + ((-t45 - t46 - m3 + t35) * z3c + zetac * (t45 - t35 + 0.3e1 + t46 - t121)) * t42 + ((-t126 + t127 + t54 - t128 + t129 - t105) * z3c + t132 - t46 + t133 + t135 - t137 + m3 + t114 + t138 - t139) * z3 + (-t33 + t133 - t112 + m2 + t114) * z3c + t145 - t107 - zeta + t59 - t126) * z2 + t20 * t117 + (-zetac * (t35 - 0.3e1 + m3) * z3c + 0.3e1 * t2 * (t12 + t44 - 0.2e1)) * t42 + ((t157 + t114 - t137 - t44) * z3c + zetac * (-t97 + t160 + 0.3e1 - t121)) * z3 + t35 - zeta * z3c;
	complex<double>	t181 = 0.2e1 * t106;
	complex<double>	t184 = 0.4e1 * t111;
	complex<double>	t185 = 0.2e1 * t2;
	complex<double>	t196 = 0.12e2 * t55;
	complex<double>	t197 = 0.3e1 * t127;
	complex<double>	t198 = 0.3e1 * t106;
	complex<double>	t199 = 0.9e1 * zeta;
	complex<double>	t202 = 0.18e2 * t65;
	complex<double>	t203 = 0.18e2 * t12;
	complex<double>	t204 = 0.3e1 * t113;
	complex<double>	t205 = 0.6e1 * t136;
	complex<double>	t206 = 0.6e1 * t111;
	complex<double>	t207 = (t196 + t126 - t197 - t54 - t198 - t129 + t199) * z3c - t202 - t203 - t132 - t204 + t95 + t139 - t133 - t33 - t135 + t205 + t206;
	complex<double>	t209 = 0.2e1 * t133;
	complex<double>	t212 = 0.6e1 * t144;
	complex<double>	t213 = zetac * t113;
	complex<double>	t214 = 0.3e1 * t213;
	complex<double>	t215 = 0.8e1 * t55;
	complex<double>	t216 = t28 * t117 + ((-t97 + t33 + t95 + t45) * z3c - 0.3e1 * zetac * (t13 - t35 + 0.3e1 + m2 - t121)) * t42 + t207 * z3 + (-t209 + t33 - t34 + t184 + t113 - t114) * z3c - t212 - t214 + t105 - t181 + t215 + t126 - t59;
	complex<double>	t224 = 0.2e1 * t127;
	complex<double>	t227 = 0.4e1 * t136;
	complex<double>	t228 = 0.18e2 * t134;
	complex<double>	t229 = (t126 + t128 + t105 - t224 + t54 - t129) * z3c - t114 - t132 + t139 + t157 + t46 - t209 + t227 - t44 - t228 - t138;
	complex<double>	t231 = 0.3e1 * t157;
	complex<double>	t234 = t57 * m3;
	complex<double>	t235 = 0.9e1 * t234;
	complex<double>	t236 = 0.12e2 * t53;
	complex<double>	t237 = t55 * m3;
	complex<double>	t238 = 0.6e1 * t237;
	complex<double>	t239 = 0.9e1 * t144;
	complex<double>	t240 = 0.27e2 * t57;
	complex<double>	a6 = (t28 * z3 + (m3 - t12) * z3c + zetac * (t12 - 0.1e1 + m2)) * t84 + (t31 * t42 + ((-t97 + t13 + t96 + t46) * z3c + 0.3e1 * zetac * (-t2 + t35 - 0.3e1 + m2 + t44)) * z3 + (t57 - t53 - t181 + t56 - zetac + t105) * z3c - t63 - t114 - t113 + m2 - t66 + t184 + t185) * t29 + t216 * t9 + (((t13 - t12 + m2) * z3c - zetac * (t13 - t12 + m2 - t44 + 0.1e1)) * t117 + t229 * t42 + ((-t203 - t231 + t206 + t205 - t133 + t96 - t46) * z3c - t235 - t224 + t236 - t238 + t105 - t239 + t240 - t129 + t196) * z3 - zeta * (-0.3e1 + t34) * z3c + t206 - m2 - t97) * z2 + (zetac * (t12 - 0.1e1 + m3) * z3c - t2 * (t12 - 0.2e1 + t33)) * t117 + ((m3 + t227 - t157 - t114) * z3c - zetac * (-t97 + t205 + t231 + 0.3e1 - 0.8e1 * m3)) * t42 + (-zeta * (t44 - 0.3e1) * z3c - t97 + t205 - m3) * z3 + zeta;
	complex<double>	t274 = -0.1e1 + m2;
	complex<double>	t275 = zeta * t274;
	complex<double>	t289 = 0.6e1 * t106;
	complex<double>	t292 = 0.12e2 * t111;
	complex<double>	t293 = 0.3e1 * t134;
	complex<double>	t294 = (-t128 - t58 + t197 + t53 + t289 + t59 - t199) * z3c - t292 - t46 + t204 + t138 + t64 - t67 + t293 - t205 - t133 + t33 + t203;
	complex<double>	t299 = 0.3e1 * t55;
	complex<double>	t310 = 0.6e1 * t127;
	complex<double>	t313 = 0.12e2 * t136;
	complex<double>	t314 = (t129 - t126 - t199 + t198 + t310 - t54 - t196) * z3c + t202 + t204 + t203 - t206 - t231 - t139 + t132 + t96 - t313 + t228 + t209 - t95;
	complex<double>	t318 = 0.24e2 * t55;
	complex<double>	t319 = 0.9e1 * t213;
	complex<double>	t320 = 0.18e2 * t144;
	complex<double>	t321 = t127 * m2;
	complex<double>	t322 = 0.2e1 * t321;
	complex<double>	t323 = 0.12e2 * t237;
	complex<double>	t324 = (t95 - t205 - t96 + t209 + t231 - t292 + t203 - t204) * z3c - t318 + t319 - t240 + t289 + t320 + t129 - t199 + t235 - t322 + t323 + t310 - t236;
	complex<double>	t329 = 0.3e1 * t111;
	complex<double>	t336 = (-t54 + t127 - t58 + t59 - t56 - zeta) * z3c + t133 + t135 - t137 + m3 + t64 + t62 + t66 - m2 - t67 - t157;
	complex<double>	t340 = zeta * t157;
	complex<double>	t341 = zetac * t157;
	complex<double>	t342 = 0.9e1 * t341;
	complex<double>	t344 = 0.18e2 * t234;
	complex<double>	t345 = 0.24e2 * t53;
	complex<double>	t346 = (t203 + t46 + t231 - t133 - t33 - t313 - t206) * z3c - t340 + t342 + 0.4e1 * t127 + t344 + t323 - t345 - t105 + t129 - t240 - t196 + t239;
	complex<double>	t351 = 0.18e2 * t111;
	complex<double>	t352 = t12 * t133;
	complex<double>	t354 = 0.18e2 * t136;
	complex<double>	t355 = 0.4e1 * t133;
	complex<double>	t356 = 0.27e2 * t12;
	complex<double>	t362 = -0.1e1 + m3;
	complex<double>	t371 = zeta * t362;
	complex<double>	a5 = (t8 * t42 + ((-t44 + t35) * z3c - zetac * (t35 + t46 - 0.3e1 + t44)) * z3 + t275 * z3c - 0.2e1 * t12 * t274) * t84 + (t8 * t117 + ((-t46 - t13 + t97 - t33) * z3c - 0.3e1 * zetac * (-t2 + t35 + t121 - 0.3e1 + m2)) * t42 + t294 * z3 - t274 * (m2 - m3 + t62) * z3c - t274 * (-t299 + t106 - t58 - t105 + zetac)) * t29 + (((-t34 + t35 - t13) * z3c + zetac * (t13 - t35 - t96 + 0.3e1 + m2)) * t117 + t314 * t42 + t324 * z3 - t275 * (-0.3e1 + m2) * z3c + t274 * (-t34 + t329 - t97)) * t9 + (t336 * t117 + t346 * t42 + (-zeta * (-t96 + t209 + 0.9e1 - t95) * z3c + t46 - t351 + 0.6e1 * t352 + t33 - t354 - t355 + t356) * z3 + 0.3e1 * t275) * z2 + (-0.2e1 * t12 * t362 * z3c + t362 * zetac * (t33 + t35 - 0.1e1)) * t117 + (-t371 * (-0.3e1 + m3) * z3c + t362 * (-t44 + t160 - t97)) * t42 + 0.3e1 * t371 * z3;
	complex<double>	t388 = t46 + m3 - 0.3e1;
	complex<double>	t393 = t274 * t274;
	complex<double>	t394 = zeta * t393;
	complex<double>	t406 = (t53 - t310 + t58 - t289 + t128 + t199 - t59) * z3c + t292 - t138 + t209 - t96 - 0.6e1 * t134 + t313 - t64 + t46 - t204 + t67 - t203 + t231;
	complex<double>	t410 = 0.12e2 * t106;
	complex<double>	t411 = 0.4e1 * t321;
	complex<double>	t415 = 0.3e1 * zeta * t113;
	complex<double>	t416 = (t206 + t204 - t46 - t133 + t137 + t44 - t114 - t157) * z3c - t410 + t196 + t126 - t59 - t239 - t319 + t199 + t411 - t238 + 0.4e1 * t53 - t310 - 0.3e1 * t234 + t415;
	complex<double>	t425 = (t105 + t58 - t197 + t54 - t59 - t106 + t107) * z3c + t34 - t114 + t231 + t112 + t67 - t110 + t205 - t33 - t135 - t133 - t113 - t64;
	complex<double>	t429 = 0.12e2 * t127;
	complex<double>	t431 = 0.3e1 * t340;
	complex<double>	t432 = (-t95 + t33 + t313 + t209 + t292 - t203 - t231 + t204) * z3c + t199 + t240 - t320 - t289 + t411 - t429 - t344 + t345 - 0.24e2 * t237 + t318 - t319 - t129 - t342 + t431;
	complex<double>	t434 = 0.12e2 * m2;
	complex<double>	t438 = 0.8e1 * t133;
	complex<double>	t439 = 0.12e2 * t352;
	complex<double>	t440 = t113 * m3;
	complex<double>	t442 = 0.36e2 * t111;
	complex<double>	t444 = 0.9e1 * t12 * t113;
	complex<double>	t445 = 0.6e1 * t113;
	complex<double>	t446 = zeta * (-t96 + t355 - t434 + t204 + 0.9e1) * z3c - t33 + t354 + t438 - t439 - 0.3e1 * t440 + t442 - t444 - t95 - t356 + t445;
	complex<double>	t453 = (-m2 + t205 + t133 - t114 + t112) * z3c - t145 + t340 - t342 + t126 + t236 - t224 - t235 - t238 - t59 + t107 + zeta;
	complex<double>	t455 = 0.12e2 * m3;
	complex<double>	t459 = 0.36e2 * t136;
	complex<double>	t460 = m2 * t157;
	complex<double>	t462 = 0.6e1 * t157;
	complex<double>	t464 = 0.9e1 * t12 * t157;
	complex<double>	t465 = zeta * (-t95 + 0.9e1 - t455 + t355 + t231) * z3c + t351 - t356 - t96 - t439 + t459 + t438 - 0.3e1 * t460 + t462 - t464 - t46;
	complex<double>	t467 = 0.3e1 - t46 + t209 - t33;
	complex<double>	t473 = t362 * t362;
	complex<double>	t474 = zeta * t473;
	complex<double>	a4 = (t3 * t117 + ((-t35 + m3) * z3c + zetac * (t35 + t121 + t46 - 0.3e1)) * t42 + (-zeta * t388 * z3c - t114 + t206 + t133 - m3 + t137) * z3 + t394) * t84 + (((t2 - t35 + m2) * z3c + zetac * (-t2 + t35 + t96 + m2 - 0.3e1)) * t117 + t406 * t42 + t416 * z3 + t394 * z3c - t393 * (m2 + t35)) * t29 + (t425 * t117 + t432 * t42 + t446 * z3 + 0.3e1 * t394) * t9 + (t453 * t117 + t465 * t42 + 0.3e1 * zeta * t467 * z3) * z2 + (t474 * z3c - t473 * (m3 + t35)) * t117 + 0.3e1 * t474 * t42;
	complex<double>	t487 = t44 + t46 - 0.3e1;
	complex<double>	t498 = (-t53 + t197 + t181 - t57 + zetac - t56 - t105) * z3c + t113 + t66 - m2 + t63 + t114 + t293 - t205 - t133 + t33 - t231 - t184 - t185;
	complex<double>	t506 = (-m3 - t133 - t227 - t206 + t46 + t157 - t204 + t114) * z3c - t199 - t196 - t126 + t319 + t239 - t415 - 0.8e1 * t321 + t323 + t429 + 0.6e1 * t234 - 0.8e1 * t53 - t431 + 0.3e1 * t341 + t410 + t59;
	complex<double>	t510 = 0.3e1 * t133;
	complex<double>	t511 = 0.9e1 * t111;
	complex<double>	t522 = (-t184 + t114 - t209 - t205 + t34 - t113) * z3c + t214 - t105 + t235 - t322 + t323 + t310 - t236 - t126 + t181 + t342 - t431 + t212 - t215 + t59;
	complex<double>	t531 = -zeta * (-t434 + 0.9e1 + t231 - t455 + t438 + t204) * z3c - t445 + t356 + t95 + 0.24e2 * t352 - t459 - 0.16e2 * t133 + t96 + 0.6e1 * t440 + 0.6e1 * t460 - t462 + t464 - t442 + t444;
	complex<double>	t533 = -t46 + t133 + 0.3e1 - t33;
	complex<double>	t542 = 0.9e1 * t136;
	complex<double>	a3 = ((t12 * z3c - zetac * (t12 - 0.1e1 + t44 + m2)) * t117 + (zeta * t487 * z3c + t114 - t157 - t209 - t227 + t44 - t206) * t42 - t275 * t487 * z3) * t84 + (t498 * t117 + t506 * t42 + (-t275 * t487 * z3c + t274 * (t204 + t510 - t46 + t511 - t97 + t205 - m3)) * z3 + zeta * t393 * t274) * t29 + (t522 * t117 + t531 * t42 + 0.3e1 * t275 * t533 * z3) * t9 + ((-t371 * (t33 + t34 - 0.3e1) * z3c + t362 * (t231 + t510 - t33 + t542 + t206 - m2 - t97)) * t117 + 0.3e1 * t371 * t533 * t42) * z2 + zeta * t473 * t362 * t117;
	complex<double>	t555 = -0.1e1 + m2 + m3;
	complex<double>	t556 = zeta * t555;
	complex<double>	a2 = ((-t556 * z3c + (m3 + t62) * t555) * t117 + t556 * t388 * t42) * t84 + (((m2 + t62) * t555 * z3c + t555 * (-t299 + t106 - t58 - t105 + zetac + t197 - t54)) * t117 + (t556 * t388 * z3c - t555 * (t204 + t510 - t46 + t511 - t97 + t160 - t44)) * t42 - 0.3e1 * t394 * t555 * z3) * t29 + ((zeta * (m2 - 0.3e1 + t33) * t555 * z3c - t555 * (-t34 + t510 + t329 + t231 - t33 - t97 + t542)) * t117 - 0.3e1 * t556 * t467 * t42) * t9 - 0.3e1 * t474 * t555 * t117 * z2;
	complex<double>	t601 = t555 * t555;
	complex<double>	t602 = zeta * t601;
	complex<double>	a1 = -t602 * t117 * t84 + ((-t602 * z3c + t601 * (m2 + m3 + t35)) * t117 + 0.3e1 * t275 * t601 * t42) * t29 + 0.3e1 * t371 * t601 * t117 * t9;
	complex<double>	a0 = -zeta * t29 * t117 * t601 * t555;  

  return vector<complex<double>>{a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10};  
}



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

  // In order to access the data in python, 
  // we copy them to array of complex<double>
  void copy_img(ImgPoint*         img,
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


