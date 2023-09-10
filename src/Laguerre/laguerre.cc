/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 **************************************************************************/

#include "laguerre.h"

// P are the polynom coefficients
Laguerre::Laguerre(vector<complex<double> > poly)
{
  _polyCoeffs = poly;
}

bool Laguerre::laguerre(const vector<complex<double> > poly, complex<double> &x)
{
  // Fractional steps to be used after breaking the iteration cycle.
  // Supposedly, it is becoming relevant for polynomials for order 20 or higher.
  vector<double> frac={0.5, 0.23, 0.73, 0.13, 0.37, 0.61, 0.89, 1.};
  const int m = poly.size()-1;
  const double md = static_cast<double>(m); // having double version for complex library
  complex<double> dx, x1, g, gp, gm, h, b, f, s;
  double error;

  // For a linear equation just return the result.
  if (m == 1)
  {
    if(std::abs(poly[m]) != 0.0 )
    {
      x = -poly[0]/poly[m];
      return true;
    }
    else
    {
      // 0.0*x + a = 0 where abs(a) > 0 does not have a finite solution, there is no point in iterating further 
      x = complex<double>(999.0,999.0);
      return false;
    }
  }

  for (unsigned int i = 1; i <= MAX_IT;i++)
  {
    b = poly[m];
    error = std::abs(b);
    s = complex<double>(0,0);
    f = complex<double>(0,0);
    for (int j=m-1; j>=0; j--)
    {
      s = x*s + f;//Second derivative
      f = x*f + b;//First derivative
      b = x*b + poly[j];//Polynom evaluation
      error = std::abs(b) + std::abs(x)*error;//Error term
    }
    error *= EPS;
    //x is already a root; <= is there for the case when error == 0
    if (std::abs(b) <= error)
      return true;

    g = f/b;
    h = g*g-2.*s/b;
    gp = g+sqrt((md-1.)*(md*h-g*g));
    gm = g-sqrt((md-1.)*(md*h-g*g));
    // make sure gp has the higher abs. value
    if (std::abs(gp) < std::abs(gm)) gp = gm;
    // if the higher value is nonzero put it into the denominator
    if (std::abs(gp) > 0.)
    {
      // dx denotes our distance from our root
      dx = md/gp;
    }
    else
    {
      // In the rare case that denominator would be zero
      dx = (1.+std::abs(x))*complex<double>(cos(1.*i), sin(1.*i));
    }

    x1 = x - dx;

    // Repeated until it converges or exceeds the maximum number of iterations
    if ((x.real() == x1.real()) && (x.imag() == x1.imag())) 
      return true;

    if (i%MT)
    {
      x = x1;
    }
    else
    {
      // Once every MT steps shift the iteration little bit further
      x -= complex<double>(frac[i/MT],0)*dx;
    }
  }

  //std::cout << "Laguerre: maximum number of iterations exceeded for pol. order " << poly.size()-1 << "\n";
  return false;
}


// Find all the roots via Laguerre method
vector<complex<double> > Laguerre::solveRoots(
                                       const vector<complex<double>>& initRoots)
{
  vector<complex<double>> tempPoly, roots;
  complex<double> x(0,0), b, c;
  int size = _polyCoeffs.size();

  tempPoly = _polyCoeffs;
  const bool isInitialised = initRoots.size()+1 == size;
  
  // Using function assignemt of lambda fuctions in a ternary operator
  std::function<complex<double>(int i)> getInit =
    (isInitialised
     ? std::function<complex<double>(int i)>( [&](int i)
       {
         return initRoots[i];
       })
     : std::function<complex<double>(int i)>( [&](int i)
       {
         return complex<double>(EPS,EPS); 
       })  
    );

  for (int j = size-1; j >= 1; j--)
  {
    x=getInit(j);
    
    if(laguerre(tempPoly, x))
    {
      roots.push_back(x); 
    }
    else
    {
      // if it wasn't initialised, try again starting from zero
      if(isInitialised)
      { x = complex<double>(EPS,EPS);
        //if(!laguerre(tempPoly,x))
        //  std::cout << "WARNING: deflating with unconverged root!!!";
      }
      roots.push_back(x); 
    }

    b = tempPoly[j];

    // Deflating the polynomial by removing the root.
    for (int jj = j-1; jj >= 0; jj--)
    {
      c = tempPoly[jj];
      tempPoly[jj] = b;
      b = x*b + c;
    }
    
    // Removing the highest order coefficiend from the polynomial.
    tempPoly.erase(tempPoly.end()-1);

  }
  // Polish the root.
  for(unsigned int j=0; j < _polyCoeffs.size()-1; j++)
    if(!laguerre(_polyCoeffs, roots[j]))
    {
      std::cout << "WARNING: Destroyed roots while polishing!!!";
    }

  return roots;
}

// In the case we have a good estimate of the root, we can run polish only routine.
vector<complex<double>> Laguerre::polishRoots(
                                         vector<complex<double>> roots
                                             )
{
  vector<complex<double>> polishedRoots;
  for (auto root: roots)
  {
    if (!laguerre(_polyCoeffs, root))
      std::cout << "WARNING: failed to polish a root";
    polishedRoots.push_back(root);
  }

  return polishedRoots;
}

// A simple test on the roost.
// Tests whether the multiple of all the roots is close enough to zeroth coeff.
// This is mainly to avoid one root added twice by polisRoots method.
bool Laguerre::checkRoots(const vector<complex<double>>& roots)
{
  // initialise to identity
  complex<double> c0 = _polyCoeffs.back(); 
  for (auto root: roots)
    c0 *= -root;

  return (std::abs(_polyCoeffs[0]) > 1.0e-8) ? 
                                        (1.0-std::abs(c0/_polyCoeffs[0]) < 1.0e-5):
                                        (std::abs(c0-_polyCoeffs[0]) < 1.0e-10);
}

// For each root check whether it satisfies the original equation
bool Laguerre::checkRootsOneByOne(const vector<complex<double>>& roots)
{
  bool successful = true;
  complex<double> b;
  const int m = _polyCoeffs.size();                
  int failedRootsCount = 0;

  for (auto root: roots)
  {
    b = _polyCoeffs[m];
    for(int j = m-1; j >= 0; j--)
    {
      b = b*root + _polyCoeffs[j];
    }

    if(std::abs(b) > 100*EPS)
    {
      successful = false;
      failedRootsCount++;
    }      

  }

  if(!successful)
    std::cout << "WARNING: " << failedRootsCount << " roots failed the full check \n";

  return successful;
}

// Reconstruct polynomial coefficients from the roots.
// This is a way how to check whole set of the roots at once
bool Laguerre::checkRootsAllAtOnce(const vector<complex<double>>& roots)
{
  bool successful = true;
  const unsigned int maxOrder = roots.size();

  if(maxOrder+1 != _polyCoeffs.size())
  {
    std::cout<< "WARNING: Trying to check all roots at once with too little roots\n";
    return false;
  }

  // Initialize coeffs to 1 to avoid assignment if highest coeff later
  vector<complex<double>> coeff(maxOrder+1, complex<double>(1.0,0.0));

  // outter loop for reconstruction coeffs for polynomials of different order
  for(unsigned int order = 1; order <= maxOrder; order++)
  {
    // thanks to the initialisation I can omit coeff[order] = complex<double>(1.0,0.0);

    // formula for the coeffs 'in the middle'
    for(unsigned int i = order-1; i >= 1; i--)
    {
      coeff[i] = coeff[i-1]-coeff[i]*roots[order-1];   
    }
    coeff[0] *= -roots[order-1];
  }

  for(unsigned int i = 0; i <= maxOrder; i++)
  {
    double distance = std::abs(coeff[i]*_polyCoeffs[maxOrder]-_polyCoeffs[i]);

    if(distance>1e4*EPS)
    {
      successful = false;
      std::cout << "WARNING: Reconstructed coeff " << i << " is " << distance << " off original value\n";
    }

  }

  return successful;
}

