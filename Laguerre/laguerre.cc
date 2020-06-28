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

complex<double> Laguerre::laguerre(vector<complex<double> > poly, complex<double> x)
{

  vector<double> frac={0,0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.};//Fractional steps to be used after breaking the iteration cycle.
  const int m = poly.size()-1.0;
  const double md = static_cast<double>(m); // having double version for complex library
  complex<double> dx, x1, g, gp, gm, h, b, d, f;
  double error;
  for (unsigned int it = 1; it <= MAX_IT; it++)
  {
    b = poly[m];
    error = abs(b);
    d = complex<double>(0,0);
    f = complex<double>(0,0);
    for ( int j=m-1; j>=0; j--)
    {
      f = x*f + d;//Second derivative
      d = x*d + b;//First derivative
      b = x*b + poly[j];//Polynom evaluation
      error = abs(b) + abs(x)*error;//Error term
    }
    error *= EPSS;
    if (abs(b) < error) return x;//x is already a root

    //
    g = d/b;
    h = g*g-2.*f/b;
    gp = g+sqrt((md-1.)*(md*h-g*g));
    gm = g-sqrt((md-1.)*(md*h-g*g));
    if (abs(gp) < abs(gm)) gp = gm;
    if (std::max(abs(gp), abs(gm)) > 0.) dx = md/gp;
    else dx = (1.+abs(x))*complex<double>(cos(1.*it), sin(1.*it));
    x1 = x - dx;
    // Reapeated until it converges or exceeds the maximum number of iterations
    if ((x.real() == x1.real()) && (x.imag() == x1.imag())) return x;
    if (it%MT) x = x1;
    else x -= complex<double>(frac[it/MT],0)*dx;
  }
  return x;
}


// Routine for the Laguerre method
vector<complex<double> > Laguerre::solveRoots()
{
  
  vector<complex<double>> tempPoly, roots;
  complex<double> x(0,0), b, c;
 
  tempPoly = _polyCoeffs;
 
  for (int j = _polyCoeffs.size()-1; j >= 1; j--)
  {
    
    x=complex<double>(0.0,0.0);
    x = laguerre(tempPoly, x);
    if (abs(x.imag()) < EPS*abs(x.real())) 
      x=complex<double>(x.real(),0.0);// If imaginary part is small enough, make it a real root
    
    roots.push_back(x);// 

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
  for(unsigned int j=0; j < _polyCoeffs.size()-1; j++) roots[j] = laguerre(_polyCoeffs, roots[j]);

  return roots;
}

// In the case we have a good estimate of the root, we can run polish only routine.
vector<complex<double>> Laguerre::polishRoots(
                                         vector<complex<double>> roots
                                                  )
{
  vector<complex<double>> polishedRoots;
  for (auto root: roots)
    polishedRoots.push_back(laguerre(_polyCoeffs, root));

  return polishedRoots;
}

// A simple test on the roost.
// Tests whether the multiple of all the roots is close enough to zeroth coeff.
// This is mainly to avoid one root added twice by polisRoots method.
bool Laguerre::checkRoots(vector<complex<double>> roots)
{
  complex<double> c0 = _polyCoeffs[6]; 
  for (auto root: roots)
    c0 *= root;

  return (1.0-abs(c0/_polyCoeffs[0]) < 1.0e-3);
}




