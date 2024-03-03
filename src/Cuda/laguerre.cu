/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 **************************************************************************/

#include "laguerre.cuh"

LaguerreCUDA::LaguerreCUDA(const std::vector<complex<double>>& coeffs) :
  _size(coeffs.size())
{
  // Check if vector size exceeds maximum size
  if (_size > MAX_SIZE) {
      throw std::runtime_error("Vector size exceeds maximum size");
  }
  
  // Initialize data array using elements of vec
  std::copy(coeffs.begin(), coeffs.end(), _polyCoeffs);
}


__device__
bool laguerre(thrust::complex<double>&             x,
              const thrust::complex<double>* const coeffs,
              const size_t                         size,
              const size_t                         maxIt)
{
  const int m = size-1;
  const double md = static_cast<double>(m); // having double version for complex library
  thrust::complex<double> dx, x1, g, gp, gm, h, b, f, s;
  double error;

  x = -thrust::complex<double>(coeffs[0]/coeffs[m]);

  for (unsigned int i = 1; i <= maxIt;i++)
  {
    b = coeffs[m];
    error = thrust::abs(b);
    s = thrust::complex<double>(0.0,0.0);
    f = thrust::complex<double>(0.0,0.0);
    for (int j=m-1; j>=0; j--)
    {
      s = x*s + f;//Second derivative
      f = x*f + b;//First derivative
      b = x*b + coeffs[j];//Polynom evaluation
      error = thrust::abs(b) + thrust::abs(x)*error;//Error term
    }
    error *= 1.0e-5;
    //x is already a root; <= is there for the case when error == 0
    //if (std::abs(b) <= error)
    //  return true;

    g = f/b;
    h = g*g-2.*s/b;
    gp = g+sqrt((md-1.)*(md*h-g*g));
    gm = g-sqrt((md-1.)*(md*h-g*g));
    // make sure gp has the higher abs. value
    if (thrust::abs(gp) < thrust::abs(gm)) gp = gm;
    // if the higher value is nonzero put it into the denominator
    if (thrust::abs(gp) > 0.)
    {
      // dx denotes our distance from our root
      dx = md/gp;
    }
    else
    {
      // In the rare case that denominator would be zero
      dx = (1.+thrust::abs(x))*complex<double>(cos(1.*i), sin(1.*i));
    }

    x1 = x - dx;

    // Repeated until it converges or exceeds the maximum number of iterations
    if ((x.real() == x1.real()) && (x.imag() == x1.imag())) 
      return true;
 
    x = x1;

  }

  //std::cout << "Laguerre: maximum number of iterations exceeded for pol. order " << poly.size()-1 << "\n";
  return false;
}



