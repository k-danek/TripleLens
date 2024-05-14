/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 **************************************************************************/

#include "laguerre.cuh"

const int MAX_POL_ORDER = 10; // Example maximum polynomial order

LaguerreCUDA::LaguerreCUDA(const std::vector<complex<double>>& coeffs) :
  _size(coeffs.size())
{
  // Check if vector size exceeds maximum size
  if (_size > MAX_COEFF_SIZE) {
      throw std::runtime_error("Vector size exceeds maximum size");
  }
  
  // Initialize data array using elements of vec
  //std::copy(coeffs.begin(), coeffs.end(), _polyCoeffs);

  // Copy and convert each element
  for (size_t i = 0; i < coeffs.size(); ++i) {
      _polyCoeffs[i] = thrust::complex<double>(coeffs[i].real(), coeffs[i].imag());
  }

}


__device__
bool laguerreCUDA(thrust::complex<double>&             x,
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
      dx = (1.0+thrust::abs(x))*thrust::complex<double>(cos(1.0*i), sin(1.0*i));
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



__device__
bool solveRootsCUDA(thrust::complex<double>*             roots,
                    const thrust::complex<double>* const coeffs,
                    const size_t                         polOrder,
                    const size_t                         maxIt)
{
  const size_t coeffSize(polOrder+1);
  
  // indicates whether all the laguerre call have succeeded
  bool success = true;
  
  thrust::complex<double> x, b, c;
  
  // BEWARE: missing check coeffSize <= MAX_POL_ORDER
  thrust::complex<double> tempPoly[MAX_POL_ORDER+1];
  for(uint32_t i=0; i <= coeffSize; i++)
  {
    tempPoly[i] = coeffs[i];
  } 
  
  //thrust::complex<double> roots[11];
  //int size = _polyCoeffs.size();

  for (int j = polOrder; j >= 1; j--)
  {
    x = thrust::complex<double>(1.0e-8,1.0e-8);
    
    success = success && laguerreCUDA(x, tempPoly, polOrder, maxIt);
    roots[polOrder-j] = x; 

    b = tempPoly[j];

    // Deflating the polynomial by removing the root.
    for (int jj = j-1; jj >= 0; jj--)
    {
      c = tempPoly[jj];
      tempPoly[jj] = b;
      b = x*b + c;
    }
    
    // Zero-ing the highest order coefficiend from the polynomial.
    tempPoly[j] = thrust::complex<double>(0.0, 0.0);

  }

  // Polishing omitted as it is done with different function

  return success;
}

__device__
void polishRootsCUDA(thrust::complex<double>*            roots,
                     const thrust::complex<double>* const coeffs,
                     const size_t                         polOrder,
                     const size_t                         maxIt)
{
  bool success = true; 
  thrust::complex<double> tempRoot;

  for (uint32_t i = 0; i < polOrder; i++)
  {
    tempRoot = roots[i];
    success = laguerreCUDA(tempRoot, coeffs, polOrder, maxIt);
    roots[i] = tempRoot*double(success)+roots[i]*double(!success); 
  };
}

