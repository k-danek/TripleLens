/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/***************************************************************************
* Following code is greatly inspired by version of the Laguerre's method published in book:
* Press, W.H.; Teukolsky, S.A.; Vetterling, W.T.; Flannery, B.P. (2007). "Numerical Recipes: The Art of Scientific Computing (3rd ed.)". New York: Cambridge University Press. [1]
* 
* For explanation of the basic idea behind the algorithm, wikipedia page is sufficient:
* https://en.wikipedia.org/wiki/Laguerre%27s_method
*
* Also see https://github.com/andresmmera/Laguerre-Root-Solver
* for an alternative implementation.
****************************************************************************/

#ifndef LAGUERRE_CUH
#define LAGUERRE_CUH

#include <iostream>
#include <vector>
#include <complex>
#include <math.h>
#include <unordered_map>
#include <limits>
#include <list>
#include <typeinfo>
#include <functional>

#include <thrust/complex.h>
#include <thrust/device_reference.h>


using std::vector;

template<class T>
using complex = std::complex<T>;

#define MAX_IT 300
#define MT 8

class LaguerreCUDA
{
public:
  LaguerreCUDA(const std::vector<complex<double>>& coeffs);//Class constructor

  static constexpr size_t MAX_SIZE = 11; // Maximum number of coeffitients
  thrust::complex<double> data[MAX_SIZE];

  //// Checks each root individually by substituing it into original polynomial
  //bool checkRootsOneByOne(const vector<complex<double>>& roots);

  //// Reconstructs polynomial coeffs using roots and checks agains their
  //// original values 
  //bool checkRootsAllAtOnce(const vector<complex<double>>& roots);

  //// Makes simplified check testing reconstructing zeroth and n-1 coeffitient
  //// of the polynomial using the found roots
  //bool checkRoots(const vector<complex<double>>& roots);

  constexpr static double EPS = 1.0e-10;  

private:
  size_t _size;
  const thrust::complex<double>* _polyCoeffs;// Polynomial coefficients
};

// Solver with initial estimate on start of the iteration
__device__
bool solveRootsCUDA(thrust::complex<double>*             roots,
                    const thrust::complex<double>* const coeffs,
                    const size_t                         polOrder,
                    const size_t                         maxIt);

// A simple root polisher - can in principle merge two roots into one
__device__
void polishRootsCUDA(const thrust::complex<double>*       roots,
                     const thrust::complex<double>* const coeffs); 

// The main routine setting x to a root of a polynomial
__device__
bool laguerreCUDA(thrust::complex<double>&             x,
                  const thrust::complex<double>* const coeffs,
                  const size_t                         size,
                  const size_t                         maxIt)




#endif
