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


#ifndef LAGUERRE_H
#define LAGUERRE_H

#include <iostream>
#include <complex>
#include <vector>
#include <limits>
#include <typeinfo>
#include <functional>

using std::vector;

template<class T>
using complex = std::complex<T>;

#define MAX_IT 300
#define MT 8

class Laguerre
{
public:
  Laguerre(vector<complex<double>>);//Class constructor
 
  // Solver with initial estimate on start of the iteration
  vector<complex<double>> solveRoots(const vector<complex<double>>& roots = {}); 
  // A simple root polisher - can in principle merge two roots into one
  vector<complex<double>> polishRoots(vector<complex<double>> roots);

  // Checks each root individually by substituing it into original polynomial
  bool checkRootsOneByOne(const vector<complex<double>>& roots);

  // Reconstructs polynomial coeffs using roots and checks agains their
  // original values 
  bool checkRootsAllAtOnce(const vector<complex<double>>& roots);

  // Makes simplified check testing reconstructing zeroth and n-1 coeffitient
  // of the polynomial using the found roots
  bool checkRoots(const vector<complex<double>>& roots);

  // The main routine setting x to a root of a polynomial
  bool laguerre(const vector<complex<double>> poly, complex<double> &x);

  constexpr static double EPS = std::numeric_limits<double>::epsilon();  

private:
  vector<complex<double>> _polyCoeffs;// Polynomial coefficients
};

#endif
