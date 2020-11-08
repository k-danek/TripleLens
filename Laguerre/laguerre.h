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
****************************************************************************/


#ifndef LAGUERRE_H
#define LAGUERRE_H

#include <iostream>
#include <complex>
#include <vector>
#include <limits>
#include <typeinfo>

using std::vector;

template<class T>
using complex = std::complex<T>;

#define MAX_IT 100
#define MT 8

class Laguerre
{
public:
  Laguerre(vector<complex<double>>);//Class constructor
  vector<complex<double>> solveRoots();
  vector<complex<double>> polishRoots(vector<complex<double>> roots);

  bool checkRoots(const vector<complex<double>>& roots);
  bool checkRootsOneByOne(const vector<complex<double>>& roots);
  bool checkRootsAllAtOnce(const vector<complex<double>>& roots);

  bool laguerre(vector<complex<double>> poly, complex<double> &x);
  
  constexpr static double EPS = std::numeric_limits<double>::epsilon();  

private:
  vector<complex<double>> _polyCoeffs;// Polynomial coefficients
};

#endif
