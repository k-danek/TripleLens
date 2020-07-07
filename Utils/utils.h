#ifndef UTILS_H
#define UTILS_H

#include<vector>
#include<complex>


template<class T>
using complex = std::complex<T>;

// Template over complex nubers
// to give bottom-left-most (minpoint)
// and top-right-most (maxpoint) points
// of enveloping rectangle for points 
// in the container
template<class T>
void minmax(const T          &container,
            complex<double>  &min,
            complex<double>  &max
           ) 
{
  for(auto el: container)
  {
    if(el.real() < min.real()) min.real(el.real());
    if(el.imag() < min.imag()) min.imag(el.imag());
    if(el.real() > max.real()) max.real(el.real());
    if(el.imag() > max.imag()) max.imag(el.imag());
  }
}

static void makeSquare(complex<double>  &min,
                       complex<double>  &max
                      ) 
{
  complex<double> centre = (max+min)/2.0;

  double halfEdge = max.real() - min.real();
  if(max.imag()-min.imag() > halfEdge)
    halfEdge = max.imag() - min.imag();

  max = centre + 1.05*complex<double>{halfEdge, halfEdge};
  min = centre - 1.05*complex<double>{halfEdge, halfEdge};
}

#endif
