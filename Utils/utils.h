#include<vector>
#include<complex>


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
