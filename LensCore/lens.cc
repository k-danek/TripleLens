#include <lens.h>

// constructor
Lens::Lens(double aa,
	         double bb,
	         double theta,
	         double mass2,
	         double mass3)
{
  a = aa;
  b = bb;
  th = theta;
  m2 = mass2;
  m3 = mass3;
  _setPos();
};

// Sets lens positions in the complex plane
void Lens::_setPos()
{
  z1 = {0.0,0.0};
  z2 = {a,0.0};
  z3 = {b*cos(th), b*sin(th)};
};

// Moves lenses so that their geometrical center is in coordinate origin
void Lens::moveToCenter()
{
  complex<double> zC = (z1+z2+z3)/3.0;  
  z1 -= zC;
  z2 -= zC;
  z3 -= zC;
  _isCentered = true;
};

void Lens::moveToZ1()
{
  z2 -= z2 - z1;
  z3 -= z3 - z1; 
  z1 = {0.0,0.0};
  _isCentered = false;  
};

