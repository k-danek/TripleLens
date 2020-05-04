/*****************************************************************************
*
* A simple class to hold basic parameters to define a triple lens.
* Chose a class instead of struct as I might need some methods later.
* 
* It represents three lens points in vertices of a triangle.
* The first lens point is in the origin of the coordinate system.
*
* Meaning of the parameters is as follows:
* a,b - length of sides of a triangle
* th - angle in radians between them
* m2 - mass fraction at the end of edge of length a
* m3 - mass fraction at the end of edge of length b
*
* Mass fraction sum up to 1, i.e., m1=1.0-m2-m3                         
*          
*****************************************************************************/

#ifndef LENS_H
#define LENS_H

#include <complex>
#include <iostream>

using namespace std;

class Lens {
  public:
	  double a; 
	  double b;
	  double th;
	  double m2;
	  double m3;
    complex<double> z1;
    complex<double> z2;
    complex<double> z3;

	  Lens(double aa,
	       double bb,
	       double theta,
	       double mass2,
	       double mass3
	      );

    void moveToCenter(); // move lens so that z1+z2+z3=0
    void moveToZ1(); // move z1 to coordinate origin

  private:
    void _setPos();
    bool _isCentered = false;    
};

#endif
