/*****************************************************************************
*
* Main class for obtaining and printing out of Critical Curves and Caustics.
*
* Note: Includes a wrapper readable by Python ctypes module.
*
*****************************************************************************/

#ifndef CRITICALCURVECAUSTICS_H
#define CRITICALCURVECAUSTICS_H

#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <algorithm> 

#include <lens.h> 
#include <laguerre.h>

class CriticalCurveCaustic: public Lens 
{
  public:
    
    CriticalCurveCaustic(double       a,
                         double       b,
                         double       th,
                         double       m2,
                         double       m3,
                         unsigned int cccLength
                        );

    vector<vector<complex<double>>> ccVec; 
    vector<vector<complex<double>>> caVec;

    // Critical curve and caustic caulculation
    void getCC(); 
    void getCa();

    void printCCC(std::string fileName);

  private:
    unsigned int _length = 500;  
    bool _ccAvailable = false;
    bool _caAvailable = false;
    complex<double> z2c;
    complex<double> z3c;
    vector<complex<double>> _tempRoots;
};


#endif
