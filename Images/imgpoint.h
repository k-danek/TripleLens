/*****************************************************************************
*
* Class for getting a point source images and amplification.
*
* Note: Includes a wrapper readable by Python ctypes module.
*
*****************************************************************************/

#ifndef IMGPOINT_H
#define IMGPOINT_H

#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <algorithm> 

#include <lens.h> 
#include <laguerre.h>

class ImgPoint: public Lens 
{
  public:
    
    ImgPoint(double       a,
             double       b,
             double       th,
             double       m2,
             double       m3,
             double       posX,
             double       posY
             );

    vector<complex<double>> roots; // As not all of the roots are images  
    vector<complex<double>> imgs;  // Roots verified to be images
    vector<complex<bool>>   isImg; // Vectors of bools to more 

    vector<complex<double>> getCoeffs();
    vector<complex<double>> getCoeffsOpt();

    // checks if an roos is an image
    bool imgCheck(complex<double> img,
                  double sourceRadius
                 );

    // Critical curve and caustic caulculation
    void getRoots(); 
    void getImages();

    // update position of the source
    void setPos(double       posX,
                double       posY
               );

    // update position of the source
    void setPos(complex<double> pos);

  private:
    bool _rootsAvailable = false;
    bool _imgsAvailable = false;
    complex<double> _sourcePos;
    vector<complex<double>> _tempRoots;
};

#endif
