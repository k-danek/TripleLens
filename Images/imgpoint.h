/*****************************************************************************
*
* Class for getting a point source images.
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
#include <iomanip>

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

    ImgPoint(const LensPar &lensParam);


    vector<complex<double>> roots; // As not all of the roots are images  
    vector<complex<double>> imgs;  // Roots verified to be images
    vector<bool>            isImg; // Vectors of bools to more 

    vector<complex<double>> getCoeffs();
    vector<complex<double>> getCoeffsOpt();

    // checks if an roos is an image
    bool imgCheck(complex<double> img,
                  double sourceRadius
                 );

    // Image position caulculation
    void getRoots(bool forceNewRoots); 
    void getImages();
    
    // update and return images in one functional call
    vector<complex<double>> getImages(complex<double> pos);

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
