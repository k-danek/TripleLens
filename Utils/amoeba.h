#ifndef AMOEBA_H
#define AMOEBA_H

#include<vector>

/* xright - leftmost x coordinate occupied by the lensed image
   xleft - rightmost x coordinate occupied by the lensed image
*/
struct XRange
{
  long int xleft;
  long int xright;
  
  XRange(long int xLeftMost=0, long int xRightMost=0)
  {
    xleft = xLeftMost;
    xright= xRightMost;
  }
};

/* main data structure - amoebas as a vector of lists of XRanges */
//typedef std::vector< vector<XRange> > amoeba;

class Amoeba
{
  public:
    Amoeba(long int size);
    
    void addNode(long int nxL,
                 long int nxR,
                 long int ny
                );

    bool checkLine(long int ny,
                   long int n 
                  );

    // Image plane grid filled with amoebae
    vector< vector<XRange> > amoebae;   

  private:
    long int _size;
};

#endif  
