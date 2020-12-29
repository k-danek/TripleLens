//#ifndef AMOEBA_H
//#define AMOEBA_H

#include <vector>
#include <list>
#include <algorithm>
#include <unordered_map>
#include <utility>
#include <iostream>

using std::vector;
using std::list;

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

  // useful for sorting non-overlapping ranges
  bool operator<(const XRange& xRange) const
  {
      return (xleft < xRange.xleft);
  }

};

/* main data structure - amoebas as a vector of lists of XRanges */

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

    // checks the line but with each already filled 
    // segment jumps at its end
    bool checkLineShift(long int  ny,
                        long int& n
                       );    

    void resize(long int size
               );

    // Image plane grid filled with amoebae
    //vector< list<XRange> > amoebae;   

    // Unordered map is substantially more effective for extreme grid sizes
    std::unordered_map<long int, list<XRange>> amoebae;

  private:
    long int _size;
};

//#endif  
