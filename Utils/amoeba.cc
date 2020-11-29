#include <amoeba.h>
#include <vector>

Amoeba::Amoeba(long int size)
{
  _size = size;
  amoebae.resize(_size);
}


/* check if the querried range of pixels coincides with already known image
   if not, add
*/
// Checks if range of pixels coincides with already known image
// If not, adds another node
void Amoeba::addNode (long int nL, long int nR, long int ny)
{
  vector<XRange>::iterator it;
  for (it = amoebae[ny].begin(); it != amoebae[ny].end(); it++) {
      if ( (*it).xright == nR && (*it).xleft == nL ) return;
  }
  amoebae[ny].push_back(XRange(nL, nR));
}

// check if range is not within known image
// returns 1 if not in existing image 
bool Amoeba::checkLine(long int ny, long int nx)
{
  vector<XRange>::const_iterator it;
   
  if ((ny < 0) || (ny >= _size)) return 0;

  for (it = amoebae[ny].begin(); it != amoebae[ny].end(); it++) {
      if ( ( nx <= (*it).xright) && ((*it).xleft <= nx)) {
           return 0;
      }
  }
  return 1;       
}

// check if range is not within known image
// returns 1 if not in existing image 
bool Amoeba::checkLineShift(long int ny, long int& nx)
{
  vector<XRange>::const_iterator it;
   
  if ((ny < 0) || (ny >= _size))
  {  
    // Got out of the range
    // Make sure to get beyond the line
    nx += _size;
    return 0;
  }

  for (it = amoebae[ny].begin(); it != amoebae[ny].end(); it++) {
      if ( ( nx <= (*it).xright) && ((*it).xleft <= nx)) {
        nx = (*it).xright+2;  
        return 0;
      }
  }
  return 1;       
}


// Returns vector of intersections of input interval with amoeba on the other line 
vector<XRange> Amoeba::getUnfilled(long int  nxL,
                                   long int  nxR,
                                   long int  ny
                                  )
{
  vector<XRange>::const_iterator it;
  vector<XRange> unfilledIntervals;
  long int l = nxL;
  long int r;
  bool intervalClosed = true;

  if ((ny < 0)  || (ny >= _size) || 
      (nxL < 0) || (nxR >= _size) 
     ) return unfilledIntervals;

  // sort the ranges
  std::sort(amoebae[ny].begin(), amoebae[ny].end());

  for (it = amoebae[ny].begin(); it != amoebae[ny].end(); it++)
  {
      // Check if there is an overlap with nxL to nxR
      if ( ( l <= (*it).xright) && ((*it).xleft <= nxR))
      {
        if (l < (*it).xleft)
        {
          // xright should be one position to the right from the filled range
          unfilledIntervals.push_back(XRange(l, (*it).xleft - 2));
          intervalClosed = true;
        }
        
        // sets l for next step
        if((*it).xright + 2 > nxR)
        {
          break;
        }
        else
        {
          l = (*it).xright + 2;
          intervalClosed = false;
        }
      }
  }

  // all last one by end of the range
  if(!intervalClosed)
  {
    unfilledIntervals.push_back(XRange(l, nxR));
  }

  return unfilledIntervals;
}


void Amoeba::resize(long int size)
{
  amoebae.resize(size);
  _size = size;
}

