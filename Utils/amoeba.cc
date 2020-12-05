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

void Amoeba::resize(long int size)
{
  amoebae.resize(size);
  _size = size;
}

