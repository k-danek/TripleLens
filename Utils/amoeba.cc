#include <amoeba.h>
#include <vector>

Amoeba::Amoeba(long int size)
{
  _size = size;
}


/* check if the querried range of pixels coincides with already known image
   if not, add
*/
// Checks if range of pixels coincides with already known image
// If not, adds another node
void Amoeba::addNode (long int nL, long int nR, long int ny)
{
  list<XRange>::iterator it;
  
  // Insert a new line into the map
  if (amoebae.find(ny) == amoebae.end())
  {
    amoebae.insert({ny,{XRange(nL,nR)}});
    return; 
  }

  for (it = amoebae[ny].begin(); it != amoebae[ny].end(); it++) {
      if ( (*it).xright == nR && (*it).xleft == nL ) return;
  }
  amoebae[ny].push_back(XRange(nL, nR));
}

// check if range is not within known image
// returns 1 if not in existing image 
bool Amoeba::checkLine(long int ny, long int nx)
{
  list<XRange>::const_iterator it;
  if (amoebae.find(ny) == amoebae.end())
  {
    return 1;
  }

  if ((ny < 0) || (ny >= _size)) return 0;

  for (it = amoebae[ny].begin(); it != amoebae[ny].end(); it++)
  {
      if ( ( nx <= (*it).xright) && ((*it).xleft <= nx))
      {
           return 0;
      }
  }
  return 1;       
}

// check if range is not within known image
// returns 1 if not in existing image 
bool Amoeba::checkLineShift(long int ny, long int& nx)
{
  list<XRange>::const_iterator it;

  //std::cout << "check line shift entered \n";

  if ((ny < 0) || (ny >= _size))
  {  
    // Got out of the range
    // Make sure to get beyond the line
    nx += _size;
    return 0;
  }

  if (amoebae.find(ny) == amoebae.end())
  {
    return 1;
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
  _size = size;
}

