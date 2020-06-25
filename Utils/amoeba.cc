#include<amoeba.h>


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
  for (it = amoebae[ny].begin(); it != amoebae[ny].end(); i++) {
      if ( (*it).xright == nR && (*it).xleft == nL ) return;
  }
  amoebas[ny].push_back(amoeba(nL, nR));
}

// check if range is no within known image
// returns 1 if in not in existing image 
bool Amoeba::checkLine(long int ny, long int nx)
{
  vector<XRange>::const_iterator it;
   
  if ((ny < 0) || (ny >= _size)) return 0;
    for (i = amoebae[ny].begin(); i != amoebae[ny].end(); i++) {
      if ( ( nx <= (*i).xright) && ((*i).xleft <= nx)) {
           return 0;
      }
    }
  return 1;       
}

  
