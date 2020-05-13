/**************************************************************
 *
 * This is an example of use of the relevant classes,
 * also used for testing.
 *
 *************************************************************/ 

#include <iostream>
#include <complex>
#include <vector>
#include <ctime>

#include <ccc.h>
#include <imgpoint.h>


int main()
{
  using std::cout;
  double a  = 1.5;
  double b  = 1.5;
  double m2 = 0.25;
  double m3 = 0.25;
  double th = 3.14159;
 

  // Critical Curves and Caustic test
  CriticalCurveCaustic ccc(a,b,th,m2,m3,500);
  
  cout << ccc.a;
  ccc.getCC();
  ccc.getCa();
  cout << "ccVec :" << ccc.ccVec.size() << "\n";
  cout << "ccVec :" << ccc.ccVec[0].size() << "\n";
  ccc.printCCC("./bin/CCC.dat");

  // Critical Curves and Caustic Time benchmark
  clock_t begin = clock();  
  for(unsigned int i = 0; i < 1000; i++)
  {
    CriticalCurveCaustic ccc(a,b,th,m2,m3,500);
    ccc.getCC();
    ccc.getCa();
  }
  clock_t end = clock();
  cout << "1k CCC computations:"
       << double(end - begin) / CLOCKS_PER_SEC
       << "s\n";


  // Image Point Test
  ImgPoint img(a,b,th,m2,m3,0.0,0.0);
  img.getImages();
  cout << "ImgVec size : " << img.imgs.size() << "\n";
  cout << "RootVec size : " << img.roots.size() << "\n";
  cout << "Booleans for completeness size : ";
  for(auto tf: img.isImg) cout << tf; 
  cout  << "\n";

  // Image Point Time benchmark
  begin = clock();  
  for(unsigned int i = 0; i < 1000; i++)
  {
    ImgPoint img(a,b,th,m2,m3,0.0,0.0);
    img.getImages();
  }
  end = clock();
  cout << "1k Point image calculations:" 
       << double(end - begin) / CLOCKS_PER_SEC 
       << "s\n";

  begin = clock();  
  for(unsigned int i = 0; i < 1000; i++)
  {
    // A spiral souce trajectory
    complex<double> source = (double)i/1000.0*(sin((double)i/100.0), cos((double)i/100.0));
    img.setPos(source);
    img.getImages();
  }
  end = clock();
  cout << "1k Point image calculations (moving case):" 
       << double(end - begin) / CLOCKS_PER_SEC 
       << "s\n";

  // Laguerre Test
  vector<complex<double>> secCoef = {{1.0,0.0},{2.0,0.0},{1.0,0.0}};
  vector<complex<double>> secSol;
  Laguerre second(secCoef);
  secSol = second.solveRoots();
  cout << "second order roots:" << "\n";
  for(auto sol: secSol)
  {
    cout << sol << "\n";  
  }

  vector<complex<double>> thiCoef = {{0.0,1.0},{0.0,0.0},{0.0,0.0},{1.0,0.0}};
  vector<complex<double>> thiSol;
  Laguerre third(thiCoef);
  thiSol = third.solveRoots();
  cout << "third order roots:" << "\n";
  for(auto sol: thiSol)
  {
    cout << sol << "\n";  
  }

  // (x-1)*(x-2)*(x-3)
  // x^3 - 6 x^2 + 11 x - 6
  vector<complex<double>> thi2Coef = {{-6.0,0.0},{11.0,0.0},{-6.0,0.0},{1.0,0.0}};
  vector<complex<double>> thi2Sol;
  Laguerre third2(thi2Coef);
  thi2Sol = third2.solveRoots();
  cout << "third order roots (1,2,3) :" << "\n";
  for(auto sol: thi2Sol)
  {
    cout << sol << "\n";  
  }
  
  // x^4 - 10 x^3 + 35 x^2 - 50 x + 24
  // (x-1)*(x-2)*(x-3)*(x-4)
  vector<complex<double>> fouCoef = {{24.0,0.0},{-50.0,0.0},{35.0,0.0},{-10.0,0.0},{1.0,0.0}};
  vector<complex<double>> fouSol;
  Laguerre fourth(fouCoef);
  fouSol = fourth.solveRoots();
  cout << "fourth order roots (1,2,3,4):" << "\n";
  for(auto sol: fouSol)
  {
    cout << sol << "\n";  
  }

  // x^4 - 10 x^3 + 35 x^2 - 50 x + 24
  // (x-1)*(x-2)*(x-3)*(x-4)
  vector<complex<double>> fou2Sol;
  fou2Sol.push_back( fourth.laguerre(fouCoef, complex<double>(1.001,0.002)));
  fou2Sol.push_back( fourth.laguerre(fouCoef, complex<double>(2.001,0.002)));
  fou2Sol.push_back( fourth.laguerre(fouCoef, complex<double>(3.001,0.002)));
  fou2Sol.push_back( fourth.laguerre(fouCoef, complex<double>(4.001,0.002)));
  fou2Sol.push_back( fourth.laguerre(fouCoef, complex<double>(-1.001,0.002)));
  fou2Sol.push_back( fourth.laguerre(fouCoef, complex<double>(4.000,0.0)));
  cout << "fourth order roots (1,2,3,4):" << "\n";
  for(auto sol: fou2Sol)
  {
    cout << sol << "\n";  
  }

  return 0;
}
