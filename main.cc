/**************************************************************
 *
 * This is an example of use of the relevant classes,
 * also used for testing.
 * It mostly demonstrates that the code does not crash on execution.
 *
 *************************************************************/ 

#include <iostream>
#include <complex>
#include <vector>
#include <ctime>

#include <ccc.h>
#include <imgpoint.h>
#include <lcbase.h>
#include <lcirs.h>

int main()
{
  using std::cout;
  double a  = 1.5;
  double b  = 1.5;
  double m2 = 0.25;
  double m3 = 0.25;
  double th = 3.14159;

  // Critical Curves and Caustic test
  cout << "Testing Critical Curve and Caustic generation: \n";
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
       << "s\n\n";


  // Image Point Test

  a = 1.0;
  b = 1.0;
  th = 1.047197551;
  m2 = 1.0/3.0;
  m3 = 1.0/3.0;
 
  cout << "Testing the point images: \n";
  ImgPoint img(a,b,th,m2,m3,0.0,0.0);
  img.getImages();
  cout << "ImgVec size : " << img.imgs.size() << "\n";
  cout << "RootVec size : " << img.roots.size() << "\n";
  cout << "Booleans for completeness size : ";
  for(auto tf: img.isImg) cout << tf; 
  cout  << "\n\n";

  // Image Point Time benchmark
  cout << "Startnig image point calculation: \n";
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
  for(unsigned int i = 0; i < 10000; i++)
  {
    // A spiral souce trajectory
    complex<double> source = {(double)i/10000.0*cos((double)i/100.0), (double)i/10000.0*sin((double)i/100.0)};
    //cout << "Img step " << i << "\n";
    //complex<double> source = {(double)i*0.1, (double)i*0.0577};     
    img.setPos(source);
    img.getImages();
  }
  end = clock();
  cout << "10k Point image calculations (moving case):" 
       << double(end - begin) / CLOCKS_PER_SEC 
       << "s\n\n";

    // Light Curve IRS test 
  cout << "Starting LC test BASE: \n";
  LightCurveBase lcBase(a,b,th,m2,m3,100);
  complex<double> startPoint;
  complex<double> endPoint;
  vector<double>  lightCurve;
  double angle = 0.0;

  // Point source
  begin = clock();  
  for(unsigned int i = 0; i <= 2; i++)
  {
    angle = (double)i/20.0*3.14159;
    endPoint = {cos(angle), sin(angle)};  
    startPoint = -endPoint;
    lcBase.getLC(startPoint, endPoint);
    lightCurve = lcBase.lcVec;
    cout << "\n Priting Light curve" << i << "\n";
    for(auto lcElem: lightCurve)
    {
      cout << lcElem << "\n";
    }

  }
  cout << "\n Light Curve Printed. \n"; 
  
  end = clock();  
  cout << "100 positions point-amp lightcurve:" 
       << double(end - begin) / CLOCKS_PER_SEC 
       << "s\n\n";

  // Extended source IRS
  LightCurveIRS lcIRS(a,b,th,m2,m3, 0.01, 500, 10000);

  begin = clock();  
  for(unsigned int i = 0; i < 100; i++)
  {
    angle = (double)i/200.0*3.14159;
    endPoint = {cos(angle), sin(angle)};  
    startPoint = -endPoint;
    lcIRS.getLCIRS(startPoint, endPoint);
    lightCurve = lcIRS.lcVec; 
  }
  end = clock();
  cout << "Img plane size was: " << lcIRS.amoebae.amoebae.size() << "\n";  
  cout << "100 positions IRS lightcurve:" 
       << double(end - begin) / CLOCKS_PER_SEC 
       << "s\n\n";

  // Laguerre Test
  cout << "Laguerre's Method tests: \n";
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
  complex<double> rootToPolish = {1.001,0.002};
  fourth.laguerre(fouCoef, rootToPolish);
  fou2Sol.push_back(rootToPolish);

  rootToPolish = {2.001,0.002};
  fourth.laguerre(fouCoef, rootToPolish);
  fou2Sol.push_back(rootToPolish);  

  rootToPolish = {3.001,0.002};
  fourth.laguerre(fouCoef, rootToPolish);
  fou2Sol.push_back(rootToPolish); 

  rootToPolish = {4.001,0.002};
  fourth.laguerre(fouCoef, rootToPolish);
  fou2Sol.push_back(rootToPolish);  

  rootToPolish = {-1.001,0.002};
  fourth.laguerre(fouCoef, rootToPolish);
  fou2Sol.push_back(rootToPolish);  

  rootToPolish = {4.00,0.0};
  fourth.laguerre(fouCoef, rootToPolish);
  fou2Sol.push_back(rootToPolish);  

  cout << "fourth order roots (1,2,3,4):" << "\n";
  for(auto sol: fou2Sol)
  {
    cout << sol << "\n";  
  }

  return 0;
}
