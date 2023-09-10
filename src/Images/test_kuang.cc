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

#include <imgpoint.h>
#include <imgpointcoeffkuang.h>

int main()
{
  using std::cout;
  double a  = 0.5;
  double b  = 1.5;
  double m2 = 0.25;
  double m3 = 0.55;
  double th = 3.14159;
  clock_t begin; 
  clock_t end; 

  std::complex<double> z1 = {0.0,0.0};
  std::complex<double> z2 = {a,0.0};
  std::complex<double> z3 = {b*cos(th), b*sin(th)};
  std::vector<std::complex<double>> coeffs;
  ImgPoint img(a,b,th,m2,m3,0.0,0.0);

  double coeffSum = 0.0;

  unsigned int numOfSteps = 10000;

  begin = clock();  
  for(unsigned int i = 0; i < numOfSteps; i++)
  {
    // A spiral souce trajectory
    complex<double> source = {(double)i/(double)numOfSteps*cos((double)i/100.0),
                              (double)i/(double)numOfSteps*sin((double)i/100.0)};
    //cout << "Img step " << i << "\n";
    //complex<double> source = {(double)i*0.1, (double)i*0.0577};     
    img.setPos(source);
    coeffs = img.getCoeffs();
    for(unsigned j = 0; j < coeffs.size(); j++)
    {
      coeffSum += std::abs(coeffs[j]); 
    }
  }
  end = clock();
  std::cout << "10k my coeff calculations:" 
            << double(end - begin) / CLOCKS_PER_SEC 
            << "s; Coeff sum = "
            << coeffSum
            << "\n\n";


  coeffSum = 0.0;

  begin = clock(); 
  for(unsigned int i = 0; i < numOfSteps; i++)
  {
    // A spiral souce trajectory
    complex<double> source = {(double)i/(double)numOfSteps*cos((double)i/100.0),
                              (double)i/(double)numOfSteps*sin((double)i/100.0)};
    //cout << "Img step " << i << "\n";
    //complex<double> source = {(double)i*0.1, (double)i*0.0577};     
    img.setPos(source);
    coeffs = img.getCoeffsOpt();
    for(unsigned j = 0; j < coeffs.size(); j++)
    {
      coeffSum += std::abs(coeffs[j]); 
    }
  }
  end = clock();
  std::cout << "10k my coeff opt calculations:" 
            << double(end - begin) / CLOCKS_PER_SEC 
            << "s; Coeff sum = "
            << coeffSum
            << "\n\n";

  //KuangCoeffCalculator(double mlens[], std::complex<double> zlens[]);
  double mlens[NLENS] = {1.0-m2-m3, m2, m3};
  std::complex<double> zlens[NLENS] = {z1, z2, z3};
  KuangCoeffCalculator kuangCalc(mlens, zlens);
  std::complex<double> kuangCoeffs[11];

  coeffSum = 0.0;

  begin = clock(); 
  for(unsigned int i = 0; i < numOfSteps; i++)
  {
    // A spiral souce trajectory
    complex<double> source = {(double)i/(double)numOfSteps*cos((double)i/100.0),
                              (double)i/(double)numOfSteps*sin((double)i/100.0)};
    kuangCalc.polynomialCoefficients(source.real(), source.imag(), kuangCoeffs);

    for(unsigned j = 0; j < 11; j++)
    {
      coeffSum += std::abs(kuangCoeffs[j]); 
    }

  }
  end = clock();
  std::cout << "10k kuang coeff calculations:" 
            << double(end - begin) / CLOCKS_PER_SEC 
            << "s; Coeff sum = "
            << coeffSum
            << "\n\n";



  return 0;
}
