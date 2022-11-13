/**************************************************************
 *
 * Test is meant to compare R.K.Kuang's https://github.com/rkkuang/triplelens
 * approach to obtaining of triple lens equating coefficients to the one used in this code.
 * It was verified that the code by Bozza etal. http://www.fisica.unisa.it/GravitationAstrophysics/VBBinaryLensing.htm
 * does no contain the same routine for the triple lenses.
 * The result of the benchmark is that R.K.Kuang's code is superior in performance unless -funsafe-math-optimizations is specified.
 * If the flag is specified, it improves overal performance of Kuang's code by 15% but also improves performance of our routine by
 * more than 100%, in effect outperforming Kuang's code. Keep in mind however, that the unsafe flags do introduce noise that 
 * can be interpreted as an inaccuracy.
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
  double a  = 1.5;
  double b  = 1.5;
  double m2 = 0.25;
  double m3 = 0.25;
  double th = 3.14159;
  clock_t begin; 
  clock_t end; 

  std::complex<double> z1 = {0.0,0.0};
  std::complex<double> z2 = {a,0.0};
  std::complex<double> z3 = {b*cos(th), b*sin(th)};
  std::vector<std::complex<double>> coeffs(11);
  ImgPoint img(a,b,th,m2,m3,0.0,0.0);
  std::vector<double> coeffSums(11);

  double coeffSum = 0.0;

  unsigned int numOfSteps = 100000;
  unsigned int precision = 16;

  // Zero Coeff summs
  std::fill(coeffSums.begin(), coeffSums.end(), 0.0);

  //////////////////////////////////////////////
  // Calculating coeffs without optimization
  //////////////////////////////////////////////
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
      coeffSums[j] += std::abs(coeffs[j]);
      coeffSum += std::abs(coeffs[j]); 
    }
  }
  end = clock();
  std::cout << "10k my coeff calculations:" 
            << double(end - begin) / CLOCKS_PER_SEC 
            << std::setprecision(precision) 
            << "s; Coeff sum = "
            << coeffSum
            << "\n\n";

  // Zero Coeff summs
  std::fill(coeffSums.begin(), coeffSums.end(), 0.0);


  ///////////////////////////////////////////
  // Calculating coeffs with optimalization
  /////////////////////////////////////////
  coeffSum = 0.0;
  begin = clock(); 
  for(unsigned int i = 0; i < numOfSteps; i++)
  {
    // A spiral souce trajectory
    complex<double> source = {(double)i/(double)numOfSteps*cos((double)i/100.0),
                              (double)i/(double)numOfSteps*sin((double)i/100.0)};
    img.setPos(source);
    coeffs = img.getCoeffsOpt();
    for(unsigned j = 0; j < coeffs.size(); j++)
    {
      coeffSums[j] += std::abs(coeffs[j]);
      coeffSum += std::abs(coeffs[j]); 
    }
  }
  end = clock();
  std::cout << "10k my coeff opt calculations:" 
            << std::setprecision(precision) 
            << double(end - begin) / CLOCKS_PER_SEC 
            << "s; Coeff sum = "
            << coeffSum
            << "\n\n";

  for(unsigned j = 0; j < coeffs.size(); j++)
  {
      std::cout << "z-not-separated n:" << j
                << std::setprecision(precision) 
                << ", valueL " << coeffSums[j]
                //<< ", valueL " << std::abs(coeffs[j])
                << "\n"; 
  }

  // Zero Coeff summs
  std::fill(coeffSums.begin(), coeffSums.end(), 0.0);

  /////////////////////////////
  // Split zeta and zeta-independent coefficients
  /////////////////////////////
  coeffSum = 0.0;
  begin = clock();
  vector<complex<double>> coeff0 = img.getCoeffsOptNoZ(); 
  for(unsigned int i = 0; i < numOfSteps; i++)
  {
    // A spiral souce trajectory
    complex<double> source = {(double)i/(double)numOfSteps*cos((double)i/100.0),
                              (double)i/(double)numOfSteps*sin((double)i/100.0)};
    img.setPos(source);
    coeffs = img.getCoeffsOptJustZ();
    for(unsigned j = 0; j < coeffs.size(); j++)
    {
      coeffSums[j] += std::abs(coeffs[j]+coeff0[j]);
      coeffSum += std::abs(coeffs[j]+coeff0[j]); 
    }
  }
  end = clock();
  std::cout << "10k my coeff opt separated calculations:" 
            << std::setprecision(precision) 
            << double(end - begin) / CLOCKS_PER_SEC 
            << "s; Coeff sum = "
            << coeffSum
            << "\n\n";
  
  for(unsigned j = 0; j < coeffs.size(); j++)
  {
      std::cout << "z-separated n:" << j
                << std::setprecision(precision) 
                << ", valueL " << coeffSums[j]
                << "\n"; 
  }

  //KuangCoeffCalculator(double mlens[], std::complex<double> zlens[]);
  double mlens[NLENS] = {1.0-m2-m3, m2, m3};
  std::complex<double> zlens[NLENS] = {z1, z2, z3};
  KuangCoeffCalculator kuangCalc(mlens, zlens);
  std::complex<double> kuangCoeffs[11];

  coeffSum = 0.0;

  // Zero Coeff summs
  std::fill(coeffSums.begin(), coeffSums.end(), 0.0);

  begin = clock(); 
  for(unsigned int i = 0; i < numOfSteps; i++)
  {
    // A spiral souce trajectory
    complex<double> source = {(double)i/(double)numOfSteps*cos((double)i/100.0),
                              (double)i/(double)numOfSteps*sin((double)i/100.0)};
    kuangCalc.polynomialCoefficients(source.real(), source.imag(), kuangCoeffs);

    for(unsigned j = 0; j < 11; j++)
    {
      coeffSums[j] += std::abs(kuangCoeffs[j]);
      coeffSum += std::abs(kuangCoeffs[j]); 
    }

  }
  end = clock();
  std::cout << "10k kuang coeff calculations:"
            << std::setprecision(precision) 
            << double(end - begin) / CLOCKS_PER_SEC 
            << "s; Coeff sum = "
            << coeffSum
            << "\n\n";

  for(unsigned j = 0; j < coeffs.size(); j++)
  {
      std::cout << "z-separated n:" << j
                << std::setprecision(precision) 
                << ", valueL " << coeffSums[j]
                << "\n"; 
  }

  return 0;
}


