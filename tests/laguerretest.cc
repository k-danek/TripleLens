#include "gtest/gtest.h"
#include "ccc.h"
#include "laguerre.h"

#include <iostream>
#include <complex>
#include <vector>
#include <ctime>
#include <algorithm>
#include <random>

const double eps = 1.0e-10;

bool lowerRealPart(const complex<double>& c1,
                   const complex<double>& c2)
{
  return (c1.real() < c2.real());
}


// Gets a complex number of fixed magnitude and random angle.
complex<double> getRandomComplex(double mag)
{
  // see https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 100.0);
  return std::polar(mag, dis(gen));
}


TEST(findRootsTest, secondOrder)
{
  vector<complex<double>> secCoef = {{1.0,0.0},{2.0,0.0},{1.0,0.0}};
  vector<complex<double>> secSol;
  Laguerre second(secCoef);
  secSol = second.solveRoots();
  ASSERT_EQ(secSol.size(), 2);
  // Check real parts of the roots.
  EXPECT_DOUBLE_EQ(secSol[0].real(), -1.0);
  EXPECT_DOUBLE_EQ(secSol[1].real(), -1.0);

  // Check absolute values of the roots.
  EXPECT_TRUE(abs(secSol[0]) < 1.0+eps);
  EXPECT_TRUE(abs(secSol[0]) < 1.0+eps);


}

TEST(findRootsTest, thirdOrder)
{
  // (x-1)*(x-2)*(x-3)
  // x^3 - 6 x^2 + 11 x - 6
  vector<complex<double>> thiCoef = {{-6.0,0.0},{11.0,0.0},{-6.0,0.0},{1.0,0.0}};
 
  // Expected results with real parts in non-descending order.
  vector<complex<double>> expRes = {{1.0,0.0}, {2.0,0.0}, {3.0,0.0}};
  vector<complex<double>> thiSol;
  Laguerre third(thiCoef);

  thiSol = third.solveRoots();
  ASSERT_EQ(thiSol.size(), 3);

  // Sort results to have real parts in non-descending order.
  std::sort(thiSol.begin(), thiSol.end(), lowerRealPart);

  std::vector<complex<double>>::iterator itExp;
  std::vector<complex<double>>::iterator itSol;

  for(itExp = expRes.begin(), itSol = thiSol.begin();
      itExp != expRes.end(), itSol != thiSol.end();
      ++itExp, ++itSol
     )
  {
    // Check solutions agains expectations
    EXPECT_TRUE(abs(*itExp-*itSol) < eps) << "Expected: " << *itExp 
                                          << ", Solution:" << *itSol;

    // Check imag part of the real roots.
    EXPECT_TRUE(abs(itSol->imag()) < eps);
  }

}

TEST(findRootsTest, fourthOrder)
{
  // x^4 - 10 x^3 + 35 x^2 - 50 x + 24
  // (x-1)*(x-2)*(x-3)*(x-4)
  vector<complex<double>> fouCoef = {{24.0,0.0},{-50.0,0.0},{35.0,0.0},{-10.0,0.0},{1.0,0.0}};
  vector<complex<double>> fouSol;
  Laguerre fourth(fouCoef);
  fouSol = fourth.solveRoots();

  // Expected results with real parts in non-descending order.
  vector<complex<double>> expRes = {{1.0,0.0}, {2.0,0.0}, {3.0,0.0}, {4.0,0.0}};

  fouSol = fourth.solveRoots();
  ASSERT_EQ(fouSol.size(), 4);

  // Sort results to have real parts in non-descending order.
  std::sort(fouSol.begin(), fouSol.end(), lowerRealPart);

  std::vector<complex<double>>::iterator itExp;
  std::vector<complex<double>>::iterator itSol;

  for(itExp = expRes.begin(), itSol = fouSol.begin();
      itExp != expRes.end(), itSol != fouSol.end();
      ++itExp, ++itSol
     )
  {
    // Check solutions agains expectations
    EXPECT_TRUE(abs(*itExp-*itSol) < eps) << "Expected: " << *itExp 
                                          << ", Solution:" << *itSol;

    // Check imag part of the real roots.
    EXPECT_TRUE(abs(itSol->imag()) < eps);
  }
}

// Checks that the polishing works
TEST(polishRootsTest, fourthOrder)
{
  // x^4 - 10 x^3 + 35 x^2 - 50 x + 24
  // (x-1)*(x-2)*(x-3)*(x-4)
  vector<complex<double>> fouCoef = {{24.0,0.0},{-50.0,0.0},{35.0,0.0},{-10.0,0.0},{1.0,0.0}};
  Laguerre fourth(fouCoef);

  // Expected results with real parts in non-descending order.
  vector<complex<double>> expRes = {{1.0,0.0}, {2.0,0.0}, {3.0,0.0}, {4.0,0.0}};

  const double magnitude = 1.0e-2;

  for(complex<double> root: expRes)
  {
    complex<double> polishedRoot = root + getRandomComplex(magnitude);
    fourth.laguerre(fouCoef, polishedRoot);

    // Check solutions agains expectations
    EXPECT_TRUE(abs(root - polishedRoot) < eps) << "Correct root: "   << root 
                                                << ", Polished root:" << polishedRoot;
  }
}


