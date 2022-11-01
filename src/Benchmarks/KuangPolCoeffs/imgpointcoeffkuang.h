#include <iostream>
#include <complex>
#include <vector>
#include <fstream>
#include <algorithm> 
#include <iomanip>

// output an image
#define TRUE_IMAGE 1
#define FALSE_IMAGE -1
// #define VERBOSE
#define _PRINT_ERRORS2

// #define verbose
// #define combine_verbose
#define combine_verbose_v 0

// #define CALCULATE_ANGLE

#define NLENS 3
#define DEGREE (NLENS*NLENS+1)
#define EPS 1.0e-5 // segment close threshold, 1.0e-5 is ok


#define SOLEPS 1.0e-5 // true or false solution of lens equation solving, 1.0e-5 is a bit too strict at some cases
#define MISS_SOL_THRESHOLD 1

#define JUMP_SEARCH_DEPTH 5 // search depth for the situation that we can not jump

#define one_24 (double)1.0/24

class KuangCoeffCalculator
{
  public:
    int nphi, secnum, basenum, distype, maxmuidx, flag,  pntnum, CQ, finalNPS, ifFinite, area_quality=1;
    // area_quality = 1 means this magnification is ok, otherwise (mainly caused by insufficient samples around the source edge)
    // if area_quality == 0, mean parity in areaFunc might be wrong
    // if area_quality == 2, means return from a looser threshold
    // if area_quality == 3, means return muPS

    double quad_err, quaderr_Tol, area, mu0, mu, relerr_mag, lambda1, lambda2, thetaJ, eacherr, muTotal, xs, ys, phi, SD, MD, CD, muPS, ds, dJ2, RelTolLimb = 1e-3, AbsTolLimb = 1e-4;
    unsigned short int ifjump;
    std::complex<double> *zlens;
    std::complex<double> zr[DEGREE];
    std::complex<double> coefficients[DEGREE + 1];
    std::complex<double> J1, J2, J3, dJ, firstterm, secondterm, thirdterm, J1c2, dJ5, dy, dz;

    // constructor
    KuangCoeffCalculator(double mlens[], std::complex<double> zlens[]);

    void polynomialCoefficients(double xs, double ys, std::complex<double> c[]);
    void multiply_z(std::complex<double> c[], std::complex<double> a, int n);
    void multiply(std::complex<double> a[], int na, std::complex<double> b[], int nb, std::complex<double> c[]);
    void multiply_z_v2(std::complex<double> c[][NLENS + 1], std::complex<double> a, int n, int firstdim);

  private:
    int finalnphi, nimages, ftime, FSflag, nsolution, trackcnt, parity, degrees, adanp0 = 2, MAXNPS = 2000, secnum_priv, basenum_priv, clpairparity1, clpairparity2, special_flag = 0;
    double TINY, ph1, ph2, ph3, adaerrTol = 0.0005, timerrTol = 1e-3, tempdis, mindis, M_PI2 = 2 * M_PI, absdzs;
    // absdzs is used to store the abs(dzs) in trueSolution()
    double r2_1, r2_2, x, y, dx_db, dy_db, dx_db2, dy_db2 , Jxx, Jyy, Jxy, rho2, areaSource, phi0, x1, x2, y1, y2, phi1, phi2, dphi, dphi3, ds1, ds2, subArea, pos_area, rep, imp, x_xj, y_yj, xy_j2, relerr_priv; // trueSolution
    std::complex<double> zs, dzs, dzsdz, zsc, zc[NLENS], z1, z2, z3, z1bar, z2bar, z3bar, zsbar, z13, z23, z33, z12, z22, z32, zsmod, z1barz2barzs, z1barz3barzs, z2barz3barzs, z1barzs, z2barzs, z3barzs, z1barzsmod, z2barzsmod, z3barzsmod, z1barzsbar, z2barzsbar, z3barzsbar, z1barz2barzsmod, z1barz3barzsmod, z2barz3barzsmod, z1barz2barzsbar, z1barz3barzsbar, z2barz3barzsbar;
    std::complex<double> tempzc, tempzc2, tempzc3, tempzc4, J1c, dz1, dz2;
    std::complex<double> p[NLENS + 1][NLENS], q[NLENS + 1][NLENS], p_const[NLENS + 1][NLENS];
    std::complex<double> temp[NLENS + 1], temp_const1[NLENS + 1][NLENS + 1], temp_const2[NLENS + 1][NLENS + 1][NLENS + 1], temp_const22[NLENS + 1];
    std::complex<double> ctemp[DEGREE + 1];
    std::complex<double> qtemp[DEGREE + 1], qtemp2[NLENS + 1];
    std::complex<double> ptemp[DEGREE + 1];
    double therr;
    int NPS;
    double *mlens, *Zlens;
};
