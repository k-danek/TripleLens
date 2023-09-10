#include <imgpointcoeffkuang.h>

#define NLENS 3

KuangCoeffCalculator::KuangCoeffCalculator(double mlens[], std::complex<double> zlens[]) {
    nphi = 32;
    secnum = 45;
    basenum = 1;
    quad_err = 0;
    quaderr_Tol = 1e-4;
    CQ = 1;
    relerr_mag = 1e-3;
    TINY = 1.0e-20;


    this -> mlens = mlens;
    this -> zlens = zlens;

    //used in polynomialCoefficients
    for (int i = 0; i < NLENS; i++) {
        zc[i] = conj(zlens[i]);
    }

    temp_const1[0][0] = 1.0;
    temp_const2[0][0][0] = 1.0;
    for (int i = 0; i < NLENS; i++) {

        for (int j = 0; j < NLENS; j++) {
            multiply_z_v2(temp_const1, zlens[j], j, i);
        }

        for (int i2 = 0; i2 <= NLENS; i2++) temp_const1[i + 1][i2] = temp_const1[i][i2];
        temp_const1[i + 1][0] = 1.0;


        //  //numerator, m[i] * (z-z1)*(z-z2) ... (z-zn)
        for (int j = 0; j <= NLENS; j++) {
            p_const[j][i] = mlens[i] * temp_const1[i][j];
        }

        for (int j = 0; j < NLENS; j++) {
            /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */

            degrees = 0;
            for (int k = 0; k < NLENS; k++) {
                if (k == j) continue;
                multiply_z_v2(temp_const2[i], zlens[k], degrees, j);
                degrees++;
            }

            for (int i2 = 0; i2 <= NLENS; i2++) temp_const2[i][j + 1][i2] = temp_const2[i][j][i2];
            temp_const2[i][j + 1][0] = 1.0;

        }

        for (int i2 = 0; i2 <= NLENS; i2++) {
            temp_const2[i + 1][0][i2] = temp_const2[i][NLENS][i2];
        }
        temp_const2[i + 1][0][0] = 1.0;
    }


}


// obtain the polynomical coefficients
// input: lens mass and positions, mlens[], zlens[], source positions: xs, ys
// output: c
void KuangCoeffCalculator::polynomialCoefficients(double xs, double ys, std::complex<double> c[])
{
    int i, j, k;
    // complex zs;
    zs = std::complex<double>(xs, ys);
    zsc = conj(zs);

    for (i = 0; i < NLENS; i++) {
        /* denominator */
        for (j = 0; j <= NLENS; j++) { /* (zsc-conjugate(z[i])) * Product_{j=1}^N (z-z[j]) */
            q[j][i] = (zsc - zc[i]) * temp_const1[i][j];
        }

        /* Sum_{j=1}^n Product (z-z_k), k=1, n, but k !=j. */
        for (j = 0; j < NLENS; j++) {
            /* coefficient for  Product (z-z_k), k=1, n, but k !=j. This is a polynomial of DEGREE n-1 */

            /* doing the sum */
            for (k = 0; k < NLENS; k++) {
                q[k][i] = q[k][i] + mlens[j] * temp_const2[i][j][k];
                // q[k][i] = q[k][i] + mlens[j] * temp[k];
            }
        }
    }

    /* now clear the fractions, z-zs - Sum_i p_i/q_i */
    /* get the polynomial Product q_i, i=1, n */

    /* first term */
    qtemp[0] = 1.0;
    degrees = 0;
    for(i = 0; i < NLENS; i++) {
        for (j = 0; j <= NLENS; j++) {
            qtemp2[j] = q[j][i];
        }

        multiply(qtemp, degrees, qtemp2, NLENS, ctemp);

        degrees += NLENS;

        for (j = 0; j <= degrees; j++) {
            qtemp[j] = ctemp[j];
        }
    }

    /* get coefficients of (z-zs) Product_i=1^n q_i */
    multiply_z(ctemp, zs, degrees);

    /* copy the coefficients */
    for (i = 0; i < DEGREE + 1; i++) {
        c[i] = ctemp[i];
    }

    /* second term */
    for (i = 0; i < NLENS; i++) {
        degrees = 0;
        qtemp[0] = 1.0;
        for (j = 0; j < NLENS; j++) {
            if (j == i) continue;

            for (k = 0; k <= NLENS; k++) {
                qtemp2[k] = q[k][j];
            }

            multiply(qtemp, degrees, qtemp2, NLENS, ctemp);

            degrees += NLENS;

            for (k = 0; k <= degrees; k++) {
                qtemp[k] = ctemp[k];
            }
        }

        for(k = 0; k <= NLENS; k++) {
            ptemp[k] = p_const[k][i];
        }

        multiply(qtemp, degrees, ptemp, NLENS, ctemp);

        for(k = 0; k < DEGREE; k++) {
            c[k] = c[k] - ctemp[k];
        }
    }
}

/* multiply a polynomial of degree n by (z-a), returns a polynomial of degree n+1 */
void KuangCoeffCalculator::multiply_z(std::complex<double> c[], std::complex<double> a, int n)
{
    int j;

    c[n + 1] = c[n];
    for (j = n; j >= 1; j--) c[j] = c[j - 1] - c[j] * a;
    c[0] = c[0] * (-a);
};

/* multiply two polynomials of degree of na and nb, input as a(na), b(nb) */
void KuangCoeffCalculator::multiply(std::complex<double> a[],
                                    int na,
                                    std::complex<double> b[],
                                    int nb,
                                    std::complex<double> c[])
{
    int i, j;

    /* zero the array */
    for (i = 0; i <= na + nb; i++) {
        c[i] = 0.0;
    }

    /* now do the product */
    for (i = 0; i <= na; i++) {
        for (j = 0; j <= nb; j++) {
            c[i + j] = c[i + j] + a[i] * b[j];
        }
    }
}

void KuangCoeffCalculator::multiply_z_v2(std::complex<double> c[][NLENS + 1], std::complex<double> a, int n, int firstdim)
{
    int j;

    c[firstdim + 1][n + 1] = c[firstdim][n];
    if (firstdim == 0) {
        c[0][n + 1] = c[0][n];
        for (j = n; j >= 1; j--) {
            c[0][j] = c[0][j - 1] - c[0][j] * a;
        }
        c[0][0] = c[0][0] * (-a);
    } else
    {
        // for (int i = 0; i<=NLENS; i++) c[firstdim][i] = c[firstdim-1][i];
        // c[firstdim][0] = 1.0;

        c[firstdim][n + 1] = c[firstdim][n];
        for (j = n; j >= 1; j--) {
            c[firstdim][j] = c[firstdim][j - 1] - c[firstdim][j] * a;
        }
        c[firstdim][0] = c[firstdim][0] * (-a);

    }

}