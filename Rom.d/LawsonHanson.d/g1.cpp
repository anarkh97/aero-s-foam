#include <cmath>

double d_sign(double a, double b)
{
  double x;
  x = (a >= 0 ? a : -a);
  return( b >= 0 ? x : -x);
}

int g1(double a, double b, double &cterm, 
       double &sterm, double &sig)
{
//     COMPUTE ORTHOGONAL ROTATION MATRIX.. 

//  The original version of this code was developed by
//  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
//  1973 JUN 12, and published in the book
//  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
//  Revised FEB 1995 to accompany reprinting of the book by SIAM.

//     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))
//                        (-S,C)         (-S,C)(B)   (   0          )
//     COMPUTE SIG = SQRT(A**2+B**2)
//        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT
//        SIG MAY BE IN THE SAME LOCATION AS A OR B .
//     ------------------------------------------------------------------
//     ------------------------------------------------------------------

    // Builtin functions
    using std::sqrt;
    using std::abs;

    // Local variables
    double xr, yr;

    // Function body
    if (abs(a) > abs(b)) {
        xr = b / a;
        yr = sqrt(xr*xr + 1.);
        cterm = d_sign(1/yr, a);
        sterm = cterm * xr;
        sig = abs(a) * yr;
        return 0;
    }
    if (b != 0.) {
        xr = a / b;
        yr = sqrt(xr*xr + 1.);
        sterm = d_sign(1/yr, b);
        cterm = sterm * xr;
        sig = abs(b) * yr;
        return 0;
    }
    sig = 0.;
    cterm = 0.;
    sterm = 1.;
    return 0;
}
