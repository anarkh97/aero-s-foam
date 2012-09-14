#include <cmath>
#include <iostream>
#include <algorithm>
#ifdef USE_EIGEN3
#include <Eigen/Core>

long h12(long mode, long lpivot, long l1, 
        long m, double *u, long iue, double &up, double *
        c, long ice, long icv, long ncv)
{
//  CONSTRUCTION AND/OR APPLICATION OF A SINGLE
//  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B

//  The original version of this code was developed by
//  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
//  1973 JUN 12, and published in the book
//  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
//  Revised FEB 1995 to accompany reprinting of the book by SIAM.
//     ------------------------------------------------------------------
//                     Subroutine Arguments

//     mode   = 1 OR 2   Selects Algorithm H1 to construct and apply a
//            Householder transformation, or Algorithm H2 to apply a
//            previously constructed transformation.
//     lpivot IS THE INDEX OF THE PIVOT ELEMENT.
//     l1,m   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO
//            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M
//            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
//     u(),iue,up    On entry with MODE = 1, U() contains the pivot
//            vector.  IUE is the storage increment between elements.
//            On exit when MODE = 1, U() and UP contain quantities
//            defining the vector U of the Householder transformation.
//            on entry with MODE = 2, U() and UP should contain
//            quantities previously computed with MODE = 1.  These will
//            not be modified during the entry with MODE = 2.
//     c()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH
//            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE
//            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED.
//            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS.
//     ice    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().
//     icv    STORAGE INCREMENT BETWEEN VECTORS IN C().
//     ncv    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0
//            NO OPERATIONS WILL BE DONE ON C().
//     ------------------------------------------------------------------

    // Builtin functions
    using std::sqrt;
    using std::abs;
    using std::max;
    using std::pow;

    // Local variables
    double b;
    long i, j, i2, i3, i4;
    double cl, sm;
    long incr;
    double clinv;
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>, Eigen::RowMajor > u_ref(u, iue, m);

    // Function Body
    if (0 >= lpivot || lpivot >= l1 || l1 > m) {
        return 0;
    }
    cl = abs(u_ref(0, lpivot-1));
    if (mode == 2) {
        if (cl <= 0.) {
           return 0;
        }
    }
    else {
//                            ****** CONSTRUCT THE TRANSFORMATION. ******
        for (long j = l1-1; j < m; ++j) {
               cl = max(abs(u_ref(0, j)),cl);
        }
        if (cl <= 0.) {
            return 0;
        }

        clinv = 1. / cl;
        sm = pow(u_ref(0, lpivot-1) * clinv, 2);
        for (long j = l1-1; j < m; ++j) {
              sm += pow(u_ref(0, j) * clinv, 2);
        }
        cl *= sqrt(sm);
        if (u_ref(0, lpivot-1) > 0.) {
            cl = -cl;
        }

        up = u_ref(0, lpivot-1) - cl;
        u_ref(0, lpivot-1) = cl;
    }
//            ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******

    if (ncv <= 0) {
        return 0;
    }
    b = up * u_ref(0, lpivot-1);
//                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.

    if (b >= 0.) {
        return 0;
    }

    b = 1. / b;
    i2 = 1 - icv + ice * (lpivot - 1);
    incr = ice * (l1 - lpivot);
    for (long j = 0; j < ncv; ++j) {
        i2 += icv;
        i3 = i2 + incr;
        i4 = i3;
        sm = c[i2-1] * up;
        for (long i = l1-1; i < m; ++i) {
            sm += c[i3-1] * u_ref(0, i);
            i3 += ice;
        }
        if (sm != 0.) {
           sm *= b;
           c[i2-1] += sm * up;
           for (long i = l1-1; i < m; ++i) {
               c[i4-1] += sm * u_ref(0, i);
               i4 += ice;
           }
        }
    }

    return 0;
}
#endif
