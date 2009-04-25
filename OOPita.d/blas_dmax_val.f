       SUBROUTINE BLAS_DMAX_VAL( N, X, INCX, K, R )
*
*      BLAS_DMAX_VAL finds the largest component of x, r, and
*      determines the smallest index, k, such that x(k) = r.
*      Craig Lucas, University of Manchester. June, 2003
*
*      Modified from the BLAS function IDAMAX:
*      Jack dongarra, LINPACK, 3/11/78.
*
*      .. Scalar Arguments ..
*
       DOUBLE PRECISION    R
       INTEGER             INCX, K, N
*      ..
*      .. Array Arguments ..
       DOUBLE PRECISION    X( * )
*      ..
*      .. Local Scalars ..
*
       INTEGER             I, IX
*      ..
       K=0
       IF( N.LT.1 .OR. INCX.LE.0 )
     $    RETURN
       K=1
       IF( N.EQ.1 )
     $    RETURN
       IF( INCX.EQ.1 )
     $    GO TO 30
*
*      Code for increment not equal to 1
*
       IX = 1
       R = X( 1 )
       IX = IX + INCX
       DO 20 I = 2, N
          IF( X( IX ).LE.R )
     $        GO TO 10
          K=I
          R = X( IX )
  10      CONTINUE
          IX = IX + INCX
  20   CONTINUE
       RETURN
*
*      Code for increment equal to 1
*
  30   CONTINUE
       R = X( 1 )
       DO 40 I = 2, N
          IF( X( I ).LE.R )
     $        GO TO 40
          K=I
          R = X( I )
  40  CONTINUE
      RETURN
      END
