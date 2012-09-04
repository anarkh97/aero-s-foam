C=BLOCK ABSTRACT
C
C     TRIG3SHAPE computes the value of the shape functions for a
C     three-noded isoparametric triangle and its
C     x-y derivatives, at a sample point given by its triangle
C     coordinates (zeta1, zeta2, zeta3)
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL  TRIG3SHAPE (HF, ZETA1, ZETA2, ZETA3, X, Y, S, SX, SY, DET)
C
C     Input arguments:
C
C       HF         A dummy character argument.
C      ZETA1,ZETA2,ZETA3  Triangular coordinates of given point
C       X         (3 x 1) array of x coordinates of triangle corners
C       Y         (3 x 1) array of y coordinates of triangle corners
C
C     Outputs arguments:
C
C       S         (3 x 1) array of shape function values
C       SX        (3 x 1) array of shape function x-derivatives
C       SY        (3 x 1) array of shape function y-derivatives
C       DET        Value of Jacobian determinant
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    TRIG3SHAPE
     $       (hf, zeta1, zeta2, zeta3, x, y, s, sx, sy, det)
C
C                   A R G U M E N T S
C
      character*(*)     hf
      double precision  zeta1, zeta2, zeta3, x(3), y(3)
      double precision  s(3), sx(3), sy(3), det
C
C                   L O C A L   V A R I A B L E S
C  
      double precision   x23, y23
      double precision   cdet
C
C                   L O G I C
C
      s(1) = zeta1
      s(2) = zeta2
      s(3) = zeta3

      y23 = y(2) - y(3)
      x23 = x(2) - x(3)
      
      det = x(1)*y23 - y(1)*x23 + x(2)*y(3) - x(3)*y(2)

      if (det .le. 0.0)       then
        return
      end if

      cdet =      1.0 /det

      sx(1) =   y23 * cdet
      sy(1) = - x23 * cdet

      sx(2) =   cdet * ( y(3) - y(1) )
      sy(2) =   cdet * ( x(1) - x(3) )
    
      sx(3) =   cdet * ( y(1) - y(2) )  
      sy(3) =   cdet * ( x(2) - x(1) ) 
 
      return

      end
