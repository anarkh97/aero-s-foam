C=DECK TRIGGAUSSQ
C=PURPOSE Get abscissas and weight factors for Gauss triangle rule
C=AUTHOR C. A. Felippa, April 1967
C=VERSION May 1982 (Fortran 77)
C=EQUIPMENT Machine independent
C=KEYWORDS triangle Gauss integration rule absissae weigth
C=BLOCK ABSTRACT
C
C     TRIGGAUSSQ returns the triangular coordinates of
C     sample points and weights for a Gauss-type integration
C     rule over a triangle. Rules are identified by total
C     number of points.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       call      TRIGGAUSSQ  (P, I, ZETA1, ZETA2, ZETA3, WEIGHT)
C
C     Input arguments:
C
C       ABS(P)    Total number of Gauss points in rule:
C                 1, 3 or 7.  If none of these, assume 1.
C                 Sign of P is used to discriminate between two
C                 rules with same number of points.  For ABS(P)<8
C                 the only conflict occurs at 3 points:
C                   P=-3  use midpoint (1/2,1/2,0) rule
C                   P=+3  use (2/3,1/6,1/6) rule.
C       I         Index of sample point in  rule.
C
C     Outputs arguments:
C
C       ZETA1,ZETA2,ZETA3   Triangular coordinates of sample point
C       WEIGHT    Weight factor
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    TRIGGAUSSQ
     $             (p, i, zeta1, zeta2, zeta3, weight)
C
C                   A R G U M E N T S
C
      integer           p, i
      double precision  zeta1, zeta2, zeta3, weight
C
C                   L O C A L   V A R I A B L E S
C
      double precision  zeta(3), sqrt15
      integer           pp
C
C                   L O G I C
C
      pp =   abs(p)
      if (pp .eq. 1)             then
        zeta(1) =   1.0D0/3.
        zeta(2) =   zeta(1)
        zeta(3) =   zeta(2)
        weight =  1.0D0
      else if (pp .eq. 3)        then
        if (p .lt. 0)            then
          zeta(1) =  0.5D0
          zeta(2) =  zeta(1)
          zeta(3) =  zeta(2)
          zeta(i) =  0.0
          weight =   1.0D0/3.
        else
          zeta(1) =  1.D0/6.
          zeta(2) =  zeta(1)
          zeta(3) =  zeta(2)
          zeta(i) =  2.D0/3.
          weight =   1.0D0/3.
        end if
      else if (pp .eq. 7)        then
        sqrt15 =   sqrt(15.D0)
        if (i .eq. 1)            then
          zeta(1) =  1.D0/3.
          zeta(2) =  zeta(1)
          zeta(3) =  zeta(2)
          weight =   9.D0/40.
        else if (i .le. 4)       then
          zeta(1) =  (6.D0-sqrt15)/21.
          zeta(2) =  zeta(1)
          zeta(3) =  zeta(2)
          zeta(i-1) = (9.D0+2.D0*sqrt15)/21.
          weight =   (155.D0-sqrt15)/1200.
        else
          zeta(1) =  (6.D0+sqrt15)/21.
          zeta(2) =  zeta(1)
          zeta(3) =  zeta(2)
          zeta(i-4) = (9.D0-2.D0*sqrt15)/21.
          weight =   (31.D0/120.)-((155.D0-sqrt15)/1200.)
        end if
      else
        zeta(1) =   1.0D0/3.
        zeta(2) =   zeta(1)
        zeta(3) =   zeta(2)
        weight =  1.0D0
      end if
      zeta1 =  zeta(1)
      zeta2 =  zeta(2)
      zeta3 =  zeta(3)
      return
      end
C=END FORTRAN
