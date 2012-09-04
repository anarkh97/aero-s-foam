C                           DISCLAIMER
C
C   This file was generated on 10/14/98 by the version of
C   ADIFOR compiled on Aug 21 1995.
C
C   ADIFOR was prepared as an account of work sponsored by an
C   agency of the United States Government, Rice University, and
C   the University of Chicago.  NEITHER THE AUTHOR(S), THE UNITED
C   STATES GOVERNMENT NOR ANY AGENCY THEREOF, NOR RICE UNIVERSITY,
C   NOR THE UNIVERSITY OF CHICAGO, INCLUDING ANY OF THEIR EMPLOYEES
C   OR OFFICERS, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES
C   ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETE-
C   NESS, OR USEFULNESS OF ANY INFORMATION OR PROCESS DISCLOSED, OR
C   REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
C
C=DECK SM3MB
C=PURPOSE Form basic membrane stiffness of 9-dof triangle
C=AUTHOR C. A. Felippa, June 1984
C=VERSION June 1984
C=EQUIPMENT Machine independent
C=KEYWORDS finite element membrane plane stress
C=KEYWORDS basic material stiffness matrix
C=BLOCK ABSTRACT
C
C     SM3MB forms the material element stiffness matrix associated with
C     the basic displacement modes (rigid modes + constant strain
C     modes) of a 9-dof plane-stress triangle based on the
C     free formulation of Bergan and Nygard.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL      SM3MB (X, Y, DM, ALPHA, F, LS, SM, M, STATUS)
C
C     where the input arguments are
C
C       X         (3 x 1) array of x coordinates of triangle nodes.
C       Y         (3 x 1) array of y coordinates of triangle nodes.
C       DM        (3 x 3) matrix relating in-plane forces to strains.
C       ALPHA     Rotational lumping factor; if zero form CST.
C       F         Factor by which stiffness entries will be multiplied.
C       LS        (9 x 1) array of stiffness location pointers
C                 (see Output SM).
C       SM        Incoming material stiffness array.
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        Output stiffness array with basic stiffness
C                 coefficients added in.  The (i,j)-th entry of the
C                 basic element stiffness is added to SM(K,L),
C                 where K=LS(I) and L=LS(J).
C       STATUS    Status character variable.  Blank if no error
C                 detected.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine g_sm3mb(x, g_x, y, g_y, dm, g_dm, 
     * alpha, f, ls, sm, g_sm,  m, status)
C
C                   T Y P E   &   D I M E N S I O N
C
C=BLOCK VAX
C      implicit      none
C=END VAX
        character*(*) status
        integer m, ls(9)
C=BLOCK DOUBLE
        double precision x(3), y(3), dm(3, 3), alpha, f, p(9, 3),
     * sm(m, m)
        double precision area2, c, coef1, coef2
        double precision d11, d12, d13, d22, d23, d33
        double precision x21, x32, x13, y21, y32, y13
        double precision x12, x23, x31, y12, y23, y31
        double precision s1, s2, s3
C=ELSE
C=END DOUBLE
        integer i, j, k, l, n
C
C                   L O G I C
C
        double precision d1_w, d8_b, d4_v, d7_b, d6_b, d5_b, d4_b,
     * d2_b, d2_v, d3_b
        double precision g_x21, g_x(3), g_x12, 
     *g_x32, g_x23, g_x13, g_x31, g_y21
     *, g_y(3), g_y12
        double precision g_y32, g_y23, g_y13,
     * g_y31, g_area2, g_p(9, 3), g_c
     *, g_d11, g_dm(3, 3), g_d22
        double precision g_d33, g_d12, g_d13,
     * g_d23, g_d1_w, g_s1, g_s2, g_s3
     *, g_sm(m, m)
        save g_d23, g_d1_w, g_s1, g_s2, g_s3
        save g_y13, g_y31, g_area2, g_p, g_c, g_d11, g_d22, g_d33,
     * g_d12, g_d13
        save g_x21, g_x12, g_x32, g_x23, g_x13, g_x31, g_y21, g_y12,
     * g_y32, g_y23
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'sm3mb','g_sm3mb.f')
C
        status = ' '
          g_x21 = -g_x(1) + g_x(2)
        x21 = x(2) - x(1)
C--------
          g_x12= -g_x21
        x12 = -x21
C--------
          g_x32 = -g_x(2) + g_x(3)
        x32 = x(3) - x(2)
C--------
          g_x23 = -g_x32
        x23 = -x32
C--------
          g_x13 = -g_x(3) + g_x(1)
        x13 = x(1) - x(3)
C--------
          g_x31 = -g_x13
        x31 = -x13
C--------
          g_y21 = -g_y(1) + g_y(2)
        y21 = y(2) - y(1)
C--------
          g_y12= -g_y21
        y12 = -y21
C--------
          g_y32 = -g_y(2) + g_y(3)
        y32 = y(3) - y(2)
C--------
          g_y23 = -g_y32
        y23 = -y32
C--------
          g_y13 = -g_y(3) + g_y(1)
        y13 = y(1) - y(3)
C--------
          g_y31 = -g_y13
        y31 = -y13
C--------
          g_area2 = -x21 * g_y13 + (-y13) * g_x21 + y21
     * * g_x13 + x13 * g_y21
        area2 = y21 * x13 - x21 * y13
C--------
        if (area2 .le. 0.0) then
          status = 'SM3MB: Zero area'
          if (area2 .eq. 0.0) then
            status = 'SM3MB: Zero area'
          endif
          return
        endif
          g_p(1, 1) = g_y23
        p(1, 1) = y23
C--------
          g_p(2, 1) = 0.0d0
        p(2, 1) = 0.0d0
C--------
          g_p(3, 1) = g_y31
        p(3, 1) = y31
C--------
          g_p(4, 1) = 0.0d0
        p(4, 1) = 0.0d0
C--------
          g_p(5, 1) = g_y12
        p(5, 1) = y12
C--------
          g_p(6, 1) = 0.0d0
        p(6, 1) = 0.0d0
C--------
          g_p(1, 2) = 0.0d0
        p(1, 2) = 0.0d0
C--------
          g_p(2, 2) = g_x32
        p(2, 2) = x32
C--------
          g_p(3, 2) = 0.0d0
        p(3, 2) = 0.0d0
C--------
          g_p(4, 2) = g_x13
        p(4, 2) = x13
C--------
          g_p(5, 2) = 0.0d0
        p(5, 2) = 0.0d0
C--------
          g_p(6, 2) = g_x21
        p(6, 2) = x21
C--------
          g_p(1, 3) = g_x32
        p(1, 3) = x32
C--------
          g_p(2, 3) = g_y23
        p(2, 3) = y23
C--------
          g_p(3, 3) = g_x13
        p(3, 3) = x13
C--------
          g_p(4, 3) = g_y31
        p(4, 3) = y31
C--------
          g_p(5, 3) = g_x21
        p(5, 3) = x21
C--------
          g_p(6, 3) = g_y12
        p(6, 3) = y12
C--------
        n = 6
        if (alpha .ne. 0.0) then
          coef1 = alpha / 6.0
          coef2 = alpha / 3.0
          d4_v = y13 - y21
          d2_b = dble(coef1)
          d3_b = d2_b * d4_v
          d4_b = d2_b * y23
            g_p(7, 1) = -d4_b * g_y21 + d4_b * g_y13 +
     * d3_b * g_y23
          p(7, 1) = y23 * d4_v * dble(coef1)
C--------
          d4_v = x31 - x12
          d2_b = dble(coef1)
          d3_b = d2_b * d4_v
          d4_b = d2_b * x32
            g_p(7, 2) = -d4_b * g_x12+ d4_b * g_x31 +
     * d3_b * g_x32
          p(7, 2) = x32 * d4_v * dble(coef1)
C--------
          d2_b = dble(coef2)
          d5_b = -d2_b * y21
          d6_b = -d2_b * x12
          d7_b = d2_b * y13
          d8_b = d2_b * x31
            g_p(7, 3) = d6_b * g_y21 + d5_b * g_x12+ 
     *d8_b * g_y13 + d7_b * g_x31
          p(7, 3) = (x31 * y13 - x12 * y21) * dble(coef2)
C--------
          d4_v = y21 - y32
          d2_b = dble(coef1)
          d3_b = d2_b * d4_v
          d4_b = d2_b * y31
            g_p(8, 1) = -d4_b * g_y32 + d4_b * g_y21 +
     * d3_b * g_y31
          p(8, 1) = y31 * d4_v * dble(coef1)
C--------
          d4_v = x12 - x23
          d2_b = dble(coef1)
          d3_b = d2_b * d4_v
          d4_b = d2_b * x13
            g_p(8, 2) = -d4_b * g_x23 + d4_b * g_x12+
     * d3_b * g_x13
          p(8, 2) = x13 * d4_v * dble(coef1)
C--------
          d2_b = dble(coef2)
          d5_b = -d2_b * y32
          d6_b = -d2_b * x23
          d7_b = d2_b * y21
          d8_b = d2_b * x12
            g_p(8, 3) = d6_b * g_y32 + d5_b * g_x23 + 
     *d8_b * g_y21 + d7_b * g_x12
          p(8, 3) = (x12 * y21 - x23 * y32) * dble(coef2)
C--------
          d4_v = y32 - y13
          d2_b = dble(coef1)
          d3_b = d2_b * d4_v
          d4_b = d2_b * y12
            g_p(9, 1) = -d4_b * g_y13 + d4_b * g_y32 +
     * d3_b * g_y12
          p(9, 1) = y12 * d4_v * dble(coef1)
C--------
          d4_v = x23 - x31
          d2_b = dble(coef1)
          d3_b = d2_b * d4_v
          d4_b = d2_b * x21
            g_p(9, 2) = -d4_b * g_x31 + d4_b * g_x23 +
     * d3_b * g_x21
          p(9, 2) = x21 * d4_v * dble(coef1)
C--------
          d2_b = dble(coef2)
          d5_b = -d2_b * y13
          d6_b = -d2_b * x31
          d7_b = d2_b * y32
          d8_b = d2_b * x23
            g_p(9, 3) = d6_b * g_y13 + d5_b * g_x31 + 
     *d8_b * g_y32 + d7_b * g_x23
          p(9, 3) = (x23 * y32 - x31 * y13) * dble(coef2)
C--------
          n = 9
        endif
        d2_v = dble(0.5) * f / area2
        d2_b = -d2_v / area2
          g_c = d2_b * g_area2
        c = d2_v
C--------
          g_d11 = c * g_dm(1, 1) + dm(1, 1) * g_c
        d11 = c * dm(1, 1)
C--------
          g_d22 = c * g_dm(2, 2) + dm(2, 2) * g_c
        d22 = c * dm(2, 2)
C--------
          g_d33 = c * g_dm(3, 3) + dm(3, 3) * g_c
        d33 = c * dm(3, 3)
C--------
          g_d12 = c * g_dm(1, 2) + dm(1, 2) * g_c
        d12 = c * dm(1, 2)
C--------
          g_d13 = c * g_dm(1, 3) + dm(1, 3) * g_c
        d13 = c * dm(1, 3)
C--------
          g_d23 = c * g_dm(2, 3) + dm(2, 3) * g_c
        d23 = c * dm(2, 3)
C--------
        do 99998 j = 1, n
          l = ls(j)
            g_d1_w = d12 * g_p(j, 2) + p(j, 2) * g_d12
     * + d11 * g_p(j, 1) + p(j, 1) * g_d11
          d1_w = d11 * p(j, 1) + d12 * p(j, 2)
            g_s1 = d13 * g_p(j, 3) + p(j, 3) * g_d13 +
     * g_d1_w
          s1 = d1_w + d13 * p(j, 3)
C--------
            g_d1_w = d22 * g_p(j, 2) + p(j, 2) * g_d22
     * + d12 * g_p(j, 1) + p(j, 1) * g_d12
          d1_w = d12 * p(j, 1) + d22 * p(j, 2)
            g_s2 = d23 * g_p(j, 3) + p(j, 3) * g_d23 +
     * g_d1_w
          s2 = d1_w + d23 * p(j, 3)
C--------
            g_d1_w = d23 * g_p(j, 2) + p(j, 2) * g_d23
     * + d13 * g_p(j, 1) + p(j, 1) * g_d13
          d1_w = d13 * p(j, 1) + d23 * p(j, 2)
            g_s3 = d33 * g_p(j, 3) + p(j, 3) * g_d33 +
     * g_d1_w
          s3 = d1_w + d33 * p(j, 3)
C--------
          do 99999 i = 1, j
            k = ls(i)
              g_d1_w = s2 * g_p(i, 2) + p(i, 2) * g_s2
     * + s1 * g_p(i, 1) + p(i, 1) * g_s1
            d1_w = s1 * p(i, 1) + s2 * p(i, 2)
              g_sm(k, l) = s3 * g_p(i, 3) + p(i, 3) * g_s3
     * + g_d1_w + g_sm(k, l)
            sm(k, l) = sm(k, l) + (d1_w + s3 * p(i, 3))
C--------
              g_sm(l, k) = g_sm(k, l)
            sm(l, k) = sm(k, l)
C--------
2500        continue
99999     continue
3000      continue
99998   continue
        return
      end
C=END FORTRAN
