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
C=DECK SM3MHE
C=PURPOSE Form high-order material stiffness of 9-dof EFF triangle
C=AUTHOR C. A. Felippa
C=VERSION June 1991
C=EQUIPMENT Machine independent
C=KEYWORDS finite element
C=KEYWORDS material stiffness matrix
C=KEYWORDS triangle membrane high-order extended free formulation
C=BLOCK ABSTRACT
C
C     SM3ME forms the higher order stiffness matrix of a 9-dof
C     membrane triangle based on the extended free formulation.
C     This implementation has alphah=5/4 hardwired, and is
C     optimized for maximum formation speed.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL      SM3MHE (X, Y, DM, F, LS, SM, M, STATUS)
C
C     The inputs are:
C
C       X         (3 x 1) array of x coordinates of triangle nodes
C       Y         (3 x 1) array of y coordinates of triangle nodes
C       DM        (3 x 3) matrix constitutive matrix already
C                 integrated through the thickness
C       F         Factor by which all stiffness entries will be multiplied.
C                 It is beta or 0.5*beta
C       SM        Incoming material stiffness array.
C       LS        (9 x 1) array of stiffness location pointers
C                 (see examples in SM3MB).
C                 three rotational DOF will appear at the end.
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        Output stiffness array with higher order stiffness
C                 coefficients added in.
C                 The (i,j)-th entry of the basic element stiffness is added
C                 to SM(K,L), where K=LS(I) and L=LS(J).
C                 (Drilling freedoms are internally 7,8,9)
C
C       STATUS    Status character variable.  Blank if no error
C                 detected.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine g_sm3mhe(x, g_x, y, g_y,dm, g_dm, 
     *f, ls, sm, g_sm, m, status)
C
C                   A R G U M E N T S
C
        integer ls(9), m
        double precision x(3), y(3), dm(3, 3), f, sm(m, m)
        character*(*) status
C
C                   T Y P E   &   D I M E N S I O N
C
        double precision x0, y0, x10, x20, x30, y10, y20, y30
        double precision x12, x21, x23, x32, x31, x13
        double precision y12, y21, y23, y32, y31, y13
        double precision aa12, aa23, aa31, ss12, ss23, ss31, ss1,
     * ss2, ss3
        double precision caa12, caa23, caa31, sum
        double precision ca, cax10, cax20, cax30, cay10, cay20, cay30
        double precision area, area2, kfac
        double precision kqh(6, 6), hmt(6, 3), hqt(6, 3), kth(3, 3)
        double precision s(3), w(6), xyij(6)
        double precision e11, e22, e33, e12, e13, e23
        integer i, j, k, l
C
C                   L O G I C
C
        double precision d12_b, d2_p, d2_w, d1_w, d10_b, d9_b, d8_b,
     * d8_v, d2_v, d3_v
        double precision d9_v, d2_b, d3_b, d4_v, d5_v, d6_v, d7_v,
     * d4_b, d5_b, d6_b
        double precision d7_b, d1_p, g_x12, g_x(1, 3), g_x21
     *, g_x23, g_x32, g_x31, g_x13
     *, g_y12
        double precision g_y(3), g_y21, g_y23, 
     *g_y32, g_y31, g_y13, g_area2, 
     *g_area, g_x0, g_y0
        double precision g_x10, g_x20, g_x30,
     * g_y10, g_y20, g_y30, g_aa12, 
     *g_aa23, g_aa31, g_caa12
        double precision g_caa23, g_caa31, g_ss12
     *, g_ss23, g_ss31, g_ss1, g_ss2
     *, g_ss3, g_cay10, g_cay20
        double precision g_cay30, g_cax10, g_cax20
     *, g_cax30, g_d1_w, g_hmt(6, 3), 
     *g_sum, g_hqt(6, 3), g_kfac, g_e11
     *
        double precision g_dm(3, 3), g_e22, g_e33
     *, g_e12, g_e13, g_e23, g_kqh(6, 6)
     *, g_d2_w, g_kth(3, 3), g_w(6)
        double precision g_s(3), g_ca, g_xyij(6)
     *, g_sm(m, m)
        save g_e23, g_kqh, g_d2_w, g_kth, g_w, g_s, g_ca, g_xyij
        save g_d1_w, g_hmt, g_sum, g_hqt, g_kfac, g_e11, g_e22, g_e33,
     * g_e12, g_e13
        save g_ss31, g_ss1, g_ss2, g_ss3, g_cay10, g_cay20, g_cay30,
     * g_cax10, g_cax20, g_cax30
        save g_y20, g_y30, g_aa12, g_aa23, g_aa31, g_caa12, g_caa23,
     * g_caa31, g_ss12, g_ss23
        save g_y31, g_y13, g_area2, g_area, g_x0, g_y0, g_x10, g_x20,
     * g_x30, g_y10
        save g_x12, g_x21, g_x23, g_x32, g_x31, g_x13, g_y12, g_y21,
     * g_y23, g_y32
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'sm3mhe','g_sm3mhe.f')
C
        status = ' '
        if (f .eq. 0.0) then
          return
        endif
          g_x12 = -g_x(1, 2) + g_x(1, 1)
        x12 = x(1) - x(2)
C--------
          g_x21 = -g_x12
        x21 = -x12
C--------
          g_x23 = -g_x(1, 3) + g_x(1, 2)
        x23 = x(2) - x(3)
C--------
          g_x32 = -g_x23
        x32 = -x23
C--------
          g_x31 = -g_x(1, 1) + g_x(1, 3)
        x31 = x(3) - x(1)
C--------
          g_x13 = -g_x31
        x13 = -x31
C--------
          g_y12 = -g_y(2) + g_y(1)
        y12 = y(1) - y(2)
C--------
          g_y21 = -g_y12
        y21 = -y12
C--------
          g_y23 = -g_y(3) + g_y(2)
        y23 = y(2) - y(3)
C--------
          g_y32 = -g_y23
        y32 = -y23
C--------
          g_y31 = -g_y(1) + g_y(3)
        y31 = y(3) - y(1)
C--------
          g_y13 = -g_y31
        y13 = -y31
C--------
          g_area2 = -x31 * g_y21 + (-y21) * g_x31 + x21
     * * g_y31 + y31 * g_x21
        area2 = x21 * y31 - x31 * y21
C--------
        if (area2 .le. 0.0) then
          status = 'SM3MBE: Negative area'
          if (area2 .eq. 0.0) then
            status = 'SM3MBE: Zero area'
          endif
          return
        endif
          g_area = 0.5d0 * g_area2
        area = 0.5d0 * area2
C--------
        d2_b = 1.0d0 / dble(3.0)
          g_x0 = d2_b * g_x(1, 3) + d2_b * g_x(1, 2) + d2_b 
     ** g_x(1, 1)
        x0 = (x(1) + x(2) + x(3)) / dble(3.0)
C--------
        d2_b = 1.0d0 / dble(3.0)
          g_y0 = d2_b * g_y(3) + d2_b * g_y(2) + d2_b 
     ** g_y(1)
        y0 = (y(1) + y(2) + y(3)) / dble(3.0)
C--------
          g_x10 = -g_x0 + g_x(1, 1)
        x10 = x(1) - x0
C--------
          g_x20 = -g_x0 + g_x(1, 2)
        x20 = x(2) - x0
C--------
          g_x30 = -g_x0 + g_x(1, 3)
        x30 = x(3) - x0
C--------
          g_y10 = -g_y0 + g_y(1)
        y10 = y(1) - y0
C--------
          g_y20 = -g_y0 + g_y(2)
        y20 = y(2) - y0
C--------
          g_y30 = -g_y0 + g_y(3)
        y30 = y(3) - y0
C--------
        d2_v = x30 * x30
        d2_p = 2.0d0 * x30
        d4_v = y30 * y30
        d1_p = 2.0d0 * y30
        d5_b = 2.25d0 * d1_p
        d6_b = 2.25d0 * d2_p
          g_aa12 = d5_b * g_y30 + d6_b * g_x30
        aa12 = 2.25d0 * (d2_v + d4_v)
C--------
        d2_v = x10 * x10
        d2_p = 2.0d0 * x10
        d4_v = y10 * y10
        d1_p = 2.0d0 * y10
        d5_b = 2.25d0 * d1_p
        d6_b = 2.25d0 * d2_p
          g_aa23 = d5_b * g_y10 + d6_b * g_x10
        aa23 = 2.25d0 * (d2_v + d4_v)
C--------
        d2_v = x20 * x20
        d2_p = 2.0d0 * x20
        d4_v = y20 * y20
        d1_p = 2.0d0 * y20
        d5_b = 2.25d0 * d1_p
        d6_b = 2.25d0 * d2_p
          g_aa31 = d5_b * g_y20 + d6_b * g_x20
        aa31 = 2.25d0 * (d2_v + d4_v)
C--------
        d2_v = dble(32.) * aa12
        d3_v = 15.d0 / d2_v
        d3_b = -d3_v / d2_v * dble(32.)       
          g_caa12 = d3_b * g_aa12
        caa12 = d3_v
C--------
        d2_v = dble(32.) * aa23
        d3_v = 15.d0 / d2_v
        d3_b = -d3_v / d2_v * dble(32.)
          g_caa23 = d3_b * g_aa23
        caa23 = d3_v
C--------
        d2_v = dble(32.) * aa31
        d3_v = 15.d0 / d2_v
        d3_b = -d3_v / d2_v * dble(32.)
          g_caa31 = d3_b * g_aa31
        caa31 = d3_v
C--------
        d2_v = x12 * x12
        d2_p = 2.0d0 * x12
        d4_v = y12 * y12
        d1_p = 2.0d0 * y12
          g_ss12 = d1_p * g_y12 + d2_p * g_x12
        ss12 = d2_v + d4_v
C--------
        d2_v = x23 * x23
        d2_p = 2.0d0 * x23
        d4_v = y23 * y23
        d1_p = 2.0d0 * y23
          g_ss23 = d1_p * g_y23 + d2_p * g_x23
        ss23 = d2_v + d4_v
C--------
        d2_v = x31 * x31
        d2_p = 2.0d0 * x31
        d4_v = y31 * y31
        d1_p = 2.0d0 * y31
          g_ss31 = d1_p * g_y31 + d2_p * g_x31
        ss31 = d2_v + d4_v
C--------
          g_ss1 = -0.25d0 * g_ss31 + 0.25d0 * g_ss12
        ss1 = 0.25d0 * (ss12 - ss31)
C--------
          g_ss2 = -0.25d0 * g_ss12 + 0.25d0 * g_ss23
        ss2 = 0.25d0 * (ss23 - ss12)
C--------
          g_ss3 = -0.25d0 * g_ss23 + 0.25d0 * g_ss31
        ss3 = 0.25d0 * (ss31 - ss23)
C--------
          g_cay10 = 0.1875d0 * g_y10
        cay10 = 0.1875d0 * y10
C--------
          g_cay20 = 0.1875d0 * g_y20
        cay20 = 0.1875d0 * y20
C--------
          g_cay30 = 0.1875d0 * g_y30
        cay30 = 0.1875d0 * y30
C--------
          g_cax10 = 0.1875d0 * g_x10
        cax10 = 0.1875d0 * x10
C--------
          g_cax20 = 0.1875d0 * g_x20
        cax20 = 0.1875d0 * x20
C--------
          g_cax30 = 0.1875d0 * g_x30
        cax30 = 0.1875d0 * x30
C--------
        d5_v = -ss3 + 0.6d0 * aa12
        d10_b = y30 * 0.6d0
          g_d1_w = area * g_x30 + x30 * g_area + d5_v 
     ** g_y30 + d10_b * g_aa12 + (-y30) * g_ss3
        d1_w = d5_v * y30 + area * x30
          g_hmt(1, 1) = caa12 * g_d1_w + d1_w * g_caa12
        hmt(1, 1) = caa12 * d1_w
C--------
        d4_b = dble(3.)
          g_hmt(1, 2) = -g_hmt(1, 1) + d4_b * g_cay30
        hmt(1, 2) = dble(3.) * cay30 - hmt(1, 1)
C--------
          g_hmt(1, 3) = g_cay30
        hmt(1, 3) = cay30
C--------
          g_hmt(2, 1) = g_cay10
        hmt(2, 1) = cay10
C--------
        d5_v = -ss1 + 0.6d0 * aa23
        d10_b = y10 * 0.6d0
          g_d1_w = area * g_x10 + x10 * g_area + d5_v 
     ** g_y10 + d10_b * g_aa23 + (-y10) * g_ss1
        d1_w = d5_v * y10 + area * x10
          g_hmt(2, 2) = caa23 * g_d1_w + d1_w * g_caa23
        hmt(2, 2) = caa23 * d1_w
C--------
        d4_b = dble(3.)
          g_hmt(2, 3) = -g_hmt(2, 2) + d4_b * g_cay10
        hmt(2, 3) = dble(3.) * cay10 - hmt(2, 2)
C--------
        d4_v = ss2 + 0.6d0 * aa31
        d10_b = y20 * 0.6d0
          g_d1_w = -area * g_x20 + (-x20) * g_area + d4_v
     * * g_y20 + d10_b * g_aa31 + y20 * g_ss2
        d1_w = d4_v * y20 - area * x20
          g_hmt(3, 1) = caa31 * g_d1_w + d1_w * g_caa31
        hmt(3, 1) = caa31 * d1_w
C--------
          g_hmt(3, 2) = g_cay20
        hmt(3, 2) = cay20
C--------
        d4_b = dble(3.)
          g_hmt(3, 3) = -g_hmt(3, 1) + d4_b * g_cay20
        hmt(3, 3) = dble(3.) * cay20 - hmt(3, 1)
C--------
        d4_v = ss3 - 0.6d0 * aa12
        d10_b = -x30 * 0.6d0
          g_d1_w = area * g_y30 + y30 * g_area + d4_v 
     ** g_x30 + d10_b * g_aa12 + x30 * g_ss3
        d1_w = d4_v * x30 + area * y30
          g_hmt(4, 1) = caa12 * g_d1_w + d1_w * g_caa12
        hmt(4, 1) = caa12 * d1_w
C--------
        d4_b = dble(-3.)
          g_hmt(4, 2) = -g_hmt(4, 1) + d4_b * g_cax30
        hmt(4, 2) = dble(-3.) * cax30 - hmt(4, 1)
C--------
          g_hmt(4, 3) = -g_cax30
        hmt(4, 3) = -cax30
C--------
          g_hmt(5, 1) = -g_cax10
        hmt(5, 1) = -cax10
C--------
        d4_v = ss1 - 0.6d0 * aa23
        d10_b = -x10 * 0.6d0
          g_d1_w = area * g_y10 + y10 * g_area + d4_v 
     ** g_x10 + d10_b * g_aa23 + x10 * g_ss1
        d1_w = d4_v * x10 + area * y10
          g_hmt(5, 2) = caa23 * g_d1_w + d1_w * g_caa23
        hmt(5, 2) = caa23 * d1_w
C--------
        d4_b = dble(-3.)
          g_hmt(5, 3) = -g_hmt(5, 2) + d4_b * g_cax10
        hmt(5, 3) = dble(-3.) * cax10 - hmt(5, 2)
C--------
        d5_v = -ss2 - 0.6d0 * aa31
        d10_b = -x20 * 0.6d0
          g_d1_w = -area * g_y20 + (-y20) * g_area + d5_v
     * * g_x20 + d10_b * g_aa31 + (-x20) * g_ss2
        d1_w = d5_v * x20 - area * y20
          g_hmt(6, 1) = caa31 * g_d1_w + d1_w * g_caa31
        hmt(6, 1) = caa31 * d1_w
C--------
          g_hmt(6, 2) = -g_cax20
        hmt(6, 2) = -cax20
C--------
        d4_b = dble(-3.)
          g_hmt(6, 3) = -g_hmt(6, 1) + d4_b * g_cax20
        hmt(6, 3) = dble(-3.) * cax20 - hmt(6, 1)
C--------
        do 99999 j = 1, 3
          d2_b = 2.d0 / dble(9.)
            g_sum= d2_b * g_hmt(3, j) + d2_b * g_hmt(2, j)
     * + d2_b * g_hmt(1, j)
          sum = 2.d0 / dble(9.) * (hmt(1, j) + hmt(2, j) + hmt(3, j))
C--------
          d4_b = -(4.d0 / dble(3.))
            g_hqt(1, j) = d4_b * g_hmt(1, j) + g_sum
          hqt(1, j) = sum - 4.d0 / dble(3.) * hmt(1, j)
C--------
          d4_b = -(4.d0 / dble(3.))
            g_hqt(2, j) = d4_b * g_hmt(2, j) + g_sum
          hqt(2, j) = sum - 4.d0 / dble(3.) * hmt(2, j)
C--------
          d4_b = -(4.d0 / dble(3.))
            g_hqt(3, j) = d4_b * g_hmt(3, j) + g_sum
          hqt(3, j) = sum - 4.d0 / dble(3.) * hmt(3, j)
C--------
          d2_b = 2.d0 / dble(9.)
            g_sum= d2_b * g_hmt(6, j) + d2_b * g_hmt(5, j)
     * + d2_b * g_hmt(4, j)
          sum = 2.d0 / dble(9.) * (hmt(4, j) + hmt(5, j) + hmt(6, j))
C--------
          d4_b = -(4.d0 / dble(3.))
            g_hqt(4, j) = d4_b * g_hmt(4, j) + g_sum
          hqt(4, j) = sum - 4.d0 / dble(3.) * hmt(4, j)
C--------
          d4_b = -(4.d0 / dble(3.))
            g_hqt(5, j) = d4_b * g_hmt(5, j) + g_sum
          hqt(5, j) = sum - 4.d0 / dble(3.) * hmt(5, j)
C--------
          d4_b = -(4.d0 / dble(3.))
            g_hqt(6, j) = d4_b * g_hmt(6, j) + g_sum
          hqt(6, j) = sum - 4.d0 / dble(3.) * hmt(6, j)
C--------
2000      continue
99999   continue
        d2_v = 1.5d0 * f / area2
        d2_b = -d2_v / area2
          g_kfac = d2_b * g_area2
        kfac = d2_v
C--------
          g_e11 = kfac * g_dm(1, 1) + dm(1, 1) * g_kfac
        e11 = kfac * dm(1, 1)
C--------
          g_e22 = kfac * g_dm(2, 2) + dm(2, 2) * g_kfac
        e22 = kfac * dm(2, 2)
C--------
          g_e33 = kfac * g_dm(3, 3) + dm(3, 3) * g_kfac
        e33 = kfac * dm(3, 3)
C--------
          g_e12 = kfac * g_dm(1, 2) + dm(1, 2) * g_kfac
        e12 = kfac * dm(1, 2)
C--------
          g_e13 = kfac * g_dm(1, 3) + dm(1, 3) * g_kfac
        e13 = kfac * dm(1, 3)
C--------
          g_e23 = kfac * g_dm(2, 3) + dm(2, 3) * g_kfac
        e23 = kfac * dm(2, 3)
C--------
        d3_v = y30 * y30
        d1_p = 2.0d0 * y30
        d6_v = dble(2) * e13
        d8_v = d6_v * x30
        d7_b = -y30 * d6_v
        d8_b = -y30 * x30 * dble(2)
        d5_b = -d8_v + e11 * d1_p
          g_d1_w = d7_b * g_x30 + d8_b * g_e13 + d5_b 
     ** g_y30 + d3_v * g_e11
        d1_w = e11 * d3_v - d8_v * y30
        d4_v = x30 * x30
        d1_p = 2.0d0 * x30
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e33 * d1_p
          g_kqh(1, 1) = d7_b * g_x30 + d5_b * g_e33 + 
     *d2_b * g_d1_w
        kqh(1, 1) = dble(2) * (d1_w + e33 * d4_v)
C--------
        d7_v = e13 * x10 - e11 * y10
        d6_b = -y30 * y10
        d7_b = -y30 * e11
        d8_b = y30 * x10
        d9_b = y30 * e13
          g_d1_w = d7_v * g_y30 + d7_b * g_y10 + d6_b 
     ** g_e11 + d9_b * g_x10 + d8_b * g_e13
        d1_w = d7_v * y30
        d7_v = e13 * y10 - e33 * x10
        d6_b = -x30 * x10
        d7_b = -x30 * e33
        d8_b = x30 * y10
        d9_b = x30 * e13
          g_d2_w = d7_v * g_x30 + d7_b * g_x10 + d6_b 
     ** g_e33 + d9_b * g_y10 + d8_b * g_e13
        d2_w = d7_v * x30
          g_kqh(1, 2) = g_d2_w + g_d1_w
        kqh(1, 2) = d1_w + d2_w
C--------
        d7_v = e13 * x20 - e11 * y20
        d6_b = -y30 * y20
        d7_b = -y30 * e11
        d8_b = y30 * x20
        d9_b = y30 * e13
          g_d1_w = d7_v * g_y30 + d7_b * g_y20 + d6_b 
     ** g_e11 + d9_b * g_x20 + d8_b * g_e13
        d1_w = d7_v * y30
        d7_v = e13 * y20 - e33 * x20
        d6_b = -x30 * x20
        d7_b = -x30 * e33
        d8_b = x30 * y20
        d9_b = x30 * e13
          g_d2_w = d7_v * g_x30 + d7_b * g_x20 + d6_b 
     ** g_e33 + d9_b * g_y20 + d8_b * g_e13
        d2_w = d7_v * x30
          g_kqh(1, 3) = g_d2_w + g_d1_w
        kqh(1, 3) = d1_w + d2_w
C--------
        d3_v = e33 + e12
        d5_v = d3_v * x30
        d4_b = y30 * x30
        d5_b = y30 * d3_v
          g_d1_w = d5_v * g_y30 + d5_b * g_x30 + d4_b 
     ** g_e12 + d4_b * g_e33
        d1_w = d5_v * y30
        d3_v = y30 * y30
        d2_p = 2.0d0 * y30
        d9_v = x30 * x30
        d1_p = 2.0d0 * x30
        d2_b = dble(2)
        d5_b = d2_b * d9_v
        d7_b = d2_b * e23 * d1_p
        d10_b = d2_b * d3_v
        d12_b = d2_b * e13 * d2_p
          g_kqh(1, 4) = d7_b * g_x30 + d5_b * g_e23 + 
     *(-d2_b) * g_d1_w + d12_b * g_y30 + d10_b * g_e13
        kqh(1, 4) = dble(2) * (e13 * d3_v - d1_w + e23 * d9_v)
C--------
        d7_v = e12 * x10 - e13 * y10
        d6_b = -y30 * y10
        d7_b = -y30 * e13
        d8_b = y30 * x10
        d9_b = y30 * e12
          g_d1_w = d7_v * g_y30 + d7_b * g_y10 + d6_b 
     ** g_e13 + d9_b * g_x10 + d8_b * g_e12
        d1_w = d7_v * y30
        d7_v = e33 * y10 - e23 * x10
        d6_b = -x30 * x10
        d7_b = -x30 * e23
        d8_b = x30 * y10
        d9_b = x30 * e33
          g_d2_w = d7_v * g_x30 + d7_b * g_x10 + d6_b 
     ** g_e23 + d9_b * g_y10 + d8_b * g_e33
        d2_w = d7_v * x30
          g_kqh(1, 5) = g_d2_w + g_d1_w
        kqh(1, 5) = d1_w + d2_w
C--------
        d7_v = e12 * x20 - e13 * y20
        d6_b = -y30 * y20
        d7_b = -y30 * e13
        d8_b = y30 * x20
        d9_b = y30 * e12
          g_d1_w = d7_v * g_y30 + d7_b * g_y20 + d6_b 
     ** g_e13 + d9_b * g_x20 + d8_b * g_e12
        d1_w = d7_v * y30
        d7_v = e33 * y20 - e23 * x20
        d6_b = -x30 * x20
        d7_b = -x30 * e23
        d8_b = x30 * y20
        d9_b = x30 * e33
          g_d2_w = d7_v * g_x30 + d7_b * g_x20 + d6_b 
     ** g_e23 + d9_b * g_y20 + d8_b * g_e33
        d2_w = d7_v * x30
          g_kqh(1, 6) = g_d2_w + g_d1_w
        kqh(1, 6) = d1_w + d2_w
C--------
          g_kqh(2, 1) = g_kqh(1, 2)
        kqh(2, 1) = kqh(1, 2)
C--------
        d3_v = y10 * y10
        d1_p = 2.0d0 * y10
        d6_v = dble(2) * e13
        d8_v = d6_v * x10
        d7_b = -y10 * d6_v
        d8_b = -y10 * x10 * dble(2)
        d5_b = -d8_v + e11 * d1_p
          g_d1_w = d7_b * g_x10 + d8_b * g_e13 + d5_b 
     ** g_y10 + d3_v * g_e11
        d1_w = e11 * d3_v - d8_v * y10
        d4_v = x10 * x10
        d1_p = 2.0d0 * x10
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e33 * d1_p
          g_kqh(2, 2) = d7_b * g_x10 + d5_b * g_e33 + 
     *d2_b * g_d1_w
        kqh(2, 2) = dble(2) * (d1_w + e33 * d4_v)
C--------
        d7_v = e13 * x10 - e11 * y10
        d6_b = -y20 * y10
        d7_b = -y20 * e11
        d8_b = y20 * x10
        d9_b = y20 * e13
          g_d1_w = d7_v * g_y20 + d7_b * g_y10 + d6_b 
     ** g_e11 + d9_b * g_x10 + d8_b * g_e13
        d1_w = d7_v * y20
        d7_v = e13 * y10 - e33 * x10
        d6_b = -x20 * x10
        d7_b = -x20 * e33
        d8_b = x20 * y10
        d9_b = x20 * e13
          g_d2_w = d7_v * g_x20 + d7_b * g_x10 + d6_b 
     ** g_e33 + d9_b * g_y10 + d8_b * g_e13
        d2_w = d7_v * x20
          g_kqh(2, 3) = g_d2_w + g_d1_w
        kqh(2, 3) = d1_w + d2_w
C--------
        d7_v = e33 * x10 - e13 * y10
        d6_b = -y30 * y10
        d7_b = -y30 * e13
        d8_b = y30 * x10
        d9_b = y30 * e33
          g_d1_w = d7_v * g_y30 + d7_b * g_y10 + d6_b 
     ** g_e13 + d9_b * g_x10 + d8_b * g_e33
        d1_w = d7_v * y30
        d7_v = e12 * y10 - e23 * x10
        d6_b = -x30 * x10
        d7_b = -x30 * e23
        d8_b = x30 * y10
        d9_b = x30 * e12
          g_d2_w = d7_v * g_x30 + d7_b * g_x10 + d6_b 
     ** g_e23 + d9_b * g_y10 + d8_b * g_e12
        d2_w = d7_v * x30
          g_kqh(2, 4) = g_d2_w + g_d1_w
        kqh(2, 4) = d1_w + d2_w
C--------
        d3_v = e33 + e12
        d5_v = d3_v * x10
        d4_b = y10 * x10
        d5_b = y10 * d3_v
          g_d1_w = d5_v * g_y10 + d5_b * g_x10 + d4_b 
     ** g_e12 + d4_b * g_e33
        d1_w = d5_v * y10
        d3_v = y10 * y10
        d2_p = 2.0d0 * y10
        d9_v = x10 * x10
        d1_p = 2.0d0 * x10
        d2_b = dble(2)
        d5_b = d2_b * d9_v
        d7_b = d2_b * e23 * d1_p
        d10_b = d2_b * d3_v
        d12_b = d2_b * e13 * d2_p
          g_kqh(2, 5) = d7_b * g_x10 + d5_b * g_e23 + 
     *(-d2_b) * g_d1_w + d12_b * g_y10 + d10_b * g_e13
        kqh(2, 5) = dble(2) * (e13 * d3_v - d1_w + e23 * d9_v)
C--------
        d7_v = e33 * x10 - e13 * y10
        d6_b = -y20 * y10
        d7_b = -y20 * e13
        d8_b = y20 * x10
        d9_b = y20 * e33
          g_d1_w = d7_v * g_y20 + d7_b * g_y10 + d6_b 
     ** g_e13 + d9_b * g_x10 + d8_b * g_e33
        d1_w = d7_v * y20
        d7_v = e12 * y10 - e23 * x10
        d6_b = -x20 * x10
        d7_b = -x20 * e23
        d8_b = x20 * y10
        d9_b = x20 * e12
          g_d2_w = d7_v * g_x20 + d7_b * g_x10 + d6_b 
     ** g_e23 + d9_b * g_y10 + d8_b * g_e12
        d2_w = d7_v * x20
          g_kqh(2, 6) = g_d2_w + g_d1_w
        kqh(2, 6) = d1_w + d2_w
C--------
          g_kqh(3, 1) = g_kqh(1, 3)
        kqh(3, 1) = kqh(1, 3)
C--------
          g_kqh(3, 2) = g_kqh(2, 3)
        kqh(3, 2) = kqh(2, 3)
C--------
        d3_v = y20 * y20
        d1_p = 2.0d0 * y20
        d6_v = dble(2) * e13
        d8_v = d6_v * x20
        d7_b = -y20 * d6_v
        d8_b = -y20 * x20 * dble(2)
        d5_b = -d8_v + e11 * d1_p
          g_d1_w = d7_b * g_x20 + d8_b * g_e13 + d5_b 
     ** g_y20 + d3_v * g_e11
        d1_w = e11 * d3_v - d8_v * y20
        d4_v = x20 * x20
        d1_p = 2.0d0 * x20
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e33 * d1_p
          g_kqh(3, 3) = d7_b * g_x20 + d5_b * g_e33 + 
     *d2_b * g_d1_w
        kqh(3, 3) = dble(2) * (d1_w + e33 * d4_v)
C--------
        d7_v = e33 * x20 - e13 * y20
        d6_b = -y30 * y20
        d7_b = -y30 * e13
        d8_b = y30 * x20
        d9_b = y30 * e33
          g_d1_w = d7_v * g_y30 + d7_b * g_y20 + d6_b 
     ** g_e13 + d9_b * g_x20 + d8_b * g_e33
        d1_w = d7_v * y30
        d7_v = e12 * y20 - e23 * x20
        d6_b = -x30 * x20
        d7_b = -x30 * e23
        d8_b = x30 * y20
        d9_b = x30 * e12
          g_d2_w = d7_v * g_x30 + d7_b * g_x20 + d6_b 
     ** g_e23 + d9_b * g_y20 + d8_b * g_e12
        d2_w = d7_v * x30
          g_kqh(3, 4) = g_d2_w + g_d1_w
        kqh(3, 4) = d1_w + d2_w
C--------
        d7_v = e12 * x10 - e13 * y10
        d6_b = -y20 * y10
        d7_b = -y20 * e13
        d8_b = y20 * x10
        d9_b = y20 * e12
          g_d1_w = d7_v * g_y20 + d7_b * g_y10 + d6_b 
     ** g_e13 + d9_b * g_x10 + d8_b * g_e12
        d1_w = d7_v * y20
        d7_v = e33 * y10 - e23 * x10
        d6_b = -x20 * x10
        d7_b = -x20 * e23
        d8_b = x20 * y10
        d9_b = x20 * e33
          g_d2_w = d7_v * g_x20 + d7_b * g_x10 + d6_b 
     ** g_e23 + d9_b * g_y10 + d8_b * g_e33
        d2_w = d7_v * x20
          g_kqh(3, 5) = g_d2_w + g_d1_w
        kqh(3, 5) = d1_w + d2_w
C--------
        d3_v = e33 + e12
        d5_v = d3_v * x20
        d4_b = y20 * x20
        d5_b = y20 * d3_v
          g_d1_w = d5_v * g_y20 + d5_b * g_x20 + d4_b 
     ** g_e12 + d4_b * g_e33
        d1_w = d5_v * y20
        d3_v = y20 * y20
        d2_p = 2.0d0 * y20
        d9_v = x20 * x20
        d1_p = 2.0d0 * x20
        d2_b = dble(2)
        d5_b = d2_b * d9_v
        d7_b = d2_b * e23 * d1_p
        d10_b = d2_b * d3_v
        d12_b = d2_b * e13 * d2_p
          g_kqh(3, 6) = d7_b * g_x20 + d5_b * g_e23 + 
     *(-d2_b) * g_d1_w + d12_b * g_y20 + d10_b * g_e13
        kqh(3, 6) = dble(2) * (e13 * d3_v - d1_w + e23 * d9_v)
C--------
          g_kqh(4, 1) = g_kqh(1, 4)
        kqh(4, 1) = kqh(1, 4)
C--------
          g_kqh(4, 2) = g_kqh(2, 4)
        kqh(4, 2) = kqh(2, 4)
C--------
          g_kqh(4, 3) = g_kqh(3, 4)
        kqh(4, 3) = kqh(3, 4)
C--------
        d3_v = y30 * y30
        d1_p = 2.0d0 * y30
        d6_v = dble(2) * e23
        d8_v = d6_v * x30
        d7_b = -y30 * d6_v
        d8_b = -y30 * x30 * dble(2)
        d5_b = -d8_v + e33 * d1_p
          g_d1_w = d7_b * g_x30 + d8_b * g_e23 + d5_b 
     ** g_y30 + d3_v * g_e33
        d1_w = e33 * d3_v - d8_v * y30
        d4_v = x30 * x30
        d1_p = 2.0d0 * x30
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e22 * d1_p
          g_kqh(4, 4) = d7_b * g_x30 + d5_b * g_e22 + 
     *d2_b * g_d1_w
        kqh(4, 4) = dble(2) * (d1_w + e22 * d4_v)
C--------
        d7_v = e23 * x10 - e33 * y10
        d6_b = -y30 * y10
        d7_b = -y30 * e33
        d8_b = y30 * x10
        d9_b = y30 * e23
          g_d1_w = d7_v * g_y30 + d7_b * g_y10 + d6_b 
     ** g_e33 + d9_b * g_x10 + d8_b * g_e23
        d1_w = d7_v * y30
        d7_v = e23 * y10 - e22 * x10
        d6_b = -x30 * x10
        d7_b = -x30 * e22
        d8_b = x30 * y10
        d9_b = x30 * e23
          g_d2_w = d7_v * g_x30 + d7_b * g_x10 + d6_b 
     ** g_e22 + d9_b * g_y10 + d8_b * g_e23
        d2_w = d7_v * x30
          g_kqh(4, 5) = g_d2_w + g_d1_w
        kqh(4, 5) = d1_w + d2_w
C--------
        d7_v = e23 * x20 - e33 * y20
        d6_b = -y30 * y20
        d7_b = -y30 * e33
        d8_b = y30 * x20
        d9_b = y30 * e23
          g_d1_w = d7_v * g_y30 + d7_b * g_y20 + d6_b 
     ** g_e33 + d9_b * g_x20 + d8_b * g_e23
        d1_w = d7_v * y30
        d7_v = e23 * y20 - e22 * x20
        d6_b = -x30 * x20
        d7_b = -x30 * e22
        d8_b = x30 * y20
        d9_b = x30 * e23
          g_d2_w = d7_v * g_x30 + d7_b * g_x20 + d6_b 
     ** g_e22 + d9_b * g_y20 + d8_b * g_e23
        d2_w = d7_v * x30
          g_kqh(4, 6) = g_d2_w + g_d1_w
        kqh(4, 6) = d1_w + d2_w
C--------
          g_kqh(5, 1) = g_kqh(1, 5)
        kqh(5, 1) = kqh(1, 5)
C--------
          g_kqh(5, 2) = g_kqh(2, 5)
        kqh(5, 2) = kqh(2, 5)
C--------
          g_kqh(5, 3) = g_kqh(3, 5)
        kqh(5, 3) = kqh(3, 5)
C--------
          g_kqh(5, 4) = g_kqh(4, 5)       
        kqh(5, 4) = kqh(4, 5)
C--------
        d3_v = y10 * y10
        d1_p = 2.0d0 * y10
        d6_v = dble(2) * e23
        d8_v = d6_v * x10
        d7_b = -y10 * d6_v
        d8_b = -y10 * x10 * dble(2)
        d5_b = -d8_v + e33 * d1_p
          g_d1_w = d7_b * g_x10 + d8_b * g_e23 + d5_b 
     ** g_y10 + d3_v * g_e33
        d1_w = e33 * d3_v - d8_v * y10
        d4_v = x10 * x10
        d1_p = 2.0d0 * x10
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e22 * d1_p
          g_kqh(5, 5) = d7_b * g_x10 + d5_b * g_e22 + 
     *d2_b * g_d1_w
        kqh(5, 5) = dble(2) * (d1_w + e22 * d4_v)
C--------
        d7_v = e23 * x10 - e33 * y10
        d6_b = -y20 * y10
        d7_b = -y20 * e33
        d8_b = y20 * x10
        d9_b = y20 * e23
          g_d1_w = d7_v * g_y20 + d7_b * g_y10 + d6_b 
     ** g_e33 + d9_b * g_x10 + d8_b * g_e23
        d1_w = d7_v * y20
        d7_v = e23 * y10 - e22 * x10
        d6_b = -x20 * x10
        d7_b = -x20 * e22
        d8_b = x20 * y10
        d9_b = x20 * e23
          g_d2_w = d7_v * g_x20 + d7_b * g_x10 + d6_b 
     ** g_e22 + d9_b * g_y10 + d8_b * g_e23
        d2_w = d7_v * x20
          g_kqh(5, 6) = g_d2_w + g_d1_w
        kqh(5, 6) = d1_w + d2_w
C--------
          g_kqh(6, 1) = g_kqh(1, 6)
        kqh(6, 1) = kqh(1, 6)
C--------
          g_kqh(6, 2) = g_kqh(2, 6)
        kqh(6, 2) = kqh(2, 6)
C--------
          g_kqh(6, 3) = g_kqh(3, 6)        
        kqh(6, 3) = kqh(3, 6)
C--------
          g_kqh(6, 4) = g_kqh(4, 6)
        kqh(6, 4) = kqh(4, 6)
C--------
          g_kqh(6, 5) = g_kqh(5, 6)
        kqh(6, 5) = kqh(5, 6)
C--------
        d3_v = y20 * y20
        d1_p = 2.0d0 * y20
        d6_v = dble(2) * e23
        d8_v = d6_v * x20
        d7_b = -y20 * d6_v
        d8_b = -y20 * x20 * dble(2)
        d5_b = -d8_v + e33 * d1_p
          g_d1_w = d7_b * g_x20 + d8_b * g_e23 + d5_b 
     ** g_y20 + d3_v * g_e33
        d1_w = e33 * d3_v - d8_v * y20
        d4_v = x20 * x20
        d1_p = 2.0d0 * x20
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e22 * d1_p
          g_kqh(6, 6) = d7_b * g_x20 + d5_b * g_e22 + 
     *d2_b * g_d1_w
        kqh(6, 6) = dble(2) * (d1_w + e22 * d4_v)
C--------
          g_kth(1, 1) = 0.0d0
        kth(1, 1) = 0.0d0
C--------
          g_kth(1, 2) = 0.0d0
        kth(1, 2) = 0.0d0
C--------
          g_kth(2, 2) = 0.0d0
        kth(2, 2) = 0.0d0
C--------
          g_kth(1, 3) = 0.0d0
        kth(1, 3) = 0.0d0
C--------
          g_kth(2, 3) = 0.0d0
        kth(2, 3) = 0.0d0
C--------
          g_kth(3, 3) = 0.0d0
        kth(3, 3) = 0.0d0
C--------
        do 99996 j = 1, 3
          do 99998 i = 1, 6
              g_d1_w = kqh(i, 2) * g_hqt(2, j) + hqt(2, j) *
     * g_kqh(i, 2) + kqh(i, 1) * g_hqt(1, j) + hqt(1, j) * 
     *g_kqh(i, 1)
            d1_w = kqh(i, 1) * hqt(1, j) + kqh(i, 2) * hqt(2, j)
              g_d2_w = kqh(i, 4) * g_hqt(4, j) + hqt(4, j) *
     * g_kqh(i, 4) + kqh(i, 3) * g_hqt(3, j) + hqt(3, j) * 
     *g_kqh(i, 3) + g_d1_w
            d2_w = d1_w + kqh(i, 3) * hqt(3, j) + kqh(i, 4) * hqt(4, j)
              g_w(i) = kqh(i, 6) * g_hqt(6, j) + hqt(6, j) *
     * g_kqh(i, 6) + kqh(i, 5) * g_hqt(5, j) + hqt(5, j) * 
     *g_kqh(i, 5) + g_d2_w
            w(i) = d2_w + kqh(i, 5) * hqt(5, j) + kqh(i, 6) * hqt(6, j)
C--------
3200        continue
99998     continue
          do 99997 i = 1, j
              g_d1_w = hqt(2, i) * g_w(2) + w(2) * 
     *g_hqt(2, i) + hqt(1, i) * g_w(1) + w(1) * g_hqt(1, i) + 
     *g_kth(i, j)
            d1_w = kth(i, j) + hqt(1, i) * w(1) + hqt(2, i) * w(2)
              g_d2_w = hqt(4, i) * g_w(4) + w(4) * 
     *g_hqt(4, i) + hqt(3, i) * g_w(3) + w(3) * g_hqt(3, i) + 
     *g_d1_w
            d2_w = d1_w + hqt(3, i) * w(3) + hqt(4, i) * w(4)
              g_kth(i, j) = hqt(6, i) * g_w(6) + w(6) * 
     *g_hqt(6, i) + hqt(5, i) * g_w(5) + w(5) * g_hqt(5, i) 
     *+ g_d2_w
            kth(i, j) = d2_w + hqt(5, i) * w(5) + hqt(6, i) * w(6)
C--------
              g_kth(j, i) = g_kth(i, j)
            kth(j, i) = kth(i, j)
C--------
3300        continue
99997     continue
3500      continue
99996   continue
          g_s(1) = g_kth(1, 3) + g_kth(1, 2) + 
     *g_kth(1, 1)
        s(1) = kth(1, 1) + kth(1, 2) + kth(1, 3)
C--------
          g_s(2) = g_kth(2, 3) + g_kth(2, 2) + 
     *g_kth(2, 1)
        s(2) = kth(2, 1) + kth(2, 2) + kth(2, 3)
C--------
          g_s(3) = g_kth(3, 3) + g_kth(3, 2) + 
     *g_kth(3, 1)
        s(3) = kth(3, 1) + kth(3, 2) + kth(3, 3)
C--------
        d2_v = 0.25d0 / area
        d2_b = -d2_v / area
          g_ca = d2_b * g_area
        ca = d2_v
C--------
          g_xyij(1) = ca * g_x32 + x32 * g_ca
        xyij(1) = ca * x32
C--------
          g_xyij(2) = ca * g_y32 + y32 * g_ca
        xyij(2) = ca * y32
C--------
          g_xyij(3) = ca * g_x13 + x13 * g_ca
        xyij(3) = ca * x13
C--------
          g_xyij(4) = ca * g_y13 + y13 * g_ca
        xyij(4) = ca * y13
C--------
          g_xyij(5) = ca * g_x21 + x21 * g_ca
        xyij(5) = ca * x21
C--------
          g_xyij(6) = ca * g_y21 + y21 * g_ca
        xyij(6) = ca * y21
C--------
        do 99993 j = 1, 9
          l = ls(j)
          do 99995 i = 1, 3
            if (j .le. 6) then
                g_w(i) = s(i) * g_xyij(j) + xyij(j) * 
     *g_s(i)
              w(i) = s(i) * xyij(j)
C--------
            else
                g_w(i) = g_kth(i, j - 6)
              w(i) = kth(i, j - 6)
C--------
            endif
3600        continue
99995     continue
            g_sum= g_w(3) + g_w(2) + g_w(1)
          sum = w(1) + w(2) + w(3)
C--------
          do 99994 i = 1, j
            k = ls(i)
            if (i .le. 6) then
                g_sm(k, l) = sum * g_xyij(i) + xyij(i) * 
     *g_sum+ g_sm(k, l)
              sm(k, l) = sm(k, l) + sum * xyij(i)
C--------
            else
                g_sm(k, l) = g_w(i - 6) + g_sm(k, l)
              sm(k, l) = sm(k, l) + w(i - 6)
C--------
            endif
              g_sm(l, k) = g_sm(k, l)
            sm(l, k) = sm(k, l)
C--------
3700        continue
99994     continue
4000      continue
99993   continue
        return
      end
C=END FORTRAN
C
