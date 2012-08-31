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
      subroutine g_sm3mhe(g_p_, x, g_x, ldg_x, y, g_y, ldg_y, dm, g_dm, 
     *ldg_dm, f, ls, sm, g_sm, ldg_sm, m, status)
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
        double precision aa12, aa23, aa31, ss12, ss23, ss31, ss1, ss2, s
     *s3
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
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_x, ldg_y, ldg_dm, ldg_sm
        double precision d12_b, d2_p, d2_w, d1_w, d10_b, d9_b, d8_b, d8_
     *v, d2_v, d3_v
        double precision d9_v, d2_b, d3_b, d4_v, d5_v, d6_v, d7_v, d4_b,
     * d5_b, d6_b
        double precision d7_b, d1_p, g_x12(g_pmax_), g_x(ldg_x, 3), g_x2
     *1(g_pmax_), g_x23(g_pmax_), g_x32(g_pmax_), g_x31(g_pmax_), g_x13(
     *g_pmax_), g_y12(g_pmax_)
        double precision g_y(ldg_y, 3), g_y21(g_pmax_), g_y23(g_pmax_), 
     *g_y32(g_pmax_), g_y31(g_pmax_), g_y13(g_pmax_), g_area2(g_pmax_), 
     *g_area(g_pmax_), g_x0(g_pmax_), g_y0(g_pmax_)
        double precision g_x10(g_pmax_), g_x20(g_pmax_), g_x30(g_pmax_),
     * g_y10(g_pmax_), g_y20(g_pmax_), g_y30(g_pmax_), g_aa12(g_pmax_), 
     *g_aa23(g_pmax_), g_aa31(g_pmax_), g_caa12(g_pmax_)
        double precision g_caa23(g_pmax_), g_caa31(g_pmax_), g_ss12(g_pm
     *ax_), g_ss23(g_pmax_), g_ss31(g_pmax_), g_ss1(g_pmax_), g_ss2(g_pm
     *ax_), g_ss3(g_pmax_), g_cay10(g_pmax_), g_cay20(g_pmax_)
        double precision g_cay30(g_pmax_), g_cax10(g_pmax_), g_cax20(g_p
     *max_), g_cax30(g_pmax_), g_d1_w(g_pmax_), g_hmt(g_pmax_, 6, 3), g_
     *sum(g_pmax_), g_hqt(g_pmax_, 6, 3), g_kfac(g_pmax_), g_e11(g_pmax_
     *)
        double precision g_dm(ldg_dm, 3, 3), g_e22(g_pmax_), g_e33(g_pma
     *x_), g_e12(g_pmax_), g_e13(g_pmax_), g_e23(g_pmax_), g_kqh(g_pmax_
     *, 6, 6), g_d2_w(g_pmax_), g_kth(g_pmax_, 3, 3), g_w(g_pmax_, 6)
        double precision g_s(g_pmax_, 3), g_ca(g_pmax_), g_xyij(g_pmax_,
     * 6), g_sm(ldg_sm, m, m)
        save g_e23, g_kqh, g_d2_w, g_kth, g_w, g_s, g_ca, g_xyij
        save g_d1_w, g_hmt, g_sum, g_hqt, g_kfac, g_e11, g_e22, g_e33, g
     *_e12, g_e13
        save g_ss31, g_ss1, g_ss2, g_ss3, g_cay10, g_cay20, g_cay30, g_c
     *ax10, g_cax20, g_cax30
        save g_y20, g_y30, g_aa12, g_aa23, g_aa31, g_caa12, g_caa23, g_c
     *aa31, g_ss12, g_ss23
        save g_y31, g_y13, g_area2, g_area, g_x0, g_y0, g_x10, g_x20, g_
     *x30, g_y10
        save g_x12, g_x21, g_x23, g_x32, g_x31, g_x13, g_y12, g_y21, g_y
     *23, g_y32
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'sm3mhe','g_sm3mhe.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        status = ' '
        if (f .eq. 0.0) then
          return
        endif
        do g_i_ = 1, g_p_
          g_x12(g_i_) = -g_x(g_i_, 2) + g_x(g_i_, 1)
        enddo
        x12 = x(1) - x(2)
C--------
        do g_i_ = 1, g_p_
          g_x21(g_i_) = -g_x12(g_i_)
        enddo
        x21 = -x12
C--------
        do g_i_ = 1, g_p_
          g_x23(g_i_) = -g_x(g_i_, 3) + g_x(g_i_, 2)
        enddo
        x23 = x(2) - x(3)
C--------
        do g_i_ = 1, g_p_
          g_x32(g_i_) = -g_x23(g_i_)
        enddo
        x32 = -x23
C--------
        do g_i_ = 1, g_p_
          g_x31(g_i_) = -g_x(g_i_, 1) + g_x(g_i_, 3)
        enddo
        x31 = x(3) - x(1)
C--------
        do g_i_ = 1, g_p_
          g_x13(g_i_) = -g_x31(g_i_)
        enddo
        x13 = -x31
C--------
        do g_i_ = 1, g_p_
          g_y12(g_i_) = -g_y(g_i_, 2) + g_y(g_i_, 1)
        enddo
        y12 = y(1) - y(2)
C--------
        do g_i_ = 1, g_p_
          g_y21(g_i_) = -g_y12(g_i_)
        enddo
        y21 = -y12
C--------
        do g_i_ = 1, g_p_
          g_y23(g_i_) = -g_y(g_i_, 3) + g_y(g_i_, 2)
        enddo
        y23 = y(2) - y(3)
C--------
        do g_i_ = 1, g_p_
          g_y32(g_i_) = -g_y23(g_i_)
        enddo
        y32 = -y23
C--------
        do g_i_ = 1, g_p_
          g_y31(g_i_) = -g_y(g_i_, 1) + g_y(g_i_, 3)
        enddo
        y31 = y(3) - y(1)
C--------
        do g_i_ = 1, g_p_
          g_y13(g_i_) = -g_y31(g_i_)
        enddo
        y13 = -y31
C--------
        do g_i_ = 1, g_p_
          g_area2(g_i_) = -x31 * g_y21(g_i_) + (-y21) * g_x31(g_i_) + x2
     *1 * g_y31(g_i_) + y31 * g_x21(g_i_)
        enddo
        area2 = x21 * y31 - x31 * y21
C--------
        if (area2 .le. 0.0) then
          status = 'SM3MBE: Negative area'
          if (area2 .eq. 0.0) then
            status = 'SM3MBE: Zero area'
          endif
          return
        endif
        do g_i_ = 1, g_p_
          g_area(g_i_) = 0.5d0 * g_area2(g_i_)
        enddo
        area = 0.5d0 * area2
C--------
        d2_b = 1.0d0 / dble(3.0)
        do g_i_ = 1, g_p_
          g_x0(g_i_) = d2_b * g_x(g_i_, 3) + d2_b * g_x(g_i_, 2) + d2_b 
     ** g_x(g_i_, 1)
        enddo
        x0 = (x(1) + x(2) + x(3)) / dble(3.0)
C--------
        d2_b = 1.0d0 / dble(3.0)
        do g_i_ = 1, g_p_
          g_y0(g_i_) = d2_b * g_y(g_i_, 3) + d2_b * g_y(g_i_, 2) + d2_b 
     ** g_y(g_i_, 1)
        enddo
        y0 = (y(1) + y(2) + y(3)) / dble(3.0)
C--------
        do g_i_ = 1, g_p_
          g_x10(g_i_) = -g_x0(g_i_) + g_x(g_i_, 1)
        enddo
        x10 = x(1) - x0
C--------
        do g_i_ = 1, g_p_
          g_x20(g_i_) = -g_x0(g_i_) + g_x(g_i_, 2)
        enddo
        x20 = x(2) - x0
C--------
        do g_i_ = 1, g_p_
          g_x30(g_i_) = -g_x0(g_i_) + g_x(g_i_, 3)
        enddo
        x30 = x(3) - x0
C--------
        do g_i_ = 1, g_p_
          g_y10(g_i_) = -g_y0(g_i_) + g_y(g_i_, 1)
        enddo
        y10 = y(1) - y0
C--------
        do g_i_ = 1, g_p_
          g_y20(g_i_) = -g_y0(g_i_) + g_y(g_i_, 2)
        enddo
        y20 = y(2) - y0
C--------
        do g_i_ = 1, g_p_
          g_y30(g_i_) = -g_y0(g_i_) + g_y(g_i_, 3)
        enddo
        y30 = y(3) - y0
C--------
        d2_v = x30 * x30
        d2_p = 2.0d0 * x30
        d4_v = y30 * y30
        d1_p = 2.0d0 * y30
        d5_b = 2.25d0 * d1_p
        d6_b = 2.25d0 * d2_p
        do g_i_ = 1, g_p_
          g_aa12(g_i_) = d5_b * g_y30(g_i_) + d6_b * g_x30(g_i_)
        enddo
        aa12 = 2.25d0 * (d2_v + d4_v)
C--------
        d2_v = x10 * x10
        d2_p = 2.0d0 * x10
        d4_v = y10 * y10
        d1_p = 2.0d0 * y10
        d5_b = 2.25d0 * d1_p
        d6_b = 2.25d0 * d2_p
        do g_i_ = 1, g_p_
          g_aa23(g_i_) = d5_b * g_y10(g_i_) + d6_b * g_x10(g_i_)
        enddo
        aa23 = 2.25d0 * (d2_v + d4_v)
C--------
        d2_v = x20 * x20
        d2_p = 2.0d0 * x20
        d4_v = y20 * y20
        d1_p = 2.0d0 * y20
        d5_b = 2.25d0 * d1_p
        d6_b = 2.25d0 * d2_p
        do g_i_ = 1, g_p_
          g_aa31(g_i_) = d5_b * g_y20(g_i_) + d6_b * g_x20(g_i_)
        enddo
        aa31 = 2.25d0 * (d2_v + d4_v)
C--------
        d2_v = dble(32.) * aa12
        d3_v = 15.d0 / d2_v
        d3_b = -d3_v / d2_v * dble(32.)
        do g_i_ = 1, g_p_
          g_caa12(g_i_) = d3_b * g_aa12(g_i_)
        enddo
        caa12 = d3_v
C--------
        d2_v = dble(32.) * aa23
        d3_v = 15.d0 / d2_v
        d3_b = -d3_v / d2_v * dble(32.)
        do g_i_ = 1, g_p_
          g_caa23(g_i_) = d3_b * g_aa23(g_i_)
        enddo
        caa23 = d3_v
C--------
        d2_v = dble(32.) * aa31
        d3_v = 15.d0 / d2_v
        d3_b = -d3_v / d2_v * dble(32.)
        do g_i_ = 1, g_p_
          g_caa31(g_i_) = d3_b * g_aa31(g_i_)
        enddo
        caa31 = d3_v
C--------
        d2_v = x12 * x12
        d2_p = 2.0d0 * x12
        d4_v = y12 * y12
        d1_p = 2.0d0 * y12
        do g_i_ = 1, g_p_
          g_ss12(g_i_) = d1_p * g_y12(g_i_) + d2_p * g_x12(g_i_)
        enddo
        ss12 = d2_v + d4_v
C--------
        d2_v = x23 * x23
        d2_p = 2.0d0 * x23
        d4_v = y23 * y23
        d1_p = 2.0d0 * y23
        do g_i_ = 1, g_p_
          g_ss23(g_i_) = d1_p * g_y23(g_i_) + d2_p * g_x23(g_i_)
        enddo
        ss23 = d2_v + d4_v
C--------
        d2_v = x31 * x31
        d2_p = 2.0d0 * x31
        d4_v = y31 * y31
        d1_p = 2.0d0 * y31
        do g_i_ = 1, g_p_
          g_ss31(g_i_) = d1_p * g_y31(g_i_) + d2_p * g_x31(g_i_)
        enddo
        ss31 = d2_v + d4_v
C--------
        do g_i_ = 1, g_p_
          g_ss1(g_i_) = -0.25d0 * g_ss31(g_i_) + 0.25d0 * g_ss12(g_i_)
        enddo
        ss1 = 0.25d0 * (ss12 - ss31)
C--------
        do g_i_ = 1, g_p_
          g_ss2(g_i_) = -0.25d0 * g_ss12(g_i_) + 0.25d0 * g_ss23(g_i_)
        enddo
        ss2 = 0.25d0 * (ss23 - ss12)
C--------
        do g_i_ = 1, g_p_
          g_ss3(g_i_) = -0.25d0 * g_ss23(g_i_) + 0.25d0 * g_ss31(g_i_)
        enddo
        ss3 = 0.25d0 * (ss31 - ss23)
C--------
        do g_i_ = 1, g_p_
          g_cay10(g_i_) = 0.1875d0 * g_y10(g_i_)
        enddo
        cay10 = 0.1875d0 * y10
C--------
        do g_i_ = 1, g_p_
          g_cay20(g_i_) = 0.1875d0 * g_y20(g_i_)
        enddo
        cay20 = 0.1875d0 * y20
C--------
        do g_i_ = 1, g_p_
          g_cay30(g_i_) = 0.1875d0 * g_y30(g_i_)
        enddo
        cay30 = 0.1875d0 * y30
C--------
        do g_i_ = 1, g_p_
          g_cax10(g_i_) = 0.1875d0 * g_x10(g_i_)
        enddo
        cax10 = 0.1875d0 * x10
C--------
        do g_i_ = 1, g_p_
          g_cax20(g_i_) = 0.1875d0 * g_x20(g_i_)
        enddo
        cax20 = 0.1875d0 * x20
C--------
        do g_i_ = 1, g_p_
          g_cax30(g_i_) = 0.1875d0 * g_x30(g_i_)
        enddo
        cax30 = 0.1875d0 * x30
C--------
        d5_v = -ss3 + 0.6d0 * aa12
        d10_b = y30 * 0.6d0
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = area * g_x30(g_i_) + x30 * g_area(g_i_) + d5_v 
     ** g_y30(g_i_) + d10_b * g_aa12(g_i_) + (-y30) * g_ss3(g_i_)
        enddo
        d1_w = d5_v * y30 + area * x30
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 1, 1) = caa12 * g_d1_w(g_i_) + d1_w * g_caa12(g_i_
     *)
        enddo
        hmt(1, 1) = caa12 * d1_w
C--------
        d4_b = dble(3.)
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 1, 2) = -g_hmt(g_i_, 1, 1) + d4_b * g_cay30(g_i_)
        enddo
        hmt(1, 2) = dble(3.) * cay30 - hmt(1, 1)
C--------
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 1, 3) = g_cay30(g_i_)
        enddo
        hmt(1, 3) = cay30
C--------
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 2, 1) = g_cay10(g_i_)
        enddo
        hmt(2, 1) = cay10
C--------
        d5_v = -ss1 + 0.6d0 * aa23
        d10_b = y10 * 0.6d0
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = area * g_x10(g_i_) + x10 * g_area(g_i_) + d5_v 
     ** g_y10(g_i_) + d10_b * g_aa23(g_i_) + (-y10) * g_ss1(g_i_)
        enddo
        d1_w = d5_v * y10 + area * x10
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 2, 2) = caa23 * g_d1_w(g_i_) + d1_w * g_caa23(g_i_
     *)
        enddo
        hmt(2, 2) = caa23 * d1_w
C--------
        d4_b = dble(3.)
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 2, 3) = -g_hmt(g_i_, 2, 2) + d4_b * g_cay10(g_i_)
        enddo
        hmt(2, 3) = dble(3.) * cay10 - hmt(2, 2)
C--------
        d4_v = ss2 + 0.6d0 * aa31
        d10_b = y20 * 0.6d0
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = -area * g_x20(g_i_) + (-x20) * g_area(g_i_) + d
     *4_v * g_y20(g_i_) + d10_b * g_aa31(g_i_) + y20 * g_ss2(g_i_)
        enddo
        d1_w = d4_v * y20 - area * x20
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 3, 1) = caa31 * g_d1_w(g_i_) + d1_w * g_caa31(g_i_
     *)
        enddo
        hmt(3, 1) = caa31 * d1_w
C--------
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 3, 2) = g_cay20(g_i_)
        enddo
        hmt(3, 2) = cay20
C--------
        d4_b = dble(3.)
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 3, 3) = -g_hmt(g_i_, 3, 1) + d4_b * g_cay20(g_i_)
        enddo
        hmt(3, 3) = dble(3.) * cay20 - hmt(3, 1)
C--------
        d4_v = ss3 - 0.6d0 * aa12
        d10_b = -x30 * 0.6d0
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = area * g_y30(g_i_) + y30 * g_area(g_i_) + d4_v 
     ** g_x30(g_i_) + d10_b * g_aa12(g_i_) + x30 * g_ss3(g_i_)
        enddo
        d1_w = d4_v * x30 + area * y30
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 4, 1) = caa12 * g_d1_w(g_i_) + d1_w * g_caa12(g_i_
     *)
        enddo
        hmt(4, 1) = caa12 * d1_w
C--------
        d4_b = dble(-3.)
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 4, 2) = -g_hmt(g_i_, 4, 1) + d4_b * g_cax30(g_i_)
        enddo
        hmt(4, 2) = dble(-3.) * cax30 - hmt(4, 1)
C--------
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 4, 3) = -g_cax30(g_i_)
        enddo
        hmt(4, 3) = -cax30
C--------
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 5, 1) = -g_cax10(g_i_)
        enddo
        hmt(5, 1) = -cax10
C--------
        d4_v = ss1 - 0.6d0 * aa23
        d10_b = -x10 * 0.6d0
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = area * g_y10(g_i_) + y10 * g_area(g_i_) + d4_v 
     ** g_x10(g_i_) + d10_b * g_aa23(g_i_) + x10 * g_ss1(g_i_)
        enddo
        d1_w = d4_v * x10 + area * y10
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 5, 2) = caa23 * g_d1_w(g_i_) + d1_w * g_caa23(g_i_
     *)
        enddo
        hmt(5, 2) = caa23 * d1_w
C--------
        d4_b = dble(-3.)
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 5, 3) = -g_hmt(g_i_, 5, 2) + d4_b * g_cax10(g_i_)
        enddo
        hmt(5, 3) = dble(-3.) * cax10 - hmt(5, 2)
C--------
        d5_v = -ss2 - 0.6d0 * aa31
        d10_b = -x20 * 0.6d0
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = -area * g_y20(g_i_) + (-y20) * g_area(g_i_) + d
     *5_v * g_x20(g_i_) + d10_b * g_aa31(g_i_) + (-x20) * g_ss2(g_i_)
        enddo
        d1_w = d5_v * x20 - area * y20
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 6, 1) = caa31 * g_d1_w(g_i_) + d1_w * g_caa31(g_i_
     *)
        enddo
        hmt(6, 1) = caa31 * d1_w
C--------
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 6, 2) = -g_cax20(g_i_)
        enddo
        hmt(6, 2) = -cax20
C--------
        d4_b = dble(-3.)
        do g_i_ = 1, g_p_
          g_hmt(g_i_, 6, 3) = -g_hmt(g_i_, 6, 1) + d4_b * g_cax20(g_i_)
        enddo
        hmt(6, 3) = dble(-3.) * cax20 - hmt(6, 1)
C--------
        do 99999 j = 1, 3
          d2_b = 2.d0 / dble(9.)
          do g_i_ = 1, g_p_
            g_sum(g_i_) = d2_b * g_hmt(g_i_, 3, j) + d2_b * g_hmt(g_i_, 
     *2, j) + d2_b * g_hmt(g_i_, 1, j)
          enddo
          sum = 2.d0 / dble(9.) * (hmt(1, j) + hmt(2, j) + hmt(3, j))
C--------
          d4_b = -(4.d0 / dble(3.))
          do g_i_ = 1, g_p_
            g_hqt(g_i_, 1, j) = d4_b * g_hmt(g_i_, 1, j) + g_sum(g_i_)
          enddo
          hqt(1, j) = sum - 4.d0 / dble(3.) * hmt(1, j)
C--------
          d4_b = -(4.d0 / dble(3.))
          do g_i_ = 1, g_p_
            g_hqt(g_i_, 2, j) = d4_b * g_hmt(g_i_, 2, j) + g_sum(g_i_)
          enddo
          hqt(2, j) = sum - 4.d0 / dble(3.) * hmt(2, j)
C--------
          d4_b = -(4.d0 / dble(3.))
          do g_i_ = 1, g_p_
            g_hqt(g_i_, 3, j) = d4_b * g_hmt(g_i_, 3, j) + g_sum(g_i_)
          enddo
          hqt(3, j) = sum - 4.d0 / dble(3.) * hmt(3, j)
C--------
          d2_b = 2.d0 / dble(9.)
          do g_i_ = 1, g_p_
            g_sum(g_i_) = d2_b * g_hmt(g_i_, 6, j) + d2_b * g_hmt(g_i_, 
     *5, j) + d2_b * g_hmt(g_i_, 4, j)
          enddo
          sum = 2.d0 / dble(9.) * (hmt(4, j) + hmt(5, j) + hmt(6, j))
C--------
          d4_b = -(4.d0 / dble(3.))
          do g_i_ = 1, g_p_
            g_hqt(g_i_, 4, j) = d4_b * g_hmt(g_i_, 4, j) + g_sum(g_i_)
          enddo
          hqt(4, j) = sum - 4.d0 / dble(3.) * hmt(4, j)
C--------
          d4_b = -(4.d0 / dble(3.))
          do g_i_ = 1, g_p_
            g_hqt(g_i_, 5, j) = d4_b * g_hmt(g_i_, 5, j) + g_sum(g_i_)
          enddo
          hqt(5, j) = sum - 4.d0 / dble(3.) * hmt(5, j)
C--------
          d4_b = -(4.d0 / dble(3.))
          do g_i_ = 1, g_p_
            g_hqt(g_i_, 6, j) = d4_b * g_hmt(g_i_, 6, j) + g_sum(g_i_)
          enddo
          hqt(6, j) = sum - 4.d0 / dble(3.) * hmt(6, j)
C--------
2000      continue
99999   continue
        d2_v = 1.5d0 * f / area2
        d2_b = -d2_v / area2
        do g_i_ = 1, g_p_
          g_kfac(g_i_) = d2_b * g_area2(g_i_)
        enddo
        kfac = d2_v
C--------
        do g_i_ = 1, g_p_
          g_e11(g_i_) = kfac * g_dm(g_i_, 1, 1) + dm(1, 1) * g_kfac(g_i_
     *)
        enddo
        e11 = kfac * dm(1, 1)
C--------
        do g_i_ = 1, g_p_
          g_e22(g_i_) = kfac * g_dm(g_i_, 2, 2) + dm(2, 2) * g_kfac(g_i_
     *)
        enddo
        e22 = kfac * dm(2, 2)
C--------
        do g_i_ = 1, g_p_
          g_e33(g_i_) = kfac * g_dm(g_i_, 3, 3) + dm(3, 3) * g_kfac(g_i_
     *)
        enddo
        e33 = kfac * dm(3, 3)
C--------
        do g_i_ = 1, g_p_
          g_e12(g_i_) = kfac * g_dm(g_i_, 1, 2) + dm(1, 2) * g_kfac(g_i_
     *)
        enddo
        e12 = kfac * dm(1, 2)
C--------
        do g_i_ = 1, g_p_
          g_e13(g_i_) = kfac * g_dm(g_i_, 1, 3) + dm(1, 3) * g_kfac(g_i_
     *)
        enddo
        e13 = kfac * dm(1, 3)
C--------
        do g_i_ = 1, g_p_
          g_e23(g_i_) = kfac * g_dm(g_i_, 2, 3) + dm(2, 3) * g_kfac(g_i_
     *)
        enddo
        e23 = kfac * dm(2, 3)
C--------
        d3_v = y30 * y30
        d1_p = 2.0d0 * y30
        d6_v = dble(2) * e13
        d8_v = d6_v * x30
        d7_b = -y30 * d6_v
        d8_b = -y30 * x30 * dble(2)
        d5_b = -d8_v + e11 * d1_p
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_b * g_x30(g_i_) + d8_b * g_e13(g_i_) + d5_b 
     ** g_y30(g_i_) + d3_v * g_e11(g_i_)
        enddo
        d1_w = e11 * d3_v - d8_v * y30
        d4_v = x30 * x30
        d1_p = 2.0d0 * x30
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e33 * d1_p
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 1, 1) = d7_b * g_x30(g_i_) + d5_b * g_e33(g_i_) + 
     *d2_b * g_d1_w(g_i_)
        enddo
        kqh(1, 1) = dble(2) * (d1_w + e33 * d4_v)
C--------
        d7_v = e13 * x10 - e11 * y10
        d6_b = -y30 * y10
        d7_b = -y30 * e11
        d8_b = y30 * x10
        d9_b = y30 * e13
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y30(g_i_) + d7_b * g_y10(g_i_) + d6_b 
     ** g_e11(g_i_) + d9_b * g_x10(g_i_) + d8_b * g_e13(g_i_)
        enddo
        d1_w = d7_v * y30
        d7_v = e13 * y10 - e33 * x10
        d6_b = -x30 * x10
        d7_b = -x30 * e33
        d8_b = x30 * y10
        d9_b = x30 * e13
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x30(g_i_) + d7_b * g_x10(g_i_) + d6_b 
     ** g_e33(g_i_) + d9_b * g_y10(g_i_) + d8_b * g_e13(g_i_)
        enddo
        d2_w = d7_v * x30
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 1, 2) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(1, 2) = d1_w + d2_w
C--------
        d7_v = e13 * x20 - e11 * y20
        d6_b = -y30 * y20
        d7_b = -y30 * e11
        d8_b = y30 * x20
        d9_b = y30 * e13
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y30(g_i_) + d7_b * g_y20(g_i_) + d6_b 
     ** g_e11(g_i_) + d9_b * g_x20(g_i_) + d8_b * g_e13(g_i_)
        enddo
        d1_w = d7_v * y30
        d7_v = e13 * y20 - e33 * x20
        d6_b = -x30 * x20
        d7_b = -x30 * e33
        d8_b = x30 * y20
        d9_b = x30 * e13
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x30(g_i_) + d7_b * g_x20(g_i_) + d6_b 
     ** g_e33(g_i_) + d9_b * g_y20(g_i_) + d8_b * g_e13(g_i_)
        enddo
        d2_w = d7_v * x30
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 1, 3) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(1, 3) = d1_w + d2_w
C--------
        d3_v = e33 + e12
        d5_v = d3_v * x30
        d4_b = y30 * x30
        d5_b = y30 * d3_v
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d5_v * g_y30(g_i_) + d5_b * g_x30(g_i_) + d4_b 
     ** g_e12(g_i_) + d4_b * g_e33(g_i_)
        enddo
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
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 1, 4) = d7_b * g_x30(g_i_) + d5_b * g_e23(g_i_) + 
     *(-d2_b) * g_d1_w(g_i_) + d12_b * g_y30(g_i_) + d10_b * g_e13(g_i_)
        enddo
        kqh(1, 4) = dble(2) * (e13 * d3_v - d1_w + e23 * d9_v)
C--------
        d7_v = e12 * x10 - e13 * y10
        d6_b = -y30 * y10
        d7_b = -y30 * e13
        d8_b = y30 * x10
        d9_b = y30 * e12
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y30(g_i_) + d7_b * g_y10(g_i_) + d6_b 
     ** g_e13(g_i_) + d9_b * g_x10(g_i_) + d8_b * g_e12(g_i_)
        enddo
        d1_w = d7_v * y30
        d7_v = e33 * y10 - e23 * x10
        d6_b = -x30 * x10
        d7_b = -x30 * e23
        d8_b = x30 * y10
        d9_b = x30 * e33
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x30(g_i_) + d7_b * g_x10(g_i_) + d6_b 
     ** g_e23(g_i_) + d9_b * g_y10(g_i_) + d8_b * g_e33(g_i_)
        enddo
        d2_w = d7_v * x30
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 1, 5) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(1, 5) = d1_w + d2_w
C--------
        d7_v = e12 * x20 - e13 * y20
        d6_b = -y30 * y20
        d7_b = -y30 * e13
        d8_b = y30 * x20
        d9_b = y30 * e12
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y30(g_i_) + d7_b * g_y20(g_i_) + d6_b 
     ** g_e13(g_i_) + d9_b * g_x20(g_i_) + d8_b * g_e12(g_i_)
        enddo
        d1_w = d7_v * y30
        d7_v = e33 * y20 - e23 * x20
        d6_b = -x30 * x20
        d7_b = -x30 * e23
        d8_b = x30 * y20
        d9_b = x30 * e33
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x30(g_i_) + d7_b * g_x20(g_i_) + d6_b 
     ** g_e23(g_i_) + d9_b * g_y20(g_i_) + d8_b * g_e33(g_i_)
        enddo
        d2_w = d7_v * x30
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 1, 6) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(1, 6) = d1_w + d2_w
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 2, 1) = g_kqh(g_i_, 1, 2)
        enddo
        kqh(2, 1) = kqh(1, 2)
C--------
        d3_v = y10 * y10
        d1_p = 2.0d0 * y10
        d6_v = dble(2) * e13
        d8_v = d6_v * x10
        d7_b = -y10 * d6_v
        d8_b = -y10 * x10 * dble(2)
        d5_b = -d8_v + e11 * d1_p
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_b * g_x10(g_i_) + d8_b * g_e13(g_i_) + d5_b 
     ** g_y10(g_i_) + d3_v * g_e11(g_i_)
        enddo
        d1_w = e11 * d3_v - d8_v * y10
        d4_v = x10 * x10
        d1_p = 2.0d0 * x10
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e33 * d1_p
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 2, 2) = d7_b * g_x10(g_i_) + d5_b * g_e33(g_i_) + 
     *d2_b * g_d1_w(g_i_)
        enddo
        kqh(2, 2) = dble(2) * (d1_w + e33 * d4_v)
C--------
        d7_v = e13 * x10 - e11 * y10
        d6_b = -y20 * y10
        d7_b = -y20 * e11
        d8_b = y20 * x10
        d9_b = y20 * e13
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y20(g_i_) + d7_b * g_y10(g_i_) + d6_b 
     ** g_e11(g_i_) + d9_b * g_x10(g_i_) + d8_b * g_e13(g_i_)
        enddo
        d1_w = d7_v * y20
        d7_v = e13 * y10 - e33 * x10
        d6_b = -x20 * x10
        d7_b = -x20 * e33
        d8_b = x20 * y10
        d9_b = x20 * e13
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x20(g_i_) + d7_b * g_x10(g_i_) + d6_b 
     ** g_e33(g_i_) + d9_b * g_y10(g_i_) + d8_b * g_e13(g_i_)
        enddo
        d2_w = d7_v * x20
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 2, 3) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(2, 3) = d1_w + d2_w
C--------
        d7_v = e33 * x10 - e13 * y10
        d6_b = -y30 * y10
        d7_b = -y30 * e13
        d8_b = y30 * x10
        d9_b = y30 * e33
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y30(g_i_) + d7_b * g_y10(g_i_) + d6_b 
     ** g_e13(g_i_) + d9_b * g_x10(g_i_) + d8_b * g_e33(g_i_)
        enddo
        d1_w = d7_v * y30
        d7_v = e12 * y10 - e23 * x10
        d6_b = -x30 * x10
        d7_b = -x30 * e23
        d8_b = x30 * y10
        d9_b = x30 * e12
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x30(g_i_) + d7_b * g_x10(g_i_) + d6_b 
     ** g_e23(g_i_) + d9_b * g_y10(g_i_) + d8_b * g_e12(g_i_)
        enddo
        d2_w = d7_v * x30
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 2, 4) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(2, 4) = d1_w + d2_w
C--------
        d3_v = e33 + e12
        d5_v = d3_v * x10
        d4_b = y10 * x10
        d5_b = y10 * d3_v
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d5_v * g_y10(g_i_) + d5_b * g_x10(g_i_) + d4_b 
     ** g_e12(g_i_) + d4_b * g_e33(g_i_)
        enddo
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
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 2, 5) = d7_b * g_x10(g_i_) + d5_b * g_e23(g_i_) + 
     *(-d2_b) * g_d1_w(g_i_) + d12_b * g_y10(g_i_) + d10_b * g_e13(g_i_)
        enddo
        kqh(2, 5) = dble(2) * (e13 * d3_v - d1_w + e23 * d9_v)
C--------
        d7_v = e33 * x10 - e13 * y10
        d6_b = -y20 * y10
        d7_b = -y20 * e13
        d8_b = y20 * x10
        d9_b = y20 * e33
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y20(g_i_) + d7_b * g_y10(g_i_) + d6_b 
     ** g_e13(g_i_) + d9_b * g_x10(g_i_) + d8_b * g_e33(g_i_)
        enddo
        d1_w = d7_v * y20
        d7_v = e12 * y10 - e23 * x10
        d6_b = -x20 * x10
        d7_b = -x20 * e23
        d8_b = x20 * y10
        d9_b = x20 * e12
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x20(g_i_) + d7_b * g_x10(g_i_) + d6_b 
     ** g_e23(g_i_) + d9_b * g_y10(g_i_) + d8_b * g_e12(g_i_)
        enddo
        d2_w = d7_v * x20
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 2, 6) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(2, 6) = d1_w + d2_w
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 3, 1) = g_kqh(g_i_, 1, 3)
        enddo
        kqh(3, 1) = kqh(1, 3)
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 3, 2) = g_kqh(g_i_, 2, 3)
        enddo
        kqh(3, 2) = kqh(2, 3)
C--------
        d3_v = y20 * y20
        d1_p = 2.0d0 * y20
        d6_v = dble(2) * e13
        d8_v = d6_v * x20
        d7_b = -y20 * d6_v
        d8_b = -y20 * x20 * dble(2)
        d5_b = -d8_v + e11 * d1_p
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_b * g_x20(g_i_) + d8_b * g_e13(g_i_) + d5_b 
     ** g_y20(g_i_) + d3_v * g_e11(g_i_)
        enddo
        d1_w = e11 * d3_v - d8_v * y20
        d4_v = x20 * x20
        d1_p = 2.0d0 * x20
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e33 * d1_p
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 3, 3) = d7_b * g_x20(g_i_) + d5_b * g_e33(g_i_) + 
     *d2_b * g_d1_w(g_i_)
        enddo
        kqh(3, 3) = dble(2) * (d1_w + e33 * d4_v)
C--------
        d7_v = e33 * x20 - e13 * y20
        d6_b = -y30 * y20
        d7_b = -y30 * e13
        d8_b = y30 * x20
        d9_b = y30 * e33
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y30(g_i_) + d7_b * g_y20(g_i_) + d6_b 
     ** g_e13(g_i_) + d9_b * g_x20(g_i_) + d8_b * g_e33(g_i_)
        enddo
        d1_w = d7_v * y30
        d7_v = e12 * y20 - e23 * x20
        d6_b = -x30 * x20
        d7_b = -x30 * e23
        d8_b = x30 * y20
        d9_b = x30 * e12
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x30(g_i_) + d7_b * g_x20(g_i_) + d6_b 
     ** g_e23(g_i_) + d9_b * g_y20(g_i_) + d8_b * g_e12(g_i_)
        enddo
        d2_w = d7_v * x30
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 3, 4) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(3, 4) = d1_w + d2_w
C--------
        d7_v = e12 * x10 - e13 * y10
        d6_b = -y20 * y10
        d7_b = -y20 * e13
        d8_b = y20 * x10
        d9_b = y20 * e12
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y20(g_i_) + d7_b * g_y10(g_i_) + d6_b 
     ** g_e13(g_i_) + d9_b * g_x10(g_i_) + d8_b * g_e12(g_i_)
        enddo
        d1_w = d7_v * y20
        d7_v = e33 * y10 - e23 * x10
        d6_b = -x20 * x10
        d7_b = -x20 * e23
        d8_b = x20 * y10
        d9_b = x20 * e33
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x20(g_i_) + d7_b * g_x10(g_i_) + d6_b 
     ** g_e23(g_i_) + d9_b * g_y10(g_i_) + d8_b * g_e33(g_i_)
        enddo
        d2_w = d7_v * x20
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 3, 5) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(3, 5) = d1_w + d2_w
C--------
        d3_v = e33 + e12
        d5_v = d3_v * x20
        d4_b = y20 * x20
        d5_b = y20 * d3_v
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d5_v * g_y20(g_i_) + d5_b * g_x20(g_i_) + d4_b 
     ** g_e12(g_i_) + d4_b * g_e33(g_i_)
        enddo
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
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 3, 6) = d7_b * g_x20(g_i_) + d5_b * g_e23(g_i_) + 
     *(-d2_b) * g_d1_w(g_i_) + d12_b * g_y20(g_i_) + d10_b * g_e13(g_i_)
        enddo
        kqh(3, 6) = dble(2) * (e13 * d3_v - d1_w + e23 * d9_v)
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 4, 1) = g_kqh(g_i_, 1, 4)
        enddo
        kqh(4, 1) = kqh(1, 4)
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 4, 2) = g_kqh(g_i_, 2, 4)
        enddo
        kqh(4, 2) = kqh(2, 4)
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 4, 3) = g_kqh(g_i_, 3, 4)
        enddo
        kqh(4, 3) = kqh(3, 4)
C--------
        d3_v = y30 * y30
        d1_p = 2.0d0 * y30
        d6_v = dble(2) * e23
        d8_v = d6_v * x30
        d7_b = -y30 * d6_v
        d8_b = -y30 * x30 * dble(2)
        d5_b = -d8_v + e33 * d1_p
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_b * g_x30(g_i_) + d8_b * g_e23(g_i_) + d5_b 
     ** g_y30(g_i_) + d3_v * g_e33(g_i_)
        enddo
        d1_w = e33 * d3_v - d8_v * y30
        d4_v = x30 * x30
        d1_p = 2.0d0 * x30
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e22 * d1_p
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 4, 4) = d7_b * g_x30(g_i_) + d5_b * g_e22(g_i_) + 
     *d2_b * g_d1_w(g_i_)
        enddo
        kqh(4, 4) = dble(2) * (d1_w + e22 * d4_v)
C--------
        d7_v = e23 * x10 - e33 * y10
        d6_b = -y30 * y10
        d7_b = -y30 * e33
        d8_b = y30 * x10
        d9_b = y30 * e23
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y30(g_i_) + d7_b * g_y10(g_i_) + d6_b 
     ** g_e33(g_i_) + d9_b * g_x10(g_i_) + d8_b * g_e23(g_i_)
        enddo
        d1_w = d7_v * y30
        d7_v = e23 * y10 - e22 * x10
        d6_b = -x30 * x10
        d7_b = -x30 * e22
        d8_b = x30 * y10
        d9_b = x30 * e23
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x30(g_i_) + d7_b * g_x10(g_i_) + d6_b 
     ** g_e22(g_i_) + d9_b * g_y10(g_i_) + d8_b * g_e23(g_i_)
        enddo
        d2_w = d7_v * x30
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 4, 5) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(4, 5) = d1_w + d2_w
C--------
        d7_v = e23 * x20 - e33 * y20
        d6_b = -y30 * y20
        d7_b = -y30 * e33
        d8_b = y30 * x20
        d9_b = y30 * e23
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y30(g_i_) + d7_b * g_y20(g_i_) + d6_b 
     ** g_e33(g_i_) + d9_b * g_x20(g_i_) + d8_b * g_e23(g_i_)
        enddo
        d1_w = d7_v * y30
        d7_v = e23 * y20 - e22 * x20
        d6_b = -x30 * x20
        d7_b = -x30 * e22
        d8_b = x30 * y20
        d9_b = x30 * e23
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x30(g_i_) + d7_b * g_x20(g_i_) + d6_b 
     ** g_e22(g_i_) + d9_b * g_y20(g_i_) + d8_b * g_e23(g_i_)
        enddo
        d2_w = d7_v * x30
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 4, 6) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(4, 6) = d1_w + d2_w
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 5, 1) = g_kqh(g_i_, 1, 5)
        enddo
        kqh(5, 1) = kqh(1, 5)
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 5, 2) = g_kqh(g_i_, 2, 5)
        enddo
        kqh(5, 2) = kqh(2, 5)
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 5, 3) = g_kqh(g_i_, 3, 5)
        enddo
        kqh(5, 3) = kqh(3, 5)
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 5, 4) = g_kqh(g_i_, 4, 5)
        enddo
        kqh(5, 4) = kqh(4, 5)
C--------
        d3_v = y10 * y10
        d1_p = 2.0d0 * y10
        d6_v = dble(2) * e23
        d8_v = d6_v * x10
        d7_b = -y10 * d6_v
        d8_b = -y10 * x10 * dble(2)
        d5_b = -d8_v + e33 * d1_p
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_b * g_x10(g_i_) + d8_b * g_e23(g_i_) + d5_b 
     ** g_y10(g_i_) + d3_v * g_e33(g_i_)
        enddo
        d1_w = e33 * d3_v - d8_v * y10
        d4_v = x10 * x10
        d1_p = 2.0d0 * x10
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e22 * d1_p
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 5, 5) = d7_b * g_x10(g_i_) + d5_b * g_e22(g_i_) + 
     *d2_b * g_d1_w(g_i_)
        enddo
        kqh(5, 5) = dble(2) * (d1_w + e22 * d4_v)
C--------
        d7_v = e23 * x10 - e33 * y10
        d6_b = -y20 * y10
        d7_b = -y20 * e33
        d8_b = y20 * x10
        d9_b = y20 * e23
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_v * g_y20(g_i_) + d7_b * g_y10(g_i_) + d6_b 
     ** g_e33(g_i_) + d9_b * g_x10(g_i_) + d8_b * g_e23(g_i_)
        enddo
        d1_w = d7_v * y20
        d7_v = e23 * y10 - e22 * x10
        d6_b = -x20 * x10
        d7_b = -x20 * e22
        d8_b = x20 * y10
        d9_b = x20 * e23
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d7_v * g_x20(g_i_) + d7_b * g_x10(g_i_) + d6_b 
     ** g_e22(g_i_) + d9_b * g_y10(g_i_) + d8_b * g_e23(g_i_)
        enddo
        d2_w = d7_v * x20
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 5, 6) = g_d2_w(g_i_) + g_d1_w(g_i_)
        enddo
        kqh(5, 6) = d1_w + d2_w
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 6, 1) = g_kqh(g_i_, 1, 6)
        enddo
        kqh(6, 1) = kqh(1, 6)
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 6, 2) = g_kqh(g_i_, 2, 6)
        enddo
        kqh(6, 2) = kqh(2, 6)
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 6, 3) = g_kqh(g_i_, 3, 6)
        enddo
        kqh(6, 3) = kqh(3, 6)
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 6, 4) = g_kqh(g_i_, 4, 6)
        enddo
        kqh(6, 4) = kqh(4, 6)
C--------
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 6, 5) = g_kqh(g_i_, 5, 6)
        enddo
        kqh(6, 5) = kqh(5, 6)
C--------
        d3_v = y20 * y20
        d1_p = 2.0d0 * y20
        d6_v = dble(2) * e23
        d8_v = d6_v * x20
        d7_b = -y20 * d6_v
        d8_b = -y20 * x20 * dble(2)
        d5_b = -d8_v + e33 * d1_p
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d7_b * g_x20(g_i_) + d8_b * g_e23(g_i_) + d5_b 
     ** g_y20(g_i_) + d3_v * g_e33(g_i_)
        enddo
        d1_w = e33 * d3_v - d8_v * y20
        d4_v = x20 * x20
        d1_p = 2.0d0 * x20
        d2_b = dble(2)
        d5_b = d2_b * d4_v
        d7_b = d2_b * e22 * d1_p
        do g_i_ = 1, g_p_
          g_kqh(g_i_, 6, 6) = d7_b * g_x20(g_i_) + d5_b * g_e22(g_i_) + 
     *d2_b * g_d1_w(g_i_)
        enddo
        kqh(6, 6) = dble(2) * (d1_w + e22 * d4_v)
C--------
        do g_i_ = 1, g_p_
          g_kth(g_i_, 1, 1) = 0.0d0
        enddo
        kth(1, 1) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_kth(g_i_, 1, 2) = 0.0d0
        enddo
        kth(1, 2) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_kth(g_i_, 2, 2) = 0.0d0
        enddo
        kth(2, 2) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_kth(g_i_, 1, 3) = 0.0d0
        enddo
        kth(1, 3) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_kth(g_i_, 2, 3) = 0.0d0
        enddo
        kth(2, 3) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_kth(g_i_, 3, 3) = 0.0d0
        enddo
        kth(3, 3) = 0.0d0
C--------
        do 99996 j = 1, 3
          do 99998 i = 1, 6
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = kqh(i, 2) * g_hqt(g_i_, 2, j) + hqt(2, j) *
     * g_kqh(g_i_, i, 2) + kqh(i, 1) * g_hqt(g_i_, 1, j) + hqt(1, j) * g
     *_kqh(g_i_, i, 1)
            enddo
            d1_w = kqh(i, 1) * hqt(1, j) + kqh(i, 2) * hqt(2, j)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = kqh(i, 4) * g_hqt(g_i_, 4, j) + hqt(4, j) *
     * g_kqh(g_i_, i, 4) + kqh(i, 3) * g_hqt(g_i_, 3, j) + hqt(3, j) * g
     *_kqh(g_i_, i, 3) + g_d1_w(g_i_)
            enddo
            d2_w = d1_w + kqh(i, 3) * hqt(3, j) + kqh(i, 4) * hqt(4, j)
            do g_i_ = 1, g_p_
              g_w(g_i_, i) = kqh(i, 6) * g_hqt(g_i_, 6, j) + hqt(6, j) *
     * g_kqh(g_i_, i, 6) + kqh(i, 5) * g_hqt(g_i_, 5, j) + hqt(5, j) * g
     *_kqh(g_i_, i, 5) + g_d2_w(g_i_)
            enddo
            w(i) = d2_w + kqh(i, 5) * hqt(5, j) + kqh(i, 6) * hqt(6, j)
C--------
3200        continue
99998     continue
          do 99997 i = 1, j
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = hqt(2, i) * g_w(g_i_, 2) + w(2) * g_hqt(g_i
     *_, 2, i) + hqt(1, i) * g_w(g_i_, 1) + w(1) * g_hqt(g_i_, 1, i) + g
     *_kth(g_i_, i, j)
            enddo
            d1_w = kth(i, j) + hqt(1, i) * w(1) + hqt(2, i) * w(2)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = hqt(4, i) * g_w(g_i_, 4) + w(4) * g_hqt(g_i
     *_, 4, i) + hqt(3, i) * g_w(g_i_, 3) + w(3) * g_hqt(g_i_, 3, i) + g
     *_d1_w(g_i_)
            enddo
            d2_w = d1_w + hqt(3, i) * w(3) + hqt(4, i) * w(4)
            do g_i_ = 1, g_p_
              g_kth(g_i_, i, j) = hqt(6, i) * g_w(g_i_, 6) + w(6) * g_hq
     *t(g_i_, 6, i) + hqt(5, i) * g_w(g_i_, 5) + w(5) * g_hqt(g_i_, 5, i
     *) + g_d2_w(g_i_)
            enddo
            kth(i, j) = d2_w + hqt(5, i) * w(5) + hqt(6, i) * w(6)
C--------
            do g_i_ = 1, g_p_
              g_kth(g_i_, j, i) = g_kth(g_i_, i, j)
            enddo
            kth(j, i) = kth(i, j)
C--------
3300        continue
99997     continue
3500      continue
99996   continue
        do g_i_ = 1, g_p_
          g_s(g_i_, 1) = g_kth(g_i_, 1, 3) + g_kth(g_i_, 1, 2) + g_kth(g
     *_i_, 1, 1)
        enddo
        s(1) = kth(1, 1) + kth(1, 2) + kth(1, 3)
C--------
        do g_i_ = 1, g_p_
          g_s(g_i_, 2) = g_kth(g_i_, 2, 3) + g_kth(g_i_, 2, 2) + g_kth(g
     *_i_, 2, 1)
        enddo
        s(2) = kth(2, 1) + kth(2, 2) + kth(2, 3)
C--------
        do g_i_ = 1, g_p_
          g_s(g_i_, 3) = g_kth(g_i_, 3, 3) + g_kth(g_i_, 3, 2) + g_kth(g
     *_i_, 3, 1)
        enddo
        s(3) = kth(3, 1) + kth(3, 2) + kth(3, 3)
C--------
        d2_v = 0.25d0 / area
        d2_b = -d2_v / area
        do g_i_ = 1, g_p_
          g_ca(g_i_) = d2_b * g_area(g_i_)
        enddo
        ca = d2_v
C--------
        do g_i_ = 1, g_p_
          g_xyij(g_i_, 1) = ca * g_x32(g_i_) + x32 * g_ca(g_i_)
        enddo
        xyij(1) = ca * x32
C--------
        do g_i_ = 1, g_p_
          g_xyij(g_i_, 2) = ca * g_y32(g_i_) + y32 * g_ca(g_i_)
        enddo
        xyij(2) = ca * y32
C--------
        do g_i_ = 1, g_p_
          g_xyij(g_i_, 3) = ca * g_x13(g_i_) + x13 * g_ca(g_i_)
        enddo
        xyij(3) = ca * x13
C--------
        do g_i_ = 1, g_p_
          g_xyij(g_i_, 4) = ca * g_y13(g_i_) + y13 * g_ca(g_i_)
        enddo
        xyij(4) = ca * y13
C--------
        do g_i_ = 1, g_p_
          g_xyij(g_i_, 5) = ca * g_x21(g_i_) + x21 * g_ca(g_i_)
        enddo
        xyij(5) = ca * x21
C--------
        do g_i_ = 1, g_p_
          g_xyij(g_i_, 6) = ca * g_y21(g_i_) + y21 * g_ca(g_i_)
        enddo
        xyij(6) = ca * y21
C--------
        do 99993 j = 1, 9
          l = ls(j)
          do 99995 i = 1, 3
            if (j .le. 6) then
              do g_i_ = 1, g_p_
                g_w(g_i_, i) = s(i) * g_xyij(g_i_, j) + xyij(j) * g_s(g_
     *i_, i)
              enddo
              w(i) = s(i) * xyij(j)
C--------
            else
              do g_i_ = 1, g_p_
                g_w(g_i_, i) = g_kth(g_i_, i, j - 6)
              enddo
              w(i) = kth(i, j - 6)
C--------
            endif
3600        continue
99995     continue
          do g_i_ = 1, g_p_
            g_sum(g_i_) = g_w(g_i_, 3) + g_w(g_i_, 2) + g_w(g_i_, 1)
          enddo
          sum = w(1) + w(2) + w(3)
C--------
          do 99994 i = 1, j
            k = ls(i)
            if (i .le. 6) then
              do g_i_ = 1, g_p_
                g_sm(g_i_, k, l) = sum * g_xyij(g_i_, i) + xyij(i) * g_s
     *um(g_i_) + g_sm(g_i_, k, l)
              enddo
              sm(k, l) = sm(k, l) + sum * xyij(i)
C--------
            else
              do g_i_ = 1, g_p_
                g_sm(g_i_, k, l) = g_w(g_i_, i - 6) + g_sm(g_i_, k, l)
              enddo
              sm(k, l) = sm(k, l) + w(i - 6)
C--------
            endif
            do g_i_ = 1, g_p_
              g_sm(g_i_, l, k) = g_sm(g_i_, k, l)
            enddo
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
