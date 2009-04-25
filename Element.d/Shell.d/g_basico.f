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
C=deck basico basico fortran
C=purpose form basico bending stiffness of c1 triangle
C=author c. a. felippa, may 1984
C=version september 1989
C=equipment machine independent
C=keywords thin plate bending
C=keywords finite element triangle basico stiffness matrix
C=block abstract
C
C     basico forms the basico material element stiffness matrix
C     of a c1 triangle constructed with the generalized
C     ff/ans formulation
C
C=end abstract
C=block usage
C
C     the calling sequence is
C
C       call      basico (x,y,db,f,clr,cqr,ls,sm,m,status)
C
C     where the input arguments are
C
C       x         (3 x 1) array of x coordinates of triangle nodes
C       y         (3 x 1) array of y coordinates of triangle nodes
C       db        (5 x 5) stress resultant-strain matrix
C       f         factor by which stiffness entries will be multiplied.
C       clr,cqr   use clr*llr+cqr*lqr for l
C                 thus clr+cqr  must add up to 1.
C                 llr=linear rotation, lqr=quadratic rotation (kirchhoff)
C       ls        (9 x 1) array of stiffness location pointers
C                 (see output sm).
C       sm        incoming material stiffness array.
C       m         first dimension of sm in calling program.
C
C     the outputs are:
C
C       sm        output stiffness array with basico stiffness
C                 coefficients added in.  the (i,j)-th entry of the
C                 (9 by 9) element bending stiffness is added to
C                 sm(k,l), where k=ls(i) and l=ls(j).
C       status    status character variable.  blank if no error
C                 detected.
C
C=end usage
C=block fortran
      subroutine g_basico(g_p_, x, g_x, ldg_x, y, g_y, ldg_y, db, g_db, 
     *ldg_db, f, clr, cqr, ls, sm, g_sm, ldg_sm, m, status)
C
C                   a r g u m e n t s
C
        real*8 x(3), y(3), db(3, 3), f, sm(m, m)
        real*8 clr, cqr
        integer m, ls(9)
        character*(*) status
C
C                   t y p e   &   d i m e n s i o n
C
        real*8 llr(9, 3), lqr(9, 3), l(9, 3)
        real*8 db11, db12, db13, db22, db23, db33
        real*8 x0, y0, cab, a1, a2, a3, b1, b2, b3
        real*8 x21, x32, x13, y21, y32, y13
        real*8 x12, x23, x31, y12, y23, y31
        real*8 xl12, xl23, xl31, c12, c23, c31, s12, s23, s31
        real*8 cc12, cc23, cc31, ss12, ss23, ss31
        real*8 cs12, cs23, cs31, s1, s2, s3
        real*8 area2, c
        integer i, j, ii, jj
C
C                   l o g i c
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_x, ldg_y, ldg_db, ldg_sm
        double precision d9_b, d3_b, d8_b, d4_b, d1_p, d1_w, d7_b, d6_b,
     * d2_v, d3_v
        double precision d5_b, d2_b, g_x21(g_pmax_), g_x(ldg_x, 3), g_x1
     *2(g_pmax_), g_x32(g_pmax_), g_x23(g_pmax_), g_x13(g_pmax_), g_x31(
     *g_pmax_), g_y21(g_pmax_)
        double precision g_y(ldg_y, 3), g_y12(g_pmax_), g_y32(g_pmax_), 
     *g_y23(g_pmax_), g_y13(g_pmax_), g_y31(g_pmax_), g_area2(g_pmax_), 
     *g_d1_w(g_pmax_), g_xl12(g_pmax_), g_xl23(g_pmax_)
        double precision g_xl31(g_pmax_), g_llr(g_pmax_, 9, 3), g_lqr(g_
     *pmax_, 9, 3), g_c12(g_pmax_), g_s12(g_pmax_), g_c23(g_pmax_), g_s2
     *3(g_pmax_), g_c31(g_pmax_), g_s31(g_pmax_), g_cc12(g_pmax_)
        double precision g_cc23(g_pmax_), g_cc31(g_pmax_), g_ss12(g_pmax
     *_), g_ss23(g_pmax_), g_ss31(g_pmax_), g_cs12(g_pmax_), g_cs23(g_pm
     *ax_), g_cs31(g_pmax_), g_l(g_pmax_, 9, 3), g_c(g_pmax_)
        double precision g_db11(g_pmax_), g_db(ldg_db, 3, 3), g_db22(g_p
     *max_), g_db33(g_pmax_), g_db12(g_pmax_), g_db13(g_pmax_), g_db23(g
     *_pmax_), g_s1(g_pmax_), g_s2(g_pmax_), g_s3(g_pmax_)
        double precision g_sm(ldg_sm, m, m)
        save g_db13, g_db23, g_s1, g_s2, g_s3
        save g_ss31, g_cs12, g_cs23, g_cs31, g_l, g_c, g_db11, g_db22, g
     *_db33, g_db12
        save g_s12, g_c23, g_s23, g_c31, g_s31, g_cc12, g_cc23, g_cc31, 
     *g_ss12, g_ss23
        save g_y13, g_y31, g_area2, g_d1_w, g_xl12, g_xl23, g_xl31, g_ll
     *r, g_lqr, g_c12
        save g_x21, g_x12, g_x32, g_x23, g_x13, g_x31, g_y21, g_y12, g_y
     *32, g_y23
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'basico','g_basico.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        status = ' '
        do g_i_ = 1, g_p_
          g_x21(g_i_) = -g_x(g_i_, 1) + g_x(g_i_, 2)
        enddo
        x21 = x(2) - x(1)
C--------
        do g_i_ = 1, g_p_
          g_x12(g_i_) = -g_x21(g_i_)
        enddo
        x12 = -x21
C--------
        do g_i_ = 1, g_p_
          g_x32(g_i_) = -g_x(g_i_, 2) + g_x(g_i_, 3)
        enddo
        x32 = x(3) - x(2)
C--------
        do g_i_ = 1, g_p_
          g_x23(g_i_) = -g_x32(g_i_)
        enddo
        x23 = -x32
C--------
        do g_i_ = 1, g_p_
          g_x13(g_i_) = -g_x(g_i_, 3) + g_x(g_i_, 1)
        enddo
        x13 = x(1) - x(3)
C--------
        do g_i_ = 1, g_p_
          g_x31(g_i_) = -g_x13(g_i_)
        enddo
        x31 = -x13
C--------
        do g_i_ = 1, g_p_
          g_y21(g_i_) = -g_y(g_i_, 1) + g_y(g_i_, 2)
        enddo
        y21 = y(2) - y(1)
C--------
        do g_i_ = 1, g_p_
          g_y12(g_i_) = -g_y21(g_i_)
        enddo
        y12 = -y21
C--------
        do g_i_ = 1, g_p_
          g_y32(g_i_) = -g_y(g_i_, 2) + g_y(g_i_, 3)
        enddo
        y32 = y(3) - y(2)
C--------
        do g_i_ = 1, g_p_
          g_y23(g_i_) = -g_y32(g_i_)
        enddo
        y23 = -y32
C--------
        do g_i_ = 1, g_p_
          g_y13(g_i_) = -g_y(g_i_, 3) + g_y(g_i_, 1)
        enddo
        y13 = y(1) - y(3)
C--------
        do g_i_ = 1, g_p_
          g_y31(g_i_) = -g_y13(g_i_)
        enddo
        y31 = -y13
C--------
        do g_i_ = 1, g_p_
          g_area2(g_i_) = -x21 * g_y13(g_i_) + (-y13) * g_x21(g_i_) + y2
     *1 * g_x13(g_i_) + x13 * g_y21(g_i_)
        enddo
        area2 = y21 * x13 - x21 * y13
C--------
        if (area2 .le. 0.0) then
          status = 'basico: negative area'
          if (area2 .eq. 0.0) then
            status = 'basico: zero area'
          endif
          return
        endif
        x0 = (x(1) + x(2) + x(3)) / 3.
        y0 = (y(1) + y(2) + y(3)) / 3.
        cab = 3.0 / area2
        a1 = -cab * (y(3) - y0)
        a2 = -cab * (y(1) - y0)
        a3 = -cab * (y(2) - y0)
        b1 = cab * (x(3) - x0)
        b2 = cab * (x(1) - x0)
        b3 = cab * (x(2) - x0)
C
        d4_b = y12 + y12
        d5_b = x12 + x12
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y12(g_i_) + d5_b * g_x12(g_i_)
        enddo
        d1_w = x12 * x12 + y12 * y12
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +222)
        endif
        do g_i_ = 1, g_p_
          g_xl12(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        xl12 = d2_v
C--------
        d4_b = y23 + y23
        d5_b = x23 + x23
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y23(g_i_) + d5_b * g_x23(g_i_)
        enddo
        d1_w = x23 * x23 + y23 * y23
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +241)
        endif
        do g_i_ = 1, g_p_
          g_xl23(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        xl23 = d2_v
C--------
        d4_b = y31 + y31
        d5_b = x31 + x31
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y31(g_i_) + d5_b * g_x31(g_i_)
        enddo
        d1_w = x31 * x31 + y31 * y31
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +260)
        endif
        do g_i_ = 1, g_p_
          g_xl31(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        xl31 = d2_v
C--------
C
        do 99998 j = 1, 3
          do 99999 i = 1, 9
            do g_i_ = 1, g_p_
              g_llr(g_i_, i, j) = 0.0d0
            enddo
            llr(i, j) = 0.0d0
C--------
            do g_i_ = 1, g_p_
              g_lqr(g_i_, i, j) = 0.0d0
            enddo
            lqr(i, j) = 0.0d0
C--------
1200        continue
99999     continue
99998   continue
C
        if (clr .ne. 0.0) then
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 3, 1) = d2_b * g_y32(g_i_)
          enddo
          llr(3, 1) = y32 * dble(.5)
C--------
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 6, 1) = d2_b * g_y13(g_i_)
          enddo
          llr(6, 1) = y13 * dble(.5)
C--------
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 9, 1) = d2_b * g_y21(g_i_)
          enddo
          llr(9, 1) = y21 * dble(.5)
C--------
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 2, 2) = d2_b * g_x32(g_i_)
          enddo
          llr(2, 2) = x32 * dble(.5)
C--------
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 5, 2) = d2_b * g_x13(g_i_)
          enddo
          llr(5, 2) = x13 * dble(.5)
C--------
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 8, 2) = d2_b * g_x21(g_i_)
          enddo
          llr(8, 2) = x21 * dble(.5)
C--------
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 2, 3) = -d2_b * g_y32(g_i_)
          enddo
          llr(2, 3) = -y32 * dble(.5)
C--------
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 3, 3) = -d2_b * g_x32(g_i_)
          enddo
          llr(3, 3) = -x32 * dble(.5)
C--------
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 5, 3) = -d2_b * g_y13(g_i_)
          enddo
          llr(5, 3) = -y13 * dble(.5)
C--------
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 6, 3) = -d2_b * g_x13(g_i_)
          enddo
          llr(6, 3) = -x13 * dble(.5)
C--------
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 8, 3) = -d2_b * g_y21(g_i_)
          enddo
          llr(8, 3) = -y21 * dble(.5)
C--------
          d2_b = dble(.5)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 9, 3) = -d2_b * g_x21(g_i_)
          enddo
          llr(9, 3) = -x21 * dble(.5)
C--------
        endif
C
        if (cqr .ne. 0.0) then
          d3_v = y21 / xl12
          d2_b = 1.0d0 / xl12
          d3_b = -d3_v / xl12
          do g_i_ = 1, g_p_
            g_c12(g_i_) = d3_b * g_xl12(g_i_) + d2_b * g_y21(g_i_)
          enddo
          c12 = d3_v
C--------
          d3_v = x12 / xl12
          d2_b = 1.0d0 / xl12
          d3_b = -d3_v / xl12
          do g_i_ = 1, g_p_
            g_s12(g_i_) = d3_b * g_xl12(g_i_) + d2_b * g_x12(g_i_)
          enddo
          s12 = d3_v
C--------
          d3_v = y32 / xl23
          d2_b = 1.0d0 / xl23
          d3_b = -d3_v / xl23
          do g_i_ = 1, g_p_
            g_c23(g_i_) = d3_b * g_xl23(g_i_) + d2_b * g_y32(g_i_)
          enddo
          c23 = d3_v
C--------
          d3_v = x23 / xl23
          d2_b = 1.0d0 / xl23
          d3_b = -d3_v / xl23
          do g_i_ = 1, g_p_
            g_s23(g_i_) = d3_b * g_xl23(g_i_) + d2_b * g_x23(g_i_)
          enddo
          s23 = d3_v
C--------
          d3_v = y13 / xl31
          d2_b = 1.0d0 / xl31
          d3_b = -d3_v / xl31
          do g_i_ = 1, g_p_
            g_c31(g_i_) = d3_b * g_xl31(g_i_) + d2_b * g_y13(g_i_)
          enddo
          c31 = d3_v
C--------
          d3_v = x31 / xl31
          d2_b = 1.0d0 / xl31
          d3_b = -d3_v / xl31
          do g_i_ = 1, g_p_
            g_s31(g_i_) = d3_b * g_xl31(g_i_) + d2_b * g_x31(g_i_)
          enddo
          s31 = d3_v
C--------
          d2_b = c12 + c12
          do g_i_ = 1, g_p_
            g_cc12(g_i_) = d2_b * g_c12(g_i_)
          enddo
          cc12 = c12 * c12
C--------
          d2_b = c23 + c23
          do g_i_ = 1, g_p_
            g_cc23(g_i_) = d2_b * g_c23(g_i_)
          enddo
          cc23 = c23 * c23
C--------
          d2_b = c31 + c31
          do g_i_ = 1, g_p_
            g_cc31(g_i_) = d2_b * g_c31(g_i_)
          enddo
          cc31 = c31 * c31
C--------
          d2_b = s12 + s12
          do g_i_ = 1, g_p_
            g_ss12(g_i_) = d2_b * g_s12(g_i_)
          enddo
          ss12 = s12 * s12
C--------
          d2_b = s23 + s23
          do g_i_ = 1, g_p_
            g_ss23(g_i_) = d2_b * g_s23(g_i_)
          enddo
          ss23 = s23 * s23
C--------
          d2_b = s31 + s31
          do g_i_ = 1, g_p_
            g_ss31(g_i_) = d2_b * g_s31(g_i_)
          enddo
          ss31 = s31 * s31
C--------
          do g_i_ = 1, g_p_
            g_cs12(g_i_) = c12 * g_s12(g_i_) + s12 * g_c12(g_i_)
          enddo
          cs12 = c12 * s12
C--------
          do g_i_ = 1, g_p_
            g_cs23(g_i_) = c23 * g_s23(g_i_) + s23 * g_c23(g_i_)
          enddo
          cs23 = c23 * s23
C--------
          do g_i_ = 1, g_p_
            g_cs31(g_i_) = c31 * g_s31(g_i_) + s31 * g_c31(g_i_)
          enddo
          cs31 = c31 * s31
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 1, 1) = -g_cs31(g_i_) + g_cs12(g_i_)
          enddo
          lqr(1, 1) = cs12 - cs31
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 1, 2) = -g_lqr(g_i_, 1, 1)
          enddo
          lqr(1, 2) = -lqr(1, 1)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 1, 3) = g_ss12(g_i_) + (-g_cc12(g_i_)) + (-g_ss3
     *1(g_i_)) + g_cc31(g_i_)
          enddo
          lqr(1, 3) = cc31 - ss31 - (cc12 - ss12)
C--------
          d2_b = dble(.5)
          d5_b = d2_b * x31
          d6_b = d2_b * cc31
          d7_b = d2_b * x12
          d8_b = d2_b * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 2, 1) = d6_b * g_x31(g_i_) + d5_b * g_cc31(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_cc12(g_i_)
          enddo
          lqr(2, 1) = (cc12 * x12 + cc31 * x31) * dble(.5)
C--------
          d2_b = dble(.5)
          d5_b = d2_b * x31
          d6_b = d2_b * ss31
          d7_b = d2_b * x12
          d8_b = d2_b * ss12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 2, 2) = d6_b * g_x31(g_i_) + d5_b * g_ss31(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_ss12(g_i_)
          enddo
          lqr(2, 2) = (ss12 * x12 + ss31 * x31) * dble(.5)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 2, 3) = ss31 * g_y13(g_i_) + y13 * g_ss31(g_i_) 
     *+ ss12 * g_y21(g_i_) + y21 * g_ss12(g_i_)
          enddo
          lqr(2, 3) = ss12 * y21 + ss31 * y13
C--------
          d2_b = dble(.5)
          d6_b = -d2_b * y13
          d7_b = -d2_b * cc31
          d8_b = -d2_b * y21
          d9_b = -d2_b * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 3, 1) = d7_b * g_y13(g_i_) + d6_b * g_cc31(g_i_)
     * + d9_b * g_y21(g_i_) + d8_b * g_cc12(g_i_)
          enddo
          lqr(3, 1) = -(cc12 * y21 + cc31 * y13) * dble(.5)
C--------
          d2_b = dble(-0.5)
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 3, 2) = d2_b * g_lqr(g_i_, 2, 3)
          enddo
          lqr(3, 2) = dble(-0.5) * lqr(2, 3)
C--------
          d2_b = dble(-2.)
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 3, 3) = d2_b * g_lqr(g_i_, 2, 1)
          enddo
          lqr(3, 3) = dble(-2.) * lqr(2, 1)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 4, 1) = -g_cs12(g_i_) + g_cs23(g_i_)
          enddo
          lqr(4, 1) = cs23 - cs12
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 4, 2) = -g_lqr(g_i_, 4, 1)
          enddo
          lqr(4, 2) = -lqr(4, 1)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 4, 3) = g_ss23(g_i_) + (-g_cc23(g_i_)) + (-g_ss1
     *2(g_i_)) + g_cc12(g_i_)
          enddo
          lqr(4, 3) = cc12 - ss12 - (cc23 - ss23)
C--------
          d2_b = dble(.5)
          d5_b = d2_b * x23
          d6_b = d2_b * cc23
          d7_b = d2_b * x12
          d8_b = d2_b * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 5, 1) = d6_b * g_x23(g_i_) + d5_b * g_cc23(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_cc12(g_i_)
          enddo
          lqr(5, 1) = (cc12 * x12 + cc23 * x23) * dble(.5)
C--------
          d2_b = dble(.5)
          d5_b = d2_b * x23
          d6_b = d2_b * ss23
          d7_b = d2_b * x12
          d8_b = d2_b * ss12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 5, 2) = d6_b * g_x23(g_i_) + d5_b * g_ss23(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_ss12(g_i_)
          enddo
          lqr(5, 2) = (ss12 * x12 + ss23 * x23) * dble(.5)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 5, 3) = ss23 * g_y32(g_i_) + y32 * g_ss23(g_i_) 
     *+ ss12 * g_y21(g_i_) + y21 * g_ss12(g_i_)
          enddo
          lqr(5, 3) = ss12 * y21 + ss23 * y32
C--------
          d2_b = dble(.5)
          d6_b = -d2_b * y32
          d7_b = -d2_b * cc23
          d8_b = -d2_b * y21
          d9_b = -d2_b * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 6, 1) = d7_b * g_y32(g_i_) + d6_b * g_cc23(g_i_)
     * + d9_b * g_y21(g_i_) + d8_b * g_cc12(g_i_)
          enddo
          lqr(6, 1) = -(cc12 * y21 + cc23 * y32) * dble(.5)
C--------
          d2_b = dble(-0.5)
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 6, 2) = d2_b * g_lqr(g_i_, 5, 3)
          enddo
          lqr(6, 2) = dble(-0.5) * lqr(5, 3)
C--------
          d2_b = dble(-2.)
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 6, 3) = d2_b * g_lqr(g_i_, 5, 1)
          enddo
          lqr(6, 3) = dble(-2.) * lqr(5, 1)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 7, 1) = -g_cs23(g_i_) + g_cs31(g_i_)
          enddo
          lqr(7, 1) = cs31 - cs23
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 7, 2) = -g_lqr(g_i_, 7, 1)
          enddo
          lqr(7, 2) = -lqr(7, 1)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 7, 3) = g_ss31(g_i_) + (-g_cc31(g_i_)) + (-g_ss2
     *3(g_i_)) + g_cc23(g_i_)
          enddo
          lqr(7, 3) = cc23 - ss23 - (cc31 - ss31)
C--------
          d2_b = dble(.5)
          d5_b = d2_b * x31
          d6_b = d2_b * cc31
          d7_b = d2_b * x23
          d8_b = d2_b * cc23
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 8, 1) = d6_b * g_x31(g_i_) + d5_b * g_cc31(g_i_)
     * + d8_b * g_x23(g_i_) + d7_b * g_cc23(g_i_)
          enddo
          lqr(8, 1) = (cc23 * x23 + cc31 * x31) * dble(.5)
C--------
          d2_b = dble(.5)
          d5_b = d2_b * x31
          d6_b = d2_b * ss31
          d7_b = d2_b * x23
          d8_b = d2_b * ss23
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 8, 2) = d6_b * g_x31(g_i_) + d5_b * g_ss31(g_i_)
     * + d8_b * g_x23(g_i_) + d7_b * g_ss23(g_i_)
          enddo
          lqr(8, 2) = (ss23 * x23 + ss31 * x31) * dble(.5)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 8, 3) = ss31 * g_y13(g_i_) + y13 * g_ss31(g_i_) 
     *+ ss23 * g_y32(g_i_) + y32 * g_ss23(g_i_)
          enddo
          lqr(8, 3) = ss23 * y32 + ss31 * y13
C--------
          d2_b = dble(.5)
          d6_b = -d2_b * y13
          d7_b = -d2_b * cc31
          d8_b = -d2_b * y32
          d9_b = -d2_b * cc23
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 9, 1) = d7_b * g_y13(g_i_) + d6_b * g_cc31(g_i_)
     * + d9_b * g_y32(g_i_) + d8_b * g_cc23(g_i_)
          enddo
          lqr(9, 1) = -(cc23 * y32 + cc31 * y13) * dble(.5)
C--------
          d2_b = dble(-0.5)
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 9, 2) = d2_b * g_lqr(g_i_, 8, 3)
          enddo
          lqr(9, 2) = dble(-0.5) * lqr(8, 3)
C--------
          d2_b = dble(-2.)
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 9, 3) = d2_b * g_lqr(g_i_, 8, 1)
          enddo
          lqr(9, 3) = dble(-2.) * lqr(8, 1)
C--------
        endif
C
C     write(*,*) ((lqr(i,j),j=1,3),i=1,9)
C
        do 99997 j = 1, 9
          do g_i_ = 1, g_p_
            g_l(g_i_, j, 1) = cqr * g_lqr(g_i_, j, 1) + clr * g_llr(g_i_
     *, j, 1)
          enddo
          l(j, 1) = clr * llr(j, 1) + cqr * lqr(j, 1)
C--------
          do g_i_ = 1, g_p_
            g_l(g_i_, j, 2) = cqr * g_lqr(g_i_, j, 2) + clr * g_llr(g_i_
     *, j, 2)
          enddo
          l(j, 2) = clr * llr(j, 2) + cqr * lqr(j, 2)
C--------
          do g_i_ = 1, g_p_
            g_l(g_i_, j, 3) = cqr * g_lqr(g_i_, j, 3) + clr * g_llr(g_i_
     *, j, 3)
          enddo
          l(j, 3) = clr * llr(j, 3) + cqr * lqr(j, 3)
C--------
1600      continue
99997   continue
C
        d2_v = dble(2.0) * f / area2
        d2_b = -d2_v / area2
        do g_i_ = 1, g_p_
          g_c(g_i_) = d2_b * g_area2(g_i_)
        enddo
        c = d2_v
C--------
        do g_i_ = 1, g_p_
          g_db11(g_i_) = c * g_db(g_i_, 1, 1) + db(1, 1) * g_c(g_i_)
        enddo
        db11 = c * db(1, 1)
C--------
        do g_i_ = 1, g_p_
          g_db22(g_i_) = c * g_db(g_i_, 2, 2) + db(2, 2) * g_c(g_i_)
        enddo
        db22 = c * db(2, 2)
C--------
        do g_i_ = 1, g_p_
          g_db33(g_i_) = c * g_db(g_i_, 3, 3) + db(3, 3) * g_c(g_i_)
        enddo
        db33 = c * db(3, 3)
C--------
        do g_i_ = 1, g_p_
          g_db12(g_i_) = c * g_db(g_i_, 1, 2) + db(1, 2) * g_c(g_i_)
        enddo
        db12 = c * db(1, 2)
C--------
        do g_i_ = 1, g_p_
          g_db13(g_i_) = c * g_db(g_i_, 1, 3) + db(1, 3) * g_c(g_i_)
        enddo
        db13 = c * db(1, 3)
C--------
        do g_i_ = 1, g_p_
          g_db23(g_i_) = c * g_db(g_i_, 2, 3) + db(2, 3) * g_c(g_i_)
        enddo
        db23 = c * db(2, 3)
C--------
        do 99995 j = 1, 9
          jj = ls(j)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db12 * g_l(g_i_, j, 2) + l(j, 2) * g_db12(g_i
     *_) + db11 * g_l(g_i_, j, 1) + l(j, 1) * g_db11(g_i_)
          enddo
          d1_w = db11 * l(j, 1) + db12 * l(j, 2)
          do g_i_ = 1, g_p_
            g_s1(g_i_) = db13 * g_l(g_i_, j, 3) + l(j, 3) * g_db13(g_i_)
     * + g_d1_w(g_i_)
          enddo
          s1 = d1_w + db13 * l(j, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db22 * g_l(g_i_, j, 2) + l(j, 2) * g_db22(g_i
     *_) + db12 * g_l(g_i_, j, 1) + l(j, 1) * g_db12(g_i_)
          enddo
          d1_w = db12 * l(j, 1) + db22 * l(j, 2)
          do g_i_ = 1, g_p_
            g_s2(g_i_) = db23 * g_l(g_i_, j, 3) + l(j, 3) * g_db23(g_i_)
     * + g_d1_w(g_i_)
          enddo
          s2 = d1_w + db23 * l(j, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db23 * g_l(g_i_, j, 2) + l(j, 2) * g_db23(g_i
     *_) + db13 * g_l(g_i_, j, 1) + l(j, 1) * g_db13(g_i_)
          enddo
          d1_w = db13 * l(j, 1) + db23 * l(j, 2)
          do g_i_ = 1, g_p_
            g_s3(g_i_) = db33 * g_l(g_i_, j, 3) + l(j, 3) * g_db33(g_i_)
     * + g_d1_w(g_i_)
          enddo
          s3 = d1_w + db33 * l(j, 3)
C--------
          do 99996 i = 1, j
            ii = ls(i)
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = s2 * g_l(g_i_, i, 2) + l(i, 2) * g_s2(g_i_)
     * + s1 * g_l(g_i_, i, 1) + l(i, 1) * g_s1(g_i_)
            enddo
            d1_w = s1 * l(i, 1) + s2 * l(i, 2)
            do g_i_ = 1, g_p_
              g_sm(g_i_, jj, ii) = s3 * g_l(g_i_, i, 3) + l(i, 3) * g_s3
     *(g_i_) + g_d1_w(g_i_) + g_sm(g_i_, jj, ii)
            enddo
            sm(jj, ii) = sm(jj, ii) + (d1_w + s3 * l(i, 3))
C--------
            do g_i_ = 1, g_p_
              g_sm(g_i_, ii, jj) = g_sm(g_i_, jj, ii)
            enddo
            sm(ii, jj) = sm(jj, ii)
C--------
3500        continue
99996     continue
4000      continue
99995   continue
        return
      end
C=end fortran
