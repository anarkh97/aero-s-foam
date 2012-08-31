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
C=deck smcbh smcbh fortran
C=purpose form higher-order bending stiffness obtained
C          from curvatures over the sides
C=author c. militello, may 1989
C=version may 1989
C=equipment machine independent
C=keywords thin plate bending
C=keywords finite element triangle higher order stiffness matrix
C=block abstract
C
C     smcbh forms the higher order material stiffness matrix of a
C     9-dof thin-plate-bending triangle obtained by using linear
C     curvatures over the sides
C
C=end abstract
C=block usage
C
C     the calling sequence is
C
C       call      smcbh (x, y, db, f, ,ls, sm, m, status)
C
C     where the input arguments are
C
C       x         (3 x 1) array of x coordinates of triangle nodes
C       y         (3 x 1) array of y coordinates of triangle nodes
C       db        (3 x 3) moment-curvature matrix.
C       f         factor by which stiffness entries will be multiplied.
C       ls        (9 x 1) array of stiffness location pointers
C                 (see output sm).
C       sm        incoming material stiffness array.
C       m         first dimension of sm in calling program.
C
C     the outputs are:
C
C       sm        output stiffness array with higher order stiffness
C                 coefficients added in.  the (i,j)-th entry of the
C                 (9 by 9) element bending stiffness is added to
C                 sm(k,l), where k=ls(i) and l=ls(j).
C       status    status character variable.  blank if no error
C                 detected.
C
C=end usage
C=block fortran
      subroutine g_smcbh(g_p_, x, g_x, ldg_x, y, g_y, ldg_y, db, g_db, l
     *dg_db, f, ls, sm, g_sm, ldg_sm, m, status)
C
C                   t y p e   &   d i m e n s i o n
C
C=block vax
C     implicit      none
C=end vax
        integer m
        character*(*) status
        double precision x(3), y(3), db(3, 3), pg(3, 3), rsd(6, 6)
        double precision sm(m, m), sq(3, 3), sds(3, 3), q(6, 9)
        double precision rm(3, 6), l1, l2, l3
        double precision f, x0, y0, x1, x2, x3, y1, y2, y3
        double precision x21, x32, x13, y21, y32, y13, area, area2
        double precision l21, l32, l13, bl2, al2, bl3, al3, bl1, al1
        double precision cc, d11, d22, d33, d12, d13, d23
        double precision s1, s2, s3, s4, s5, s6
        integer ls(9)
        integer i, j, k, l
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_x, ldg_y, ldg_db, ldg_sm
        double precision d3_w, d14_b, d14_v, d4_p, d3_p, d13_b, d12_b, d
     *10_b, d2_v, d3_v
        double precision d4_v, d5_v, d6_v, d13_v, d2_b, d3_b, d4_b, d5_b
     *, d6_b, d7_v
        double precision d7_b, d1_w, d1_p, d2_p, d2_w, d8_v, d8_b, d9_v,
     * d10_v, d9_b
        double precision g_rm(g_pmax_, 3, 6), g_q(g_pmax_, 6, 9), g_rsd(
     *g_pmax_, 6, 6), g_sds(g_pmax_, 3, 3), g_x0(g_pmax_), g_x(ldg_x, 3)
     *, g_y0(g_pmax_), g_y(ldg_y, 3), g_x1(g_pmax_), g_x2(g_pmax_)
        double precision g_x3(g_pmax_), g_y1(g_pmax_), g_y2(g_pmax_), g_
     *y3(g_pmax_), g_x21(g_pmax_), g_x32(g_pmax_), g_x13(g_pmax_), g_y21
     *(g_pmax_), g_y32(g_pmax_), g_y13(g_pmax_)
        double precision g_area2(g_pmax_), g_d1_w(g_pmax_), g_l21(g_pmax
     *_), g_l32(g_pmax_), g_l13(g_pmax_), g_d2_w(g_pmax_), g_bl2(g_pmax_
     *), g_al2(g_pmax_), g_bl3(g_pmax_), g_al3(g_pmax_)
        double precision g_bl1(g_pmax_), g_al1(g_pmax_), g_cc(g_pmax_), 
     *g_sq(g_pmax_, 3, 3), g_d11(g_pmax_), g_db(ldg_db, 3, 3), g_d22(g_p
     *max_), g_d33(g_pmax_), g_d12(g_pmax_), g_d13(g_pmax_)
        double precision g_d23(g_pmax_), g_area(g_pmax_), g_s1(g_pmax_),
     * g_s2(g_pmax_), g_s3(g_pmax_), g_s4(g_pmax_), g_s5(g_pmax_), g_s6(
     *g_pmax_), g_d3_w(g_pmax_), g_sm(ldg_sm, m, m)
        save g_s2, g_s3, g_s4, g_s5, g_s6, g_d3_w
        save g_cc, g_sq, g_d11, g_d22, g_d33, g_d12, g_d13, g_d23, g_are
     *a, g_s1
        save g_l21, g_l32, g_l13, g_d2_w, g_bl2, g_al2, g_bl3, g_al3, g_
     *bl1, g_al1
        save g_y2, g_y3, g_x21, g_x32, g_x13, g_y21, g_y32, g_y13, g_are
     *a2, g_d1_w
        save g_rm, g_q, g_rsd, g_sds, g_x0, g_y0, g_x1, g_x2, g_x3, g_y1
        intrinsic dble
        data pg /0.d0, 0.5d0, 0.5d0, 0.5d0, 0.d0, 0.5d0, 0.5d0, 0.5d0, 0
     *.d0/
C
C                   l o g i c
C
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'smcbh','g_smcbh.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        status = ' '
C     cleaning
C
        do 99998 i = 1, 3
          do 99999 j = 1, 6
            do g_i_ = 1, g_p_
              g_rm(g_i_, i, j) = 0.0d0
            enddo
            rm(i, j) = 0.0d0
C--------
100         continue
99999     continue
99998   continue
        do 99996 i = 1, 6
          do 99997 j = 1, 9
            do g_i_ = 1, g_p_
              g_q(g_i_, i, j) = 0.0d0
            enddo
            q(i, j) = 0.0d0
C--------
200         continue
99997     continue
99996   continue
        do 99994 i = 1, 6
          do 99995 j = 1, 6
300         do g_i_ = 1, g_p_
              g_rsd(g_i_, i, j) = 0.0d0
            enddo
            rsd(i, j) = 0.0d0
C--------
99995     continue
99994   continue
        do 99992 i = 1, 3
          do 99993 j = 1, 3
400         do g_i_ = 1, g_p_
              g_sds(g_i_, i, j) = 0.0d0
            enddo
            sds(i, j) = 0.0d0
C--------
99993     continue
99992   continue
C     coordinates
        d2_b = 1.0d0 / dble(3.)
        do g_i_ = 1, g_p_
          g_x0(g_i_) = d2_b * g_x(g_i_, 3) + d2_b * g_x(g_i_, 2) + d2_b 
     ** g_x(g_i_, 1)
        enddo
        x0 = (x(1) + x(2) + x(3)) / dble(3.)
C--------
        d2_b = 1.0d0 / dble(3.)
        do g_i_ = 1, g_p_
          g_y0(g_i_) = d2_b * g_y(g_i_, 3) + d2_b * g_y(g_i_, 2) + d2_b 
     ** g_y(g_i_, 1)
        enddo
        y0 = (y(1) + y(2) + y(3)) / dble(3.)
C--------
        do g_i_ = 1, g_p_
          g_x1(g_i_) = -g_x0(g_i_) + g_x(g_i_, 1)
        enddo
        x1 = x(1) - x0
C--------
        do g_i_ = 1, g_p_
          g_x2(g_i_) = -g_x0(g_i_) + g_x(g_i_, 2)
        enddo
        x2 = x(2) - x0
C--------
        do g_i_ = 1, g_p_
          g_x3(g_i_) = -g_x0(g_i_) + g_x(g_i_, 3)
        enddo
        x3 = x(3) - x0
C--------
        do g_i_ = 1, g_p_
          g_y1(g_i_) = -g_y0(g_i_) + g_y(g_i_, 1)
        enddo
        y1 = y(1) - y0
C--------
        do g_i_ = 1, g_p_
          g_y2(g_i_) = -g_y0(g_i_) + g_y(g_i_, 2)
        enddo
        y2 = y(2) - y0
C--------
        do g_i_ = 1, g_p_
          g_y3(g_i_) = -g_y0(g_i_) + g_y(g_i_, 3)
        enddo
        y3 = y(3) - y0
C--------
        do g_i_ = 1, g_p_
          g_x21(g_i_) = -g_x1(g_i_) + g_x2(g_i_)
        enddo
        x21 = x2 - x1
C--------
        do g_i_ = 1, g_p_
          g_x32(g_i_) = -g_x2(g_i_) + g_x3(g_i_)
        enddo
        x32 = x3 - x2
C--------
        do g_i_ = 1, g_p_
          g_x13(g_i_) = -g_x3(g_i_) + g_x1(g_i_)
        enddo
        x13 = x1 - x3
C--------
        do g_i_ = 1, g_p_
          g_y21(g_i_) = -g_y1(g_i_) + g_y2(g_i_)
        enddo
        y21 = y2 - y1
C--------
        do g_i_ = 1, g_p_
          g_y32(g_i_) = -g_y2(g_i_) + g_y3(g_i_)
        enddo
        y32 = y3 - y2
C--------
        do g_i_ = 1, g_p_
          g_y13(g_i_) = -g_y3(g_i_) + g_y1(g_i_)
        enddo
        y13 = y1 - y3
C--------
        do g_i_ = 1, g_p_
          g_area2(g_i_) = -x21 * g_y13(g_i_) + (-y13) * g_x21(g_i_) + y2
     *1 * g_x13(g_i_) + x13 * g_y21(g_i_)
        enddo
        area2 = y21 * x13 - x21 * y13
C--------
        if (area2 .le. 0.0) then
          status = 'nega_area'
          if (area2 .eq. 0.0) then
            status = 'zero_area'
          endif
          return
        endif
C          side lenghts
        d2_v = x21 * x21
        d2_p = 2.0d0 * x21
        d4_v = y21 * y21
        d1_p = 2.0d0 * y21
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d1_p * g_y21(g_i_) + d2_p * g_x21(g_i_)
        enddo
        d1_w = d2_v + d4_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +270)
        endif
        do g_i_ = 1, g_p_
          g_l21(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        l21 = d2_v
C--------
        d2_v = x32 * x32
        d2_p = 2.0d0 * x32
        d4_v = y32 * y32
        d1_p = 2.0d0 * y32
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d1_p * g_y32(g_i_) + d2_p * g_x32(g_i_)
        enddo
        d1_w = d2_v + d4_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +291)
        endif
        do g_i_ = 1, g_p_
          g_l32(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        l32 = d2_v
C--------
        d2_v = x13 * x13
        d2_p = 2.0d0 * x13
        d4_v = y13 * y13
        d1_p = 2.0d0 * y13
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d1_p * g_y13(g_i_) + d2_p * g_x13(g_i_)
        enddo
        d1_w = d2_v + d4_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +312)
        endif
        do g_i_ = 1, g_p_
          g_l13(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        l13 = d2_v
C--------
C          side proyections
        d3_v = x3 - x1
        d5_v = x2 - x1
        d5_b = -d3_v + (-d5_v)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_v * g_x2(g_i_) + d5_b * g_x1(g_i_) + d5_v * 
     *g_x3(g_i_)
        enddo
        d1_w = d3_v * d5_v
        d4_v = y3 - y1
        d6_v = y2 - y1
        d7_b = -d4_v + (-d6_v)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_v * g_y2(g_i_) + d7_b * g_y1(g_i_) + d6_v * 
     *g_y3(g_i_) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + d4_v * d6_v
        d3_v = l21 * l21
        d1_p = 2.0d0 * l21
        d4_v = d2_w / d3_v
        d2_b = 1.0d0 / d3_v
        d4_b = -d4_v / d3_v * d1_p
        do g_i_ = 1, g_p_
          g_bl2(g_i_) = d4_b * g_l21(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        bl2 = d4_v
C--------
        do g_i_ = 1, g_p_
          g_al2(g_i_) = -g_bl2(g_i_)
        enddo
        al2 = 1.0d0 - bl2
C--------
        d3_v = x1 - x2
        d5_v = x3 - x2
        d5_b = -d3_v + (-d5_v)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_v * g_x3(g_i_) + d5_b * g_x2(g_i_) + d5_v * 
     *g_x1(g_i_)
        enddo
        d1_w = d3_v * d5_v
        d4_v = y1 - y2
        d6_v = y3 - y2
        d7_b = -d4_v + (-d6_v)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_v * g_y3(g_i_) + d7_b * g_y2(g_i_) + d6_v * 
     *g_y1(g_i_) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + d4_v * d6_v
        d3_v = l32 * l32
        d1_p = 2.0d0 * l32
        d4_v = d2_w / d3_v
        d2_b = 1.0d0 / d3_v
        d4_b = -d4_v / d3_v * d1_p
        do g_i_ = 1, g_p_
          g_bl3(g_i_) = d4_b * g_l32(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        bl3 = d4_v
C--------
        do g_i_ = 1, g_p_
          g_al3(g_i_) = -g_bl3(g_i_)
        enddo
        al3 = 1.0d0 - bl3
C--------
        d3_v = x2 - x3
        d5_v = x1 - x3
        d5_b = -d3_v + (-d5_v)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_v * g_x1(g_i_) + d5_b * g_x3(g_i_) + d5_v * 
     *g_x2(g_i_)
        enddo
        d1_w = d3_v * d5_v
        d4_v = y2 - y3
        d6_v = y1 - y3
        d7_b = -d4_v + (-d6_v)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_v * g_y1(g_i_) + d7_b * g_y3(g_i_) + d6_v * 
     *g_y2(g_i_) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + d4_v * d6_v
        d3_v = l13 * l13
        d1_p = 2.0d0 * l13
        d4_v = d2_w / d3_v
        d2_b = 1.0d0 / d3_v
        d4_b = -d4_v / d3_v * d1_p
        do g_i_ = 1, g_p_
          g_bl1(g_i_) = d4_b * g_l13(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        bl1 = d4_v
C--------
        do g_i_ = 1, g_p_
          g_al1(g_i_) = -g_bl1(g_i_)
        enddo
        al1 = 1.0d0 - bl1
C--------
C          inverse of the matrix relating inside curvatures
C          xx,yy,xy with boundary curvatures
C
        d2_v = area2 ** ( 3 - 2)
        d2_v =  d2_v * area2
        d1_p =  3 *  d2_v
        d2_v =  d2_v * area2
        do g_i_ = 1, g_p_
          g_cc(g_i_) = d1_p * g_area2(g_i_)
        enddo
        cc = d2_v
C--------
        d4_v = -x21 * y21
        d6_v = y32 * y32
        d1_p = 2.0d0 * y32
        d4_b = d4_v * d1_p
        d6_b = d6_v * (-x21)
        d7_b = -(d6_v * y21)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y32(g_i_) + d6_b * g_y21(g_i_) + d7_b 
     ** g_x21(g_i_)
        enddo
        d1_w = d4_v * d6_v
        d4_v = y21 * y21
        d1_p = 2.0d0 * y21
        d5_v = x32 * d4_v
        d10_v = (d1_w + d5_v * y32) / cc
        d2_b = 1.0d0 / cc
        d3_b = -d10_v / cc
        d6_b = d2_b * y32
        d7_b = d2_b * d5_v
        d8_b = d6_b * d4_v
        d10_b = d6_b * x32 * d1_p
        do g_i_ = 1, g_p_
          g_sq(g_i_, 1, 1) = d3_b * g_cc(g_i_) + d7_b * g_y32(g_i_) + d1
     *0_b * g_y21(g_i_) + d8_b * g_x32(g_i_) + d2_b * g_d1_w(g_i_)
        enddo
        sq(1, 1) = d10_v
C--------
        d3_v = x13 * y13
        d5_v = y32 * y32
        d1_p = 2.0d0 * y32
        d4_b = d3_v * d1_p
        d5_b = d5_v * y13
        d6_b = d5_v * x13
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y32(g_i_) + d6_b * g_y13(g_i_) + d5_b 
     ** g_x13(g_i_)
        enddo
        d1_w = d3_v * d5_v
        d4_v = y13 * y13
        d1_p = 2.0d0 * y13
        d5_v = x32 * d4_v
        d10_v = (d1_w - d5_v * y32) / cc
        d2_b = 1.0d0 / cc
        d3_b = -d10_v / cc
        d6_b = -d2_b * y32
        d7_b = -d2_b * d5_v
        d8_b = d6_b * d4_v
        d10_b = d6_b * x32 * d1_p
        do g_i_ = 1, g_p_
          g_sq(g_i_, 1, 2) = d3_b * g_cc(g_i_) + d7_b * g_y32(g_i_) + d1
     *0_b * g_y13(g_i_) + d8_b * g_x32(g_i_) + d2_b * g_d1_w(g_i_)
        enddo
        sq(1, 2) = d10_v
C--------
        d3_v = x21 * y21
        d5_v = y13 * y13
        d1_p = 2.0d0 * y13
        d4_b = d3_v * d1_p
        d5_b = d5_v * y21
        d6_b = d5_v * x21
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y13(g_i_) + d6_b * g_y21(g_i_) + d5_b 
     ** g_x21(g_i_)
        enddo
        d1_w = d3_v * d5_v
        d4_v = y21 * y21
        d1_p = 2.0d0 * y21
        d5_v = x13 * d4_v
        d10_v = (d1_w - d5_v * y13) / cc
        d2_b = 1.0d0 / cc
        d3_b = -d10_v / cc
        d6_b = -d2_b * y13
        d7_b = -d2_b * d5_v
        d8_b = d6_b * d4_v
        d10_b = d6_b * x13 * d1_p
        do g_i_ = 1, g_p_
          g_sq(g_i_, 1, 3) = d3_b * g_cc(g_i_) + d7_b * g_y13(g_i_) + d1
     *0_b * g_y21(g_i_) + d8_b * g_x13(g_i_) + d2_b * g_d1_w(g_i_)
        enddo
        sq(1, 3) = d10_v
C--------
        d3_v = x32 * x32
        d1_p = 2.0d0 * x32
        d4_v = x21 * d3_v
        d4_b = y21 * d3_v
        d6_b = y21 * x21 * d1_p
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_v * g_y21(g_i_) + d6_b * g_x32(g_i_) + d4_b 
     ** g_x21(g_i_)
        enddo
        d1_w = d4_v * y21
        d3_v = x21 * x21
        d1_p = 2.0d0 * x21
        d5_v = d3_v * x32
        d10_v = (d1_w - d5_v * y32) / cc
        d2_b = 1.0d0 / cc
        d3_b = -d10_v / cc
        d6_b = -d2_b * y32
        d7_b = -d2_b * d5_v
        d9_b = d6_b * d3_v
        d10_b = d6_b * x32 * d1_p
        do g_i_ = 1, g_p_
          g_sq(g_i_, 2, 1) = d3_b * g_cc(g_i_) + d7_b * g_y32(g_i_) + d9
     *_b * g_x32(g_i_) + d10_b * g_x21(g_i_) + d2_b * g_d1_w(g_i_)
        enddo
        sq(2, 1) = d10_v
C--------
        d4_v = x32 * x32
        d1_p = 2.0d0 * x32
        d5_v = -x13 * d4_v
        d6_b = y13 * (-x13) * d1_p
        d7_b = -(y13 * d4_v)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d5_v * g_y13(g_i_) + d6_b * g_x32(g_i_) + d7_b 
     ** g_x13(g_i_)
        enddo
        d1_w = d5_v * y13
        d3_v = x13 * x13
        d1_p = 2.0d0 * x13
        d5_v = d3_v * x32
        d10_v = (d1_w + d5_v * y32) / cc
        d2_b = 1.0d0 / cc
        d3_b = -d10_v / cc
        d6_b = d2_b * y32
        d7_b = d2_b * d5_v
        d9_b = d6_b * d3_v
        d10_b = d6_b * x32 * d1_p
        do g_i_ = 1, g_p_
          g_sq(g_i_, 2, 2) = d3_b * g_cc(g_i_) + d7_b * g_y32(g_i_) + d9
     *_b * g_x32(g_i_) + d10_b * g_x13(g_i_) + d2_b * g_d1_w(g_i_)
        enddo
        sq(2, 2) = d10_v
C--------
        d4_v = x13 * x13
        d1_p = 2.0d0 * x13
        d5_v = -x21 * d4_v
        d6_b = y21 * (-x21) * d1_p
        d7_b = -(y21 * d4_v)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d5_v * g_y21(g_i_) + d6_b * g_x13(g_i_) + d7_b 
     ** g_x21(g_i_)
        enddo
        d1_w = d5_v * y21
        d3_v = x21 * x21
        d1_p = 2.0d0 * x21
        d5_v = d3_v * x13
        d10_v = (d1_w + d5_v * y13) / cc
        d2_b = 1.0d0 / cc
        d3_b = -d10_v / cc
        d6_b = d2_b * y13
        d7_b = d2_b * d5_v
        d9_b = d6_b * d3_v
        d10_b = d6_b * x13 * d1_p
        do g_i_ = 1, g_p_
          g_sq(g_i_, 2, 3) = d3_b * g_cc(g_i_) + d7_b * g_y13(g_i_) + d9
     *_b * g_x13(g_i_) + d10_b * g_x21(g_i_) + d2_b * g_d1_w(g_i_)
        enddo
        sq(2, 3) = d10_v
C--------
        d2_v = x21 * x21
        d4_p = 2.0d0 * x21
        d4_v = y32 * y32
        d3_p = 2.0d0 * y32
        d7_v = x32 * x32
        d2_p = 2.0d0 * x32
        d9_v = y21 * y21
        d1_p = 2.0d0 * y21
        d13_v = (d2_v * d4_v - d7_v * d9_v) / cc
        d2_b = 1.0d0 / cc
        d3_b = -d13_v / cc
        d8_b = -d2_b * d7_v * d1_p
        d9_b = -d2_b * d9_v * d2_p
        d12_b = d2_b * d2_v * d3_p
        d13_b = d2_b * d4_v * d4_p
        do g_i_ = 1, g_p_
          g_sq(g_i_, 3, 1) = d3_b * g_cc(g_i_) + d8_b * g_y21(g_i_) + d9
     *_b * g_x32(g_i_) + d12_b * g_y32(g_i_) + d13_b * g_x21(g_i_)
        enddo
        sq(3, 1) = d13_v
C--------
        d2_v = x13 * x13
        d4_p = 2.0d0 * x13
        d5_v = y32 * y32
        d3_p = 2.0d0 * y32
        d8_v = x32 * x32
        d2_p = 2.0d0 * x32
        d10_v = y13 * y13
        d1_p = 2.0d0 * y13
        d14_v = (-d2_v * d5_v + d8_v * d10_v) / cc
        d2_b = 1.0d0 / cc
        d3_b = -d14_v / cc
        d8_b = d2_b * d8_v * d1_p
        d9_b = d2_b * d10_v * d2_p
        d12_b = d2_b * (-d2_v) * d3_p
        d14_b = -(d2_b * d5_v) * d4_p
        do g_i_ = 1, g_p_
          g_sq(g_i_, 3, 2) = d3_b * g_cc(g_i_) + d8_b * g_y13(g_i_) + d9
     *_b * g_x32(g_i_) + d12_b * g_y32(g_i_) + d14_b * g_x13(g_i_)
        enddo
        sq(3, 2) = d14_v
C--------
        d2_v = x21 * x21
        d4_p = 2.0d0 * x21
        d5_v = y13 * y13
        d3_p = 2.0d0 * y13
        d8_v = x13 * x13
        d2_p = 2.0d0 * x13
        d10_v = y21 * y21
        d1_p = 2.0d0 * y21
        d14_v = (-d2_v * d5_v + d8_v * d10_v) / cc
        d2_b = 1.0d0 / cc
        d3_b = -d14_v / cc
        d8_b = d2_b * d8_v * d1_p
        d9_b = d2_b * d10_v * d2_p
        d12_b = d2_b * (-d2_v) * d3_p
        d14_b = -(d2_b * d5_v) * d4_p
        do g_i_ = 1, g_p_
          g_sq(g_i_, 3, 3) = d3_b * g_cc(g_i_) + d8_b * g_y21(g_i_) + d9
     *_b * g_x13(g_i_) + d12_b * g_y13(g_i_) + d14_b * g_x21(g_i_)
        enddo
        sq(3, 3) = d14_v
C--------
C     print '(''area'',f8.3)',area
C     print '(''inver'',3f8.3)',((sq(i,j),j=1,3),i=1,3)
        do g_i_ = 1, g_p_
          g_d11(g_i_) = g_db(g_i_, 1, 1)
        enddo
        d11 = db(1, 1)
C--------
        do g_i_ = 1, g_p_
          g_d22(g_i_) = g_db(g_i_, 2, 2)
        enddo
        d22 = db(2, 2)
C--------
        do g_i_ = 1, g_p_
          g_d33(g_i_) = g_db(g_i_, 3, 3)
        enddo
        d33 = db(3, 3)
C--------
        do g_i_ = 1, g_p_
          g_d12(g_i_) = g_db(g_i_, 1, 2)
        enddo
        d12 = db(1, 2)
C--------
        do g_i_ = 1, g_p_
          g_d13(g_i_) = g_db(g_i_, 1, 3)
        enddo
        d13 = db(1, 3)
C--------
        do g_i_ = 1, g_p_
          g_d23(g_i_) = g_db(g_i_, 2, 3)
        enddo
        d23 = db(2, 3)
C--------
        d2_b = dble(0.5)
        do g_i_ = 1, g_p_
          g_area(g_i_) = d2_b * g_area2(g_i_)
        enddo
        area = dble(0.5) * area2
C--------
        do 99990 j = 1, 3
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d12 * g_sq(g_i_, 2, j) + sq(2, j) * g_d12(g_i
     *_) + d11 * g_sq(g_i_, 1, j) + sq(1, j) * g_d11(g_i_)
          enddo
          d1_w = d11 * sq(1, j) + d12 * sq(2, j)
          do g_i_ = 1, g_p_
            g_s1(g_i_) = d13 * g_sq(g_i_, 3, j) + sq(3, j) * g_d13(g_i_)
     * + g_d1_w(g_i_)
          enddo
          s1 = d1_w + d13 * sq(3, j)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d22 * g_sq(g_i_, 2, j) + sq(2, j) * g_d22(g_i
     *_) + d12 * g_sq(g_i_, 1, j) + sq(1, j) * g_d12(g_i_)
          enddo
          d1_w = d12 * sq(1, j) + d22 * sq(2, j)
          do g_i_ = 1, g_p_
            g_s2(g_i_) = d23 * g_sq(g_i_, 3, j) + sq(3, j) * g_d23(g_i_)
     * + g_d1_w(g_i_)
          enddo
          s2 = d1_w + d23 * sq(3, j)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d23 * g_sq(g_i_, 2, j) + sq(2, j) * g_d23(g_i
     *_) + d13 * g_sq(g_i_, 1, j) + sq(1, j) * g_d13(g_i_)
          enddo
          d1_w = d13 * sq(1, j) + d23 * sq(2, j)
          do g_i_ = 1, g_p_
            g_s3(g_i_) = d33 * g_sq(g_i_, 3, j) + sq(3, j) * g_d33(g_i_)
     * + g_d1_w(g_i_)
          enddo
          s3 = d1_w + d33 * sq(3, j)
C--------
          do 99991 i = 1, j
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = s2 * g_sq(g_i_, 2, i) + sq(2, i) * g_s2(g_i
     *_) + s1 * g_sq(g_i_, 1, i) + sq(1, i) * g_s1(g_i_)
            enddo
            d1_w = s1 * sq(1, i) + s2 * sq(2, i)
            do g_i_ = 1, g_p_
              g_sds(g_i_, i, j) = s3 * g_sq(g_i_, 3, i) + sq(3, i) * g_s
     *3(g_i_) + g_d1_w(g_i_) + g_sds(g_i_, i, j)
            enddo
            sds(i, j) = sds(i, j) + (d1_w + s3 * sq(3, i))
C--------
            do g_i_ = 1, g_p_
              g_sds(g_i_, j, i) = g_sds(g_i_, i, j)
            enddo
            sds(j, i) = sds(i, j)
C--------
1500        continue
99991     continue
2000      continue
99990   continue
        do 99988 j = 1, 3
          do 99989 i = 1, 3
2100        d3_b = f * (1.0d0 / dble(3.0))
            d4_b = d3_b * area
            d5_b = d3_b * sds(i, j)
            do g_i_ = 1, g_p_
              g_sds(g_i_, i, j) = d5_b * g_area(g_i_) + d4_b * g_sds(g_i
     *_, i, j)
            enddo
            sds(i, j) = sds(i, j) * area / dble(3.0) * f
C--------
99989     continue
99988   continue
C     print '(''sds'',3f8.3)',((sds(i,j),j=1,3),i=1,3)
C
C    
C         matrix q
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 1) = 0.0d0
        enddo
        q(1, 1) = 6.0d0
C--------
        d2_b = dble(-2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 2) = d2_b * g_y13(g_i_)
        enddo
        q(1, 2) = dble(-2.0) * y13
C--------
        d2_b = dble(2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 3) = d2_b * g_x13(g_i_)
        enddo
        q(1, 3) = dble(2.0) * x13
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 7) = 0.0d0
        enddo
        q(1, 7) = -6.0d0
C--------
        d2_b = dble(-4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 8) = d2_b * g_y13(g_i_)
        enddo
        q(1, 8) = dble(-4.0) * y13
C--------
        d2_b = dble(4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 9) = d2_b * g_x13(g_i_)
        enddo
        q(1, 9) = dble(4.0) * x13
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 1) = 0.0d0
        enddo
        q(2, 1) = -6.0d0
C--------
        d2_b = dble(4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 2) = d2_b * g_y13(g_i_)
        enddo
        q(2, 2) = dble(4.0) * y13
C--------
        d2_b = dble(-4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 3) = d2_b * g_x13(g_i_)
        enddo
        q(2, 3) = dble(-4.0) * x13
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 7) = 0.0d0
        enddo
        q(2, 7) = 6.0d0
C--------
        d2_b = dble(2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 8) = d2_b * g_y13(g_i_)
        enddo
        q(2, 8) = dble(2.0) * y13
C--------
        d2_b = dble(-2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 9) = d2_b * g_x13(g_i_)
        enddo
        q(2, 9) = dble(-2.0) * x13
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 1) = 0.0d0
        enddo
        q(3, 1) = -6.0d0
C--------
        d2_b = dble(-4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 2) = d2_b * g_y21(g_i_)
        enddo
        q(3, 2) = dble(-4.0) * y21
C--------
        d2_b = dble(4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 3) = d2_b * g_x21(g_i_)
        enddo
        q(3, 3) = dble(4.0) * x21
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 4) = 0.0d0
        enddo
        q(3, 4) = 6.0d0
C--------
        d2_b = dble(-2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 5) = d2_b * g_y21(g_i_)
        enddo
        q(3, 5) = dble(-2.0) * y21
C--------
        d2_b = dble(2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 6) = d2_b * g_x21(g_i_)
        enddo
        q(3, 6) = dble(2.0) * x21
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 1) = 0.0d0
        enddo
        q(4, 1) = 6.0d0
C--------
        d2_b = dble(2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 2) = d2_b * g_y21(g_i_)
        enddo
        q(4, 2) = dble(2.0) * y21
C--------
        d2_b = dble(-2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 3) = d2_b * g_x21(g_i_)
        enddo
        q(4, 3) = dble(-2.0) * x21
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 4) = 0.0d0
        enddo
        q(4, 4) = -6.0d0
C--------
        d2_b = dble(4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 5) = d2_b * g_y21(g_i_)
        enddo
        q(4, 5) = dble(4.0) * y21
C--------
        d2_b = dble(-4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 6) = d2_b * g_x21(g_i_)
        enddo
        q(4, 6) = dble(-4.0) * x21
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 4) = 0.0d0
        enddo
        q(5, 4) = -6.0d0
C--------
        d2_b = dble(-4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 5) = d2_b * g_y32(g_i_)
        enddo
        q(5, 5) = dble(-4.0) * y32
C--------
        d2_b = dble(4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 6) = d2_b * g_x32(g_i_)
        enddo
        q(5, 6) = dble(4.0) * x32
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 7) = 0.0d0
        enddo
        q(5, 7) = 6.0d0
C--------
        d2_b = dble(-2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 8) = d2_b * g_y32(g_i_)
        enddo
        q(5, 8) = dble(-2.0) * y32
C--------
        d2_b = dble(2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 9) = d2_b * g_x32(g_i_)
        enddo
        q(5, 9) = dble(2.0) * x32
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 4) = 0.0d0
        enddo
        q(6, 4) = 6.0d0
C--------
        d2_b = dble(2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 5) = d2_b * g_y32(g_i_)
        enddo
        q(6, 5) = dble(2.0) * y32
C--------
        d2_b = dble(-2.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 6) = d2_b * g_x32(g_i_)
        enddo
        q(6, 6) = dble(-2.0) * x32
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 7) = 0.0d0
        enddo
        q(6, 7) = -6.0d0
C--------
        d2_b = dble(4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 8) = d2_b * g_y32(g_i_)
        enddo
        q(6, 8) = dble(4.0) * y32
C--------
        d2_b = dble(-4.0)
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 9) = d2_b * g_x32(g_i_)
        enddo
        q(6, 9) = dble(-4.0) * x32
C--------
C     print '(''q'',9f5.1)',((q(i,j),j=1,9),i=1,6)
C
C    numerical integration
C
        do 99985 k = 1, 3
          l1 = pg(k, 1)
          l2 = pg(k, 2)
          l3 = pg(k, 3)
C
C    compute rm in the integration point
C
          d5_b = -(1.0d0 / dble(3.)) + l2
          do g_i_ = 1, g_p_
            g_rm(g_i_, 1, 1) = d5_b * g_al1(g_i_)
          enddo
          rm(1, 1) = l3 + al1 * l2 - (1.0d0 + al1) / dble(3.)
C--------
          d5_b = -(1.0d0 / dble(3.)) + l2
          do g_i_ = 1, g_p_
            g_rm(g_i_, 1, 2) = d5_b * g_bl1(g_i_)
          enddo
          rm(1, 2) = l1 + bl1 * l2 - (1.0d0 + bl1) / dble(3.)
C--------
          d5_b = -(1.0d0 / dble(3.)) + l3
          do g_i_ = 1, g_p_
            g_rm(g_i_, 2, 3) = d5_b * g_al2(g_i_)
          enddo
          rm(2, 3) = l1 + al2 * l3 - (1.0d0 + al2) / dble(3.)
C--------
          d5_b = -(1.0d0 / dble(3.)) + l3
          do g_i_ = 1, g_p_
            g_rm(g_i_, 2, 4) = d5_b * g_bl2(g_i_)
          enddo
          rm(2, 4) = l2 + bl2 * l3 - (1.0d0 + bl2) / dble(3.)
C--------
          d5_b = -(1.0d0 / dble(3.)) + l1
          do g_i_ = 1, g_p_
            g_rm(g_i_, 3, 5) = d5_b * g_al3(g_i_)
          enddo
          rm(3, 5) = l2 + al3 * l1 - (1.0d0 + al3) / dble(3.)
C--------
          d5_b = -(1.0d0 / dble(3.)) + l1
          do g_i_ = 1, g_p_
            g_rm(g_i_, 3, 6) = d5_b * g_bl3(g_i_)
          enddo
          rm(3, 6) = l3 + bl3 * l1 - (1.0d0 + bl3) / dble(3.)
C--------
          do 99986 j = 1, 6
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = sds(1, 2) * g_rm(g_i_, 2, j) + rm(2, j) * g
     *_sds(g_i_, 1, 2) + sds(1, 1) * g_rm(g_i_, 1, j) + rm(1, j) * g_sds
     *(g_i_, 1, 1)
            enddo
            d1_w = sds(1, 1) * rm(1, j) + sds(1, 2) * rm(2, j)
            do g_i_ = 1, g_p_
              g_s1(g_i_) = sds(1, 3) * g_rm(g_i_, 3, j) + rm(3, j) * g_s
     *ds(g_i_, 1, 3) + g_d1_w(g_i_)
            enddo
            s1 = d1_w + sds(1, 3) * rm(3, j)
C--------
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = sds(2, 2) * g_rm(g_i_, 2, j) + rm(2, j) * g
     *_sds(g_i_, 2, 2) + sds(2, 1) * g_rm(g_i_, 1, j) + rm(1, j) * g_sds
     *(g_i_, 2, 1)
            enddo
            d1_w = sds(2, 1) * rm(1, j) + sds(2, 2) * rm(2, j)
            do g_i_ = 1, g_p_
              g_s2(g_i_) = sds(2, 3) * g_rm(g_i_, 3, j) + rm(3, j) * g_s
     *ds(g_i_, 2, 3) + g_d1_w(g_i_)
            enddo
            s2 = d1_w + sds(2, 3) * rm(3, j)
C--------
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = sds(3, 2) * g_rm(g_i_, 2, j) + rm(2, j) * g
     *_sds(g_i_, 3, 2) + sds(3, 1) * g_rm(g_i_, 1, j) + rm(1, j) * g_sds
     *(g_i_, 3, 1)
            enddo
            d1_w = sds(3, 1) * rm(1, j) + sds(3, 2) * rm(2, j)
            do g_i_ = 1, g_p_
              g_s3(g_i_) = sds(3, 3) * g_rm(g_i_, 3, j) + rm(3, j) * g_s
     *ds(g_i_, 3, 3) + g_d1_w(g_i_)
            enddo
            s3 = d1_w + sds(3, 3) * rm(3, j)
C--------
            do 99987 i = 1, j
              do g_i_ = 1, g_p_
                g_d1_w(g_i_) = s2 * g_rm(g_i_, 2, i) + rm(2, i) * g_s2(g
     *_i_) + s1 * g_rm(g_i_, 1, i) + rm(1, i) * g_s1(g_i_)
              enddo
              d1_w = s1 * rm(1, i) + s2 * rm(2, i)
              do g_i_ = 1, g_p_
                g_rsd(g_i_, i, j) = s3 * g_rm(g_i_, 3, i) + rm(3, i) * g
     *_s3(g_i_) + g_d1_w(g_i_) + g_rsd(g_i_, i, j)
              enddo
              rsd(i, j) = rsd(i, j) + (d1_w + s3 * rm(3, i))
C--------
              do g_i_ = 1, g_p_
                g_rsd(g_i_, j, i) = g_rsd(g_i_, i, j)
              enddo
              rsd(j, i) = rsd(i, j)
C--------
2300          continue
99987       continue
2400        continue
99986     continue
2500      continue
99985   continue
C     print '(''int'',6f7.3)',((rsd(i,j),j=1,6),i=1,6)
C
C    computing the kh matrix
C
        do 99983 j = 1, 9
          k = ls(j)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = rsd(1, 2) * g_q(g_i_, 2, j) + q(2, j) * g_rsd
     *(g_i_, 1, 2) + rsd(1, 1) * g_q(g_i_, 1, j) + q(1, j) * g_rsd(g_i_,
     * 1, 1)
          enddo
          d1_w = rsd(1, 1) * q(1, j) + rsd(1, 2) * q(2, j)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = rsd(1, 4) * g_q(g_i_, 4, j) + q(4, j) * g_rsd
     *(g_i_, 1, 4) + rsd(1, 3) * g_q(g_i_, 3, j) + q(3, j) * g_rsd(g_i_,
     * 1, 3) + g_d1_w(g_i_)
          enddo
          d2_w = d1_w + rsd(1, 3) * q(3, j) + rsd(1, 4) * q(4, j)
          do g_i_ = 1, g_p_
            g_s1(g_i_) = rsd(1, 6) * g_q(g_i_, 6, j) + q(6, j) * g_rsd(g
     *_i_, 1, 6) + rsd(1, 5) * g_q(g_i_, 5, j) + q(5, j) * g_rsd(g_i_, 1
     *, 5) + g_d2_w(g_i_)
          enddo
          s1 = d2_w + rsd(1, 5) * q(5, j) + rsd(1, 6) * q(6, j)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = rsd(2, 2) * g_q(g_i_, 2, j) + q(2, j) * g_rsd
     *(g_i_, 2, 2) + rsd(2, 1) * g_q(g_i_, 1, j) + q(1, j) * g_rsd(g_i_,
     * 2, 1)
          enddo
          d1_w = rsd(2, 1) * q(1, j) + rsd(2, 2) * q(2, j)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = rsd(2, 4) * g_q(g_i_, 4, j) + q(4, j) * g_rsd
     *(g_i_, 2, 4) + rsd(2, 3) * g_q(g_i_, 3, j) + q(3, j) * g_rsd(g_i_,
     * 2, 3) + g_d1_w(g_i_)
          enddo
          d2_w = d1_w + rsd(2, 3) * q(3, j) + rsd(2, 4) * q(4, j)
          do g_i_ = 1, g_p_
            g_s2(g_i_) = rsd(2, 6) * g_q(g_i_, 6, j) + q(6, j) * g_rsd(g
     *_i_, 2, 6) + rsd(2, 5) * g_q(g_i_, 5, j) + q(5, j) * g_rsd(g_i_, 2
     *, 5) + g_d2_w(g_i_)
          enddo
          s2 = d2_w + rsd(2, 5) * q(5, j) + rsd(2, 6) * q(6, j)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = rsd(3, 2) * g_q(g_i_, 2, j) + q(2, j) * g_rsd
     *(g_i_, 3, 2) + rsd(3, 1) * g_q(g_i_, 1, j) + q(1, j) * g_rsd(g_i_,
     * 3, 1)
          enddo
          d1_w = rsd(3, 1) * q(1, j) + rsd(3, 2) * q(2, j)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = rsd(3, 4) * g_q(g_i_, 4, j) + q(4, j) * g_rsd
     *(g_i_, 3, 4) + rsd(3, 3) * g_q(g_i_, 3, j) + q(3, j) * g_rsd(g_i_,
     * 3, 3) + g_d1_w(g_i_)
          enddo
          d2_w = d1_w + rsd(3, 3) * q(3, j) + rsd(3, 4) * q(4, j)
          do g_i_ = 1, g_p_
            g_s3(g_i_) = rsd(3, 6) * g_q(g_i_, 6, j) + q(6, j) * g_rsd(g
     *_i_, 3, 6) + rsd(3, 5) * g_q(g_i_, 5, j) + q(5, j) * g_rsd(g_i_, 3
     *, 5) + g_d2_w(g_i_)
          enddo
          s3 = d2_w + rsd(3, 5) * q(5, j) + rsd(3, 6) * q(6, j)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = rsd(4, 2) * g_q(g_i_, 2, j) + q(2, j) * g_rsd
     *(g_i_, 4, 2) + rsd(4, 1) * g_q(g_i_, 1, j) + q(1, j) * g_rsd(g_i_,
     * 4, 1)
          enddo
          d1_w = rsd(4, 1) * q(1, j) + rsd(4, 2) * q(2, j)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = rsd(4, 4) * g_q(g_i_, 4, j) + q(4, j) * g_rsd
     *(g_i_, 4, 4) + rsd(4, 3) * g_q(g_i_, 3, j) + q(3, j) * g_rsd(g_i_,
     * 4, 3) + g_d1_w(g_i_)
          enddo
          d2_w = d1_w + rsd(4, 3) * q(3, j) + rsd(4, 4) * q(4, j)
          do g_i_ = 1, g_p_
            g_s4(g_i_) = rsd(4, 6) * g_q(g_i_, 6, j) + q(6, j) * g_rsd(g
     *_i_, 4, 6) + rsd(4, 5) * g_q(g_i_, 5, j) + q(5, j) * g_rsd(g_i_, 4
     *, 5) + g_d2_w(g_i_)
          enddo
          s4 = d2_w + rsd(4, 5) * q(5, j) + rsd(4, 6) * q(6, j)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = rsd(5, 2) * g_q(g_i_, 2, j) + q(2, j) * g_rsd
     *(g_i_, 5, 2) + rsd(5, 1) * g_q(g_i_, 1, j) + q(1, j) * g_rsd(g_i_,
     * 5, 1)
          enddo
          d1_w = rsd(5, 1) * q(1, j) + rsd(5, 2) * q(2, j)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = rsd(5, 4) * g_q(g_i_, 4, j) + q(4, j) * g_rsd
     *(g_i_, 5, 4) + rsd(5, 3) * g_q(g_i_, 3, j) + q(3, j) * g_rsd(g_i_,
     * 5, 3) + g_d1_w(g_i_)
          enddo
          d2_w = d1_w + rsd(5, 3) * q(3, j) + rsd(5, 4) * q(4, j)
          do g_i_ = 1, g_p_
            g_s5(g_i_) = rsd(5, 6) * g_q(g_i_, 6, j) + q(6, j) * g_rsd(g
     *_i_, 5, 6) + rsd(5, 5) * g_q(g_i_, 5, j) + q(5, j) * g_rsd(g_i_, 5
     *, 5) + g_d2_w(g_i_)
          enddo
          s5 = d2_w + rsd(5, 5) * q(5, j) + rsd(5, 6) * q(6, j)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = rsd(6, 2) * g_q(g_i_, 2, j) + q(2, j) * g_rsd
     *(g_i_, 6, 2) + rsd(6, 1) * g_q(g_i_, 1, j) + q(1, j) * g_rsd(g_i_,
     * 6, 1)
          enddo
          d1_w = rsd(6, 1) * q(1, j) + rsd(6, 2) * q(2, j)
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = rsd(6, 4) * g_q(g_i_, 4, j) + q(4, j) * g_rsd
     *(g_i_, 6, 4) + rsd(6, 3) * g_q(g_i_, 3, j) + q(3, j) * g_rsd(g_i_,
     * 6, 3) + g_d1_w(g_i_)
          enddo
          d2_w = d1_w + rsd(6, 3) * q(3, j) + rsd(6, 4) * q(4, j)
          do g_i_ = 1, g_p_
            g_s6(g_i_) = rsd(6, 6) * g_q(g_i_, 6, j) + q(6, j) * g_rsd(g
     *_i_, 6, 6) + rsd(6, 5) * g_q(g_i_, 5, j) + q(5, j) * g_rsd(g_i_, 6
     *, 5) + g_d2_w(g_i_)
          enddo
          s6 = d2_w + rsd(6, 5) * q(5, j) + rsd(6, 6) * q(6, j)
C--------
          do 99984 i = 1, j
            l = ls(i)
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = s2 * g_q(g_i_, 2, i) + q(2, i) * g_s2(g_i_)
     * + s1 * g_q(g_i_, 1, i) + q(1, i) * g_s1(g_i_)
            enddo
            d1_w = s1 * q(1, i) + s2 * q(2, i)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = s4 * g_q(g_i_, 4, i) + q(4, i) * g_s4(g_i_)
     * + s3 * g_q(g_i_, 3, i) + q(3, i) * g_s3(g_i_) + g_d1_w(g_i_)
            enddo
            d2_w = d1_w + s3 * q(3, i) + s4 * q(4, i)
            do g_i_ = 1, g_p_
              g_d3_w(g_i_) = s6 * g_q(g_i_, 6, i) + q(6, i) * g_s6(g_i_)
     * + s5 * g_q(g_i_, 5, i) + q(5, i) * g_s5(g_i_) + g_d2_w(g_i_)
            enddo
            d3_w = d2_w + s5 * q(5, i) + s6 * q(6, i)
            do g_i_ = 1, g_p_
              g_sm(g_i_, l, k) = g_d3_w(g_i_) + g_sm(g_i_, l, k)
            enddo
            sm(l, k) = sm(l, k) + d3_w
C--------
            do g_i_ = 1, g_p_
              g_sm(g_i_, k, l) = g_sm(g_i_, l, k)
            enddo
            sm(k, l) = sm(l, k)
C--------
2800        continue
99984     continue
2600      continue
99983   continue
C      write(46,5000) ((sm(i,j),i=1,9),j=1,9)
C5000  format(9f10.3)
        return
      end
C=end fortran
