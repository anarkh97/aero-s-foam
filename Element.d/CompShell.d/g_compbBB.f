C                           DISCLAIMER
C
C   This file was generated on 12/17/98 by the version of
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
C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine g_compbbb(g_p_, elm, type, x, g_x, ldg_x, y, g_y, ldg_y
     *, db, g_db, ldg_db, f, clr, cqr, rowb, colb, rot, g_rot, ldg_rot, 
     *l, g_l, ldg_l, kbbb, g_kbbb, ldg_kbbb)
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
        integer elm, type, rowb(9), colb(9)
        real*8 x(3), y(3), db(3, 3), kbbb(18, 18), rot(6, 6)
        real*8 f, clr, cqr, l(9, 3)
C
C.....Local Variables
C
        integer i, j, row, col
        real*8 zero, twicearea, llr(9, 3), lqr(9, 3)
        real*8 x21, x32, x13, y21, y32, y13
        real*8 x12, x23, x31, y12, y23, y31
        real*8 x0, y0, dist12, dist23, dist31
        real*8 c12, c23, c31, s12, s23, s31
        real*8 cc12, cc23, cc31, ss12, ss23, ss31
        real*8 one, cs12, cs23, cs31, area
        real*8 dblt1, dblt2, dblt3, factor
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_kbbb, ldg_l, ldg_x, ldg_y, ldg_db, ldg_r
     *ot
        double precision d2_b, d3_b, d8_b, d4_b, d1_p, d1_w, d7_b, d6_b,
     * d2_v, d3_v
        double precision d5_b, g_kbbb(ldg_kbbb, 18, 18), g_llr(g_pmax_, 
     *9, 3), g_lqr(g_pmax_, 9, 3), g_l(ldg_l, 9, 3), g_x21(g_pmax_), g_x
     *(ldg_x, 3), g_x12(g_pmax_), g_x32(g_pmax_), g_x23(g_pmax_)
        double precision g_x13(g_pmax_), g_x31(g_pmax_), g_y21(g_pmax_),
     * g_y(ldg_y, 3), g_y12(g_pmax_), g_y32(g_pmax_), g_y23(g_pmax_), g_
     *y13(g_pmax_), g_y31(g_pmax_), g_twicearea(g_pmax_)
        double precision g_area(g_pmax_), g_d1_w(g_pmax_), g_dist12(g_pm
     *ax_), g_dist23(g_pmax_), g_dist31(g_pmax_), g_c12(g_pmax_), g_s12(
     *g_pmax_), g_c23(g_pmax_), g_s23(g_pmax_), g_c31(g_pmax_)
        double precision g_s31(g_pmax_), g_cc12(g_pmax_), g_cc23(g_pmax_
     *), g_cc31(g_pmax_), g_ss12(g_pmax_), g_ss23(g_pmax_), g_ss31(g_pma
     *x_), g_cs12(g_pmax_), g_cs23(g_pmax_), g_cs31(g_pmax_)
        double precision g_factor(g_pmax_), g_db(ldg_db, 3, 3), g_dblt1(
     *g_pmax_), g_dblt2(g_pmax_), g_dblt3(g_pmax_), g_rot(ldg_rot, 6, 6)
        save g_ss23, g_ss31, g_cs12, g_cs23, g_cs31, g_factor, g_dblt1, 
     *g_dblt2, g_dblt3
        save g_c12, g_s12, g_c23, g_s23, g_c31, g_s31, g_cc12, g_cc23, g
     *_cc31, g_ss12
        save g_y32, g_y23, g_y13, g_y31, g_twicearea, g_area, g_d1_w, g_
     *dist12, g_dist23, g_dist31
        save g_llr, g_lqr, g_x21, g_x12, g_x32, g_x23, g_x13, g_x31, g_y
     *21, g_y12
        data zero /0.000000d+00/
        data one /1.000000d+00/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE OUTPUT STIFFNESS MATRIX
C
        integer g_ehfid
        data g_ehfid /0/
C

C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do 99998 j = 1, 18
          do 99999 i = 1, 18
            do g_i_ = 1, g_p_
              g_kbbb(g_i_, i, j) = 0.0d0
            enddo
            kbbb(i, j) = zero
C--------
1002        continue
99999     continue
1001      continue
99998   continue
C
C.....CLEAR THE LOCAL MATRICES
C
        do 99994 j = 1, 3
          do 99997 i = 1, 9
            do g_i_ = 1, g_p_
              g_llr(g_i_, i, j) = 0.0d0
            enddo
            llr(i, j) = zero
C--------
1004        continue
99997     continue
          do 99996 i = 1, 9
            do g_i_ = 1, g_p_
              g_lqr(g_i_, i, j) = 0.0d0
            enddo
            lqr(i, j) = zero
C--------
1005        continue
99996     continue
          do 99995 i = 1, 9
            do g_i_ = 1, g_p_
              g_l(g_i_, i, j) = 0.0d0
            enddo
            l(i, j) = zero
C--------
1006        continue
99995     continue
1003      continue
99994   continue
C
C.....RETURN IF THE STIFFNESS FACTOR IS ZERO
C
        if (f .eq. zero) then
          return
        endif
C
C.....CHECK IF [CLR] AND [CQR] SATISFY THE CONSTRAINT [CLR]+[CQR]=1
C
        if ((clr + cqr) .ne. one) then
          goto 100
        endif
C
C.....CHECK IF THE STIFFNESS FACTOR [F] IS POSITIVE
C
        if (f .lt. zero) then
          goto 200
        endif
C
C.....GET THE DISTANCES BETWEEN NODAL POINT X- AND Y- COORDINATES
C
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
C
C.....CALCULATE TWICE THE AREA OF THE TRIANGLE
C
        do g_i_ = 1, g_p_
          g_twicearea(g_i_) = -x21 * g_y13(g_i_) + (-y13) * g_x21(g_i_) 
     *+ y21 * g_x13(g_i_) + x13 * g_y21(g_i_)
        enddo
        twicearea = y21 * x13 - x21 * y13
C--------
C
C.....CALCULATE THE AREA OF THE TRIANGLE
C
        do g_i_ = 1, g_p_
          g_area(g_i_) = 0.50d+00 * g_twicearea(g_i_)
        enddo
        area = 0.50d+00 * twicearea
C--------
C
C.....CHECK THE AREA (ERROR IF NOT POSITIVE)
C
        if (twicearea .le. zero) then
          goto 300
        endif
C
C.....GET THE COORDINATES OF THE CENTROID OF THE TRIANGLE
C
        x0 = (x(1) + x(2) + x(3)) / 3.00d+00
        y0 = (y(1) + y(2) + y(3)) / 3.00d+00
C
C.....GET THE DISTANCES BETWEEN NODES 1-2, 2-3 AND 3-1
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
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compbBB.f',
     +262)
        endif
        do g_i_ = 1, g_p_
          g_dist12(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        dist12 = d2_v
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
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compbBB.f',
     +281)
        endif
        do g_i_ = 1, g_p_
          g_dist23(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        dist23 = d2_v
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
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compbBB.f',
     +300)
        endif
        do g_i_ = 1, g_p_
          g_dist31(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        dist31 = d2_v
C--------
C
C.....ASSEMBLE THE LOCAL MATRIX [LLR] W/ SHAPE FUNCTION DERIVATIVES
C
        if (clr .ne. zero) then
C
          do g_i_ = 1, g_p_
            g_llr(g_i_, 3, 1) = 0.50d+00 * g_y32(g_i_)
          enddo
          llr(3, 1) = 0.50d+00 * y32
C--------
          do g_i_ = 1, g_p_
            g_llr(g_i_, 6, 1) = 0.50d+00 * g_y13(g_i_)
          enddo
          llr(6, 1) = 0.50d+00 * y13
C--------
          do g_i_ = 1, g_p_
            g_llr(g_i_, 9, 1) = 0.50d+00 * g_y21(g_i_)
          enddo
          llr(9, 1) = 0.50d+00 * y21
C--------
C
          do g_i_ = 1, g_p_
            g_llr(g_i_, 2, 2) = 0.50d+00 * g_x32(g_i_)
          enddo
          llr(2, 2) = 0.50d+00 * x32
C--------
          do g_i_ = 1, g_p_
            g_llr(g_i_, 5, 2) = 0.50d+00 * g_x13(g_i_)
          enddo
          llr(5, 2) = 0.50d+00 * x13
C--------
          do g_i_ = 1, g_p_
            g_llr(g_i_, 8, 2) = 0.50d+00 * g_x21(g_i_)
          enddo
          llr(8, 2) = 0.50d+00 * x21
C--------
C
          do g_i_ = 1, g_p_
            g_llr(g_i_, 2, 3) = -0.50d+00 * g_y32(g_i_)
          enddo
          llr(2, 3) = -0.50d+00 * y32
C--------
          do g_i_ = 1, g_p_
            g_llr(g_i_, 3, 3) = -0.50d+00 * g_x32(g_i_)
          enddo
          llr(3, 3) = -0.50d+00 * x32
C--------
          do g_i_ = 1, g_p_
            g_llr(g_i_, 5, 3) = -0.50d+00 * g_y13(g_i_)
          enddo
          llr(5, 3) = -0.50d+00 * y13
C--------
          do g_i_ = 1, g_p_
            g_llr(g_i_, 6, 3) = -0.50d+00 * g_x13(g_i_)
          enddo
          llr(6, 3) = -0.50d+00 * x13
C--------
          do g_i_ = 1, g_p_
            g_llr(g_i_, 8, 3) = -0.50d+00 * g_y21(g_i_)
          enddo
          llr(8, 3) = -0.50d+00 * y21
C--------
          do g_i_ = 1, g_p_
            g_llr(g_i_, 9, 3) = -0.50d+00 * g_x21(g_i_)
          enddo
          llr(9, 3) = -0.50d+00 * x21
C--------
C
        endif
C
C.....ASSEMBLE THE LOCAL MATRIX [LQR] W/ SHAPE FUNCTION DERIVATIVES
C
        if (cqr .ne. zero) then
C
          d3_v = y21 / dist12
          d2_b = 1.0d0 / dist12
          d3_b = -d3_v / dist12
          do g_i_ = 1, g_p_
            g_c12(g_i_) = d3_b * g_dist12(g_i_) + d2_b * g_y21(g_i_)
          enddo
          c12 = d3_v
C--------
          d3_v = x12 / dist12
          d2_b = 1.0d0 / dist12
          d3_b = -d3_v / dist12
          do g_i_ = 1, g_p_
            g_s12(g_i_) = d3_b * g_dist12(g_i_) + d2_b * g_x12(g_i_)
          enddo
          s12 = d3_v
C--------
          d3_v = y32 / dist23
          d2_b = 1.0d0 / dist23
          d3_b = -d3_v / dist23
          do g_i_ = 1, g_p_
            g_c23(g_i_) = d3_b * g_dist23(g_i_) + d2_b * g_y32(g_i_)
          enddo
          c23 = d3_v
C--------
          d3_v = x23 / dist23
          d2_b = 1.0d0 / dist23
          d3_b = -d3_v / dist23
          do g_i_ = 1, g_p_
            g_s23(g_i_) = d3_b * g_dist23(g_i_) + d2_b * g_x23(g_i_)
          enddo
          s23 = d3_v
C--------
          d3_v = y13 / dist31
          d2_b = 1.0d0 / dist31
          d3_b = -d3_v / dist31
          do g_i_ = 1, g_p_
            g_c31(g_i_) = d3_b * g_dist31(g_i_) + d2_b * g_y13(g_i_)
          enddo
          c31 = d3_v
C--------
          d3_v = x31 / dist31
          d2_b = 1.0d0 / dist31
          d3_b = -d3_v / dist31
          do g_i_ = 1, g_p_
            g_s31(g_i_) = d3_b * g_dist31(g_i_) + d2_b * g_x31(g_i_)
          enddo
          s31 = d3_v
C--------
C
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
C
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
C
          d5_b = 0.50d+00 * x31
          d6_b = 0.50d+00 * cc31
          d7_b = 0.50d+00 * x12
          d8_b = 0.50d+00 * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 2, 1) = d6_b * g_x31(g_i_) + d5_b * g_cc31(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_cc12(g_i_)
          enddo
          lqr(2, 1) = 0.50d+00 * (cc12 * x12 + cc31 * x31)
C--------
          d5_b = 0.50d+00 * x31
          d6_b = 0.50d+00 * ss31
          d7_b = 0.50d+00 * x12
          d8_b = 0.50d+00 * ss12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 2, 2) = d6_b * g_x31(g_i_) + d5_b * g_ss31(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_ss12(g_i_)
          enddo
          lqr(2, 2) = 0.50d+00 * (ss12 * x12 + ss31 * x31)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 2, 3) = ss31 * g_y13(g_i_) + y13 * g_ss31(g_i_) 
     *+ ss12 * g_y21(g_i_) + y21 * g_ss12(g_i_)
          enddo
          lqr(2, 3) = ss12 * y21 + ss31 * y13
C--------
C
          d5_b = -0.50d+00 * y13
          d6_b = -0.50d+00 * cc31
          d7_b = -0.50d+00 * y21
          d8_b = -0.50d+00 * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 3, 1) = d6_b * g_y13(g_i_) + d5_b * g_cc31(g_i_)
     * + d8_b * g_y21(g_i_) + d7_b * g_cc12(g_i_)
          enddo
          lqr(3, 1) = -0.50d+00 * (cc12 * y21 + cc31 * y13)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 3, 2) = -0.50d+00 * g_lqr(g_i_, 2, 3)
          enddo
          lqr(3, 2) = -0.50d+00 * lqr(2, 3)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 3, 3) = -2.00d+00 * g_lqr(g_i_, 2, 1)
          enddo
          lqr(3, 3) = -2.00d+00 * lqr(2, 1)
C--------
C
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
C
          d5_b = 0.50d+00 * x23
          d6_b = 0.50d+00 * cc23
          d7_b = 0.50d+00 * x12
          d8_b = 0.50d+00 * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 5, 1) = d6_b * g_x23(g_i_) + d5_b * g_cc23(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_cc12(g_i_)
          enddo
          lqr(5, 1) = 0.50d+00 * (cc12 * x12 + cc23 * x23)
C--------
          d5_b = 0.50d+00 * x23
          d6_b = 0.50d+00 * ss23
          d7_b = 0.50d+00 * x12
          d8_b = 0.50d+00 * ss12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 5, 2) = d6_b * g_x23(g_i_) + d5_b * g_ss23(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_ss12(g_i_)
          enddo
          lqr(5, 2) = 0.50d+00 * (ss12 * x12 + ss23 * x23)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 5, 3) = ss23 * g_y32(g_i_) + y32 * g_ss23(g_i_) 
     *+ ss12 * g_y21(g_i_) + y21 * g_ss12(g_i_)
          enddo
          lqr(5, 3) = ss12 * y21 + ss23 * y32
C--------
C
          d5_b = -0.50d+00 * y32
          d6_b = -0.50d+00 * cc23
          d7_b = -0.50d+00 * y21
          d8_b = -0.50d+00 * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 6, 1) = d6_b * g_y32(g_i_) + d5_b * g_cc23(g_i_)
     * + d8_b * g_y21(g_i_) + d7_b * g_cc12(g_i_)
          enddo
          lqr(6, 1) = -0.50d+00 * (cc12 * y21 + cc23 * y32)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 6, 2) = -0.50d+00 * g_lqr(g_i_, 5, 3)
          enddo
          lqr(6, 2) = -0.50d+00 * lqr(5, 3)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 6, 3) = -2.00d+00 * g_lqr(g_i_, 5, 1)
          enddo
          lqr(6, 3) = -2.00d+00 * lqr(5, 1)
C--------
C
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
C
          d5_b = 0.50d+00 * x31
          d6_b = 0.50d+00 * cc31
          d7_b = 0.50d+00 * x23
          d8_b = 0.50d+00 * cc23
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 8, 1) = d6_b * g_x31(g_i_) + d5_b * g_cc31(g_i_)
     * + d8_b * g_x23(g_i_) + d7_b * g_cc23(g_i_)
          enddo
          lqr(8, 1) = 0.50d+00 * (cc23 * x23 + cc31 * x31)
C--------
          d5_b = 0.50d+00 * x31
          d6_b = 0.50d+00 * ss31
          d7_b = 0.50d+00 * x23
          d8_b = 0.50d+00 * ss23
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 8, 2) = d6_b * g_x31(g_i_) + d5_b * g_ss31(g_i_)
     * + d8_b * g_x23(g_i_) + d7_b * g_ss23(g_i_)
          enddo
          lqr(8, 2) = 0.50d+00 * (ss23 * x23 + ss31 * x31)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 8, 3) = ss31 * g_y13(g_i_) + y13 * g_ss31(g_i_) 
     *+ ss23 * g_y32(g_i_) + y32 * g_ss23(g_i_)
          enddo
          lqr(8, 3) = ss23 * y32 + ss31 * y13
C--------
C
          d5_b = -0.50d+00 * y13
          d6_b = -0.50d+00 * cc31
          d7_b = -0.50d+00 * y32
          d8_b = -0.50d+00 * cc23
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 9, 1) = d6_b * g_y13(g_i_) + d5_b * g_cc31(g_i_)
     * + d8_b * g_y32(g_i_) + d7_b * g_cc23(g_i_)
          enddo
          lqr(9, 1) = -0.50d+00 * (cc23 * y32 + cc31 * y13)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 9, 2) = -0.50d+00 * g_lqr(g_i_, 8, 3)
          enddo
          lqr(9, 2) = -0.50d+00 * lqr(8, 3)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 9, 3) = -2.00d+00 * g_lqr(g_i_, 8, 1)
          enddo
          lqr(9, 3) = -2.00d+00 * lqr(8, 1)
C--------
C
        endif
C
C.....ASSEMBLE THE LOCAL MATRIX [L] AS [CLR]*[LLR] + [CQR]*[LQR]
C
        if (clr .eq. zero) then
          do 99992 j = 1, 3
            do 99993 i = 1, 9
              do g_i_ = 1, g_p_
                g_l(g_i_, i, j) = g_lqr(g_i_, i, j)
              enddo
              l(i, j) = lqr(i, j)
C--------
2002          continue
99993       continue
2001        continue
99992     continue
        endif
C
        if (cqr .eq. zero) then
          do 99990 j = 1, 3
            do 99991 i = 1, 9
              do g_i_ = 1, g_p_
                g_l(g_i_, i, j) = g_llr(g_i_, i, j)
              enddo
              l(i, j) = llr(i, j)
C--------
2004          continue
99991       continue
2003        continue
99990     continue
        endif
C
        if ((clr .ne. zero) .and. (cqr .ne. zero)) then
          do 99988 j = 1, 3
            do 99989 i = 1, 9
              do g_i_ = 1, g_p_
                g_l(g_i_, i, j) = cqr * g_lqr(g_i_, i, j) + clr * g_llr(
     *g_i_, i, j)
              enddo
              l(i, j) = clr * llr(i, j) + cqr * lqr(i, j)
C--------
2006          continue
99989       continue
2005        continue
99988     continue
        endif
C
C.....MULTIPLY MATRIX [L] BY THE SQUARE ROOT OF THE STIFFNESS FACTOR
C.....AND DIVIDE BY THE SQUARE ROOT OF THE AREA
C
        d2_v = f / area
        d2_b = -d2_v / area
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d2_b * g_area(g_i_)
        enddo
        d1_w = d2_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compbBB.f',
     +739)
        endif
        do g_i_ = 1, g_p_
          g_factor(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        factor = d2_v
C--------
C
        do 99986 j = 1, 3
          do 99987 i = 1, 9
            do g_i_ = 1, g_p_
              g_l(g_i_, i, j) = factor * g_l(g_i_, i, j) + l(i, j) * g_f
     *actor(g_i_)
            enddo
            l(i, j) = factor * l(i, j)
C--------
3002        continue
99987     continue
3001      continue
99986   continue
C
C.....ASSEMBLE THE OUTPUT STIFFNESS SUCH THAT [kbBB] = [L]*[db]*[L]^T
C
        do 99984 i = 1, 9
C
          col = colb(i)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(1, 2) * g_l(g_i_, i, 2) + l(i, 2) * g_db(g
     *_i_, 1, 2) + db(1, 1) * g_l(g_i_, i, 1) + l(i, 1) * g_db(g_i_, 1, 
     *1)
          enddo
          d1_w = db(1, 1) * l(i, 1) + db(1, 2) * l(i, 2)
          do g_i_ = 1, g_p_
            g_dblt1(g_i_) = db(1, 3) * g_l(g_i_, i, 3) + l(i, 3) * g_db(
     *g_i_, 1, 3) + g_d1_w(g_i_)
          enddo
          dblt1 = d1_w + db(1, 3) * l(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(2, 2) * g_l(g_i_, i, 2) + l(i, 2) * g_db(g
     *_i_, 2, 2) + db(2, 1) * g_l(g_i_, i, 1) + l(i, 1) * g_db(g_i_, 2, 
     *1)
          enddo
          d1_w = db(2, 1) * l(i, 1) + db(2, 2) * l(i, 2)
          do g_i_ = 1, g_p_
            g_dblt2(g_i_) = db(2, 3) * g_l(g_i_, i, 3) + l(i, 3) * g_db(
     *g_i_, 2, 3) + g_d1_w(g_i_)
          enddo
          dblt2 = d1_w + db(2, 3) * l(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(3, 2) * g_l(g_i_, i, 2) + l(i, 2) * g_db(g
     *_i_, 3, 2) + db(3, 1) * g_l(g_i_, i, 1) + l(i, 1) * g_db(g_i_, 3, 
     *1)
          enddo
          d1_w = db(3, 1) * l(i, 1) + db(3, 2) * l(i, 2)
          do g_i_ = 1, g_p_
            g_dblt3(g_i_) = db(3, 3) * g_l(g_i_, i, 3) + l(i, 3) * g_db(
     *g_i_, 3, 3) + g_d1_w(g_i_)
          enddo
          dblt3 = d1_w + db(3, 3) * l(i, 3)
C--------
C
          do 99985 j = 1, i
            row = rowb(j)
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = dblt2 * g_l(g_i_, j, 2) + l(j, 2) * g_dblt2
     *(g_i_) + dblt1 * g_l(g_i_, j, 1) + l(j, 1) * g_dblt1(g_i_) + g_kbb
     *b(g_i_, row, col)
            enddo
            d1_w = kbbb(row, col) + dblt1 * l(j, 1) + dblt2 * l(j, 2)
            do g_i_ = 1, g_p_
              g_kbbb(g_i_, row, col) = dblt3 * g_l(g_i_, j, 3) + l(j, 3)
     * * g_dblt3(g_i_) + g_d1_w(g_i_)
            enddo
            kbbb(row, col) = d1_w + dblt3 * l(j, 3)
C--------
            do g_i_ = 1, g_p_
              g_kbbb(g_i_, col, row) = g_kbbb(g_i_, row, col)
            enddo
            kbbb(col, row) = kbbb(row, col)
C--------
4002        continue
99985     continue
C
4001      continue
99984   continue
C
C.....OUTPUT THE MATRIX PRIOR TO ROTATION (FOR DEBUGGING ONLY)
C
C     open(unit=90,file="kbBB.m")
C     write(90,*) "kbBB=["
C     do 991 i=1,18
C 991 write(90,9) (kbBB(i,j),j=1,18)
C     write(90,*) "     ];"
C     close(90)
C   9 format(18(1x,E16.9))
C
C.....ROTATE THE OUTPUT STIFFNESS MATRIX
C
C      call compmrot( kbBB , rot , rot , rot )
C
C     ------
C     RETURN
C     ------
C
        return
C
C     ---------------
C     ERROR TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF [CLR]+[CQR] IS DIFFERENT FROM ONE
C
100     continue
        write (*, *) '*** FATAL ERROR in routine COMPBBB      ***'
        write (*, *) '*** The Factors [clr] and [cqr] Violate ***'
        write (*, *) '*** the Constraint [clr]+[cqr]=1:       ***'
        write (*, *) '*** Check the Calling Sequence          ***'
        write (*, *) '*** EXECUTION TERNINATED RIGHT HERE     ***'
        stop
C
C.....ERROR-MESSAGE IF THE STIFFNESS FACTOR [F] IS NEGATIVE
C
200     continue
        write (*, *) '*** FATAL ERROR in routine COMPBBB       ***'
        write (*, *) '*** The Stiffness Factor [f] is Negative ***'
        write (*, *) '*** Check the Calling Sequence:          ***'
        write (*, *) '*** Factor [f] Must be Positive or Zero  ***'
        write (*, *) '*** EXECUTION TERNINATED RIGHT HERE      ***'
        stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
300     continue
        write (*, *) '*** FATAL ERROR in routine COMPBBB         ***'
        write (*, *) '*** The Triangle Area is Found Negative or ***'
        write (*, *) '*** Zero: Check the Nodal Point Numbering  ***'
        write (*, *) '*** ... Counterclockwise?                  ***'
        write (*, *) '*** EXECUTION TERNINATED RIGHT HERE        ***'
        stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPBBB"
C========================C
