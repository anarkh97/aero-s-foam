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
      subroutine g_compbmb(g_p_, elm, type, x, g_x, ldg_x, y, g_y, ldg_y
     *, dmb, g_dmb, ldg_dmb, fb, clr, cqr, alpha, fm, rowm, colb, rot, g
     *_rot, ldg_rot, l, g_l, ldg_l, p, g_p, ldg_p, fastcal, kbmb, g_kbmb
     *, ldg_kbmb)
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
        integer elm, type, rowm(9), colb(9)
        real*8 x(3), y(3), dmb(3, 3), kbmb(18, 18), rot(6, 6)
        real*8 fb, clr, cqr, alpha, fm, l(9, 3), p(9, 3)
        logical fastcal
C
C.....Local Variables
C
        integer i, j, row, col, dimp
        real*8 zero, twicearea, llr(9, 3), lqr(9, 3)
        real*8 x21, x32, x13, y21, y32, y13
        real*8 x12, x23, x31, y12, y23, y31
        real*8 x0, y0, dist12, dist23, dist31
        real*8 c12, c23, c31, s12, s23, s31
        real*8 cc12, cc23, cc31, ss12, ss23, ss31
        real*8 one, cs12, cs23, cs31, area, factorm
        real*8 dmblt1, dmblt2, dmblt3, factorb
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_kbmb, ldg_l, ldg_p, ldg_x, ldg_y, ldg_dm
     *b, ldg_rot
        double precision d9_b, d4_v, d8_b, d4_b, d1_p, d1_w, d7_b, d6_b,
     * d2_v, d3_v
        double precision d5_b, d2_b, d3_b, g_kbmb(ldg_kbmb, 18, 18), g_l
     *lr(g_pmax_, 9, 3), g_lqr(g_pmax_, 9, 3), g_l(ldg_l, 9, 3), g_p(ldg
     *_p, 9, 3), g_x21(g_pmax_), g_x(ldg_x, 3)
        double precision g_x12(g_pmax_), g_x32(g_pmax_), g_x23(g_pmax_),
     * g_x13(g_pmax_), g_x31(g_pmax_), g_y21(g_pmax_), g_y(ldg_y, 3), g_
     *y12(g_pmax_), g_y32(g_pmax_), g_y23(g_pmax_)
        double precision g_y13(g_pmax_), g_y31(g_pmax_), g_twicearea(g_p
     *max_), g_area(g_pmax_), g_d1_w(g_pmax_), g_dist12(g_pmax_), g_dist
     *23(g_pmax_), g_dist31(g_pmax_), g_c12(g_pmax_), g_s12(g_pmax_)
        double precision g_c23(g_pmax_), g_s23(g_pmax_), g_c31(g_pmax_),
     * g_s31(g_pmax_), g_cc12(g_pmax_), g_cc23(g_pmax_), g_cc31(g_pmax_)
     *, g_ss12(g_pmax_), g_ss23(g_pmax_), g_ss31(g_pmax_)
        double precision g_cs12(g_pmax_), g_cs23(g_pmax_), g_cs31(g_pmax
     *_), g_factorb(g_pmax_), g_factorm(g_pmax_), g_dmb(ldg_dmb, 3, 3), 
     *g_dmblt1(g_pmax_), g_dmblt2(g_pmax_), g_dmblt3(g_pmax_), g_rot(ldg
     *_rot, 6, 6)
        save g_ss23, g_ss31, g_cs12, g_cs23, g_cs31, g_factorb, g_factor
     *m, g_dmblt1, g_dmblt2, g_dmblt3
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
C     --------------------------------
C
C     BASIC CHECKS AND INITIALIZATIONS
C
C     --------------------------------
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
              g_kbmb(g_i_, i, j) = 0.0d0
            enddo
            kbmb(i, j) = zero
C--------
1002        continue
99999     continue
1001      continue
99998   continue
C
C.....CLEAR THE LOCAL MATRICES FOR BENDING
C
        do 99995 j = 1, 3
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
1003      continue
99995   continue
C
C.....INITIALIZE THE INTEGRATED STRAIN-TO-DISPLACEMENT AND
C.....CURVATURE-TO-DISPLACEMENT MATRICES TO ZERO OR SKIP
C.....THEIR ASSEMBLY IF THEY ARE AVAILABLE ALREADY
C
        if (fastcal) then
C
          goto 900
C
        else
C
          do 99992 j = 1, 3
            do 99994 i = 1, 9
              do g_i_ = 1, g_p_
                g_l(g_i_, i, j) = 0.0d0
              enddo
              l(i, j) = zero
C--------
1007          continue
99994       continue
            do 99993 i = 1, 9
              do g_i_ = 1, g_p_
                g_p(g_i_, i, j) = 0.0d0
              enddo
              p(i, j) = zero
C--------
1008          continue
99993       continue
1006        continue
99992     continue
C
        endif
C
C.....RETURN IF THE STIFFNESS FACTOR FOR BENDING IS ZERO
C
        if (fb .eq. zero) then
          return
        endif
C
C.....RETURN IF THE STIFFNESS FACTOR FOR MEMBRANE IS ZERO
C
        if (fm .eq. zero) then
          return
        endif
C
C.....CHECK IF [CLR] AND [CQR] SATISFY THE CONSTRAINT [CLR]+[CQR]=1
C
        if ((clr + cqr) .ne. one) then
          goto 100
        endif
C
C.....CHECK IF THE STIFFNESS FACTOR [FB] FOR BENDING IS POSITIVE
C
        if (fb .lt. zero) then
          goto 200
        endif
C
C.....CHECK IF THE STIFFNESS FACTOR [FM] FOR MEMBRANE IS POSITIVE
C
        if (fm .lt. zero) then
          goto 300
        endif
C
C     ------------------------------------------------------
C
C     FORM THE LOCAL MOMENT-CURVATURE MATRIX [L] FOR BENDING
C
C     ------------------------------------------------------
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
          goto 400
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
     +'g_compbMB.f',
     +314)
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
     +'g_compbMB.f',
     +333)
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
     +'g_compbMB.f',
     +352)
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
          do 99990 j = 1, 3
            do 99991 i = 1, 9
              do g_i_ = 1, g_p_
                g_l(g_i_, i, j) = g_lqr(g_i_, i, j)
              enddo
              l(i, j) = lqr(i, j)
C--------
2002          continue
99991       continue
2001        continue
99990     continue
        endif
C
        if (cqr .eq. zero) then
          do 99988 j = 1, 3
            do 99989 i = 1, 9
              do g_i_ = 1, g_p_
                g_l(g_i_, i, j) = g_llr(g_i_, i, j)
              enddo
              l(i, j) = llr(i, j)
C--------
2004          continue
99989       continue
2003        continue
99988     continue
        endif
C
        if ((clr .ne. zero) .and. (cqr .ne. zero)) then
          do 99986 j = 1, 3
            do 99987 i = 1, 9
              do g_i_ = 1, g_p_
                g_l(g_i_, i, j) = cqr * g_lqr(g_i_, i, j) + clr * g_llr(
     *g_i_, i, j)
              enddo
              l(i, j) = clr * llr(i, j) + cqr * lqr(i, j)
C--------
2006          continue
99987       continue
2005        continue
99986     continue
        endif
C
C.....MULTIPLY MATRIX [L] BY THE SQUARE ROOT OF THE STIFFNESS FACTOR
C.....(FOR BENDING) AND DIVIDE BY THE SQUARE ROOT OF THE AREA
C
        d2_v = fb / area
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
     +'g_compbMB.f',
     +791)
        endif
        do g_i_ = 1, g_p_
          g_factorb(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        factorb = d2_v
C--------
C
        do 99984 j = 1, 3
          do 99985 i = 1, 9
            do g_i_ = 1, g_p_
              g_l(g_i_, i, j) = factorb * g_l(g_i_, i, j) + l(i, j) * g_
     *factorb(g_i_)
            enddo
            l(i, j) = factorb * l(i, j)
C--------
3002        continue
99985     continue
3001      continue
99984   continue
C
C     ----------------------------------------------------------
C
C     FORM THE LOCAL NORMAL FORCE-STRAIN MATRIX [P] FOR MEMBRANE
C
C     ----------------------------------------------------------
C
C.....ASSEMBLE THE MATRIX [P] W/ SHAPE FUNCTION DERIVATIVES
C
        do g_i_ = 1, g_p_
          g_p(g_i_, 1, 1) = g_y23(g_i_)
        enddo
        p(1, 1) = y23
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 2, 1) = 0.0d0
        enddo
        p(2, 1) = zero
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 3, 1) = g_y31(g_i_)
        enddo
        p(3, 1) = y31
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 4, 1) = 0.0d0
        enddo
        p(4, 1) = zero
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 5, 1) = g_y12(g_i_)
        enddo
        p(5, 1) = y12
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 6, 1) = 0.0d0
        enddo
        p(6, 1) = zero
C--------
C
        do g_i_ = 1, g_p_
          g_p(g_i_, 1, 2) = 0.0d0
        enddo
        p(1, 2) = zero
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 2, 2) = g_x32(g_i_)
        enddo
        p(2, 2) = x32
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 3, 2) = 0.0d0
        enddo
        p(3, 2) = zero
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 4, 2) = g_x13(g_i_)
        enddo
        p(4, 2) = x13
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 5, 2) = 0.0d0
        enddo
        p(5, 2) = zero
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 6, 2) = g_x21(g_i_)
        enddo
        p(6, 2) = x21
C--------
C
        do g_i_ = 1, g_p_
          g_p(g_i_, 1, 3) = g_x32(g_i_)
        enddo
        p(1, 3) = x32
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 2, 3) = g_y23(g_i_)
        enddo
        p(2, 3) = y23
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 3, 3) = g_x13(g_i_)
        enddo
        p(3, 3) = x13
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 4, 3) = g_y31(g_i_)
        enddo
        p(4, 3) = y31
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 5, 3) = g_x21(g_i_)
        enddo
        p(5, 3) = x21
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 6, 3) = g_y12(g_i_)
        enddo
        p(6, 3) = y12
C--------
C
        dimp = 6
C
        if (alpha .ne. zero) then
C
          d4_v = y13 - y21
          d3_b = 1.0d0 / 6.00d+00 * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * y23
          do g_i_ = 1, g_p_
            g_p(g_i_, 7, 1) = -d5_b * g_y21(g_i_) + d5_b * g_y13(g_i_) +
     * d4_b * g_y23(g_i_)
          enddo
          p(7, 1) = y23 * d4_v * alpha / 6.00d+00
C--------
          d4_v = x31 - x12
          d3_b = 1.0d0 / 6.00d+00 * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x32
          do g_i_ = 1, g_p_
            g_p(g_i_, 7, 2) = -d5_b * g_x12(g_i_) + d5_b * g_x31(g_i_) +
     * d4_b * g_x32(g_i_)
          enddo
          p(7, 2) = x32 * d4_v * alpha / 6.00d+00
C--------
          d3_b = 1.0d0 / 3.00d+00 * alpha
          d6_b = -d3_b * y21
          d7_b = -d3_b * x12
          d8_b = d3_b * y13
          d9_b = d3_b * x31
          do g_i_ = 1, g_p_
            g_p(g_i_, 7, 3) = d7_b * g_y21(g_i_) + d6_b * g_x12(g_i_) + 
     *d9_b * g_y13(g_i_) + d8_b * g_x31(g_i_)
          enddo
          p(7, 3) = (x31 * y13 - x12 * y21) * alpha / 3.00d+00
C--------
C
          d4_v = y21 - y32
          d3_b = 1.0d0 / 6.00d+00 * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * y31
          do g_i_ = 1, g_p_
            g_p(g_i_, 8, 1) = -d5_b * g_y32(g_i_) + d5_b * g_y21(g_i_) +
     * d4_b * g_y31(g_i_)
          enddo
          p(8, 1) = y31 * d4_v * alpha / 6.00d+00
C--------
          d4_v = x12 - x23
          d3_b = 1.0d0 / 6.00d+00 * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x13
          do g_i_ = 1, g_p_
            g_p(g_i_, 8, 2) = -d5_b * g_x23(g_i_) + d5_b * g_x12(g_i_) +
     * d4_b * g_x13(g_i_)
          enddo
          p(8, 2) = x13 * d4_v * alpha / 6.00d+00
C--------
          d3_b = 1.0d0 / 3.00d+00 * alpha
          d6_b = -d3_b * y32
          d7_b = -d3_b * x23
          d8_b = d3_b * y21
          d9_b = d3_b * x12
          do g_i_ = 1, g_p_
            g_p(g_i_, 8, 3) = d7_b * g_y32(g_i_) + d6_b * g_x23(g_i_) + 
     *d9_b * g_y21(g_i_) + d8_b * g_x12(g_i_)
          enddo
          p(8, 3) = (x12 * y21 - x23 * y32) * alpha / 3.00d+00
C--------
C
          d4_v = y32 - y13
          d3_b = 1.0d0 / 6.00d+00 * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * y12
          do g_i_ = 1, g_p_
            g_p(g_i_, 9, 1) = -d5_b * g_y13(g_i_) + d5_b * g_y32(g_i_) +
     * d4_b * g_y12(g_i_)
          enddo
          p(9, 1) = y12 * d4_v * alpha / 6.00d+00
C--------
          d4_v = x23 - x31
          d3_b = 1.0d0 / 6.00d+00 * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x21
          do g_i_ = 1, g_p_
            g_p(g_i_, 9, 2) = -d5_b * g_x31(g_i_) + d5_b * g_x23(g_i_) +
     * d4_b * g_x21(g_i_)
          enddo
          p(9, 2) = x21 * d4_v * alpha / 6.00d+00
C--------
          d3_b = 1.0d0 / 3.00d+00 * alpha
          d6_b = -d3_b * y13
          d7_b = -d3_b * x31
          d8_b = d3_b * y32
          d9_b = d3_b * x23
          do g_i_ = 1, g_p_
            g_p(g_i_, 9, 3) = d7_b * g_y13(g_i_) + d6_b * g_x31(g_i_) + 
     *d9_b * g_y32(g_i_) + d8_b * g_x23(g_i_)
          enddo
          p(9, 3) = (x23 * y32 - x31 * y13) * alpha / 3.00d+00
C--------
C
          dimp = 9
C
        endif
C
C.....MULTIPLY MATRIX [P] BY THE SQUARE ROOT OF THE STIFFNESS FACTOR
C.....(FOR MEMBRANE) AND DIVIDE BY THE SQUARE ROOT OF FOUR TIMES THE AREA
C
C
        d2_v = 2.00d+00 * twicearea
        d3_v = fm / d2_v
        d3_b = -d3_v / d2_v * 2.00d+00
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_twicearea(g_i_)
        enddo
        d1_w = d3_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compbMB.f',
     +1034)
        endif
        do g_i_ = 1, g_p_
          g_factorm(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        factorm = d2_v
C--------
C
        do 99982 j = 1, 3
          do 99983 i = 1, dimp
            do g_i_ = 1, g_p_
              g_p(g_i_, i, j) = factorm * g_p(g_i_, i, j) + p(i, j) * g_
     *factorm(g_i_)
            enddo
            p(i, j) = factorm * p(i, j)
C--------
4002        continue
99983     continue
4001      continue
99982   continue
C
C     -------------------------------------------------------------
C
C     STIFFNESS MATRIX ASSEMBLY FOR BASIC MEMBRANE-BENDING COUPLING
C
C     -------------------------------------------------------------
C
900     continue
C
C.....ASSEMBLE THE OUTPUT STIFFNESS SUCH THAT [kbMB] = [P]*[dmb]*[L]^T
C
        do 99980 i = 1, 9
C
          col = colb(i)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dmb(1, 2) * g_l(g_i_, i, 2) + l(i, 2) * g_dmb
     *(g_i_, 1, 2) + dmb(1, 1) * g_l(g_i_, i, 1) + l(i, 1) * g_dmb(g_i_,
     * 1, 1)
          enddo
          d1_w = dmb(1, 1) * l(i, 1) + dmb(1, 2) * l(i, 2)
          do g_i_ = 1, g_p_
            g_dmblt1(g_i_) = dmb(1, 3) * g_l(g_i_, i, 3) + l(i, 3) * g_d
     *mb(g_i_, 1, 3) + g_d1_w(g_i_)
          enddo
          dmblt1 = d1_w + dmb(1, 3) * l(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dmb(2, 2) * g_l(g_i_, i, 2) + l(i, 2) * g_dmb
     *(g_i_, 2, 2) + dmb(2, 1) * g_l(g_i_, i, 1) + l(i, 1) * g_dmb(g_i_,
     * 2, 1)
          enddo
          d1_w = dmb(2, 1) * l(i, 1) + dmb(2, 2) * l(i, 2)
          do g_i_ = 1, g_p_
            g_dmblt2(g_i_) = dmb(2, 3) * g_l(g_i_, i, 3) + l(i, 3) * g_d
     *mb(g_i_, 2, 3) + g_d1_w(g_i_)
          enddo
          dmblt2 = d1_w + dmb(2, 3) * l(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dmb(3, 2) * g_l(g_i_, i, 2) + l(i, 2) * g_dmb
     *(g_i_, 3, 2) + dmb(3, 1) * g_l(g_i_, i, 1) + l(i, 1) * g_dmb(g_i_,
     * 3, 1)
          enddo
          d1_w = dmb(3, 1) * l(i, 1) + dmb(3, 2) * l(i, 2)
          do g_i_ = 1, g_p_
            g_dmblt3(g_i_) = dmb(3, 3) * g_l(g_i_, i, 3) + l(i, 3) * g_d
     *mb(g_i_, 3, 3) + g_d1_w(g_i_)
          enddo
          dmblt3 = d1_w + dmb(3, 3) * l(i, 3)
C--------
C
          do 99981 j = 1, 9
            row = rowm(j)
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = dmblt2 * g_p(g_i_, j, 2) + p(j, 2) * g_dmbl
     *t2(g_i_) + dmblt1 * g_p(g_i_, j, 1) + p(j, 1) * g_dmblt1(g_i_) + g
     *_kbmb(g_i_, row, col)
            enddo
            d1_w = kbmb(row, col) + dmblt1 * p(j, 1) + dmblt2 * p(j, 2)
            do g_i_ = 1, g_p_
              g_kbmb(g_i_, row, col) = dmblt3 * g_p(g_i_, j, 3) + p(j, 3
     *) * g_dmblt3(g_i_) + g_d1_w(g_i_)
            enddo
            kbmb(row, col) = d1_w + dmblt3 * p(j, 3)
C--------
5002        continue
99981     continue
C
5001      continue
99980   continue
C
C.....OUTPUT THE MATRIX PRIOR TO ROTATION (FOR DEBUGGING ONLY)
C
C     open(unit=90,file="kbMB.m")
C     write(90,*) "kbMB=["
C     do 991 i=1,18
C 991 write(90,9) (kbMB(i,j),j=1,18)
C     write(90,*) "     ];"
C     close(90)
C   9 format(18(1x,E16.9))
C
C.....ROTATE THE OUTPUT STIFFNESS MATRIX
C
C      call compmrot( kbMB , rot , rot , rot )
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
        write (*, *) '*** FATAL ERROR in routine COMPBMB      ***'
        write (*, *) '*** The Factors [clr] and [cqr] Violate ***'
        write (*, *) '*** the Constraint [clr]+[cqr]=1:       ***'
        write (*, *) '*** Check the Calling Sequence          ***'
        write (*, *) '*** EXECUTION TERNINATED RIGHT HERE     ***'
        stop
C
C.....ERROR-MESSAGE IF THE STIFFNESS FACTOR [FB] (BENDING) IS NEGATIVE
C
200     continue
        write (*, *) '*** FATAL ERROR in routine COMPBMB          ***'
        write (*, *) '*** The Stiffness Factor [fb] for Bending   ***'
        write (*, *) '*** is Negative: Check the Calling Sequence ***'
        write (*, *) '*** Factor [fb] Must be Positive or Zero    ***'
        write (*, *) '*** EXECUTION TERNINATED RIGHT HERE         ***'
        stop
C
C.....ERROR-MESSAGE IF THE STIFFNESS FACTOR [FM] (MEMBRANE) IS NEGATIVE
C
300     continue
        write (*, *) '*** FATAL ERROR in routine COMPBMB          ***'
        write (*, *) '*** The Stiffness Factor [fm] for Membrane  ***'
        write (*, *) '*** is Negative: Check the Calling Sequence ***'
        write (*, *) '*** Factor [fm] Must be Positive or Zero    ***'
        write (*, *) '*** EXECUTION TERNINATED RIGHT HERE         ***'
        stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
400     continue
        write (*, *) '*** FATAL ERROR in routine COMPBMB         ***'
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
C=end of routine "COMPBMB"
C========================C
