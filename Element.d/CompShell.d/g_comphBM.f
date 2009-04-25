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
      subroutine g_comphbm(g_p_, elm, type, x, g_x, ldg_x, y, g_y, ldg_y
     *, dbm, g_dbm, ldg_dbm, fb, fm, rowb, colm, rot, g_rot, ldg_rot, lh
     *1, g_lh1, ldg_lh1, lh2, g_lh2, ldg_lh2, lh3, g_lh3, ldg_lh3, ph1, 
     *g_ph1, ldg_ph1, ph2, g_ph2, ldg_ph2, ph3, g_ph3, ldg_ph3, fastcal,
     * khbm, g_khbm, ldg_khbm)
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
        integer elm, type, rowb(9), colm(9)
        real*8 x(3), y(3), dbm(3, 3), khbm(18, 18), rot(6, 6)
        real*8 lh1(9, 3), lh2(9, 3), lh3(9, 3)
        real*8 fb, fm, ph1(9, 3), ph2(9, 3), ph3(9, 3)
        logical fastcal
C
C.....Local Variables
C
        integer i, j, col, row
        real*8 zero, x0, y0, x1, x2, x3, y1, y2, y3
        real*8 x21, x32, x13, y21, y32, y13
        real*8 area, twicearea, x2ap3, dist21, dist32, dist13
        real*8 al1, al2, al3, bl1, bl2, bl3, factor
        real*8 s1, s2, s3, q1, q2, q3, q4, q5, q6
        real*8 q(6, 9), sq(3, 3), rm1(3, 6), rm2(3, 6), rm3(3, 6)
        real*8 rm1tsqt(6, 3), rm2tsqt(6, 3), rm3tsqt(6, 3)
        real*8 gauss11, gauss21, gauss31, invweight1
        real*8 gauss12, gauss22, gauss32, invweight2
        real*8 gauss13, gauss23, gauss33, invweight3
        real*8 mult11, mult12, mult13
        real*8 mult21, mult22, mult23
        real*8 mult31, mult32, mult33
        real*8 x12, x23, x31, y12, y23, y31
        real*8 x10, x20, x30, y10, y20, y30
        real*8 aa12, aa23, aa31, ss12, ss23, ss31
        real*8 ca, caa12, caa23, caa31, sum123, sum456
        real*8 cax10, cax20, cax30, cay10, cay20, cay30
        real*8 h1, h2, h3, h4, h5, h6, factor1, factor2
        real*8 hmv(6, 9), hqv(6, 9), b(3, 6), factor3
        real*8 z1(6), z2(6), z3(6)
        real*8 z1bt(6, 3), z2bt(6, 3), z3bt(6, 3)
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_khbm, ldg_lh1, ldg_lh2, ldg_lh3, ldg_ph1
     *, ldg_ph2, ldg_ph3, ldg_x
        integer ldg_y, ldg_dbm, ldg_rot
        double precision d4_w, d3_w, d6_b, d1_p, d7_b, d1_w, d8_b, d8_v,
     * d2_v, d3_v
        double precision d4_v, d5_v, d6_v, d2_w, d2_b, d3_b, d4_b, d5_b,
     * g_khbm(ldg_khbm, 18, 18), g_q(g_pmax_, 6, 9)
        double precision g_sq(g_pmax_, 3, 3), g_rm1(g_pmax_, 3, 6), g_rm
     *2(g_pmax_, 3, 6), g_rm3(g_pmax_, 3, 6), g_rm1tsqt(g_pmax_, 6, 3), 
     *g_rm2tsqt(g_pmax_, 6, 3), g_rm3tsqt(g_pmax_, 6, 3), g_hmv(g_pmax_,
     * 6, 9), g_hqv(g_pmax_, 6, 9), g_b(g_pmax_, 3, 6)
        double precision g_z1bt(g_pmax_, 6, 3), g_z2bt(g_pmax_, 6, 3), g
     *_z3bt(g_pmax_, 6, 3), g_lh1(ldg_lh1, 9, 3), g_lh2(ldg_lh2, 9, 3), 
     *g_lh3(ldg_lh3, 9, 3), g_ph1(ldg_ph1, 9, 3), g_ph2(ldg_ph2, 9, 3), 
     *g_ph3(ldg_ph3, 9, 3), g_x0(g_pmax_)
        double precision g_x(ldg_x, 3), g_y0(g_pmax_), g_y(ldg_y, 3), g_
     *x1(g_pmax_), g_x2(g_pmax_), g_x3(g_pmax_), g_y1(g_pmax_), g_y2(g_p
     *max_), g_y3(g_pmax_), g_x21(g_pmax_)
        double precision g_x32(g_pmax_), g_x13(g_pmax_), g_y21(g_pmax_),
     * g_y32(g_pmax_), g_y13(g_pmax_), g_x12(g_pmax_), g_x23(g_pmax_), g
     *_x31(g_pmax_), g_y12(g_pmax_), g_y23(g_pmax_)
        double precision g_y31(g_pmax_), g_twicearea(g_pmax_), g_area(g_
     *pmax_), g_d1_w(g_pmax_), g_dist21(g_pmax_), g_dist32(g_pmax_), g_d
     *ist13(g_pmax_), g_d2_w(g_pmax_), g_bl1(g_pmax_), g_bl2(g_pmax_)
        double precision g_bl3(g_pmax_), g_al1(g_pmax_), g_al2(g_pmax_),
     * g_al3(g_pmax_), g_x2ap3(g_pmax_), g_factor(g_pmax_), g_s1(g_pmax_
     *), g_s2(g_pmax_), g_s3(g_pmax_), g_q1(g_pmax_)
        double precision g_q2(g_pmax_), g_q3(g_pmax_), g_q4(g_pmax_), g_
     *q5(g_pmax_), g_q6(g_pmax_), g_x10(g_pmax_), g_x20(g_pmax_), g_x30(
     *g_pmax_), g_y10(g_pmax_), g_y20(g_pmax_)
        double precision g_y30(g_pmax_), g_aa12(g_pmax_), g_aa23(g_pmax_
     *), g_aa31(g_pmax_), g_caa12(g_pmax_), g_caa23(g_pmax_), g_caa31(g_
     *pmax_), g_ss12(g_pmax_), g_ss23(g_pmax_), g_ss31(g_pmax_)
        double precision g_ca(g_pmax_), g_cax10(g_pmax_), g_cax20(g_pmax
     *_), g_cax30(g_pmax_), g_cay10(g_pmax_), g_cay20(g_pmax_), g_cay30(
     *g_pmax_), g_sum123(g_pmax_), g_sum456(g_pmax_), g_factor1(g_pmax_)
        double precision g_factor2(g_pmax_), g_factor3(g_pmax_), g_h1(g_
     *pmax_), g_h2(g_pmax_), g_h3(g_pmax_), g_h4(g_pmax_), g_h5(g_pmax_)
     *, g_h6(g_pmax_), g_dbm(ldg_dbm, 3, 3), g_mult11(g_pmax_)
        double precision g_mult12(g_pmax_), g_mult13(g_pmax_), g_mult21(
     *g_pmax_), g_mult22(g_pmax_), g_mult23(g_pmax_), g_mult31(g_pmax_),
     * g_mult32(g_pmax_), g_mult33(g_pmax_), g_d3_w(g_pmax_), g_d4_w(g_p
     *max_)
        double precision g_rot(ldg_rot, 6, 6)
        save g_d3_w, g_d4_w
        save g_h6, g_mult11, g_mult12, g_mult13, g_mult21, g_mult22, g_m
     *ult23, g_mult31, g_mult32, g_mult33
        save g_sum123, g_sum456, g_factor1, g_factor2, g_factor3, g_h1, 
     *g_h2, g_h3, g_h4, g_h5
        save g_ss12, g_ss23, g_ss31, g_ca, g_cax10, g_cax20, g_cax30, g_
     *cay10, g_cay20, g_cay30
        save g_x30, g_y10, g_y20, g_y30, g_aa12, g_aa23, g_aa31, g_caa12
     *, g_caa23, g_caa31
        save g_s2, g_s3, g_q1, g_q2, g_q3, g_q4, g_q5, g_q6, g_x10, g_x2
     *0
        save g_d2_w, g_bl1, g_bl2, g_bl3, g_al1, g_al2, g_al3, g_x2ap3, 
     *g_factor, g_s1
        save g_x31, g_y12, g_y23, g_y31, g_twicearea, g_area, g_d1_w, g_
     *dist21, g_dist32, g_dist13
        save g_y2, g_y3, g_x21, g_x32, g_x13, g_y21, g_y32, g_y13, g_x12
     *, g_x23
        save g_b, g_z1bt, g_z2bt, g_z3bt, g_x0, g_y0, g_x1, g_x2, g_x3, 
     *g_y1
        save g_q, g_sq, g_rm1, g_rm2, g_rm3, g_rm1tsqt, g_rm2tsqt, g_rm3
     *tsqt, g_hmv, g_hqv
        data zero /0.000000d+00/
C
C.....DEFINE THE FIRST SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
        data gauss11 /0.000000d+00/
        data gauss21 /0.500000d+00/
        data gauss31 /0.500000d+00/
        data invweight1 /3.000000d+00/
C
C.....DEFINE THE SECOND SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
        data gauss12 /0.500000d+00/
        data gauss22 /0.000000d+00/
        data gauss32 /0.500000d+00/
        data invweight2 /3.000000d+00/
C
C.....DEFINE THE THIRD SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
        data gauss13 /0.500000d+00/
        data gauss23 /0.500000d+00/
        data gauss33 /0.000000d+00/
        data invweight3 /3.000000d+00/
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
              g_khbm(g_i_, i, j) = 0.0d0
            enddo
            khbm(i, j) = zero
C--------
1002        continue
99999     continue
1001      continue
99998   continue
C
C.....CLEAR THE LOCAL MATRICES FOR BENDING
C
        do 99996 j = 1, 9
          do 99997 i = 1, 6
            do g_i_ = 1, g_p_
              g_q(g_i_, i, j) = 0.0d0
            enddo
            q(i, j) = zero
C--------
1004        continue
99997     continue
1003      continue
99996   continue
C
        do 99994 j = 1, 3
          do 99995 i = 1, 3
            do g_i_ = 1, g_p_
              g_sq(g_i_, i, j) = 0.0d0
            enddo
            sq(i, j) = zero
C--------
1006        continue
99995     continue
1005      continue
99994   continue
C
        do 99990 j = 1, 6
          do 99993 i = 1, 3
            do g_i_ = 1, g_p_
              g_rm1(g_i_, i, j) = 0.0d0
            enddo
            rm1(i, j) = zero
C--------
1008        continue
99993     continue
          do 99992 i = 1, 3
            do g_i_ = 1, g_p_
              g_rm2(g_i_, i, j) = 0.0d0
            enddo
            rm2(i, j) = zero
C--------
1009        continue
99992     continue
          do 99991 i = 1, 3
            do g_i_ = 1, g_p_
              g_rm3(g_i_, i, j) = 0.0d0
            enddo
            rm3(i, j) = zero
C--------
1010        continue
99991     continue
1007      continue
99990   continue
C
        do 99986 j = 1, 3
          do 99989 i = 1, 6
            do g_i_ = 1, g_p_
              g_rm1tsqt(g_i_, i, j) = 0.0d0
            enddo
            rm1tsqt(i, j) = zero
C--------
1012        continue
99989     continue
          do 99988 i = 1, 6
            do g_i_ = 1, g_p_
              g_rm2tsqt(g_i_, i, j) = 0.0d0
            enddo
            rm2tsqt(i, j) = zero
C--------
1013        continue
99988     continue
          do 99987 i = 1, 6
            do g_i_ = 1, g_p_
              g_rm3tsqt(g_i_, i, j) = 0.0d0
            enddo
            rm3tsqt(i, j) = zero
C--------
1014        continue
99987     continue
1011      continue
99986   continue
C
C.....CLEAR THE LOCAL MATRICES FOR MEMBRANE
C
        do 99983 j = 1, 9
          do 99985 i = 1, 6
            do g_i_ = 1, g_p_
              g_hmv(g_i_, i, j) = 0.0d0
            enddo
            hmv(i, j) = zero
C--------
1016        continue
99985     continue
          do 99984 i = 1, 6
            do g_i_ = 1, g_p_
              g_hqv(g_i_, i, j) = 0.0d0
            enddo
            hqv(i, j) = zero
C--------
1017        continue
99984     continue
1015      continue
99983   continue
C
        do 99981 j = 1, 6
          do 99982 i = 1, 3
            do g_i_ = 1, g_p_
              g_b(g_i_, i, j) = 0.0d0
            enddo
            b(i, j) = zero
C--------
1019        continue
99982     continue
1018      continue
99981   continue
C
        do 99980 i = 1, 6
          z1(i) = zero
1020      continue
99980   continue
C
        do 99979 i = 1, 6
          z2(i) = zero
1021      continue
99979   continue
C
        do 99978 i = 1, 6
          z3(i) = zero
1022      continue
99978   continue
C
        do 99974 j = 1, 3
          do 99977 i = 1, 6
            do g_i_ = 1, g_p_
              g_z1bt(g_i_, i, j) = 0.0d0
            enddo
            z1bt(i, j) = zero
C--------
1024        continue
99977     continue
          do 99976 i = 1, 6
            do g_i_ = 1, g_p_
              g_z2bt(g_i_, i, j) = 0.0d0
            enddo
            z2bt(i, j) = zero
C--------
1025        continue
99976     continue
          do 99975 i = 1, 6
            do g_i_ = 1, g_p_
              g_z3bt(g_i_, i, j) = 0.0d0
            enddo
            z3bt(i, j) = zero
C--------
1026        continue
99975     continue
1023      continue
99974   continue
C
C.....INITIALIZE THE INTEGRATED STRAIN-TO-DISPLACEMENT AND
C.....CURVATURE-TO-DISPLACEMENT MATRICES TI ZERO OR SKIP
C.....THEIR ASSEMBLY IF THEY ARE AVAILABLE ALREADY
C
        if (fastcal) then
C
          goto 900
C
        else
C
          do 99970 j = 1, 3
            do 99973 i = 1, 9
              do g_i_ = 1, g_p_
                g_lh1(g_i_, i, j) = 0.0d0
              enddo
              lh1(i, j) = zero
C--------
1028          continue
99973       continue
            do 99972 i = 1, 9
              do g_i_ = 1, g_p_
                g_lh2(g_i_, i, j) = 0.0d0
              enddo
              lh2(i, j) = zero
C--------
1029          continue
99972       continue
            do 99971 i = 1, 9
              do g_i_ = 1, g_p_
                g_lh3(g_i_, i, j) = 0.0d0
              enddo
              lh3(i, j) = zero
C--------
1030          continue
99971       continue
1027        continue
99970     continue
C
          do 99966 j = 1, 3
            do 99969 i = 1, 9
              do g_i_ = 1, g_p_
                g_ph1(g_i_, i, j) = 0.0d0
              enddo
              ph1(i, j) = zero
C--------
1032          continue
99969       continue
            do 99968 i = 1, 9
              do g_i_ = 1, g_p_
                g_ph2(g_i_, i, j) = 0.0d0
              enddo
              ph2(i, j) = zero
C--------
1033          continue
99968       continue
            do 99967 i = 1, 9
              do g_i_ = 1, g_p_
                g_ph3(g_i_, i, j) = 0.0d0
              enddo
              ph3(i, j) = zero
C--------
1034          continue
99967       continue
1031        continue
99966     continue
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
C.....CHECK IF THE STIFFNESS FACTOR [FB] FOR BENDING IS POSITIVE
C
        if (fb .lt. zero) then
          goto 100
        endif
C
C.....CHECK IF THE STIFFNESS FACTOR [FM] FOR MEMBRANE IS POSITIVE
C
        if (fm .lt. zero) then
          goto 200
        endif
C
C     ------------------------------------------------------
C
C     FORM THE LOCAL MOMENT-CURVATURE MATRICES FOR BENDING
C
C             AND FOR THE THREE INTEGRATION POINTS
C
C     ------------------------------------------------------
C
C.....GET THE COORDINATES OF THE CENTROID OF THE TRIANGLE
C
        d2_b = 1.0d0 / 3.00d+00
        do g_i_ = 1, g_p_
          g_x0(g_i_) = d2_b * g_x(g_i_, 3) + d2_b * g_x(g_i_, 2) + d2_b 
     ** g_x(g_i_, 1)
        enddo
        x0 = (x(1) + x(2) + x(3)) / 3.00d+00
C--------
        d2_b = 1.0d0 / 3.00d+00
        do g_i_ = 1, g_p_
          g_y0(g_i_) = d2_b * g_y(g_i_, 3) + d2_b * g_y(g_i_, 2) + d2_b 
     ** g_y(g_i_, 1)
        enddo
        y0 = (y(1) + y(2) + y(3)) / 3.00d+00
C--------
C
C.....GET THE DISTANCES BETWEEN NODES 1, 2 AND 3 AND THE CENTROID
C
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
C
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
C
C.....GET THE DISTANCES BETWEEN NODAL POINT X- AND Y- COORDINATES
C
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
C
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
C
        do g_i_ = 1, g_p_
          g_x12(g_i_) = -g_x21(g_i_)
        enddo
        x12 = -x21
C--------
        do g_i_ = 1, g_p_
          g_x23(g_i_) = -g_x32(g_i_)
        enddo
        x23 = -x32
C--------
        do g_i_ = 1, g_p_
          g_x31(g_i_) = -g_x13(g_i_)
        enddo
        x31 = -x13
C--------
C
        do g_i_ = 1, g_p_
          g_y12(g_i_) = -g_y21(g_i_)
        enddo
        y12 = -y21
C--------
        do g_i_ = 1, g_p_
          g_y23(g_i_) = -g_y32(g_i_)
        enddo
        y23 = -y32
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
C.....GET THE DISTANCES BETWEEN NODES 1-2, 2-3 AND 3-1
C
        d4_b = y21 + y21
        d5_b = x21 + x21
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y21(g_i_) + d5_b * g_x21(g_i_)
        enddo
        d1_w = x21 * x21 + y21 * y21
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_comphBM.f',
     +605)
        endif
        do g_i_ = 1, g_p_
          g_dist21(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        dist21 = d2_v
C--------
        d4_b = y32 + y32
        d5_b = x32 + x32
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y32(g_i_) + d5_b * g_x32(g_i_)
        enddo
        d1_w = x32 * x32 + y32 * y32
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_comphBM.f',
     +624)
        endif
        do g_i_ = 1, g_p_
          g_dist32(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        dist32 = d2_v
C--------
        d4_b = y13 + y13
        d5_b = x13 + x13
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y13(g_i_) + d5_b * g_x13(g_i_)
        enddo
        d1_w = x13 * x13 + y13 * y13
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_comphBM.f',
     +643)
        endif
        do g_i_ = 1, g_p_
          g_dist13(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        dist13 = d2_v
C--------
C
C.....GET THE SIDE PROJECTIONS
C
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
        d3_v = dist13 * dist13
        d4_v = d2_w / d3_v
        d2_b = 1.0d0 / d3_v
        d3_b = -d4_v / d3_v
        d4_b = d3_b * dist13 + d3_b * dist13
        do g_i_ = 1, g_p_
          g_bl1(g_i_) = d4_b * g_dist13(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        bl1 = d4_v
C--------
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
        d3_v = dist21 * dist21
        d4_v = d2_w / d3_v
        d2_b = 1.0d0 / d3_v
        d3_b = -d4_v / d3_v
        d4_b = d3_b * dist21 + d3_b * dist21
        do g_i_ = 1, g_p_
          g_bl2(g_i_) = d4_b * g_dist21(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        bl2 = d4_v
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
        d3_v = dist32 * dist32
        d4_v = d2_w / d3_v
        d2_b = 1.0d0 / d3_v
        d3_b = -d4_v / d3_v
        d4_b = d3_b * dist32 + d3_b * dist32
        do g_i_ = 1, g_p_
          g_bl3(g_i_) = d4_b * g_dist32(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        bl3 = d4_v
C--------
C
        do g_i_ = 1, g_p_
          g_al1(g_i_) = -g_bl1(g_i_)
        enddo
        al1 = 1.00d+00 - bl1
C--------
        do g_i_ = 1, g_p_
          g_al2(g_i_) = -g_bl2(g_i_)
        enddo
        al2 = 1.00d+00 - bl2
C--------
        do g_i_ = 1, g_p_
          g_al3(g_i_) = -g_bl3(g_i_)
        enddo
        al3 = 1.00d+00 - bl3
C--------
C
C.....FORM THE MATRIX [Q] W/ SHAPE FUNCTION DERIVATIVES
C
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 1) = 0.0d0
        enddo
        q(1, 1) = 6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 2) = -2.00d+00 * g_y13(g_i_)
        enddo
        q(1, 2) = -2.00d+00 * y13
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 3) = 2.00d+00 * g_x13(g_i_)
        enddo
        q(1, 3) = 2.00d+00 * x13
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 7) = 0.0d0
        enddo
        q(1, 7) = -6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 8) = -4.00d+00 * g_y13(g_i_)
        enddo
        q(1, 8) = -4.00d+00 * y13
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 9) = 4.00d+00 * g_x13(g_i_)
        enddo
        q(1, 9) = 4.00d+00 * x13
C--------
C
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 1) = 0.0d0
        enddo
        q(2, 1) = -6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 2) = 4.00d+00 * g_y13(g_i_)
        enddo
        q(2, 2) = 4.00d+00 * y13
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 3) = -4.00d+00 * g_x13(g_i_)
        enddo
        q(2, 3) = -4.00d+00 * x13
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 7) = 0.0d0
        enddo
        q(2, 7) = 6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 8) = 2.00d+00 * g_y13(g_i_)
        enddo
        q(2, 8) = 2.00d+00 * y13
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 9) = -2.00d+00 * g_x13(g_i_)
        enddo
        q(2, 9) = -2.00d+00 * x13
C--------
C
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 1) = 0.0d0
        enddo
        q(3, 1) = -6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 2) = -4.00d+00 * g_y21(g_i_)
        enddo
        q(3, 2) = -4.00d+00 * y21
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 3) = 4.00d+00 * g_x21(g_i_)
        enddo
        q(3, 3) = 4.00d+00 * x21
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 4) = 0.0d0
        enddo
        q(3, 4) = 6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 5) = -2.00d+00 * g_y21(g_i_)
        enddo
        q(3, 5) = -2.00d+00 * y21
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 6) = 2.00d+00 * g_x21(g_i_)
        enddo
        q(3, 6) = 2.00d+00 * x21
C--------
C
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 1) = 0.0d0
        enddo
        q(4, 1) = 6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 2) = 2.00d+00 * g_y21(g_i_)
        enddo
        q(4, 2) = 2.00d+00 * y21
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 3) = -2.00d+00 * g_x21(g_i_)
        enddo
        q(4, 3) = -2.00d+00 * x21
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 4) = 0.0d0
        enddo
        q(4, 4) = -6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 5) = 4.00d+00 * g_y21(g_i_)
        enddo
        q(4, 5) = 4.00d+00 * y21
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 4, 6) = -4.00d+00 * g_x21(g_i_)
        enddo
        q(4, 6) = -4.00d+00 * x21
C--------
C
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 4) = 0.0d0
        enddo
        q(5, 4) = -6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 5) = -4.00d+00 * g_y32(g_i_)
        enddo
        q(5, 5) = -4.00d+00 * y32
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 6) = 4.00d+00 * g_x32(g_i_)
        enddo
        q(5, 6) = 4.00d+00 * x32
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 7) = 0.0d0
        enddo
        q(5, 7) = 6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 8) = -2.00d+00 * g_y32(g_i_)
        enddo
        q(5, 8) = -2.00d+00 * y32
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 5, 9) = 2.00d+00 * g_x32(g_i_)
        enddo
        q(5, 9) = 2.00d+00 * x32
C--------
C
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 4) = 0.0d0
        enddo
        q(6, 4) = 6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 5) = 2.00d+00 * g_y32(g_i_)
        enddo
        q(6, 5) = 2.00d+00 * y32
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 6) = -2.00d+00 * g_x32(g_i_)
        enddo
        q(6, 6) = -2.00d+00 * x32
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 7) = 0.0d0
        enddo
        q(6, 7) = -6.00d+00
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 8) = 4.00d+00 * g_y32(g_i_)
        enddo
        q(6, 8) = 4.00d+00 * y32
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 6, 9) = -4.00d+00 * g_x32(g_i_)
        enddo
        q(6, 9) = -4.00d+00 * x32
C--------
C
C.....GET THE MATRIX [SQ] THAT REPRESENTS THE INVERSE OF THE MATRIX
C.....RELATING INSIDE CURVATURES WITH BOUNDARY CURVATURES
C
        d2_v = twicearea * twicearea
        d3_b = d2_v + twicearea * twicearea + twicearea * twicearea
        do g_i_ = 1, g_p_
          g_x2ap3(g_i_) = d3_b * g_twicearea(g_i_)
        enddo
        x2ap3 = d2_v * twicearea
C--------
C
        d4_v = -x21 * y21
        d6_v = d4_v * y32
        d4_b = y32 * y32
        d3_b = d6_v + y32 * d4_v
        d6_b = d4_b * (-x21)
        d7_b = -(d4_b * y21)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_y32(g_i_) + d6_b * g_y21(g_i_) + d7_b 
     ** g_x21(g_i_)
        enddo
        d1_w = d6_v * y32
        d4_v = x32 * y21
        d5_v = d4_v * y21
        d6_b = y32 * y21
        d8_b = d6_b * y21
        d7_b = y32 * d4_v + d6_b * x32
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d5_v * g_y32(g_i_) + d7_b * g_y21(g_i_) + d8_b 
     ** g_x32(g_i_) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + d5_v * y32
        d3_v = d2_w / x2ap3
        d2_b = 1.0d0 / x2ap3
        d3_b = -d3_v / x2ap3
        do g_i_ = 1, g_p_
          g_sq(g_i_, 1, 1) = d3_b * g_x2ap3(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        sq(1, 1) = d3_v
C--------
        d3_v = x13 * y13
        d5_v = d3_v * y32
        d4_b = y32 * y32
        d3_b = d5_v + y32 * d3_v
        d5_b = d4_b * y13
        d6_b = d4_b * x13
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_y32(g_i_) + d6_b * g_y13(g_i_) + d5_b 
     ** g_x13(g_i_)
        enddo
        d1_w = d5_v * y32
        d4_v = x32 * y13
        d5_v = d4_v * y13
        d6_b = -y32 * y13
        d8_b = d6_b * y13
        d7_b = -y32 * d4_v + d6_b * x32
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = -d5_v * g_y32(g_i_) + d7_b * g_y13(g_i_) + d8_b
     * * g_x32(g_i_) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w - d5_v * y32
        d3_v = d2_w / x2ap3
        d2_b = 1.0d0 / x2ap3
        d3_b = -d3_v / x2ap3
        do g_i_ = 1, g_p_
          g_sq(g_i_, 1, 2) = d3_b * g_x2ap3(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        sq(1, 2) = d3_v
C--------
        d3_v = x21 * y21
        d5_v = d3_v * y13
        d4_b = y13 * y13
        d3_b = d5_v + y13 * d3_v
        d5_b = d4_b * y21
        d6_b = d4_b * x21
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_y13(g_i_) + d6_b * g_y21(g_i_) + d5_b 
     ** g_x21(g_i_)
        enddo
        d1_w = d5_v * y13
        d4_v = x13 * y21
        d5_v = d4_v * y21
        d6_b = -y13 * y21
        d8_b = d6_b * y21
        d7_b = -y13 * d4_v + d6_b * x13
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = -d5_v * g_y13(g_i_) + d7_b * g_y21(g_i_) + d8_b
     * * g_x13(g_i_) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w - d5_v * y13
        d3_v = d2_w / x2ap3
        d2_b = 1.0d0 / x2ap3
        d3_b = -d3_v / x2ap3
        do g_i_ = 1, g_p_
          g_sq(g_i_, 1, 3) = d3_b * g_x2ap3(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        sq(1, 3) = d3_v
C--------
C
        d3_v = x21 * x32
        d4_v = d3_v * x32
        d4_b = y21 * x32
        d6_b = d4_b * x32
        d5_b = y21 * d3_v + d4_b * x21
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_v * g_y21(g_i_) + d5_b * g_x32(g_i_) + d6_b 
     ** g_x21(g_i_)
        enddo
        d1_w = d4_v * y21
        d3_v = x21 * x21
        d5_v = d3_v * x32
        d6_b = -y32 * x32
        d7_b = -y32 * d3_v
        d8_b = d6_b * x21 + d6_b * x21
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = -d5_v * g_y32(g_i_) + d7_b * g_x32(g_i_) + d8_b
     * * g_x21(g_i_) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w - d5_v * y32
        d3_v = d2_w / x2ap3
        d2_b = 1.0d0 / x2ap3
        d3_b = -d3_v / x2ap3
        do g_i_ = 1, g_p_
          g_sq(g_i_, 2, 1) = d3_b * g_x2ap3(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        sq(2, 1) = d3_v
C--------
        d4_v = -x13 * x32
        d5_v = d4_v * x32
        d4_b = y13 * x32
        d5_b = y13 * d4_v + d4_b * (-x13)
        d7_b = -(d4_b * x32)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d5_v * g_y13(g_i_) + d5_b * g_x32(g_i_) + d7_b 
     ** g_x13(g_i_)
        enddo
        d1_w = d5_v * y13
        d3_v = x13 * x13
        d5_v = d3_v * x32
        d6_b = y32 * x32
        d7_b = y32 * d3_v
        d8_b = d6_b * x13 + d6_b * x13
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d5_v * g_y32(g_i_) + d7_b * g_x32(g_i_) + d8_b 
     ** g_x13(g_i_) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + d5_v * y32
        d3_v = d2_w / x2ap3
        d2_b = 1.0d0 / x2ap3
        d3_b = -d3_v / x2ap3
        do g_i_ = 1, g_p_
          g_sq(g_i_, 2, 2) = d3_b * g_x2ap3(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        sq(2, 2) = d3_v
C--------
        d4_v = -x21 * x13
        d5_v = d4_v * x13
        d4_b = y21 * x13
        d5_b = y21 * d4_v + d4_b * (-x21)
        d7_b = -(d4_b * x13)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d5_v * g_y21(g_i_) + d5_b * g_x13(g_i_) + d7_b 
     ** g_x21(g_i_)
        enddo
        d1_w = d5_v * y21
        d3_v = x21 * x21
        d5_v = d3_v * x13
        d6_b = y13 * x13
        d7_b = y13 * d3_v
        d8_b = d6_b * x21 + d6_b * x21
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d5_v * g_y13(g_i_) + d7_b * g_x13(g_i_) + d8_b 
     ** g_x21(g_i_) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + d5_v * y13
        d3_v = d2_w / x2ap3
        d2_b = 1.0d0 / x2ap3
        d3_b = -d3_v / x2ap3
        do g_i_ = 1, g_p_
          g_sq(g_i_, 2, 3) = d3_b * g_x2ap3(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        sq(2, 3) = d3_v
C--------
C
        d2_v = x21 * x21
        d4_v = d2_v * y32
        d4_b = y32 * y32
        d3_b = d4_v + y32 * d2_v
        d5_b = d4_b * x21 + d4_b * x21
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_y32(g_i_) + d5_b * g_x21(g_i_)
        enddo
        d1_w = d4_v * y32
        d3_v = x32 * x32
        d5_v = d3_v * y21
        d6_b = -y21 * y21
        d5_b = -d5_v + (-y21) * d3_v
        d7_b = d6_b * x32 + d6_b * x32
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d5_b * g_y21(g_i_) + d7_b * g_x32(g_i_) + g_d1_
     *w(g_i_)
        enddo
        d2_w = d1_w - d5_v * y21
        d3_v = d2_w / x2ap3
        d2_b = 1.0d0 / x2ap3
        d3_b = -d3_v / x2ap3
        do g_i_ = 1, g_p_
          g_sq(g_i_, 3, 1) = d3_b * g_x2ap3(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        sq(3, 1) = d3_v
C--------
        d3_v = -x13 * x13
        d5_v = d3_v * y32
        d4_b = y32 * y32
        d3_b = d5_v + y32 * d3_v
        d6_b = d4_b * (-x13) + (-(d4_b * x13))
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_y32(g_i_) + d6_b * g_x13(g_i_)
        enddo
        d1_w = d5_v * y32
        d3_v = x32 * x32
        d5_v = d3_v * y13
        d6_b = y13 * y13
        d5_b = d5_v + y13 * d3_v
        d7_b = d6_b * x32 + d6_b * x32
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d5_b * g_y13(g_i_) + d7_b * g_x32(g_i_) + g_d1_
     *w(g_i_)
        enddo
        d2_w = d1_w + d5_v * y13
        d3_v = d2_w / x2ap3
        d2_b = 1.0d0 / x2ap3
        d3_b = -d3_v / x2ap3
        do g_i_ = 1, g_p_
          g_sq(g_i_, 3, 2) = d3_b * g_x2ap3(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        sq(3, 2) = d3_v
C--------
        d3_v = -x21 * x21
        d5_v = d3_v * y13
        d4_b = y13 * y13
        d3_b = d5_v + y13 * d3_v
        d6_b = d4_b * (-x21) + (-(d4_b * x21))
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_y13(g_i_) + d6_b * g_x21(g_i_)
        enddo
        d1_w = d5_v * y13
        d3_v = x13 * x13
        d5_v = d3_v * y21
        d6_b = y21 * y21
        d5_b = d5_v + y21 * d3_v
        d7_b = d6_b * x13 + d6_b * x13
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d5_b * g_y21(g_i_) + d7_b * g_x13(g_i_) + g_d1_
     *w(g_i_)
        enddo
        d2_w = d1_w + d5_v * y21
        d3_v = d2_w / x2ap3
        d2_b = 1.0d0 / x2ap3
        d3_b = -d3_v / x2ap3
        do g_i_ = 1, g_p_
          g_sq(g_i_, 3, 3) = d3_b * g_x2ap3(g_i_) + d2_b * g_d2_w(g_i_)
        enddo
        sq(3, 3) = d3_v
C--------
C
C.....MULTIPLY MATRIX [SQ] BY THE SQUARE ROOT OF THE STIFFNESS FACTOR
C.....AND THE SQUARE ROOT OF THE AREA DIVIDED BY THE WEIGHT OF THE
C.....NUMERICAL INTEGRATION (EQUAL TO THREE W/ MID-POINT RULE HERE)
C
        d3_b = 1.0d0 / 3.00d+00 * fb
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_area(g_i_)
        enddo
        d1_w = fb * area / 3.00d+00
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_comphBM.f',
     +1217)
        endif
        do g_i_ = 1, g_p_
          g_factor(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        factor = d2_v
C--------
C
        do 99964 j = 1, 3
          do 99965 i = 1, 3
            do g_i_ = 1, g_p_
              g_sq(g_i_, i, j) = factor * g_sq(g_i_, i, j) + sq(i, j) * 
     *g_factor(g_i_)
            enddo
            sq(i, j) = factor * sq(i, j)
C--------
2002        continue
99965     continue
2001      continue
99964   continue
C
C.....ESTIMATE THE MATRIX [RM] AT THE FIRST GAUSS INTEGRATION POINT
C
        d5_b = -(1.0d0 / 3.00d+00) + gauss21
        do g_i_ = 1, g_p_
          g_rm1(g_i_, 1, 1) = d5_b * g_al1(g_i_)
        enddo
        rm1(1, 1) = gauss31 + al1 * gauss21 - (1.00d+00 + al1) / 3.00d+0
     *0
C--------
        d5_b = -(1.0d0 / 3.00d+00) + gauss21
        do g_i_ = 1, g_p_
          g_rm1(g_i_, 1, 2) = d5_b * g_bl1(g_i_)
        enddo
        rm1(1, 2) = gauss11 + bl1 * gauss21 - (1.00d+00 + bl1) / 3.00d+0
     *0
C--------
C
        d5_b = -(1.0d0 / 3.00d+00) + gauss31
        do g_i_ = 1, g_p_
          g_rm1(g_i_, 2, 3) = d5_b * g_al2(g_i_)
        enddo
        rm1(2, 3) = gauss11 + al2 * gauss31 - (1.00d+00 + al2) / 3.00d+0
     *0
C--------
        d5_b = -(1.0d0 / 3.00d+00) + gauss31
        do g_i_ = 1, g_p_
          g_rm1(g_i_, 2, 4) = d5_b * g_bl2(g_i_)
        enddo
        rm1(2, 4) = gauss21 + bl2 * gauss31 - (1.00d+00 + bl2) / 3.00d+0
     *0
C--------
C
        d5_b = -(1.0d0 / 3.00d+00) + gauss11
        do g_i_ = 1, g_p_
          g_rm1(g_i_, 3, 5) = d5_b * g_al3(g_i_)
        enddo
        rm1(3, 5) = gauss21 + al3 * gauss11 - (1.00d+00 + al3) / 3.00d+0
     *0
C--------
        d5_b = -(1.0d0 / 3.00d+00) + gauss11
        do g_i_ = 1, g_p_
          g_rm1(g_i_, 3, 6) = d5_b * g_bl3(g_i_)
        enddo
        rm1(3, 6) = gauss31 + bl3 * gauss11 - (1.00d+00 + bl3) / 3.00d+0
     *0
C--------
C
C.....ESTIMATE THE MATRIX [RM] AT THE SECOND GAUSS INTEGRATION POINT
C
        d5_b = -(1.0d0 / 3.00d+00) + gauss22
        do g_i_ = 1, g_p_
          g_rm2(g_i_, 1, 1) = d5_b * g_al1(g_i_)
        enddo
        rm2(1, 1) = gauss32 + al1 * gauss22 - (1.00d+00 + al1) / 3.00d+0
     *0
C--------
        d5_b = -(1.0d0 / 3.00d+00) + gauss22
        do g_i_ = 1, g_p_
          g_rm2(g_i_, 1, 2) = d5_b * g_bl1(g_i_)
        enddo
        rm2(1, 2) = gauss12 + bl1 * gauss22 - (1.00d+00 + bl1) / 3.00d+0
     *0
C--------
C
        d5_b = -(1.0d0 / 3.00d+00) + gauss32
        do g_i_ = 1, g_p_
          g_rm2(g_i_, 2, 3) = d5_b * g_al2(g_i_)
        enddo
        rm2(2, 3) = gauss12 + al2 * gauss32 - (1.00d+00 + al2) / 3.00d+0
     *0
C--------
        d5_b = -(1.0d0 / 3.00d+00) + gauss32
        do g_i_ = 1, g_p_
          g_rm2(g_i_, 2, 4) = d5_b * g_bl2(g_i_)
        enddo
        rm2(2, 4) = gauss22 + bl2 * gauss32 - (1.00d+00 + bl2) / 3.00d+0
     *0
C--------
C
        d5_b = -(1.0d0 / 3.00d+00) + gauss12
        do g_i_ = 1, g_p_
          g_rm2(g_i_, 3, 5) = d5_b * g_al3(g_i_)
        enddo
        rm2(3, 5) = gauss22 + al3 * gauss12 - (1.00d+00 + al3) / 3.00d+0
     *0
C--------
        d5_b = -(1.0d0 / 3.00d+00) + gauss12
        do g_i_ = 1, g_p_
          g_rm2(g_i_, 3, 6) = d5_b * g_bl3(g_i_)
        enddo
        rm2(3, 6) = gauss32 + bl3 * gauss12 - (1.00d+00 + bl3) / 3.00d+0
     *0
C--------
C
C.....ESTIMATE THE MATRIX [RM] AT THE THIRD GAUSS INTEGRATION POINT
C
        d5_b = -(1.0d0 / 3.00d+00) + gauss23
        do g_i_ = 1, g_p_
          g_rm3(g_i_, 1, 1) = d5_b * g_al1(g_i_)
        enddo
        rm3(1, 1) = gauss33 + al1 * gauss23 - (1.00d+00 + al1) / 3.00d+0
     *0
C--------
        d5_b = -(1.0d0 / 3.00d+00) + gauss23
        do g_i_ = 1, g_p_
          g_rm3(g_i_, 1, 2) = d5_b * g_bl1(g_i_)
        enddo
        rm3(1, 2) = gauss13 + bl1 * gauss23 - (1.00d+00 + bl1) / 3.00d+0
     *0
C--------
C
        d5_b = -(1.0d0 / 3.00d+00) + gauss33
        do g_i_ = 1, g_p_
          g_rm3(g_i_, 2, 3) = d5_b * g_al2(g_i_)
        enddo
        rm3(2, 3) = gauss13 + al2 * gauss33 - (1.00d+00 + al2) / 3.00d+0
     *0
C--------
        d5_b = -(1.0d0 / 3.00d+00) + gauss33
        do g_i_ = 1, g_p_
          g_rm3(g_i_, 2, 4) = d5_b * g_bl2(g_i_)
        enddo
        rm3(2, 4) = gauss23 + bl2 * gauss33 - (1.00d+00 + bl2) / 3.00d+0
     *0
C--------
C
        d5_b = -(1.0d0 / 3.00d+00) + gauss13
        do g_i_ = 1, g_p_
          g_rm3(g_i_, 3, 5) = d5_b * g_al3(g_i_)
        enddo
        rm3(3, 5) = gauss23 + al3 * gauss13 - (1.00d+00 + al3) / 3.00d+0
     *0
C--------
        d5_b = -(1.0d0 / 3.00d+00) + gauss13
        do g_i_ = 1, g_p_
          g_rm3(g_i_, 3, 6) = d5_b * g_bl3(g_i_)
        enddo
        rm3(3, 6) = gauss33 + bl3 * gauss13 - (1.00d+00 + bl3) / 3.00d+0
     *0
C--------
C
C.....COMPUTE THE MATRIX-MATRIX PRODUCT [rm1]^T*[sq]^T INTO MATRIX [rm1tsqt]
C.....COMPUTE THE MATRIX-MATRIX PRODUCT [rm2]^T*[sq]^T INTO MATRIX [rm2tsqt]
C.....COMPUTE THE MATRIX-MATRIX PRODUCT [rm3]^T*[sq]^T INTO MATRIX [rm3tsqt]
C
        do 99960 j = 1, 3
C
          do g_i_ = 1, g_p_
            g_s1(g_i_) = g_sq(g_i_, j, 1)
          enddo
          s1 = sq(j, 1)
C--------
          do g_i_ = 1, g_p_
            g_s2(g_i_) = g_sq(g_i_, j, 2)
          enddo
          s2 = sq(j, 2)
C--------
          do g_i_ = 1, g_p_
            g_s3(g_i_) = g_sq(g_i_, j, 3)
          enddo
          s3 = sq(j, 3)
C--------
C
          do 99963 i = 1, 6
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = s2 * g_rm1(g_i_, 2, i) + rm1(2, i) * g_s2(g
     *_i_) + s1 * g_rm1(g_i_, 1, i) + rm1(1, i) * g_s1(g_i_)
            enddo
            d1_w = s1 * rm1(1, i) + s2 * rm1(2, i)
            do g_i_ = 1, g_p_
              g_rm1tsqt(g_i_, i, j) = s3 * g_rm1(g_i_, 3, i) + rm1(3, i)
     * * g_s3(g_i_) + g_d1_w(g_i_)
            enddo
            rm1tsqt(i, j) = d1_w + s3 * rm1(3, i)
C--------
3002        continue
99963     continue
C
          do 99962 i = 1, 6
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = s2 * g_rm2(g_i_, 2, i) + rm2(2, i) * g_s2(g
     *_i_) + s1 * g_rm2(g_i_, 1, i) + rm2(1, i) * g_s1(g_i_)
            enddo
            d1_w = s1 * rm2(1, i) + s2 * rm2(2, i)
            do g_i_ = 1, g_p_
              g_rm2tsqt(g_i_, i, j) = s3 * g_rm2(g_i_, 3, i) + rm2(3, i)
     * * g_s3(g_i_) + g_d1_w(g_i_)
            enddo
            rm2tsqt(i, j) = d1_w + s3 * rm2(3, i)
C--------
3003        continue
99962     continue
C
          do 99961 i = 1, 6
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = s2 * g_rm3(g_i_, 2, i) + rm3(2, i) * g_s2(g
     *_i_) + s1 * g_rm3(g_i_, 1, i) + rm3(1, i) * g_s1(g_i_)
            enddo
            d1_w = s1 * rm3(1, i) + s2 * rm3(2, i)
            do g_i_ = 1, g_p_
              g_rm3tsqt(g_i_, i, j) = s3 * g_rm3(g_i_, 3, i) + rm3(3, i)
     * * g_s3(g_i_) + g_d1_w(g_i_)
            enddo
            rm3tsqt(i, j) = d1_w + s3 * rm3(3, i)
C--------
3004        continue
99961     continue
C
3001      continue
99960   continue
C
C.....ASSEMBLE MATRICES [Lh1], [Lh2] AND [Lh3] FOR THE THREE GAUSS POINTS:
C.....[Lh1] = [q]^T * [rm1]^T * [sq]^T = [q]^T * [rm1tsqt]
C.....[Lh2] = [q]^T * [rm2]^T * [sq]^T = [q]^T * [rm2tsqt]
C.....[Lh3] = [q]^T * [rm3]^T * [sq]^T = [q]^T * [rm3tsqt]
C
        do 99956 i = 1, 9
C
          do g_i_ = 1, g_p_
            g_q1(g_i_) = g_q(g_i_, 1, i)
          enddo
          q1 = q(1, i)
C--------
          do g_i_ = 1, g_p_
            g_q2(g_i_) = g_q(g_i_, 2, i)
          enddo
          q2 = q(2, i)
C--------
          do g_i_ = 1, g_p_
            g_q3(g_i_) = g_q(g_i_, 3, i)
          enddo
          q3 = q(3, i)
C--------
          do g_i_ = 1, g_p_
            g_q4(g_i_) = g_q(g_i_, 4, i)
          enddo
          q4 = q(4, i)
C--------
          do g_i_ = 1, g_p_
            g_q5(g_i_) = g_q(g_i_, 5, i)
          enddo
          q5 = q(5, i)
C--------
          do g_i_ = 1, g_p_
            g_q6(g_i_) = g_q(g_i_, 6, i)
          enddo
          q6 = q(6, i)
C--------
C
          do 99959 j = 1, 3
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = q2 * g_rm1tsqt(g_i_, 2, j) + rm1tsqt(2, j) 
     ** g_q2(g_i_) + q1 * g_rm1tsqt(g_i_, 1, j) + rm1tsqt(1, j) * g_q1(g
     *_i_)
            enddo
            d1_w = q1 * rm1tsqt(1, j) + q2 * rm1tsqt(2, j)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = q4 * g_rm1tsqt(g_i_, 4, j) + rm1tsqt(4, j) 
     ** g_q4(g_i_) + q3 * g_rm1tsqt(g_i_, 3, j) + rm1tsqt(3, j) * g_q3(g
     *_i_) + g_d1_w(g_i_)
            enddo
            d2_w = d1_w + q3 * rm1tsqt(3, j) + q4 * rm1tsqt(4, j)
            do g_i_ = 1, g_p_
              g_lh1(g_i_, i, j) = q6 * g_rm1tsqt(g_i_, 6, j) + rm1tsqt(6
     *, j) * g_q6(g_i_) + q5 * g_rm1tsqt(g_i_, 5, j) + rm1tsqt(5, j) * g
     *_q5(g_i_) + g_d2_w(g_i_)
            enddo
            lh1(i, j) = d2_w + q5 * rm1tsqt(5, j) + q6 * rm1tsqt(6, j)
C--------
4002        continue
99959     continue
C
          do 99958 j = 1, 3
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = q2 * g_rm2tsqt(g_i_, 2, j) + rm2tsqt(2, j) 
     ** g_q2(g_i_) + q1 * g_rm2tsqt(g_i_, 1, j) + rm2tsqt(1, j) * g_q1(g
     *_i_)
            enddo
            d1_w = q1 * rm2tsqt(1, j) + q2 * rm2tsqt(2, j)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = q4 * g_rm2tsqt(g_i_, 4, j) + rm2tsqt(4, j) 
     ** g_q4(g_i_) + q3 * g_rm2tsqt(g_i_, 3, j) + rm2tsqt(3, j) * g_q3(g
     *_i_) + g_d1_w(g_i_)
            enddo
            d2_w = d1_w + q3 * rm2tsqt(3, j) + q4 * rm2tsqt(4, j)
            do g_i_ = 1, g_p_
              g_lh2(g_i_, i, j) = q6 * g_rm2tsqt(g_i_, 6, j) + rm2tsqt(6
     *, j) * g_q6(g_i_) + q5 * g_rm2tsqt(g_i_, 5, j) + rm2tsqt(5, j) * g
     *_q5(g_i_) + g_d2_w(g_i_)
            enddo
            lh2(i, j) = d2_w + q5 * rm2tsqt(5, j) + q6 * rm2tsqt(6, j)
C--------
4003        continue
99958     continue
C
          do 99957 j = 1, 3
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = q2 * g_rm3tsqt(g_i_, 2, j) + rm3tsqt(2, j) 
     ** g_q2(g_i_) + q1 * g_rm3tsqt(g_i_, 1, j) + rm3tsqt(1, j) * g_q1(g
     *_i_)
            enddo
            d1_w = q1 * rm3tsqt(1, j) + q2 * rm3tsqt(2, j)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = q4 * g_rm3tsqt(g_i_, 4, j) + rm3tsqt(4, j) 
     ** g_q4(g_i_) + q3 * g_rm3tsqt(g_i_, 3, j) + rm3tsqt(3, j) * g_q3(g
     *_i_) + g_d1_w(g_i_)
            enddo
            d2_w = d1_w + q3 * rm3tsqt(3, j) + q4 * rm3tsqt(4, j)
            do g_i_ = 1, g_p_
              g_lh3(g_i_, i, j) = q6 * g_rm3tsqt(g_i_, 6, j) + rm3tsqt(6
     *, j) * g_q6(g_i_) + q5 * g_rm3tsqt(g_i_, 5, j) + rm3tsqt(5, j) * g
     *_q5(g_i_) + g_d2_w(g_i_)
            enddo
            lh3(i, j) = d2_w + q5 * rm3tsqt(5, j) + q6 * rm3tsqt(6, j)
C--------
4004        continue
99957     continue
C
4001      continue
99956   continue
C
C     --------------------------------------------------------
C
C     FORM THE LOCAL NORMAL FORCE-STRAIN MATRICES FOR MEMBRANE
C
C              AND FOR THE THREE INTEGRATION POINTS
C
C     --------------------------------------------------------
C
C.....GET THE DISTANCES BETWEEN NODES 1, 2 AND 3 AND THE CENTROID
C
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
C
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
C
C.....CALCULATE BASIC COEFFICIENTS FOR THE ASSEMBLY OF MATRIX [HMV]
C
        d5_b = 2.25d+00 * y30 + 2.25d+00 * y30
        d6_b = 2.25d+00 * x30 + 2.25d+00 * x30
        do g_i_ = 1, g_p_
          g_aa12(g_i_) = d5_b * g_y30(g_i_) + d6_b * g_x30(g_i_)
        enddo
        aa12 = 2.25d+00 * (x30 * x30 + y30 * y30)
C--------
        d5_b = 2.25d+00 * y10 + 2.25d+00 * y10
        d6_b = 2.25d+00 * x10 + 2.25d+00 * x10
        do g_i_ = 1, g_p_
          g_aa23(g_i_) = d5_b * g_y10(g_i_) + d6_b * g_x10(g_i_)
        enddo
        aa23 = 2.25d+00 * (x10 * x10 + y10 * y10)
C--------
        d5_b = 2.25d+00 * y20 + 2.25d+00 * y20
        d6_b = 2.25d+00 * x20 + 2.25d+00 * x20
        do g_i_ = 1, g_p_
          g_aa31(g_i_) = d5_b * g_y20(g_i_) + d6_b * g_x20(g_i_)
        enddo
        aa31 = 2.25d+00 * (x20 * x20 + y20 * y20)
C--------
C
        d2_v = 128.00d+00 * aa12
        d3_v = 15.00d+00 / d2_v
        d3_b = -d3_v / d2_v * 128.00d+00
        do g_i_ = 1, g_p_
          g_caa12(g_i_) = d3_b * g_aa12(g_i_)
        enddo
        caa12 = d3_v
C--------
        d2_v = 128.00d+00 * aa23
        d3_v = 15.00d+00 / d2_v
        d3_b = -d3_v / d2_v * 128.00d+00
        do g_i_ = 1, g_p_
          g_caa23(g_i_) = d3_b * g_aa23(g_i_)
        enddo
        caa23 = d3_v
C--------
        d2_v = 128.00d+00 * aa31
        d3_v = 15.00d+00 / d2_v
        d3_b = -d3_v / d2_v * 128.00d+00
        do g_i_ = 1, g_p_
          g_caa31(g_i_) = d3_b * g_aa31(g_i_)
        enddo
        caa31 = d3_v
C--------
C
        d4_b = y12 + y12
        d5_b = x12 + x12
        do g_i_ = 1, g_p_
          g_ss12(g_i_) = d4_b * g_y12(g_i_) + d5_b * g_x12(g_i_)
        enddo
        ss12 = x12 * x12 + y12 * y12
C--------
        d4_b = y23 + y23
        d5_b = x23 + x23
        do g_i_ = 1, g_p_
          g_ss23(g_i_) = d4_b * g_y23(g_i_) + d5_b * g_x23(g_i_)
        enddo
        ss23 = x23 * x23 + y23 * y23
C--------
        d4_b = y31 + y31
        d5_b = x31 + x31
        do g_i_ = 1, g_p_
          g_ss31(g_i_) = d4_b * g_y31(g_i_) + d5_b * g_x31(g_i_)
        enddo
        ss31 = x31 * x31 + y31 * y31
C--------
C
        d2_v = 16.00d+00 * area
        d3_v = 3.00d+00 / d2_v
        d3_b = -d3_v / d2_v * 16.00d+00
        do g_i_ = 1, g_p_
          g_ca(g_i_) = d3_b * g_area(g_i_)
        enddo
        ca = d3_v
C--------
C
        do g_i_ = 1, g_p_
          g_cax10(g_i_) = ca * g_x10(g_i_) + x10 * g_ca(g_i_)
        enddo
        cax10 = ca * x10
C--------
        do g_i_ = 1, g_p_
          g_cax20(g_i_) = ca * g_x20(g_i_) + x20 * g_ca(g_i_)
        enddo
        cax20 = ca * x20
C--------
        do g_i_ = 1, g_p_
          g_cax30(g_i_) = ca * g_x30(g_i_) + x30 * g_ca(g_i_)
        enddo
        cax30 = ca * x30
C--------
C
        do g_i_ = 1, g_p_
          g_cay10(g_i_) = ca * g_y10(g_i_) + y10 * g_ca(g_i_)
        enddo
        cay10 = ca * y10
C--------
        do g_i_ = 1, g_p_
          g_cay20(g_i_) = ca * g_y20(g_i_) + y20 * g_ca(g_i_)
        enddo
        cay20 = ca * y20
C--------
        do g_i_ = 1, g_p_
          g_cay30(g_i_) = ca * g_y30(g_i_) + y30 * g_ca(g_i_)
        enddo
        cay30 = ca * y30
C--------
C
C.....CONSTRUCT LOCAL MATRIX [HMV] W/ SHAPE FUNCTION DERIVATIVES
C
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 1, 1) = cay30 * g_x32(g_i_) + x32 * g_cay30(g_i_)
        enddo
        hmv(1, 1) = cay30 * x32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 1, 2) = cay30 * g_y32(g_i_) + y32 * g_cay30(g_i_)
        enddo
        hmv(1, 2) = cay30 * y32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 1, 3) = cay30 * g_x13(g_i_) + x13 * g_cay30(g_i_)
        enddo
        hmv(1, 3) = cay30 * x13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 1, 4) = cay30 * g_y13(g_i_) + y13 * g_cay30(g_i_)
        enddo
        hmv(1, 4) = cay30 * y13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 1, 5) = cay30 * g_x21(g_i_) + x21 * g_cay30(g_i_)
        enddo
        hmv(1, 5) = cay30 * x21
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 1, 6) = cay30 * g_y21(g_i_) + y21 * g_cay30(g_i_)
        enddo
        hmv(1, 6) = cay30 * y21
C--------
        d6_v = ss23 - ss31 + 2.40d+00 * aa12
        d8_v = d6_v * y30
        d4_b = caa12 * y30
        d5_b = caa12 * d6_v
        d8_b = d4_b * 2.40d+00
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d8_v * g_caa12(g_i_) + d5_b * g_y30(g_i_) + d8_
     *b * g_aa12(g_i_) + (-d4_b) * g_ss31(g_i_) + d4_b * g_ss23(g_i_)
        enddo
        d1_w = d8_v * caa12
        d3_v = 4.00d+00 * area
        d5_v = d3_v * x30
        d7_b = caa12 * d3_v
        d8_b = caa12 * x30 * 4.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 1, 7) = d5_v * g_caa12(g_i_) + d7_b * g_x30(g_i_) 
     *+ d8_b * g_area(g_i_) + g_d1_w(g_i_)
        enddo
        hmv(1, 7) = d1_w + d5_v * caa12
C--------
        d4_b = 9.00d+00 / 16.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 1, 8) = -g_hmv(g_i_, 1, 7) + d4_b * g_y30(g_i_)
        enddo
        hmv(1, 8) = 9.00d+00 / 16.00d+00 * y30 - hmv(1, 7)
C--------
        d2_b = 3.00d+00 / 16.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 1, 9) = d2_b * g_y30(g_i_)
        enddo
        hmv(1, 9) = 3.00d+00 / 16.00d+00 * y30
C--------
C
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 2, 1) = cay10 * g_x32(g_i_) + x32 * g_cay10(g_i_)
        enddo
        hmv(2, 1) = cay10 * x32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 2, 2) = cay10 * g_y32(g_i_) + y32 * g_cay10(g_i_)
        enddo
        hmv(2, 2) = cay10 * y32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 2, 3) = cay10 * g_x13(g_i_) + x13 * g_cay10(g_i_)
        enddo
        hmv(2, 3) = cay10 * x13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 2, 4) = cay10 * g_y13(g_i_) + y13 * g_cay10(g_i_)
        enddo
        hmv(2, 4) = cay10 * y13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 2, 5) = cay10 * g_x21(g_i_) + x21 * g_cay10(g_i_)
        enddo
        hmv(2, 5) = cay10 * x21
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 2, 6) = cay10 * g_y21(g_i_) + y21 * g_cay10(g_i_)
        enddo
        hmv(2, 6) = cay10 * y21
C--------
        d2_b = 3.00d+00 / 16.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 2, 7) = d2_b * g_y10(g_i_)
        enddo
        hmv(2, 7) = 3.00d+00 / 16.00d+00 * y10
C--------
        d6_v = ss31 - ss12 + 2.40d+00 * aa23
        d8_v = d6_v * y10
        d4_b = caa23 * y10
        d5_b = caa23 * d6_v
        d8_b = d4_b * 2.40d+00
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d8_v * g_caa23(g_i_) + d5_b * g_y10(g_i_) + d8_
     *b * g_aa23(g_i_) + (-d4_b) * g_ss12(g_i_) + d4_b * g_ss31(g_i_)
        enddo
        d1_w = d8_v * caa23
        d3_v = 4.00d+00 * area
        d5_v = d3_v * x10
        d7_b = caa23 * d3_v
        d8_b = caa23 * x10 * 4.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 2, 8) = d5_v * g_caa23(g_i_) + d7_b * g_x10(g_i_) 
     *+ d8_b * g_area(g_i_) + g_d1_w(g_i_)
        enddo
        hmv(2, 8) = d1_w + d5_v * caa23
C--------
        d4_b = 9.00d+00 / 16.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 2, 9) = -g_hmv(g_i_, 2, 8) + d4_b * g_y10(g_i_)
        enddo
        hmv(2, 9) = 9.00d+00 / 16.00d+00 * y10 - hmv(2, 8)
C--------
C
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 3, 1) = cay20 * g_x32(g_i_) + x32 * g_cay20(g_i_)
        enddo
        hmv(3, 1) = cay20 * x32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 3, 2) = cay20 * g_y32(g_i_) + y32 * g_cay20(g_i_)
        enddo
        hmv(3, 2) = cay20 * y32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 3, 3) = cay20 * g_x13(g_i_) + x13 * g_cay20(g_i_)
        enddo
        hmv(3, 3) = cay20 * x13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 3, 4) = cay20 * g_y13(g_i_) + y13 * g_cay20(g_i_)
        enddo
        hmv(3, 4) = cay20 * y13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 3, 5) = cay20 * g_x21(g_i_) + x21 * g_cay20(g_i_)
        enddo
        hmv(3, 5) = cay20 * x21
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 3, 6) = cay20 * g_y21(g_i_) + y21 * g_cay20(g_i_)
        enddo
        hmv(3, 6) = cay20 * y21
C--------
        d6_v = ss23 - ss12 + 2.40d+00 * aa31
        d8_v = d6_v * y20
        d4_b = caa31 * y20
        d5_b = caa31 * d6_v
        d8_b = d4_b * 2.40d+00
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d8_v * g_caa31(g_i_) + d5_b * g_y20(g_i_) + d8_
     *b * g_aa31(g_i_) + (-d4_b) * g_ss12(g_i_) + d4_b * g_ss23(g_i_)
        enddo
        d1_w = d8_v * caa31
        d3_v = 4.00d+00 * area
        d5_v = d3_v * x20
        d7_b = -caa31 * d3_v
        d8_b = -caa31 * x20 * 4.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 3, 7) = -d5_v * g_caa31(g_i_) + d7_b * g_x20(g_i_)
     * + d8_b * g_area(g_i_) + g_d1_w(g_i_)
        enddo
        hmv(3, 7) = d1_w - d5_v * caa31
C--------
        d2_b = 3.00d+00 / 16.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 3, 8) = d2_b * g_y20(g_i_)
        enddo
        hmv(3, 8) = 3.00d+00 / 16.00d+00 * y20
C--------
        d4_b = 9.00d+00 / 16.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 3, 9) = -g_hmv(g_i_, 3, 7) + d4_b * g_y20(g_i_)
        enddo
        hmv(3, 9) = 9.00d+00 / 16.00d+00 * y20 - hmv(3, 7)
C--------
C
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 4, 1) = -cax30 * g_x32(g_i_) + (-x32) * g_cax30(g_
     *i_)
        enddo
        hmv(4, 1) = -cax30 * x32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 4, 2) = -cax30 * g_y32(g_i_) + (-y32) * g_cax30(g_
     *i_)
        enddo
        hmv(4, 2) = -cax30 * y32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 4, 3) = -cax30 * g_x13(g_i_) + (-x13) * g_cax30(g_
     *i_)
        enddo
        hmv(4, 3) = -cax30 * x13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 4, 4) = -cax30 * g_y13(g_i_) + (-y13) * g_cax30(g_
     *i_)
        enddo
        hmv(4, 4) = -cax30 * y13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 4, 5) = -cax30 * g_x21(g_i_) + (-x21) * g_cax30(g_
     *i_)
        enddo
        hmv(4, 5) = -cax30 * x21
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 4, 6) = -cax30 * g_y21(g_i_) + (-y21) * g_cax30(g_
     *i_)
        enddo
        hmv(4, 6) = -cax30 * y21
C--------
        d6_v = ss31 - ss23 - 2.40d+00 * aa12
        d8_v = d6_v * x30
        d4_b = caa12 * x30
        d5_b = caa12 * d6_v
        d8_b = -d4_b * 2.40d+00
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d8_v * g_caa12(g_i_) + d5_b * g_x30(g_i_) + d8_
     *b * g_aa12(g_i_) + (-d4_b) * g_ss23(g_i_) + d4_b * g_ss31(g_i_)
        enddo
        d1_w = d8_v * caa12
        d3_v = 4.00d+00 * area
        d5_v = d3_v * y30
        d7_b = caa12 * d3_v
        d8_b = caa12 * y30 * 4.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 4, 7) = d5_v * g_caa12(g_i_) + d7_b * g_y30(g_i_) 
     *+ d8_b * g_area(g_i_) + g_d1_w(g_i_)
        enddo
        hmv(4, 7) = d1_w + d5_v * caa12
C--------
        d4_b = -(9.00d+00 / 16.00d+00)
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 4, 8) = -g_hmv(g_i_, 4, 7) + d4_b * g_x30(g_i_)
        enddo
        hmv(4, 8) = -(9.00d+00 / 16.00d+00) * x30 - hmv(4, 7)
C--------
        d2_b = -(3.00d+00 / 16.00d+00)
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 4, 9) = d2_b * g_x30(g_i_)
        enddo
        hmv(4, 9) = -(3.00d+00 / 16.00d+00) * x30
C--------
C
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 5, 1) = -cax10 * g_x32(g_i_) + (-x32) * g_cax10(g_
     *i_)
        enddo
        hmv(5, 1) = -cax10 * x32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 5, 2) = -cax10 * g_y32(g_i_) + (-y32) * g_cax10(g_
     *i_)
        enddo
        hmv(5, 2) = -cax10 * y32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 5, 3) = -cax10 * g_x13(g_i_) + (-x13) * g_cax10(g_
     *i_)
        enddo
        hmv(5, 3) = -cax10 * x13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 5, 4) = -cax10 * g_y13(g_i_) + (-y13) * g_cax10(g_
     *i_)
        enddo
        hmv(5, 4) = -cax10 * y13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 5, 5) = -cax10 * g_x21(g_i_) + (-x21) * g_cax10(g_
     *i_)
        enddo
        hmv(5, 5) = -cax10 * x21
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 5, 6) = -cax10 * g_y21(g_i_) + (-y21) * g_cax10(g_
     *i_)
        enddo
        hmv(5, 6) = -cax10 * y21
C--------
        d2_b = -(3.00d+00 / 16.00d+00)
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 5, 7) = d2_b * g_x10(g_i_)
        enddo
        hmv(5, 7) = -(3.00d+00 / 16.00d+00) * x10
C--------
        d6_v = ss12 - ss31 - 2.40d+00 * aa23
        d8_v = d6_v * x10
        d4_b = caa23 * x10
        d5_b = caa23 * d6_v
        d8_b = -d4_b * 2.40d+00
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d8_v * g_caa23(g_i_) + d5_b * g_x10(g_i_) + d8_
     *b * g_aa23(g_i_) + (-d4_b) * g_ss31(g_i_) + d4_b * g_ss12(g_i_)
        enddo
        d1_w = d8_v * caa23
        d3_v = 4.00d+00 * area
        d5_v = d3_v * y10
        d7_b = caa23 * d3_v
        d8_b = caa23 * y10 * 4.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 5, 8) = d5_v * g_caa23(g_i_) + d7_b * g_y10(g_i_) 
     *+ d8_b * g_area(g_i_) + g_d1_w(g_i_)
        enddo
        hmv(5, 8) = d1_w + d5_v * caa23
C--------
        d4_b = -(9.00d+00 / 16.00d+00)
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 5, 9) = -g_hmv(g_i_, 5, 8) + d4_b * g_x10(g_i_)
        enddo
        hmv(5, 9) = -(9.00d+00 / 16.00d+00) * x10 - hmv(5, 8)
C--------
C
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 6, 1) = -cax20 * g_x32(g_i_) + (-x32) * g_cax20(g_
     *i_)
        enddo
        hmv(6, 1) = -cax20 * x32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 6, 2) = -cax20 * g_y32(g_i_) + (-y32) * g_cax20(g_
     *i_)
        enddo
        hmv(6, 2) = -cax20 * y32
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 6, 3) = -cax20 * g_x13(g_i_) + (-x13) * g_cax20(g_
     *i_)
        enddo
        hmv(6, 3) = -cax20 * x13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 6, 4) = -cax20 * g_y13(g_i_) + (-y13) * g_cax20(g_
     *i_)
        enddo
        hmv(6, 4) = -cax20 * y13
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 6, 5) = -cax20 * g_x21(g_i_) + (-x21) * g_cax20(g_
     *i_)
        enddo
        hmv(6, 5) = -cax20 * x21
C--------
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 6, 6) = -cax20 * g_y21(g_i_) + (-y21) * g_cax20(g_
     *i_)
        enddo
        hmv(6, 6) = -cax20 * y21
C--------
        d6_v = ss12 - ss23 - 2.40d+00 * aa31
        d8_v = d6_v * x20
        d4_b = caa31 * x20
        d5_b = caa31 * d6_v
        d8_b = -d4_b * 2.40d+00
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d8_v * g_caa31(g_i_) + d5_b * g_x20(g_i_) + d8_
     *b * g_aa31(g_i_) + (-d4_b) * g_ss23(g_i_) + d4_b * g_ss12(g_i_)
        enddo
        d1_w = d8_v * caa31
        d3_v = 4.00d+00 * area
        d5_v = d3_v * y20
        d7_b = -caa31 * d3_v
        d8_b = -caa31 * y20 * 4.00d+00
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 6, 7) = -d5_v * g_caa31(g_i_) + d7_b * g_y20(g_i_)
     * + d8_b * g_area(g_i_) + g_d1_w(g_i_)
        enddo
        hmv(6, 7) = d1_w - d5_v * caa31
C--------
        d2_b = -(3.00d+00 / 16.00d+00)
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 6, 8) = d2_b * g_x20(g_i_)
        enddo
        hmv(6, 8) = -(3.00d+00 / 16.00d+00) * x20
C--------
        d4_b = -(9.00d+00 / 16.00d+00)
        do g_i_ = 1, g_p_
          g_hmv(g_i_, 6, 9) = -g_hmv(g_i_, 6, 7) + d4_b * g_x20(g_i_)
        enddo
        hmv(6, 9) = -(9.00d+00 / 16.00d+00) * x20 - hmv(6, 7)
C--------
C
C.....CONSTRUCT LOCAL MATRIX [HQV] FROM [HMV]
C
        do 99955 j = 1, 9
C
          d2_b = 2.00d+00 / 9.00d+00
          do g_i_ = 1, g_p_
            g_sum123(g_i_) = d2_b * g_hmv(g_i_, 3, j) + d2_b * g_hmv(g_i
     *_, 2, j) + d2_b * g_hmv(g_i_, 1, j)
          enddo
          sum123 = 2.00d+00 / 9.00d+00 * (hmv(1, j) + hmv(2, j) + hmv(3,
     * j))
C--------
C
          d4_b = -(4.00d+00 / 3.00d+00)
          do g_i_ = 1, g_p_
            g_hqv(g_i_, 1, j) = d4_b * g_hmv(g_i_, 1, j) + g_sum123(g_i_
     *)
          enddo
          hqv(1, j) = sum123 - 4.00d+00 / 3.00d+00 * hmv(1, j)
C--------
          d4_b = -(4.00d+00 / 3.00d+00)
          do g_i_ = 1, g_p_
            g_hqv(g_i_, 2, j) = d4_b * g_hmv(g_i_, 2, j) + g_sum123(g_i_
     *)
          enddo
          hqv(2, j) = sum123 - 4.00d+00 / 3.00d+00 * hmv(2, j)
C--------
          d4_b = -(4.00d+00 / 3.00d+00)
          do g_i_ = 1, g_p_
            g_hqv(g_i_, 3, j) = d4_b * g_hmv(g_i_, 3, j) + g_sum123(g_i_
     *)
          enddo
          hqv(3, j) = sum123 - 4.00d+00 / 3.00d+00 * hmv(3, j)
C--------
C
          d2_b = 2.00d+00 / 9.00d+00
          do g_i_ = 1, g_p_
            g_sum456(g_i_) = d2_b * g_hmv(g_i_, 6, j) + d2_b * g_hmv(g_i
     *_, 5, j) + d2_b * g_hmv(g_i_, 4, j)
          enddo
          sum456 = 2.00d+00 / 9.00d+00 * (hmv(4, j) + hmv(5, j) + hmv(6,
     * j))
C--------
C
          d4_b = -(4.00d+00 / 3.00d+00)
          do g_i_ = 1, g_p_
            g_hqv(g_i_, 4, j) = d4_b * g_hmv(g_i_, 4, j) + g_sum456(g_i_
     *)
          enddo
          hqv(4, j) = sum456 - 4.00d+00 / 3.00d+00 * hmv(4, j)
C--------
          d4_b = -(4.00d+00 / 3.00d+00)
          do g_i_ = 1, g_p_
            g_hqv(g_i_, 5, j) = d4_b * g_hmv(g_i_, 5, j) + g_sum456(g_i_
     *)
          enddo
          hqv(5, j) = sum456 - 4.00d+00 / 3.00d+00 * hmv(5, j)
C--------
          d4_b = -(4.00d+00 / 3.00d+00)
          do g_i_ = 1, g_p_
            g_hqv(g_i_, 6, j) = d4_b * g_hmv(g_i_, 6, j) + g_sum456(g_i_
     *)
          enddo
          hqv(6, j) = sum456 - 4.00d+00 / 3.00d+00 * hmv(6, j)
C--------
C
5001      continue
99955   continue
C
C.....FORM THE LOCAL MATRIX [B]
C
        do g_i_ = 1, g_p_
          g_b(g_i_, 1, 1) = g_y30(g_i_)
        enddo
        b(1, 1) = y30
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 2, 1) = 0.0d0
        enddo
        b(2, 1) = zero
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 3, 1) = -g_x30(g_i_)
        enddo
        b(3, 1) = -x30
C--------
C
        do g_i_ = 1, g_p_
          g_b(g_i_, 1, 2) = g_y10(g_i_)
        enddo
        b(1, 2) = y10
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 2, 2) = 0.0d0
        enddo
        b(2, 2) = zero
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 3, 2) = -g_x10(g_i_)
        enddo
        b(3, 2) = -x10
C--------
C
        do g_i_ = 1, g_p_
          g_b(g_i_, 1, 3) = g_y20(g_i_)
        enddo
        b(1, 3) = y20
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 2, 3) = 0.0d0
        enddo
        b(2, 3) = zero
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 3, 3) = -g_x20(g_i_)
        enddo
        b(3, 3) = -x20
C--------
C
        do g_i_ = 1, g_p_
          g_b(g_i_, 1, 4) = 0.0d0
        enddo
        b(1, 4) = zero
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 2, 4) = -g_x30(g_i_)
        enddo
        b(2, 4) = -x30
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 3, 4) = g_y30(g_i_)
        enddo
        b(3, 4) = y30
C--------
C
        do g_i_ = 1, g_p_
          g_b(g_i_, 1, 5) = 0.0d0
        enddo
        b(1, 5) = zero
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 2, 5) = -g_x10(g_i_)
        enddo
        b(2, 5) = -x10
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 3, 5) = g_y10(g_i_)
        enddo
        b(3, 5) = y10
C--------
C
        do g_i_ = 1, g_p_
          g_b(g_i_, 1, 6) = 0.0d0
        enddo
        b(1, 6) = zero
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 2, 6) = -g_x20(g_i_)
        enddo
        b(2, 6) = -x20
C--------
        do g_i_ = 1, g_p_
          g_b(g_i_, 3, 6) = g_y20(g_i_)
        enddo
        b(3, 6) = y20
C--------
C
C.....CALCULATE THE DIAGONAL TRIANGULAR COORDINATE MATRIX [Z1]
C.....AT THE FIRST GAUSS INTEGRATION POINT (W/ MID-POINT RULE)
C
        z1(1) = gauss21 - gauss11
        z1(2) = gauss31 - gauss21
        z1(3) = gauss11 - gauss31
        z1(4) = gauss21 - gauss11
        z1(5) = gauss31 - gauss21
        z1(6) = gauss11 - gauss31
C
C.....CALCULATE THE DIAGONAL TRIANGULAR COORDINATE MATRIX [Z2]
C.....AT THE SECOND GAUSS INTEGRATION POINT (W/ MID-POINT RULE)
C
        z2(1) = gauss22 - gauss12
        z2(2) = gauss32 - gauss22
        z2(3) = gauss12 - gauss32
        z2(4) = gauss22 - gauss12
        z2(5) = gauss32 - gauss22
        z2(6) = gauss12 - gauss32
C
C.....CALCULATE THE DIAGONAL TRIANGULAR COORDINATE MATRIX [Z3]
C.....AT THE THIRD GAUSS INTEGRATION POINT (W/ MID-POINT RULE)
C
        z3(1) = gauss23 - gauss13
        z3(2) = gauss33 - gauss23
        z3(3) = gauss13 - gauss33
        z3(4) = gauss23 - gauss13
        z3(5) = gauss33 - gauss23
        z3(6) = gauss13 - gauss33
C
C.....FORM THE PRODUCTS [Z1]*[B]^T, [Z2]*[B]^T AND [Z3]*[B]^T
C.....IN MATRICES [z1bt], [z2bt] AND [z3bt], RESPECTIVELY
C
        do 99951 j = 1, 3
C
          do 99954 i = 1, 6
            do g_i_ = 1, g_p_
              g_z1bt(g_i_, i, j) = z1(i) * g_b(g_i_, j, i)
            enddo
            z1bt(i, j) = z1(i) * b(j, i)
C--------
5003        continue
99954     continue
C
          do 99953 i = 1, 6
            do g_i_ = 1, g_p_
              g_z2bt(g_i_, i, j) = z2(i) * g_b(g_i_, j, i)
            enddo
            z2bt(i, j) = z2(i) * b(j, i)
C--------
5004        continue
99953     continue
C
          do 99952 i = 1, 6
            do g_i_ = 1, g_p_
              g_z3bt(g_i_, i, j) = z3(i) * g_b(g_i_, j, i)
            enddo
            z3bt(i, j) = z3(i) * b(j, i)
C--------
5006        continue
99952     continue
C
5002      continue
99951   continue
C
C.....MULTIPLY LOCAL MATRICES [z1bt], [z2bt] AND [z3bt] BY THE SQUARE
C.....ROOT OF ( 9.0 * [fm] * INTEGRATION_WEIGHT / AREA )
C
        d2_v = area * invweight1
        d3_v = 9.00d+00 * fm / d2_v
        d3_b = -d3_v / d2_v * invweight1
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_area(g_i_)
        enddo
        d1_w = d3_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_comphBM.f',
     +2359)
        endif
        do g_i_ = 1, g_p_
          g_factor1(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        factor1 = d2_v
C--------
        d2_v = area * invweight2
        d3_v = 9.00d+00 * fm / d2_v
        d3_b = -d3_v / d2_v * invweight2
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_area(g_i_)
        enddo
        d1_w = d3_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_comphBM.f',
     +2379)
        endif
        do g_i_ = 1, g_p_
          g_factor2(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        factor2 = d2_v
C--------
        d2_v = area * invweight3
        d3_v = 9.00d+00 * fm / d2_v
        d3_b = -d3_v / d2_v * invweight3
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_area(g_i_)
        enddo
        d1_w = d3_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_comphBM.f',
     +2399)
        endif
        do g_i_ = 1, g_p_
          g_factor3(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        factor3 = d2_v
C--------
C
        do 99947 j = 1, 3
C
          do 99950 i = 1, 6
            do g_i_ = 1, g_p_
              g_z1bt(g_i_, i, j) = factor1 * g_z1bt(g_i_, i, j) + z1bt(i
     *, j) * g_factor1(g_i_)
            enddo
            z1bt(i, j) = factor1 * z1bt(i, j)
C--------
5008        continue
99950     continue
C
          do 99949 i = 1, 6
            do g_i_ = 1, g_p_
              g_z2bt(g_i_, i, j) = factor2 * g_z2bt(g_i_, i, j) + z2bt(i
     *, j) * g_factor2(g_i_)
            enddo
            z2bt(i, j) = factor2 * z2bt(i, j)
C--------
5009        continue
99949     continue
C
          do 99948 i = 1, 6
            do g_i_ = 1, g_p_
              g_z3bt(g_i_, i, j) = factor3 * g_z3bt(g_i_, i, j) + z3bt(i
     *, j) * g_factor3(g_i_)
            enddo
            z3bt(i, j) = factor3 * z3bt(i, j)
C--------
5010        continue
99948     continue
C
5007      continue
99947   continue
C
C.....FORM THE LOCAL PRODUCTS:
C.....[Ph1] = [HQV]^T * [Z1] * [B]^T = [HQV]^T * [z1bt]
C.....[Ph2] = [HQV]^T * [Z2] * [B]^T = [HQV]^T * [z2bt]
C.....[Ph3] = [HQV]^T * [Z3] * [B]^T = [HQV]^T * [z3bt]
C
        do 99943 i = 1, 9
C
          do g_i_ = 1, g_p_
            g_h1(g_i_) = g_hqv(g_i_, 1, i)
          enddo
          h1 = hqv(1, i)
C--------
          do g_i_ = 1, g_p_
            g_h2(g_i_) = g_hqv(g_i_, 2, i)
          enddo
          h2 = hqv(2, i)
C--------
          do g_i_ = 1, g_p_
            g_h3(g_i_) = g_hqv(g_i_, 3, i)
          enddo
          h3 = hqv(3, i)
C--------
          do g_i_ = 1, g_p_
            g_h4(g_i_) = g_hqv(g_i_, 4, i)
          enddo
          h4 = hqv(4, i)
C--------
          do g_i_ = 1, g_p_
            g_h5(g_i_) = g_hqv(g_i_, 5, i)
          enddo
          h5 = hqv(5, i)
C--------
          do g_i_ = 1, g_p_
            g_h6(g_i_) = g_hqv(g_i_, 6, i)
          enddo
          h6 = hqv(6, i)
C--------
C
          do 99946 j = 1, 3
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = h2 * g_z1bt(g_i_, 2, j) + z1bt(2, j) * g_h2
     *(g_i_) + h1 * g_z1bt(g_i_, 1, j) + z1bt(1, j) * g_h1(g_i_)
            enddo
            d1_w = h1 * z1bt(1, j) + h2 * z1bt(2, j)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = h4 * g_z1bt(g_i_, 4, j) + z1bt(4, j) * g_h4
     *(g_i_) + h3 * g_z1bt(g_i_, 3, j) + z1bt(3, j) * g_h3(g_i_) + g_d1_
     *w(g_i_)
            enddo
            d2_w = d1_w + h3 * z1bt(3, j) + h4 * z1bt(4, j)
            do g_i_ = 1, g_p_
              g_ph1(g_i_, i, j) = h6 * g_z1bt(g_i_, 6, j) + z1bt(6, j) *
     * g_h6(g_i_) + h5 * g_z1bt(g_i_, 5, j) + z1bt(5, j) * g_h5(g_i_) + 
     *g_d2_w(g_i_)
            enddo
            ph1(i, j) = d2_w + h5 * z1bt(5, j) + h6 * z1bt(6, j)
C--------
5012        continue
99946     continue
C
          do 99945 j = 1, 3
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = h2 * g_z2bt(g_i_, 2, j) + z2bt(2, j) * g_h2
     *(g_i_) + h1 * g_z2bt(g_i_, 1, j) + z2bt(1, j) * g_h1(g_i_)
            enddo
            d1_w = h1 * z2bt(1, j) + h2 * z2bt(2, j)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = h4 * g_z2bt(g_i_, 4, j) + z2bt(4, j) * g_h4
     *(g_i_) + h3 * g_z2bt(g_i_, 3, j) + z2bt(3, j) * g_h3(g_i_) + g_d1_
     *w(g_i_)
            enddo
            d2_w = d1_w + h3 * z2bt(3, j) + h4 * z2bt(4, j)
            do g_i_ = 1, g_p_
              g_ph2(g_i_, i, j) = h6 * g_z2bt(g_i_, 6, j) + z2bt(6, j) *
     * g_h6(g_i_) + h5 * g_z2bt(g_i_, 5, j) + z2bt(5, j) * g_h5(g_i_) + 
     *g_d2_w(g_i_)
            enddo
            ph2(i, j) = d2_w + h5 * z2bt(5, j) + h6 * z2bt(6, j)
C--------
5013        continue
99945     continue
C
          do 99944 j = 1, 3
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = h2 * g_z3bt(g_i_, 2, j) + z3bt(2, j) * g_h2
     *(g_i_) + h1 * g_z3bt(g_i_, 1, j) + z3bt(1, j) * g_h1(g_i_)
            enddo
            d1_w = h1 * z3bt(1, j) + h2 * z3bt(2, j)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = h4 * g_z3bt(g_i_, 4, j) + z3bt(4, j) * g_h4
     *(g_i_) + h3 * g_z3bt(g_i_, 3, j) + z3bt(3, j) * g_h3(g_i_) + g_d1_
     *w(g_i_)
            enddo
            d2_w = d1_w + h3 * z3bt(3, j) + h4 * z3bt(4, j)
            do g_i_ = 1, g_p_
              g_ph3(g_i_, i, j) = h6 * g_z3bt(g_i_, 6, j) + z3bt(6, j) *
     * g_h6(g_i_) + h5 * g_z3bt(g_i_, 5, j) + z3bt(5, j) * g_h5(g_i_) + 
     *g_d2_w(g_i_)
            enddo
            ph3(i, j) = d2_w + h5 * z3bt(5, j) + h6 * z3bt(6, j)
C--------
5014        continue
99944     continue
C
5011      continue
99943   continue
C
C     --------------------------------------------------------------------
C
C     STIFFNESS MATRIX ASSEMBLY FOR HIGHER ORDER BENDING-MEMBRANE COUPLING
C
C     --------------------------------------------------------------------
C
900     continue
C
C.....ASSEMBLE THE OUTPUT STIFFNESS SUCH THAT:
C.....[khBM] = [Lh1]*[dbm]*[Ph1]^T + [Lh2]*[dbm]*[Ph2]^T + [Lh3]*[dbm]*[Ph3]^T
C
        do 99941 i = 1, 9
C
          col = colm(i)
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dbm(1, 2) * g_ph1(g_i_, i, 2) + ph1(i, 2) * g
     *_dbm(g_i_, 1, 2) + dbm(1, 1) * g_ph1(g_i_, i, 1) + ph1(i, 1) * g_d
     *bm(g_i_, 1, 1)
          enddo
          d1_w = dbm(1, 1) * ph1(i, 1) + dbm(1, 2) * ph1(i, 2)
          do g_i_ = 1, g_p_
            g_mult11(g_i_) = dbm(1, 3) * g_ph1(g_i_, i, 3) + ph1(i, 3) *
     * g_dbm(g_i_, 1, 3) + g_d1_w(g_i_)
          enddo
          mult11 = d1_w + dbm(1, 3) * ph1(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dbm(2, 2) * g_ph1(g_i_, i, 2) + ph1(i, 2) * g
     *_dbm(g_i_, 2, 2) + dbm(2, 1) * g_ph1(g_i_, i, 1) + ph1(i, 1) * g_d
     *bm(g_i_, 2, 1)
          enddo
          d1_w = dbm(2, 1) * ph1(i, 1) + dbm(2, 2) * ph1(i, 2)
          do g_i_ = 1, g_p_
            g_mult12(g_i_) = dbm(2, 3) * g_ph1(g_i_, i, 3) + ph1(i, 3) *
     * g_dbm(g_i_, 2, 3) + g_d1_w(g_i_)
          enddo
          mult12 = d1_w + dbm(2, 3) * ph1(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dbm(3, 2) * g_ph1(g_i_, i, 2) + ph1(i, 2) * g
     *_dbm(g_i_, 3, 2) + dbm(3, 1) * g_ph1(g_i_, i, 1) + ph1(i, 1) * g_d
     *bm(g_i_, 3, 1)
          enddo
          d1_w = dbm(3, 1) * ph1(i, 1) + dbm(3, 2) * ph1(i, 2)
          do g_i_ = 1, g_p_
            g_mult13(g_i_) = dbm(3, 3) * g_ph1(g_i_, i, 3) + ph1(i, 3) *
     * g_dbm(g_i_, 3, 3) + g_d1_w(g_i_)
          enddo
          mult13 = d1_w + dbm(3, 3) * ph1(i, 3)
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dbm(1, 2) * g_ph2(g_i_, i, 2) + ph2(i, 2) * g
     *_dbm(g_i_, 1, 2) + dbm(1, 1) * g_ph2(g_i_, i, 1) + ph2(i, 1) * g_d
     *bm(g_i_, 1, 1)
          enddo
          d1_w = dbm(1, 1) * ph2(i, 1) + dbm(1, 2) * ph2(i, 2)
          do g_i_ = 1, g_p_
            g_mult21(g_i_) = dbm(1, 3) * g_ph2(g_i_, i, 3) + ph2(i, 3) *
     * g_dbm(g_i_, 1, 3) + g_d1_w(g_i_)
          enddo
          mult21 = d1_w + dbm(1, 3) * ph2(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dbm(2, 2) * g_ph2(g_i_, i, 2) + ph2(i, 2) * g
     *_dbm(g_i_, 2, 2) + dbm(2, 1) * g_ph2(g_i_, i, 1) + ph2(i, 1) * g_d
     *bm(g_i_, 2, 1)
          enddo
          d1_w = dbm(2, 1) * ph2(i, 1) + dbm(2, 2) * ph2(i, 2)
          do g_i_ = 1, g_p_
            g_mult22(g_i_) = dbm(2, 3) * g_ph2(g_i_, i, 3) + ph2(i, 3) *
     * g_dbm(g_i_, 2, 3) + g_d1_w(g_i_)
          enddo
          mult22 = d1_w + dbm(2, 3) * ph2(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dbm(3, 2) * g_ph2(g_i_, i, 2) + ph2(i, 2) * g
     *_dbm(g_i_, 3, 2) + dbm(3, 1) * g_ph2(g_i_, i, 1) + ph2(i, 1) * g_d
     *bm(g_i_, 3, 1)
          enddo
          d1_w = dbm(3, 1) * ph2(i, 1) + dbm(3, 2) * ph2(i, 2)
          do g_i_ = 1, g_p_
            g_mult23(g_i_) = dbm(3, 3) * g_ph2(g_i_, i, 3) + ph2(i, 3) *
     * g_dbm(g_i_, 3, 3) + g_d1_w(g_i_)
          enddo
          mult23 = d1_w + dbm(3, 3) * ph2(i, 3)
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dbm(1, 2) * g_ph3(g_i_, i, 2) + ph3(i, 2) * g
     *_dbm(g_i_, 1, 2) + dbm(1, 1) * g_ph3(g_i_, i, 1) + ph3(i, 1) * g_d
     *bm(g_i_, 1, 1)
          enddo
          d1_w = dbm(1, 1) * ph3(i, 1) + dbm(1, 2) * ph3(i, 2)
          do g_i_ = 1, g_p_
            g_mult31(g_i_) = dbm(1, 3) * g_ph3(g_i_, i, 3) + ph3(i, 3) *
     * g_dbm(g_i_, 1, 3) + g_d1_w(g_i_)
          enddo
          mult31 = d1_w + dbm(1, 3) * ph3(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dbm(2, 2) * g_ph3(g_i_, i, 2) + ph3(i, 2) * g
     *_dbm(g_i_, 2, 2) + dbm(2, 1) * g_ph3(g_i_, i, 1) + ph3(i, 1) * g_d
     *bm(g_i_, 2, 1)
          enddo
          d1_w = dbm(2, 1) * ph3(i, 1) + dbm(2, 2) * ph3(i, 2)
          do g_i_ = 1, g_p_
            g_mult32(g_i_) = dbm(2, 3) * g_ph3(g_i_, i, 3) + ph3(i, 3) *
     * g_dbm(g_i_, 2, 3) + g_d1_w(g_i_)
          enddo
          mult32 = d1_w + dbm(2, 3) * ph3(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dbm(3, 2) * g_ph3(g_i_, i, 2) + ph3(i, 2) * g
     *_dbm(g_i_, 3, 2) + dbm(3, 1) * g_ph3(g_i_, i, 1) + ph3(i, 1) * g_d
     *bm(g_i_, 3, 1)
          enddo
          d1_w = dbm(3, 1) * ph3(i, 1) + dbm(3, 2) * ph3(i, 2)
          do g_i_ = 1, g_p_
            g_mult33(g_i_) = dbm(3, 3) * g_ph3(g_i_, i, 3) + ph3(i, 3) *
     * g_dbm(g_i_, 3, 3) + g_d1_w(g_i_)
          enddo
          mult33 = d1_w + dbm(3, 3) * ph3(i, 3)
C--------
C
          do 99942 j = 1, 9
            row = rowb(j)
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = mult12 * g_lh1(g_i_, j, 2) + lh1(j, 2) * g_
     *mult12(g_i_) + mult11 * g_lh1(g_i_, j, 1) + lh1(j, 1) * g_mult11(g
     *_i_) + g_khbm(g_i_, row, col)
            enddo
            d1_w = khbm(row, col) + mult11 * lh1(j, 1) + mult12 * lh1(j,
     * 2)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = mult21 * g_lh2(g_i_, j, 1) + lh2(j, 1) * g_
     *mult21(g_i_) + mult13 * g_lh1(g_i_, j, 3) + lh1(j, 3) * g_mult13(g
     *_i_) + g_d1_w(g_i_)
            enddo
            d2_w = d1_w + mult13 * lh1(j, 3) + mult21 * lh2(j, 1)
            do g_i_ = 1, g_p_
              g_d3_w(g_i_) = mult23 * g_lh2(g_i_, j, 3) + lh2(j, 3) * g_
     *mult23(g_i_) + mult22 * g_lh2(g_i_, j, 2) + lh2(j, 2) * g_mult22(g
     *_i_) + g_d2_w(g_i_)
            enddo
            d3_w = d2_w + mult22 * lh2(j, 2) + mult23 * lh2(j, 3)
            do g_i_ = 1, g_p_
              g_d4_w(g_i_) = mult32 * g_lh3(g_i_, j, 2) + lh3(j, 2) * g_
     *mult32(g_i_) + mult31 * g_lh3(g_i_, j, 1) + lh3(j, 1) * g_mult31(g
     *_i_) + g_d3_w(g_i_)
            enddo
            d4_w = d3_w + mult31 * lh3(j, 1) + mult32 * lh3(j, 2)
            do g_i_ = 1, g_p_
              g_khbm(g_i_, row, col) = mult33 * g_lh3(g_i_, j, 3) + lh3(
     *j, 3) * g_mult33(g_i_) + g_d4_w(g_i_)
            enddo
            khbm(row, col) = d4_w + mult33 * lh3(j, 3)
C--------
6002        continue
99942     continue
C
6001      continue
99941   continue
C
C.....OUTPUT THE MATRIX PRIOR TO ROTATION (FOR DEBUGGING ONLY)
C
C     open(unit=90,file="khBM.m")
C     write(90,*) "khBM=["
C     do 991 i=1,18
C 991 write(90,9) (khBM(i,j),j=1,18)
C     write(90,*) "     ];"
C     close(90)
C   9 format(18(1x,E16.9))
C
C.....ROTATE THE OUTPUT STIFFNESS MATRIX
C
C      call compmrot( khBM , rot , rot , rot )
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
C.....ERROR-MESSAGE IF THE STIFFNESS FACTOR [FB] (BENDING) IS NEGATIVE
C
100     continue
        write (*, *) '*** FATAL ERROR in routine COMPHBM          ***'
        write (*, *) '*** The Stiffness Factor [fb] for Bending   ***'
        write (*, *) '*** is Negative: Check the Calling Sequence ***'
        write (*, *) '*** Factor [fb] Must be Positive or Zero    ***'
        write (*, *) '*** EXECUTION TERNINATED RIGHT HERE         ***'
        stop
C
C.....ERROR-MESSAGE IF THE STIFFNESS FACTOR [FM] (MEMBRANE) IS NEGATIVE
C
200     continue
        write (*, *) '*** FATAL ERROR in routine COMPHBM          ***'
        write (*, *) '*** The Stiffness Factor [fm] for Membrane  ***'
        write (*, *) '*** is Negative: Check the Calling Sequence ***'
        write (*, *) '*** Factor [fm] Must be Positive or Zero    ***'
        write (*, *) '*** EXECUTION TERNINATED RIGHT HERE         ***'
        stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
300     continue
        write (*, *) '*** FATAL ERROR in routine COMPHBM         ***'
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
C=end of routine "COMPHBM"
C========================C
