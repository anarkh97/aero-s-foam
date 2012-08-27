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
      subroutine g_comphbb(g_p_, elm, type, x, g_x, ldg_x, y, g_y, ldg_y
     *, db, g_db, ldg_db, f, rowb, colb, rot, g_rot, ldg_rot, lh1, g_lh1
     *, ldg_lh1, lh2, g_lh2, ldg_lh2, lh3, g_lh3, ldg_lh3, khbb, g_khbb,
     * ldg_khbb)
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
        integer elm, type, rowb(9), colb(9)
        real*8 x(3), y(3), db(3, 3), khbb(18, 18), f
        real*8 lh1(9, 3), lh2(9, 3), lh3(9, 3), rot(6, 6)
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
        real*8 gauss11, gauss21, gauss31
        real*8 gauss12, gauss22, gauss32
        real*8 gauss13, gauss23, gauss33
        real*8 dblh1t1, dblh1t2, dblh1t3
        real*8 dblh2t1, dblh2t2, dblh2t3
        real*8 dblh3t1, dblh3t2, dblh3t3
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_khbb, ldg_lh1, ldg_lh2, ldg_lh3, ldg_x, 
     *ldg_y, ldg_db, ldg_rot
        double precision d4_w, d3_w, d5_b, d6_b, d8_b, d7_b, d2_w, d1_p,
     * d2_v, d3_v
        double precision d4_v, d5_v, d6_v, d1_w, d2_b, d3_b, d4_b, g_khb
     *b(ldg_khbb, 18, 18), g_q(g_pmax_, 6, 9), g_sq(g_pmax_, 3, 3)
        double precision g_rm1(g_pmax_, 3, 6), g_rm2(g_pmax_, 3, 6), g_r
     *m3(g_pmax_, 3, 6), g_rm1tsqt(g_pmax_, 6, 3), g_rm2tsqt(g_pmax_, 6,
     * 3), g_rm3tsqt(g_pmax_, 6, 3), g_lh1(ldg_lh1, 9, 3), g_lh2(ldg_lh2
     *, 9, 3), g_lh3(ldg_lh3, 9, 3), g_x0(g_pmax_)
        double precision g_x(ldg_x, 3), g_y0(g_pmax_), g_y(ldg_y, 3), g_
     *x1(g_pmax_), g_x2(g_pmax_), g_x3(g_pmax_), g_y1(g_pmax_), g_y2(g_p
     *max_), g_y3(g_pmax_), g_x21(g_pmax_)
        double precision g_x32(g_pmax_), g_x13(g_pmax_), g_y21(g_pmax_),
     * g_y32(g_pmax_), g_y13(g_pmax_), g_twicearea(g_pmax_), g_area(g_pm
     *ax_), g_d1_w(g_pmax_), g_dist21(g_pmax_), g_dist32(g_pmax_)
        double precision g_dist13(g_pmax_), g_d2_w(g_pmax_), g_bl1(g_pma
     *x_), g_bl2(g_pmax_), g_bl3(g_pmax_), g_al1(g_pmax_), g_al2(g_pmax_
     *), g_al3(g_pmax_), g_x2ap3(g_pmax_), g_factor(g_pmax_)
        double precision g_s1(g_pmax_), g_s2(g_pmax_), g_s3(g_pmax_), g_
     *q1(g_pmax_), g_q2(g_pmax_), g_q3(g_pmax_), g_q4(g_pmax_), g_q5(g_p
     *max_), g_q6(g_pmax_), g_db(ldg_db, 3, 3)
        double precision g_dblh1t1(g_pmax_), g_dblh1t2(g_pmax_), g_dblh1
     *t3(g_pmax_), g_dblh2t1(g_pmax_), g_dblh2t2(g_pmax_), g_dblh2t3(g_p
     *max_), g_dblh3t1(g_pmax_), g_dblh3t2(g_pmax_), g_dblh3t3(g_pmax_),
     * g_d3_w(g_pmax_)
        double precision g_d4_w(g_pmax_), g_rot(ldg_rot, 6, 6)
        save g_dblh2t2, g_dblh2t3, g_dblh3t1, g_dblh3t2, g_dblh3t3, g_d3
     *_w, g_d4_w
        save g_q1, g_q2, g_q3, g_q4, g_q5, g_q6, g_dblh1t1, g_dblh1t2, g
     *_dblh1t3, g_dblh2t1
        save g_bl2, g_bl3, g_al1, g_al2, g_al3, g_x2ap3, g_factor, g_s1,
     * g_s2, g_s3
        save g_y32, g_y13, g_twicearea, g_area, g_d1_w, g_dist21, g_dist
     *32, g_dist13, g_d2_w, g_bl1
        save g_x1, g_x2, g_x3, g_y1, g_y2, g_y3, g_x21, g_x32, g_x13, g_
     *y21
        save g_q, g_sq, g_rm1, g_rm2, g_rm3, g_rm1tsqt, g_rm2tsqt, g_rm3
     *tsqt, g_x0, g_y0
        data zero /0.000000d+00/
C
C.....DEFINE THE FIRST SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
        data gauss11 /0.000000d+00/
        data gauss21 /0.500000d+00/
        data gauss31 /0.500000d+00/
C
C.....DEFINE THE SECOND SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
        data gauss12 /0.500000d+00/
        data gauss22 /0.000000d+00/
        data gauss32 /0.500000d+00/
C
C.....DEFINE THE THIRD SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
        data gauss13 /0.500000d+00/
        data gauss23 /0.500000d+00/
        data gauss33 /0.000000d+00/
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
              g_khbb(g_i_, i, j) = 0.0d0
            enddo
            khbb(i, j) = zero
C--------
1002        continue
99999     continue
1001      continue
99998   continue
C
C.....CLEAR THE LOCAL MATRICES
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
        do 99982 j = 1, 3
          do 99985 i = 1, 9
            do g_i_ = 1, g_p_
              g_lh1(g_i_, i, j) = 0.0d0
            enddo
            lh1(i, j) = zero
C--------
1016        continue
99985     continue
          do 99984 i = 1, 9
            do g_i_ = 1, g_p_
              g_lh2(g_i_, i, j) = 0.0d0
            enddo
            lh2(i, j) = zero
C--------
1017        continue
99984     continue
          do 99983 i = 1, 9
            do g_i_ = 1, g_p_
              g_lh3(g_i_, i, j) = 0.0d0
            enddo
            lh3(i, j) = zero
C--------
1018        continue
99983     continue
1015      continue
99982   continue
C
C.....RETURN IF THE STIFFNESS FACTOR IS ZERO
C
        if (f .eq. zero) then
          return
        endif
C
C.....CHECK IF THE STIFFNESS FACTOR [F] IS POSITIVE
C
        if (f .lt. zero) then
          goto 100
        endif
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
          goto 200
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
     +'g_comphBB.f',
     +389)
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
     +'g_comphBB.f',
     +408)
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
     +'g_comphBB.f',
     +427)
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
        d3_b = 1.0d0 / 3.00d+00 * f
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_area(g_i_)
        enddo
        d1_w = f * area / 3.00d+00
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_comphBB.f',
     +1001)
        endif
        do g_i_ = 1, g_p_
          g_factor(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        factor = d2_v
C--------
C
        do 99980 j = 1, 3
          do 99981 i = 1, 3
            do g_i_ = 1, g_p_
              g_sq(g_i_, i, j) = factor * g_sq(g_i_, i, j) + sq(i, j) * 
     *g_factor(g_i_)
            enddo
            sq(i, j) = factor * sq(i, j)
C--------
2002        continue
99981     continue
2001      continue
99980   continue
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
        do 99976 j = 1, 3
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
          do 99979 i = 1, 6
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
99979     continue
C
          do 99978 i = 1, 6
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
99978     continue
C
          do 99977 i = 1, 6
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
99977     continue
C
3001      continue
99976   continue
C
C.....ASSEMBLE MATRICES [Lh1], [Lh2] AND [Lh3] FOR THE THREE GAUSS POINTS:
C.....[Lh1] = [q]^T * [rm1]^T * [sq]^T = [q]^T * [rm1tsqt]
C.....[Lh2] = [q]^T * [rm2]^T * [sq]^T = [q]^T * [rm2tsqt]
C.....[Lh3] = [q]^T * [rm3]^T * [sq]^T = [q]^T * [rm3tsqt]
C
        do 99972 i = 1, 9
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
          do 99975 j = 1, 3
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
99975     continue
C
          do 99974 j = 1, 3
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
99974     continue
C
          do 99973 j = 1, 3
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
99973     continue
C
4001      continue
99972   continue
C
C.....ASSEMBLE THE OUTPUT STIFFNESS SUCH THAT:
C.....[khBB] = [Lh1]*[db]*[Lh1]^T + [Lh2]*[db]*[Lh2]^T + [Lh3]*[db]*[Lh3]^T
C
        do 99970 i = 1, 9
C
          col = colb(i)
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(1, 2) * g_lh1(g_i_, i, 2) + lh1(i, 2) * g_
     *db(g_i_, 1, 2) + db(1, 1) * g_lh1(g_i_, i, 1) + lh1(i, 1) * g_db(g
     *_i_, 1, 1)
          enddo
          d1_w = db(1, 1) * lh1(i, 1) + db(1, 2) * lh1(i, 2)
          do g_i_ = 1, g_p_
            g_dblh1t1(g_i_) = db(1, 3) * g_lh1(g_i_, i, 3) + lh1(i, 3) *
     * g_db(g_i_, 1, 3) + g_d1_w(g_i_)
          enddo
          dblh1t1 = d1_w + db(1, 3) * lh1(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(2, 2) * g_lh1(g_i_, i, 2) + lh1(i, 2) * g_
     *db(g_i_, 2, 2) + db(2, 1) * g_lh1(g_i_, i, 1) + lh1(i, 1) * g_db(g
     *_i_, 2, 1)
          enddo
          d1_w = db(2, 1) * lh1(i, 1) + db(2, 2) * lh1(i, 2)
          do g_i_ = 1, g_p_
            g_dblh1t2(g_i_) = db(2, 3) * g_lh1(g_i_, i, 3) + lh1(i, 3) *
     * g_db(g_i_, 2, 3) + g_d1_w(g_i_)
          enddo
          dblh1t2 = d1_w + db(2, 3) * lh1(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(3, 2) * g_lh1(g_i_, i, 2) + lh1(i, 2) * g_
     *db(g_i_, 3, 2) + db(3, 1) * g_lh1(g_i_, i, 1) + lh1(i, 1) * g_db(g
     *_i_, 3, 1)
          enddo
          d1_w = db(3, 1) * lh1(i, 1) + db(3, 2) * lh1(i, 2)
          do g_i_ = 1, g_p_
            g_dblh1t3(g_i_) = db(3, 3) * g_lh1(g_i_, i, 3) + lh1(i, 3) *
     * g_db(g_i_, 3, 3) + g_d1_w(g_i_)
          enddo
          dblh1t3 = d1_w + db(3, 3) * lh1(i, 3)
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(1, 2) * g_lh2(g_i_, i, 2) + lh2(i, 2) * g_
     *db(g_i_, 1, 2) + db(1, 1) * g_lh2(g_i_, i, 1) + lh2(i, 1) * g_db(g
     *_i_, 1, 1)
          enddo
          d1_w = db(1, 1) * lh2(i, 1) + db(1, 2) * lh2(i, 2)
          do g_i_ = 1, g_p_
            g_dblh2t1(g_i_) = db(1, 3) * g_lh2(g_i_, i, 3) + lh2(i, 3) *
     * g_db(g_i_, 1, 3) + g_d1_w(g_i_)
          enddo
          dblh2t1 = d1_w + db(1, 3) * lh2(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(2, 2) * g_lh2(g_i_, i, 2) + lh2(i, 2) * g_
     *db(g_i_, 2, 2) + db(2, 1) * g_lh2(g_i_, i, 1) + lh2(i, 1) * g_db(g
     *_i_, 2, 1)
          enddo
          d1_w = db(2, 1) * lh2(i, 1) + db(2, 2) * lh2(i, 2)
          do g_i_ = 1, g_p_
            g_dblh2t2(g_i_) = db(2, 3) * g_lh2(g_i_, i, 3) + lh2(i, 3) *
     * g_db(g_i_, 2, 3) + g_d1_w(g_i_)
          enddo
          dblh2t2 = d1_w + db(2, 3) * lh2(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(3, 2) * g_lh2(g_i_, i, 2) + lh2(i, 2) * g_
     *db(g_i_, 3, 2) + db(3, 1) * g_lh2(g_i_, i, 1) + lh2(i, 1) * g_db(g
     *_i_, 3, 1)
          enddo
          d1_w = db(3, 1) * lh2(i, 1) + db(3, 2) * lh2(i, 2)
          do g_i_ = 1, g_p_
            g_dblh2t3(g_i_) = db(3, 3) * g_lh2(g_i_, i, 3) + lh2(i, 3) *
     * g_db(g_i_, 3, 3) + g_d1_w(g_i_)
          enddo
          dblh2t3 = d1_w + db(3, 3) * lh2(i, 3)
C--------
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(1, 2) * g_lh3(g_i_, i, 2) + lh3(i, 2) * g_
     *db(g_i_, 1, 2) + db(1, 1) * g_lh3(g_i_, i, 1) + lh3(i, 1) * g_db(g
     *_i_, 1, 1)
          enddo
          d1_w = db(1, 1) * lh3(i, 1) + db(1, 2) * lh3(i, 2)
          do g_i_ = 1, g_p_
            g_dblh3t1(g_i_) = db(1, 3) * g_lh3(g_i_, i, 3) + lh3(i, 3) *
     * g_db(g_i_, 1, 3) + g_d1_w(g_i_)
          enddo
          dblh3t1 = d1_w + db(1, 3) * lh3(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(2, 2) * g_lh3(g_i_, i, 2) + lh3(i, 2) * g_
     *db(g_i_, 2, 2) + db(2, 1) * g_lh3(g_i_, i, 1) + lh3(i, 1) * g_db(g
     *_i_, 2, 1)
          enddo
          d1_w = db(2, 1) * lh3(i, 1) + db(2, 2) * lh3(i, 2)
          do g_i_ = 1, g_p_
            g_dblh3t2(g_i_) = db(2, 3) * g_lh3(g_i_, i, 3) + lh3(i, 3) *
     * g_db(g_i_, 2, 3) + g_d1_w(g_i_)
          enddo
          dblh3t2 = d1_w + db(2, 3) * lh3(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = db(3, 2) * g_lh3(g_i_, i, 2) + lh3(i, 2) * g_
     *db(g_i_, 3, 2) + db(3, 1) * g_lh3(g_i_, i, 1) + lh3(i, 1) * g_db(g
     *_i_, 3, 1)
          enddo
          d1_w = db(3, 1) * lh3(i, 1) + db(3, 2) * lh3(i, 2)
          do g_i_ = 1, g_p_
            g_dblh3t3(g_i_) = db(3, 3) * g_lh3(g_i_, i, 3) + lh3(i, 3) *
     * g_db(g_i_, 3, 3) + g_d1_w(g_i_)
          enddo
          dblh3t3 = d1_w + db(3, 3) * lh3(i, 3)
C--------
C
          do 99971 j = 1, i
            row = rowb(j)
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = dblh1t2 * g_lh1(g_i_, j, 2) + lh1(j, 2) * g
     *_dblh1t2(g_i_) + dblh1t1 * g_lh1(g_i_, j, 1) + lh1(j, 1) * g_dblh1
     *t1(g_i_) + g_khbb(g_i_, row, col)
            enddo
            d1_w = khbb(row, col) + dblh1t1 * lh1(j, 1) + dblh1t2 * lh1(
     *j, 2)
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = dblh2t1 * g_lh2(g_i_, j, 1) + lh2(j, 1) * g
     *_dblh2t1(g_i_) + dblh1t3 * g_lh1(g_i_, j, 3) + lh1(j, 3) * g_dblh1
     *t3(g_i_) + g_d1_w(g_i_)
            enddo
            d2_w = d1_w + dblh1t3 * lh1(j, 3) + dblh2t1 * lh2(j, 1)
            do g_i_ = 1, g_p_
              g_d3_w(g_i_) = dblh2t3 * g_lh2(g_i_, j, 3) + lh2(j, 3) * g
     *_dblh2t3(g_i_) + dblh2t2 * g_lh2(g_i_, j, 2) + lh2(j, 2) * g_dblh2
     *t2(g_i_) + g_d2_w(g_i_)
            enddo
            d3_w = d2_w + dblh2t2 * lh2(j, 2) + dblh2t3 * lh2(j, 3)
            do g_i_ = 1, g_p_
              g_d4_w(g_i_) = dblh3t2 * g_lh3(g_i_, j, 2) + lh3(j, 2) * g
     *_dblh3t2(g_i_) + dblh3t1 * g_lh3(g_i_, j, 1) + lh3(j, 1) * g_dblh3
     *t1(g_i_) + g_d3_w(g_i_)
            enddo
            d4_w = d3_w + dblh3t1 * lh3(j, 1) + dblh3t2 * lh3(j, 2)
            do g_i_ = 1, g_p_
              g_khbb(g_i_, row, col) = dblh3t3 * g_lh3(g_i_, j, 3) + lh3
     *(j, 3) * g_dblh3t3(g_i_) + g_d4_w(g_i_)
            enddo
            khbb(row, col) = d4_w + dblh3t3 * lh3(j, 3)
C--------
            do g_i_ = 1, g_p_
              g_khbb(g_i_, col, row) = g_khbb(g_i_, row, col)
            enddo
            khbb(col, row) = khbb(row, col)
C--------
5002        continue
99971     continue
C
5001      continue
99970   continue
C
C.....OUTPUT THE MATRIX PRIOR TO ROTATION (FOR DEBUGGING ONLY)
C
C     open(unit=90,file="khBB.m")
C     write(90,*) "khBB=["
C     do 991 i=1,18
C 991 write(90,9) (khBB(i,j),j=1,18)
C     write(90,*) "     ];"
C     close(90)
C   9 format(18(1x,E16.9))
C
C.....ROTATE THE OUTPUT STIFFNESS MATRIX
C
C      call compmrot( khBB , rot , rot , rot )
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
C.....ERROR-MESSAGE IF THE STIFFNESS FACTOR [F] IS NEGATIVE
C
100     continue
        write (*, *) '*** FATAL ERROR in routine COMPHBB       ***'
        write (*, *) '*** The Stiffness Factor [f] is Negative ***'
        write (*, *) '*** Check the Calling Sequence:          ***'
        write (*, *) '*** Factor [f] Must be Positive or Zero  ***'
        write (*, *) '*** EXECUTION TERNINATED RIGHT HERE      ***'
        stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
200     continue
        write (*, *) '*** FATAL ERROR in routine COMPHBB         ***'
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
C=end of routine "COMPHBB"
C========================C
