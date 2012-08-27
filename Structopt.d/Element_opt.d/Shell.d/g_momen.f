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
C=block fortran
      subroutine g_momen(x, g_x, y, g_y, lb, l, g_l,
     * status)
C
C                   a r g u m e n t s
C
C       implicit none
       integer lb(*)
        real*8 x(*), y(*), l(18, *)
        character*(*) status
C
C                   t y p e   &   d i m e n s i o n
C
        real*8 lqr(9, 3)
        real*8 x0, y0, cab, a1, a2, a3, b1, b2, b3
        real*8 x21, x32, x13, y21, y32, y13
        real*8 x12, x23, x31, y12, y23, y31
        real*8 xl12, xl23, xl31, c12, c23, c31, s12, s23, s31
        real*8 cc12, cc23, cc31, ss12, ss23, ss31
        real*8 cs12, cs23, cs31
        real*8 area2
        integer i, j, kk
C
C                   l o g i c
C
        double precision d9_b, d5_b, d8_b, d4_b, d2_p, d1_p, d1_w,
     * d7_b, d2_v, d3_v
        double precision d6_b, d2_b, d3_b, d4_v, g_x21, 
     *g_x(*), g_x12, g_x32, g_x23, g_x13
        double precision g_x31, g_y21, g_y(*), 
     *g_y12, g_y32, g_y23, g_y13, 
     *g_y31, g_area2, g_d1_w
        double precision g_xl12, g_xl23, g_xl31
     *, g_lqr(9, 3), g_c12, g_s12, 
     *g_c23, g_s23, g_c31, g_s31
        double precision g_cc12, g_cc23, g_cc31
     *, g_ss12, g_ss23, g_ss31, g_cs12
     *, g_cs23, g_cs31, g_l(18, *)
        save g_cs12, g_cs23, g_cs31
        save g_c23, g_s23, g_c31, g_s31, g_cc12, g_cc23, g_cc31,
     * g_ss12, g_ss23, g_ss31
        save g_y13, g_y31, g_area2, g_d1_w, g_xl12, g_xl23, g_xl31,
     * g_lqr, g_c12, g_s12
        save g_x21, g_x12, g_x32, g_x23, g_x13, g_x31, g_y21, g_y12,
     * g_y32, g_y23
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'momen','g_momen.f')
C
        status = ' '
          g_x21 = -g_x(1) + g_x(2)
        x21 = x(2) - x(1)
C--------
          g_x12 = -g_x21
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
          g_y12 = -g_y21
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
          g_area2 = -x21 * g_y13 + (-y13) * g_x21 + y2
     *1 * g_x13 + x13 * g_y21
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
        d2_v = x12 * x12
        d2_p = 2.0d0 * x12
        d4_v = y12 * y12
        d1_p = 2.0d0 * y12
          g_d1_w = d1_p * g_y12 + d2_p * g_x12
        d1_w = d2_v + d4_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +171)
        endif
          g_xl12 = d1_p * g_d1_w
        xl12 = d2_v
C--------
        d2_v = x23 * x23
        d2_p = 2.0d0 * x23
        d4_v = y23 * y23
        d1_p = 2.0d0 * y23
          g_d1_w = d1_p * g_y23 + d2_p * g_x23
        d1_w = d2_v + d4_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +192)
        endif
          g_xl23 = d1_p * g_d1_w
        xl23 = d2_v
C--------
        d2_v = x31 * x31
        d2_p = 2.0d0 * x31
        d4_v = y31 * y31
        d1_p = 2.0d0 * y31
          g_d1_w = d1_p * g_y31 + d2_p * g_x31
        d1_w = d2_v + d4_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +213)
        endif
          g_xl31 = d1_p * g_d1_w
        xl31 = d2_v
C--------
        do 99998 j = 1, 3
          do 99999 i = 1, 9
              g_lqr(i, j) = 0.0d0
            lqr(i, j) = 0.0d0
C--------
1200        continue
99999     continue
99998   continue
C
C
        d3_v = y21 / xl12
        d2_b = 1.0d0 / xl12
        d3_b = -d3_v / xl12
          g_c12 = d3_b * g_xl12 + d2_b * g_y21
        c12 = d3_v
C--------
        d3_v = x12 / xl12
        d2_b = 1.0d0 / xl12
        d3_b = -d3_v / xl12
          g_s12 = d3_b * g_xl12 + d2_b * g_x12
        s12 = d3_v
C--------
        d3_v = y32 / xl23
        d2_b = 1.0d0 / xl23
        d3_b = -d3_v / xl23
          g_c23 = d3_b * g_xl23 + d2_b * g_y32
        c23 = d3_v
C--------
        d3_v = x23 / xl23
        d2_b = 1.0d0 / xl23
        d3_b = -d3_v / xl23
          g_s23 = d3_b * g_xl23 + d2_b * g_x23
        s23 = d3_v
C--------
        d3_v = y13 / xl31
        d2_b = 1.0d0 / xl31
        d3_b = -d3_v / xl31
          g_c31 = d3_b * g_xl31 + d2_b * g_y13
        c31 = d3_v
C--------
        d3_v = x31 / xl31
        d2_b = 1.0d0 / xl31
        d3_b = -d3_v / xl31
          g_s31 = d3_b * g_xl31 + d2_b * g_x31
        s31 = d3_v
C--------
        d2_b = c12 + c12
          g_cc12 = d2_b * g_c12
        cc12 = c12 * c12
C--------
        d2_b = c23 + c23
          g_cc23 = d2_b * g_c23
        cc23 = c23 * c23
C--------
        d2_b = c31 + c31
          g_cc31 = d2_b * g_c31
        cc31 = c31 * c31
C--------
        d2_b = s12 + s12
          g_ss12 = d2_b * g_s12 
        ss12 = s12 * s12
C--------
        d2_b = s23 + s23
          g_ss23 = d2_b * g_s23
        ss23 = s23 * s23
C--------
        d2_b = s31 + s31
          g_ss31 = d2_b * g_s31
        ss31 = s31 * s31
C--------
          g_cs12 = c12 * g_s12 + s12 * g_c12
        cs12 = c12 * s12
C--------
          g_cs23 = c23 * g_s23 + s23 * g_c23
        cs23 = c23 * s23
C--------
          g_cs31 = c31 * g_s31 + s31 * g_c31
        cs31 = c31 * s31
C--------
          g_lqr(1, 1) = -g_cs31 + g_cs12
        lqr(1, 1) = cs12 - cs31
C--------
          g_lqr(1, 2) = -g_lqr(1, 1)
        lqr(1, 2) = -lqr(1, 1)
C--------
          g_lqr(1, 3) = g_ss12 + (-g_cc12) + 
     *(-g_ss31) + g_cc31
        lqr(1, 3) = cc31 - ss31 - (cc12 - ss12)
C--------
        d2_b = dble(.5)
        d5_b = d2_b * x31
        d6_b = d2_b * cc31
        d7_b = d2_b * x12
        d8_b = d2_b * cc12
          g_lqr(2, 1) = d6_b * g_x31 + d5_b * g_cc31 +
     * d8_b * g_x12 + d7_b * g_cc12
        lqr(2, 1) = (cc12 * x12 + cc31 * x31) * dble(.5)
C--------
        d2_b = dble(.5)
        d5_b = d2_b * x31
        d6_b = d2_b * ss31
        d7_b = d2_b * x12
        d8_b = d2_b * ss12
          g_lqr(2, 2) = d6_b * g_x31 + d5_b * g_ss31 +
     * d8_b * g_x12 + d7_b * g_ss12
        lqr(2, 2) = (ss12 * x12 + ss31 * x31) * dble(.5)
C--------
          g_lqr(2, 3) = ss31 * g_y13 + y13 * g_ss31 + 
     *ss12 * g_y21 + y21 * g_ss12
        lqr(2, 3) = ss12 * y21 + ss31 * y13
C--------
        d2_b = dble(.5)
        d6_b = -d2_b * y13
        d7_b = -d2_b * cc31
        d8_b = -d2_b * y21
        d9_b = -d2_b * cc12
          g_lqr(3, 1) = d7_b * g_y13 + d6_b * g_cc31 +
     * d9_b * g_y21 + d8_b * g_cc12
        lqr(3, 1) = -(cc12 * y21 + cc31 * y13) * dble(.5)
C--------
        d2_b = dble(-0.5)
          g_lqr(3, 2) = d2_b * g_lqr(2, 3)
        lqr(3, 2) = dble(-0.5) * lqr(2, 3)
C--------
        d2_b = dble(-2.)       
          g_lqr(3, 3) = d2_b * g_lqr(2, 1)
        lqr(3, 3) = dble(-2.) * lqr(2, 1)
C--------
          g_lqr(4, 1) = -g_cs12 + g_cs23
        lqr(4, 1) = cs23 - cs12
C--------
          g_lqr(4, 2) = -g_lqr(4, 1)       
        lqr(4, 2) = -lqr(4, 1)
C--------
          g_lqr(4, 3) = g_ss23 + (-g_cc23) + 
     *(-g_ss12) + g_cc12
        lqr(4, 3) = cc12 - ss12 - (cc23 - ss23)
C--------
        d2_b = dble(.5)
        d5_b = d2_b * x23
        d6_b = d2_b * cc23
        d7_b = d2_b * x12
        d8_b = d2_b * cc12
          g_lqr(5, 1) = d6_b * g_x23 + d5_b * g_cc23 +
     * d8_b * g_x12 + d7_b * g_cc12
        lqr(5, 1) = (cc12 * x12 + cc23 * x23) * dble(.5)
C--------
        d2_b = dble(.5)
        d5_b = d2_b * x23
        d6_b = d2_b * ss23
        d7_b = d2_b * x12
        d8_b = d2_b * ss12
          g_lqr(5, 2) = d6_b * g_x23 + d5_b * g_ss23 +
     * d8_b * g_x12 + d7_b * g_ss12
        lqr(5, 2) = (ss12 * x12 + ss23 * x23) * dble(.5)
C--------
          g_lqr(5, 3) = ss23 * g_y32 + y32 * g_ss23 + 
     *ss12 * g_y21 + y21 * g_ss12
        lqr(5, 3) = ss12 * y21 + ss23 * y32
C--------
        d2_b = dble(.5)
        d6_b = -d2_b * y32
        d7_b = -d2_b * cc23
        d8_b = -d2_b * y21
        d9_b = -d2_b * cc12
          g_lqr(6, 1) = d7_b * g_y32 + d6_b * g_cc23 +
     * d9_b * g_y21 + d8_b * g_cc12
        lqr(6, 1) = -(cc12 * y21 + cc23 * y32) * dble(.5)
C--------
        d2_b = dble(-0.5)
          g_lqr(6, 2) = d2_b * g_lqr(5, 3)
        lqr(6, 2) = dble(-0.5) * lqr(5, 3)
C--------
        d2_b = dble(-2.)
          g_lqr(6, 3) = d2_b * g_lqr(5, 1)
        lqr(6, 3) = dble(-2.) * lqr(5, 1)
C--------
          g_lqr(7, 1) = -g_cs23 + g_cs31
        lqr(7, 1) = cs31 - cs23
C--------
          g_lqr(7, 2) = -g_lqr(7, 1)
        lqr(7, 2) = -lqr(7, 1)
C--------
          g_lqr(7, 3) = g_ss31 + (-g_cc31) + 
     *(-g_ss23) + g_cc23
        lqr(7, 3) = cc23 - ss23 - (cc31 - ss31)
C--------
        d2_b = dble(.5)
        d5_b = d2_b * x31
        d6_b = d2_b * cc31
        d7_b = d2_b * x23
        d8_b = d2_b * cc23
          g_lqr(8, 1) = d6_b * g_x31 + d5_b * g_cc31 +
     * d8_b * g_x23 + d7_b * g_cc23
        lqr(8, 1) = (cc23 * x23 + cc31 * x31) * dble(.5)
C--------
        d2_b = dble(.5)
        d5_b = d2_b * x31
        d6_b = d2_b * ss31
        d7_b = d2_b * x23
        d8_b = d2_b * ss23
          g_lqr(8, 2) = d6_b * g_x31 + d5_b * g_ss31 +
     * d8_b * g_x23 + d7_b * g_ss23
        lqr(8, 2) = (ss23 * x23 + ss31 * x31) * dble(.5)
C--------
          g_lqr(8, 3) = ss31 * g_y13 + y13 * g_ss31 + 
     *ss23 * g_y32 + y32 * g_ss23
        lqr(8, 3) = ss23 * y32 + ss31 * y13
C--------
        d2_b = dble(.5)
        d6_b = -d2_b * y13
        d7_b = -d2_b * cc31
        d8_b = -d2_b * y32
        d9_b = -d2_b * cc23
          g_lqr(9, 1) = d7_b * g_y13 + d6_b * g_cc31 +
     * d9_b * g_y32 + d8_b * g_cc23
        lqr(9, 1) = -(cc23 * y32 + cc31 * y13) * dble(.5)
C--------
        d2_b = dble(-0.5)
          g_lqr(9, 2) = d2_b * g_lqr(8, 3)
        lqr(9, 2) = dble(-0.5) * lqr(8, 3)
C--------
        d2_b = dble(-2.)
          g_lqr(9, 3) = d2_b * g_lqr(8, 1)
        lqr(9, 3) = dble(-2.) * lqr(8, 1)
C--------
C
        do 99997 j = 1, 9
          kk = lb(j)
C       l(kk,1) = lqr(j,1)*2./area2
C       l(kk,2) = lqr(j,2)*2./area2
C       l(kk,3) = lqr(j,3)*2./area2
          d2_v = 2.00d+00 / area2
          d2_b = -d2_v / area2
            g_d1_w = d2_b * g_area2
          d1_w = d2_v
          d3_v = sqrt(d1_w)
          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d3_v)
          else
             call ehufDO (9,d1_w, d3_v, d1_p,
     +g_ehfid,
     +550)
          endif
          d4_b = lqr(j, 1) * d1_p
            g_l(kk, 1) = d4_b * g_d1_w + d3_v * 
     *g_lqr(j, 1)
          l(kk, 1) = lqr(j, 1) * d3_v
C--------
          d2_v = 2.00d+00 / area2
          d2_b = -d2_v / area2
            g_d1_w = d2_b * g_area2
          d1_w = d2_v
          d3_v = sqrt(d1_w)
          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d3_v)
          else
             call ehufDO (9,d1_w, d3_v, d1_p,
     +g_ehfid,
     +571)
          endif
          d4_b = lqr(j, 2) * d1_p
            g_l(kk, 2) = d4_b * g_d1_w + d3_v * 
     *g_lqr(j, 2)
          l(kk, 2) = lqr(j, 2) * d3_v
C--------
          d2_v = 2.00d+00 / area2
          d2_b = -d2_v / area2
            g_d1_w = d2_b * g_area2
          d1_w = d2_v
          d3_v = sqrt(d1_w)
          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d3_v)
          else
             call ehufDO (9,d1_w, d3_v, d1_p,
     +g_ehfid,
     +592)
          endif
          d4_b = lqr(j, 3) * d1_p
            g_l(kk, 3) = d4_b * g_d1_w + d3_v * 
     *g_lqr(j, 3)
          l(kk, 3) = lqr(j, 3) * d3_v
C--------
1600      continue
99997   continue
        return
      end
C=end fortran
