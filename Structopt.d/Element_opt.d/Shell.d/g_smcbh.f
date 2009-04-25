      subroutine g_smcbh( x, g_x, y, g_y, db, g_db, f, ls, sm, g_sm, 
     &                    m, status)
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
        double precision d3_w, d14_b, d14_v, d4_p, d3_p, d13_b, d12_b, 
     &  d10_b, d2_v, d3_v
        double precision d4_v, d5_v, d6_v, d13_v, d2_b, d3_b, d4_b, 
     &  d5_b, d6_b, d7_v
        double precision d7_b, d1_w, d1_p, d2_p, d2_w, d8_v, d8_b, 
     &  d9_v, d10_v, d9_b
        double precision g_rm( 3, 6), g_q( 6, 9), g_rsd( 6, 6), 
     &  g_sds( 3, 3), g_x0, g_x( 3), g_y0, g_y(  3), g_x1, g_x2
        double precision g_x3, g_y1, g_y2, g_y3, g_x21, g_x32, g_x13, 
     &  g_y21, g_y32, g_y13
        double precision g_area2, g_d1_w, g_l21, g_l32, g_l13, 
     &  g_d2_w, g_bl2, g_al2, g_bl3, g_al3
        double precision g_bl1, g_al1, g_cc, g_sq( 3, 3), g_d11, 
     &  g_db( 3, 3), g_d22, g_d33, g_d12, g_d13
        double precision g_d23, g_area, g_s1, g_s2, g_s3, g_s4, 
     &  g_s5, g_s6, g_d3_w, g_sm( m, m)
        save g_s2, g_s3, g_s4, g_s5, g_s6, g_d3_w
        save g_cc, g_sq, g_d11, g_d22, g_d33, g_d12, g_d13, 
     &  g_d23, g_area, g_s1
        save g_l21, g_l32, g_l13, g_d2_w, g_bl2, g_al2, g_bl3, 
     &  g_al3, g_bl1, g_al1
        save g_y2, g_y3, g_x21, g_x32, g_x13, g_y21, g_y32, 
     &  g_y13, g_area2, g_d1_w
        save g_rm, g_q, g_rsd, g_sds, g_x0, g_y0, g_x1, g_x2, 
     &  g_x3, g_y1
        intrinsic dble
        data pg /0.d0, 0.5d0, 0.5d0, 0.5d0, 0.d0, 0.5d0, 
     &  0.5d0, 0.5d0, 0.d0/
C
C                   l o g i c
C
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'smcbh','g_smcbh.f')
C
        status = ' '
C     cleaning

        do 99998 i = 1, 3
          do 99999 j = 1, 6
              g_rm( i, j) = 0.0d0
            rm(i, j) = 0.0d0
C--------
100         continue
99999     continue
99998   continue
        do 99996 i = 1, 6
          do 99997 j = 1, 9
              g_q( i, j) = 0.0d0
            q(i, j) = 0.0d0
C--------
200         continue
99997     continue
99996   continue
        do 99994 i = 1, 6
          do 99995 j = 1, 6
              g_rsd( i, j) = 0.0d0
            rsd(i, j) = 0.0d0
C--------
99995     continue
99994   continue
        do 99992 i = 1, 3
          do 99993 j = 1, 3
              g_sds( i, j) = 0.0d0
            sds(i, j) = 0.0d0
C--------
99993     continue
99992   continue
C     coordinates
        d2_b = 1.0d0 / dble(3.)
          g_x0 = d2_b * g_x( 3) + d2_b * g_x( 2) + d2_b * g_x( 1)
        x0 = (x(1) + x(2) + x(3)) / dble(3.)
C--------
        d2_b = 1.0d0 / dble(3.)
          g_y0 = d2_b * g_y(  3) + d2_b * g_y(  2) + d2_b * g_y(  1)
        y0 = (y(1) + y(2) + y(3)) / dble(3.)
C--------
          g_x1 = -g_x0 + g_x( 1)
        x1 = x(1) - x0
C--------
          g_x2 = -g_x0 + g_x( 2)
        x2 = x(2) - x0
C--------
          g_x3 = -g_x0 + g_x( 3)
        x3 = x(3) - x0
C--------
          g_y1 = -g_y0 + g_y(  1)
        y1 = y(1) - y0
C--------
          g_y2 = -g_y0 + g_y(  2)
        y2 = y(2) - y0
C--------
          g_y3 = -g_y0 + g_y(  3)
        y3 = y(3) - y0
C--------
          g_x21 = -g_x1 + g_x2
        x21 = x2 - x1
C--------
          g_x32 = -g_x2 + g_x3
        x32 = x3 - x2
C--------
          g_x13 = -g_x3 + g_x1
        x13 = x1 - x3
C--------
          g_y21 = -g_y1 + g_y2
        y21 = y2 - y1
C--------
          g_y32 = -g_y2 + g_y3
        y32 = y3 - y2
C--------
          g_y13 = -g_y3 + g_y1
        y13 = y1 - y3
C--------
          g_area2 = -x21 * g_y13 + (-y13) * g_x21 + y21 * g_x13 + 
     &  x13 * g_y21
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
          g_d1_w = d1_p * g_y21 + d2_p * g_x21
        d1_w = d2_v + d4_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +270)
        endif
          g_l21 = d1_p * g_d1_w
        l21 = d2_v
C--------
        d2_v = x32 * x32
        d2_p = 2.0d0 * x32
        d4_v = y32 * y32
        d1_p = 2.0d0 * y32
          g_d1_w = d1_p * g_y32 + d2_p * g_x32
        d1_w = d2_v + d4_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +291)
        endif
          g_l32 = d1_p * g_d1_w
        l32 = d2_v
C--------
        d2_v = x13 * x13
        d2_p = 2.0d0 * x13
        d4_v = y13 * y13
        d1_p = 2.0d0 * y13
          g_d1_w = d1_p * g_y13 + d2_p * g_x13
        d1_w = d2_v + d4_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +312)
        endif
          g_l13 = d1_p * g_d1_w
        l13 = d2_v
C--------
C          side proyections
        d3_v = x3 - x1
        d5_v = x2 - x1
        d5_b = -d3_v + (-d5_v)
          g_d1_w = d3_v * g_x2 + d5_b * g_x1 + d5_v * g_x3
        d1_w = d3_v * d5_v
        d4_v = y3 - y1
        d6_v = y2 - y1
        d7_b = -d4_v + (-d6_v)
          g_d2_w = d4_v * g_y2 + d7_b * g_y1 + d6_v * g_y3 + g_d1_w
        d2_w = d1_w + d4_v * d6_v
        d3_v = l21 * l21
        d1_p = 2.0d0 * l21
        d4_v = d2_w / d3_v
        d2_b = 1.0d0 / d3_v
        d4_b = -d4_v / d3_v * d1_p
          g_bl2 = d4_b * g_l21 + d2_b * g_d2_w
        bl2 = d4_v
C--------
          g_al2 = -g_bl2
        al2 = 1.0d0 - bl2
C--------
        d3_v = x1 - x2
        d5_v = x3 - x2
        d5_b = -d3_v + (-d5_v)
          g_d1_w = d3_v * g_x3 + d5_b * g_x2 + d5_v * g_x1
        d1_w = d3_v * d5_v
        d4_v = y1 - y2
        d6_v = y3 - y2
        d7_b = -d4_v + (-d6_v)
          g_d2_w = d4_v * g_y3 + d7_b * g_y2 + d6_v * g_y1 + g_d1_w
        d2_w = d1_w + d4_v * d6_v
        d3_v = l32 * l32
        d1_p = 2.0d0 * l32
        d4_v = d2_w / d3_v
        d2_b = 1.0d0 / d3_v
        d4_b = -d4_v / d3_v * d1_p
          g_bl3 = d4_b * g_l32 + d2_b * g_d2_w
        bl3 = d4_v
C--------
          g_al3 = -g_bl3
        al3 = 1.0d0 - bl3
C--------
        d3_v = x2 - x3
        d5_v = x1 - x3
        d5_b = -d3_v + (-d5_v)
          g_d1_w = d3_v * g_x1 + d5_b * g_x3 + d5_v * g_x2
        d1_w = d3_v * d5_v
        d4_v = y2 - y3
        d6_v = y1 - y3
        d7_b = -d4_v + (-d6_v)
          g_d2_w = d4_v * g_y1 + d7_b * g_y3 + d6_v * g_y2 + g_d1_w
        d2_w = d1_w + d4_v * d6_v
        d3_v = l13 * l13
        d1_p = 2.0d0 * l13
        d4_v = d2_w / d3_v
        d2_b = 1.0d0 / d3_v
        d4_b = -d4_v / d3_v * d1_p
          g_bl1 = d4_b * g_l13 + d2_b * g_d2_w
        bl1 = d4_v
C--------
          g_al1 = -g_bl1
        al1 = 1.0d0 - bl1
C--------
C          inverse of the matrix relating inside curvatures
C          xx,yy,xy with boundary curvatures
C
        d2_v = area2 ** ( 3 - 2)
        d2_v =  d2_v * area2
        d1_p =  3 *  d2_v
        d2_v =  d2_v * area2
          g_cc = d1_p * g_area2
        cc = d2_v
C--------
        d4_v = -x21 * y21
        d6_v = y32 * y32
        d1_p = 2.0d0 * y32
        d4_b = d4_v * d1_p
        d6_b = d6_v * (-x21)
        d7_b = -(d6_v * y21)
          g_d1_w = d4_b * g_y32 + d6_b * g_y21 + d7_b * g_x21
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
          g_sq( 1, 1) = d3_b * g_cc + d7_b * g_y32 + d10_b * g_y21 + 
     &  d8_b * g_x32 + d2_b * g_d1_w
        sq(1, 1) = d10_v
C--------
        d3_v = x13 * y13
        d5_v = y32 * y32
        d1_p = 2.0d0 * y32
        d4_b = d3_v * d1_p
        d5_b = d5_v * y13
        d6_b = d5_v * x13
          g_d1_w = d4_b * g_y32 + d6_b * g_y13 + d5_b * g_x13
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
          g_sq( 1, 2) = d3_b * g_cc + d7_b * g_y32 + d10_b * g_y13 + 
     &  d8_b * g_x32 + d2_b * g_d1_w
        sq(1, 2) = d10_v
C--------
        d3_v = x21 * y21
        d5_v = y13 * y13
        d1_p = 2.0d0 * y13
        d4_b = d3_v * d1_p
        d5_b = d5_v * y21
        d6_b = d5_v * x21
          g_d1_w = d4_b * g_y13 + d6_b * g_y21 + d5_b * g_x21
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
          g_sq( 1, 3) = d3_b * g_cc + d7_b * g_y13 + d10_b * g_y21 + 
     &  d8_b * g_x13 + d2_b * g_d1_w
        sq(1, 3) = d10_v
C--------
        d3_v = x32 * x32
        d1_p = 2.0d0 * x32
        d4_v = x21 * d3_v
        d4_b = y21 * d3_v
        d6_b = y21 * x21 * d1_p
          g_d1_w = d4_v * g_y21 + d6_b * g_x32 + d4_b * g_x21
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
          g_sq( 2, 1) = d3_b * g_cc + d7_b * g_y32 + d9_b * g_x32 + 
     &    d10_b * g_x21 + d2_b * g_d1_w
        sq(2, 1) = d10_v
C--------
        d4_v = x32 * x32
        d1_p = 2.0d0 * x32
        d5_v = -x13 * d4_v
        d6_b = y13 * (-x13) * d1_p
        d7_b = -(y13 * d4_v)
          g_d1_w = d5_v * g_y13 + d6_b * g_x32 + d7_b * g_x13
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
          g_sq( 2, 2) = d3_b * g_cc + d7_b * g_y32 + d9_b * g_x32 + 
     &    d10_b * g_x13 + d2_b * g_d1_w
        sq(2, 2) = d10_v
C--------
        d4_v = x13 * x13
        d1_p = 2.0d0 * x13
        d5_v = -x21 * d4_v
        d6_b = y21 * (-x21) * d1_p
        d7_b = -(y21 * d4_v)
          g_d1_w = d5_v * g_y21 + d6_b * g_x13 + d7_b * g_x21
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
          g_sq( 2, 3) = d3_b * g_cc + d7_b * g_y13 + d9_b * g_x13 +
     &     d10_b * g_x21 + d2_b * g_d1_w
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
          g_sq( 3, 1) = d3_b * g_cc + d8_b * g_y21 + d9_b * g_x32 + 
     &    d12_b * g_y32 + d13_b * g_x21
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
          g_sq( 3, 2) = d3_b * g_cc + d8_b * g_y13 + d9_b * g_x32 + 
     &    d12_b * g_y32 + d14_b * g_x13
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
          g_sq( 3, 3) = d3_b * g_cc + d8_b * g_y21 + d9_b * g_x13 + 
     &    d12_b * g_y13 + d14_b * g_x21
        sq(3, 3) = d14_v
C--------
C     print '(''area'',f8.3)',area
C     print '(''inver'',3f8.3)',((sq(i,j),j=1,3),i=1,3)
          g_d11 = g_db( 1, 1)
        d11 = db(1, 1)
C--------
          g_d22 = g_db( 2, 2)
        d22 = db(2, 2)
C--------
          g_d33 = g_db( 3, 3)
        d33 = db(3, 3)
C--------
          g_d12 = g_db( 1, 2)
        d12 = db(1, 2)
C--------
          g_d13 = g_db( 1, 3)
        d13 = db(1, 3)
C--------
          g_d23 = g_db( 2, 3)
        d23 = db(2, 3)
C--------
        d2_b = dble(0.5)
          g_area = d2_b * g_area2
        area = dble(0.5) * area2
C--------
        do 99990 j = 1, 3
            g_d1_w = d12 * g_sq( 2, j) + sq(2, j) * g_d12 + d11 * 
     &    g_sq( 1, j) + sq(1, j) * g_d11
          d1_w = d11 * sq(1, j) + d12 * sq(2, j)
            g_s1 = d13 * g_sq( 3, j) + sq(3, j) * g_d13 + g_d1_w
          s1 = d1_w + d13 * sq(3, j)
C--------
            g_d1_w = d22 * g_sq( 2, j) + sq(2, j) * g_d22 + d12 * 
     &    g_sq( 1, j) + sq(1, j) * g_d12
          d1_w = d12 * sq(1, j) + d22 * sq(2, j)
            g_s2 = d23 * g_sq( 3, j) + sq(3, j) * g_d23 + g_d1_w
          s2 = d1_w + d23 * sq(3, j)
C--------
            g_d1_w = d23 * g_sq( 2, j) + sq(2, j) * g_d23 + d13 * 
     &      g_sq( 1, j) + sq(1, j) * g_d13
          d1_w = d13 * sq(1, j) + d23 * sq(2, j)
            g_s3 = d33 * g_sq( 3, j) + sq(3, j) * g_d33 + g_d1_w
          s3 = d1_w + d33 * sq(3, j)
C--------
          do 99991 i = 1, j
              g_d1_w = s2 * g_sq( 2, i) + sq(2, i) * g_s2 + s1 * 
     &    g_sq( 1, i) + sq(1, i) * g_s1
            d1_w = s1 * sq(1, i) + s2 * sq(2, i)
              g_sds( i, j) = s3 * g_sq( 3, i) + sq(3, i) * g_s3 + 
     &    g_d1_w + g_sds( i, j)
            sds(i, j) = sds(i, j) + (d1_w + s3 * sq(3, i))
C--------
              g_sds( j, i) = g_sds( i, j)
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
              g_sds( i, j) = d5_b * g_area + d4_b * g_sds( i, j)
            sds(i, j) = sds(i, j) * area / dble(3.0) * f
C--------
99989     continue
99988   continue
C     print '(''sds'',3f8.3)',((sds(i,j),j=1,3),i=1,3)
C
C    
C         matrix q
          g_q( 1, 1) = 0.0d0
        q(1, 1) = 6.0d0
C--------
        d2_b = dble(-2.0)
          g_q( 1, 2) = d2_b * g_y13
        q(1, 2) = dble(-2.0) * y13
C--------
        d2_b = dble(2.0)
          g_q( 1, 3) = d2_b * g_x13
        q(1, 3) = dble(2.0) * x13
C--------
          g_q( 1, 7) = 0.0d0
        q(1, 7) = -6.0d0
C--------
        d2_b = dble(-4.0)
          g_q( 1, 8) = d2_b * g_y13
        q(1, 8) = dble(-4.0) * y13
C--------
        d2_b = dble(4.0)
          g_q( 1, 9) = d2_b * g_x13
        q(1, 9) = dble(4.0) * x13
C--------
          g_q( 2, 1) = 0.0d0
        q(2, 1) = -6.0d0
C--------
        d2_b = dble(4.0)
          g_q( 2, 2) = d2_b * g_y13
        q(2, 2) = dble(4.0) * y13
C--------
        d2_b = dble(-4.0)
          g_q( 2, 3) = d2_b * g_x13
        q(2, 3) = dble(-4.0) * x13
C--------
          g_q( 2, 7) = 0.0d0
        q(2, 7) = 6.0d0
C--------
        d2_b = dble(2.0)
          g_q( 2, 8) = d2_b * g_y13
        q(2, 8) = dble(2.0) * y13
C--------
        d2_b = dble(-2.0)
          g_q( 2, 9) = d2_b * g_x13
        q(2, 9) = dble(-2.0) * x13
C--------
          g_q( 3, 1) = 0.0d0
        q(3, 1) = -6.0d0
C--------
        d2_b = dble(-4.0)
          g_q( 3, 2) = d2_b * g_y21
        q(3, 2) = dble(-4.0) * y21
C--------
        d2_b = dble(4.0)
          g_q( 3, 3) = d2_b * g_x21
        q(3, 3) = dble(4.0) * x21
C--------
          g_q( 3, 4) = 0.0d0
        q(3, 4) = 6.0d0
C--------
        d2_b = dble(-2.0)
          g_q( 3, 5) = d2_b * g_y21
        q(3, 5) = dble(-2.0) * y21
C--------
        d2_b = dble(2.0)
          g_q( 3, 6) = d2_b * g_x21
        q(3, 6) = dble(2.0) * x21
C--------
          g_q( 4, 1) = 0.0d0
        q(4, 1) = 6.0d0
C--------
        d2_b = dble(2.0)
          g_q( 4, 2) = d2_b * g_y21
        q(4, 2) = dble(2.0) * y21
C--------
        d2_b = dble(-2.0)
          g_q( 4, 3) = d2_b * g_x21
        q(4, 3) = dble(-2.0) * x21
C--------
          g_q( 4, 4) = 0.0d0
        q(4, 4) = -6.0d0
C--------
        d2_b = dble(4.0)
          g_q( 4, 5) = d2_b * g_y21
        q(4, 5) = dble(4.0) * y21
C--------
        d2_b = dble(-4.0)
          g_q( 4, 6) = d2_b * g_x21
        q(4, 6) = dble(-4.0) * x21
C--------
          g_q( 5, 4) = 0.0d0
        q(5, 4) = -6.0d0
C--------
        d2_b = dble(-4.0)
          g_q( 5, 5) = d2_b * g_y32
        q(5, 5) = dble(-4.0) * y32
C--------
        d2_b = dble(4.0)
          g_q( 5, 6) = d2_b * g_x32
        q(5, 6) = dble(4.0) * x32
C--------
          g_q( 5, 7) = 0.0d0
        q(5, 7) = 6.0d0
C--------
        d2_b = dble(-2.0)
          g_q( 5, 8) = d2_b * g_y32
        q(5, 8) = dble(-2.0) * y32
C--------
        d2_b = dble(2.0)
          g_q( 5, 9) = d2_b * g_x32
        q(5, 9) = dble(2.0) * x32
C--------
          g_q( 6, 4) = 0.0d0
        q(6, 4) = 6.0d0
C--------
        d2_b = dble(2.0)
          g_q( 6, 5) = d2_b * g_y32
        q(6, 5) = dble(2.0) * y32
C--------
        d2_b = dble(-2.0)
          g_q( 6, 6) = d2_b * g_x32
        q(6, 6) = dble(-2.0) * x32
C--------
          g_q( 6, 7) = 0.0d0
        q(6, 7) = -6.0d0
C--------
        d2_b = dble(4.0)
          g_q( 6, 8) = d2_b * g_y32
        q(6, 8) = dble(4.0) * y32
C--------
        d2_b = dble(-4.0)
          g_q( 6, 9) = d2_b * g_x32
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
            g_rm( 1, 1) = d5_b * g_al1
          rm(1, 1) = l3 + al1 * l2 - (1.0d0 + al1) / dble(3.)
C--------
          d5_b = -(1.0d0 / dble(3.)) + l2
            g_rm( 1, 2) = d5_b * g_bl1
          rm(1, 2) = l1 + bl1 * l2 - (1.0d0 + bl1) / dble(3.)
C--------
          d5_b = -(1.0d0 / dble(3.)) + l3
            g_rm( 2, 3) = d5_b * g_al2
          rm(2, 3) = l1 + al2 * l3 - (1.0d0 + al2) / dble(3.)
C--------
          d5_b = -(1.0d0 / dble(3.)) + l3
            g_rm( 2, 4) = d5_b * g_bl2
          rm(2, 4) = l2 + bl2 * l3 - (1.0d0 + bl2) / dble(3.)
C--------
          d5_b = -(1.0d0 / dble(3.)) + l1
            g_rm( 3, 5) = d5_b * g_al3
          rm(3, 5) = l2 + al3 * l1 - (1.0d0 + al3) / dble(3.)
C--------
          d5_b = -(1.0d0 / dble(3.)) + l1
            g_rm( 3, 6) = d5_b * g_bl3
          rm(3, 6) = l3 + bl3 * l1 - (1.0d0 + bl3) / dble(3.)
C--------
          do 99986 j = 1, 6
              g_d1_w = sds(1, 2) * g_rm( 2, j) + rm(2, j) * 
     &    g_sds( 1, 2) + sds(1, 1) * g_rm( 1, j) + rm(1, j) * 
     &    g_sds( 1, 1)
            d1_w = sds(1, 1) * rm(1, j) + sds(1, 2) * rm(2, j)
              g_s1 = sds(1, 3) * g_rm( 3, j) + rm(3, j) * 
     &    g_sds( 1, 3) + g_d1_w
            s1 = d1_w + sds(1, 3) * rm(3, j)
C--------
              g_d1_w = sds(2, 2) * g_rm( 2, j) + rm(2, j) *
     &        g_sds( 2, 2) + sds(2, 1) * g_rm( 1, j) + rm(1, j) * 
     &        g_sds( 2, 1)
            d1_w = sds(2, 1) * rm(1, j) + sds(2, 2) * rm(2, j)
              g_s2 = sds(2, 3) * g_rm( 3, j) + rm(3, j) * 
     &        g_sds( 2, 3) + g_d1_w
            s2 = d1_w + sds(2, 3) * rm(3, j)
C--------
              g_d1_w = sds(3, 2) * g_rm( 2, j) + rm(2, j) *
     &         g_sds( 3, 2) + sds(3, 1) * g_rm( 1, j) + rm(1, j) *
     &         g_sds( 3, 1)
            d1_w = sds(3, 1) * rm(1, j) + sds(3, 2) * rm(2, j)
              g_s3 = sds(3, 3) * g_rm( 3, j) + rm(3, j) * 
     &        g_sds( 3, 3) + g_d1_w
            s3 = d1_w + sds(3, 3) * rm(3, j)
C--------
            do 99987 i = 1, j
                g_d1_w = s2 * g_rm( 2, i) + rm(2, i) * g_s2 + 
     &        s1 * g_rm( 1, i) + rm(1, i) * g_s1
              d1_w = s1 * rm(1, i) + s2 * rm(2, i)
                g_rsd( i, j) = s3 * g_rm( 3, i) + rm(3, i) *
     &         g_s3 + g_d1_w + g_rsd( i, j)
              rsd(i, j) = rsd(i, j) + (d1_w + s3 * rm(3, i))
C--------
                g_rsd( j, i) = g_rsd( i, j)
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
            g_d1_w = rsd(1, 2) * g_q( 2, j) + q(2, j) * g_rsd( 1, 2) + 
     &        rsd(1, 1) * g_q( 1, j) + q(1, j) * g_rsd( 1, 1)
          d1_w = rsd(1, 1) * q(1, j) + rsd(1, 2) * q(2, j)
            g_d2_w = rsd(1, 4) * g_q( 4, j) + q(4, j) * g_rsd( 1, 4) + 
     &        rsd(1, 3) * g_q( 3, j) + q(3, j) * g_rsd( 1, 3) + g_d1_w
          d2_w = d1_w + rsd(1, 3) * q(3, j) + rsd(1, 4) * q(4, j)
            g_s1 = rsd(1, 6) * g_q( 6, j) + q(6, j) * g_rsd( 1, 6) +
     &         rsd(1, 5) * g_q( 5, j) + q(5, j) * g_rsd( 1, 5) + g_d2_w
          s1 = d2_w + rsd(1, 5) * q(5, j) + rsd(1, 6) * q(6, j)
C--------
            g_d1_w = rsd(2, 2) * g_q( 2, j) + q(2, j) * g_rsd( 2, 2) + 
     &        rsd(2, 1) * g_q( 1, j) + q(1, j) * g_rsd( 2, 1)
          d1_w = rsd(2, 1) * q(1, j) + rsd(2, 2) * q(2, j)
            g_d2_w = rsd(2, 4) * g_q( 4, j) + q(4, j) * g_rsd( 2, 4) + 
     &        rsd(2, 3) * g_q( 3, j) + q(3, j) * g_rsd( 2, 3) + g_d1_w
          d2_w = d1_w + rsd(2, 3) * q(3, j) + rsd(2, 4) * q(4, j)
            g_s2 = rsd(2, 6) * g_q( 6, j) + q(6, j) * g_rsd( 2, 6) + 
     &        rsd(2, 5) * g_q( 5, j) + q(5, j) * g_rsd( 2, 5) + g_d2_w
          s2 = d2_w + rsd(2, 5) * q(5, j) + rsd(2, 6) * q(6, j)
C--------
            g_d1_w = rsd(3, 2) * g_q( 2, j) + q(2, j) * g_rsd( 3, 2) +
     &         rsd(3, 1) * g_q( 1, j) + q(1, j) * g_rsd( 3, 1)
          d1_w = rsd(3, 1) * q(1, j) + rsd(3, 2) * q(2, j)
            g_d2_w = rsd(3, 4) * g_q( 4, j) + q(4, j) * g_rsd( 3, 4) + 
     &        rsd(3, 3) * g_q( 3, j) + q(3, j) * g_rsd( 3, 3) + g_d1_w
          d2_w = d1_w + rsd(3, 3) * q(3, j) + rsd(3, 4) * q(4, j)
            g_s3 = rsd(3, 6) * g_q( 6, j) + q(6, j) * g_rsd( 3, 6) + 
     &        rsd(3, 5) * g_q( 5, j) + q(5, j) * g_rsd( 3, 5) + g_d2_w
          s3 = d2_w + rsd(3, 5) * q(5, j) + rsd(3, 6) * q(6, j)
C--------
            g_d1_w = rsd(4, 2) * g_q( 2, j) + q(2, j) * g_rsd( 4, 2) + 
     &        rsd(4, 1) * g_q( 1, j) + q(1, j) * g_rsd( 4, 1)
          d1_w = rsd(4, 1) * q(1, j) + rsd(4, 2) * q(2, j)
            g_d2_w = rsd(4, 4) * g_q( 4, j) + q(4, j) * g_rsd( 4, 4) + 
     &        rsd(4, 3) * g_q( 3, j) + q(3, j) * g_rsd( 4, 3) + g_d1_w
          d2_w = d1_w + rsd(4, 3) * q(3, j) + rsd(4, 4) * q(4, j)
            g_s4 = rsd(4, 6) * g_q( 6, j) + q(6, j) * g_rsd( 4, 6) + 
     &        rsd(4, 5) * g_q( 5, j) + q(5, j) * g_rsd( 4, 5) + g_d2_w
          s4 = d2_w + rsd(4, 5) * q(5, j) + rsd(4, 6) * q(6, j)
C--------
            g_d1_w = rsd(5, 2) * g_q( 2, j) + q(2, j) * g_rsd( 5, 2) + 
     &        rsd(5, 1) * g_q( 1, j) + q(1, j) * g_rsd( 5, 1)
          d1_w = rsd(5, 1) * q(1, j) + rsd(5, 2) * q(2, j)
            g_d2_w = rsd(5, 4) * g_q( 4, j) + q(4, j) * g_rsd( 5, 4) + 
     &        rsd(5, 3) * g_q( 3, j) + q(3, j) * g_rsd( 5, 3) + g_d1_w
          d2_w = d1_w + rsd(5, 3) * q(3, j) + rsd(5, 4) * q(4, j)
            g_s5 = rsd(5, 6) * g_q( 6, j) + q(6, j) * g_rsd( 5, 6) + 
     &        rsd(5, 5) * g_q( 5, j) + q(5, j) * g_rsd( 5, 5) + g_d2_w
          s5 = d2_w + rsd(5, 5) * q(5, j) + rsd(5, 6) * q(6, j)
C--------
            g_d1_w = rsd(6, 2) * g_q( 2, j) + q(2, j) * g_rsd( 6, 2) + 
     &        rsd(6, 1) * g_q( 1, j) + q(1, j) * g_rsd( 6, 1)
          d1_w = rsd(6, 1) * q(1, j) + rsd(6, 2) * q(2, j)
            g_d2_w = rsd(6, 4) * g_q( 4, j) + q(4, j) * g_rsd( 6, 4) +
     &         rsd(6, 3) * g_q( 3, j) + q(3, j) * g_rsd( 6, 3) + g_d1_w
          d2_w = d1_w + rsd(6, 3) * q(3, j) + rsd(6, 4) * q(4, j)
            g_s6 = rsd(6, 6) * g_q( 6, j) + q(6, j) * g_rsd( 6, 6) + 
     &        rsd(6, 5) * g_q( 5, j) + q(5, j) * g_rsd( 6, 5) + g_d2_w
          s6 = d2_w + rsd(6, 5) * q(5, j) + rsd(6, 6) * q(6, j)
C--------
          do 99984 i = 1, j
            l = ls(i)
              g_d1_w = s2 * g_q( 2, i) + q(2, i) * g_s2 + s1 *
     &         g_q( 1, i) + q(1, i) * g_s1
            d1_w = s1 * q(1, i) + s2 * q(2, i)
              g_d2_w = s4 * g_q( 4, i) + q(4, i) * g_s4 + s3 * 
     &        g_q( 3, i) + q(3, i) * g_s3 + g_d1_w
            d2_w = d1_w + s3 * q(3, i) + s4 * q(4, i)
              g_d3_w = s6 * g_q( 6, i) + q(6, i) * g_s6
     &      + s5 * g_q( 5, i) + q(5, i) * g_s5 + g_d2_w
            d3_w = d2_w + s5 * q(5, i) + s6 * q(6, i)
              g_sm( l, k) = g_d3_w + g_sm( l, k)
            sm(l, k) = sm(l, k) + d3_w
C--------
              g_sm( k, l) = g_sm( l, k)
            sm(k, l) = sm(l, k)
C--------
2800        continue
99984     continue
2600      continue
99983   continue
        return
      end
C=end fortran
