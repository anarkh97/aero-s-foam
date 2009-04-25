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
C=BLOCK FORTRAN
      subroutine g_membra(x, g_x, y, g_y, alpha, le,
     * q, g_q, status)
C
C                   T Y P E   &   D I M E N S I O N
C
C=BLOCK VAX
C     implicit      none
C=END VAX
        character*(*) status
C=BLOCK REAL*8
        real*8 x(3), y(3), p(9, 3), q(18, 3)
        real*8 area2, c, alpha
        real*8 x21, x32, x13, y21, y32, y13
        real*8 x12, x23, x31, y12, y23, y31
C=ELSE
C=END REAL*8
        integer i, j, kk
        integer le(9), n
C
C                   L O G I C
C
        double precision d9_b, d8_b, d4_v, d2_b, d7_b, d6_b, d5_b, d4_b,
     * d2_v, d3_b
        double precision g_x21, g_x(3), g_x12, 
     *g_x32, g_x23, g_x13, g_x31, 
     *g_y21, g_y(3), g_y12
        double precision g_y32, g_y23, g_y13,
     * g_y31, g_area2, g_p(9, 3), 
     *g_c, g_q(18, 3)
        save g_y13, g_y31, g_area2, g_p, g_c
        save g_x21, g_x12, g_x32, g_x23, g_x13, g_x31, g_y21, g_y12, 
     *g_y32, g_y23
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'membra','g_membra.f')
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
          g_area2 = -x21 * g_y13 + (-y13) * g_x21 + 
     *y21 * g_x13 + x13 * g_y21
        area2 = y21 * x13 - x21 * y13
C--------
        if (area2 .le. 0.0) then
          status = 'SM3MB: Zero area'
          if (area2 .eq. 0.0) then
            status = 'SM3MB: Zero area'
          endif
          return
        endif
          g_p(1, 1) = g_y23
        p(1, 1) = y23
C--------
          g_p(2, 1) = 0.0d0
        p(2, 1) = 0.0d0
C--------
          g_p(3, 1) = g_y31
        p(3, 1) = y31
C--------
          g_p(4, 1) = 0.0d0
        p(4, 1) = 0.0d0
C--------
          g_p(5, 1) = g_y12
        p(5, 1) = y12
C--------
          g_p(6, 1) = 0.0d0
        p(6, 1) = 0.0d0
C--------
          g_p(1, 2) = 0.0d0     
        p(1, 2) = 0.0d0
C--------
          g_p(2, 2) = g_x32
        p(2, 2) = x32
C--------
          g_p(3, 2) = 0.0d0
        p(3, 2) = 0.0d0
C--------
          g_p(4, 2) = g_x13
        p(4, 2) = x13
C--------
          g_p(5, 2) = 0.0d0
        p(5, 2) = 0.0d0
C--------
          g_p(6, 2) = g_x21
        p(6, 2) = x21
C--------
          g_p(1, 3) = g_x32
        p(1, 3) = x32
C--------
          g_p(2, 3) = g_y23
        p(2, 3) = y23
C--------
          g_p(3, 3) = g_x13
        p(3, 3) = x13
C--------
          g_p(4, 3) = g_y31
        p(4, 3) = y31
C--------
          g_p(5, 3) = g_x21
        p(5, 3) = x21
C--------
          g_p(6, 3) = g_y12     
        p(6, 3) = y12
C--------
        n = 6
        if (alpha .ne. 0.0) then
          d4_v = y13 - y21
          d3_b = 1.0d0 / dble(6.) * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * y23
            g_p(7, 1) = -d5_b * g_y21 + d5_b * g_y13 +
     * d4_b * g_y23
          p(7, 1) = y23 * d4_v * alpha / dble(6.)
C--------
          d4_v = x31 - x12
          d3_b = 1.0d0 / dble(6.) * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x32
            g_p(7, 2) = -d5_b * g_x12 + d5_b * g_x31 +
     * d4_b * g_x32
          p(7, 2) = x32 * d4_v * alpha / dble(6.)
C--------
          d3_b = 1.0d0 / dble(3.) * alpha
          d6_b = -d3_b * y21
          d7_b = -d3_b * x12
          d8_b = d3_b * y13
          d9_b = d3_b * x31
            g_p(7, 3) = d7_b * g_y21 + d6_b * g_x12 + 
     *d9_b * g_y13 + d8_b * g_x31
          p(7, 3) = (x31 * y13 - x12 * y21) * alpha / dble(3.)
C--------
          d4_v = y21 - y32
          d3_b = 1.0d0 / dble(6.) * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * y31
            g_p(8, 1) = -d5_b * g_y32 + d5_b * g_y21 +
     * d4_b * g_y31
          p(8, 1) = y31 * d4_v * alpha / dble(6.)
C--------
          d4_v = x12 - x23
          d3_b = 1.0d0 / dble(6.) * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x13
            g_p(8, 2) = -d5_b * g_x23 + d5_b * g_x12 +
     * d4_b * g_x13
          p(8, 2) = x13 * d4_v * alpha / dble(6.)
C--------
          d3_b = 1.0d0 / dble(3.) * alpha
          d6_b = -d3_b * y32
          d7_b = -d3_b * x23
          d8_b = d3_b * y21
          d9_b = d3_b * x12
            g_p(8, 3) = d7_b * g_y32 + d6_b * g_x23 + 
     *d9_b * g_y21 + d8_b * g_x12
          p(8, 3) = (x12 * y21 - x23 * y32) * alpha / dble(3.)
C--------
          d4_v = y32 - y13
          d3_b = 1.0d0 / dble(6.) * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * y12
            g_p(9, 1) = -d5_b * g_y13 + d5_b * g_y32 +
     * d4_b * g_y12
          p(9, 1) = y12 * d4_v * alpha / dble(6.)
C--------
          d4_v = x23 - x31
          d3_b = 1.0d0 / dble(6.) * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x21
            g_p(9, 2) = -d5_b * g_x31 + d5_b * g_x23 +
     * d4_b * g_x21
          p(9, 2) = x21 * d4_v * alpha / dble(6.)
C--------
          d3_b = 1.0d0 / dble(3.) * alpha
          d6_b = -d3_b * y13
          d7_b = -d3_b * x31
          d8_b = d3_b * y32
          d9_b = d3_b * x23
            g_p(9, 3) = d7_b * g_y13 + d6_b * g_x31 + 
     *d9_b * g_y32 + d8_b * g_x23
          p(9, 3) = (x23 * y32 - x31 * y13) * alpha / dble(3.)
C--------
          n = 9
        endif
        d2_v = 1.0d0 / area2
        d2_b = -d2_v / area2
          g_c = d2_b * g_area2
        c = d2_v
C--------
        do 99998 i = 1, n
          kk = le(i)
          do 99999 j = 1, 3
              g_q(kk, j) = p(i, j) * g_c + c *  
     *g_p(i,j)
            q(kk, j) = p(i, j) * c
C--------
300         continue
99999     continue
99998   continue
        return
      end
C=END FORTRAN
