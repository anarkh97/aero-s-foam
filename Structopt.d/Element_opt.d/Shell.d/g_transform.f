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
      subroutine g_transform(xl, g_xl, yl, g_yl, zl
     *, g_zl, xg, yg, zg, str, g_str)
C
C**********************************************************************C
C
C Purpose: to form a transformation matrix from local coordinates to
C          global coordinates
C
C input variables:
C      xl = x local unit vector
C      yl = y local unit vector
C      zl = z local unit vector
C      xg = x global unit vector
C      yg = y global unit vector
C      zg = z global unit vector
C      str = stress/strain 6x1 vector
C            sigmaxx, sigmayy, sigmazz, sigma12, sigma23, sigma13
C
C local variables:
C      l1 = direction cosine between xl and xg
C      l2 = direction cosine between xl and yg
C      l3 = direction cosine between xl and zg
C      m1 = direction cosine between yl and xg
C      m2 = direction cosine between yl and yg
C      m3 = direction cosine between yl and zg
C      n1 = direction cosine between zl and xg
C      n2 = direction cosine between zl and yg
C      n3 = direction cosine between zl and zg
C      t  = transformation matrix from local to global
C
C**********************************************************************C
C
C Declarations
C
        double precision xl(3), yl(3), zl(3)
        double precision xg(3), yg(3), zg(3)
        double precision str(6)
C
C Local Declarations
C
        double precision l1, l2, l3
        double precision m1, m2, m3
        double precision n1, n2, n3
        double precision t(6, 6), s(6)
C
C Copy stress/strain values to a temporary array
C
        double precision d2_v, d4_b, d2_w, d1_w, d2_b, g_s(6), 
     *g_str(6), g_l1, g_xl(3), g_l2
        double precision g_l3, g_m1, g_yl(3), 
     *g_m2, g_m3, g_n1, g_zl(3), g_n2,
     * g_n3, g_t(6, 6)
        double precision g_d1_w, g_d2_w
        save g_t, g_d1_w, g_d2_w
        save g_s, g_l1, g_l2, g_l3, g_m1, g_m2, g_m3, g_n1, g_n2, g_n3
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'transform','g_transform.f')
C
          g_s(1) = g_str(1)
        s(1) = str(1)
C--------
          g_s(2) = g_str(2)
        s(2) = str(2)
C--------
          g_s(3) = g_str(3)
        s(3) = str(3)
C--------
          g_s(4) = g_str(4)
        s(4) = str(4)
C--------
          g_s(5) = g_str(5)
        s(5) = str(5)
C--------       
          g_s(6) = g_str(6)
        s(6) = str(6)
C--------
C
C Compute direction cosines
C     
          g_l1 = xg(3) * g_xl(3) + xg(2) * g_xl(2) + 
     *xg(1) * g_xl(1)
        l1 = xg(1) * xl(1) + xg(2) * xl(2) + xg(3) * xl(3)
C--------
          g_l2 = yg(3) * g_xl(3) + yg(2) * g_xl(2) + 
     *yg(1) * g_xl(1)
        l2 = yg(1) * xl(1) + yg(2) * xl(2) + yg(3) * xl(3)
C--------
          g_l3 = zg(3) * g_xl(3) + zg(2) * g_xl(2) + 
     *zg(1) * g_xl(1)        
        l3 = zg(1) * xl(1) + zg(2) * xl(2) + zg(3) * xl(3)
C--------
C
          g_m1 = xg(3) * g_yl(3) + xg(2) * g_yl(2) + 
     *xg(1) * g_yl(1)
        m1 = xg(1) * yl(1) + xg(2) * yl(2) + xg(3) * yl(3)
C--------
          g_m2 = yg(3) * g_yl(3) + yg(2) * g_yl(2) + 
     *yg(1) * g_yl(1)
        m2 = yg(1) * yl(1) + yg(2) * yl(2) + yg(3) * yl(3)
C--------
          g_m3 = zg(3) * g_yl(3) + zg(2) * g_yl(2) + 
     *zg(1) * g_yl(1)
        m3 = zg(1) * yl(1) + zg(2) * yl(2) + zg(3) * yl(3)
C--------
C
          g_n1 = xg(3) * g_zl(3) + xg(2) * g_zl(2) + 
     *xg(1) * g_zl(1)
        n1 = xg(1) * zl(1) + xg(2) * zl(2) + xg(3) * zl(3)
C--------        do 1 = 1, 1
          g_n2 = yg(3) * g_zl(3) + yg(2) * g_zl(2) + 
     *yg(1) * g_zl(1)
        n2 = yg(1) * zl(1) + yg(2) * zl(2) + yg(3) * zl(3)
C--------
          g_n3 = zg(3) * g_zl(3) + zg(2) * g_zl(2) + 
     *zg(1) * g_zl(1)
        n3 = zg(1) * zl(1) + zg(2) * zl(2) + zg(3) * zl(3)
C--------
C
C Construct the 6x6 transformation matrix
C     
        d2_b = l1 + l1
          g_t(1, 1) = d2_b * g_l1
        t(1, 1) = l1 * l1
C--------
        d2_b = m1 + m1
          g_t(1, 2) = d2_b * g_m1
        t(1, 2) = m1 * m1
C--------
        d2_b = n1 + n1
          g_t(1, 3) = d2_b * g_n1
        t(1, 3) = n1 * n1
C--------
        d2_v = dble(2.0) * l1
        d4_b = m1 * dble(2.0)
          g_t(1, 4) = d2_v * g_m1 + d4_b * g_l1
        t(1, 4) = d2_v * m1
C--------
        d2_v = dble(2.0) * m1
        d4_b = n1 * dble(2.0)
          g_t(1, 5) = d2_v * g_n1 + d4_b * g_m1
        t(1, 5) = d2_v * n1
C--------
        d2_v = dble(2.0) * n1
        d4_b = l1 * dble(2.0)
          g_t(1, 6) = d2_v * g_l1 + d4_b * g_n1
        t(1, 6) = d2_v * l1
C--------
C     
        d2_b = l2 + l2
          g_t(2, 1) = d2_b * g_l2
        t(2, 1) = l2 * l2
C--------
        d2_b = m2 + m2
          g_t(2, 2) = d2_b * g_m2
        t(2, 2) = m2 * m2
C--------
        d2_b = n2 + n2
          g_t(2, 3) = d2_b * g_n2
        t(2, 3) = n2 * n2
C--------
        d2_v = dble(2.0) * l2
        d4_b = m2 * dble(2.0)
          g_t(2, 4) = d2_v * g_m2 + d4_b * g_l2
        t(2, 4) = d2_v * m2
C--------
        d2_v = dble(2.0) * m2
        d4_b = n2 * dble(2.0)
          g_t(2, 5) = d2_v * g_n2 + d4_b * g_m2
        t(2, 5) = d2_v * n2
C--------
        d2_v = dble(2.0) * n2
        d4_b = l2 * dble(2.0)
          g_t(2, 6) = d2_v * g_l2 + d4_b * g_n2
        t(2, 6) = d2_v * l2
C--------
C     
        d2_b = l3 + l3
          g_t(3, 1) = d2_b * g_l3
        t(3, 1) = l3 * l3
C--------
        d2_b = m3 + m3
          g_t(3, 2) = d2_b * g_m3
        t(3, 2) = m3 * m3
C--------
        d2_b = n3 + n3
          g_t(3, 3) = d2_b * g_n3
        t(3, 3) = n3 * n3
C--------
        d2_v = dble(2.0) * l3
        d4_b = m3 * dble(2.0)
          g_t(3, 4) = d2_v * g_m3 + d4_b * g_l3
        t(3, 4) = d2_v * m3
C--------
        d2_v = dble(2.0) * m3
        d4_b = n3 * dble(2.0)
          g_t(3, 5) = d2_v * g_n3 + d4_b * g_m3
        t(3, 5) = d2_v * n3
C--------
        d2_v = dble(2.0) * n3
        d4_b = l3 * dble(2.0)
          g_t(3, 6) = d2_v * g_l3 + d4_b * g_n3  
        t(3, 6) = d2_v * l3
C--------
C     
          g_t(4, 1) = l1 * g_l2 + l2 * g_l1
        t(4, 1) = l1 * l2
C--------
          g_t(4, 2) = m1 * g_m2 + m2 * g_m1
        t(4, 2) = m1 * m2
C--------
          g_t(4, 3) = n1 * g_n2 + n2 * g_n1
        t(4, 3) = n1 * n2
C--------
          g_t(4, 4) = l2 * g_m1 + m1 * g_l2 + l1 * 
     *g_m2 + m2 * g_l1
        t(4, 4) = l1 * m2 + l2 * m1
C--------
          g_t(4, 5) = m2 * g_n1 + n1 * g_m2 + m1 * 
     *g_n2 + n2 * g_m1  
        t(4, 5) = m1 * n2 + m2 * n1
C--------
          g_t(4, 6) = n2 * g_l1 + l1 * g_n2 + n1 * 
     *g_l2 + l2 * g_n1  
        t(4, 6) = n1 * l2 + n2 * l1
C--------
C     
          g_t(5, 1) = l2 * g_l3 + l3 * g_l2
        t(5, 1) = l2 * l3
C--------
          g_t(5, 2) = m2 * g_m3 + m3 * g_m2
        t(5, 2) = m2 * m3
C--------
          g_t(5, 3) = n2 * g_n3 + n3 * g_n2
        t(5, 3) = n2 * n3
C--------
          g_t(5, 4) = l3 * g_m2 + m2 * g_l3 + l2 * 
     *g_m3 + m3 * g_l2
        t(5, 4) = l2 * m3 + l3 * m2
C--------
          g_t(5, 5) = m3 * g_n2 + n2 * g_m3 + m2 * 
     *g_n3 + n3 * g_m2
        t(5, 5) = m2 * n3 + m3 * n2
C--------
          g_t(5, 6) = n3 * g_l2 + l2 * g_n3 + n2 * 
     *g_l3 + l3 * g_n2
        t(5, 6) = n2 * l3 + n3 * l2
C--------
C     
          g_t(6, 1) = l3 * g_l1 + l1 * g_l3
        t(6, 1) = l3 * l1
C--------
          g_t(6, 2) = m3 * g_m1 + m1 * g_m3
        t(6, 2) = m3 * m1
C--------
          g_t(6, 3) = n3 * g_n1 + n1 * g_n3
        t(6, 3) = n3 * n1
C--------
          g_t(6, 4) = l1 * g_m3 + m3 * g_l1 + l3 * 
     *g_m1 + m1 * g_l3
        t(6, 4) = l3 * m1 + l1 * m3
C--------
          g_t(6, 5) = m1 * g_n3 + n3 * g_m1 + m3 * 
     *g_n1 + n1 * g_m3
        t(6, 5) = m3 * n1 + m1 * n3
C--------
          g_t(6, 6) = n1 * g_l3 + l3 * g_n1 + n3 * 
     *g_l1 + l1 * g_n3
        t(6, 6) = n3 * l1 + n1 * l3
C--------
C
C Perform the multiplication {str'} = T{str}
C     
          g_d1_w = t(1, 2) * g_s(2) + s(2) * g_t(1, 2)
     * + t(1, 1) * g_s(1) + s(1) * g_t(1, 1)  
        d1_w = t(1, 1) * s(1) + t(1, 2) * s(2)
          g_d2_w = t(1, 4) * g_s(4) + s(4) * g_t(1, 4)
     * + t(1, 3) * g_s(3) + s(3) * g_t(1, 3) + g_d1_w
        d2_w = d1_w + t(1, 3) * s(3) + t(1, 4) * s(4)
          g_str(1) = t(1, 6) * g_s(6) + s(6) * g_t(1, 
     *6) + t(1, 5) * g_s(5) + s(5) * g_t(1, 5) + g_d2_w
        str(1) = d2_w + t(1, 5) * s(5) + t(1, 6) * s(6)
C--------
C
          g_d1_w = t(2, 2) * g_s(2) + s(2) * g_t(2, 2)
     * + t(2, 1) * g_s(1) + s(1) * g_t(2, 1)
        d1_w = t(2, 1) * s(1) + t(2, 2) * s(2)
          g_d2_w = t(2, 4) * g_s(4) + s(4) * g_t(2, 4)
     * + t(2, 3) * g_s(3) + s(3) * g_t(2, 3) + g_d1_w 
        d2_w = d1_w + t(2, 3) * s(3) + t(2, 4) * s(4)
          g_str(2) = t(2, 6) * g_s(6) + s(6) * g_t(2, 
     *6) + t(2, 5) * g_s(5) + s(5) * g_t(2, 5) + g_d2_w
        str(2) = d2_w + t(2, 5) * s(5) + t(2, 6) * s(6)
C--------
C
          g_d1_w = t(3, 2) * g_s(2) + s(2) * g_t(3, 2)
     * + t(3, 1) * g_s(1) + s(1) * g_t(3, 1)
        d1_w = t(3, 1) * s(1) + t(3, 2) * s(2)
          g_d2_w = t(3, 4) * g_s(4) + s(4) * g_t(3, 4)
     * + t(3, 3) * g_s(3) + s(3) * g_t(3, 3) + g_d1_w
        d2_w = d1_w + t(3, 3) * s(3) + t(3, 4) * s(4)
          g_str(3) = t(3, 6) * g_s(6) + s(6) * g_t(3, 
     *6) + t(3, 5) * g_s(5) + s(5) * g_t(3, 5) + g_d2_w
        str(3) = d2_w + t(3, 5) * s(5) + t(3, 6) * s(6)
C--------
C
          g_d1_w = t(4, 2) * g_s(2) + s(2) * g_t(4, 2)
     * + t(4, 1) * g_s(1) + s(1) * g_t(4, 1)
        d1_w = t(4, 1) * s(1) + t(4, 2) * s(2)
          g_d2_w = t(4, 4) * g_s(4) + s(4) * g_t(4, 4)
     * + t(4, 3) * g_s(3) + s(3) * g_t(4, 3) + g_d1_w
        d2_w = d1_w + t(4, 3) * s(3) + t(4, 4) * s(4)
          g_str(4) = t(4, 6) * g_s(6) + s(6) * g_t(4, 
     *6) + t(4, 5) * g_s(5) + s(5) * g_t(4, 5) + g_d2_w
        str(4) = d2_w + t(4, 5) * s(5) + t(4, 6) * s(6)
C--------
C
          g_d1_w = t(5, 2) * g_s(2) + s(2) * g_t(5, 2)
     * + t(5, 1) * g_s(1) + s(1) * g_t(5, 1)
        d1_w = t(5, 1) * s(1) + t(5, 2) * s(2)
          g_d2_w = t(5, 4) * g_s(4) + s(4) * g_t(5, 4)
     * + t(5, 3) * g_s(3) + s(3) * g_t(5, 3) + g_d1_w
        d2_w = d1_w + t(5, 3) * s(3) + t(5, 4) * s(4)
          g_str(5) = t(5, 6) * g_s(6) + s(6) * g_t(5, 
     *6) + t(5, 5) * g_s(5) + s(5) * g_t(5, 5) + g_d2_w
        str(5) = d2_w + t(5, 5) * s(5) + t(5, 6) * s(6)
C--------
C
          g_d1_w = t(6, 2) * g_s(2) + s(2) * g_t(6, 2)
     * + t(6, 1) * g_s(1) + s(1) * g_t(6, 1)
        d1_w = t(6, 1) * s(1) + t(6, 2) * s(2)
          g_d2_w = t(6, 4) * g_s(4) + s(4) * g_t(6, 4)
     * + t(6, 3) * g_s(3) + s(3) * g_t(6, 3) + g_d1_w 
        d2_w = d1_w + t(6, 3) * s(3) + t(6, 4) * s(4)
          g_str(6) = t(6, 6) * g_s(6) + s(6) * g_t(6, 
     *6) + t(6, 5) * g_s(5) + s(5) * g_t(6, 5) + g_d2_w
        str(6) = d2_w + t(6, 5) * s(5) + t(6, 6) * s(6)
C--------
C
C
        return
      end
