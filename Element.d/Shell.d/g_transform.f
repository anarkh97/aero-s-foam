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
      subroutine g_transform(g_p_, xl, g_xl, ldg_xl, yl, g_yl, ldg_yl, z
     *l, g_zl, ldg_zl, xg, yg, zg, str, g_str, ldg_str)
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
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_str, ldg_xl, ldg_yl, ldg_zl
        double precision d2_v, d4_b, d2_w, d1_w, d2_b, g_s(g_pmax_, 6), 
     *g_str(ldg_str, 6), g_l1(g_pmax_), g_xl(ldg_xl, 3), g_l2(g_pmax_)
        double precision g_l3(g_pmax_), g_m1(g_pmax_), g_yl(ldg_yl, 3), 
     *g_m2(g_pmax_), g_m3(g_pmax_), g_n1(g_pmax_), g_zl(ldg_zl, 3), g_n2
     *(g_pmax_), g_n3(g_pmax_), g_t(g_pmax_, 6, 6)
        double precision g_d1_w(g_pmax_), g_d2_w(g_pmax_)
        save g_t, g_d1_w, g_d2_w
        save g_s, g_l1, g_l2, g_l3, g_m1, g_m2, g_m3, g_n1, g_n2, g_n3
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'transform','g_transform.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do g_i_ = 1, g_p_
          g_s(g_i_, 1) = g_str(g_i_, 1)
        enddo
        s(1) = str(1)
C--------
        do g_i_ = 1, g_p_
          g_s(g_i_, 2) = g_str(g_i_, 2)
        enddo
        s(2) = str(2)
C--------
        do g_i_ = 1, g_p_
          g_s(g_i_, 3) = g_str(g_i_, 3)
        enddo
        s(3) = str(3)
C--------
        do g_i_ = 1, g_p_
          g_s(g_i_, 4) = g_str(g_i_, 4)
        enddo
        s(4) = str(4)
C--------
        do g_i_ = 1, g_p_
          g_s(g_i_, 5) = g_str(g_i_, 5)
        enddo
        s(5) = str(5)
C--------
        do g_i_ = 1, g_p_
          g_s(g_i_, 6) = g_str(g_i_, 6)
        enddo
        s(6) = str(6)
C--------
C
C Compute direction cosines
C     
        do g_i_ = 1, g_p_
          g_l1(g_i_) = xg(3) * g_xl(g_i_, 3) + xg(2) * g_xl(g_i_, 2) + x
     *g(1) * g_xl(g_i_, 1)
        enddo
        l1 = xg(1) * xl(1) + xg(2) * xl(2) + xg(3) * xl(3)
C--------
        do g_i_ = 1, g_p_
          g_l2(g_i_) = yg(3) * g_xl(g_i_, 3) + yg(2) * g_xl(g_i_, 2) + y
     *g(1) * g_xl(g_i_, 1)
        enddo
        l2 = yg(1) * xl(1) + yg(2) * xl(2) + yg(3) * xl(3)
C--------
        do g_i_ = 1, g_p_
          g_l3(g_i_) = zg(3) * g_xl(g_i_, 3) + zg(2) * g_xl(g_i_, 2) + z
     *g(1) * g_xl(g_i_, 1)
        enddo
        l3 = zg(1) * xl(1) + zg(2) * xl(2) + zg(3) * xl(3)
C--------
C
        do g_i_ = 1, g_p_
          g_m1(g_i_) = xg(3) * g_yl(g_i_, 3) + xg(2) * g_yl(g_i_, 2) + x
     *g(1) * g_yl(g_i_, 1)
        enddo
        m1 = xg(1) * yl(1) + xg(2) * yl(2) + xg(3) * yl(3)
C--------
        do g_i_ = 1, g_p_
          g_m2(g_i_) = yg(3) * g_yl(g_i_, 3) + yg(2) * g_yl(g_i_, 2) + y
     *g(1) * g_yl(g_i_, 1)
        enddo
        m2 = yg(1) * yl(1) + yg(2) * yl(2) + yg(3) * yl(3)
C--------
        do g_i_ = 1, g_p_
          g_m3(g_i_) = zg(3) * g_yl(g_i_, 3) + zg(2) * g_yl(g_i_, 2) + z
     *g(1) * g_yl(g_i_, 1)
        enddo
        m3 = zg(1) * yl(1) + zg(2) * yl(2) + zg(3) * yl(3)
C--------
C
        do g_i_ = 1, g_p_
          g_n1(g_i_) = xg(3) * g_zl(g_i_, 3) + xg(2) * g_zl(g_i_, 2) + x
     *g(1) * g_zl(g_i_, 1)
        enddo
        n1 = xg(1) * zl(1) + xg(2) * zl(2) + xg(3) * zl(3)
C--------
        do g_i_ = 1, g_p_
          g_n2(g_i_) = yg(3) * g_zl(g_i_, 3) + yg(2) * g_zl(g_i_, 2) + y
     *g(1) * g_zl(g_i_, 1)
        enddo
        n2 = yg(1) * zl(1) + yg(2) * zl(2) + yg(3) * zl(3)
C--------
        do g_i_ = 1, g_p_
          g_n3(g_i_) = zg(3) * g_zl(g_i_, 3) + zg(2) * g_zl(g_i_, 2) + z
     *g(1) * g_zl(g_i_, 1)
        enddo
        n3 = zg(1) * zl(1) + zg(2) * zl(2) + zg(3) * zl(3)
C--------
C
C Construct the 6x6 transformation matrix
C     
        d2_b = l1 + l1
        do g_i_ = 1, g_p_
          g_t(g_i_, 1, 1) = d2_b * g_l1(g_i_)
        enddo
        t(1, 1) = l1 * l1
C--------
        d2_b = m1 + m1
        do g_i_ = 1, g_p_
          g_t(g_i_, 1, 2) = d2_b * g_m1(g_i_)
        enddo
        t(1, 2) = m1 * m1
C--------
        d2_b = n1 + n1
        do g_i_ = 1, g_p_
          g_t(g_i_, 1, 3) = d2_b * g_n1(g_i_)
        enddo
        t(1, 3) = n1 * n1
C--------
        d2_v = dble(2.0) * l1
        d4_b = m1 * dble(2.0)
        do g_i_ = 1, g_p_
          g_t(g_i_, 1, 4) = d2_v * g_m1(g_i_) + d4_b * g_l1(g_i_)
        enddo
        t(1, 4) = d2_v * m1
C--------
        d2_v = dble(2.0) * m1
        d4_b = n1 * dble(2.0)
        do g_i_ = 1, g_p_
          g_t(g_i_, 1, 5) = d2_v * g_n1(g_i_) + d4_b * g_m1(g_i_)
        enddo
        t(1, 5) = d2_v * n1
C--------
        d2_v = dble(2.0) * n1
        d4_b = l1 * dble(2.0)
        do g_i_ = 1, g_p_
          g_t(g_i_, 1, 6) = d2_v * g_l1(g_i_) + d4_b * g_n1(g_i_)
        enddo
        t(1, 6) = d2_v * l1
C--------
C     
        d2_b = l2 + l2
        do g_i_ = 1, g_p_
          g_t(g_i_, 2, 1) = d2_b * g_l2(g_i_)
        enddo
        t(2, 1) = l2 * l2
C--------
        d2_b = m2 + m2
        do g_i_ = 1, g_p_
          g_t(g_i_, 2, 2) = d2_b * g_m2(g_i_)
        enddo
        t(2, 2) = m2 * m2
C--------
        d2_b = n2 + n2
        do g_i_ = 1, g_p_
          g_t(g_i_, 2, 3) = d2_b * g_n2(g_i_)
        enddo
        t(2, 3) = n2 * n2
C--------
        d2_v = dble(2.0) * l2
        d4_b = m2 * dble(2.0)
        do g_i_ = 1, g_p_
          g_t(g_i_, 2, 4) = d2_v * g_m2(g_i_) + d4_b * g_l2(g_i_)
        enddo
        t(2, 4) = d2_v * m2
C--------
        d2_v = dble(2.0) * m2
        d4_b = n2 * dble(2.0)
        do g_i_ = 1, g_p_
          g_t(g_i_, 2, 5) = d2_v * g_n2(g_i_) + d4_b * g_m2(g_i_)
        enddo
        t(2, 5) = d2_v * n2
C--------
        d2_v = dble(2.0) * n2
        d4_b = l2 * dble(2.0)
        do g_i_ = 1, g_p_
          g_t(g_i_, 2, 6) = d2_v * g_l2(g_i_) + d4_b * g_n2(g_i_)
        enddo
        t(2, 6) = d2_v * l2
C--------
C     
        d2_b = l3 + l3
        do g_i_ = 1, g_p_
          g_t(g_i_, 3, 1) = d2_b * g_l3(g_i_)
        enddo
        t(3, 1) = l3 * l3
C--------
        d2_b = m3 + m3
        do g_i_ = 1, g_p_
          g_t(g_i_, 3, 2) = d2_b * g_m3(g_i_)
        enddo
        t(3, 2) = m3 * m3
C--------
        d2_b = n3 + n3
        do g_i_ = 1, g_p_
          g_t(g_i_, 3, 3) = d2_b * g_n3(g_i_)
        enddo
        t(3, 3) = n3 * n3
C--------
        d2_v = dble(2.0) * l3
        d4_b = m3 * dble(2.0)
        do g_i_ = 1, g_p_
          g_t(g_i_, 3, 4) = d2_v * g_m3(g_i_) + d4_b * g_l3(g_i_)
        enddo
        t(3, 4) = d2_v * m3
C--------
        d2_v = dble(2.0) * m3
        d4_b = n3 * dble(2.0)
        do g_i_ = 1, g_p_
          g_t(g_i_, 3, 5) = d2_v * g_n3(g_i_) + d4_b * g_m3(g_i_)
        enddo
        t(3, 5) = d2_v * n3
C--------
        d2_v = dble(2.0) * n3
        d4_b = l3 * dble(2.0)
        do g_i_ = 1, g_p_
          g_t(g_i_, 3, 6) = d2_v * g_l3(g_i_) + d4_b * g_n3(g_i_)
        enddo
        t(3, 6) = d2_v * l3
C--------
C     
        do g_i_ = 1, g_p_
          g_t(g_i_, 4, 1) = l1 * g_l2(g_i_) + l2 * g_l1(g_i_)
        enddo
        t(4, 1) = l1 * l2
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 4, 2) = m1 * g_m2(g_i_) + m2 * g_m1(g_i_)
        enddo
        t(4, 2) = m1 * m2
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 4, 3) = n1 * g_n2(g_i_) + n2 * g_n1(g_i_)
        enddo
        t(4, 3) = n1 * n2
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 4, 4) = l2 * g_m1(g_i_) + m1 * g_l2(g_i_) + l1 * g_m
     *2(g_i_) + m2 * g_l1(g_i_)
        enddo
        t(4, 4) = l1 * m2 + l2 * m1
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 4, 5) = m2 * g_n1(g_i_) + n1 * g_m2(g_i_) + m1 * g_n
     *2(g_i_) + n2 * g_m1(g_i_)
        enddo
        t(4, 5) = m1 * n2 + m2 * n1
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 4, 6) = n2 * g_l1(g_i_) + l1 * g_n2(g_i_) + n1 * g_l
     *2(g_i_) + l2 * g_n1(g_i_)
        enddo
        t(4, 6) = n1 * l2 + n2 * l1
C--------
C     
        do g_i_ = 1, g_p_
          g_t(g_i_, 5, 1) = l2 * g_l3(g_i_) + l3 * g_l2(g_i_)
        enddo
        t(5, 1) = l2 * l3
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 5, 2) = m2 * g_m3(g_i_) + m3 * g_m2(g_i_)
        enddo
        t(5, 2) = m2 * m3
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 5, 3) = n2 * g_n3(g_i_) + n3 * g_n2(g_i_)
        enddo
        t(5, 3) = n2 * n3
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 5, 4) = l3 * g_m2(g_i_) + m2 * g_l3(g_i_) + l2 * g_m
     *3(g_i_) + m3 * g_l2(g_i_)
        enddo
        t(5, 4) = l2 * m3 + l3 * m2
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 5, 5) = m3 * g_n2(g_i_) + n2 * g_m3(g_i_) + m2 * g_n
     *3(g_i_) + n3 * g_m2(g_i_)
        enddo
        t(5, 5) = m2 * n3 + m3 * n2
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 5, 6) = n3 * g_l2(g_i_) + l2 * g_n3(g_i_) + n2 * g_l
     *3(g_i_) + l3 * g_n2(g_i_)
        enddo
        t(5, 6) = n2 * l3 + n3 * l2
C--------
C     
        do g_i_ = 1, g_p_
          g_t(g_i_, 6, 1) = l3 * g_l1(g_i_) + l1 * g_l3(g_i_)
        enddo
        t(6, 1) = l3 * l1
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 6, 2) = m3 * g_m1(g_i_) + m1 * g_m3(g_i_)
        enddo
        t(6, 2) = m3 * m1
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 6, 3) = n3 * g_n1(g_i_) + n1 * g_n3(g_i_)
        enddo
        t(6, 3) = n3 * n1
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 6, 4) = l1 * g_m3(g_i_) + m3 * g_l1(g_i_) + l3 * g_m
     *1(g_i_) + m1 * g_l3(g_i_)
        enddo
        t(6, 4) = l3 * m1 + l1 * m3
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 6, 5) = m1 * g_n3(g_i_) + n3 * g_m1(g_i_) + m3 * g_n
     *1(g_i_) + n1 * g_m3(g_i_)
        enddo
        t(6, 5) = m3 * n1 + m1 * n3
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 6, 6) = n1 * g_l3(g_i_) + l3 * g_n1(g_i_) + n3 * g_l
     *1(g_i_) + l1 * g_n3(g_i_)
        enddo
        t(6, 6) = n3 * l1 + n1 * l3
C--------
C
C Perform the multiplication {str'} = T{str}
C     
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = t(1, 2) * g_s(g_i_, 2) + s(2) * g_t(g_i_, 1, 2)
     * + t(1, 1) * g_s(g_i_, 1) + s(1) * g_t(g_i_, 1, 1)
        enddo
        d1_w = t(1, 1) * s(1) + t(1, 2) * s(2)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = t(1, 4) * g_s(g_i_, 4) + s(4) * g_t(g_i_, 1, 4)
     * + t(1, 3) * g_s(g_i_, 3) + s(3) * g_t(g_i_, 1, 3) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + t(1, 3) * s(3) + t(1, 4) * s(4)
        do g_i_ = 1, g_p_
          g_str(g_i_, 1) = t(1, 6) * g_s(g_i_, 6) + s(6) * g_t(g_i_, 1, 
     *6) + t(1, 5) * g_s(g_i_, 5) + s(5) * g_t(g_i_, 1, 5) + g_d2_w(g_i_
     *)
        enddo
        str(1) = d2_w + t(1, 5) * s(5) + t(1, 6) * s(6)
C--------
C
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = t(2, 2) * g_s(g_i_, 2) + s(2) * g_t(g_i_, 2, 2)
     * + t(2, 1) * g_s(g_i_, 1) + s(1) * g_t(g_i_, 2, 1)
        enddo
        d1_w = t(2, 1) * s(1) + t(2, 2) * s(2)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = t(2, 4) * g_s(g_i_, 4) + s(4) * g_t(g_i_, 2, 4)
     * + t(2, 3) * g_s(g_i_, 3) + s(3) * g_t(g_i_, 2, 3) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + t(2, 3) * s(3) + t(2, 4) * s(4)
        do g_i_ = 1, g_p_
          g_str(g_i_, 2) = t(2, 6) * g_s(g_i_, 6) + s(6) * g_t(g_i_, 2, 
     *6) + t(2, 5) * g_s(g_i_, 5) + s(5) * g_t(g_i_, 2, 5) + g_d2_w(g_i_
     *)
        enddo
        str(2) = d2_w + t(2, 5) * s(5) + t(2, 6) * s(6)
C--------
C
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = t(3, 2) * g_s(g_i_, 2) + s(2) * g_t(g_i_, 3, 2)
     * + t(3, 1) * g_s(g_i_, 1) + s(1) * g_t(g_i_, 3, 1)
        enddo
        d1_w = t(3, 1) * s(1) + t(3, 2) * s(2)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = t(3, 4) * g_s(g_i_, 4) + s(4) * g_t(g_i_, 3, 4)
     * + t(3, 3) * g_s(g_i_, 3) + s(3) * g_t(g_i_, 3, 3) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + t(3, 3) * s(3) + t(3, 4) * s(4)
        do g_i_ = 1, g_p_
          g_str(g_i_, 3) = t(3, 6) * g_s(g_i_, 6) + s(6) * g_t(g_i_, 3, 
     *6) + t(3, 5) * g_s(g_i_, 5) + s(5) * g_t(g_i_, 3, 5) + g_d2_w(g_i_
     *)
        enddo
        str(3) = d2_w + t(3, 5) * s(5) + t(3, 6) * s(6)
C--------
C
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = t(4, 2) * g_s(g_i_, 2) + s(2) * g_t(g_i_, 4, 2)
     * + t(4, 1) * g_s(g_i_, 1) + s(1) * g_t(g_i_, 4, 1)
        enddo
        d1_w = t(4, 1) * s(1) + t(4, 2) * s(2)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = t(4, 4) * g_s(g_i_, 4) + s(4) * g_t(g_i_, 4, 4)
     * + t(4, 3) * g_s(g_i_, 3) + s(3) * g_t(g_i_, 4, 3) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + t(4, 3) * s(3) + t(4, 4) * s(4)
        do g_i_ = 1, g_p_
          g_str(g_i_, 4) = t(4, 6) * g_s(g_i_, 6) + s(6) * g_t(g_i_, 4, 
     *6) + t(4, 5) * g_s(g_i_, 5) + s(5) * g_t(g_i_, 4, 5) + g_d2_w(g_i_
     *)
        enddo
        str(4) = d2_w + t(4, 5) * s(5) + t(4, 6) * s(6)
C--------
C
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = t(5, 2) * g_s(g_i_, 2) + s(2) * g_t(g_i_, 5, 2)
     * + t(5, 1) * g_s(g_i_, 1) + s(1) * g_t(g_i_, 5, 1)
        enddo
        d1_w = t(5, 1) * s(1) + t(5, 2) * s(2)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = t(5, 4) * g_s(g_i_, 4) + s(4) * g_t(g_i_, 5, 4)
     * + t(5, 3) * g_s(g_i_, 3) + s(3) * g_t(g_i_, 5, 3) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + t(5, 3) * s(3) + t(5, 4) * s(4)
        do g_i_ = 1, g_p_
          g_str(g_i_, 5) = t(5, 6) * g_s(g_i_, 6) + s(6) * g_t(g_i_, 5, 
     *6) + t(5, 5) * g_s(g_i_, 5) + s(5) * g_t(g_i_, 5, 5) + g_d2_w(g_i_
     *)
        enddo
        str(5) = d2_w + t(5, 5) * s(5) + t(5, 6) * s(6)
C--------
C
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = t(6, 2) * g_s(g_i_, 2) + s(2) * g_t(g_i_, 6, 2)
     * + t(6, 1) * g_s(g_i_, 1) + s(1) * g_t(g_i_, 6, 1)
        enddo
        d1_w = t(6, 1) * s(1) + t(6, 2) * s(2)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = t(6, 4) * g_s(g_i_, 4) + s(4) * g_t(g_i_, 6, 4)
     * + t(6, 3) * g_s(g_i_, 3) + s(3) * g_t(g_i_, 6, 3) + g_d1_w(g_i_)
        enddo
        d2_w = d1_w + t(6, 3) * s(3) + t(6, 4) * s(4)
        do g_i_ = 1, g_p_
          g_str(g_i_, 6) = t(6, 6) * g_s(g_i_, 6) + s(6) * g_t(g_i_, 6, 
     *6) + t(6, 5) * g_s(g_i_, 5) + s(5) * g_t(g_i_, 6, 5) + g_d2_w(g_i_
     *)
        enddo
        str(6) = d2_w + t(6, 5) * s(5) + t(6, 6) * s(6)
C--------
C
C
        return
      end
