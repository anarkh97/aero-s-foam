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
      subroutine gxtria3d(flag,xl,g_xl,yl,g_yl,zl,g_zl,e,g_e,nu,
     *                    h,g_h,rk,g_rk)
C
C-------------------------------------------------------------------*
C This subroutine evaluates the stiffness matrix and mass vector
C for the spatial 18 d.o.f tree node triangle obtained as
C a combination of the aqr bending triangle plus the membrane
C with drilling d.o.f. developed by Felippa et al.
C
C        rkb will be used for the basic stiffness 
C      rkm will be used for the higher order stiffness
C
C input variables:
C     xl = x coordinates
C     yl = y coordinates
C     zl = z coordinates
C     e  = elastic modulus
C     nu = poisson's ratio
C     h  = thickness
C
C output variables:
C     rk = stiffness matrix
C 
C-------------------------------------------------------------------*
C
C  SUBROUTINES CALLED: BASICO
C                      SM3MB
C                      SMCBH
C                      SM3MHEFF
C                      ROTATION
C                      TRIROT
C-------------------------------------------------------------------*
C
        real*8 e, nu, alpha
        real*8 xl(3), yl(3), zl(3), h(3), rk(18, 18)
        real*8 db(3, 3), dm(3, 3)
        real*8 xp(3), yp(3), zp(3), xlp(3), ylp(3), zlp(3)
        real*8 r1(3, 3), rkm(18, 18), rkb(18, 18)
        real*8 v1n(3), v2n(3), v3n(3)
        real*8 cb, x21, y21, z21, x32, y32, z32, x13, y13, z13
        real*8 rx, ry, rz, bx, by, bz, rlr, rlb, bpr, area, ylr, zlr
        real*8 ycg, xcg, zcg, xlcg, ylcg, zlcg, f, clr, cqr, esp
        integer lb(9), le(9)
        integer i, j, flag
        character*10 status
C
        parameter (f = 1.0d+00, clr = 0.0d+00, cqr = 1.0d+00)
        parameter (alpha = 1.5d+00)
C
        double precision d1
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_h, ldg_e, ldg_xl, ldg_yl, ldg_zl, ldg_rk
c
c manually inserted - begin	     
c

      parameter (g_p_=1,ldg_xl=1,ldg_yl=1,ldg_zl=1,ldg_e=1,
     *           ldg_h=1,ldg_rk=1)
c
c manually inserted - end      
c
        double precision d4_b, d5_b, d2_v, d3_v, d2_w, d1_w, d1_p, d6_b,
     * d2_b, d3_b
        double precision g_esp(g_pmax_), g_h(ldg_h, 3), g_cb(g_pmax_), g
     *_e(ldg_e), g_db(g_pmax_, 3, 3), g_dm(g_pmax_, 3, 3), g_x21(g_pmax_
     *), g_xl(ldg_xl, 3), g_y21(g_pmax_), g_yl(ldg_yl, 3)
        double precision g_z21(g_pmax_), g_zl(ldg_zl, 3), g_x32(g_pmax_)
     *, g_y32(g_pmax_), g_z32(g_pmax_), g_d2_w(g_pmax_), g_d1_w(g_pmax_)
     *, g_rlr(g_pmax_), g_xp(g_pmax_, 3), g_zp(g_pmax_, 3)
        double precision g_zlr(g_pmax_), g_yp(g_pmax_, 3), g_ylr(g_pmax_
     *), g_xcg(g_pmax_), g_ycg(g_pmax_), g_zcg(g_pmax_), g_xlcg(g_pmax_)
     *, g_ylcg(g_pmax_), g_zlcg(g_pmax_), g_xlp(g_pmax_, 3)
        double precision g_ylp(g_pmax_, 3), g_rk(ldg_rk, 18, 18), g_rkm(
     *g_pmax_, 18, 18), g_rkb(g_pmax_, 18, 18), g_r1(g_pmax_, 3, 3)
        save g_zcg, g_xlcg, g_ylcg, g_zlcg, g_xlp, g_ylp, g_rkm, g_rkb, 
     *g_r1
        save g_d2_w, g_d1_w, g_rlr, g_xp, g_zp, g_zlr, g_yp, g_ylr, g_xc
     *g, g_ycg
        save g_esp, g_cb, g_db, g_dm, g_x21, g_y21, g_z21, g_x32, g_y32,
     * g_z32
        external g_trirotation
        external g_rotation
        external g_sm3mhe
        external g_smcbh
        external g_sm3mb
        external g_basico
        intrinsic dble
        data lb /3, 4, 5, 9, 10, 11, 15, 16, 17/
        data le /1, 2, 7, 8, 13, 14, 6, 12, 18/
        data v1n /1.0, 0.0, 0.0/
        data v2n /0.0, 1.0, 0.0/
        data v3n /0.0, 0.0, 1.0/
C
C
C h   = shell thickness
C e   = elastic modulus
C nu  = poisson's ratio
C
C flag = 0 do NOT perform transform from local to global
C flag = 1 perform transformation from local to global
C
C set thickness
C
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'tria3d','g_tria3d.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do g_i_ = 1, g_p_
          g_esp(g_i_) = g_h(g_i_, 1)
        enddo
        esp = h(1)
C--------
C
C membrane elastic matrix
C
        d3_v = esp ** ( 3 - 2)
        d3_v =  d3_v * esp
        d1_p =  3 *  d3_v
        d3_v =  d3_v * esp
        d3_b = 1.0d0 / (1.0d0 - nu * nu) * (1.0d0 / dble(12.0))
        d4_b = d3_b * d3_v
        d6_b = d3_b * e * d1_p
        do g_i_ = 1, g_p_
          g_cb(g_i_) = d6_b * g_esp(g_i_) + d4_b * g_e(g_i_)
        enddo
        cb = e * d3_v / dble(12.0) / (1.0d0 - nu * nu)
C--------
        do g_i_ = 1, g_p_
          g_db(g_i_, 1, 1) = g_cb(g_i_)
        enddo
        db(1, 1) = cb
C--------
        do g_i_ = 1, g_p_
          g_db(g_i_, 1, 2) = nu * g_cb(g_i_)
        enddo
        db(1, 2) = nu * cb
C--------
        do g_i_ = 1, g_p_
          g_db(g_i_, 1, 3) = 0.0d0
        enddo
        db(1, 3) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_db(g_i_, 2, 1) = g_db(g_i_, 1, 2)
        enddo
        db(2, 1) = db(1, 2)
C--------
        do g_i_ = 1, g_p_
          g_db(g_i_, 2, 2) = g_cb(g_i_)
        enddo
        db(2, 2) = cb
C--------
        do g_i_ = 1, g_p_
          g_db(g_i_, 2, 3) = 0.0d0
        enddo
        db(2, 3) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_db(g_i_, 3, 1) = 0.0d0
        enddo
        db(3, 1) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_db(g_i_, 3, 2) = 0.0d0
        enddo
        db(3, 2) = 0.0d0
C--------
        d2_b = (1.0d0 - nu) / dble(2.0)
        do g_i_ = 1, g_p_
          g_db(g_i_, 3, 3) = d2_b * g_cb(g_i_)
        enddo
        db(3, 3) = (1.0d0 - nu) / dble(2.0) * cb
C--------
C
C bending elastic matrix
C
        d2_b = 1.0d0 / (1.0d0 - nu * nu)
        d3_b = d2_b * esp
        d4_b = d2_b * e
        do g_i_ = 1, g_p_
          g_cb(g_i_) = d4_b * g_esp(g_i_) + d3_b * g_e(g_i_)
        enddo
        cb = e * esp / (1.0d0 - nu * nu)
C--------
        do g_i_ = 1, g_p_
          g_dm(g_i_, 1, 1) = g_cb(g_i_)
        enddo
        dm(1, 1) = cb
C--------
        do g_i_ = 1, g_p_
          g_dm(g_i_, 1, 2) = nu * g_cb(g_i_)
        enddo
        dm(1, 2) = nu * cb
C--------
        do g_i_ = 1, g_p_
          g_dm(g_i_, 1, 3) = 0.0d0
        enddo
        dm(1, 3) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_dm(g_i_, 2, 1) = g_dm(g_i_, 1, 2)
        enddo
        dm(2, 1) = dm(1, 2)
C--------
        do g_i_ = 1, g_p_
          g_dm(g_i_, 2, 2) = g_cb(g_i_)
        enddo
        dm(2, 2) = cb
C--------
        do g_i_ = 1, g_p_
          g_dm(g_i_, 2, 3) = 0.0d0
        enddo
        dm(2, 3) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_dm(g_i_, 3, 1) = 0.0d0
        enddo
        dm(3, 1) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_dm(g_i_, 3, 2) = 0.0d0
        enddo
        dm(3, 2) = 0.0d0
C--------
        d2_b = (1.0d0 - nu) / dble(2.0)
        do g_i_ = 1, g_p_
          g_dm(g_i_, 3, 3) = d2_b * g_cb(g_i_)
        enddo
        dm(3, 3) = (1.0d0 - nu) / dble(2.0) * cb
C--------
C
C dimension variables
C
        do g_i_ = 1, g_p_
          g_x21(g_i_) = -g_xl(g_i_, 1) + g_xl(g_i_, 2)
        enddo
        x21 = xl(2) - xl(1)
C--------
        do g_i_ = 1, g_p_
          g_y21(g_i_) = -g_yl(g_i_, 1) + g_yl(g_i_, 2)
        enddo
        y21 = yl(2) - yl(1)
C--------
        do g_i_ = 1, g_p_
          g_z21(g_i_) = -g_zl(g_i_, 1) + g_zl(g_i_, 2)
        enddo
        z21 = zl(2) - zl(1)
C--------
        do g_i_ = 1, g_p_
          g_x32(g_i_) = -g_xl(g_i_, 2) + g_xl(g_i_, 3)
        enddo
        x32 = xl(3) - xl(2)
C--------
        do g_i_ = 1, g_p_
          g_y32(g_i_) = -g_yl(g_i_, 2) + g_yl(g_i_, 3)
        enddo
        y32 = yl(3) - yl(2)
C--------
        do g_i_ = 1, g_p_
          g_z32(g_i_) = -g_zl(g_i_, 2) + g_zl(g_i_, 3)
        enddo
        z32 = zl(3) - zl(2)
C--------
        x13 = xl(1) - xl(3)
        y13 = yl(1) - yl(3)
        z13 = zl(1) - zl(3)
C triangle in space : we compute the length of one side and the distance of the
C opposing node to that side to compute the area
        d4_b = y21 + y21
        d5_b = x21 + x21
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_y21(g_i_) + d5_b * g_x21(g_i_)
        enddo
        d2_w = x21 * x21 + y21 * y21
        d4_b = z21 + z21
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_z21(g_i_) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + z21 * z21
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +299)
        endif
        do g_i_ = 1, g_p_
          g_rlr(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        rlr = d2_v
C--------
        rlb = dsqrt(x32 * x32 + y32 * y32 + z32 * z32)
        bpr = dsqrt((x21 * x32 + y21 * y32 + z21 * z32) ** 2) / rlr
        area = rlr * dsqrt(rlb ** 2 - bpr ** 2) / 2.0d+00
C direction cosines of the local system . X' is directed parallel 
C to the 2-1 side
C Z' is the external normal (counterclockwise). Y' computed as Z' x X'
        d3_v = x21 / rlr
        d2_b = 1.0d0 / rlr
        d3_b = -d3_v / rlr
        do g_i_ = 1, g_p_
          g_xp(g_i_, 1) = d3_b * g_rlr(g_i_) + d2_b * g_x21(g_i_)
        enddo
        xp(1) = d3_v
C--------
        d3_v = y21 / rlr
        d2_b = 1.0d0 / rlr
        d3_b = -d3_v / rlr
        do g_i_ = 1, g_p_
          g_xp(g_i_, 2) = d3_b * g_rlr(g_i_) + d2_b * g_y21(g_i_)
        enddo
        xp(2) = d3_v
C--------
        d3_v = z21 / rlr
        d2_b = 1.0d0 / rlr
        d3_b = -d3_v / rlr
        do g_i_ = 1, g_p_
          g_xp(g_i_, 3) = d3_b * g_rlr(g_i_) + d2_b * g_z21(g_i_)
        enddo
        xp(3) = d3_v
C--------
C Z'
        do g_i_ = 1, g_p_
          g_zp(g_i_, 1) = -z21 * g_y32(g_i_) + (-y32) * g_z21(g_i_) + y2
     *1 * g_z32(g_i_) + z32 * g_y21(g_i_)
        enddo
        zp(1) = y21 * z32 - z21 * y32
C--------
        do g_i_ = 1, g_p_
          g_zp(g_i_, 2) = -x21 * g_z32(g_i_) + (-z32) * g_x21(g_i_) + z2
     *1 * g_x32(g_i_) + x32 * g_z21(g_i_)
        enddo
        zp(2) = z21 * x32 - x21 * z32
C--------
        do g_i_ = 1, g_p_
          g_zp(g_i_, 3) = -y21 * g_x32(g_i_) + (-x32) * g_y21(g_i_) + x2
     *1 * g_y32(g_i_) + y32 * g_x21(g_i_)
        enddo
        zp(3) = x21 * y32 - y21 * x32
C--------
        d4_b = zp(2) + zp(2)
        d5_b = zp(1) + zp(1)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_zp(g_i_, 2) + d5_b * g_zp(g_i_, 1)
        enddo
        d2_w = zp(1) * zp(1) + zp(2) * zp(2)
        d4_b = zp(3) + zp(3)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_zp(g_i_, 3) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + zp(3) * zp(3)
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +372)
        endif
        do g_i_ = 1, g_p_
          g_zlr(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        zlr = d2_v
C--------
        d3_v = zp(1) / zlr
        d2_b = 1.0d0 / zlr
        d3_b = -d3_v / zlr
        do g_i_ = 1, g_p_
          g_zp(g_i_, 1) = d3_b * g_zlr(g_i_) + d2_b * g_zp(g_i_, 1)
        enddo
        zp(1) = d3_v
C--------
        d3_v = zp(2) / zlr
        d2_b = 1.0d0 / zlr
        d3_b = -d3_v / zlr
        do g_i_ = 1, g_p_
          g_zp(g_i_, 2) = d3_b * g_zlr(g_i_) + d2_b * g_zp(g_i_, 2)
        enddo
        zp(2) = d3_v
C--------
        d3_v = zp(3) / zlr
        d2_b = 1.0d0 / zlr
        d3_b = -d3_v / zlr
        do g_i_ = 1, g_p_
          g_zp(g_i_, 3) = d3_b * g_zlr(g_i_) + d2_b * g_zp(g_i_, 3)
        enddo
        zp(3) = d3_v
C--------
C Y'
        do g_i_ = 1, g_p_
          g_yp(g_i_, 1) = -zp(3) * g_xp(g_i_, 2) + (-xp(2)) * g_zp(g_i_,
     * 3) + zp(2) * g_xp(g_i_, 3) + xp(3) * g_zp(g_i_, 2)
        enddo
        yp(1) = zp(2) * xp(3) - zp(3) * xp(2)
C--------
        do g_i_ = 1, g_p_
          g_yp(g_i_, 2) = -zp(1) * g_xp(g_i_, 3) + (-xp(3)) * g_zp(g_i_,
     * 1) + zp(3) * g_xp(g_i_, 1) + xp(1) * g_zp(g_i_, 3)
        enddo
        yp(2) = zp(3) * xp(1) - zp(1) * xp(3)
C--------
        do g_i_ = 1, g_p_
          g_yp(g_i_, 3) = -zp(2) * g_xp(g_i_, 1) + (-xp(1)) * g_zp(g_i_,
     * 2) + zp(1) * g_xp(g_i_, 2) + xp(2) * g_zp(g_i_, 1)
        enddo
        yp(3) = zp(1) * xp(2) - zp(2) * xp(1)
C--------
        d4_b = yp(2) + yp(2)
        d5_b = yp(1) + yp(1)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_yp(g_i_, 2) + d5_b * g_yp(g_i_, 1)
        enddo
        d2_w = yp(1) * yp(1) + yp(2) * yp(2)
        d4_b = yp(3) + yp(3)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_yp(g_i_, 3) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + yp(3) * yp(3)
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +439)
        endif
        do g_i_ = 1, g_p_
          g_ylr(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        ylr = d2_v
C--------
        d3_v = yp(1) / ylr
        d2_b = 1.0d0 / ylr
        d3_b = -d3_v / ylr
        do g_i_ = 1, g_p_
          g_yp(g_i_, 1) = d3_b * g_ylr(g_i_) + d2_b * g_yp(g_i_, 1)
        enddo
        yp(1) = d3_v
C--------
        d3_v = yp(2) / ylr
        d2_b = 1.0d0 / ylr
        d3_b = -d3_v / ylr
        do g_i_ = 1, g_p_
          g_yp(g_i_, 2) = d3_b * g_ylr(g_i_) + d2_b * g_yp(g_i_, 2)
        enddo
        yp(2) = d3_v
C--------
        d3_v = yp(3) / ylr
        d2_b = 1.0d0 / ylr
        d3_b = -d3_v / ylr
        do g_i_ = 1, g_p_
          g_yp(g_i_, 3) = d3_b * g_ylr(g_i_) + d2_b * g_yp(g_i_, 3)
        enddo
        yp(3) = d3_v
C--------
C compute center of gravity
        d2_b = 1.0d0 / 3.0d+00
        do g_i_ = 1, g_p_
          g_xcg(g_i_) = d2_b * g_xl(g_i_, 3) + d2_b * g_xl(g_i_, 2) + d2
     *_b * g_xl(g_i_, 1)
        enddo
        xcg = (xl(1) + xl(2) + xl(3)) / 3.0d+00
C--------
        d2_b = 1.0d0 / 3.0d+00
        do g_i_ = 1, g_p_
          g_ycg(g_i_) = d2_b * g_yl(g_i_, 3) + d2_b * g_yl(g_i_, 2) + d2
     *_b * g_yl(g_i_, 1)
        enddo
        ycg = (yl(1) + yl(2) + yl(3)) / 3.0d+00
C--------
        d2_b = 1.0d0 / 3.0d+00
        do g_i_ = 1, g_p_
          g_zcg(g_i_) = d2_b * g_zl(g_i_, 3) + d2_b * g_zl(g_i_, 2) + d2
     *_b * g_zl(g_i_, 1)
        enddo
        zcg = (zl(1) + zl(2) + zl(3)) / 3.0d+00
C--------
C compute local coordinates 
        do 99999 i = 1, 3
          do g_i_ = 1, g_p_
            g_xlcg(g_i_) = -g_xcg(g_i_) + g_xl(g_i_, i)
          enddo
          xlcg = xl(i) - xcg
C--------
          do g_i_ = 1, g_p_
            g_ylcg(g_i_) = -g_ycg(g_i_) + g_yl(g_i_, i)
          enddo
          ylcg = yl(i) - ycg
C--------
          do g_i_ = 1, g_p_
            g_zlcg(g_i_) = -g_zcg(g_i_) + g_zl(g_i_, i)
          enddo
          zlcg = zl(i) - zcg
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = xp(2) * g_ylcg(g_i_) + ylcg * g_xp(g_i_, 2) +
     * xp(1) * g_xlcg(g_i_) + xlcg * g_xp(g_i_, 1)
          enddo
          d1_w = xp(1) * xlcg + xp(2) * ylcg
          do g_i_ = 1, g_p_
            g_xlp(g_i_, i) = xp(3) * g_zlcg(g_i_) + zlcg * g_xp(g_i_, 3)
     * + g_d1_w(g_i_)
          enddo
          xlp(i) = d1_w + xp(3) * zlcg
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = yp(2) * g_ylcg(g_i_) + ylcg * g_yp(g_i_, 2) +
     * yp(1) * g_xlcg(g_i_) + xlcg * g_yp(g_i_, 1)
          enddo
          d1_w = yp(1) * xlcg + yp(2) * ylcg
          do g_i_ = 1, g_p_
            g_ylp(g_i_, i) = yp(3) * g_zlcg(g_i_) + zlcg * g_yp(g_i_, 3)
     * + g_d1_w(g_i_)
          enddo
          ylp(i) = d1_w + yp(3) * zlcg
C--------
          zlp(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg
43        continue
99999   continue
C
C zero stiffness matrices
C
        do 99997 i = 1, 18
          do 99998 j = 1, 18
            do g_i_ = 1, g_p_
              g_rk(g_i_, i, j) = 0.0d0
            enddo
            rk(i, j) = 0.0d+00
C--------
            do g_i_ = 1, g_p_
              g_rkm(g_i_, i, j) = 0.0d0
            enddo
            rkm(i, j) = 0.0d+00
C--------
15          do g_i_ = 1, g_p_
              g_rkb(g_i_, i, j) = 0.0d0
            enddo
            rkb(i, j) = 0.0d+00
C--------
99998     continue
99997   continue
C
C form local basic bending stiffness
C
        d1 = 1.0d+00
        call g_basico(g_p_, xlp, g_xlp, g_pmax_, ylp, g_ylp, g_pmax_, db
     *, g_db, g_pmax_, d1, clr, cqr, lb, rkb, g_rkb, g_pmax_, 18, status
     *)
C
C form local basic membrane stiffness
C
        d1 = 1.0d+00
        call g_sm3mb(g_p_, xlp, g_xlp, g_pmax_, ylp, g_ylp, g_pmax_, dm,
     * g_dm, g_pmax_, alpha, d1, le, rkm, g_rkm, g_pmax_, 18, status)
C
C form local higher order bending stiffness
C
        d1 = 1.0d+00
        call g_smcbh(g_p_, xlp, g_xlp, g_pmax_, ylp, g_ylp, g_pmax_, db,
     * g_db, g_pmax_, d1, lb, rkb, g_rkb, g_pmax_, 18, status)
C
C form local higher order membrane stiffness
C
        d1 = 0.32d+00
        call g_sm3mhe(g_p_, xlp, g_xlp, g_pmax_, ylp, g_ylp, g_pmax_, dm
     *, g_dm, g_pmax_, d1, le, rkm, g_rkm, g_pmax_, 18, status)
C
C add bending stiffness and membrane stiffness
C
        do 99995 i = 1, 18
          do 99996 j = 1, 18
            do g_i_ = 1, g_p_
              g_rk(g_i_, i, j) = g_rkm(g_i_, i, j) + g_rkb(g_i_, i, j)
            enddo
            rk(i, j) = rkb(i, j) + rkm(i, j)
C--------
16          continue
99996     continue
99995   continue
C
C rotate stiffness matrix from local coordinate system to global
C coordinate system in the case of linear FEM. In the case of
C nonlinear FEM with corotational method, do not perform this 
C transformation as the corotational routines expect a stiffness
C matrix in local coordinates.
C
        if (flag .eq. 1) then
C
C     compute local to global rotation matrix
C
          call g_rotation(g_p_, xp, g_xp, g_pmax_, yp, g_yp, g_pmax_, zp
     *, g_zp, g_pmax_, v1n, v2n, v3n, r1, g_r1, g_pmax_)
C
          call g_trirotation(g_p_, rk, g_rk, ldg_rk, r1, g_r1, g_pmax_)
C
        endif
C
        return
      end
C
