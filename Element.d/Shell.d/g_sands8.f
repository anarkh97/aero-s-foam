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
      subroutine gxsands8(xl,g_xl,yl,g_yl,zl,g_zl,
     *                    e,g_e,nu,h,g_h,
     *                    v,g_v,stress,g_stress,
     *                    strainflg, maxsze, maxstr, 
     *                    maxgus,elm,surface)
C-----------------------------------------------------------------
C This subroutine evaluates the von mises stress at the centroid
C of the triangle. it computes the values at the top and bottom
C and picks the maximum as output
C for the spacial 18 d.o.f tree node triangle obtained as
C a combination of the aqr bending triangle plus the membrane
C with driling d.o.f. developed by Felippa et al.
C
C MODIFIED: 10-07-97  K. H. Pierson
C Added stress calculations for: sigmaxx, sigmayy, sigmaxy,
C epsilonxx, epsilonyy, epsilonzz, epsilonxy and an equivalent
C strain similar to the vonmises stress. These stresses and strains
C can be calculated at the top, median or bottom surfaces.
C Stresses and strains are computed locally and then transformed
C to the global coordinate system.
C 
C-------------------------------------------------------------------*
C  CALLED BY : ThreeNodeShell.C 
C
C  SUBROUTINES CALLED:
C                      ROTATION
C                      TRIROT
C                      MEMBRA
C                      MOMEN
C         TRANSFORMATION
C                      VONMISES
C-------------------------------------------------------------------* 
C
C t = thickness of triangle
C
        integer maxsze, maxstr, maxgus
        integer elm, surface
        double precision e, nu
        double precision epsxx, epsyy, epszz, epsxy
        double precision rnxt, rnyt, rnxyt, rmxt, rmyt, rmxyt
        double precision t2, sx, sy, sxy
        double precision xl(*), yl(*), zl(*), h(*), v(*)
        double precision stress(maxsze, maxstr, maxgus)
        double precision db(3, 3), dm(3, 3)
        double precision xp(3), yp(3), zp(3), xlp(3), ylp(3), zlp(3)
        double precision r1(3, 3), dll(18)
        double precision xg(3), yg(3), zg(3), str(6)
        double precision cb, x21, y21, z21, x32, y32, z32, x13, y13, z13
        double precision rx, ry, rz, bx, by, bz, rlr, rlb, bpr, area, yl
     *r, zlr
        double precision ycg, xcg, zcg, xlcg, ylcg, zlcg, f, t
        double precision rmom(18, 3), rmem(18, 3)
        double precision rmx, rmy, rmxy, rnx, rny, rnxy, clr, cqr
        double precision rmmx, rmmy, rmmxy, rnnx, rnny, rnnxy, sbf
        double precision ebar
        integer lb(9), le(9), strainflg
        character*10 status
        integer i, j
        double precision d1
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_h, ldg_e, ldg_xl, ldg_yl, ldg_zl, ldg_v,
     * ldg_stress
c
c manually inserted - begin	     
c
        parameter (g_p_=1,ldg_xl=1,ldg_yl=1,ldg_zl=1,ldg_e=1,
     *             ldg_h=1,ldg_v=1,ldg_stress=1)
c
c manually inserted - end      
c
        double precision d1_p, d2_w, d2_v, d3_v, d4_v, d1_w, d7_b, d6_b,
     * d5_b, d2_b
        double precision d3_b, d4_b, g_rmom(g_pmax_, 18, 3), g_rmem(g_pm
     *ax_, 18, 3), g_t(g_pmax_), g_h(ldg_h, *), g_cb(g_pmax_), g_e(ldg_e
     *), g_db(g_pmax_, 3, 3), g_dm(g_pmax_, 3, 3)
        double precision g_x21(g_pmax_), g_xl(ldg_xl, *), g_y21(g_pmax_)
     *, g_yl(ldg_yl, *), g_z21(g_pmax_), g_zl(ldg_zl, *), g_x32(g_pmax_)
     *, g_y32(g_pmax_), g_z32(g_pmax_), g_rx(g_pmax_)
        double precision g_ry(g_pmax_), g_rz(g_pmax_), g_d2_w(g_pmax_), 
     *g_d1_w(g_pmax_), g_rlr(g_pmax_), g_xp(g_pmax_, 3), g_zp(g_pmax_, 3
     *), g_zlr(g_pmax_), g_yp(g_pmax_, 3), g_ylr(g_pmax_)
        double precision g_xcg(g_pmax_), g_ycg(g_pmax_), g_zcg(g_pmax_),
     * g_xlcg(g_pmax_), g_ylcg(g_pmax_), g_zlcg(g_pmax_), g_xlp(g_pmax_,
     * 3), g_ylp(g_pmax_, 3), g_dll(g_pmax_, 18), g_r1(g_pmax_, 3, 3)
        double precision g_v(ldg_v, *), g_rmx(g_pmax_), g_rmy(g_pmax_), 
     *g_rmxy(g_pmax_), g_rnx(g_pmax_), g_rny(g_pmax_), g_rnxy(g_pmax_), 
     *g_t2(g_pmax_), g_epsxx(g_pmax_), g_epsyy(g_pmax_)
        double precision g_epsxy(g_pmax_), g_epszz(g_pmax_), g_str(g_pma
     *x_, 6), g_stress(ldg_stress, maxsze, maxstr, maxgus), g_ebar(g_pma
     *x_), g_rmmx(g_pmax_), g_rmmy(g_pmax_), g_rmmxy(g_pmax_), g_rnnx(g_
     *pmax_), g_rnny(g_pmax_)
        double precision g_rnnxy(g_pmax_), g_rnxt(g_pmax_), g_rnyt(g_pma
     *x_), g_rnxyt(g_pmax_), g_rmxt(g_pmax_), g_rmyt(g_pmax_), g_rmxyt(g
     *_pmax_), g_sx(g_pmax_), g_sy(g_pmax_), g_sxy(g_pmax_)
        double precision g_sbf(g_pmax_)
        save g_sxy, g_sbf
        save g_rnny, g_rnnxy, g_rnxt, g_rnyt, g_rnxyt, g_rmxt, g_rmyt, g
     *_rmxyt, g_sx, g_sy
        save g_epsxx, g_epsyy, g_epsxy, g_epszz, g_str, g_ebar, g_rmmx, 
     *g_rmmy, g_rmmxy, g_rnnx
        save g_ylp, g_dll, g_r1, g_rmx, g_rmy, g_rmxy, g_rnx, g_rny, g_r
     *nxy, g_t2
        save g_zlr, g_yp, g_ylr, g_xcg, g_ycg, g_zcg, g_xlcg, g_ylcg, g_
     *zlcg, g_xlp
        save g_y32, g_z32, g_rx, g_ry, g_rz, g_d2_w, g_d1_w, g_rlr, g_xp
     *, g_zp
        save g_rmom, g_rmem, g_t, g_cb, g_db, g_dm, g_x21, g_y21, g_z21,
     * g_x32
        external g_vonmis
        external g_straineq
        external g_transform
        external g_momen
        external g_membra
        external g_rotation
        intrinsic dble
        data lb /3, 4, 5, 9, 10, 11, 15, 16, 17/
        data le /1, 2, 7, 8, 13, 14, 6, 12, 18/
C
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'sands8','g_sands8.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do 99999 i = 1, 18
          do g_i_ = 1, g_p_
            g_rmom(g_i_, i, 1) = 0.0d0
          enddo
          rmom(i, 1) = 0.0d0
C--------
          do g_i_ = 1, g_p_
            g_rmom(g_i_, i, 2) = 0.0d0
          enddo
          rmom(i, 2) = 0.0d0
C--------
          do g_i_ = 1, g_p_
            g_rmom(g_i_, i, 3) = 0.0d0
          enddo
          rmom(i, 3) = 0.0d0
C--------
          do g_i_ = 1, g_p_
            g_rmem(g_i_, i, 1) = 0.0d0
          enddo
          rmem(i, 1) = 0.0d0
C--------
          do g_i_ = 1, g_p_
            g_rmem(g_i_, i, 2) = 0.0d0
          enddo
          rmem(i, 2) = 0.0d0
C--------
          do g_i_ = 1, g_p_
            g_rmem(g_i_, i, 3) = 0.0d0
          enddo
          rmem(i, 3) = 0.0d0
C--------
5         continue
99999   continue
C
        f = 1.0d+00
        clr = 0.0d+00
        cqr = 1.0d+00
        do g_i_ = 1, g_p_
          g_t(g_i_) = g_h(g_i_, 1)
        enddo
        t = h(1)
C--------
C
C set the bending constitutive matrix
C
        d3_v = t * t
        d4_v = d3_v * t
        d3_b = 1.0d0 / (1.0d0 - nu * nu) * (1.0d0 / dble(12.0))
        d4_b = d3_b * d4_v
        d5_b = d3_b * e
        d6_b = d5_b * t
        d7_b = d5_b * d3_v + d6_b * t + d6_b * t
        do g_i_ = 1, g_p_
          g_cb(g_i_) = d7_b * g_t(g_i_) + d4_b * g_e(g_i_)
        enddo
        cb = e * d4_v / dble(12.0) / (1.0d0 - nu * nu)
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
        d2_b = dble(0.5) * (1.0d0 - nu)
        do g_i_ = 1, g_p_
          g_db(g_i_, 3, 3) = d2_b * g_cb(g_i_)
        enddo
        db(3, 3) = dble(0.5) * (1.0d0 - nu) * cb
C--------
C
C set the membrane constitutive matrix
C
        d2_b = 1.0d0 / (1.0d0 - nu * nu)
        d3_b = d2_b * t
        d4_b = d2_b * e
        do g_i_ = 1, g_p_
          g_cb(g_i_) = d4_b * g_t(g_i_) + d3_b * g_e(g_i_)
        enddo
        cb = e * t / (1.0d0 - nu * nu)
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
          g_dm(g_i_, 2, 1) = g_db(g_i_, 1, 2)
        enddo
        dm(2, 1) = db(1, 2)
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
        d2_b = dble(0.5) * (1.0d0 - nu)
        do g_i_ = 1, g_p_
          g_dm(g_i_, 3, 3) = d2_b * g_cb(g_i_)
        enddo
        dm(3, 3) = dble(0.5) * (1.0d0 - nu) * cb
C--------
C
C triangular dimension variables
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
C  triangle in space : we compute the length of one side and the distance of the
C  opposing node to that side to compute the area
        do g_i_ = 1, g_p_
          g_rx(g_i_) = g_x21(g_i_)
        enddo
        rx = x21
C--------
        do g_i_ = 1, g_p_
          g_ry(g_i_) = g_y21(g_i_)
        enddo
        ry = y21
C--------
        do g_i_ = 1, g_p_
          g_rz(g_i_) = g_z21(g_i_)
        enddo
        rz = z21
C--------
        bx = x32
        by = y32
        bz = z32
        d4_b = ry + ry
        d5_b = rx + rx
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_ry(g_i_) + d5_b * g_rx(g_i_)
        enddo
        d2_w = rx * rx + ry * ry
        d4_b = rz + rz
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_rz(g_i_) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + rz * rz
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +368)
        endif
        do g_i_ = 1, g_p_
          g_rlr(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        rlr = d2_v
C--------
        rlb = dsqrt(bx * bx + by * by + bz * bz)
        bpr = dsqrt((rx * bx + ry * by + rz * bz) ** 2) / rlr
        area = 0.5 * rlr * dsqrt(rlb * rlb - bpr * bpr)
C
C Direction cosines of the local system . X' is directed parallel to the 
C 2-1 side Z' is the external normal (counterclockwise). 
C Y' computed as Z' x X'
C
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
C
C Z' local axis
C
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
     +445)
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
C
C Y' local axis
C
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
     +514)
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
C
C center of gravity
C
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
C
C computing local coordinates
C
        do 99998 i = 1, 3
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
99998   continue
C
C Set Global axes
C
        xg(1) = 1.0
        xg(2) = 0.0
        xg(3) = 0.0
        yg(1) = 0.0
        yg(2) = 1.0
        yg(3) = 0.0
        zg(1) = 0.0
        zg(2) = 0.0
        zg(3) = 1.0
C
C computing nodal rotation matrix
C
        call g_rotation(g_p_, xp, g_xp, g_pmax_, yp, g_yp, g_pmax_, zp, 
     *g_zp, g_pmax_, xg, yg, zg, r1, g_r1, g_pmax_)
C
C compute the von mises stress
C rotate nodal displacements to local system
C
        do 99997 i = 1, 18
          do g_i_ = 1, g_p_
            g_dll(g_i_, i) = 0.0d0
          enddo
          dll(i) = 0.0d+00
C--------
270       continue
99997   continue
        do 99995 i = 1, 3
          do 99996 j = 1, 3
            do g_i_ = 1, g_p_
              g_dll(g_i_, i) = r1(j, i) * g_v(g_i_, j) + v(j) * g_r1(g_i
     *_, j, i) + g_dll(g_i_, i)
            enddo
            dll(i) = dll(i) + r1(j, i) * v(j)
C--------
            do g_i_ = 1, g_p_
              g_dll(g_i_, i + 3) = r1(j, i) * g_v(g_i_, j + 3) + v(j + 3
     *) * g_r1(g_i_, j, i) + g_dll(g_i_, i + 3)
            enddo
            dll(i + 3) = dll(i + 3) + r1(j, i) * v(j + 3)
C--------
            do g_i_ = 1, g_p_
              g_dll(g_i_, i + 6) = r1(j, i) * g_v(g_i_, j + 6) + v(j + 6
     *) * g_r1(g_i_, j, i) + g_dll(g_i_, i + 6)
            enddo
            dll(i + 6) = dll(i + 6) + r1(j, i) * v(j + 6)
C--------
            do g_i_ = 1, g_p_
              g_dll(g_i_, i + 9) = r1(j, i) * g_v(g_i_, j + 9) + v(j + 9
     *) * g_r1(g_i_, j, i) + g_dll(g_i_, i + 9)
            enddo
            dll(i + 9) = dll(i + 9) + r1(j, i) * v(j + 9)
C--------
            do g_i_ = 1, g_p_
              g_dll(g_i_, i + 12) = r1(j, i) * g_v(g_i_, j + 12) + v(j +
     * 12) * g_r1(g_i_, j, i) + g_dll(g_i_, i + 12)
            enddo
            dll(i + 12) = dll(i + 12) + r1(j, i) * v(j + 12)
C--------
            do g_i_ = 1, g_p_
              g_dll(g_i_, i + 15) = r1(j, i) * g_v(g_i_, j + 15) + v(j +
     * 15) * g_r1(g_i_, j, i) + g_dll(g_i_, i + 15)
            enddo
            dll(i + 15) = dll(i + 15) + r1(j, i) * v(j + 15)
C--------
280         continue
99996     continue
99995   continue
C
C      compute centroidal membrane strains
        d1 = 1.5d+00
        call g_membra(g_p_, xlp, g_xlp, g_pmax_, ylp, g_ylp, g_pmax_, d1
     *, le, rmem, g_rmem, g_pmax_, status)
C
C      compute centroidal bending strains (curvatures (1/radius))
        call g_momen(g_p_, xlp, g_xlp, g_pmax_, ylp, g_ylp, g_pmax_, lb,
     * rmom, g_rmom, g_pmax_, status)
C
        do g_i_ = 1, g_p_
          g_rmx(g_i_) = 0.0d0
        enddo
        rmx = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_rmy(g_i_) = 0.0d0
        enddo
        rmy = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_rmxy(g_i_) = 0.0d0
        enddo
        rmxy = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_rnx(g_i_) = 0.0d0
        enddo
        rnx = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_rny(g_i_) = 0.0d0
        enddo
        rny = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_rnxy(g_i_) = 0.0d0
        enddo
        rnxy = 0.0d0
C--------
        do 99994 j = 1, 18
          do g_i_ = 1, g_p_
            g_rmx(g_i_) = rmom(j, 1) * g_dll(g_i_, j) + dll(j) * g_rmom(
     *g_i_, j, 1) + g_rmx(g_i_)
          enddo
          rmx = rmx + rmom(j, 1) * dll(j)
C--------
          do g_i_ = 1, g_p_
            g_rmy(g_i_) = rmom(j, 2) * g_dll(g_i_, j) + dll(j) * g_rmom(
     *g_i_, j, 2) + g_rmy(g_i_)
          enddo
          rmy = rmy + rmom(j, 2) * dll(j)
C--------
          do g_i_ = 1, g_p_
            g_rmxy(g_i_) = rmom(j, 3) * g_dll(g_i_, j) + dll(j) * g_rmom
     *(g_i_, j, 3) + g_rmxy(g_i_)
          enddo
          rmxy = rmxy + rmom(j, 3) * dll(j)
C--------
          do g_i_ = 1, g_p_
            g_rnx(g_i_) = rmem(j, 1) * g_dll(g_i_, j) + dll(j) * g_rmem(
     *g_i_, j, 1) + g_rnx(g_i_)
          enddo
          rnx = rnx + rmem(j, 1) * dll(j)
C--------
          do g_i_ = 1, g_p_
            g_rny(g_i_) = rmem(j, 2) * g_dll(g_i_, j) + dll(j) * g_rmem(
     *g_i_, j, 2) + g_rny(g_i_)
          enddo
          rny = rny + rmem(j, 2) * dll(j)
C--------
          do g_i_ = 1, g_p_
            g_rnxy(g_i_) = rmem(j, 3) * g_dll(g_i_, j) + dll(j) * g_rmem
     *(g_i_, j, 3) + g_rnxy(g_i_)
          enddo
          rnxy = rnxy + rmem(j, 3) * dll(j)
C--------
290       continue
99994   continue
C
C      COMPUTE VON MISES STRAIN RESULTANT
C
        if (strainflg .eq. 1) then
C
          d2_b = dble(0.5)
          do g_i_ = 1, g_p_
            g_t2(g_i_) = d2_b * g_t(g_i_)
          enddo
          t2 = dble(0.5) * t
C--------
C
          do g_i_ = 1, g_p_
            g_rmx(g_i_) = t2 * g_rmx(g_i_) + rmx * g_t2(g_i_)
          enddo
          rmx = t2 * rmx
C--------
          do g_i_ = 1, g_p_
            g_rmy(g_i_) = t2 * g_rmy(g_i_) + rmy * g_t2(g_i_)
          enddo
          rmy = t2 * rmy
C--------
          do g_i_ = 1, g_p_
            g_rmxy(g_i_) = t2 * g_rmxy(g_i_) + rmxy * g_t2(g_i_)
          enddo
          rmxy = t2 * rmxy
C--------
C
C ... COMPUTE STRAINS AT MEDIAN SURFACE
          if (surface .eq. 2) then
            do g_i_ = 1, g_p_
              g_epsxx(g_i_) = g_rnx(g_i_)
            enddo
            epsxx = rnx
C--------
            do g_i_ = 1, g_p_
              g_epsyy(g_i_) = g_rny(g_i_)
            enddo
            epsyy = rny
C--------
            d2_b = dble(0.5)
            do g_i_ = 1, g_p_
              g_epsxy(g_i_) = d2_b * g_rnxy(g_i_)
            enddo
            epsxy = dble(0.5) * rnxy
C--------
C ... COMPUTE STRAINS AT BOTTOM SURFACE
          else
            if (surface .eq. 3) then
              do g_i_ = 1, g_p_
                g_epsxx(g_i_) = -g_rmx(g_i_) + g_rnx(g_i_)
              enddo
              epsxx = rnx - rmx
C--------
              do g_i_ = 1, g_p_
                g_epsyy(g_i_) = -g_rmy(g_i_) + g_rny(g_i_)
              enddo
              epsyy = rny - rmy
C--------
              d2_b = dble(0.5)
              do g_i_ = 1, g_p_
                g_epsxy(g_i_) = -d2_b * g_rmxy(g_i_) + d2_b * g_rnxy(g_i
     *_)
              enddo
              epsxy = dble(0.5) * (rnxy - rmxy)
C--------
C ... COMPUTE STRAINS AT TOP SURFACE
            else
              do g_i_ = 1, g_p_
                g_epsxx(g_i_) = g_rmx(g_i_) + g_rnx(g_i_)
              enddo
              epsxx = rnx + rmx
C--------
              do g_i_ = 1, g_p_
                g_epsyy(g_i_) = g_rmy(g_i_) + g_rny(g_i_)
              enddo
              epsyy = rny + rmy
C--------
              d2_b = dble(0.5)
              do g_i_ = 1, g_p_
                g_epsxy(g_i_) = d2_b * g_rmxy(g_i_) + d2_b * g_rnxy(g_i_
     *)
              enddo
              epsxy = dble(0.5) * (rnxy + rmxy)
C--------
            endif
          endif
C
          d2_b = -nu / (1.0d0 - nu)
          do g_i_ = 1, g_p_
            g_epszz(g_i_) = d2_b * g_epsyy(g_i_) + d2_b * g_epsxx(g_i_)
          enddo
          epszz = -nu / (1.0d0 - nu) * (epsxx + epsyy)
C--------
C
          do g_i_ = 1, g_p_
            g_str(g_i_, 1) = g_epsxx(g_i_)
          enddo
          str(1) = epsxx
C--------
          do g_i_ = 1, g_p_
            g_str(g_i_, 2) = g_epsyy(g_i_)
          enddo
          str(2) = epsyy
C--------
          do g_i_ = 1, g_p_
            g_str(g_i_, 3) = g_epszz(g_i_)
          enddo
          str(3) = epszz
C--------
          do g_i_ = 1, g_p_
            g_str(g_i_, 4) = g_epsxy(g_i_)
          enddo
          str(4) = epsxy
C--------
          do g_i_ = 1, g_p_
            g_str(g_i_, 5) = 0.0d0
          enddo
          str(5) = 0.0d0
C--------
          do g_i_ = 1, g_p_
            g_str(g_i_, 6) = 0.0d0
          enddo
          str(6) = 0.0d0
C--------
C
          call g_transform(g_p_, xp, g_xp, g_pmax_, yp, g_yp, g_pmax_, z
     *p, g_zp, g_pmax_, xg, yg, zg, str, g_str, g_pmax_)
C
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 1, 1) = g_str(g_i_, 1)
          enddo
          stress(elm, 1, 1) = str(1)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 1, 2) = g_str(g_i_, 1)
          enddo
          stress(elm, 1, 2) = str(1)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 1, 3) = g_str(g_i_, 1)
          enddo
          stress(elm, 1, 3) = str(1)
C--------
C
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 2, 1) = g_str(g_i_, 2)
          enddo
          stress(elm, 2, 1) = str(2)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 2, 2) = g_str(g_i_, 2)
          enddo
          stress(elm, 2, 2) = str(2)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 2, 3) = g_str(g_i_, 2)
          enddo
          stress(elm, 2, 3) = str(2)
C--------
C
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 3, 1) = g_str(g_i_, 3)
          enddo
          stress(elm, 3, 1) = str(3)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 3, 2) = g_str(g_i_, 3)
          enddo
          stress(elm, 3, 2) = str(3)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 3, 3) = g_str(g_i_, 3)
          enddo
          stress(elm, 3, 3) = str(3)
C--------
C
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 4, 1) = g_str(g_i_, 4)
          enddo
          stress(elm, 4, 1) = str(4)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 4, 2) = g_str(g_i_, 4)
          enddo
          stress(elm, 4, 2) = str(4)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 4, 3) = g_str(g_i_, 4)
          enddo
          stress(elm, 4, 3) = str(4)
C--------
C
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 5, 1) = g_str(g_i_, 5)
          enddo
          stress(elm, 5, 1) = str(5)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 5, 2) = g_str(g_i_, 5)
          enddo
          stress(elm, 5, 2) = str(5)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 5, 3) = g_str(g_i_, 5)
          enddo
          stress(elm, 5, 3) = str(5)
C--------
C
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 6, 1) = g_str(g_i_, 6)
          enddo
          stress(elm, 6, 1) = str(6)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 6, 2) = g_str(g_i_, 6)
          enddo
          stress(elm, 6, 2) = str(6)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 6, 3) = g_str(g_i_, 6)
          enddo
          stress(elm, 6, 3) = str(6)
C--------
C
          call g_straineq(g_p_, rmx, g_rmx, g_pmax_, rmy, g_rmy, g_pmax_
     *, rmxy, g_rmxy, g_pmax_, rnx, g_rnx, g_pmax_, rny, g_rny, g_pmax_,
     * rnxy, g_rnxy, g_pmax_, t, g_t, g_pmax_, surface, ebar, g_ebar, g_
     *pmax_)
C
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 7, 1) = g_ebar(g_i_)
          enddo
          stress(elm, 7, 1) = ebar
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 7, 2) = g_ebar(g_i_)
          enddo
          stress(elm, 7, 2) = ebar
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 7, 3) = g_ebar(g_i_)
          enddo
          stress(elm, 7, 3) = ebar
C--------
C
          return
        endif
C
C     compute centroidal stress resultants
C
C     bending resultants
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = db(1, 2) * g_rmy(g_i_) + rmy * g_db(g_i_, 1, 2)
     * + db(1, 1) * g_rmx(g_i_) + rmx * g_db(g_i_, 1, 1)
        enddo
        d1_w = db(1, 1) * rmx + db(1, 2) * rmy
        do g_i_ = 1, g_p_
          g_rmmx(g_i_) = db(1, 3) * g_rmxy(g_i_) + rmxy * g_db(g_i_, 1, 
     *3) + g_d1_w(g_i_)
        enddo
        rmmx = d1_w + db(1, 3) * rmxy
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = db(2, 2) * g_rmy(g_i_) + rmy * g_db(g_i_, 2, 2)
     * + db(2, 1) * g_rmx(g_i_) + rmx * g_db(g_i_, 2, 1)
        enddo
        d1_w = db(2, 1) * rmx + db(2, 2) * rmy
        do g_i_ = 1, g_p_
          g_rmmy(g_i_) = db(2, 3) * g_rmxy(g_i_) + rmxy * g_db(g_i_, 2, 
     *3) + g_d1_w(g_i_)
        enddo
        rmmy = d1_w + db(2, 3) * rmxy
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = db(3, 2) * g_rmy(g_i_) + rmy * g_db(g_i_, 3, 2)
     * + db(3, 1) * g_rmx(g_i_) + rmx * g_db(g_i_, 3, 1)
        enddo
        d1_w = db(3, 1) * rmx + db(3, 2) * rmy
        do g_i_ = 1, g_p_
          g_rmmxy(g_i_) = db(3, 3) * g_rmxy(g_i_) + rmxy * g_db(g_i_, 3,
     * 3) + g_d1_w(g_i_)
        enddo
        rmmxy = d1_w + db(3, 3) * rmxy
C--------
C
C     membrane resultants
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = dm(1, 2) * g_rny(g_i_) + rny * g_dm(g_i_, 1, 2)
     * + dm(1, 1) * g_rnx(g_i_) + rnx * g_dm(g_i_, 1, 1)
        enddo
        d1_w = dm(1, 1) * rnx + dm(1, 2) * rny
        do g_i_ = 1, g_p_
          g_rnnx(g_i_) = dm(1, 3) * g_rnxy(g_i_) + rnxy * g_dm(g_i_, 1, 
     *3) + g_d1_w(g_i_)
        enddo
        rnnx = d1_w + dm(1, 3) * rnxy
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = dm(2, 2) * g_rny(g_i_) + rny * g_dm(g_i_, 2, 2)
     * + dm(2, 1) * g_rnx(g_i_) + rnx * g_dm(g_i_, 2, 1)
        enddo
        d1_w = dm(2, 1) * rnx + dm(2, 2) * rny
        do g_i_ = 1, g_p_
          g_rnny(g_i_) = dm(2, 3) * g_rnxy(g_i_) + rnxy * g_dm(g_i_, 2, 
     *3) + g_d1_w(g_i_)
        enddo
        rnny = d1_w + dm(2, 3) * rnxy
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = dm(3, 2) * g_rny(g_i_) + rny * g_dm(g_i_, 3, 2)
     * + dm(3, 1) * g_rnx(g_i_) + rnx * g_dm(g_i_, 3, 1)
        enddo
        d1_w = dm(3, 1) * rnx + dm(3, 2) * rny
        do g_i_ = 1, g_p_
          g_rnnxy(g_i_) = dm(3, 3) * g_rnxy(g_i_) + rnxy * g_dm(g_i_, 3,
     * 3) + g_d1_w(g_i_)
        enddo
        rnnxy = d1_w + dm(3, 3) * rnxy
C--------
C
        d2_v = abs(t)
        if (t .gt. 0.0d0) then
           d1_p =  1.0d0
        else if (t .lt. 0.0d0) then
           d1_p = -1.0d0
        else
           call ehufDO (3,t, d2_v, d1_p,
     +g_ehfid,
     +1091)
        endif
        do g_i_ = 1, g_p_
          g_t(g_i_) = d1_p * g_t(g_i_)
        enddo
        t = d2_v
C--------
        d2_b = t + t
        do g_i_ = 1, g_p_
          g_t2(g_i_) = d2_b * g_t(g_i_)
        enddo
        t2 = t * t
C--------
C
        d3_v = rnnx / t
        d2_b = 1.0d0 / t
        d3_b = -d3_v / t
        do g_i_ = 1, g_p_
          g_rnxt(g_i_) = d3_b * g_t(g_i_) + d2_b * g_rnnx(g_i_)
        enddo
        rnxt = d3_v
C--------
        d3_v = rnny / t
        d2_b = 1.0d0 / t
        d3_b = -d3_v / t
        do g_i_ = 1, g_p_
          g_rnyt(g_i_) = d3_b * g_t(g_i_) + d2_b * g_rnny(g_i_)
        enddo
        rnyt = d3_v
C--------
        d3_v = rnnxy / t
        d2_b = 1.0d0 / t
        d3_b = -d3_v / t
        do g_i_ = 1, g_p_
          g_rnxyt(g_i_) = d3_b * g_t(g_i_) + d2_b * g_rnnxy(g_i_)
        enddo
        rnxyt = d3_v
C--------
C
        d4_v = dble(6.0) * rmmx / t2
        d3_b = -d4_v / t2
        d4_b = 1.0d0 / t2 * dble(6.0)
        do g_i_ = 1, g_p_
          g_rmxt(g_i_) = d3_b * g_t2(g_i_) + d4_b * g_rmmx(g_i_)
        enddo
        rmxt = d4_v
C--------
        d4_v = dble(6.0) * rmmy / t2
        d3_b = -d4_v / t2
        d4_b = 1.0d0 / t2 * dble(6.0)
        do g_i_ = 1, g_p_
          g_rmyt(g_i_) = d3_b * g_t2(g_i_) + d4_b * g_rmmy(g_i_)
        enddo
        rmyt = d4_v
C--------
        d4_v = dble(6.0) * rmmxy / t2
        d3_b = -d4_v / t2
        d4_b = 1.0d0 / t2 * dble(6.0)
        do g_i_ = 1, g_p_
          g_rmxyt(g_i_) = d3_b * g_t2(g_i_) + d4_b * g_rmmxy(g_i_)
        enddo
        rmxyt = d4_v
C--------
C
C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT TOP SURFACE (DEFAULT)
        do g_i_ = 1, g_p_
          g_sx(g_i_) = g_rmxt(g_i_) + g_rnxt(g_i_)
        enddo
        sx = rnxt + rmxt
C--------
        do g_i_ = 1, g_p_
          g_sy(g_i_) = g_rmyt(g_i_) + g_rnyt(g_i_)
        enddo
        sy = rnyt + rmyt
C--------
        do g_i_ = 1, g_p_
          g_sxy(g_i_) = g_rmxyt(g_i_) + g_rnxyt(g_i_)
        enddo
        sxy = rnxyt + rmxyt
C--------
C
C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT MEDIAN SURFACE
        if (surface .eq. 2) then
          do g_i_ = 1, g_p_
            g_sx(g_i_) = g_rnxt(g_i_)
          enddo
          sx = rnxt
C--------
          do g_i_ = 1, g_p_
            g_sy(g_i_) = g_rnyt(g_i_)
          enddo
          sy = rnyt
C--------
          do g_i_ = 1, g_p_
            g_sxy(g_i_) = g_rnxyt(g_i_)
          enddo
          sxy = rnxyt
C--------
        endif
C
C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT BOTTOM SURFACE
        if (surface .eq. 3) then
          do g_i_ = 1, g_p_
            g_sx(g_i_) = -g_rmxt(g_i_) + g_rnxt(g_i_)
          enddo
          sx = rnxt - rmxt
C--------
          do g_i_ = 1, g_p_
            g_sy(g_i_) = -g_rmyt(g_i_) + g_rnyt(g_i_)
          enddo
          sy = rnyt - rmyt
C--------
          do g_i_ = 1, g_p_
            g_sxy(g_i_) = -g_rmxyt(g_i_) + g_rnxyt(g_i_)
          enddo
          sxy = rnxyt - rmxyt
C--------
        endif
C
        do g_i_ = 1, g_p_
          g_str(g_i_, 1) = g_sx(g_i_)
        enddo
        str(1) = sx
C--------
        do g_i_ = 1, g_p_
          g_str(g_i_, 2) = g_sy(g_i_)
        enddo
        str(2) = sy
C--------
        do g_i_ = 1, g_p_
          g_str(g_i_, 3) = 0.0d0
        enddo
        str(3) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_str(g_i_, 4) = g_sxy(g_i_)
        enddo
        str(4) = sxy
C--------
        do g_i_ = 1, g_p_
          g_str(g_i_, 5) = 0.0d0
        enddo
        str(5) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_str(g_i_, 6) = 0.0d0
        enddo
        str(6) = 0.0d0
C--------
C
        call g_transform(g_p_, xp, g_xp, g_pmax_, yp, g_yp, g_pmax_, zp,
     * g_zp, g_pmax_, xg, yg, zg, str, g_str, g_pmax_)
C
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 1, 1) = g_str(g_i_, 1)
        enddo
        stress(elm, 1, 1) = str(1)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 1, 2) = g_str(g_i_, 1)
        enddo
        stress(elm, 1, 2) = str(1)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 1, 3) = g_str(g_i_, 1)
        enddo
        stress(elm, 1, 3) = str(1)
C--------
C
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 2, 1) = g_str(g_i_, 2)
        enddo
        stress(elm, 2, 1) = str(2)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 2, 2) = g_str(g_i_, 2)
        enddo
        stress(elm, 2, 2) = str(2)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 2, 3) = g_str(g_i_, 2)
        enddo
        stress(elm, 2, 3) = str(2)
C--------
C
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 3, 1) = g_str(g_i_, 3)
        enddo
        stress(elm, 3, 1) = str(3)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 3, 2) = g_str(g_i_, 3)
        enddo
        stress(elm, 3, 2) = str(3)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 3, 3) = g_str(g_i_, 3)
        enddo
        stress(elm, 3, 3) = str(3)
C--------
C
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 4, 1) = g_str(g_i_, 4)
        enddo
        stress(elm, 4, 1) = str(4)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 4, 2) = g_str(g_i_, 4)
        enddo
        stress(elm, 4, 2) = str(4)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 4, 3) = g_str(g_i_, 4)
        enddo
        stress(elm, 4, 3) = str(4)
C--------
C
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 5, 1) = g_str(g_i_, 5)
        enddo
        stress(elm, 5, 1) = str(5)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 5, 2) = g_str(g_i_, 5)
        enddo
        stress(elm, 5, 2) = str(5)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 5, 3) = g_str(g_i_, 5)
        enddo
        stress(elm, 5, 3) = str(5)
C--------
C
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 6, 1) = g_str(g_i_, 6)
        enddo
        stress(elm, 6, 1) = str(6)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 6, 2) = g_str(g_i_, 6)
        enddo
        stress(elm, 6, 2) = str(6)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 6, 3) = g_str(g_i_, 6)
        enddo
        stress(elm, 6, 3) = str(6)
C--------
C
C     compute von mises stress resultant
C
        call g_vonmis(g_p_, rmmx, g_rmmx, g_pmax_, rmmy, g_rmmy, g_pmax_
     *, rmmxy, g_rmmxy, g_pmax_, rnnx, g_rnnx, g_pmax_, rnny, g_rnny, g_
     *pmax_, rnnxy, g_rnnxy, g_pmax_, t, g_t, g_pmax_, sbf, g_sbf, g_pma
     *x_, surface)
C
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 7, 1) = g_sbf(g_i_)
        enddo
        stress(elm, 7, 1) = sbf
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 7, 2) = g_sbf(g_i_)
        enddo
        stress(elm, 7, 2) = sbf
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 7, 3) = g_sbf(g_i_)
        enddo
        stress(elm, 7, 3) = sbf
C--------
C
        return
      end
