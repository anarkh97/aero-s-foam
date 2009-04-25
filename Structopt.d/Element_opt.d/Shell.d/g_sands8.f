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
     *                    maxgus,elm,surface,alpha,beta,thrmStr,
     *                    dthrmStr)
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
C       implicit none
        integer maxsze, maxstr, maxgus
        integer elm, surface
        double precision e, nu,alpha,beta,thrmStr,dthrmStr
        double precision epsxx, epsyy, epszz, epsxy
        double precision rnxt, rnyt, rnxyt, rmxt, rmyt, rmxyt
        double precision t2, sx, sy, sxy
        double precision xl(3), yl(3), zl(3), h(3), v(18)
        double precision stress(maxsze, maxstr, maxgus)
        double precision db(3, 3), dm(3, 3)
        double precision xp(3), yp(3), zp(3), xlp(3), ylp(3), zlp(3)
        double precision r1(3, 3), dll(18)
        double precision xg(3), yg(3), zg(3), str(6)
        double precision cb, x21, y21, z21, x32, y32, z32, x13, y13, 
     *z13
        double precision rx, ry, rz, bx, by, bz, rlr, rlb, bpr, 
     *area, ylr, zlr
        double precision ycg, xcg, zcg, xlcg, ylcg, zlcg, f, t
        double precision rmom(18, 3), rmem(18, 3)
        double precision rmx, rmy, rmxy, rnx, rny, rnxy, clr, cqr
        double precision rmmx, rmmy, rmmxy, rnnx, rnny, rnnxy, sbf
        double precision ebar
        integer lb(9), le(9), strainflg
        character*10 status
        integer i, j
        double precision d1
c
c manually inserted - begin	     
c
c
c manually inserted - end      
c
        double precision d1_p, d2_w, d2_v, d3_v, d4_v, d1_w, d7_b,
     * d6_b, d5_b, d2_b
        double precision d3_b, d4_b, g_rmom(18, 3), 
     *g_rmem(18, 3), g_t, g_h(3), g_cb, 
     *g_e, g_db(3, 3), g_dm(3, 3)
        double precision g_x21, g_xl(3), g_y21
     *, g_yl(3), g_z21, g_zl(3), g_x32
     *, g_y32, g_z32, g_rx
        double precision g_ry, g_rz, g_d2_w, 
     *g_d1_w, g_rlr, g_xp(3), g_zp(3)
     *, g_zlr, g_yp(3), g_ylr
        double precision g_xcg, g_ycg, g_zcg,
     * g_xlcg, g_ylcg, g_zlcg, g_xlp(3)
     *, g_ylp(3), g_dll(18), g_r1(3, 3)
        double precision g_v(18), g_rmx, g_rmy, 
     *g_rmxy, g_rnx, g_rny, g_rnxy, 
     *g_t2, g_epsxx, g_epsyy
        double precision g_epsxy, g_epszz, g_str(6)
     *, g_stress(maxsze, maxstr, maxgus), g_ebar
     *, g_rmmx, g_rmmy, g_rmmxy, g_rnnx
     *, g_rnny
        double precision g_rnnxy, g_rnxt, g_rnyt
     *, g_rnxyt, g_rmxt, g_rmyt, g_rmxyt
     *, g_sx, g_sy, g_sxy
        double precision g_sbf
        save g_sxy, g_sbf
        save g_rnny, g_rnnxy, g_rnxt, g_rnyt, g_rnxyt, g_rmxt, g_rmyt,
     * g_rmxyt, g_sx, g_sy
        save g_epsxx, g_epsyy, g_epsxy, g_epszz, g_str, g_ebar,  
     *g_rmmx,g_rmmy, g_rmmxy, g_rnnx
        save g_ylp, g_dll, g_r1, g_rmx, g_rmy, g_rmxy, g_rnx, g_rny, 
     *g_rnxy, g_t2
        save g_zlr, g_yp, g_ylr, g_xcg, g_ycg, g_zcg, g_xlcg, g_ylcg,
     * g_zlcg, g_xlp
        save g_y32, g_z32, g_rx, g_ry, g_rz, g_d2_w, g_d1_w, g_rlr, 
     *g_xp, g_zp
        save g_rmom, g_rmem, g_t, g_cb, g_db, g_dm, g_x21, g_y21, 
     *g_z21, g_x32
        external g_vonmis
        external g_straineq
        external g_transform
        external g_momen
        external g_membra
        external g_rotation
        intrinsic dble
        data lb /3, 4, 5, 9, 10, 11, 15, 16, 17/
        data le /1, 2, 7, 8, 13, 14, 6, 12, 18/
c
c	 PARAMETER(alpha = 1.50d+00)
c	 PARAMETER(beta  = 0.32d+00)
c       
c	PARAMETER(alpha = 0.0d+00)
c	PARAMETER(beta  = 0.0d+00)
c
C
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'sands8','g_sands8.f')
C
        do 99999 i = 1, 18
            g_rmom(i, 1) = 0.0d0
          rmom(i, 1) = 0.0d0
C--------
            g_rmom(i, 2) = 0.0d0
          rmom(i, 2) = 0.0d0
C--------
            g_rmom(i, 3) = 0.0d0
          rmom(i, 3) = 0.0d0
C--------
            g_rmem(i, 1) = 0.0d0
          rmem(i, 1) = 0.0d0
C--------
            g_rmem(i, 2) = 0.0d0
          rmem(i, 2) = 0.0d0
C--------
            g_rmem(i, 3) = 0.0d0
          rmem(i, 3) = 0.0d0
C--------
5         continue
99999   continue
C
        f = 1.0d+00
        clr = 0.0d+00
        cqr = 1.0d+00
          g_t = g_h(1)
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
          g_cb = d7_b * g_t + d4_b * g_e
        cb = e * d4_v / dble(12.0) / (1.0d0 - nu * nu)
C--------
          g_db(1, 1) = g_cb
        db(1, 1) = cb
C--------
          g_db(1, 2) = nu * g_cb
        db(1, 2) = nu * cb
C--------
          g_db(1, 3) = 0.0d0
        db(1, 3) = 0.0d0
C--------
          g_db(2, 1) = g_db(1, 2)
        db(2, 1) = db(1, 2)
C--------
          g_db(2, 2) = g_cb
        db(2, 2) = cb
C--------
          g_db(2, 3) = 0.0d0
        db(2, 3) = 0.0d0
C--------
          g_db(3, 1) = 0.0d0
        db(3, 1) = 0.0d0
C--------
          g_db(3, 2) = 0.0d0
        db(3, 2) = 0.0d0
C--------
        d2_b = dble(0.5) * (1.0d0 - nu)
          g_db(3, 3) = d2_b * g_cb
        db(3, 3) = dble(0.5) * (1.0d0 - nu) * cb
C--------
C
C set the membrane constitutive matrix
C
        d2_b = 1.0d0 / (1.0d0 - nu * nu)
        d3_b = d2_b * t
        d4_b = d2_b * e
          g_cb = d4_b * g_t + d3_b * g_e
        cb = e * t / (1.0d0 - nu * nu)
C--------
          g_dm(1, 1) = g_cb
        dm(1, 1) = cb
C--------
          g_dm(1, 2) = nu * g_cb
        dm(1, 2) = nu * cb
C--------
          g_dm(1, 3) = 0.0d0
        dm(1, 3) = 0.0d0
C--------
          g_dm(2, 1) = g_db(1, 2)
        dm(2, 1) = db(1, 2)
C--------
          g_dm(2, 2) = g_cb
        dm(2, 2) = cb
C--------
          g_dm(2, 3) = 0.0d0
        dm(2, 3) = 0.0d0
C--------
          g_dm(3, 1) = 0.0d0
        dm(3, 1) = 0.0d0
C--------
          g_dm(3, 2) = 0.0d0
        dm(3, 2) = 0.0d0
C--------
        d2_b = dble(0.5) * (1.0d0 - nu)
          g_dm(3, 3) = d2_b * g_cb
        dm(3, 3) = dble(0.5) * (1.0d0 - nu) * cb
C--------
C
C triangular dimension variables
C
          g_x21 = -g_xl(1) + g_xl(2)
        x21 = xl(2) - xl(1)
C--------
          g_y21 = -g_yl(1) + g_yl(2)
        y21 = yl(2) - yl(1)
C--------
          g_z21 = -g_zl(1) + g_zl(2)
        z21 = zl(2) - zl(1)
C--------
          g_x32 = -g_xl(2) + g_xl(3)
        x32 = xl(3) - xl(2)
C--------
          g_y32 = -g_yl(2) + g_yl(3)
        y32 = yl(3) - yl(2)
C--------
          g_z32 = -g_zl(2) + g_zl(3)
        z32 = zl(3) - zl(2)
C--------
        x13 = xl(1) - xl(3)
        y13 = yl(1) - yl(3)
        z13 = zl(1) - zl(3)
C  triangle in space : we compute the length of one side and the distance of the
C  opposing node to that side to compute the area
          g_rx = g_x21
        rx = x21
C--------
          g_ry = g_y21
        ry = y21
C--------
          g_rz = g_z21
        rz = z21
C--------
        bx = x32
        by = y32
        bz = z32
        d4_b = ry + ry
        d5_b = rx + rx
          g_d2_w = d4_b * g_ry + d5_b * g_rx
        d2_w = rx * rx + ry * ry
        d4_b = rz + rz
          g_d1_w = d4_b * g_rz + g_d2_w
        d1_w = d2_w + rz * rz
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +368)
        endif
          g_rlr = d1_p * g_d1_w
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
          g_xp(1) = d3_b * g_rlr + d2_b * g_x21
        xp(1) = d3_v
C--------
        d3_v = y21 / rlr
        d2_b = 1.0d0 / rlr
        d3_b = -d3_v / rlr
          g_xp(2) = d3_b * g_rlr + d2_b * g_y21
        xp(2) = d3_v
C--------
        d3_v = z21 / rlr
        d2_b = 1.0d0 / rlr
        d3_b = -d3_v / rlr
          g_xp(3) = d3_b * g_rlr + d2_b * g_z21
        xp(3) = d3_v
C--------
C
C Z' local axis
C
          g_zp(1) = -z21 * g_y32 + (-y32) * g_z21 + 
     *y21 * g_z32 + z32 * g_y21
        zp(1) = y21 * z32 - z21 * y32
C--------
          g_zp(2) = -x21 * g_z32 + (-z32) * g_x21 + 
     *z21 * g_x32 + x32 * g_z21
        zp(2) = z21 * x32 - x21 * z32
C--------
          g_zp(3) = -y21 * g_x32 + (-x32) * g_y21 + 
     *x21 * g_y32 + y32 * g_x21
        zp(3) = x21 * y32 - y21 * x32
C--------
        d4_b = zp(2) + zp(2)
        d5_b = zp(1) + zp(1)
          g_d2_w = d4_b * g_zp(2) + d5_b * g_zp(1)
        d2_w = zp(1) * zp(1) + zp(2) * zp(2)
        d4_b = zp(3) + zp(3)
          g_d1_w = d4_b * g_zp(3) + g_d2_w
        d1_w = d2_w + zp(3) * zp(3)
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +445)
        endif
          g_zlr = d1_p * g_d1_w
        zlr = d2_v
C--------
        d3_v = zp(1) / zlr
        d2_b = 1.0d0 / zlr
        d3_b = -d3_v / zlr
          g_zp(1) = d3_b * g_zlr + d2_b * g_zp(1)
        zp(1) = d3_v
C--------
        d3_v = zp(2) / zlr
        d2_b = 1.0d0 / zlr
        d3_b = -d3_v / zlr
          g_zp(2) = d3_b * g_zlr + d2_b * g_zp(2)
        zp(2) = d3_v
C--------
        d3_v = zp(3) / zlr
        d2_b = 1.0d0 / zlr
        d3_b = -d3_v / zlr
          g_zp(3) = d3_b * g_zlr + d2_b * g_zp(3)
        zp(3) = d3_v
C--------
C
C Y' local axis
C
          g_yp(1) = -zp(3) * g_xp(2) + (-xp(2)) * 
     *g_zp(3) + zp(2) * g_xp(3) + xp(3) * g_zp(2)
        yp(1) = zp(2) * xp(3) - zp(3) * xp(2)
C--------
          g_yp(2) = -zp(1) * g_xp(3) + (-xp(3)) * 
     *g_zp(1) + zp(3) * g_xp(1) + xp(1) * g_zp(3)
        yp(2) = zp(3) * xp(1) - zp(1) * xp(3)
C--------
          g_yp(3) = -zp(2) * g_xp(1) + (-xp(1)) * 
     *g_zp(2) + zp(1) * g_xp(2) + xp(2) * g_zp(1)
        yp(3) = zp(1) * xp(2) - zp(2) * xp(1)
C--------
        d4_b = yp(2) + yp(2)
        d5_b = yp(1) + yp(1)
          g_d2_w = d4_b * g_yp(2) + d5_b * g_yp(1)
        d2_w = yp(1) * yp(1) + yp(2) * yp(2)
        d4_b = yp(3) + yp(3)
          g_d1_w = d4_b * g_yp(3) + g_d2_w
        d1_w = d2_w + yp(3) * yp(3)
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +514)
        endif
          g_ylr = d1_p * g_d1_w
        ylr = d2_v
C--------
        d3_v = yp(1) / ylr
        d2_b = 1.0d0 / ylr
        d3_b = -d3_v / ylr
          g_yp(1) = d3_b * g_ylr + d2_b * g_yp(1)
        yp(1) = d3_v
C--------
        d3_v = yp(2) / ylr
        d2_b = 1.0d0 / ylr
        d3_b = -d3_v / ylr
          g_yp(2) = d3_b * g_ylr + d2_b * g_yp(2)
        yp(2) = d3_v
C--------
        d3_v = yp(3) / ylr
        d2_b = 1.0d0 / ylr
        d3_b = -d3_v / ylr
          g_yp(3) = d3_b * g_ylr + d2_b * g_yp(3)
        yp(3) = d3_v
C--------
C
C center of gravity
C
        d2_b = 1.0d0 / 3.0d+00
          g_xcg = d2_b * g_xl(3) + d2_b * g_xl(2) + 
     *d2_b * g_xl(1)
        xcg = (xl(1) + xl(2) + xl(3)) / 3.0d+00
C--------
        d2_b = 1.0d0 / 3.0d+00
          g_ycg = d2_b * g_yl(3) + d2_b * g_yl(2) + 
     *d2_b * g_yl(1)
        ycg = (yl(1) + yl(2) + yl(3)) / 3.0d+00
C--------
        d2_b = 1.0d0 / 3.0d+00
          g_zcg = d2_b * g_zl(3) + d2_b * g_zl(2) + 
     *d2_b * g_zl(1)
        zcg = (zl(1) + zl(2) + zl(3)) / 3.0d+00
C--------
C
C computing local coordinates
C
        do 99998 i = 1, 3
            g_xlcg = -g_xcg + g_xl(i)
          xlcg = xl(i) - xcg
C--------
            g_ylcg = -g_ycg + g_yl(i)
          ylcg = yl(i) - ycg
C--------
            g_zlcg = -g_zcg + g_zl(i)
          zlcg = zl(i) - zcg
C--------
            g_d1_w = xp(2) * g_ylcg + ylcg * g_xp(2) +
     * xp(1) * g_xlcg + xlcg * g_xp(1)
          d1_w = xp(1) * xlcg + xp(2) * ylcg
            g_xlp(i) = xp(3) * g_zlcg + zlcg * g_xp(3)
     * + g_d1_w
          xlp(i) = d1_w + xp(3) * zlcg
C--------
            g_d1_w = yp(2) * g_ylcg + ylcg * g_yp(2) +
     * yp(1) * g_xlcg + xlcg * g_yp(1)
          d1_w = yp(1) * xlcg + yp(2) * ylcg
            g_ylp(i) = yp(3) * g_zlcg + zlcg * g_yp(3)
     * + g_d1_w
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
        call g_rotation(xp, g_xp, yp, g_yp, zp, 
     *g_zp, xg, yg, zg, r1, g_r1)
C
C compute the von mises stress
C rotate nodal displacements to local system
C
        do 99997 i = 1, 18
            g_dll(i) = 0.0d0
          dll(i) = 0.0d+00
C--------
270       continue
99997   continue
        do 99995 i = 1, 3
          do 99996 j = 1, 3
              g_dll(i) = r1(j, i) * g_v( j) + v(j) * 
     *g_r1(j, i) + g_dll(i)
            dll(i) = dll(i) + r1(j, i) * v(j)
C--------
              g_dll(i + 3) = r1(j, i) * g_v( j + 3) + 
     *v(j + 3) * g_r1(j, i) + g_dll(i + 3)
            dll(i + 3) = dll(i + 3) + r1(j, i) * v(j + 3)
C--------
              g_dll(i + 6) = r1(j, i) * g_v( j + 6) + v(j + 6
     *) * g_r1(j, i) + g_dll(i + 6)
            dll(i + 6) = dll(i + 6) + r1(j, i) * v(j + 6)
C--------
              g_dll(i + 9) = r1(j, i) * g_v( j + 9) + v(j + 9
     *) * g_r1(j, i) + g_dll(i + 9)
            dll(i + 9) = dll(i + 9) + r1(j, i) * v(j + 9)
C--------
              g_dll(i + 12) = r1(j, i) * g_v( j + 12) + v(j +
     * 12) * g_r1(j, i) + g_dll(i + 12)
            dll(i + 12) = dll(i + 12) + r1(j, i) * v(j + 12)
C--------
              g_dll(i + 15) = r1(j, i) * g_v( j + 15) + v(j +
     * 15) * g_r1(j, i) + g_dll(i + 15)
            dll(i + 15) = dll(i + 15) + r1(j, i) * v(j + 15)
C--------
280         continue
99996     continue
99995   continue
C
C      compute centroidal membrane strains
        d1 = alpha
        call g_membra(xlp, g_xlp, ylp, g_ylp, d1
     *, le, rmem, g_rmem, status)
C
C      compute centroidal bending strains (curvatures (1/radius))
        call g_momen(xlp, g_xlp, ylp, g_ylp, lb,
     * rmom, g_rmom, status)
C
          g_rmx = 0.0d0
        rmx = 0.0d0
C--------
          g_rmy = 0.0d0
        rmy = 0.0d0
C--------
          g_rmxy = 0.0d0
        rmxy = 0.0d0
C--------
          g_rnx = 0.0d0
        rnx = 0.0d0
C--------
          g_rny = 0.0d0
        rny = 0.0d0
C--------
          g_rnxy = 0.0d0
        rnxy = 0.0d0
C--------
        do 99994 j = 1, 18
            g_rmx = rmom(j, 1) * g_dll(j) + dll(j) * 
     *g_rmom(j, 1) + g_rmx
          rmx = rmx + rmom(j, 1) * dll(j)
C--------
            g_rmy = rmom(j, 2) * g_dll(j) + dll(j) * 
     *g_rmom(j, 2) + g_rmy   
          rmy = rmy + rmom(j, 2) * dll(j)
C--------
            g_rmxy = rmom(j, 3) * g_dll(j) + dll(j) * 
     *g_rmom(j, 3) + g_rmxy
          rmxy = rmxy + rmom(j, 3) * dll(j)
C--------
            g_rnx = rmem(j, 1) * g_dll(j) + dll(j) * 
     *g_rmem(j, 1) + g_rnx
          rnx = rnx + rmem(j, 1) * dll(j)
C--------
            g_rny = rmem(j, 2) * g_dll(j) + dll(j) * 
     *g_rmem(j, 2) + g_rny
          rny = rny + rmem(j, 2) * dll(j)
C--------
            g_rnxy = rmem(j, 3) * g_dll(j) + dll(j) * 
     *g_rmem(j, 3) + g_rnxy
          rnxy = rnxy + rmem(j, 3) * dll(j)
C--------
290       continue
99994   continue
C
C      subtract off thermal strain
       rnx = rnx - thrmStr
       rny = rny - thrmStr
C
       g_rnx = g_rnx - dthrmStr
       g_rny = g_rny - dthrmStr
C
C      COMPUTE VON MISES STRAIN RESULTANT
C
        if (strainflg .eq. 1) then
C
          d2_b = dble(0.5)
            g_t2 = d2_b * g_t
          t2 = dble(0.5) * t
C--------
C
            g_rmx = t2 * g_rmx + rmx * g_t2
          rmx = t2 * rmx
C--------
            g_rmy = t2 * g_rmy + rmy * g_t2
          rmy = t2 * rmy
C--------
            g_rmxy = t2 * g_rmxy + rmxy * g_t2
          rmxy = t2 * rmxy
C--------
C
C ... COMPUTE STRAINS AT MEDIAN SURFACE
          if (surface .eq. 2) then
              g_epsxx = g_rnx
            epsxx = rnx
C--------
              g_epsyy = g_rny
            epsyy = rny
C--------
            d2_b = dble(0.5)
              g_epsxy = d2_b * g_rnxy
            epsxy = dble(0.5) * rnxy
C--------
C ... COMPUTE STRAINS AT BOTTOM SURFACE
          else
            if (surface .eq. 3) then
                g_epsxx = -g_rmx + g_rnx
              epsxx = rnx - rmx
C--------
                g_epsyy = -g_rmy + g_rny
              epsyy = rny - rmy
C--------
              d2_b = dble(0.5)
                g_epsxy = -d2_b * g_rmxy + d2_b * 
     *g_rnxy
              epsxy = dble(0.5) * (rnxy - rmxy)
C--------
C ... COMPUTE STRAINS AT TOP SURFACE
            else
                g_epsxx = g_rmx + g_rnx
              epsxx = rnx + rmx
C--------
                g_epsyy = g_rmy + g_rny
              epsyy = rny + rmy
C--------
              d2_b = dble(0.5)
                g_epsxy = d2_b * g_rmxy + d2_b * 
     *g_rnxy
              epsxy = dble(0.5) * (rnxy + rmxy)
C--------
            endif
          endif
C
          d2_b = -nu / (1.0d0 - nu)
            g_epszz = d2_b * g_epsyy + d2_b * g_epsxx
          epszz = -nu / (1.0d0 - nu) * (epsxx + epsyy)
C--------
C
            g_str(1) = g_epsxx
          str(1) = epsxx
C--------
            g_str(2) = g_epsyy
          str(2) = epsyy
C--------
            g_str(3) = g_epszz
          str(3) = epszz
C--------
            g_str(4) = g_epsxy
          str(4) = epsxy
C--------
            g_str(5) = 0.0d0
          str(5) = 0.0d0
C--------
            g_str(6) = 0.0d0
          str(6) = 0.0d0
C--------
C
          call g_transform( xp, g_xp, yp, g_yp, z
     *p, g_zp, xg, yg, zg, str, g_str)
C
            g_stress(elm, 1, 1) = g_str(1)
          stress(elm, 1, 1) = str(1)
C--------
            g_stress(elm, 1, 2) = g_str(1)
          stress(elm, 1, 2) = str(1)
C--------
            g_stress(elm, 1, 3) = g_str(1)
          stress(elm, 1, 3) = str(1)
C--------
C
            g_stress(elm, 2, 1) = g_str(2)
          stress(elm, 2, 1) = str(2)
C--------
            g_stress(elm, 2, 2) = g_str(2)
          stress(elm, 2, 2) = str(2)
C--------
            g_stress(elm, 2, 3) = g_str(2)
          stress(elm, 2, 3) = str(2)
C--------
C
            g_stress(elm, 3, 1) = g_str(3)
          stress(elm, 3, 1) = str(3)
C--------
            g_stress(elm, 3, 2) = g_str(3)
          stress(elm, 3, 2) = str(3)
C--------
            g_stress(elm, 3, 3) = g_str(3)
          stress(elm, 3, 3) = str(3)
C--------
C
            g_stress(elm, 4, 1) = 2.0D0*g_str(4)
          stress(elm, 4, 1) = 2.0D0*str(4)
C--------
            g_stress(elm, 4, 2) = 2.0D0*g_str(4)
          stress(elm, 4, 2) = 2.0D0*str(4)
C--------
            g_stress(elm, 4, 3) = 2.0D0*g_str(4)
          stress(elm, 4, 3) = 2.0D0*str(4)
C--------
C
            g_stress(elm, 5, 1) = 2.0D0*g_str(5)
          stress(elm, 5, 1) = 2.0D0*str(5)
C--------
            g_stress(elm, 5, 2) = 2.0D0*g_str(5)
          stress(elm, 5, 2) = 2.0D0*str(5)
C--------
            g_stress(elm, 5, 3) = 2.0D0*g_str(5)
          stress(elm, 5, 3) = 2.0D0*str(5)
C--------
C
            g_stress(elm, 6, 1) = 2.0D0*g_str(6)
          stress(elm, 6, 1) = 2.0D0*str(6)
C--------
            g_stress(elm, 6, 2) = 2.0D0*g_str(6)
          stress(elm, 6, 2) = 2.0D0*str(6)
C--------
            g_stress(elm, 6, 3) = 2.0D0*g_str(6)
          stress(elm, 6, 3) = 2.0D0*str(6)
C--------
C
          call g_straineq(rmx, g_rmx, rmy, g_rmy
     *, rmxy, g_rmxy, rnx, g_rnx, rny, g_rny,
     * rnxy, g_rnxy, t, g_t, surface, ebar, g_ebar 
     *)
C
            g_stress(elm, 7, 1) = g_ebar
          stress(elm, 7, 1) = ebar
C--------
            g_stress(elm, 7, 2) = g_ebar
          stress(elm, 7, 2) = ebar
C--------
            g_stress(elm, 7, 3) = g_ebar
          stress(elm, 7, 3) = ebar
C--------
C
          return
        endif
C
C     compute centroidal stress resultants
C
C     bending resultants
          g_d1_w = db(1, 2) * g_rmy + rmy * g_db(1, 2)
     * + db(1, 1) * g_rmx + rmx * g_db(1, 1)
        d1_w = db(1, 1) * rmx + db(1, 2) * rmy   
          g_rmmx = db(1, 3) * g_rmxy + rmxy *  
     *g_db(1, 3) + g_d1_w
        rmmx = d1_w + db(1, 3) * rmxy
C--------
          g_d1_w = db(2, 2) * g_rmy + rmy * g_db(2, 2)
     * + db(2, 1) * g_rmx + rmx * g_db(2, 1)
        d1_w = db(2, 1) * rmx + db(2, 2) * rmy
          g_rmmy = db(2, 3) * g_rmxy + rmxy *
     * g_db(2, 3) + g_d1_w
        rmmy = d1_w + db(2, 3) * rmxy
C--------
          g_d1_w = db(3, 2) * g_rmy + rmy * g_db(3, 2)
     * + db(3, 1) * g_rmx + rmx * g_db(3, 1)
        d1_w = db(3, 1) * rmx + db(3, 2) * rmy
          g_rmmxy = db(3, 3) * g_rmxy + rmxy * 
     *g_db(3, 3) + g_d1_w
        rmmxy = d1_w + db(3, 3) * rmxy
C--------
C
C     membrane resultants
          g_d1_w = dm(1, 2) * g_rny + rny * g_dm(1, 2)
     * + dm(1, 1) * g_rnx + rnx * g_dm(1, 1)
        d1_w = dm(1, 1) * rnx + dm(1, 2) * rny
          g_rnnx = dm(1, 3) * g_rnxy + rnxy * 
     *g_dm(1, 3) + g_d1_w
        rnnx = d1_w + dm(1, 3) * rnxy
C--------
          g_d1_w = dm(2, 2) * g_rny + rny * g_dm(2, 2)
     * + dm(2, 1) * g_rnx + rnx * g_dm(2, 1)
        d1_w = dm(2, 1) * rnx + dm(2, 2) * rny
          g_rnny = dm(2, 3) * g_rnxy + rnxy * 
     *g_dm(2, 3) + g_d1_w
        rnny = d1_w + dm(2, 3) * rnxy
C--------
          g_d1_w = dm(3, 2) * g_rny + rny * g_dm(3, 2)
     * + dm(3, 1) * g_rnx + rnx * g_dm(3, 1)
        d1_w = dm(3, 1) * rnx + dm(3, 2) * rny
          g_rnnxy = dm(3, 3) * g_rnxy + rnxy * 
     *g_dm(3, 3) + g_d1_w
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
          g_t = d1_p * g_t
        t = d2_v
C--------
        d2_b = t + t
          g_t2 = d2_b * g_t
        t2 = t * t
C--------
C
        d3_v = rnnx / t
        d2_b = 1.0d0 / t
        d3_b = -d3_v / t
          g_rnxt = d3_b * g_t + d2_b * g_rnnx
        rnxt = d3_v
C--------
        d3_v = rnny / t
        d2_b = 1.0d0 / t
        d3_b = -d3_v / t
          g_rnyt = d3_b * g_t + d2_b * g_rnny
        rnyt = d3_v
C--------
        d3_v = rnnxy / t
        d2_b = 1.0d0 / t
        d3_b = -d3_v / t
          g_rnxyt = d3_b * g_t + d2_b * g_rnnxy
        rnxyt = d3_v
C--------
C
        d4_v = dble(6.0) * rmmx / t2
        d3_b = -d4_v / t2
        d4_b = 1.0d0 / t2 * dble(6.0)
          g_rmxt = d3_b * g_t2 + d4_b * g_rmmx
        rmxt = d4_v
C--------
        d4_v = dble(6.0) * rmmy / t2
        d3_b = -d4_v / t2
        d4_b = 1.0d0 / t2 * dble(6.0)
          g_rmyt = d3_b * g_t2 + d4_b * g_rmmy
        rmyt = d4_v
C--------
        d4_v = dble(6.0) * rmmxy / t2
        d3_b = -d4_v / t2
        d4_b = 1.0d0 / t2 * dble(6.0)
          g_rmxyt = d3_b * g_t2 + d4_b * g_rmmxy
        rmxyt = d4_v
C--------
C
C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT TOP SURFACE (DEFAULT)
          g_sx = g_rmxt + g_rnxt
        sx = rnxt + rmxt
C--------
          g_sy = g_rmyt + g_rnyt
        sy = rnyt + rmyt
C--------
          g_sxy = g_rmxyt + g_rnxyt
        sxy = rnxyt + rmxyt
C--------
C
C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT MEDIAN SURFACE
        if (surface .eq. 2) then
            g_sx = g_rnxt
          sx = rnxt
C--------
            g_sy = g_rnyt
          sy = rnyt
C--------
            g_sxy = g_rnxyt
          sxy = rnxyt
C--------
        endif
C
C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT BOTTOM SURFACE
        if (surface .eq. 3) then
            g_sx = -g_rmxt + g_rnxt
          sx = rnxt - rmxt
C--------
            g_sy = -g_rmyt + g_rnyt
          sy = rnyt - rmyt
C--------
            g_sxy = -g_rmxyt + g_rnxyt
          sxy = rnxyt - rmxyt
C--------
        endif
C
          g_str(1) = g_sx
        str(1) = sx
C--------
          g_str(2) = g_sy
        str(2) = sy
C--------
          g_str(3) = 0.0d0
        str(3) = 0.0d0
C--------
          g_str(4) = g_sxy
        str(4) = sxy
C--------
          g_str(5) = 0.0d0
        str(5) = 0.0d0
C--------
          g_str(6) = 0.0d0
        str(6) = 0.0d0
C--------
C
        call g_transform(xp, g_xp, yp, g_yp, zp,
     * g_zp, xg, yg, zg, str, g_str)
C
          g_stress(elm, 1, 1) = g_str(1)
        stress(elm, 1, 1) = str(1)
C--------
          g_stress(elm, 1, 2) = g_str(1)
        stress(elm, 1, 2) = str(1)
C--------
          g_stress(elm, 1, 3) = g_str(1)
        stress(elm, 1, 3) = str(1)
C--------
C
          g_stress(elm, 2, 1) = g_str(2)
        stress(elm, 2, 1) = str(2)
C--------
          g_stress(elm, 2, 2) = g_str(2)
        stress(elm, 2, 2) = str(2)
C--------
          g_stress(elm, 2, 3) = g_str(2)
        stress(elm, 2, 3) = str(2)
C--------
C
          g_stress(elm, 3, 1) = g_str(3)
        stress(elm, 3, 1) = str(3)
C--------
          g_stress(elm, 3, 2) = g_str(3)
        stress(elm, 3, 2) = str(3)
C--------
          g_stress(elm, 3, 3) = g_str(3)
        stress(elm, 3, 3) = str(3)
C--------
C
          g_stress(elm, 4, 1) = g_str(4)
        stress(elm, 4, 1) = str(4)
C--------
          g_stress(elm, 4, 2) = g_str(4)
        stress(elm, 4, 2) = str(4)
C--------
          g_stress(elm, 4, 3) = g_str(4)
        stress(elm, 4, 3) = str(4)
C--------
C
          g_stress(elm, 5, 1) = g_str(5)
        stress(elm, 5, 1) = str(5)
C--------
          g_stress(elm, 5, 2) = g_str(5)
        stress(elm, 5, 2) = str(5)
C--------
          g_stress(elm, 5, 3) = g_str(5)
        stress(elm, 5, 3) = str(5)
C--------
C
          g_stress(elm, 6, 1) = g_str(6)
        stress(elm, 6, 1) = str(6)
C--------
          g_stress(elm, 6, 2) = g_str(6)
        stress(elm, 6, 2) = str(6)
C--------
          g_stress(elm, 6, 3) = g_str(6)
        stress(elm, 6, 3) = str(6)
C--------
C
C     compute von mises stress resultant
C
        call g_vonmis(rmmx, g_rmmx, rmmy, g_rmmy
     *, rmmxy, g_rmmxy, rnnx, g_rnnx, rnny, g_rnny, 
     *rnnxy, g_rnnxy, t, g_t, sbf, g_sbf, 
     *surface)
C
          g_stress(elm, 7, 1) = g_sbf
        stress(elm, 7, 1) = sbf
C--------
          g_stress(elm, 7, 2) = g_sbf
        stress(elm, 7, 2) = sbf
C--------
          g_stress(elm, 7, 3) = g_sbf
        stress(elm, 7, 3) = sbf
C--------
C
        return
      end
