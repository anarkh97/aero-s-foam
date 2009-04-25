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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This subroutine computes the mass matrix for the
C AQR element, shell version. The main idea is to distri-
C bute the complete mass of the element over three
C fictitious beams allocated over the sides. The big
C problem is that we must construct the consistent mass
C matrix in the global system and then lump it.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
      subroutine gxmass8(xl,g_xl,yl,g_yl,zl,g_zl,h,g_h,dens,g_dens,
     *                   fv,g_fv,eltdof,gamma,
     *                    grvfor,grvflg,totmas,masflg)
c
c
c
        integer i, j, i1, i2, i3, eltdof, nd
        real*8 xl(3), yl(3), zl(3), h(3)
        real*8 fv(eltdof, eltdof), rl(3)
        real*8 x21, y21, z21, x32, y32, z32, x13, y13, z13
        real*8 rx, ry, rz, bx, by, bz, rlr, rlb, bpr, area, rmas
        real*8 dens, esp, totmas,coef
        real*8 gamma(*), grvfor(*)
        logical grvflg, masflg
c
c manually inserted - begin	     
c
        
c
c manually inserted - end      
c
        double precision d2_w, d3_p, d2_p, d1_p, d8_b, d7_b, d6_b, d5_b,
     * d2_v, d3_v
        double precision d4_b, d2_b, d3_b, d1_w, d4_v, d7_v, d6_v, g_esp
     *, g_h(3), g_fv(eltdof, eltdof)
        double precision g_x21, g_xl(3), g_y21
     *, g_yl(3), g_z21, g_zl(3), g_x32
     *, g_y32, g_z32, g_x13
        double precision g_y13, g_z13, g_d1_w
     *, g_rl(3), g_rx, g_ry, g_rz, g_bx
     *, g_by, g_bz
        double precision g_d2_w, g_rlr, g_rlb
     *, g_bpr, g_area, g_rmas, g_dens
     *
        save g_rlb, g_bpr, g_area, g_rmas
        save g_d1_w, g_rl, g_rx, g_ry, g_rz, g_bx, g_by, g_bz, g_d2_w, 
     *g_rlr
        save g_esp, g_x21, g_y21, g_z21, g_x32, g_y32, g_z32, g_x13, 
     *g_y13, g_z13
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'mass8','g_mass8.f')
C
        nd = 18
C initialize mass matrix to zero
          g_esp = g_h(1)
        esp = h(1)
C--------
        do 99998 i = 1, eltdof
          do 99999 j = 1, eltdof
              g_fv(i, j) = 0.0d0
            fv(i, j) = 0.0d0
C--------
10          continue
99999     continue
99998   continue
C dimension variables
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
          g_x13 = -g_xl(3) + g_xl(1)
        x13 = xl(1) - xl(3)
C--------
          g_y13 = -g_yl(3) + g_yl(1)
        y13 = yl(1) - yl(3)
C--------
          g_z13 = -g_zl(3) + g_zl(1)
        z13 = zl(1) - zl(3)
C--------
        d2_v = x21 * x21
        d3_p = 2.0d0 * x21
        d4_v = y21 * y21
        d2_p = 2.0d0 * y21
        d7_v = z21 * z21
        d1_p = 2.0d0 * z21
          g_d1_w = d1_p * g_z21 + d2_p * g_y21 + d3_p 
     ** g_x21
        d1_w = d2_v + d4_v + d7_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +148)
        endif
          g_rl(1) = d1_p * g_d1_w
        rl(1) = d2_v
C--------
        d2_v = x32 * x32
        d3_p = 2.0d0 * x32
        d4_v = y32 * y32
        d2_p = 2.0d0 * y32
        d7_v = z32 * z32
        d1_p = 2.0d0 * z32
          g_d1_w = d1_p * g_z32 + d2_p * g_y32 + d3_p 
     ** g_x32
        d1_w = d2_v + d4_v + d7_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +172)
        endif
          g_rl(2) = d1_p * g_d1_w
        rl(2) = d2_v
C--------
        d2_v = x13 * x13
        d3_p = 2.0d0 * x13
        d4_v = y13 * y13
        d2_p = 2.0d0 * y13
        d7_v = z13 * z13
        d1_p = 2.0d0 * z13
          g_d1_w = d1_p * g_z13 + d2_p * g_y13 + d3_p 
     ** g_x13
        d1_w = d2_v + d4_v + d7_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +196)
        endif
          g_rl(3) = d1_p * g_d1_w
        rl(3) = d2_v
C--------
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
          g_bx = g_x32
        bx = x32
C--------
          g_by = g_y32
        by = y32
C--------
          g_bz = g_z32
        bz = z32
C--------
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
     +252)
        endif
          g_rlr = d1_p * g_d1_w
        rlr = d2_v
C--------
        d4_b = by + by
        d5_b = bx + bx
          g_d2_w = d4_b * g_by + d5_b * g_bx
        d2_w = bx * bx + by * by
        d4_b = bz + bz
          g_d1_w = d4_b * g_bz + g_d2_w
        d1_w = d2_w + bz * bz
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +276)
        endif
          g_rlb = d1_p * g_d1_w
        rlb = d2_v
C--------
          g_d2_w = ry * g_by + by * g_ry + rx * 
     *g_bx + bx * g_rx
        d2_w = rx * bx + ry * by
        d6_v = (d2_w + rz * bz) * (d2_w + rz * bz)
        d1_p = 2.0d0 * (d2_w + rz * bz)
        d5_b = d1_p * bz
        d6_b = d1_p * rz
          g_d1_w = d6_b * g_bz + d5_b * g_rz + d1_p * 
     *g_d2_w
        d1_w = d6_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +303)
        endif
        d4_v = d2_v / rlr
        d3_b = -d4_v / rlr
        d4_b = 1.0d0 / rlr * d1_p
          g_bpr = d3_b * g_rlr + d4_b * g_d1_w
        bpr = d4_v
C--------
        d2_v = rlb * rlb
        d2_p = 2.0d0 * rlb
        d4_v = bpr * bpr
        d1_p = 2.0d0 * bpr
          g_d1_w = -d1_p * g_bpr + d2_p * g_rlb
        d1_w = d2_v - d4_v
        d3_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d3_v)
        else
           call ehufDO (9,d1_w, d3_v, d1_p,
     +g_ehfid,
     +327)
        endif
        d2_b = 1.0d0 / 2.0d+00
        d3_b = d2_b * d3_v
        d5_b = d2_b * rlr * d1_p
          g_area = d5_b * g_d1_w + d3_b * g_rlr
        area = rlr * d3_v / 2.0d+00
C--------
        d3_v = dens * area
        d2_b = 1.0d0 / dble(3.0)
        d3_b = d2_b * esp
        d4_b = d2_b * d3_v
        d5_b = d3_b * area
        d6_b = d3_b * dens
          g_rmas = d4_b * g_esp + d6_b * g_area + d5_b
     * * g_dens
        rmas = d3_v * esp / dble(3.0)
C--------
C     
C constructing mass matrix
        do 99997 i = 1, 3
          i2 = 6 + i
          i3 = 12 + i
            g_fv(i, i) = g_rmas
          fv(i, i) = rmas
C--------
            g_fv(i2, i2) = g_rmas
          fv(i2, i2) = rmas
C--------
            g_fv(i3, i3) = g_rmas
          fv(i3, i3) = rmas
C--------
20        continue
99997   continue
        do 99996 i = 1, 3
          i1 = i + 3
          i2 = i + 9
          i3 = i + 15
          d2_v = rl(1) * rl(1)
          d2_p = 2.0d0 * rl(1)
          d4_v = rl(3) * rl(3)
          d1_p = 2.0d0 * rl(3)
          d6_v = (d2_v + d4_v) / dble(420.)
          d4_b = rmas * (1.0d0 / dble(420.))
          d7_b = d4_b * d1_p
          d8_b = d4_b * d2_p
            g_fv(i1, i1) = d6_v * g_rmas + d7_b * 
     *g_rl(3) + d8_b * g_rl(1)
          fv(i1, i1) = d6_v * rmas
C--------
          d2_v = rl(2) * rl(2)
          d2_p = 2.0d0 * rl(2)
          d4_v = rl(1) * rl(1)
          d1_p = 2.0d0 * rl(1)
          d6_v = (d2_v + d4_v) / dble(420.)
          d4_b = rmas * (1.0d0 / dble(420.))
          d7_b = d4_b * d1_p
          d8_b = d4_b * d2_p
            g_fv(i2, i2) = d6_v * g_rmas + d7_b * 
     *g_rl(1) + d8_b * g_rl(2)
          fv(i2, i2) = d6_v * rmas
C--------
          d2_v = rl(3) * rl(3)
          d2_p = 2.0d0 * rl(3)
          d4_v = rl(2) * rl(2)
          d1_p = 2.0d0 * rl(2)
          d6_v = (d2_v + d4_v) / dble(420.)
          d4_b = rmas * (1.0d0 / dble(420.))
          d7_b = d4_b * d1_p
          d8_b = d4_b * d2_p
            g_fv(i3, i3) = d6_v * g_rmas + d7_b * 
     *g_rl(2) + d8_b * g_rl(3)
          fv(i3, i3) = d6_v * rmas
C--------
30        continue
99996   continue
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
        if (grvflg) then
          coef = 3.0d0 * rmas
          grvfor(1) = coef * gamma(1)
          grvfor(2) = coef * gamma(2)
          grvfor(3) = coef * gamma(3)
        endif
C
C
C.... ACCUMULATE THE SUBDOMAIN MASS
C
        if (masflg) then
          totmas = totmas + 3.0d0 * rmas
        endif
C
        return
      end
