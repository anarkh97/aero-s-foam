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
        real*8 dens, esp, totmas
        real*8 gamma(*), grvfor(*)
        logical grvflg, masflg
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_h, ldg_fv, ldg_xl, ldg_yl, ldg_zl, ldg_d
     *ens
c
c manually inserted - begin	     
c
        parameter (g_p_=1,ldg_xl=1,ldg_yl=1,ldg_zl=1,ldg_h=1,ldg_dens=1,
     *             ldg_fv=1)
c
c manually inserted - end      
c
        double precision d2_w, d3_p, d2_p, d1_p, d8_b, d7_b, d6_b, d5_b,
     * d2_v, d3_v
        double precision d4_b, d2_b, d3_b, d1_w, d4_v, d7_v, d6_v, g_esp
     *(g_pmax_), g_h(ldg_h, 3), g_fv(ldg_fv, eltdof, eltdof)
        double precision g_x21(g_pmax_), g_xl(ldg_xl, 3), g_y21(g_pmax_)
     *, g_yl(ldg_yl, 3), g_z21(g_pmax_), g_zl(ldg_zl, 3), g_x32(g_pmax_)
     *, g_y32(g_pmax_), g_z32(g_pmax_), g_x13(g_pmax_)
        double precision g_y13(g_pmax_), g_z13(g_pmax_), g_d1_w(g_pmax_)
     *, g_rl(g_pmax_, 3), g_rx(g_pmax_), g_ry(g_pmax_), g_rz(g_pmax_), g
     *_bx(g_pmax_), g_by(g_pmax_), g_bz(g_pmax_)
        double precision g_d2_w(g_pmax_), g_rlr(g_pmax_), g_rlb(g_pmax_)
     *, g_bpr(g_pmax_), g_area(g_pmax_), g_rmas(g_pmax_), g_dens(ldg_den
     *s)
        save g_rlb, g_bpr, g_area, g_rmas
        save g_d1_w, g_rl, g_rx, g_ry, g_rz, g_bx, g_by, g_bz, g_d2_w, g
     *_rlr
        save g_esp, g_x21, g_y21, g_z21, g_x32, g_y32, g_z32, g_x13, g_y
     *13, g_z13
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'mass8','g_mass8.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        nd = 18
C initialize mass matrix to zero
        do g_i_ = 1, g_p_
          g_esp(g_i_) = g_h(g_i_, 1)
        enddo
        esp = h(1)
C--------
        do 99998 i = 1, eltdof
          do 99999 j = 1, eltdof
            do g_i_ = 1, g_p_
              g_fv(g_i_, i, j) = 0.0d0
            enddo
            fv(i, j) = 0.0d0
C--------
10          continue
99999     continue
99998   continue
C dimension variables
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
        do g_i_ = 1, g_p_
          g_x13(g_i_) = -g_xl(g_i_, 3) + g_xl(g_i_, 1)
        enddo
        x13 = xl(1) - xl(3)
C--------
        do g_i_ = 1, g_p_
          g_y13(g_i_) = -g_yl(g_i_, 3) + g_yl(g_i_, 1)
        enddo
        y13 = yl(1) - yl(3)
C--------
        do g_i_ = 1, g_p_
          g_z13(g_i_) = -g_zl(g_i_, 3) + g_zl(g_i_, 1)
        enddo
        z13 = zl(1) - zl(3)
C--------
        d2_v = x21 * x21
        d3_p = 2.0d0 * x21
        d4_v = y21 * y21
        d2_p = 2.0d0 * y21
        d7_v = z21 * z21
        d1_p = 2.0d0 * z21
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d1_p * g_z21(g_i_) + d2_p * g_y21(g_i_) + d3_p 
     ** g_x21(g_i_)
        enddo
        d1_w = d2_v + d4_v + d7_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +148)
        endif
        do g_i_ = 1, g_p_
          g_rl(g_i_, 1) = d1_p * g_d1_w(g_i_)
        enddo
        rl(1) = d2_v
C--------
        d2_v = x32 * x32
        d3_p = 2.0d0 * x32
        d4_v = y32 * y32
        d2_p = 2.0d0 * y32
        d7_v = z32 * z32
        d1_p = 2.0d0 * z32
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d1_p * g_z32(g_i_) + d2_p * g_y32(g_i_) + d3_p 
     ** g_x32(g_i_)
        enddo
        d1_w = d2_v + d4_v + d7_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +172)
        endif
        do g_i_ = 1, g_p_
          g_rl(g_i_, 2) = d1_p * g_d1_w(g_i_)
        enddo
        rl(2) = d2_v
C--------
        d2_v = x13 * x13
        d3_p = 2.0d0 * x13
        d4_v = y13 * y13
        d2_p = 2.0d0 * y13
        d7_v = z13 * z13
        d1_p = 2.0d0 * z13
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d1_p * g_z13(g_i_) + d2_p * g_y13(g_i_) + d3_p 
     ** g_x13(g_i_)
        enddo
        d1_w = d2_v + d4_v + d7_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +196)
        endif
        do g_i_ = 1, g_p_
          g_rl(g_i_, 3) = d1_p * g_d1_w(g_i_)
        enddo
        rl(3) = d2_v
C--------
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
        do g_i_ = 1, g_p_
          g_bx(g_i_) = g_x32(g_i_)
        enddo
        bx = x32
C--------
        do g_i_ = 1, g_p_
          g_by(g_i_) = g_y32(g_i_)
        enddo
        by = y32
C--------
        do g_i_ = 1, g_p_
          g_bz(g_i_) = g_z32(g_i_)
        enddo
        bz = z32
C--------
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
     +252)
        endif
        do g_i_ = 1, g_p_
          g_rlr(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        rlr = d2_v
C--------
        d4_b = by + by
        d5_b = bx + bx
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_by(g_i_) + d5_b * g_bx(g_i_)
        enddo
        d2_w = bx * bx + by * by
        d4_b = bz + bz
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_bz(g_i_) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + bz * bz
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +276)
        endif
        do g_i_ = 1, g_p_
          g_rlb(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        rlb = d2_v
C--------
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = ry * g_by(g_i_) + by * g_ry(g_i_) + rx * g_bx(g
     *_i_) + bx * g_rx(g_i_)
        enddo
        d2_w = rx * bx + ry * by
        d6_v = (d2_w + rz * bz) * (d2_w + rz * bz)
        d1_p = 2.0d0 * (d2_w + rz * bz)
        d5_b = d1_p * bz
        d6_b = d1_p * rz
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d6_b * g_bz(g_i_) + d5_b * g_rz(g_i_) + d1_p * 
     *g_d2_w(g_i_)
        enddo
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
        do g_i_ = 1, g_p_
          g_bpr(g_i_) = d3_b * g_rlr(g_i_) + d4_b * g_d1_w(g_i_)
        enddo
        bpr = d4_v
C--------
        d2_v = rlb * rlb
        d2_p = 2.0d0 * rlb
        d4_v = bpr * bpr
        d1_p = 2.0d0 * bpr
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = -d1_p * g_bpr(g_i_) + d2_p * g_rlb(g_i_)
        enddo
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
        do g_i_ = 1, g_p_
          g_area(g_i_) = d5_b * g_d1_w(g_i_) + d3_b * g_rlr(g_i_)
        enddo
        area = rlr * d3_v / 2.0d+00
C--------
        d3_v = dens * area
        d2_b = 1.0d0 / dble(3.0)
        d3_b = d2_b * esp
        d4_b = d2_b * d3_v
        d5_b = d3_b * area
        d6_b = d3_b * dens
        do g_i_ = 1, g_p_
          g_rmas(g_i_) = d4_b * g_esp(g_i_) + d6_b * g_area(g_i_) + d5_b
     * * g_dens(g_i_)
        enddo
        rmas = d3_v * esp / dble(3.0)
C--------
C     
C constructing mass matrix
        do 99997 i = 1, 3
          i2 = 6 + i
          i3 = 12 + i
          do g_i_ = 1, g_p_
            g_fv(g_i_, i, i) = g_rmas(g_i_)
          enddo
          fv(i, i) = rmas
C--------
          do g_i_ = 1, g_p_
            g_fv(g_i_, i2, i2) = g_rmas(g_i_)
          enddo
          fv(i2, i2) = rmas
C--------
          do g_i_ = 1, g_p_
            g_fv(g_i_, i3, i3) = g_rmas(g_i_)
          enddo
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
          do g_i_ = 1, g_p_
            g_fv(g_i_, i1, i1) = d6_v * g_rmas(g_i_) + d7_b * g_rl(g_i_,
     * 3) + d8_b * g_rl(g_i_, 1)
          enddo
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
          do g_i_ = 1, g_p_
            g_fv(g_i_, i2, i2) = d6_v * g_rmas(g_i_) + d7_b * g_rl(g_i_,
     * 1) + d8_b * g_rl(g_i_, 2)
          enddo
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
          do g_i_ = 1, g_p_
            g_fv(g_i_, i3, i3) = d6_v * g_rmas(g_i_) + d7_b * g_rl(g_i_,
     * 2) + d8_b * g_rl(g_i_, 3)
          enddo
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
