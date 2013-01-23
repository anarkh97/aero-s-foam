C                           DISCLAIMER
C
C   This file was generated on 12/17/98 by the version of
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
      subroutine g_straineq(g_p_, rmx, g_rmx, ldg_rmx, rmy, g_rmy, ldg_r
     *my, rmxy, g_rmxy, ldg_rmxy, rnx, g_rnx, ldg_rnx, rny, g_rny, ldg_r
     *ny, rnxy, g_rnxy, ldg_rnxy, t, g_t, ldg_t, surface, ebar, g_ebar, 
     *ldg_ebar)
C**************************************************************
C       THIS ROUTINE CALCULATES THE EQUIVLAENT VON MISES      *
C       STRAIN FOR THE 3 NODE SHELL ELEMENT                   *
C**************************************************************
C                                                             *
C       AUTHOR  :       K.H. PIERSON                          *
C       DATE    :       MARCH 1997                            *
C       VERSION :       FEM-C++ 1.00                          *
C                                                             *
C**************************************************************
C                                                             *
C     rmx  = centroidal bending curvature (kxxc)            *
C     rmy  = centroidal bending curvature (kyyc)            *
C     rmxy = centroidal bending curvature (kxyc)            *
C     rnx  = centroidal membrane strain   (exxc)            *
C     rny  = centroidal membrane strain   (eyyc)            *
C     rnxy = centroidal membrane strain   (exyc)            *
C     t    = element thickness                              *
C     ebar = equivalent strain                              *
C                                                             *
C**************************************************************
C                                                             *
C       CALLED BY :  SANDS8.F                                 *
C                                                             *
C**************************************************************
C
C ... ARGUMENTS
        integer surface
        double precision rmx, rmy, rmxy, rnx, rny, rnxy, t, ebar
C
C ... LOCAL VARIABLES
        double precision ex, ey, exy, etop, ebot, emid, th
C
C ... RETURN IF MEDIAN SURFACE IS REQUESTED
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_ebar, ldg_rnx, ldg_rmx, ldg_rny, ldg_rmy
     *, ldg_rmxy, ldg_rnxy, ldg_t
        double precision d2_p, d1_p, d3_v, g_ebar(ldg_ebar), g_emid(g_pm
     *ax_), g_ex(g_pmax_), g_rnx(ldg_rnx), g_rmx(ldg_rmx), g_ey(g_pmax_)
     *, g_rny(ldg_rny)
        double precision g_rmy(ldg_rmy), g_etop(g_pmax_), g_ebot(g_pmax_
     *), g_rmxy(ldg_rmxy), g_rnxy(ldg_rnxy), g_t(ldg_t), g_exy(g_pmax_)
        save g_emid, g_ex, g_ey, g_etop, g_ebot, g_exy
        external g_equiv
        integer g_ehfid
        data g_ehfid /0/
C

C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        if (surface .eq. 2) then
          call g_equiv(g_p_, rmx, g_rmx, ldg_rmx, rmy, g_rmy, ldg_rmy, r
     *mxy, g_rmxy, ldg_rmxy, emid, g_emid, g_pmax_)
          do g_i_ = 1, g_p_
            g_ebar(g_i_) = g_emid(g_i_)
          enddo
          ebar = emid
C--------
          return
        endif
C
C ... DIVIDE THICKNESS BY 2
        th = 0.5 * t
C
C ... CONVERT CENTROIDAL BENDING CURVATURES TO BENDING STRAINS
C        rmx  = th/rmx
C        rmy  = th/rmy
C        rmxy = th/rmxy
C
C ... CALCULATE STRAINS AT TOP SURFACE
        do g_i_ = 1, g_p_
          g_ex(g_i_) = g_rmx(g_i_) + g_rnx(g_i_)
        enddo
        ex = rnx + rmx
C--------
        do g_i_ = 1, g_p_
          g_ey(g_i_) = g_rmy(g_i_) + g_rny(g_i_)
        enddo
        ey = rny + rmy
C--------
        exy = rnxy + rmxy
C
C ... COMPUTE EQUIVALENT STRAIN AT TOP SURFACE
        call g_equiv(g_p_, ex, g_ex, g_pmax_, ey, g_ey, g_pmax_, exy, g_
     *exy, g_pmax_, etop, g_etop, g_pmax_)
C
C ... RETURN IF TOP SURFACE VALUE IS REQUESTED
        if (surface .eq. 1) then
          do g_i_ = 1, g_p_
            g_ebar(g_i_) = g_etop(g_i_)
          enddo
          ebar = etop
C--------
          return
        endif
C
C ... CALCULATE STRAINS AT BOTTOM SURFACE
        do g_i_ = 1, g_p_
          g_ex(g_i_) = -g_rmx(g_i_) + g_rnx(g_i_)
        enddo
        ex = rnx - rmx
C--------
        do g_i_ = 1, g_p_
          g_ey(g_i_) = -g_rmy(g_i_) + g_rny(g_i_)
        enddo
        ey = rny - rmy
C--------
        exy = rnxy - rmxy
C
C ... COMPUTE EQUIVALENT STRAIN AT BOTTOM SURFACE
        call g_equiv(g_p_, ex, g_ex, g_pmax_, ey, g_ey, g_pmax_, exy, g_
     *exy, g_pmax_, ebot, g_ebot, g_pmax_)
C
C ... RETURN IF BOTTOM SURFACE VALUE IS REQUESTED
        if (surface .eq. 3) then
          do g_i_ = 1, g_p_
            g_ebar(g_i_) = g_ebot(g_i_)
          enddo
          ebar = ebot
C--------
          return
        endif
C
C
C ... RETURN THE MAXIMUM EQUIVALENT STRAIN
        d3_v = max (etop, ebot)
        if (etop .gt.  ebot) then
           d1_p = 1.0d0
           d2_p = 0.0d0
        else if (etop .lt.  ebot) then
           d1_p = 0.0d0
           d2_p = 1.0d0
        else
           call ehbfDV (7,etop, ebot, d3_v, d1_p, d2_p,
     +'g_straineq.f',
     +159)
           d2_p = 1.0d0 -  d1_p
        endif
        do g_i_ = 1, g_p_
          g_ebar(g_i_) = d2_p * g_ebot(g_i_) + d1_p * g_etop(g_i_)
        enddo
        ebar = d3_v
C--------
C
        return
      end
C
C
C ... SUBROUTINE TO CALCULATE EQUIVALENT STRAIN
C
      subroutine g_equiv(g_p_, ex, g_ex, ldg_ex, ey, g_ey, ldg_ey, exy, 
     *g_exy, ldg_exy, eq, g_eq, ldg_eq)
C
C ... ARGUMENTS
        double precision ex, ey, exy, eq
C
C ... LOCAL VARIABLES
        double precision e0, dex, dey, dez
C
C ... COMPUTE MEAN HYDROSTATIC STRAIN
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_ex, ldg_ey, ldg_eq, ldg_exy
        double precision d1_p, d6_b, d4_b, d2_v, d2_b, d2_w, d1_w, g_e0(
     *g_pmax_), g_ex(ldg_ex), g_ey(ldg_ey)
        double precision g_dex(g_pmax_), g_dey(g_pmax_), g_dez(g_pmax_),
     * g_d2_w(g_pmax_), g_d1_w(g_pmax_), g_eq(ldg_eq), g_exy(ldg_exy)
        save g_e0, g_dex, g_dey, g_dez, g_d2_w, g_d1_w
        integer g_ehfid
        data g_ehfid /0/
C

C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        d2_b = 1.0d0 / 3.0d0
        do g_i_ = 1, g_p_
          g_e0(g_i_) = d2_b * g_ey(g_i_) + d2_b * g_ex(g_i_)
        enddo
        e0 = (ex + ey) / 3.0d0
C--------
C
C ... COMPUTE DEVIATORIC STRAINS
        do g_i_ = 1, g_p_
          g_dex(g_i_) = -g_e0(g_i_) + g_ex(g_i_)
        enddo
        dex = ex - e0
C--------
        do g_i_ = 1, g_p_
          g_dey(g_i_) = -g_e0(g_i_) + g_ey(g_i_)
        enddo
        dey = ey - e0
C--------
        do g_i_ = 1, g_p_
          g_dez(g_i_) = -g_e0(g_i_)
        enddo
        dez = -e0
C--------
C
C ... COMPUTE EQUIVALENT STRAIN
        d2_v = 2.0d0 / 3.0d0 * dex
        d4_b = dey + dey
        d6_b = d2_v + dex * (2.0d0 / 3.0d0)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_dey(g_i_) + d6_b * g_dex(g_i_)
        enddo
        d2_w = d2_v * dex + dey * dey
        d4_b = dez + dez
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_dez(g_i_) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + dez * dez
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_straineq.f',
     +244)
        endif
        do g_i_ = 1, g_p_
          g_eq(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        eq = d2_v
C--------
C
        return
      end
