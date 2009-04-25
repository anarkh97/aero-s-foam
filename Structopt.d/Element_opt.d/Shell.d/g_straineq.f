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
      subroutine g_straineq(rmx, g_rmx, rmy, g_rmy
     *, rmxy, g_rmxy, rnx, g_rnx, rny, g_rny
     *, rnxy, g_rnxy, t, g_t, surface, ebar, g_ebar 
     *)
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
        double precision d2_p, d1_p, d3_v, g_ebar, g_emid
     *, g_ex, g_rnx, g_rmx, g_ey
     *, g_rny
        double precision g_rmy, g_etop, g_ebot
     *, g_rmxy, g_rnxy, g_t, g_exy
        save g_emid, g_ex, g_ey, g_etop, g_ebot, g_exy
        external g_equiv
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'straineq','g_straineq.f')
C
        if (surface .eq. 2) then
          call g_equiv(rmx, g_rmx, rmy, g_rmy, rmxy,
     * g_rmxy, emid, g_emid)
            g_ebar = g_emid
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
          g_ex = g_rmx + g_rnx
        ex = rnx + rmx
C--------
          g_ey = g_rmy + g_rny
        ey = rny + rmy
C--------
        exy = rnxy + rmxy
C
C ... COMPUTE EQUIVALENT STRAIN AT TOP SURFACE
        call g_equiv(ex, g_ex, ey, g_ey, exy, 
     *g_exy, etop, g_etop)
C
C ... RETURN IF TOP SURFACE VALUE IS REQUESTED
        if (surface .eq. 1) then
            g_ebar = g_etop
          ebar = etop
C--------
          return
        endif
C
C ... CALCULATE STRAINS AT BOTTOM SURFACE
          g_ex = -g_rmx + g_rnx
        ex = rnx - rmx
C--------
          g_ey = -g_rmy + g_rny
        ey = rny - rmy
C--------
        exy = rnxy - rmxy
C
C ... COMPUTE EQUIVALENT STRAIN AT BOTTOM SURFACE
        call g_equiv(ex, g_ex, ey, g_ey, exy, 
     *g_exy, ebot, g_ebot)
C
C ... RETURN IF BOTTOM SURFACE VALUE IS REQUESTED
        if (surface .eq. 3) then
            g_ebar = g_ebot
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
           call ehbfDO (7,etop, ebot, d3_v, d1_p, d2_p,
     +g_ehfid,
     +159)
           d2_p = 1.0d0 -  d1_p
        endif
          g_ebar = d2_p * g_ebot + d1_p * g_etop
        ebar = d3_v
C--------
C
        return
      end
C
C
C ... SUBROUTINE TO CALCULATE EQUIVALENT STRAIN
C
      subroutine g_equiv(ex, g_ex, ey, g_ey, exy, 
     *g_exy, eq, g_eq)
C
C ... ARGUMENTS
        double precision ex, ey, exy, eq
C
C ... LOCAL VARIABLES
        double precision e0, dex, dey, dez
C
C ... COMPUTE MEAN HYDROSTATIC STRAIN
        double precision d1_p, d6_b, d4_b, d2_v, d2_b, d2_w, d1_w, 
     *g_e0, g_ex, g_ey
        double precision g_dex, g_dey, g_dez,
     * g_d2_w, g_d1_w, g_eq, g_exy
        save g_e0, g_dex, g_dey, g_dez, g_d2_w, g_d1_w
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'equiv','g_straineq.f')
C
        d2_b = 1.0d0 / 3.0d0
          g_e0 = d2_b * g_ey + d2_b * g_ex
        e0 = (ex + ey) / 3.0d0
C--------
C
C ... COMPUTE DEVIATORIC STRAINS
          g_dex = -g_e0 + g_ex
        dex = ex - e0
C--------
          g_dey = -g_e0 + g_ey
        dey = ey - e0
C--------
          g_dez = -g_e0
        dez = -e0
C--------
C
C ... COMPUTE EQUIVALENT STRAIN
        d2_v = 2.0d0 / 3.0d0 * dex
        d4_b = dey + dey
        d6_b = d2_v + dex * (2.0d0 / 3.0d0)
          g_d2_w = d4_b * g_dey + d6_b * g_dex
        d2_w = d2_v * dex + dey * dey
        d4_b = dez + dez
          g_d1_w = d4_b * g_dez + g_d2_w
        d1_w = d2_w + dez * dez
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +244)
        endif
          g_eq = d1_p * g_d1_w
        eq = d2_v
C--------
C
        return
      end
