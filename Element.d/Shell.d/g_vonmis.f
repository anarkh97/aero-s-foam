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
C -------------------------------------------------------------- C
C subroutine to compute the equivalent von mises stress          C
C in a plate or shell in the top surface (z=+t/2)                C
C and in the bottom surface (z=-t/2)                             C
C returns maximum value of the two stresses                      C
C                                                                C
C MODIFIED: March 3, 1997                                        C 
C BY: K. H. PIERSON                                              C
C WHY: the z component deviatoric stress was ignored previously  C
C                                                                C
C -------------------------------------------------------------- C
      subroutine g_vonmis(g_p_, rmx, g_rmx, ldg_rmx, rmy, g_rmy, ldg_rmy
     *, rmxy, g_rmxy, ldg_rmxy, rnx, g_rnx, ldg_rnx, rny, g_rny, ldg_rny
     *, rnxy, g_rnxy, ldg_rnxy, t, g_t, ldg_t, sv, g_sv, ldg_sv, surface
     *)
C
C.... GLOBAL VARIABLES
        integer surface
        double precision rmx, rmy, rmxy, rnx, rny, rnxy, t, sv
C
C.... LOCAL VARIABLES
C st = von mises stress in top surface
C sm = von mises stress in median surface
C sb = von mises stress in bottom surface
C
        double precision sx, sy, sxy, st, sb, sm, t2
        double precision rnxt, rnyt, rnxyt, rmxt, rmyt, rmxyt
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_t, ldg_rnx, ldg_rny, ldg_rnxy, ldg_rmx, 
     *ldg_rmy, ldg_rmxy, ldg_sv
        double precision d2_p, d4_b, d4_v, d3_b, d3_v, d1_p, d2_b, d2_v,
     * g_t(ldg_t), g_t2(g_pmax_)
        double precision g_rnxt(g_pmax_), g_rnx(ldg_rnx), g_rnyt(g_pmax_
     *), g_rny(ldg_rny), g_rnxyt(g_pmax_), g_rnxy(ldg_rnxy), g_rmxt(g_pm
     *ax_), g_rmx(ldg_rmx), g_rmyt(g_pmax_), g_rmy(ldg_rmy)
        double precision g_rmxyt(g_pmax_), g_rmxy(ldg_rmxy), g_sx(g_pmax
     *_), g_sy(g_pmax_), g_sxy(g_pmax_), g_sv(ldg_sv), g_sb(g_pmax_), g_
     *sm(g_pmax_), g_st(g_pmax_)
        save g_sb, g_sm, g_st
        save g_t2, g_rnxt, g_rnyt, g_rnxyt, g_rmxt, g_rmyt, g_rmxyt, g_s
     *x, g_sy, g_sxy
        external g_compj2
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'vonmis','g_vonmis.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        d2_v = abs(t)
        if (t .gt. 0.0d0) then
           d1_p =  1.0d0
        else if (t .lt. 0.0d0) then
           d1_p = -1.0d0
        else
           call ehufDO (3,t, d2_v, d1_p,
     +g_ehfid,
     +78)
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
        d3_v = rnx / t
        d2_b = 1.0d0 / t
        d3_b = -d3_v / t
        do g_i_ = 1, g_p_
          g_rnxt(g_i_) = d3_b * g_t(g_i_) + d2_b * g_rnx(g_i_)
        enddo
        rnxt = d3_v
C--------
        d3_v = rny / t
        d2_b = 1.0d0 / t
        d3_b = -d3_v / t
        do g_i_ = 1, g_p_
          g_rnyt(g_i_) = d3_b * g_t(g_i_) + d2_b * g_rny(g_i_)
        enddo
        rnyt = d3_v
C--------
        d3_v = rnxy / t
        d2_b = 1.0d0 / t
        d3_b = -d3_v / t
        do g_i_ = 1, g_p_
          g_rnxyt(g_i_) = d3_b * g_t(g_i_) + d2_b * g_rnxy(g_i_)
        enddo
        rnxyt = d3_v
C--------
C
        d4_v = dble(6.0) * rmx / t2
        d3_b = -d4_v / t2
        d4_b = 1.0d0 / t2 * dble(6.0)
        do g_i_ = 1, g_p_
          g_rmxt(g_i_) = d3_b * g_t2(g_i_) + d4_b * g_rmx(g_i_)
        enddo
        rmxt = d4_v
C--------
        d4_v = dble(6.0) * rmy / t2
        d3_b = -d4_v / t2
        d4_b = 1.0d0 / t2 * dble(6.0)
        do g_i_ = 1, g_p_
          g_rmyt(g_i_) = d3_b * g_t2(g_i_) + d4_b * g_rmy(g_i_)
        enddo
        rmyt = d4_v
C--------
        d4_v = dble(6.0) * rmxy / t2
        d3_b = -d4_v / t2
        d4_b = 1.0d0 / t2 * dble(6.0)
        do g_i_ = 1, g_p_
          g_rmxyt(g_i_) = d3_b * g_t2(g_i_) + d4_b * g_rmxy(g_i_)
        enddo
        rmxyt = d4_v
C--------
C
C ... COMPUTE VON MISES STRESS IN BOTTOM SURFACE
C
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
C
        call g_compj2(g_p_, sx, g_sx, g_pmax_, sy, g_sy, g_pmax_, sxy, g
     *_sxy, g_pmax_, sb, g_sb, g_pmax_)
C
        if (surface .eq. 3) then
          do g_i_ = 1, g_p_
            g_sv(g_i_) = g_sb(g_i_)
          enddo
          sv = sb
C--------
          return
        endif
C
C ... COMPUTE VON MISES STRESS IN MEDIAN SURFACE
C
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
C
          call g_compj2(g_p_, sx, g_sx, g_pmax_, sy, g_sy, g_pmax_, sxy,
     * g_sxy, g_pmax_, sm, g_sm, g_pmax_)
C
          do g_i_ = 1, g_p_
            g_sv(g_i_) = g_sm(g_i_)
          enddo
          sv = sm
C--------
          return
        endif
C
C ... COMPUTE VON MISES STRESS IN TOP SURFACE
C
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
        call g_compj2(g_p_, sx, g_sx, g_pmax_, sy, g_sy, g_pmax_, sxy, g
     *_sxy, g_pmax_, st, g_st, g_pmax_)
C
        if (surface .eq. 1) then
          do g_i_ = 1, g_p_
            g_sv(g_i_) = g_st(g_i_)
          enddo
          sv = st
C--------
          return
        endif
C
C ... VON MISES STRESS = max(sb,st)
C
        d3_v = max (sb, st)
        if (sb .gt.  st) then
           d1_p = 1.0d0
           d2_p = 0.0d0
        else if (sb .lt.  st) then
           d1_p = 0.0d0
           d2_p = 1.0d0
        else
           call ehbfDO (7,sb, st, d3_v, d1_p, d2_p,
     +g_ehfid,
     +244)
           d2_p = 1.0d0 -  d1_p
        endif
        do g_i_ = 1, g_p_
          g_sv(g_i_) = d2_p * g_st(g_i_) + d1_p * g_sb(g_i_)
        enddo
        sv = d3_v
C--------
C
        return
      end
C
C ... SUBROUTINE TO CALCULATE J2
C
      subroutine g_compj2(g_p_, sx, g_sx, ldg_sx, sy, g_sy, ldg_sy, sxy,
     * g_sxy, ldg_sxy, svm, g_svm, ldg_svm)
C
C ... GLOBAL VARIABLES
        double precision sx, sy, sxy, svm
C
C ... LOCAL VARIABLES
        double precision sz, s0, dsx, dsy, dsz, j2
C
C ... SET sz = 0 TO REMIND USER OF THIS ASSUMPTION
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_sx, ldg_sy, ldg_sxy, ldg_svm
        double precision d1_p, d8_b, d2_v, d5_b, d4_b, d3_b, d2_b, d1_w,
     * g_s0(g_pmax_), g_sx(ldg_sx)
        double precision g_sy(ldg_sy), g_dsx(g_pmax_), g_dsy(g_pmax_), g
     *_dsz(g_pmax_), g_d1_w(g_pmax_), g_j2(g_pmax_), g_sxy(ldg_sxy), g_s
     *vm(ldg_svm)
        save g_s0, g_dsx, g_dsy, g_dsz, g_d1_w, g_j2
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'compj2','g_vonmis.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        sz = 0.0d0
C
C ... COMPUTE AVERAGE HYDROSTATIC STRESS
        d3_b = 1.0d0 / 3.0d0
        do g_i_ = 1, g_p_
          g_s0(g_i_) = d3_b * g_sy(g_i_) + d3_b * g_sx(g_i_)
        enddo
        s0 = (sx + sy + sz) / 3.0d0
C--------
C
C ... COMPUTE DEVIATORIC STRESSES
        do g_i_ = 1, g_p_
          g_dsx(g_i_) = -g_s0(g_i_) + g_sx(g_i_)
        enddo
        dsx = sx - s0
C--------
        do g_i_ = 1, g_p_
          g_dsy(g_i_) = -g_s0(g_i_) + g_sy(g_i_)
        enddo
        dsy = sy - s0
C--------
        do g_i_ = 1, g_p_
          g_dsz(g_i_) = -g_s0(g_i_)
        enddo
        dsz = sz - s0
C--------
C
C ... RCFEM IGNORED THE DEVIATORIC Z COMPONENT: 
C     dsz  = 0.0
C
C ... COMPUTE J2
        d4_b = dsy + dsy
        d5_b = dsx + dsx
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_dsy(g_i_) + d5_b * g_dsx(g_i_)
        enddo
        d1_w = dsx * dsx + dsy * dsy
        d2_b = dble(0.5)
        d5_b = d2_b * sxy + d2_b * sxy
        d8_b = d2_b * dsz + d2_b * dsz
        do g_i_ = 1, g_p_
          g_j2(g_i_) = d5_b * g_sxy(g_i_) + d8_b * g_dsz(g_i_) + d2_b * 
     *g_d1_w(g_i_)
        enddo
        j2 = dble(0.5) * (d1_w + dsz * dsz + sxy * sxy)
C--------
C
C ... COMPUTE VON MISES STRESS
        d2_v = sqrt(j2)
        if ( j2 .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,j2, d2_v, d1_p,
     +g_ehfid,
     +341)
        endif
        do g_i_ = 1, g_p_
          g_svm(g_i_) = d1_p * g_j2(g_i_)
        enddo
        svm = d2_v
C--------
C
        return
      end
