C                           DISCLAIMER
C
C   This file was generated on 05/26/00 by the version of
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
C=====================================================================C
      subroutine gxmass7(elm, g_emass,a, g_a, rho, g_rho,  
     *                   x, g_x, y, g_y, z, g_z, gamma, g_grvfor,  
     *                   grvflg, totmas, masflg)
C=====================================================================C
C                                                                     C
C     Assemble the Elemental Mass Matrix of a 3D Timoschenko          C
C     Beam Finite Element. Lumping is Assumed Here.                   C
C                                                                     C
C     The Output Mass Matrix is a Block 12 by 12 Stored in the        C
C     Upper Left Corner of [emass] (First 12 Rows and Columns).       C
C                                                                     C
C     Francois M. Hemez - July 11th 1994 - Version 1.0                C
C                                                                     C
C     Original Parallel C Code By J.C. Chiou.                         C
C     See Directory: /mars/limbo/chiou/Beam.d/Timo.d.                 C
C     Gravity Body Force Modification November 1994 P.R. Stern        C
C                                                                     C
C=====================================================================C
C
C     ------------
C     DECLARATIONS
C     ------------
C
C.....GLOBAL VARIABLES
C
        integer elm
        real*8 a, rho, emass(12, 12), x(*), y(*), z(*)
        real*8 gamma(*), totmas
        logical grvflg, masflg
C
C.....LOCAL VARIABLES
C
        integer i, j
        real*8 zero, massdis, massrot, c1
        real*8 dx, dy, dz, length, two, twentyfour
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_emass, ldg_x, ldg_y, ldg_z, ldg_rho, ldg
     *_a, ldg_grvfor

C
      parameter (g_p_      = 1)
      parameter (ldg_emass = 1)
      parameter (ldg_a     = 1)
      parameter (ldg_rho   = 1)
      parameter (ldg_x     = 1)
      parameter (ldg_y     = 1)
      parameter (ldg_z     = 1)
      parameter (ldg_grvfor= 1) 
C
        double precision d8_b, d7_b, d2_w, d5_v, d6_b, d6_v, d1_p, d5_b,
     * d2_v, d3_v
        double precision d4_b, d2_b, d3_b, d1_w, g_emass(ldg_emass, 12, 
     *12), g_dx(g_pmax_), g_x(ldg_x, *), g_dy(g_pmax_), g_y(ldg_y, *), g
     *_dz(g_pmax_)
        double precision g_z(ldg_z, *), g_d2_w(g_pmax_), g_d1_w(g_pmax_)
     *, g_length(g_pmax_), g_massdis(g_pmax_), g_rho(ldg_rho), g_a(ldg_a
     *), g_massrot(g_pmax_), g_c1(g_pmax_), g_grvfor(ldg_grvfor, *)
        save g_dx, g_dy, g_dz, g_d2_w, g_d1_w, g_length, g_massdis, g_ma
     *ssrot, g_c1
        data zero /0.000000d+00/
        data two /0.200000d+01/
        data twentyfour /0.240000d+02/
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE OUTPUT MASS MATRIX
C
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'mass7','g_mass7.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do 99998 j = 1, 12
          do 99999 i = 1, 12
            do g_i_ = 1, g_p_
              g_emass(g_i_, i, j) = 0.0d0
            enddo
c            emass(i, j) = zero
C--------
1002        continue
99999     continue
1001      continue
99998   continue
C
C.....COMPUTE THE LENGTH OF THE BEAM ELEMENT
C
        do g_i_ = 1, g_p_
          g_dx(g_i_) = -g_x(g_i_, 1) + g_x(g_i_, 2)
        enddo
        dx = x(2) - x(1)
C--------
        do g_i_ = 1, g_p_
          g_dy(g_i_) = -g_y(g_i_, 1) + g_y(g_i_, 2)
        enddo
        dy = y(2) - y(1)
C--------
        do g_i_ = 1, g_p_
          g_dz(g_i_) = -g_z(g_i_, 1) + g_z(g_i_, 2)
        enddo
        dz = z(2) - z(1)
C--------
        d4_b = dy + dy
        d5_b = dx + dx
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_dy(g_i_) + d5_b * g_dx(g_i_)
        enddo
        d2_w = dx * dx + dy * dy
        d4_b = dz + dz
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_dz(g_i_) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + dz * dz
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDO (9,d1_w, d2_v, d1_p,
     +g_ehfid,
     +136)
        endif
        do g_i_ = 1, g_p_
          g_length(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        length = d2_v
C--------
C
        if (length .eq. zero) then
          goto 100
        endif
C
C.....INITIALIZE THE MASS COEFFICIENT FOR TRANSLATIONAL DOFs
C
        d3_v = rho * a
        d2_b = 1.0d0 / two
        d3_b = d2_b * length
        d4_b = d2_b * d3_v
        d5_b = d3_b * a
        d6_b = d3_b * rho
        do g_i_ = 1, g_p_
          g_massdis(g_i_) = d4_b * g_length(g_i_) + d6_b * g_a(g_i_) + d
     *5_b * g_rho(g_i_)
        enddo
        massdis = d3_v * length / two
C--------
C
C.....INITIALIZE THE MASS COEFFICIENT FOR ROTATIONAL DOFs
C
        d3_v = rho * a
        d5_v = d3_v * length
        d6_v = d5_v * length
        d2_b = 1.0d0 / twentyfour
        d3_b = d2_b * length
        d5_b = d3_b * length
        d6_b = d5_b * length
        d4_b = d2_b * d6_v + d3_b * d5_v + d5_b * d3_v
        d7_b = d6_b * a
        d8_b = d6_b * rho
        do g_i_ = 1, g_p_
          g_massrot(g_i_) = d4_b * g_length(g_i_) + d8_b * g_a(g_i_) + d
     *7_b * g_rho(g_i_)
        enddo
        massrot = d6_v * length / twentyfour
C--------
C
C.....ASSEMBLE THE (LUMPED) OUTPUT ELEMENTAL MASS MATRIX
C
        do g_i_ = 1, g_p_
          g_emass(g_i_, 1, 1) = g_massdis(g_i_)
        enddo
c        emass(1, 1) = massdis
C--------
        do g_i_ = 1, g_p_
          g_emass(g_i_, 2, 2) = g_massdis(g_i_)
        enddo
c        emass(2, 2) = massdis
C--------
        do g_i_ = 1, g_p_
          g_emass(g_i_, 3, 3) = g_massdis(g_i_)
        enddo
c        emass(3, 3) = massdis
C--------
C
        do g_i_ = 1, g_p_
          g_emass(g_i_, 4, 4) = g_massrot(g_i_)
        enddo
c        emass(4, 4) = massrot
C--------
        do g_i_ = 1, g_p_
          g_emass(g_i_, 5, 5) = g_massrot(g_i_)
        enddo
c        emass(5, 5) = massrot
C--------
        do g_i_ = 1, g_p_
          g_emass(g_i_, 6, 6) = g_massrot(g_i_)
        enddo
c        emass(6, 6) = massrot
C--------
C
        do g_i_ = 1, g_p_
          g_emass(g_i_, 7, 7) = g_massdis(g_i_)
        enddo
c        emass(7, 7) = massdis
C--------
        do g_i_ = 1, g_p_
          g_emass(g_i_, 8, 8) = g_massdis(g_i_)
        enddo
c        emass(8, 8) = massdis
C--------
        do g_i_ = 1, g_p_
          g_emass(g_i_, 9, 9) = g_massdis(g_i_)
        enddo
c        emass(9, 9) = massdis
C--------
C
        do g_i_ = 1, g_p_
          g_emass(g_i_, 10, 10) = g_massrot(g_i_)
        enddo
c        emass(10, 10) = massrot
C--------
        do g_i_ = 1, g_p_
          g_emass(g_i_, 11, 11) = g_massrot(g_i_)
        enddo
c        emass(11, 11) = massrot
C--------
        do g_i_ = 1, g_p_
          g_emass(g_i_, 12, 12) = g_massrot(g_i_)
        enddo
c        emass(12, 12) = massrot
C--------
C
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
        if (grvflg) then
          do g_i_ = 1, g_p_
            g_c1(g_i_) = 2.0d0 * g_massdis(g_i_)
          enddo
          c1 = 2.0d0 * massdis
C--------
          do g_i_ = 1, g_p_
            g_grvfor(g_i_, 1) = gamma(1) * g_c1(g_i_)
          enddo
c          grvfor(1) = c1 * gamma(1)
C--------
          do g_i_ = 1, g_p_
            g_grvfor(g_i_, 2) = gamma(2) * g_c1(g_i_)
          enddo
c          grvfor(2) = c1 * gamma(2)
C--------
          do g_i_ = 1, g_p_
            g_grvfor(g_i_, 3) = gamma(3) * g_c1(g_i_)
          enddo
c          grvfor(3) = c1 * gamma(3)
C--------
        endif
C
C.... ACCUMULATE THE SUBDOMAIN MASS
C
        if (masflg) then
          totmas = totmas + 2.0d0 * massdis
        endif
C     ------
C     RETURN
C     ------
C
        return
C
C     ---------------
C     ERROR-TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF THE BEAM ELEMENT HAS ZERO LENGTH
C
100     continue
        write (*, *) '*** FATAL ERROR in Routine MASS7  ***'
        write (*, *) '*** The Timoschenko Beam Element  ***'
        write (*, *) '*** Has Zero Length!              ***'
        write (*, *) '*** ... All Treatments Terminated ***'
        stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "MASS7"
C======================C
