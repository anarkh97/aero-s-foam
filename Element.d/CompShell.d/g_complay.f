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
C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine g_complay(g_p_, nlayer, idlayer, mtlayer, g_mtlayer, ld
     *g_mtlayer, x, g_x, ldg_x, y, g_y, ldg_y, z, type, ilayer, cstbb, g
     *_cstbb, ldg_cstbb, cstmm, g_cstmm, ldg_cstmm, cstbm, g_cstbm, ldg_
     *cstbm, cstmb, g_cstmb, ldg_cstmb, eframe, g_eframe, ldg_eframe, af
     *rame)
C=====================================================================C
C                                                                     C
C     Perform =   Assembles the Four 3 by 3 Constitutive Matrices for C
C     ---------   a given Layer of a Multi-layer Composite Laminate.  C
C                                                                     C
C                                                                     C
C     Input/Output =                                                  C
C     --------------                                                  C
C     NLAYER  <input>  number of layers of the composite element      C
C     IDLAYER <input>  identificators for each layer                  C
C     MTLAYER <input>  material properties of each layer              C
C     X       <input>  triangular coordinates in local x-direction    C
C     Y       <input>  triangular coordinates in local y-direction    C
C     Z       <input>  triangular coordinates in local z-direction    C
C     TYPE    <input>  type of constitutive law                       C
C     ILAYER  <input>  layer number                                   C
C     CSTBB   <output> 3 by 3  bending- bending constitutive matrix   C
C     CSTMM   <output> 3 by 3 membrane-membrane constitutive matrix   C
C     CSTBM   <output> 3 by 3  bending-membrane constitutive matrix   C
C     CSTMB   <output> 3 by 3 membrane- bending constitutive matrix   C
C     EFRAME  <input>  coordinates of the elemental frame system      C
C     AFRAME  <input>  frame system for orienting the fibers          C
C                                                                     C
C                                                                     C
C     Computations =                                                  C
C     --------------                                                  C
C     Four different matrices [D] are assembled for the given layer   C
C     number [ilayer] of the composite shell element:                 C
C                                                                     C
C                             [ [D_mm]  [D_mb] ]                      C
C     [Constitutive_Matrix] = [                ]                      C
C            6 by 6           [ [D_bm]  [D_bb] ]                      C
C                                                                     C
C     where "b" and "m" stand for bending and membrane, respectively. C
C                                                                     C
C     The constitutive matrix [D_bb] relates the element's curvatures C
C     to the bending moments and the constitutive matrix [D_mm]       C
C     relates the element's normal efforts to the strains. Similarly, C
C     the constitutive matrices [D_bm] and [D_mb] couple the bending  C
C     and membrane effects:                                           C
C                                                                     C
C     [ M_x  ]            [ k_x  ]                                    C
C     [ M_y  ] = [D_bb] * [ k_y  ]                                    C
C     [ M_xy ]            [ k_xy ]                                    C
C                                                                     C
C     [ N_x  ]            [ e_x  ]                                    C
C     [ N_y  ] = [D_mm] * [ e_y  ]                                    C
C     [ N_xy ]            [ e_xy ]                                    C
C                                                                     C
C     [ M_x  ]            [ e_x  ]                                    C
C     [ M_y  ] = [D_bm] * [ e_y  ]                                    C
C     [ M_xy ]            [ e_xy ]                                    C
C                                                                     C
C     [ N_x  ]            [ k_x  ]                                    C
C     [ N_y  ] = [D_mb] * [ k_y  ]                                    C
C     [ N_xy ]            [ k_xy ]                                    C
C                                                                     C
C     The symmetry constraints are:                                   C
C                                                                     C
C     [D_bb]^T = [D_bb]                                               C
C     [D_mm]^T = [D_mm]                                               C
C     [D_bm]^T = [D_mb]                                               C
C                                                                     C
C                                                                     C
C     Caution =   It is assumed that the layer has a constant         C
C     ---------   thickness so that no numerical interpolation is     C
C                 required.                                           C
C                                                                     C
C=====================================================================C
C=Author  = Francois M. Hemez                                         C
C=Date    = June 10th, 1995                                           C
C=Version = 2.0                                                       C
C=Comment =                                                           C
C=====================================================================C
C
C     ------------
C     DECLARATIONS
C     ------------
C
C.....Global Variables
C
        integer type, nlayer, idlayer(5, nlayer), ilayer
        real*8 cstbb(3, 3), cstmm(3, 3), cstbm(3, 3), cstmb(3, 3)
        real*8 x(3), y(3), z(3), mtlayer(8, nlayer)
        real*8 eframe(3, 3), aframe(3, 3)
C
C.....Local Variables
C
        integer i, j, layernumber
        real*8 zero, one, pi, twopi, intthick
        real*8 e1, e2, nu12, g12, mu1, mu2
        real*8 zsup, zinf, thetaf, thetad, theta
        real*8 s11, s12, s13, s22, s23, s33, dets
        real*8 q(3, 3), qbar(3, 3), t(3, 3), r(3)
        real*8 invt(3, 3), costheta, sintheta
        real*8 qt1, qt2, qt3, refvec(3)
        real*8 norm1, norm2, normref, proj1, proj2
        real*8 cosine1, cosine2, orifiber(3)
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_cstbb, ldg_cstmm, ldg_cstbm, ldg_cstmb, 
     *ldg_mtlayer, ldg_eframe, ldg_x, ldg_y
        double precision d8_b, d1_p, d2_w, d1_w, d7_b, d6_b, d5_b, d7_v,
     * d2_v, d5_v
        double precision d2_b, d3_v, d3_b, d4_v, d4_b, g_cstbb(ldg_cstbb
     *, 3, 3), g_cstmm(ldg_cstmm, 3, 3), g_cstbm(ldg_cstbm, 3, 3), g_cst
     *mb(ldg_cstmb, 3, 3), g_e1(g_pmax_)
        double precision g_mtlayer(ldg_mtlayer, 8, nlayer), g_e2(g_pmax_
     *), g_nu12(g_pmax_), g_g12(g_pmax_), g_mu1(g_pmax_), g_mu2(g_pmax_)
     *, g_thetaf(g_pmax_), g_zinf(g_pmax_), g_zsup(g_pmax_), g_orifiber(
     *g_pmax_, 3)
        double precision g_norm1(g_pmax_), g_norm2(g_pmax_), g_proj1(g_p
     *max_), g_proj2(g_pmax_), g_eframe(ldg_eframe, 3, 3), g_cosine1(g_p
     *max_), g_cosine2(g_pmax_), g_thetad(g_pmax_), g_d1_w(g_pmax_), g_t
     *heta(g_pmax_)
        double precision g_s11(g_pmax_), g_s12(g_pmax_), g_s13(g_pmax_),
     * g_s22(g_pmax_), g_s23(g_pmax_), g_s33(g_pmax_), g_d2_w(g_pmax_), 
     *g_dets(g_pmax_), g_q(g_pmax_, 3, 3), g_t(g_pmax_, 3, 3)
        double precision g_costheta(g_pmax_), g_sintheta(g_pmax_), g_inv
     *t(g_pmax_, 3, 3), g_qbar(g_pmax_, 3, 3), g_qt1(g_pmax_), g_qt2(g_p
     *max_), g_qt3(g_pmax_), g_intthick(g_pmax_), g_x(ldg_x, 3), g_y(ldg
     *_y, 3)
        save g_sintheta, g_invt, g_qbar, g_qt1, g_qt2, g_qt3, g_intthick
        save g_s12, g_s13, g_s22, g_s23, g_s33, g_d2_w, g_dets, g_q, g_t
     *, g_costheta
        save g_norm1, g_norm2, g_proj1, g_proj2, g_cosine1, g_cosine2, g
     *_thetad, g_d1_w, g_theta, g_s11
        save g_e1, g_e2, g_nu12, g_g12, g_mu1, g_mu2, g_thetaf, g_zinf, 
     *g_zsup, g_orifiber
        intrinsic dble
        data zero /0.000000d+00/
        data one /1.000000d+00/
C
C.....INITIALIZE THE MAIN DIAGONAL OF REUTER'S MATRIX
C
        data r /1.000000d+00, 1.000000d+00, 0.500000d+00/
C
C     -----
C     LOGIC
C     -----
C
        integer g_ehfid
        data g_ehfid /0/
C

C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        pi = acos(-one)
        twopi = 2.000000d+00 * pi
C
C.....CHECK IF TYPE OF CONTITUTIVE LAW IS CORRECT
C
        if ((type .ne. 2) .and. (type .ne. 3)) then
          goto 100
        endif
C
C.....CLEAR THE 3 BY 3 CONSTITUTIVE MATRICES
C
        do 99995 j = 1, 3
          do 99999 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstbb(g_i_, i, j) = 0.0d0
            enddo
            cstbb(i, j) = zero
C--------
1002        continue
99999     continue
          do 99998 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstmm(g_i_, i, j) = 0.0d0
            enddo
            cstmm(i, j) = zero
C--------
1003        continue
99998     continue
          do 99997 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstbm(g_i_, i, j) = 0.0d0
            enddo
            cstbm(i, j) = zero
C--------
1004        continue
99997     continue
          do 99996 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstmb(g_i_, i, j) = 0.0d0
            enddo
            cstmb(i, j) = zero
C--------
1005        continue
99996     continue
1001      continue
99995   continue
C
C.....INITIALIZE THE LAYER MATERIAL PROPERTIES
C
        layernumber = idlayer(3, ilayer)
C
        do g_i_ = 1, g_p_
          g_e1(g_i_) = g_mtlayer(g_i_, 1, ilayer)
        enddo
        e1 = mtlayer(1, ilayer)
C--------
        do g_i_ = 1, g_p_
          g_e2(g_i_) = g_mtlayer(g_i_, 2, ilayer)
        enddo
        e2 = mtlayer(2, ilayer)
C--------
        do g_i_ = 1, g_p_
          g_nu12(g_i_) = g_mtlayer(g_i_, 3, ilayer)
        enddo
        nu12 = mtlayer(3, ilayer)
C--------
        do g_i_ = 1, g_p_
          g_g12(g_i_) = g_mtlayer(g_i_, 4, ilayer)
        enddo
        g12 = mtlayer(4, ilayer)
C--------
        do g_i_ = 1, g_p_
          g_mu1(g_i_) = g_mtlayer(g_i_, 5, ilayer)
        enddo
        mu1 = mtlayer(5, ilayer)
C--------
        do g_i_ = 1, g_p_
          g_mu2(g_i_) = g_mtlayer(g_i_, 6, ilayer)
        enddo
        mu2 = mtlayer(6, ilayer)
C--------
        do g_i_ = 1, g_p_
          g_thetaf(g_i_) = g_mtlayer(g_i_, 8, ilayer)
        enddo
        thetaf = mtlayer(8, ilayer)
C--------
C
C.....CHECK FOR OBVIOUS ERRORS IN THE INPUT DATA
C
        if ((layernumber .lt. 1) .or. (layernumber .gt. nlayer)) then
          goto 200
        endif
C
        if (e1 .le. zero) then
          goto 300
        endif
        if (e2 .le. zero) then
          goto 300
        endif
        if (nu12 .le. zero) then
          goto 300
        endif
        if (g12 .le. zero) then
          goto 300
        endif
        if (mu1 .lt. zero) then
          goto 300
        endif
        if (mu2 .lt. zero) then
          goto 300
        endif
C
C.....TRANSFORM ANGLE IN THE RANGE BETWEEN 0-360      
C
        if ((thetaf .lt. zero) .or. (thetaf .gt. 360.00d+00)) then
          irot = thetaf / 360
          do g_i_ = 1, g_p_
            g_thetaf(g_i_) = g_thetaf(g_i_)
          enddo
          thetaf = thetaf - dble(real(irot)) * 360.0d0
C--------
          if (thetaf .lt. 0.0d0) then
            do g_i_ = 1, g_p_
              g_thetaf(g_i_) = -g_thetaf(g_i_)
            enddo
            thetaf = 360.0d0 - thetaf
C--------
          endif
        endif
C
C.....TRANSFORM FROM DEGREE TO RADIAN THE ANGLE BETWEEN THE
C.....REFERENCE ORIENTATION VECTOR AND THE ORIENTATION OF THE FIBERS
C
        d3_b = 1.0d0 / 180.00d+00 * pi
        do g_i_ = 1, g_p_
          g_thetaf(g_i_) = d3_b * g_thetaf(g_i_)
        enddo
        thetaf = pi * thetaf / 180.00d+00
C--------
C
C.....SET THE REFERENCE VECTOR FOR ORIENTING THE FIBERS OF THE LAYER
C.....(ALWAYS TAKE THE FIRST VECTOR OF THE FRAME)
C
        refvec(1) = aframe(1, 1)
        refvec(2) = aframe(2, 1)
        refvec(3) = aframe(3, 1)
C
C.....INITIALIZE THE UPPER AND LOWER [z] COORDINATES FOR THE LAYER
C
        do g_i_ = 1, g_p_
          g_zinf(g_i_) = -0.50d+00 * g_mtlayer(g_i_, 7, ilayer)
        enddo
        zinf = -0.50d+00 * mtlayer(7, ilayer)
C--------
        do g_i_ = 1, g_p_
          g_zsup(g_i_) = 0.50d+00 * g_mtlayer(g_i_, 7, ilayer)
        enddo
        zsup = 0.50d+00 * mtlayer(7, ilayer)
C--------
C
C.....PROJECT THE REFERENCE VECTOR INTO THE PLANE OF
C.....THE ELEMENT TO GET THE FIBER ORIENTATION VECTOR
C
        do g_i_ = 1, g_p_
          g_orifiber(g_i_, 1) = 0.0d0
        enddo
        orifiber(1) = zero
C--------
        do g_i_ = 1, g_p_
          g_orifiber(g_i_, 2) = 0.0d0
        enddo
        orifiber(2) = zero
C--------
        do g_i_ = 1, g_p_
          g_orifiber(g_i_, 3) = 0.0d0
        enddo
        orifiber(3) = zero
C--------
C
        do g_i_ = 1, g_p_
          g_norm1(g_i_) = 0.0d0
        enddo
        norm1 = zero
C--------
        do g_i_ = 1, g_p_
          g_norm2(g_i_) = 0.0d0
        enddo
        norm2 = zero
C--------
        normref = zero
        do g_i_ = 1, g_p_
          g_proj1(g_i_) = 0.0d0
        enddo
        proj1 = zero
C--------
        do g_i_ = 1, g_p_
          g_proj2(g_i_) = 0.0d0
        enddo
        proj2 = zero
C--------
C
        do 99994 i = 1, 3
          d4_b = eframe(i, 1) + eframe(i, 1)
          do g_i_ = 1, g_p_
            g_norm1(g_i_) = d4_b * g_eframe(g_i_, i, 1) + g_norm1(g_i_)
          enddo
          norm1 = norm1 + eframe(i, 1) * eframe(i, 1)
C--------
          d4_b = eframe(i, 2) + eframe(i, 2)
          do g_i_ = 1, g_p_
            g_norm2(g_i_) = d4_b * g_eframe(g_i_, i, 2) + g_norm2(g_i_)
          enddo
          norm2 = norm2 + eframe(i, 2) * eframe(i, 2)
C--------
          normref = normref + refvec(i) * refvec(i)
          do g_i_ = 1, g_p_
            g_proj1(g_i_) = refvec(i) * g_eframe(g_i_, i, 1) + g_proj1(g
     *_i_)
          enddo
          proj1 = proj1 + eframe(i, 1) * refvec(i)
C--------
          do g_i_ = 1, g_p_
            g_proj2(g_i_) = refvec(i) * g_eframe(g_i_, i, 2) + g_proj2(g
     *_i_)
          enddo
          proj2 = proj2 + eframe(i, 2) * refvec(i)
C--------
2001      continue
99994   continue
C
        d2_v = sqrt(norm1)
        if ( norm1 .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,norm1, d2_v, d1_p,
     +'g_complay.f',
     +416)
        endif
        do g_i_ = 1, g_p_
          g_norm1(g_i_) = d1_p * g_norm1(g_i_)
        enddo
        norm1 = d2_v
C--------
        d2_v = sqrt(norm2)
        if ( norm2 .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,norm2, d2_v, d1_p,
     +'g_complay.f',
     +429)
        endif
        do g_i_ = 1, g_p_
          g_norm2(g_i_) = d1_p * g_norm2(g_i_)
        enddo
        norm2 = d2_v
C--------
        normref = sqrt(normref)
C
        if (normref .eq. zero) then
          do g_i_ = 1, g_p_
            g_cosine1(g_i_) = 0.0d0
          enddo
          cosine1 = one
C--------
          do g_i_ = 1, g_p_
            g_cosine2(g_i_) = 0.0d0
          enddo
          cosine2 = zero
C--------
        else
          d3_v = norm1 * normref
          d4_v = proj1 / d3_v
          d2_b = 1.0d0 / d3_v
          d4_b = -d4_v / d3_v * normref
          do g_i_ = 1, g_p_
            g_cosine1(g_i_) = d4_b * g_norm1(g_i_) + d2_b * g_proj1(g_i_
     *)
          enddo
          cosine1 = d4_v
C--------
          d3_v = norm2 * normref
          d4_v = proj2 / d3_v
          d2_b = 1.0d0 / d3_v
          d4_b = -d4_v / d3_v * normref
          do g_i_ = 1, g_p_
            g_cosine2(g_i_) = d4_b * g_norm2(g_i_) + d2_b * g_proj2(g_i_
     *)
          enddo
          cosine2 = d4_v
C--------
        endif
C
        do 99993 i = 1, 3
          do g_i_ = 1, g_p_
            g_orifiber(g_i_, i) = cosine2 * g_eframe(g_i_, i, 2) + efram
     *e(i, 2) * g_cosine2(g_i_) + cosine1 * g_eframe(g_i_, i, 1) + efram
     *e(i, 1) * g_cosine1(g_i_)
          enddo
          orifiber(i) = cosine1 * eframe(i, 1) + cosine2 * eframe(i, 2)
C--------
2002      continue
99993   continue
C
C.....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL
C.....COORDINATE SYSTEM TO THE REFERENCE ORIENTATION VECTOR
C
        do g_i_ = 1, g_p_
          g_proj1(g_i_) = 0.0d0
        enddo
        proj1 = zero
C--------
        do g_i_ = 1, g_p_
          g_proj2(g_i_) = 0.0d0
        enddo
        proj2 = zero
C--------
        do g_i_ = 1, g_p_
          g_thetad(g_i_) = 0.0d0
        enddo
        thetad = zero
C--------
C
        do 99992 i = 1, 3
          do g_i_ = 1, g_p_
            g_proj1(g_i_) = eframe(i, 1) * g_orifiber(g_i_, i) + orifibe
     *r(i) * g_eframe(g_i_, i, 1) + g_proj1(g_i_)
          enddo
          proj1 = proj1 + eframe(i, 1) * orifiber(i)
C--------
          do g_i_ = 1, g_p_
            g_proj2(g_i_) = eframe(i, 2) * g_orifiber(g_i_, i) + orifibe
     *r(i) * g_eframe(g_i_, i, 2) + g_proj2(g_i_)
          enddo
          proj2 = proj2 + eframe(i, 2) * orifiber(i)
C--------
2003      continue
99992   continue
C
        if (proj1 .eq. zero) then
          if (proj2 .eq. zero) then
            goto 500
          endif
          if (proj2 .gt. zero) then
            do g_i_ = 1, g_p_
              g_thetad(g_i_) = 0.0d0
            enddo
            thetad = 0.50d+00 * pi
C--------
          endif
          if (proj2 .lt. zero) then
            do g_i_ = 1, g_p_
              g_thetad(g_i_) = 0.0d0
            enddo
            thetad = 1.50d+00 * pi
C--------
          endif
        else
          d3_v = proj2 / proj1
          d2_b = 1.0d0 / proj1
          d3_b = -d3_v / proj1
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_proj1(g_i_) + d2_b * g_proj2(g_i_)
          enddo
          d1_w = d3_v
          d2_v = atan(d1_w)
          d1_p = 1.0d0 / ( 1.0d0 + d1_w * d1_w )
          do g_i_ = 1, g_p_
            g_thetad(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          thetad = d2_v
C--------
        endif
C
        if (thetad .lt. zero) then
          do g_i_ = 1, g_p_
            g_thetad(g_i_) = g_thetad(g_i_)
          enddo
          thetad = thetad + twopi
C--------
        endif
C
        if ((thetad .lt. zero) .or. (thetad .gt. twopi)) then
          goto 600
        endif
C
C.....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL
C.....COORDINATE SYSTEM TO THE DIRECTION OF THE FIBER
C
        do g_i_ = 1, g_p_
          g_theta(g_i_) = g_thetaf(g_i_) + g_thetad(g_i_)
        enddo
        theta = thetad + thetaf
C--------
C
        if (theta .gt. twopi) then
          do g_i_ = 1, g_p_
            g_theta(g_i_) = g_theta(g_i_)
          enddo
          theta = theta - twopi
C--------
        endif
C
C.....CALCULATE THE COMPLIANCE MATRIX [S] WHICH RELATES THE STRESSES
C.....[s1], [s2] AND [s12] TO THE STRAINS [e1], [e2] AND [e12] IN
C.....THE COORDINATE SYSTEM {1;2} ASSOCIATED WITH THE FIBER ORIENTATION
C
        d2_v = one / e1
        d2_b = -d2_v / e1
        do g_i_ = 1, g_p_
          g_s11(g_i_) = d2_b * g_e1(g_i_)
        enddo
        s11 = d2_v
C--------
        d4_v = -nu12 / e1
        d3_b = -d4_v / e1
        d4_b = -(1.0d0 / e1)
        do g_i_ = 1, g_p_
          g_s12(g_i_) = d3_b * g_e1(g_i_) + d4_b * g_nu12(g_i_)
        enddo
        s12 = d4_v
C--------
        d3_v = mu1 / g12
        d2_b = 1.0d0 / g12
        d3_b = -d3_v / g12
        do g_i_ = 1, g_p_
          g_s13(g_i_) = d3_b * g_g12(g_i_) + d2_b * g_mu1(g_i_)
        enddo
        s13 = d3_v
C--------
        d2_v = one / e2
        d2_b = -d2_v / e2
        do g_i_ = 1, g_p_
          g_s22(g_i_) = d2_b * g_e2(g_i_)
        enddo
        s22 = d2_v
C--------
        d3_v = mu2 / g12
        d2_b = 1.0d0 / g12
        d3_b = -d3_v / g12
        do g_i_ = 1, g_p_
          g_s23(g_i_) = d3_b * g_g12(g_i_) + d2_b * g_mu2(g_i_)
        enddo
        s23 = d3_v
C--------
        d2_v = 1.0d0 / g12
        d2_b = -d2_v / g12
        do g_i_ = 1, g_p_
          g_s33(g_i_) = d2_b * g_g12(g_i_)
        enddo
        s33 = d2_v
C--------
C
C.....CALCULATE THE DETERMINANT OF THE COMPLIANCE MATRIX
C
        d7_v = s11 * s22 - s12 * s12
        d6_b = -s33 * s12 + (-s33) * s12
        d7_b = s33 * s22
        d8_b = s33 * s11
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d6_b * g_s12(g_i_) + d8_b * g_s22(g_i_) + d7_b 
     ** g_s11(g_i_) + d7_v * g_s33(g_i_)
        enddo
        d1_w = s33 * d7_v
        d4_v = s11 * s23
        d6_b = -s23 * s23
        d5_b = -d4_v + (-s23) * s11
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d5_b * g_s23(g_i_) + d6_b * g_s11(g_i_) + g_d1_
     *w(g_i_)
        enddo
        d2_w = d1_w - d4_v * s23
        d4_v = s22 * s13
        d6_b = -s13 * s13
        d5_b = -d4_v + (-s13) * s22
        do g_i_ = 1, g_p_
          g_dets(g_i_) = d5_b * g_s13(g_i_) + d6_b * g_s22(g_i_) + g_d2_
     *w(g_i_)
        enddo
        dets = d2_w - d4_v * s13
C--------
        d3_v = 2.00d+00 * s12
        d5_v = d3_v * s13
        d7_b = s23 * d3_v
        d8_b = s23 * s13 * 2.00d+00
        do g_i_ = 1, g_p_
          g_dets(g_i_) = d5_v * g_s23(g_i_) + d7_b * g_s13(g_i_) + d8_b 
     ** g_s12(g_i_) + g_dets(g_i_)
        enddo
        dets = dets + d5_v * s23
C--------
C
        if (dets .eq. zero) then
          goto 700
        endif
C
C.....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX WHICH RELATES
C.....THE STRAINS [e1], [e2] AND [e12] TO THE STRESSES [s1], [s2]
C.....AND [s12] IN THE COORDINATE SYSTEM {1;2} OF THE FIBER ORIENTATION
C
        do 99990 j = 1, 3
          do 99991 i = 1, 3
            do g_i_ = 1, g_p_
              g_q(g_i_, i, j) = 0.0d0
            enddo
            q(i, j) = zero
C--------
2203        continue
99991     continue
2202      continue
99990   continue
C
        d4_b = -s23 + (-s23)
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 1) = d4_b * g_s23(g_i_) + s22 * g_s33(g_i_) + s33
     * * g_s22(g_i_)
        enddo
        q(1, 1) = s22 * s33 - s23 * s23
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 2) = -s12 * g_s33(g_i_) + (-s33) * g_s12(g_i_) + 
     *s13 * g_s23(g_i_) + s23 * g_s13(g_i_)
        enddo
        q(1, 2) = s13 * s23 - s12 * s33
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 1, 3) = -s13 * g_s22(g_i_) + (-s22) * g_s13(g_i_) + 
     *s12 * g_s23(g_i_) + s23 * g_s12(g_i_)
        enddo
        q(1, 3) = s12 * s23 - s13 * s22
C--------
C
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 1) = g_q(g_i_, 1, 2)
        enddo
        q(2, 1) = q(1, 2)
C--------
        d4_b = -s13 + (-s13)
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 2) = d4_b * g_s13(g_i_) + s11 * g_s33(g_i_) + s33
     * * g_s11(g_i_)
        enddo
        q(2, 2) = s11 * s33 - s13 * s13
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 2, 3) = -s11 * g_s23(g_i_) + (-s23) * g_s11(g_i_) + 
     *s12 * g_s13(g_i_) + s13 * g_s12(g_i_)
        enddo
        q(2, 3) = s12 * s13 - s11 * s23
C--------
C
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 1) = g_q(g_i_, 1, 3)
        enddo
        q(3, 1) = q(1, 3)
C--------
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 2) = g_q(g_i_, 2, 3)
        enddo
        q(3, 2) = q(2, 3)
C--------
        d4_b = -s12 + (-s12)
        do g_i_ = 1, g_p_
          g_q(g_i_, 3, 3) = d4_b * g_s12(g_i_) + s11 * g_s22(g_i_) + s22
     * * g_s11(g_i_)
        enddo
        q(3, 3) = s11 * s22 - s12 * s12
C--------
C
        do 99988 j = 1, 3
          do 99989 i = 1, 3
            d3_v = q(i, j) / dets
            d2_b = 1.0d0 / dets
            d3_b = -d3_v / dets
            do g_i_ = 1, g_p_
              g_q(g_i_, i, j) = d3_b * g_dets(g_i_) + d2_b * g_q(g_i_, i
     *, j)
            enddo
            q(i, j) = d3_v
C--------
2205        continue
99989     continue
2204      continue
99988   continue
C
C.....INITIALIZE THE ROTATION MATRIX FROM THE FIBER COORDINATE
C.....SYSTEM {1;2} TO THE ELEMENT TRIANGULAR SYSTEM {x;y}
C
        do 99986 j = 1, 3
          do 99987 i = 1, 3
            do g_i_ = 1, g_p_
              g_t(g_i_, i, j) = 0.0d0
            enddo
            t(i, j) = zero
C--------
2207        continue
99987     continue
2206      continue
99986   continue
C
        d2_v = cos(theta)
        d1_p = -sin(theta)
        do g_i_ = 1, g_p_
          g_costheta(g_i_) = d1_p * g_theta(g_i_)
        enddo
        costheta = d2_v
C--------
        d2_v = sin(theta)
        d1_p = cos(theta)
        do g_i_ = 1, g_p_
          g_sintheta(g_i_) = d1_p * g_theta(g_i_)
        enddo
        sintheta = d2_v
C--------
C
        d2_b = costheta + costheta
        do g_i_ = 1, g_p_
          g_t(g_i_, 1, 1) = d2_b * g_costheta(g_i_)
        enddo
        t(1, 1) = costheta * costheta
C--------
        d2_b = sintheta + sintheta
        do g_i_ = 1, g_p_
          g_t(g_i_, 1, 2) = d2_b * g_sintheta(g_i_)
        enddo
        t(1, 2) = sintheta * sintheta
C--------
        d2_v = 2.00d+00 * costheta
        d4_b = sintheta * 2.00d+00
        do g_i_ = 1, g_p_
          g_t(g_i_, 1, 3) = d2_v * g_sintheta(g_i_) + d4_b * g_costheta(
     *g_i_)
        enddo
        t(1, 3) = d2_v * sintheta
C--------
C
        d2_b = sintheta + sintheta
        do g_i_ = 1, g_p_
          g_t(g_i_, 2, 1) = d2_b * g_sintheta(g_i_)
        enddo
        t(2, 1) = sintheta * sintheta
C--------
        d2_b = costheta + costheta
        do g_i_ = 1, g_p_
          g_t(g_i_, 2, 2) = d2_b * g_costheta(g_i_)
        enddo
        t(2, 2) = costheta * costheta
C--------
        d2_v = -2.00d+00 * costheta
        d4_b = sintheta * (-2.00d+00)
        do g_i_ = 1, g_p_
          g_t(g_i_, 2, 3) = d2_v * g_sintheta(g_i_) + d4_b * g_costheta(
     *g_i_)
        enddo
        t(2, 3) = d2_v * sintheta
C--------
C
        do g_i_ = 1, g_p_
          g_t(g_i_, 3, 1) = -costheta * g_sintheta(g_i_) + (-sintheta) *
     * g_costheta(g_i_)
        enddo
        t(3, 1) = -costheta * sintheta
C--------
        do g_i_ = 1, g_p_
          g_t(g_i_, 3, 2) = costheta * g_sintheta(g_i_) + sintheta * g_c
     *ostheta(g_i_)
        enddo
        t(3, 2) = costheta * sintheta
C--------
        d4_b = -sintheta + (-sintheta)
        d5_b = costheta + costheta
        do g_i_ = 1, g_p_
          g_t(g_i_, 3, 3) = d4_b * g_sintheta(g_i_) + d5_b * g_costheta(
     *g_i_)
        enddo
        t(3, 3) = costheta * costheta - sintheta * sintheta
C--------
C
C.....COMPUTE THE INVERSE OF [T]:
C.....[invT] = inverse(diag[R]) * [T]^t * diag[R]
C
        do 99984 j = 1, 3
          do 99985 i = 1, 3
            do g_i_ = 1, g_p_
              g_invt(g_i_, i, j) = 0.0d0
            enddo
            invt(i, j) = zero
C--------
2209        continue
99985     continue
2208      continue
99984   continue
C
        do 99982 j = 1, 3
          do 99983 i = 1, 3
            d2_b = r(j) / r(i)
            do g_i_ = 1, g_p_
              g_invt(g_i_, i, j) = d2_b * g_t(g_i_, j, i)
            enddo
            invt(i, j) = t(j, i) * (r(j) / r(i))
C--------
2211        continue
99983     continue
2210      continue
99982   continue
C
C.....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX WHICH RELATES
C.....THE STRAINS [ex], [ey] AND [exy] TO THE STRESSES [sx], [sy]
C.....AND [sxy] IN THE TRIANGULAR COORDINATE SYSTEM {x;y}:
C.....[Qbar] = [invT] * [Q] * [invT]^t
C
        do 99980 j = 1, 3
          do 99981 i = 1, 3
            do g_i_ = 1, g_p_
              g_qbar(g_i_, i, j) = 0.0d0
            enddo
            qbar(i, j) = zero
C--------
2213        continue
99981     continue
2212      continue
99980   continue
C
        do 99978 j = 1, 3
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = q(1, 2) * g_invt(g_i_, j, 2) + invt(j, 2) * g
     *_q(g_i_, 1, 2) + q(1, 1) * g_invt(g_i_, j, 1) + invt(j, 1) * g_q(g
     *_i_, 1, 1)
          enddo
          d1_w = q(1, 1) * invt(j, 1) + q(1, 2) * invt(j, 2)
          do g_i_ = 1, g_p_
            g_qt1(g_i_) = q(1, 3) * g_invt(g_i_, j, 3) + invt(j, 3) * g_
     *q(g_i_, 1, 3) + g_d1_w(g_i_)
          enddo
          qt1 = d1_w + q(1, 3) * invt(j, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = q(2, 2) * g_invt(g_i_, j, 2) + invt(j, 2) * g
     *_q(g_i_, 2, 2) + q(2, 1) * g_invt(g_i_, j, 1) + invt(j, 1) * g_q(g
     *_i_, 2, 1)
          enddo
          d1_w = q(2, 1) * invt(j, 1) + q(2, 2) * invt(j, 2)
          do g_i_ = 1, g_p_
            g_qt2(g_i_) = q(2, 3) * g_invt(g_i_, j, 3) + invt(j, 3) * g_
     *q(g_i_, 2, 3) + g_d1_w(g_i_)
          enddo
          qt2 = d1_w + q(2, 3) * invt(j, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = q(3, 2) * g_invt(g_i_, j, 2) + invt(j, 2) * g
     *_q(g_i_, 3, 2) + q(3, 1) * g_invt(g_i_, j, 1) + invt(j, 1) * g_q(g
     *_i_, 3, 1)
          enddo
          d1_w = q(3, 1) * invt(j, 1) + q(3, 2) * invt(j, 2)
          do g_i_ = 1, g_p_
            g_qt3(g_i_) = q(3, 3) * g_invt(g_i_, j, 3) + invt(j, 3) * g_
     *q(g_i_, 3, 3) + g_d1_w(g_i_)
          enddo
          qt3 = d1_w + q(3, 3) * invt(j, 3)
C--------
C
          do 99979 i = 1, 3
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = qt2 * g_invt(g_i_, i, 2) + invt(i, 2) * g_q
     *t2(g_i_) + qt1 * g_invt(g_i_, i, 1) + invt(i, 1) * g_qt1(g_i_)
            enddo
            d1_w = qt1 * invt(i, 1) + qt2 * invt(i, 2)
            do g_i_ = 1, g_p_
              g_qbar(g_i_, i, j) = qt3 * g_invt(g_i_, i, 3) + invt(i, 3)
     * * g_qt3(g_i_) + g_d1_w(g_i_)
            enddo
            qbar(i, j) = d1_w + qt3 * invt(i, 3)
C--------
2215        continue
99979     continue
C
2214      continue
99978   continue
C
C     ------------------------------------------------
C      COMPOSITE MATERIAL WITH KNOWN LAYER PROPERTIES
C     NO COUPLING BETWEEN BENDING AND MEMBRANE EFFECTS
C      (NUMERICAL INTEGRATION THROUGH THE THICKNESS)
C     ------------------------------------------------
C
        if (type .eq. 2) then
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
          d2_v = zsup * zsup
          d3_b = d2_v + zsup * zsup + zsup * zsup
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_zsup(g_i_)
          enddo
          d1_w = d2_v * zsup
          d3_v = zinf * zinf
          d2_b = 1.0d0 / 3.00d+00
          d5_b = -d2_b * zinf
          d6_b = -d2_b * d3_v + d5_b * zinf + d5_b * zinf
          do g_i_ = 1, g_p_
            g_intthick(g_i_) = d6_b * g_zinf(g_i_) + d2_b * g_d1_w(g_i_)
          enddo
          intthick = (d1_w - d3_v * zinf) / 3.00d+00
C--------
C
          do 99976 j = 1, 3
            do 99977 i = 1, 3
              do g_i_ = 1, g_p_
                g_cstbb(g_i_, i, j) = qbar(i, j) * g_intthick(g_i_) + in
     *tthick * g_qbar(g_i_, i, j) + g_cstbb(g_i_, i, j)
              enddo
              cstbb(i, j) = cstbb(i, j) + qbar(i, j) * intthick
C--------
3002          continue
99977       continue
3001        continue
99976     continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
          do g_i_ = 1, g_p_
            g_intthick(g_i_) = -g_zinf(g_i_) + g_zsup(g_i_)
          enddo
          intthick = zsup - zinf
C--------
C
          do 99974 j = 1, 3
            do 99975 i = 1, 3
              do g_i_ = 1, g_p_
                g_cstmm(g_i_, i, j) = qbar(i, j) * g_intthick(g_i_) + in
     *tthick * g_qbar(g_i_, i, j) + g_cstmm(g_i_, i, j)
              enddo
              cstmm(i, j) = cstmm(i, j) + qbar(i, j) * intthick
C--------
3004          continue
99975       continue
3003        continue
99974     continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C
          do 99972 j = 1, 3
            do 99973 i = 1, 3
              do g_i_ = 1, g_p_
                g_cstbm(g_i_, i, j) = 0.0d0
              enddo
              cstbm(i, j) = zero
C--------
3006          continue
99973       continue
3005        continue
99972     continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
          do 99970 j = 1, 3
            do 99971 i = 1, 3
              do g_i_ = 1, g_p_
                g_cstmb(g_i_, i, j) = 0.0d0
              enddo
              cstmb(i, j) = zero
C--------
3008          continue
99971       continue
3007        continue
99970     continue
C
C.....END OF TREATMENT FOR TYPE-2 CONSTITUTIVE LAW
C
        endif
C
C     --------------------------------------------------
C       COMPOSITE MATERIAL WITH KNOWN LAYER PROPERTIES
C     WITH COUPLING BETWEEN BENDING AND MEMBRANE EFFECTS
C       (NUMERICAL INTEGRATION THROUGH THE THICKNESS)
C     --------------------------------------------------
C
        if (type .eq. 3) then
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
          d2_v = zsup * zsup
          d3_b = d2_v + zsup * zsup + zsup * zsup
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_zsup(g_i_)
          enddo
          d1_w = d2_v * zsup
          d3_v = zinf * zinf
          d2_b = 1.0d0 / 3.00d+00
          d5_b = -d2_b * zinf
          d6_b = -d2_b * d3_v + d5_b * zinf + d5_b * zinf
          do g_i_ = 1, g_p_
            g_intthick(g_i_) = d6_b * g_zinf(g_i_) + d2_b * g_d1_w(g_i_)
          enddo
          intthick = (d1_w - d3_v * zinf) / 3.00d+00
C--------
C
          do 99968 j = 1, 3
            do 99969 i = 1, 3
              do g_i_ = 1, g_p_
                g_cstbb(g_i_, i, j) = qbar(i, j) * g_intthick(g_i_) + in
     *tthick * g_qbar(g_i_, i, j) + g_cstbb(g_i_, i, j)
              enddo
              cstbb(i, j) = cstbb(i, j) + qbar(i, j) * intthick
C--------
4002          continue
99969       continue
4001        continue
99968     continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
          do g_i_ = 1, g_p_
            g_intthick(g_i_) = -g_zinf(g_i_) + g_zsup(g_i_)
          enddo
          intthick = zsup - zinf
C--------
C
          do 99966 j = 1, 3
            do 99967 i = 1, 3
              do g_i_ = 1, g_p_
                g_cstmm(g_i_, i, j) = qbar(i, j) * g_intthick(g_i_) + in
     *tthick * g_qbar(g_i_, i, j) + g_cstmm(g_i_, i, j)
              enddo
              cstmm(i, j) = cstmm(i, j) + qbar(i, j) * intthick
C--------
4004          continue
99967       continue
4003        continue
99966     continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C.....(THE TRANSPOSED OF MATRIX [Qbar] IS TAKEN TO GET BEND-MEMB COUPLING)
C
          d5_b = -0.50d+00 * zinf + (-0.50d+00) * zinf
          d6_b = 0.50d+00 * zsup + 0.50d+00 * zsup
          do g_i_ = 1, g_p_
            g_intthick(g_i_) = d5_b * g_zinf(g_i_) + d6_b * g_zsup(g_i_)
          enddo
          intthick = 0.50d+00 * (zsup * zsup - zinf * zinf)
C--------
C
          do 99964 j = 1, 3
            do 99965 i = 1, 3
              do g_i_ = 1, g_p_
                g_cstbm(g_i_, i, j) = qbar(j, i) * g_intthick(g_i_) + in
     *tthick * g_qbar(g_i_, j, i) + g_cstbm(g_i_, i, j)
              enddo
              cstbm(i, j) = cstbm(i, j) + qbar(j, i) * intthick
C--------
4006          continue
99965       continue
4005        continue
99964     continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
          d5_b = -0.50d+00 * zinf + (-0.50d+00) * zinf
          d6_b = 0.50d+00 * zsup + 0.50d+00 * zsup
          do g_i_ = 1, g_p_
            g_intthick(g_i_) = d5_b * g_zinf(g_i_) + d6_b * g_zsup(g_i_)
          enddo
          intthick = 0.50d+00 * (zsup * zsup - zinf * zinf)
C--------
C
          do 99962 j = 1, 3
            do 99963 i = 1, 3
              do g_i_ = 1, g_p_
                g_cstmb(g_i_, i, j) = qbar(i, j) * g_intthick(g_i_) + in
     *tthick * g_qbar(g_i_, i, j) + g_cstmb(g_i_, i, j)
              enddo
              cstmb(i, j) = cstmb(i, j) + qbar(i, j) * intthick
C--------
4008          continue
99963       continue
4007        continue
99962     continue
C
C.....END OF TREATMENT FOR TYPE-3 CONSTITUTIVE LAW
C
        endif
C
C     ------
C     RETURN
C     ------
C
        return
C
C     ------
C     FORMAT
C     ------
C
91      format ('*** Type Used is: ',i10,12x,' ***')
C
C     ---------------
C     ERROR TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF THE CONSTITUTIVE LAW IS INCORRECT
C
100     continue
        write (*, *) '*** FATAL ERROR in Routine COMPLAY       ***'
        write (*, *) '*** Wrong Type of Constitutive Law       ***'
        write (*, 91) type
        write (*, *) '*** Types Allowed are:                   ***'
        write (*, *) '*** 2 = given layers properties          ***'
        write (*, *) '***     (no coupling bending/membrane)   ***'
        write (*, *) '*** 3 = given layers properties          ***'
        write (*, *) '***     (with coupling bending/membrane) ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE       ***'
        stop
C
C.....ERROR-MESSAGE IF A LAYER IDENTIFICATION NUMBER IS NOT CORRECT
C
200     continue
        write (*, *) '*** FATAL ERROR in Routine COMPLAY  ***'
        write (*, *) '*** The Local Layer Number is Not   ***'
        write (*, *) '*** Correct or Out of Bounds        ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE  ***'
        stop
C
C.....ERROR-MESSAGE IF A LAYER PROPERTY IS NOT CORRECT
C
300     continue
        write (*, *) '*** FATAL ERROR in Routine COMPLAY     ***'
        write (*, *) '*** One of the Layer Material Property ***'
        write (*, *) '*** is Negative or Zero!               ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE     ***'
        stop
C
C.....ERROR-MESSAGE IF THE FIBER ANGLE IS OUT-OF-BOUNDS
C
400     continue
        write (*, *) '*** FATAL ERROR in Routine COMPLAY    ***'
        write (*, *) '*** The Angle From the Reference      ***'
        write (*, *) '*** Direction to the Direction of     ***'
        write (*, *) '*** Fibers is Out-of-Bounds: it Must  ***'
        write (*, *) '*** be Within the Range 0-360 Degrees ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE    ***'
        stop
C
C.....ERROR-MESSAGE IF THE REFERENCE ORIENTATION IS BUGGY
C
500     continue
        write (*, *) '*** FATAL ERROR in Routine COMPLAY   ***'
        write (*, *) '*** The Reference Orientation Vector ***'
        write (*, *) '*** is Parallel to the Two Inplane   ***'
        write (*, *) '*** and Orthogonal Local Frames!     ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE   ***'
        stop
C
C.....ERROR-MESSAGE IF THE ORIENTATION ANGLE IS OUT-OF-BOUNDS
C
600     continue
        write (*, *) '*** FATAL ERROR in Routine COMPLAY    ***'
        write (*, *) '*** The Angle From the Local [x]      ***'
        write (*, *) '*** Axis of the Triangular Coordinate ***'
        write (*, *) '*** System to the Reference Direction ***'
        write (*, *) '*** is Out-of-Bounds: it Must be      ***'
        write (*, *) '*** Within the Range 0-2pi Radians    ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE    ***'
        stop
C
C.....ERROR-MESSAGE IF THE COMPLIANCE MATRIX IS SINGULAR
C
700     continue
        write (*, *) '*** FATAL ERROR in Routine COMPLAY    ***'
        write (*, *) '*** The Compliance Matrix is Singular ***'
        write (*, *) '*** ... Check Material Properties ... ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE    ***'
        stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPLAY"
C========================C
