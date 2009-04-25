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
      subroutine gxcompms2(x,g_x,y,g_y,z,g_z,h,g_h,rho,g_rho,emass,
     *                     g_emass,medof,nttly,ncmpfr,elm,idlay,mtlay,
     *                     g_mtlay,cmpfr,iatt,ctyp,catt,cfrm,gamma,
     *                     grvfor,grvflg,totmas,masflg)
C=====================================================================C
C                                                                     C
C     Performs =   This subroutine will form the elemental mass       C
C     ----------   matrix of the 3D 3-node ANDES composite shell.     C
C                  Lumping is assumed here.                           C
C                                                                     C
C                                                                     C
C     Inputs/Outputs =                                                C
C     ----------------                                                C
C     X        <input>   nodal coordinates in the X-direction         C
C     Y        <input>   nodal coordinates in the Y-direction         C
C     Z        <input>   nodal coordinates in the Z-direction         C
C     H        <input>   element thicknesses (assumed constant)       C
C     RHO      <input>   density                                      C
C     EMASS    <output>  element mass matrix                          C
C     MEDOF    <input>   maximum number of DOFs per FE                C
C     NTTLY    <input>   total number of composite layers             C
C     NCMPFR   <input>   number of composite frames                   C
C     ELM      <input>   finite element number                        C
C     IDLAY    <input>   identificators for the composite layers      C
C     MTLAY    <input>   material properties of the composite layers  C
C     IATT     <input>   attribute number of the composite shell      C
C     CTYP     <input>   type of composite law (type 0, 1, 2, or 3)   C
C     CATT     <input>   composite attribute number of the element    C
C     CFRM     <input>   frame number of the composite shell element  C
C                                                                     C
C                                                                     C
C     Computations =                                                  C
C     ---------------                                                 C
C                                                                     C
C     The lumped mass matrix [M] is equal to:                         C
C                                                                     C
C           [ mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0 ]  C
C     [M] = [ 0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3 0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3 0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3]  C
C                                                                     C
C     with the following ordering of local degrees of freedom:        C
C                                                                     C
C                                [     U_x1 ]                         C
C                                [     U_y1 ]                         C
C                                [     U_z1 ]                         C
C                                [ theta_x1 ]                         C
C                                [ theta_y1 ]                         C
C                                [ theta_z1 ]                         C
C                                [     U_x2 ]                         C
C                                [     U_y2 ]                         C
C                                [     U_z2 ]                         C
C     Ordering_of_Local_DOFs  =  [ theta_x2 ]                         C
C                                [ theta_y2 ]                         C
C                                [ theta_z2 ]                         C
C                                [     U_x3 ]                         C
C                                [     U_y3 ]                         C
C                                [     U_z3 ]                         C
C                                [ theta_x3 ]                         C
C                                [ theta_y3 ]                         C
C                                [ theta_z3 ]                         C
C                                                                     C
C     No rotation of local-to-global basis is implemented since the   C
C     mass matrix [M] is formed of 3 by 3 blocks proportional to the  C
C     identity. The lumping factors are equal to:                     C
C                                                                     C
C     [mt] = [rho] * [A] * [h]        /    3.0                        C
C     [m1] = [rho] * [A] * [h] * [Ix] / 1260.0                        C
C     [m2] = [rho] * [A] * [h] * [Iy] / 1260.0                        C
C     [m3] = [rho] * [A] * [h] * [Iz] / 1260.0                        C
C                                                                     C
C     where:                                                          C
C                                                                     C
C     [rho]                 equivalent density of the composite       C
C     [h]                   equivalent thickness of the composite     C
C     [A]                   area                                      C
C     [Ix], [Iy] and [Iz]   equivalent pseudo-moments of inertia      C
C                                                                     C
C     The computation of these various quantities is performed        C
C     according to the type of composite law prescribed by the        C
C     variable [CTYP]. Four types of composite laws are available:    C
C                                                                     C
C     type-0: isotropic element                                       C
C                                                                     C
C     type-1: constitutive coefficients are given                     C
C                                                                     C
C     type-2: properties of each layer are given and no coupling      C
C             is assumed between bending and membrane                 C
C                                                                     C
C     type-3: properties of each layer are given and coupling         C
C             between bending and membrane is assumed                 C
C                                                                     C
C     The meaning of each input argument is specified in the          C
C     following as a function of the type of constitutive law.        C
C                                                                     C
C     1. Constitutive Law of Type-0                                   C
C     - - - - - - - - - - - - - - -                                   C
C                                                                     C
C     [rho]                 density of the isotropic material         C
C     [h]                   thickness of the isotropic shell          C
C     [A]                   area of the shell                         C
C     [Ix], [Iy] and [Iz]   pseudo-moments of inertia of the shell    C
C                                                                     C
C     Quantities [A], [Ix], [Iy] and [Iz] are obtained via            C
C     numerical integration.                                          C
C                                                                     C
C     2. Constitutive Law of Type-1                                   C
C     - - - - - - - - - - - - - - -                                   C
C                                                                     C
C     The density parameter must be initialized in the input file as  C
C     the density per unit surface of the composite element. In that  C
C     case, the thickness [h] is not required whereas quantities [A], C
C     [Ix], [Iy] and [Iz] are obtained via numerical integration once C
C     again. The lumping factors are equal to:                        C
C                                                                     C
C     [mt] = [rho_surface] * [A]        /    3.0                      C
C     [m1] = [rho_surface] * [A] * [Ix] / 1260.0                      C
C     [m2] = [rho_surface] * [A] * [Iy] / 1260.0                      C
C     [m3] = [rho_surface] * [A] * [Iz] / 1260.0                      C
C                                                                     C
C     where the density per unit surface [rho_surface] is stored in   C
C     the same variable as before (type-0), that is, [rho].           C
C                                                                     C
C     3. Constitutive Law of Type-2                                   C
C     - - - - - - - - - - - - - - -                                   C
C                                                                     C
C     The mass coefficients are computed by adding together the       C
C     contributions of each layer of the composite shell:             C
C                                                                     C
C     [mt] = sum{ [rho_k] * [h_k] } * [A]        /    3.0             C
C     [m1] = sum{ [rho_k] * [h_k] } * [A] * [Ix] / 1260.0             C
C     [m2] = sum{ [rho_k] * [h_k] } * [A] * [Iy] / 1260.0             C
C     [m3] = sum{ [rho_k] * [h_k] } * [A] * [Iz] / 1260.0             C
C                                                                     C
C     where:                                                          C
C                                                                     C
C     [rho_k]               density of the layer number [k]           C
C     [h_k]                 thickness of the layer number [k]         C
C     [A]                   area of the shell                         C
C     [Ix], [Iy] and [Iz]   pseudo-moments of inertia of the shell    C
C                                                                     C
C     Quantities [A], [Ix], [Iy] and [Iz] are obtained via            C
C     numerical integration. Densities (per unit volume) [rho_k] and  C
C     thicknesses [h_k] for each layer [k] of the composite are       C
C     given in the input file and retrieved from arrays [IDLAY] and   C
C     [MTLAY] using the information stored in the column that         C
C     corresponds to the particular layer number [k]:                 C
C                                                                     C
C     [IDLAY] Row 1: attribute number as read in the input file       C
C     [IDLAY] Row 2: number of layers of the composite shell          C
C     [IDLAY] Row 3: layer number (that is, number [k])               C
C     [IDLAY] Row 4: type of constitutive law (2 or 3)                C
C     [IDLAY] Row 5: frame number for the fibers reference vector     C
C                                                                     C
C     [MTLAY] Row 1: orthotropic Young modulus E1                     C
C     [MTLAY] Row 2: orthotropic Young modulus E2                     C
C     [MTLAY] Row 3: orthotropic Poisson's ratio nu12                 C
C     [MTLAY] Row 4: orthotropic shear mogulus G12                    C
C     [MTLAY] Row 5: first coefficient of mutual influence mu1,12     C
C     [MTLAY] Row 6: second coefficient of mutual influence mu2,12    C
C     [MTLAY] Row 7: density of the layer number [k]                  C
C     [MTLAY] Row 8: thickness of the layer number [k]                C
C     [MTLAY] Row 9: angle (degrees) between fibers and reference     C
C                                                                     C
C     4. Constitutive Law of Type-3                                   C
C     - - - - - - - - - - - - - - -                                   C
C                                                                     C
C     Same as Type-2.                                                 C
C                                                                     C
C                                                                     C
C     Caution =                                                       C
C     ---------                                                       C
C     The finite element is assumed to have constant thickness so     C
C     that no numerical interpolation is required. The array storing  C
C     the thicknesses at each nodal points [h] is reduced to a scalar C
C     in this routine. Also, The outputed element mass matrix is a    C
C     (diagonal) 18 by 18 block stored in the left upper corner of    C
C     the output storage [EMASS].                                     C
C                                                                     C
C                                                                     C
C     Outputs =   no output.                                          C
C     ---------                                                       C
C                                                                     C
C=====================================================================C
C=Author  = Francois M. Hemez                                         C
C=Date    = June 9, 1995                                              C
C=Version = 2.0                                                       C
C=Comment =                                                           C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Globlal Variables
C
        integer nttly, ctyp, catt, cfrm, elm
        integer idlay(5, nttly), iatt, ncmpfr, medof
C
        real*8 rho, emass(medof, medof), mtlay(9, nttly)
        real*8 h(3), cmpfr(9, ncmpfr), x(3), y(3), z(3)
        real*8 totmas, gamma(*), grvfor(*)
C
        logical grvflg, masflg
C
C.....Local Maximum Number of Layers per Composite Element
C
        integer maxlayer
        parameter (maxlayer = 1000)
C
C.....Local Variables
C
        integer i, j, i1, i2, i3, ilayer, nlayer
C
        real*8 zero, thick, mass0, mass1, mass2, mass3
        real*8 x13, y13, z13, hlayer(maxlayer)
        real*8 x32, y32, z32, rholayer(maxlayer)
        real*8 x21, y21, z21, rhoh
        real*8 dist(3), rlr, rlb, bpr, area, twicearea2
        real*8 ix, iy, iz
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_emass, ldg_mtlay, ldg_x, ldg_y, ldg_z, l
     *dg_h, ldg_rho
c
c manually inserted - begin	     
c
        parameter (g_p_=1,ldg_x=1,ldg_y=1,ldg_z=1,ldg_h=1,
     *             ldg_rho=1,ldg_emass=1,ldg_mtlay=1)
c
c manually inserted - end      
c
        double precision d8_b, d2_w, d7_b, d6_b, d4_v, d5_v, d1_p, d5_b,
     * d2_v, d3_v
        double precision d4_b, d2_b, d3_b, d1_w, g_mass0(g_pmax_), g_mas
     *s1(g_pmax_), g_mass2(g_pmax_), g_mass3(g_pmax_), g_rholayer(g_pmax
     *_, maxlayer), g_hlayer(g_pmax_, maxlayer)
        double precision g_emass(ldg_emass, medof, medof), g_mtlay(ldg_m
     *tlay, 9, nttly), g_x21(g_pmax_), g_x(ldg_x, 3), g_y21(g_pmax_), g_
     *y(ldg_y, 3), g_z21(g_pmax_), g_z(ldg_z, 3), g_x32(g_pmax_), g_y32(
     *g_pmax_)
        double precision g_z32(g_pmax_), g_x13(g_pmax_), g_y13(g_pmax_),
     * g_z13(g_pmax_), g_d2_w(g_pmax_), g_d1_w(g_pmax_), g_dist(g_pmax_,
     * 3), g_rlr(g_pmax_), g_rlb(g_pmax_), g_bpr(g_pmax_)
        double precision g_twicearea2(g_pmax_), g_area(g_pmax_), g_ix(g_
     *pmax_), g_iy(g_pmax_), g_iz(g_pmax_), g_thick(g_pmax_), g_h(ldg_h,
     * 3), g_rho(ldg_rho), g_rhoh(g_pmax_)
        save g_bpr, g_twicearea2, g_area, g_ix, g_iy, g_iz, g_thick, g_r
     *hoh
        save g_y32, g_z32, g_x13, g_y13, g_z13, g_d2_w, g_d1_w, g_dist, 
     *g_rlr, g_rlb
        save g_mass0, g_mass1, g_mass2, g_mass3, g_rholayer, g_hlayer, g
     *_x21, g_y21, g_z21, g_x32
        data zero /0.000000d+00/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE MASS LUMPING FACTORS
C
        integer g_ehfid
        data g_ehfid /0/
C

C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do g_i_ = 1, g_p_
          g_mass0(g_i_) = 0.0d0
        enddo
        mass0 = zero
C--------
        do g_i_ = 1, g_p_
          g_mass1(g_i_) = 0.0d0
        enddo
        mass1 = zero
C--------
        do g_i_ = 1, g_p_
          g_mass2(g_i_) = 0.0d0
        enddo
        mass2 = zero
C--------
        do g_i_ = 1, g_p_
          g_mass3(g_i_) = 0.0d0
        enddo
        mass3 = zero
C--------
C
C.....CLEAR THE LOCAL STORAGES FOR DENSITY AND THICKNESS PER LAYER
C
        do 99999 ilayer = 1, maxlayer
          do g_i_ = 1, g_p_
            g_rholayer(g_i_, ilayer) = 0.0d0
          enddo
          rholayer(ilayer) = zero
C--------
1001      continue
99999   continue
C
        do 99998 ilayer = 1, maxlayer
          do g_i_ = 1, g_p_
            g_hlayer(g_i_, ilayer) = 0.0d0
          enddo
          hlayer(ilayer) = zero
C--------
1002      continue
99998   continue
C
C.....CLEAR THE NUMBER OF LAYERS OF THE ELEMENT
C
        nlayer = 0
C
C.....CLEAR THE OUTPUT MASS MATRIX
C
        do 99996 j = 1, medof
          do 99997 i = 1, medof
            do g_i_ = 1, g_p_
              g_emass(g_i_, i, j) = 0.0d0
            enddo
            emass(i, j) = zero
C--------
1004        continue
99997     continue
1003      continue
99996   continue
C
C     --------------------------
C     CHECKS AND INITIALIZATIONS
C     --------------------------
C
C.....CHECK THE TYPE OF CONSTITUTIVE LAW FOR THIS ELEMENT
C
        if ((ctyp .ne. 0) .and. (ctyp .ne. 1) .and. (ctyp .ne. 2) .and. 
     *(ctyp .ne. 3)) then
          goto 100
        endif
C
C.....CHECK THE ADDRESSING IN ARRAYS [IDLAY] AND [MTLAY]
C
        if ((ctyp .eq. 2) .or. (ctyp .eq. 3)) then
          if ((catt .lt. 1) .or. (catt .gt. nttly)) then
            goto 200
          endif
        endif
C
C.....EXTRACT THE LAYER INFORMATION FOR TYPES 2 AND 3 LAWS
C
        if ((ctyp .eq. 2) .or. (ctyp .eq. 3)) then
          nlayer = idlay(2, catt)
          if (nlayer .gt. maxlayer) then
            goto 300
          endif
          do 99995 i = catt, (catt + nlayer - 1)
            ilayer = idlay(3, i)
            if (ilayer .gt. nlayer) then
              goto 400
            endif
            do g_i_ = 1, g_p_
              g_rholayer(g_i_, ilayer) = g_mtlay(g_i_, 7, i)
            enddo
            rholayer(ilayer) = mtlay(7, i)
C--------
            do g_i_ = 1, g_p_
              g_hlayer(g_i_, ilayer) = g_mtlay(g_i_, 8, i)
            enddo
            hlayer(ilayer) = mtlay(8, i)
C--------
2001        continue
99995     continue
        endif
C
C.....COMPUTE THE DISTANCE BETWEEN X-, Y- AND Z- NODAL COORDINATES
C
        do g_i_ = 1, g_p_
          g_x21(g_i_) = -g_x(g_i_, 1) + g_x(g_i_, 2)
        enddo
        x21 = x(2) - x(1)
C--------
        do g_i_ = 1, g_p_
          g_y21(g_i_) = -g_y(g_i_, 1) + g_y(g_i_, 2)
        enddo
        y21 = y(2) - y(1)
C--------
        do g_i_ = 1, g_p_
          g_z21(g_i_) = -g_z(g_i_, 1) + g_z(g_i_, 2)
        enddo
        z21 = z(2) - z(1)
C--------
C
        do g_i_ = 1, g_p_
          g_x32(g_i_) = -g_x(g_i_, 2) + g_x(g_i_, 3)
        enddo
        x32 = x(3) - x(2)
C--------
        do g_i_ = 1, g_p_
          g_y32(g_i_) = -g_y(g_i_, 2) + g_y(g_i_, 3)
        enddo
        y32 = y(3) - y(2)
C--------
        do g_i_ = 1, g_p_
          g_z32(g_i_) = -g_z(g_i_, 2) + g_z(g_i_, 3)
        enddo
        z32 = z(3) - z(2)
C--------
C
        do g_i_ = 1, g_p_
          g_x13(g_i_) = -g_x(g_i_, 3) + g_x(g_i_, 1)
        enddo
        x13 = x(1) - x(3)
C--------
        do g_i_ = 1, g_p_
          g_y13(g_i_) = -g_y(g_i_, 3) + g_y(g_i_, 1)
        enddo
        y13 = y(1) - y(3)
C--------
        do g_i_ = 1, g_p_
          g_z13(g_i_) = -g_z(g_i_, 3) + g_z(g_i_, 1)
        enddo
        z13 = z(1) - z(3)
C--------
C
C.....COMPUTE THE DISTANCE BETWEEN NODES 1-2, 2-3 AND 3-1
C
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
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compms.f',
     +474)
        endif
        do g_i_ = 1, g_p_
          g_dist(g_i_, 1) = d1_p * g_d1_w(g_i_)
        enddo
        dist(1) = d2_v
C--------
        d4_b = y32 + y32
        d5_b = x32 + x32
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_y32(g_i_) + d5_b * g_x32(g_i_)
        enddo
        d2_w = x32 * x32 + y32 * y32
        d4_b = z32 + z32
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_z32(g_i_) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + z32 * z32
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compms.f',
     +498)
        endif
        do g_i_ = 1, g_p_
          g_dist(g_i_, 2) = d1_p * g_d1_w(g_i_)
        enddo
        dist(2) = d2_v
C--------
        d4_b = y13 + y13
        d5_b = x13 + x13
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_y13(g_i_) + d5_b * g_x13(g_i_)
        enddo
        d2_w = x13 * x13 + y13 * y13
        d4_b = z13 + z13
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_z13(g_i_) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + z13 * z13
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compms.f',
     +522)
        endif
        do g_i_ = 1, g_p_
          g_dist(g_i_, 3) = d1_p * g_d1_w(g_i_)
        enddo
        dist(3) = d2_v
C--------
C
C.....COMPUTE THE LENGTH OF SIDE 1-2
C
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
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compms.f',
     +549)
        endif
        do g_i_ = 1, g_p_
          g_rlr(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        rlr = d2_v
C--------
C
C.....CHECK FOR ZERO-SIDE LENGTH
C
        if (rlr .eq. zero) then
          goto 500
        endif
C
C.....COMPUTE THE DISTANCE OF THE OPPOSING NODE (3) TO THAT SIDE (1-2)
C
        d4_b = y32 + y32
        d5_b = x32 + x32
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_y32(g_i_) + d5_b * g_x32(g_i_)
        enddo
        d2_w = x32 * x32 + y32 * y32
        d4_b = z32 + z32
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_z32(g_i_) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + z32 * z32
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compms.f',
     +582)
        endif
        do g_i_ = 1, g_p_
          g_rlb(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        rlb = d2_v
C--------
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = y21 * g_y32(g_i_) + y32 * g_y21(g_i_) + x21 * g
     *_x32(g_i_) + x32 * g_x21(g_i_)
        enddo
        d2_w = x21 * x32 + y21 * y32
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = z21 * g_z32(g_i_) + z32 * g_z21(g_i_) + g_d2_w(
     *g_i_)
        enddo
        d1_w = d2_w + z21 * z32
        d2_v = abs(d1_w)
        if (d1_w .gt. 0.0d0) then
           d1_p =  1.0d0
        else if (d1_w .lt. 0.0d0) then
           d1_p = -1.0d0
        else
           call ehufDV (3,d1_w, d2_v, d1_p,
     +'g_compms.f',
     +607)
        endif
        d4_v = d2_v / rlr
        d3_b = -d4_v / rlr
        d4_b = 1.0d0 / rlr * d1_p
        do g_i_ = 1, g_p_
          g_bpr(g_i_) = d3_b * g_rlr(g_i_) + d4_b * g_d1_w(g_i_)
        enddo
        bpr = d4_v
C--------
C
C.....COMPUTE THE SQUARE OF TWICE THE TRIANGLE'S AREA
C
        d4_b = -bpr + (-bpr)
        d5_b = rlb + rlb
        do g_i_ = 1, g_p_
          g_twicearea2(g_i_) = d4_b * g_bpr(g_i_) + d5_b * g_rlb(g_i_)
        enddo
        twicearea2 = rlb * rlb - bpr * bpr
C--------
C
C.....CHECK IF THE TRIANGLE'S AREA IS POSITIVE
C
        if (twicearea2 .le. zero) then
          goto 600
        endif
C
C.....COMPUTE THE AREA OF THE TRIANGLE
C
        d2_v = 0.50d+00 * rlr
        d4_v = sqrt(twicearea2)
        if ( twicearea2 .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d4_v)
        else
           call ehufDV (9,twicearea2, d4_v, d1_p,
     +'g_compms.f',
     +643)
        endif
        d4_b = d2_v * d1_p
        d5_b = d4_v * 0.50d+00
        do g_i_ = 1, g_p_
          g_area(g_i_) = d4_b * g_twicearea2(g_i_) + d5_b * g_rlr(g_i_)
        enddo
        area = d2_v * d4_v
C--------
C
C.....COMPUTE THE THREE PSEUDO MOMENTS OF INERTIA
C
        d4_b = dist(3) + dist(3)
        d5_b = dist(1) + dist(1)
        do g_i_ = 1, g_p_
          g_ix(g_i_) = d4_b * g_dist(g_i_, 3) + d5_b * g_dist(g_i_, 1)
        enddo
        ix = dist(1) * dist(1) + dist(3) * dist(3)
C--------
        d4_b = dist(2) + dist(2)
        d5_b = dist(1) + dist(1)
        do g_i_ = 1, g_p_
          g_iy(g_i_) = d4_b * g_dist(g_i_, 2) + d5_b * g_dist(g_i_, 1)
        enddo
        iy = dist(1) * dist(1) + dist(2) * dist(2)
C--------
        d4_b = dist(3) + dist(3)
        d5_b = dist(2) + dist(2)
        do g_i_ = 1, g_p_
          g_iz(g_i_) = d4_b * g_dist(g_i_, 3) + d5_b * g_dist(g_i_, 2)
        enddo
        iz = dist(2) * dist(2) + dist(3) * dist(3)
C--------
C
C     -----------------------
C     TYPE-0 CONSTITUTIVE LAW
C     -----------------------
C
        if (ctyp .eq. 0) then
C
C.....INITIALIZE THE ELEMENT'S CONSTANT THICKNESS
C
          do g_i_ = 1, g_p_
            g_thick(g_i_) = g_h(g_i_, 1)
          enddo
          thick = h(1)
C--------
C
C.....FORM THE MASS COEFFICIENTS PER DEGREE OF FREEDOM
C
          d3_v = rho * area
          d2_b = 1.0d0 / 3.00d+00
          d3_b = d2_b * thick
          d4_b = d2_b * d3_v
          d5_b = d3_b * area
          d6_b = d3_b * rho
          do g_i_ = 1, g_p_
            g_mass0(g_i_) = d4_b * g_thick(g_i_) + d6_b * g_area(g_i_) +
     * d5_b * g_rho(g_i_)
          enddo
          mass0 = d3_v * thick / 3.00d+00
C--------
          d3_v = rho * area
          d5_v = d3_v * thick
          d2_b = 1.0d0 / 1260.00d+00
          d3_b = d2_b * ix
          d4_b = d2_b * d5_v
          d5_b = d3_b * thick
          d6_b = d3_b * d3_v
          d7_b = d5_b * area
          d8_b = d5_b * rho
          do g_i_ = 1, g_p_
            g_mass1(g_i_) = d4_b * g_ix(g_i_) + d6_b * g_thick(g_i_) + d
     *8_b * g_area(g_i_) + d7_b * g_rho(g_i_)
          enddo
          mass1 = d5_v * ix / 1260.00d+00
C--------
          d3_v = rho * area
          d5_v = d3_v * thick
          d2_b = 1.0d0 / 1260.00d+00
          d3_b = d2_b * iy
          d4_b = d2_b * d5_v
          d5_b = d3_b * thick
          d6_b = d3_b * d3_v
          d7_b = d5_b * area
          d8_b = d5_b * rho
          do g_i_ = 1, g_p_
            g_mass2(g_i_) = d4_b * g_iy(g_i_) + d6_b * g_thick(g_i_) + d
     *8_b * g_area(g_i_) + d7_b * g_rho(g_i_)
          enddo
          mass2 = d5_v * iy / 1260.00d+00
C--------
          d3_v = rho * area
          d5_v = d3_v * thick
          d2_b = 1.0d0 / 1260.00d+00
          d3_b = d2_b * iz
          d4_b = d2_b * d5_v
          d5_b = d3_b * thick
          d6_b = d3_b * d3_v
          d7_b = d5_b * area
          d8_b = d5_b * rho
          do g_i_ = 1, g_p_
            g_mass3(g_i_) = d4_b * g_iz(g_i_) + d6_b * g_thick(g_i_) + d
     *8_b * g_area(g_i_) + d7_b * g_rho(g_i_)
          enddo
          mass3 = d5_v * iz / 1260.00d+00
C--------
C
C.....END OF TREATMENT FOR A TYPE-0 CONSTITUTIVE LAW
C
        endif
C
C     -----------------------
C     TYPE-1 CONSTITUTIVE LAW
C     -----------------------
C
        if (ctyp .eq. 1) then
C
C.....FORM THE MASS COEFFICIENTS PER DEGREE OF FREEDOM
C.....WARNING: THE "DENSITY" PARAMETER PRESCRIBED IN THE
C.....ATTRIBUTE SECTION OF THE INPUT FILE AND USED HERE
C.....MUST BE A DENSITY PER UNIT SURFACE
C
          d2_b = 1.0d0 / 3.00d+00
          d3_b = d2_b * area
          d4_b = d2_b * rho
          do g_i_ = 1, g_p_
            g_mass0(g_i_) = d4_b * g_area(g_i_) + d3_b * g_rho(g_i_)
          enddo
          mass0 = rho * area / 3.00d+00
C--------
          d3_v = rho * area
          d2_b = 1.0d0 / 1260.00d+00
          d3_b = d2_b * ix
          d4_b = d2_b * d3_v
          d5_b = d3_b * area
          d6_b = d3_b * rho
          do g_i_ = 1, g_p_
            g_mass1(g_i_) = d4_b * g_ix(g_i_) + d6_b * g_area(g_i_) + d5
     *_b * g_rho(g_i_)
          enddo
          mass1 = d3_v * ix / 1260.00d+00
C--------
          d3_v = rho * area
          d2_b = 1.0d0 / 1260.00d+00
          d3_b = d2_b * iy
          d4_b = d2_b * d3_v
          d5_b = d3_b * area
          d6_b = d3_b * rho
          do g_i_ = 1, g_p_
            g_mass2(g_i_) = d4_b * g_iy(g_i_) + d6_b * g_area(g_i_) + d5
     *_b * g_rho(g_i_)
          enddo
          mass2 = d3_v * iy / 1260.00d+00
C--------
          d3_v = rho * area
          d2_b = 1.0d0 / 1260.00d+00
          d3_b = d2_b * iz
          d4_b = d2_b * d3_v
          d5_b = d3_b * area
          d6_b = d3_b * rho
          do g_i_ = 1, g_p_
            g_mass3(g_i_) = d4_b * g_iz(g_i_) + d6_b * g_area(g_i_) + d5
     *_b * g_rho(g_i_)
          enddo
          mass3 = d3_v * iz / 1260.00d+00
C--------
C
C.....END OF TREATMENT FOR A TYPE-1 CONSTITUTIVE LAW
C
        endif
C
C     -----------------------------------
C     TYPE-2 AND TYPE-3 CONSTITUTIVE LAWS
C     -----------------------------------
C
        if ((ctyp .eq. 2) .or. (ctyp .eq. 3)) then
C
C.....ACCUMULATE THE PRODUCT DENSITY BY THICKNESS PER LAYER
C
          do g_i_ = 1, g_p_
            g_rhoh(g_i_) = 0.0d0
          enddo
          rhoh = zero
C--------
C
          do 99994 ilayer = 1, nlayer
            do g_i_ = 1, g_p_
              g_rhoh(g_i_) = rholayer(ilayer) * g_hlayer(g_i_, ilayer) +
     * hlayer(ilayer) * g_rholayer(g_i_, ilayer) + g_rhoh(g_i_)
            enddo
            rhoh = rhoh + rholayer(ilayer) * hlayer(ilayer)
C--------
3001        continue
99994     continue
C
C.....FORM THE MASS COEFFICIENTS PER DEGREE OF FREEDOM
C
          d2_b = 1.0d0 / 3.00d+00
          d3_b = d2_b * area
          d4_b = d2_b * rhoh
          do g_i_ = 1, g_p_
            g_mass0(g_i_) = d4_b * g_area(g_i_) + d3_b * g_rhoh(g_i_)
          enddo
          mass0 = rhoh * area / 3.00d+00
C--------
          d3_v = rhoh * area
          d2_b = 1.0d0 / 1260.00d+00
          d3_b = d2_b * ix
          d4_b = d2_b * d3_v
          d5_b = d3_b * area
          d6_b = d3_b * rhoh
          do g_i_ = 1, g_p_
            g_mass1(g_i_) = d4_b * g_ix(g_i_) + d6_b * g_area(g_i_) + d5
     *_b * g_rhoh(g_i_)
          enddo
          mass1 = d3_v * ix / 1260.00d+00
C--------
          d3_v = rhoh * area
          d2_b = 1.0d0 / 1260.00d+00
          d3_b = d2_b * iy
          d4_b = d2_b * d3_v
          d5_b = d3_b * area
          d6_b = d3_b * rhoh
          do g_i_ = 1, g_p_
            g_mass2(g_i_) = d4_b * g_iy(g_i_) + d6_b * g_area(g_i_) + d5
     *_b * g_rhoh(g_i_)
          enddo
          mass2 = d3_v * iy / 1260.00d+00
C--------
          d3_v = rhoh * area
          d2_b = 1.0d0 / 1260.00d+00
          d3_b = d2_b * iz
          d4_b = d2_b * d3_v
          d5_b = d3_b * area
          d6_b = d3_b * rhoh
          do g_i_ = 1, g_p_
            g_mass3(g_i_) = d4_b * g_iz(g_i_) + d6_b * g_area(g_i_) + d5
     *_b * g_rhoh(g_i_)
          enddo
          mass3 = d3_v * iz / 1260.00d+00
C--------
C
C.....END OF TREATMENT FOR TYPE-2 AND TYPE-3 CONSTITUTIVE LAWS
C
        endif
C
C     -------------------------------------
C     ASSEMBLY OF THE ELEMENTAL MASS MATRIX
C     -------------------------------------
C
C.....FORM THE LUMPED ELEMENT MASS MATRIX
C
        do 99993 i = 1, 3
          i2 = i + 6
          i3 = i + 12
          do g_i_ = 1, g_p_
            g_emass(g_i_, i, i) = g_mass0(g_i_)
          enddo
          emass(i, i) = mass0
C--------
          do g_i_ = 1, g_p_
            g_emass(g_i_, i2, i2) = g_mass0(g_i_)
          enddo
          emass(i2, i2) = mass0
C--------
          do g_i_ = 1, g_p_
            g_emass(g_i_, i3, i3) = g_mass0(g_i_)
          enddo
          emass(i3, i3) = mass0
C--------
4001      continue
99993   continue
C
        do 99992 i = 1, 3
          i1 = i + 3
          i2 = i + 9
          i3 = i + 15
          do g_i_ = 1, g_p_
            g_emass(g_i_, i1, i1) = g_mass1(g_i_)
          enddo
          emass(i1, i1) = mass1
C--------
          do g_i_ = 1, g_p_
            g_emass(g_i_, i2, i2) = g_mass2(g_i_)
          enddo
          emass(i2, i2) = mass2
C--------
          do g_i_ = 1, g_p_
            g_emass(g_i_, i3, i3) = g_mass3(g_i_)
          enddo
          emass(i3, i3) = mass3
C--------
4002      continue
99992   continue
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
        if (grvflg) then
          grvfor(1) = 3.0d0 * mass0 * gamma(1)
          grvfor(2) = 3.0d0 * mass0 * gamma(2)
          grvfor(3) = 3.0d0 * mass0 * gamma(3)
        endif
C
C
C.... ACCUMULATE THE SUBDOMAIN MASS
C
        if (masflg) then
          totmas = totmas + 3.0d0 * mass0
        endif
C
C     ------
C     RETURN
C     ------
C
        return
C
C     ---------------
C     ERROR TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF UNKNOWN TYPE OF CONSTITUTIVE LAW
C
100     continue
        write (*, *) '*** FATAL ERROR in Routine COMPMS ***'
        write (*, *) '*** The Type of Constitutive Law  ***'
        write (*, *) '*** Must Either be 0, 1, 2, or 3  ***'
        write (*, *) '*** Execution Terminated Here     ***'
        stop
C
C.....ERROR-MESSAGE IF THE ADDRESSING IS NOT CORRECT
C
200     continue
        write (*, *) '*** FATAL ERROR in Routine COMPMS ***'
        write (*, *) '*** The Address in Arrays [IDLAY] ***'
        write (*, *) '*** and [MTLAY] is Out-of-Bounds  ***'
        write (*, *) '*** Execution Terminated Here     ***'
        stop
C
C.....ERROR-MESSAGE IF THE MAXIMUM NUMBER OF LAYERS IS EXCEEDED
C
300     continue
        write (*, *) '*** FATAL ERROR in Routine COMPMS     ***'
        write (*, *) '*** Maximum Number of Layers Exceeded ***'
        write (*, *) '*** Boost Local Parameter [MAXLAYER]  ***'
        write (*, *) '*** Execution Terminated Here         ***'
        stop
C
C.....ERROR-MESSAGE IF THE TOTAL NUMBER OF LAYERS IS NOT CONSISTENT
C
400     continue
        write (*, *) '*** FATAL ERROR in Routine COMPMS      ***'
        write (*, *) '*** A Layer Number Exceeds the Total   ***'
        write (*, *) '*** Number of Layers Stored in [IDLAY] ***'
        write (*, *) '*** Execution Terminated Here          ***'
        stop
C
C.....ERROR-MESSAGE IF A SIDE HAS ZERO-LENGTH
C
500     continue
        write (*, *) '*** FATAL ERROR in routine COMPMS     ***'
        write (*, *) '*** The Side 1-2 has Zero Length      ***'
        write (*, *) '*** Check Coordinates and FE Topology ***'
        write (*, *) '*** Execution Terminated Here         ***'
        stop
C
C.....ERROR-MESSAGE IF THE AREA IS NEGATIVE OR ZERO
C
600     continue
        write (*, *) '*** FATAL ERROR in routine COMPMS     ***'
        write (*, *) '*** The Area is Negative or Zero      ***'
        write (*, *) '*** Check Coordinates and FE Topology ***'
        write (*, *) '*** ... Counterclock Nodal Numbering? ***'
        write (*, *) '*** Execution Terminated Here         ***'
        stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPMS"
C=======================C
