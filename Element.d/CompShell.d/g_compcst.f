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
      subroutine g_compcst(g_p_, e, g_e, ldg_e, thick, g_thick, ldg_thic
     *k, nu, coef, g_coef, ldg_coef, nlayer, idlayer, mtlayer, g_mtlayer
     *, ldg_mtlayer, x, g_x, ldg_x, y, g_y, ldg_y, z, d, g_d, ldg_d, typ
     *e, eframe, g_eframe, ldg_eframe, aframe, effect)
C=====================================================================C
C                                                                     C
C     Perform =   Assembles the 3 by 3 Constitutive Matrix According  C
C     ---------   to the Type of Constitutive Law Requested.          C
C                                                                     C
C                                                                     C
C     Input/Output =                                                  C
C     --------------                                                  C
C     E       <input>  Young modulus                                  C
C     THICK   <input>  thickness (assumed constant over the element)  C
C     NU      <input>  Poisson's ratio                                C
C     COEF    <input>  coefficients of the constitutive law           C
C     NLAYER  <input>  number of layers of the composite element      C
C     IDLAYER <input>  identificators for each layer                  C
C     MTLAYER <input>  material properties of each layer              C
C     X       <input>  triangular coordinates in local x-direction    C
C     Y       <input>  triangular coordinates in local y-direction    C
C     Z       <input>  triangular coordinates in local z-direction    C
C     D       <output> 3 by 3 constitutive matrix                     C
C     TYPE    <input>  type of constitutive law                       C
C     EFRAME  <input>  element level 3x3 frame                        C
C     AFRAME  <input>  arbitrary 3x3 frame of the constitutive law    C
C     EFFECT  <input>  type of matrix [D] requested                   C
C                                                                     C
C                                                                     C
C     Computations =                                                  C
C     --------------                                                  C
C     Four different matrices [D] are assembled according to whether  C
C     the flag [effect] is equal to "BB", "MM", "BM", or "MB":        C
C                                                                     C
C                             [ [D_mm]  [D_mb] ]                      C
C     [Constitutive_Matrix] = [                ]                      C
C            6 by 6           [ [D_bm]  [D_bb] ]                      C
C                                                                     C
C     where "b" and "m" stand for bending and membrane, respectively: C
C                                                                     C
C     [effect] = "BB"  =>  Assemble [D_bb] in [D]                     C
C     [effect] = "MM"  =>  Assemble [D_mm] in [D]                     C
C     [effect] = "BM"  =>  Assemble [D_bm] in [D]                     C
C     [effect] = "MB"  =>  Assemble [D_mb] in [D]                     C
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
C     The assembly of these matrices is performed according to the    C
C     type of constitutive law available and prescribed by [type].    C
C     In the following, the assembly of each one of the four matrices C
C     is briefly summarized.                                          C
C                                                                     C
C     1.  Constitutive Law of type-0:                                 C
C     - - - - - - - - - - - - - - - -                                 C
C     The material is assumed isotropic and known via the Young       C
C     modulus [E], the Poisson's ratio [nu] and thickness [thick]:    C
C                                                                     C
C     1.1  Pure Bending ("BB"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C              [ d_11  d_12   0   ]                                   C
C     [D_bb] = [ d_12  d_11   0   ]                                   C
C              [  0     0    d_33 ]                                   C
C                                                                     C
C     with:                                                           C
C                                                                     C
C              [E]*[thick]^3                                          C
C     [d_11] = -------------                                          C
C              12*(1-[nu]^2)                                          C
C                                                                     C
C              [nu]*[E]*[thick]^3                                     C
C     [d_12] = ------------------                                     C
C                12*(1-[nu]^2)                                        C
C                                                                     C
C              [E]*[thick]^3                                          C
C     [d_33] = -------------                                          C
C               24*(1+[nu])                                           C
C                                                                     C
C     1.2 Pure Membrane ("MM"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C              [ d_44  d_45   0   ]                                   C
C     [D_mm] = [ d_45  d_55   0   ]                                   C
C              [  0     0    d_66 ]                                   C
C                                                                     C
C     with:                                                           C
C                                                                     C
C              [E]*[thick]                                            C
C     [d_44] = -----------                                            C
C              (1-[nu]^2)                                             C
C                                                                     C
C              [nu]*[E]*[thick]                                       C
C     [d_45] = ----------------                                       C
C                 (1-[nu]^2)                                          C
C                                                                     C
C              [E]*[thick]                                            C
C     [d_66] = ------------                                           C
C               2*(1+[nu])                                            C
C                                                                     C
C     1.3 Coupling Bending-Membrane ("BM"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C     [D_bm] = zero (no coupling for isotropic material)              C
C                                                                     C
C     1.4 Coupling Membrane-Bending ("MB"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C     [D_mb] = zero (no coupling for isotropic material)              C
C                                                                     C
C     2.  Constitutive Law of type-1:                                 C
C     - - - - - - - - - - - - - - - -                                 C
C     The coefficients [d_ij] for i=1 to 6 and j=1 to 6 are given in  C
C     the output. They are stored in a vector of length 36 according  C
C     to the following convention:                                    C
C                                                                     C
C     [  d_11  d_12  d_13  d_14  d_15  d_16  ]                        C
C     [  d_12  d_22  d_23  d_24  d_25  d_26  ]                        C
C     [  d_13  d_23  d_33  d_34  d_35  d_36  ]                        C
C     [  d_14  d_24  d_33  d_44  d_45  d_46  ]                        C
C     [  d_15  d_25  d_33  d_44  d_55  d_56  ]                        C
C     [  d_16  d_26  d_33  d_44  d_55  d_66  ]                        C
C                                                                     C
C     is stored in the vector [coef] of size 36 at the following      C
C     location:                                                       C
C                                                                     C
C     [   01    02    03    04    05    06   ]                        C
C     [   07    08    09    10    11    12   ]                        C
C     [   13    14    15    16    17    18   ]                        C
C     [   19    20    21    22    23    24   ]                        C
C     [   25    26    27    28    29    30   ]                        C
C     [   31    32    33    34    35    36   ]                        C
C                                                                     C
C     Therefore, the entry on the ith row and jth column is stored at C
C     position number 6*(i-1)+j in the vector [coef].                 C
C                                                                     C
C     2.1  Pure Bending ("BB"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C              [ d_11 d_12 d_13 ]                    [ 22  23  24 ]   C
C     [D_bb] = [ d_12 d_22 d_23 ] found at positions [ 28  29  30 ]   C
C              [ d_13 d_23 d_33 ]                    [ 34  35  36 ]   C
C                                                                     C
C     2.2 Pure Membrane ("MM"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C              [ d_44 d_42 d_43 ]                    [ 01  02  03 ]   C
C     [D_mm] = [ d_45 d_55 d_56 ] found at positions [ 07  08  09 ]   C
C              [ d_43 d_56 d_66 ]                    [ 13  14  15 ]   C
C                                                                     C
C     2.3 Coupling Bending-Membrane ("BM"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C                                                                     C
C              [ d_14 d_15 d_16 ]                    [ 19  20  21 ]   C
C     [D_bm] = [ d_24 d_25 d_26 ] found at positions [ 25  26  27 ]   C
C              [ d_34 d_35 d_36 ]                    [ 31  32  33 ]   C
C                                                                     C
C     2.4 Coupling Membrane-Bending ("MB"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C                                                                     C
C              [ d_41 d_42 d_43 ]                    [ 04  05  06 ]   C
C     [D_mb] = [ d_51 d_52 d_53 ] found at positions [ 10  11  12 ]   C
C              [ d_61 d_62 d_63 ]                    [ 16  17  18 ]   C
C                                                                     C
C     3.  Constitutive Law of type-2:                                 C
C     - - - - - - - - - - - - - - - -                                 C
C     The material properties of each layer are known and integration C
C     through the thickness of the composite material is performed.   C
C     It is assumed that there is NO coupling between bending and     C
C     membrane effects (even though these terms are found non-zero    C
C     when numerical integration through the thickness is performed). C
C                                                                     C
C     3.1  Pure Bending ("BB"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C     3.2 Pure Membrane ("MM"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C     3.3 Coupling Bending-Membrane ("BM"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C     [D_bm] = zero (assumed)                                         C
C                                                                     C
C     3.4 Coupling Membrane-Bending ("MB"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C     [D_mb] = zero (assumed)                                         C
C                                                                     C
C     4.  Constitutive Law of type-3:                                 C
C     - - - - - - - - - - - - - - - -                                 C
C     The material properties of each layer are known and integration C
C     through the thickness of the composite material is performed.   C
C                                                                     C
C     4.1  Pure Bending ("BB"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C     4.2 Pure Membrane ("MM"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C     4.3 Coupling Bending-Membrane ("BM"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C                                                                     C
C     4.4 Coupling Membrane-Bending ("MB"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C                                                                     C
C                                                                     C
C     Caution =   It is assumed that the element has a constant       C
C     ---------   thickness so that no numerical interpolation is     C
C                 required. It is also assumed that the symmetry of   C
C                 the [D] matrix has been checked when its 36 entries C
C                 are inputed ([type]=1). (See routines "precmp.f"    C
C                 and "reacmp.f" in directory "Input".)               C
C                                                                     C
C=====================================================================C
C=Author  = Francois M. Hemez                                         C
C=Date    = June 9th, 1995                                            C
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
        integer type, nlayer, idlayer(5, nlayer)
        real*8 e, nu, thick, coef(36), d(3, 3)
        real*8 x(3), y(3), z(3), mtlayer(8, nlayer)
        real*8 eframe(3, 3), aframe(3, 3)
        character effect*2
C
C.....Local Variables
C
        integer i, j, k, ilayer, layernumber, irot
        real*8 zero, one, pi, twopi, intthick
        real*8 e1, e2, nu12, g12, mu1, mu2
        real*8 zsup, zinf, thetaf, thetad, theta
        real*8 s11, s12, s13, s22, s23, s33, dets
        real*8 q(3, 3), qbar(3, 3), t(3, 3), r(3)
        real*8 invt(3, 3), costheta, sintheta
        real*8 z0, qt1, qt2, qt3
        real*8 rotate(3, 3), rotated(3, 3), refvec(3)
        real*8 norm1, norm2, normref, proj1, proj2
        real*8 cosine1, cosine2, orifiber(3)
        logical pureben, puremem, cbenmem, cmemben
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_d, ldg_e, ldg_thick, ldg_coef, ldg_efram
     *e, ldg_mtlayer, ldg_x, ldg_y
        double precision d8_b, d3_b, d2_w, d1_w, d1_p, d7_b, d7_v, d6_b,
     * d2_v, d3_v
        double precision d4_v, d5_v, d5_b, d4_b, d2_b, g_d(ldg_d, 3, 3),
     * g_e(ldg_e), g_thick(ldg_thick), g_coef(ldg_coef, 36), g_rotate(g_
     *pmax_, 3, 3)
        double precision g_eframe(ldg_eframe, 3, 3), g_rotated(g_pmax_, 
     *3, 3), g_z0(g_pmax_), g_mtlayer(ldg_mtlayer, 8, nlayer), g_e1(g_pm
     *ax_), g_e2(g_pmax_), g_nu12(g_pmax_), g_g12(g_pmax_), g_mu1(g_pmax
     *_), g_mu2(g_pmax_)
        double precision g_thetaf(g_pmax_), g_zinf(g_pmax_), g_zsup(g_pm
     *ax_), g_orifiber(g_pmax_, 3), g_norm1(g_pmax_), g_norm2(g_pmax_), 
     *g_proj1(g_pmax_), g_proj2(g_pmax_), g_cosine1(g_pmax_), g_cosine2(
     *g_pmax_)
        double precision g_thetad(g_pmax_), g_d1_w(g_pmax_), g_theta(g_p
     *max_), g_s11(g_pmax_), g_s12(g_pmax_), g_s13(g_pmax_), g_s22(g_pma
     *x_), g_s23(g_pmax_), g_s33(g_pmax_), g_d2_w(g_pmax_)
        double precision g_dets(g_pmax_), g_q(g_pmax_, 3, 3), g_t(g_pmax
     *_, 3, 3), g_costheta(g_pmax_), g_sintheta(g_pmax_), g_invt(g_pmax_
     *, 3, 3), g_qbar(g_pmax_, 3, 3), g_qt1(g_pmax_), g_qt2(g_pmax_), g_
     *qt3(g_pmax_)
        double precision g_intthick(g_pmax_), g_x(ldg_x, 3), g_y(ldg_y, 
     *3)
        save g_q, g_t, g_costheta, g_sintheta, g_invt, g_qbar, g_qt1, g_
     *qt2, g_qt3, g_intthick
        save g_d1_w, g_theta, g_s11, g_s12, g_s13, g_s22, g_s23, g_s33, 
     *g_d2_w, g_dets
        save g_zinf, g_zsup, g_orifiber, g_norm1, g_norm2, g_proj1, g_pr
     *oj2, g_cosine1, g_cosine2, g_thetad
        save g_rotate, g_rotated, g_z0, g_e1, g_e2, g_nu12, g_g12, g_mu1
     *, g_mu2, g_thetaf
        intrinsic dble
        data zero /0.000000d+00/
        data one /1.000000d+00/
C
C.....INITIALIZE THE MAIN DIAGONAL OF REUTER'S MATRIX
C
        data r /1.000000d+00, 1.000000d+00, 2.000000d+00/
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
C.....CHECK IF THE ELEMENT IS A COMPOSITE SHELL
C
        if (type .eq. (-1)) then
          goto 100
        endif
C
C.....CHECK IF TYPE OF CONTITUTIVE LAW IS CORRECT
C
        if ((type .ne. 0) .and. (type .ne. 1) .and. (type .ne. 2) .and. 
     *(type .ne. 3)) then
          goto 200
        endif
C
C.....CHECK THE PHYSICAL EFFECT
C
        if ((effect .ne. 'BB') .and. (effect .ne. 'MM') .and. (effect .n
     *e. 'BM') .and. (effect .ne. 'MB')) then
          goto 300
        endif
C
C.....SET THE TYPE OF PHYSICAL EFFECT
C
        pureben = (effect .eq. 'BB')
        puremem = (effect .eq. 'MM')
        cbenmem = (effect .eq. 'BM')
        cmemben = (effect .eq. 'MB')
C
C.....CLEAR THE 3 BY 3 CONSTITUTIVE MATRIX
C
        do 99998 j = 1, 3
          do 99999 i = 1, 3
            do g_i_ = 1, g_p_
              g_d(g_i_, i, j) = 0.0d0
            enddo
            d(i, j) = zero
C--------
1002        continue
99999     continue
1001      continue
99998   continue
C
C     ------------------------------------------------
C     ISOTROPIC (ESSENTIALLY PLANE) MATERIAL
C     NO COUPLING BETWEEN BENDING AND MEMBRANE EFFECTS
C     ------------------------------------------------
C
        if (type .eq. 0) then
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
          if (pureben) then
C
            d3_v = thick * thick
            d4_v = d3_v * thick
            d2_b = 1.0d0 / (12.00d+00 * (one - nu * nu))
            d3_b = d2_b * d4_v
            d4_b = d2_b * e
            d5_b = d4_b * thick
            d6_b = d4_b * d3_v + d5_b * thick + d5_b * thick
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 1) = d6_b * g_thick(g_i_) + d3_b * g_e(g_i_)
            enddo
            d(1, 1) = e * d4_v / (12.00d+00 * (one - nu * nu))
C--------
            d2_v = nu * e
            d4_v = thick * thick
            d5_v = d4_v * thick
            d2_b = 1.0d0 / (12.00d+00 * (one - nu * nu))
            d4_b = d2_b * d2_v
            d5_b = d4_b * thick
            d6_b = d4_b * d4_v + d5_b * thick + d5_b * thick
            d7_b = d2_b * d5_v * nu
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 2) = d6_b * g_thick(g_i_) + d7_b * g_e(g_i_)
            enddo
            d(1, 2) = d2_v * d5_v / (12.00d+00 * (one - nu * nu))
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 3) = 0.0d0
            enddo
            d(1, 3) = zero
C--------
            d2_v = nu * e
            d4_v = thick * thick
            d5_v = d4_v * thick
            d2_b = 1.0d0 / (12.00d+00 * (one - nu * nu))
            d4_b = d2_b * d2_v
            d5_b = d4_b * thick
            d6_b = d4_b * d4_v + d5_b * thick + d5_b * thick
            d7_b = d2_b * d5_v * nu
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 1) = d6_b * g_thick(g_i_) + d7_b * g_e(g_i_)
            enddo
            d(2, 1) = d2_v * d5_v / (12.00d+00 * (one - nu * nu))
C--------
            d3_v = thick * thick
            d4_v = d3_v * thick
            d2_b = 1.0d0 / (12.00d+00 * (one - nu * nu))
            d3_b = d2_b * d4_v
            d4_b = d2_b * e
            d5_b = d4_b * thick
            d6_b = d4_b * d3_v + d5_b * thick + d5_b * thick
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 2) = d6_b * g_thick(g_i_) + d3_b * g_e(g_i_)
            enddo
            d(2, 2) = e * d4_v / (12.00d+00 * (one - nu * nu))
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 3) = 0.0d0
            enddo
            d(2, 3) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 1) = 0.0d0
            enddo
            d(3, 1) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 2) = 0.0d0
            enddo
            d(3, 2) = zero
C--------
            d3_v = thick * thick
            d4_v = d3_v * thick
            d2_b = 1.0d0 / (24.00d+00 * (one + nu))
            d3_b = d2_b * d4_v
            d4_b = d2_b * e
            d5_b = d4_b * thick
            d6_b = d4_b * d3_v + d5_b * thick + d5_b * thick
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 3) = d6_b * g_thick(g_i_) + d3_b * g_e(g_i_)
            enddo
            d(3, 3) = e * d4_v / (24.00d+00 * (one + nu))
C--------
C
          endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
          if (puremem) then
C
            d2_b = 1.0d0 / (one - nu * nu)
            d3_b = d2_b * thick
            d4_b = d2_b * e
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 1) = d4_b * g_thick(g_i_) + d3_b * g_e(g_i_)
            enddo
            d(1, 1) = e * thick / (one - nu * nu)
C--------
            d2_v = nu * e
            d2_b = 1.0d0 / (one - nu * nu)
            d4_b = d2_b * d2_v
            d5_b = d2_b * thick * nu
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 2) = d4_b * g_thick(g_i_) + d5_b * g_e(g_i_)
            enddo
            d(1, 2) = d2_v * thick / (one - nu * nu)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 3) = 0.0d0
            enddo
            d(1, 3) = zero
C--------
            d2_v = nu * e
            d2_b = 1.0d0 / (one - nu * nu)
            d4_b = d2_b * d2_v
            d5_b = d2_b * thick * nu
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 1) = d4_b * g_thick(g_i_) + d5_b * g_e(g_i_)
            enddo
            d(2, 1) = d2_v * thick / (one - nu * nu)
C--------
            d2_b = 1.0d0 / (one - nu * nu)
            d3_b = d2_b * thick
            d4_b = d2_b * e
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 2) = d4_b * g_thick(g_i_) + d3_b * g_e(g_i_)
            enddo
            d(2, 2) = e * thick / (one - nu * nu)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 3) = 0.0d0
            enddo
            d(2, 3) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 1) = 0.0d0
            enddo
            d(3, 1) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 2) = 0.0d0
            enddo
            d(3, 2) = zero
C--------
            d2_b = 1.0d0 / (2.00d+00 * (one + nu))
            d3_b = d2_b * thick
            d4_b = d2_b * e
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 3) = d4_b * g_thick(g_i_) + d3_b * g_e(g_i_)
            enddo
            d(3, 3) = e * thick / (2.00d+00 * (one + nu))
C--------
C
          endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C
          if (cbenmem) then
C
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 1) = 0.0d0
            enddo
            d(1, 1) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 2) = 0.0d0
            enddo
            d(1, 2) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 3) = 0.0d0
            enddo
            d(1, 3) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 1) = 0.0d0
            enddo
            d(2, 1) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 2) = 0.0d0
            enddo
            d(2, 2) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 3) = 0.0d0
            enddo
            d(2, 3) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 1) = 0.0d0
            enddo
            d(3, 1) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 2) = 0.0d0
            enddo
            d(3, 2) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 3) = 0.0d0
            enddo
            d(3, 3) = zero
C--------
C
          endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
          if (cmemben) then
C
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 1) = 0.0d0
            enddo
            d(1, 1) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 2) = 0.0d0
            enddo
            d(1, 2) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 3) = 0.0d0
            enddo
            d(1, 3) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 1) = 0.0d0
            enddo
            d(2, 1) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 2) = 0.0d0
            enddo
            d(2, 2) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 3) = 0.0d0
            enddo
            d(2, 3) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 1) = 0.0d0
            enddo
            d(3, 1) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 2) = 0.0d0
            enddo
            d(3, 2) = zero
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 3) = 0.0d0
            enddo
            d(3, 3) = zero
C--------
C
          endif
C
C.....END OF TREATMENT FOR ISOTROPIC MATERIAL
C
        endif
C
C     -------------------------------------------------
C     COMPOSITE MATERIAL WITH KNOWN CONSTITUTIVE MATRIX
C     -------------------------------------------------
C
        if (type .eq. 1) then
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
          if (pureben) then
C
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 1) = g_coef(g_i_, 22)
            enddo
            d(1, 1) = coef(22)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 2) = g_coef(g_i_, 23)
            enddo
            d(1, 2) = coef(23)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 3) = g_coef(g_i_, 24)
            enddo
            d(1, 3) = coef(24)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 1) = g_coef(g_i_, 28)
            enddo
            d(2, 1) = coef(28)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 2) = g_coef(g_i_, 29)
            enddo
            d(2, 2) = coef(29)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 3) = g_coef(g_i_, 30)
            enddo
            d(2, 3) = coef(30)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 1) = g_coef(g_i_, 34)
            enddo
            d(3, 1) = coef(34)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 2) = g_coef(g_i_, 35)
            enddo
            d(3, 2) = coef(35)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 3) = g_coef(g_i_, 36)
            enddo
            d(3, 3) = coef(36)
C--------
C
          endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
          if (puremem) then
C
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 1) = g_coef(g_i_, 1)
            enddo
            d(1, 1) = coef(1)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 2) = g_coef(g_i_, 2)
            enddo
            d(1, 2) = coef(2)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 3) = g_coef(g_i_, 3)
            enddo
            d(1, 3) = coef(3)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 1) = g_coef(g_i_, 7)
            enddo
            d(2, 1) = coef(7)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 2) = g_coef(g_i_, 8)
            enddo
            d(2, 2) = coef(8)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 3) = g_coef(g_i_, 9)
            enddo
            d(2, 3) = coef(9)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 1) = g_coef(g_i_, 13)
            enddo
            d(3, 1) = coef(13)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 2) = g_coef(g_i_, 14)
            enddo
            d(3, 2) = coef(14)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 3) = g_coef(g_i_, 15)
            enddo
            d(3, 3) = coef(15)
C--------
C
          endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C
          if (cbenmem) then
C
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 1) = g_coef(g_i_, 19)
            enddo
            d(1, 1) = coef(19)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 2) = g_coef(g_i_, 20)
            enddo
            d(1, 2) = coef(20)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 3) = g_coef(g_i_, 21)
            enddo
            d(1, 3) = coef(21)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 1) = g_coef(g_i_, 25)
            enddo
            d(2, 1) = coef(25)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 2) = g_coef(g_i_, 26)
            enddo
            d(2, 2) = coef(26)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 3) = g_coef(g_i_, 27)
            enddo
            d(2, 3) = coef(27)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 1) = g_coef(g_i_, 31)
            enddo
            d(3, 1) = coef(31)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 2) = g_coef(g_i_, 32)
            enddo
            d(3, 2) = coef(32)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 3) = g_coef(g_i_, 33)
            enddo
            d(3, 3) = coef(33)
C--------
C
          endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
          if (cmemben) then
C
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 1) = g_coef(g_i_, 4)
            enddo
            d(1, 1) = coef(4)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 2) = g_coef(g_i_, 5)
            enddo
            d(1, 2) = coef(5)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 1, 3) = g_coef(g_i_, 6)
            enddo
            d(1, 3) = coef(6)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 1) = g_coef(g_i_, 10)
            enddo
            d(2, 1) = coef(10)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 2) = g_coef(g_i_, 11)
            enddo
            d(2, 2) = coef(11)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 2, 3) = g_coef(g_i_, 12)
            enddo
            d(2, 3) = coef(12)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 1) = g_coef(g_i_, 16)
            enddo
            d(3, 1) = coef(16)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 2) = g_coef(g_i_, 17)
            enddo
            d(3, 2) = coef(17)
C--------
            do g_i_ = 1, g_p_
              g_d(g_i_, 3, 3) = g_coef(g_i_, 18)
            enddo
            d(3, 3) = coef(18)
C--------
C
          endif
C
C.....COMPUTE THE ROTATION MATRIX FROM THE ARBITRARY FRAME SYSTEM
C.....TO THE ELEMENT LEVEL FRAME SYSTEM ASSUMING THAT EACH ONE OF
C.....THESE FRAMES ARE ORTHONORMAL, RIGHT-HANDED SYSTEMS (THE I-TH
C.....VECTOR, I=1,2 OR 3, OF THE ARBITRARY FRAME SYSTEM IS ROTATED
C.....ONTO THE I-TH VECTOR OF THE ELEMENTAL FRAME SYSTEM)
C
          do 99996 j = 1, 3
            do 99997 i = 1, 3
              do g_i_ = 1, g_p_
                g_rotate(g_i_, i, j) = 0.0d0
              enddo
              rotate(i, j) = zero
C--------
1102          continue
99997       continue
1101        continue
99996     continue
C
          do 99993 k = 1, 3
            do 99994 j = 1, 3
              do 99995 i = 1, 3
                do g_i_ = 1, g_p_
                  g_rotate(g_i_, i, j) = aframe(j, k) * g_eframe(g_i_, i
     *, k) + g_rotate(g_i_, i, j)
                enddo
                rotate(i, j) = rotate(i, j) + eframe(i, k) * aframe(j, k
     *)
C--------
1105            continue
99995         continue
1104          continue
99994       continue
1103        continue
99993     continue
C
C.....ROTATE THE CONSTITUTIVE MATRIX IN THE ELEMENT LEVEL FRAME
C
          do 99991 j = 1, 3
            do 99992 i = 1, 3
              do g_i_ = 1, g_p_
                g_rotated(g_i_, i, j) = 0.0d0
              enddo
              rotated(i, j) = zero
C--------
1107          continue
99992       continue
1106        continue
99991     continue
C
          do 99988 j = 1, 3
            do 99989 k = 1, 3
              do 99990 i = 1, 3
                do g_i_ = 1, g_p_
                  g_rotated(g_i_, i, j) = d(i, k) * g_rotate(g_i_, k, j)
     * + rotate(k, j) * g_d(g_i_, i, k) + g_rotated(g_i_, i, j)
                enddo
                rotated(i, j) = rotated(i, j) + d(i, k) * rotate(k, j)
C--------
1110            continue
99990         continue
1109          continue
99989       continue
1108        continue
99988     continue
C
          do 99986 j = 1, 3
            do 99987 i = 1, 3
              do g_i_ = 1, g_p_
                g_d(g_i_, i, j) = 0.0d0
              enddo
              d(i, j) = zero
C--------
1112          continue
99987       continue
1111        continue
99986     continue
C
          do 99983 j = 1, 3
            do 99984 i = 1, 3
              do 99985 k = 1, 3
                do g_i_ = 1, g_p_
                  g_d(g_i_, i, j) = rotate(k, i) * g_rotated(g_i_, k, j)
     * + rotated(k, j) * g_rotate(g_i_, k, i) + g_d(g_i_, i, j)
                enddo
                d(i, j) = d(i, j) + rotate(k, i) * rotated(k, j)
C--------
1115            continue
99985         continue
1114          continue
99984       continue
1113        continue
99983     continue
C
C.....END OF TREATMENT FOR TYPE-1 CONSTITUTIVE LAW
C
        endif
C
C     ----------------------------------------------
C     COMPOSITE MATERIAL WITH KNOWN LAYER PROPERTIES
C         (COMPUTATIONS COMMUN TO TYPES 2 AND 3)
C     ----------------------------------------------
C
        if ((type .eq. 2) .or. (type .eq. 3)) then
C
C.....EXIT (SPEED UP COMPUTATIONS) IF COUPLING IS REQUESTED FOR TYPE-2
C.....BECAUSE TYPE-2 IS PRECISELY THE NO COUPLING CASE
C
          if (type .eq. 2) then
            if (cbenmem) then
              return
            endif
            if (cmemben) then
              return
            endif
          endif
C
C.....CALCULATE THE TOTAL HALF-HEIGHT OF THE LAYER
C
          do g_i_ = 1, g_p_
            g_z0(g_i_) = 0.0d0
          enddo
          z0 = zero
C--------
C
          do 99982 ilayer = 1, nlayer
            do g_i_ = 1, g_p_
              g_z0(g_i_) = g_mtlayer(g_i_, 7, ilayer) + g_z0(g_i_)
            enddo
            z0 = z0 + mtlayer(7, ilayer)
C--------
2001        continue
99982     continue
C
          do g_i_ = 1, g_p_
            g_z0(g_i_) = -0.50d+00 * g_z0(g_i_)
          enddo
          z0 = -0.50d+00 * z0
C--------
C
C.....LOOP ON LAYERS OF THE COMPOSITE SHELL ELEMENT
C
          do 99947 ilayer = 1, nlayer
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
            if ((layernumber .lt. 1) .or. (layernumber .gt. nlayer)) the
     *n
              goto 400
            endif
C
            if (e1 .le. zero) then
              goto 500
            endif
            if (e2 .le. zero) then
              goto 500
            endif
            if (nu12 .le. zero) then
              goto 500
            endif
            if (g12 .le. zero) then
              goto 500
            endif
            if (mu1 .lt. zero) then
              goto 500
            endif
            if (mu2 .lt. zero) then
              goto 500
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
              g_zinf(g_i_) = 0.0d0
            enddo
            zinf = zero
C--------
            do g_i_ = 1, g_p_
              g_zsup(g_i_) = 0.0d0
            enddo
            zsup = zero
C--------
C
            do 99981 i = 1, nlayer
              if (idlayer(3, i) .lt. layernumber) then
                do g_i_ = 1, g_p_
                  g_zinf(g_i_) = g_mtlayer(g_i_, 7, i) + g_zinf(g_i_)
                enddo
                zinf = zinf + mtlayer(7, i)
C--------
              endif
              if (idlayer(3, i) .le. layernumber) then
                do g_i_ = 1, g_p_
                  g_zsup(g_i_) = g_mtlayer(g_i_, 7, i) + g_zsup(g_i_)
                enddo
                zsup = zsup + mtlayer(7, i)
C--------
              endif
2003          continue
99981       continue
C
            do g_i_ = 1, g_p_
              g_zinf(g_i_) = g_zinf(g_i_) + g_z0(g_i_)
            enddo
            zinf = z0 + zinf
C--------
            do g_i_ = 1, g_p_
              g_zsup(g_i_) = g_zsup(g_i_) + g_z0(g_i_)
            enddo
            zsup = z0 + zsup
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
            do 99980 i = 1, 3
              d4_b = eframe(i, 1) + eframe(i, 1)
              do g_i_ = 1, g_p_
                g_norm1(g_i_) = d4_b * g_eframe(g_i_, i, 1) + g_norm1(g_
     *i_)
              enddo
              norm1 = norm1 + eframe(i, 1) * eframe(i, 1)
C--------
              d4_b = eframe(i, 2) + eframe(i, 2)
              do g_i_ = 1, g_p_
                g_norm2(g_i_) = d4_b * g_eframe(g_i_, i, 2) + g_norm2(g_
     *i_)
              enddo
              norm2 = norm2 + eframe(i, 2) * eframe(i, 2)
C--------
              normref = normref + refvec(i) * refvec(i)
              do g_i_ = 1, g_p_
                g_proj1(g_i_) = refvec(i) * g_eframe(g_i_, i, 1) + g_pro
     *j1(g_i_)
              enddo
              proj1 = proj1 + eframe(i, 1) * refvec(i)
C--------
              do g_i_ = 1, g_p_
                g_proj2(g_i_) = refvec(i) * g_eframe(g_i_, i, 2) + g_pro
     *j2(g_i_)
              enddo
              proj2 = proj2 + eframe(i, 2) * refvec(i)
C--------
2004          continue
99980       continue
C
            d2_v = sqrt(norm1)
            if ( norm1 .gt. 0.0d0 ) then
               d1_p = 1.0d0 / (2.0d0 *  d2_v)
            else
               call ehufDV (9,norm1, d2_v, d1_p,
     +'g_compcst.f',
     +1242)
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
     +'g_compcst.f',
     +1255)
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
                g_cosine1(g_i_) = d4_b * g_norm1(g_i_) + d2_b * g_proj1(
     *g_i_)
              enddo
              cosine1 = d4_v
C--------
              d3_v = norm2 * normref
              d4_v = proj2 / d3_v
              d2_b = 1.0d0 / d3_v
              d4_b = -d4_v / d3_v * normref
              do g_i_ = 1, g_p_
                g_cosine2(g_i_) = d4_b * g_norm2(g_i_) + d2_b * g_proj2(
     *g_i_)
              enddo
              cosine2 = d4_v
C--------
            endif
C
            do 99979 i = 1, 3
              do g_i_ = 1, g_p_
                g_orifiber(g_i_, i) = cosine2 * g_eframe(g_i_, i, 2) + e
     *frame(i, 2) * g_cosine2(g_i_) + cosine1 * g_eframe(g_i_, i, 1) + e
     *frame(i, 1) * g_cosine1(g_i_)
              enddo
              orifiber(i) = cosine1 * eframe(i, 1) + cosine2 * eframe(i,
     * 2)
C--------
2005          continue
99979       continue
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
            do 99978 i = 1, 3
              do g_i_ = 1, g_p_
                g_proj1(g_i_) = eframe(i, 1) * g_orifiber(g_i_, i) + ori
     *fiber(i) * g_eframe(g_i_, i, 1) + g_proj1(g_i_)
              enddo
              proj1 = proj1 + eframe(i, 1) * orifiber(i)
C--------
              do g_i_ = 1, g_p_
                g_proj2(g_i_) = eframe(i, 2) * g_orifiber(g_i_, i) + ori
     *fiber(i) * g_eframe(g_i_, i, 2) + g_proj2(g_i_)
              enddo
              proj2 = proj2 + eframe(i, 2) * orifiber(i)
C--------
2006          continue
99978       continue
C
            if (proj1 .eq. zero) then
              if (proj2 .eq. zero) then
                goto 700
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
              if (proj2 .eq. zero) then
                if (proj1 .eq. zero) then
                  goto 700
                endif
                if (proj1 .gt. zero) then
                  do g_i_ = 1, g_p_
                    g_thetad(g_i_) = 0.0d0
                  enddo
                  thetad = zero
C--------
                endif
                if (proj1 .lt. zero) then
                  do g_i_ = 1, g_p_
                    g_thetad(g_i_) = 0.0d0
                  enddo
                  thetad = pi
C--------
                endif
              else
                d3_v = proj2 / proj1
                d2_b = 1.0d0 / proj1
                d3_b = -d3_v / proj1
                do g_i_ = 1, g_p_
                  g_d1_w(g_i_) = d3_b * g_proj1(g_i_) + d2_b * g_proj2(g
     *_i_)
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
              goto 800
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
              g_d1_w(g_i_) = d6_b * g_s12(g_i_) + d8_b * g_s22(g_i_) + d
     *7_b * g_s11(g_i_) + d7_v * g_s33(g_i_)
            enddo
            d1_w = s33 * d7_v
            d4_v = s11 * s23
            d6_b = -s23 * s23
            d5_b = -d4_v + (-s23) * s11
            do g_i_ = 1, g_p_
              g_d2_w(g_i_) = d5_b * g_s23(g_i_) + d6_b * g_s11(g_i_) + g
     *_d1_w(g_i_)
            enddo
            d2_w = d1_w - d4_v * s23
            d4_v = s22 * s13
            d6_b = -s13 * s13
            d5_b = -d4_v + (-s13) * s22
            do g_i_ = 1, g_p_
              g_dets(g_i_) = d5_b * g_s13(g_i_) + d6_b * g_s22(g_i_) + g
     *_d2_w(g_i_)
            enddo
            dets = d2_w - d4_v * s13
C--------
            d3_v = 2.00d+00 * s12
            d5_v = d3_v * s13
            d7_b = s23 * d3_v
            d8_b = s23 * s13 * 2.00d+00
            do g_i_ = 1, g_p_
              g_dets(g_i_) = d5_v * g_s23(g_i_) + d7_b * g_s13(g_i_) + d
     *8_b * g_s12(g_i_) + g_dets(g_i_)
            enddo
            dets = dets + d5_v * s23
C--------
C
            if (dets .eq. zero) then
              goto 900
            endif
C
C.....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX WHICH RELATES
C.....THE STRAINS [e1], [e2] AND [e12] TO THE STRESSES [s1], [s2]
C.....AND [s12] IN THE COORDINATE SYSTEM {1;2} OF THE FIBER ORIENTATION
C
            do 99976 j = 1, 3
              do 99977 i = 1, 3
                do g_i_ = 1, g_p_
                  g_q(g_i_, i, j) = 0.0d0
                enddo
                q(i, j) = zero
C--------
2203            continue
99977         continue
2202          continue
99976       continue
C
            d4_b = -s23 + (-s23)
            do g_i_ = 1, g_p_
              g_q(g_i_, 1, 1) = d4_b * g_s23(g_i_) + s22 * g_s33(g_i_) +
     * s33 * g_s22(g_i_)
            enddo
            q(1, 1) = s22 * s33 - s23 * s23
C--------
            do g_i_ = 1, g_p_
              g_q(g_i_, 1, 2) = -s12 * g_s33(g_i_) + (-s33) * g_s12(g_i_
     *) + s13 * g_s23(g_i_) + s23 * g_s13(g_i_)
            enddo
            q(1, 2) = s13 * s23 - s12 * s33
C--------
            do g_i_ = 1, g_p_
              g_q(g_i_, 1, 3) = -s13 * g_s22(g_i_) + (-s22) * g_s13(g_i_
     *) + s12 * g_s23(g_i_) + s23 * g_s12(g_i_)
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
              g_q(g_i_, 2, 2) = d4_b * g_s13(g_i_) + s11 * g_s33(g_i_) +
     * s33 * g_s11(g_i_)
            enddo
            q(2, 2) = s11 * s33 - s13 * s13
C--------
            do g_i_ = 1, g_p_
              g_q(g_i_, 2, 3) = -s11 * g_s23(g_i_) + (-s23) * g_s11(g_i_
     *) + s12 * g_s13(g_i_) + s13 * g_s12(g_i_)
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
              g_q(g_i_, 3, 3) = d4_b * g_s12(g_i_) + s11 * g_s22(g_i_) +
     * s22 * g_s11(g_i_)
            enddo
            q(3, 3) = s11 * s22 - s12 * s12
C--------
C
            do 99974 j = 1, 3
              do 99975 i = 1, 3
                d3_v = q(i, j) / dets
                d2_b = 1.0d0 / dets
                d3_b = -d3_v / dets
                do g_i_ = 1, g_p_
                  g_q(g_i_, i, j) = d3_b * g_dets(g_i_) + d2_b * g_q(g_i
     *_, i, j)
                enddo
                q(i, j) = d3_v
C--------
2205            continue
99975         continue
2204          continue
99974       continue
C
C.....INITIALIZE THE ROTATION MATRIX FROM THE FIBER COORDINATE
C.....SYSTEM {1;2} TO THE ELEMENT TRIANGULAR SYSTEM {x;y}
C
            do 99972 j = 1, 3
              do 99973 i = 1, 3
                do g_i_ = 1, g_p_
                  g_t(g_i_, i, j) = 0.0d0
                enddo
                t(i, j) = zero
C--------
2207            continue
99973         continue
2206          continue
99972       continue
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
              g_t(g_i_, 1, 3) = d2_v * g_sintheta(g_i_) + d4_b * g_costh
     *eta(g_i_)
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
              g_t(g_i_, 2, 3) = d2_v * g_sintheta(g_i_) + d4_b * g_costh
     *eta(g_i_)
            enddo
            t(2, 3) = d2_v * sintheta
C--------
C
            do g_i_ = 1, g_p_
              g_t(g_i_, 3, 1) = -costheta * g_sintheta(g_i_) + (-sinthet
     *a) * g_costheta(g_i_)
            enddo
            t(3, 1) = -costheta * sintheta
C--------
            do g_i_ = 1, g_p_
              g_t(g_i_, 3, 2) = costheta * g_sintheta(g_i_) + sintheta *
     * g_costheta(g_i_)
            enddo
            t(3, 2) = costheta * sintheta
C--------
            d4_b = -sintheta + (-sintheta)
            d5_b = costheta + costheta
            do g_i_ = 1, g_p_
              g_t(g_i_, 3, 3) = d4_b * g_sintheta(g_i_) + d5_b * g_costh
     *eta(g_i_)
            enddo
            t(3, 3) = costheta * costheta - sintheta * sintheta
C--------
C
C.....COMPUTE THE INVERSE OF [T]:
C.....[invT] = inverse(diag[R]) * [T]^t * diag[R]
C
            do 99970 j = 1, 3
              do 99971 i = 1, 3
                do g_i_ = 1, g_p_
                  g_invt(g_i_, i, j) = 0.0d0
                enddo
                invt(i, j) = zero
C--------
2209            continue
99971         continue
2208          continue
99970       continue
C
            do 99968 j = 1, 3
              do 99969 i = 1, 3
                d2_b = r(j) / r(i)
                do g_i_ = 1, g_p_
                  g_invt(g_i_, i, j) = d2_b * g_t(g_i_, j, i)
                enddo
                invt(i, j) = t(j, i) * (r(j) / r(i))
C--------
2211            continue
99969         continue
2210          continue
99968       continue
C
C.....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX WHICH RELATES
C.....THE STRAINS [ex], [ey] AND [exy] TO THE STRESSES [sx], [sy]
C.....AND [sxy] IN THE TRIANGULAR COORDINATE SYSTEM {x;y}:
C.....[Qbar] = [invT] * [Q] * [invT]^t
C
            do 99966 j = 1, 3
              do 99967 i = 1, 3
                do g_i_ = 1, g_p_
                  g_qbar(g_i_, i, j) = 0.0d0
                enddo
                qbar(i, j) = zero
C--------
2213            continue
99967         continue
2212          continue
99966       continue
C
            do 99964 j = 1, 3
C
              do g_i_ = 1, g_p_
                g_d1_w(g_i_) = q(1, 2) * g_invt(g_i_, j, 2) + invt(j, 2)
     * * g_q(g_i_, 1, 2) + q(1, 1) * g_invt(g_i_, j, 1) + invt(j, 1) * g
     *_q(g_i_, 1, 1)
              enddo
              d1_w = q(1, 1) * invt(j, 1) + q(1, 2) * invt(j, 2)
              do g_i_ = 1, g_p_
                g_qt1(g_i_) = q(1, 3) * g_invt(g_i_, j, 3) + invt(j, 3) 
     ** g_q(g_i_, 1, 3) + g_d1_w(g_i_)
              enddo
              qt1 = d1_w + q(1, 3) * invt(j, 3)
C--------
              do g_i_ = 1, g_p_
                g_d1_w(g_i_) = q(2, 2) * g_invt(g_i_, j, 2) + invt(j, 2)
     * * g_q(g_i_, 2, 2) + q(2, 1) * g_invt(g_i_, j, 1) + invt(j, 1) * g
     *_q(g_i_, 2, 1)
              enddo
              d1_w = q(2, 1) * invt(j, 1) + q(2, 2) * invt(j, 2)
              do g_i_ = 1, g_p_
                g_qt2(g_i_) = q(2, 3) * g_invt(g_i_, j, 3) + invt(j, 3) 
     ** g_q(g_i_, 2, 3) + g_d1_w(g_i_)
              enddo
              qt2 = d1_w + q(2, 3) * invt(j, 3)
C--------
              do g_i_ = 1, g_p_
                g_d1_w(g_i_) = q(3, 2) * g_invt(g_i_, j, 2) + invt(j, 2)
     * * g_q(g_i_, 3, 2) + q(3, 1) * g_invt(g_i_, j, 1) + invt(j, 1) * g
     *_q(g_i_, 3, 1)
              enddo
              d1_w = q(3, 1) * invt(j, 1) + q(3, 2) * invt(j, 2)
              do g_i_ = 1, g_p_
                g_qt3(g_i_) = q(3, 3) * g_invt(g_i_, j, 3) + invt(j, 3) 
     ** g_q(g_i_, 3, 3) + g_d1_w(g_i_)
              enddo
              qt3 = d1_w + q(3, 3) * invt(j, 3)
C--------
C
              do 99965 i = 1, 3
                do g_i_ = 1, g_p_
                  g_d1_w(g_i_) = qt2 * g_invt(g_i_, i, 2) + invt(i, 2) *
     * g_qt2(g_i_) + qt1 * g_invt(g_i_, i, 1) + invt(i, 1) * g_qt1(g_i_)
                enddo
                d1_w = qt1 * invt(i, 1) + qt2 * invt(i, 2)
                do g_i_ = 1, g_p_
                  g_qbar(g_i_, i, j) = qt3 * g_invt(g_i_, i, 3) + invt(i
     *, 3) * g_qt3(g_i_) + g_d1_w(g_i_)
                enddo
                qbar(i, j) = d1_w + qt3 * invt(i, 3)
C--------
2215            continue
99965         continue
C
2214          continue
99964       continue
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
              if (pureben) then
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
                  g_intthick(g_i_) = d6_b * g_zinf(g_i_) + d2_b * g_d1_w
     *(g_i_)
                enddo
                intthick = (d1_w - d3_v * zinf) / 3.00d+00
C--------
C
                do 99962 j = 1, 3
                  do 99963 i = 1, 3
                    do g_i_ = 1, g_p_
                      g_d(g_i_, i, j) = qbar(i, j) * g_intthick(g_i_) + 
     *intthick * g_qbar(g_i_, i, j) + g_d(g_i_, i, j)
                    enddo
                    d(i, j) = d(i, j) + qbar(i, j) * intthick
C--------
3002                continue
99963             continue
3001              continue
99962           continue
C
              endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
              if (puremem) then
C
                do g_i_ = 1, g_p_
                  g_intthick(g_i_) = -g_zinf(g_i_) + g_zsup(g_i_)
                enddo
                intthick = zsup - zinf
C--------
C
                do 99960 j = 1, 3
                  do 99961 i = 1, 3
                    do g_i_ = 1, g_p_
                      g_d(g_i_, i, j) = qbar(i, j) * g_intthick(g_i_) + 
     *intthick * g_qbar(g_i_, i, j) + g_d(g_i_, i, j)
                    enddo
                    d(i, j) = d(i, j) + qbar(i, j) * intthick
C--------
3004                continue
99961             continue
3003              continue
99960           continue
C
              endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C
              if (cbenmem) then
C
                do 99958 j = 1, 3
                  do 99959 i = 1, 3
                    do g_i_ = 1, g_p_
                      g_d(g_i_, i, j) = 0.0d0
                    enddo
                    d(i, j) = zero
C--------
3006                continue
99959             continue
3005              continue
99958           continue
C
              endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
              if (cmemben) then
C
                do 99956 j = 1, 3
                  do 99957 i = 1, 3
                    do g_i_ = 1, g_p_
                      g_d(g_i_, i, j) = 0.0d0
                    enddo
                    d(i, j) = zero
C--------
3008                continue
99957             continue
3007              continue
99956           continue
C
              endif
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
              if (pureben) then
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
                  g_intthick(g_i_) = d6_b * g_zinf(g_i_) + d2_b * g_d1_w
     *(g_i_)
                enddo
                intthick = (d1_w - d3_v * zinf) / 3.00d+00
C--------
C
                do 99954 j = 1, 3
                  do 99955 i = 1, 3
                    do g_i_ = 1, g_p_
                      g_d(g_i_, i, j) = qbar(i, j) * g_intthick(g_i_) + 
     *intthick * g_qbar(g_i_, i, j) + g_d(g_i_, i, j)
                    enddo
                    d(i, j) = d(i, j) + qbar(i, j) * intthick
C--------
4002                continue
99955             continue
4001              continue
99954           continue
C
              endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
              if (puremem) then
C
                do g_i_ = 1, g_p_
                  g_intthick(g_i_) = -g_zinf(g_i_) + g_zsup(g_i_)
                enddo
                intthick = zsup - zinf
C--------
C
                do 99952 j = 1, 3
                  do 99953 i = 1, 3
                    do g_i_ = 1, g_p_
                      g_d(g_i_, i, j) = qbar(i, j) * g_intthick(g_i_) + 
     *intthick * g_qbar(g_i_, i, j) + g_d(g_i_, i, j)
                    enddo
                    d(i, j) = d(i, j) + qbar(i, j) * intthick
C--------
4004                continue
99953             continue
4003              continue
99952           continue
C
              endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C.....(THE TRANSPOSED OF MATRIX [Qbar] IS TAKEN TO GET BEND-MEMB COUPLING)
C
              if (cbenmem) then
C
                d5_b = -0.50d+00 * zinf + (-0.50d+00) * zinf
                d6_b = 0.50d+00 * zsup + 0.50d+00 * zsup
                do g_i_ = 1, g_p_
                  g_intthick(g_i_) = d5_b * g_zinf(g_i_) + d6_b * g_zsup
     *(g_i_)
                enddo
                intthick = 0.50d+00 * (zsup * zsup - zinf * zinf)
C--------
C
                do 99950 j = 1, 3
                  do 99951 i = 1, 3
                    do g_i_ = 1, g_p_
                      g_d(g_i_, i, j) = qbar(j, i) * g_intthick(g_i_) + 
     *intthick * g_qbar(g_i_, j, i) + g_d(g_i_, i, j)
                    enddo
                    d(i, j) = d(i, j) + qbar(j, i) * intthick
C--------
4006                continue
99951             continue
4005              continue
99950           continue
C
              endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
              if (cmemben) then
C
                d5_b = -0.50d+00 * zinf + (-0.50d+00) * zinf
                d6_b = 0.50d+00 * zsup + 0.50d+00 * zsup
                do g_i_ = 1, g_p_
                  g_intthick(g_i_) = d5_b * g_zinf(g_i_) + d6_b * g_zsup
     *(g_i_)
                enddo
                intthick = 0.50d+00 * (zsup * zsup - zinf * zinf)
C--------
C
                do 99948 j = 1, 3
                  do 99949 i = 1, 3
                    do g_i_ = 1, g_p_
                      g_d(g_i_, i, j) = qbar(i, j) * g_intthick(g_i_) + 
     *intthick * g_qbar(g_i_, i, j) + g_d(g_i_, i, j)
                    enddo
                    d(i, j) = d(i, j) + qbar(i, j) * intthick
C--------
4008                continue
99949             continue
4007              continue
99948           continue
C
              endif
C
C.....END OF TREATMENT FOR TYPE-3 CONSTITUTIVE LAW
C
            endif
C
C.....END OF TREATMENT FOR TYPES 2 AND 3 CONSTITUTIVE LAW
C
2002        continue
99947     continue
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
92      format ('*** Type Used is: ',5x,a2,9x,' ***')
C
C     ---------------
C     ERROR TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF THE ELEMENT IS NOT A COMPOSITE SHELL
C
100     continue
        write (*, *) '*** FATAL ERROR in Routine COMPCST        ***'
        write (*, *) '*** The Finite Element is not a Composite ***'
        write (*, *) '*** Shell Element: Inconsistency Detected ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE        ***'
        stop
C
C.....ERROR-MESSAGE IF THE CONSTITUTIVE LAW IS INCORRECT
C
200     continue
        write (*, *) '*** FATAL ERROR in Routine COMPCST       ***'
        write (*, *) '*** Wrong Type of Constitutive Law       ***'
        write (*, 91) type
        write (*, *) '*** Types Allowed are:                   ***'
        write (*, *) '*** 0 = isotropic material (default)     ***'
        write (*, *) '***     (no coupling bending/membrane)   ***'
        write (*, *) '*** 1 = given constitutive law           ***'
        write (*, *) '*** 2 = given layers properties          ***'
        write (*, *) '***     (no coupling bending/membrane)   ***'
        write (*, *) '*** 3 = given layers properties          ***'
        write (*, *) '***     (with coupling bending/membrane) ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE       ***'
        stop
C
C.....ERROR-MESSAGE IF THE PHYSICAL EFFECT IS INCORRECT
C
300     continue
        write (*, *) '*** FATAL ERROR in Routine COMPCST ***'
        write (*, *) '*** Wrong Type of Physical Effect  ***'
        write (*, 92) effect
        write (*, *) '*** Types Allowed are:             ***'
        write (*, *) '*** BB = pure bending              ***'
        write (*, *) '*** MM = pure membrane             ***'
        write (*, *) '*** BM = coupling bending-membrane ***'
        write (*, *) '*** MB = coupling membrane-bending ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE ***'
        stop
C
C.....ERROR-MESSAGE IF A LAYER IDENTIFICATION NUMBER IS NOT CORRECT
C
400     continue
        write (*, *) '*** FATAL ERROR in Routine COMPCST  ***'
        write (*, *) '*** The Local Layer Number is Not   ***'
        write (*, *) '*** Correct or Out of Bounds        ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE  ***'
        stop
C
C.....ERROR-MESSAGE IF A LAYER PROPERTY IS NOT CORRECT
C
500     continue
        write (*, *) '*** FATAL ERROR in Routine COMPCST     ***'
        write (*, *) '*** One of the Layer Material Property ***'
        write (*, *) '*** is Negative or Zero!               ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE     ***'
        stop
C
C.....ERROR-MESSAGE IF THE FIBER ANGLE IS OUT-OF-BOUNDS
C
600     continue
        write (*, *) '*** FATAL ERROR in Routine COMPCST    ***'
        write (*, *) '*** The Angle From the Reference      ***'
        write (*, *) '*** Direction to the Direction of     ***'
        write (*, *) '*** Fibers is Out-of-Bounds: it Must  ***'
        write (*, *) '*** be Within the Range 0-360 Degrees ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE    ***'
        stop
C
C.....ERROR-MESSAGE IF THE REFERENCE ORIENTATION IS BUGGY
C
700     continue
        write (*, *) '*** FATAL ERROR in Routine COMPCST   ***'
        write (*, *) '*** The Reference Orientation Vector ***'
        write (*, *) '*** is Parallel to the Two Inplane   ***'
        write (*, *) '*** and Orthogonal Local Frames!     ***'
        write (*, *) '*** STOP ALL TREATMENTS RIGHT HERE   ***'
        stop
C
C.....ERROR-MESSAGE IF THE ORIENTATION ANGLE IS OUT-OF-BOUNDS
C
800     continue
        write (*, *) '*** FATAL ERROR in Routine COMPCST    ***'
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
900     continue
        write (*, *) '*** FATAL ERROR in Routine COMPCST    ***'
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
C=end of routine "COMPCST"
C========================C
