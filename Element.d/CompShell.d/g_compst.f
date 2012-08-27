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
C===================================================================C
      subroutine gxcompst(e,g_e,elm,h,g_h,estiff,g_estiff,medof,nu,
     *                    x,g_x,y,g_y,z,g_z,nttco,nttly,ncmpfr, 
     *	                  cmpco,g_cmpco,idlay,mtlay,g_mtlay,cmpfr,
     *                    iatt,ctyp,catt,cfrm,flag)
C===================================================================C
C                                                                   C
C     Perform =    This subroutine will form the element stiffness  C
C     ---------    matrix for the 3D 3-node composite shell element C
C                  derived from the ANDES formulation (ANS shell    C
C                  element with Assumed Quadratic Rotations).       C
C                                                                   C
C                                                                   C
C     Inputs/Outputs =                                              C
C     ----------------                                              C
C     E       <input>  Young modulus                                C
C     ELM     <input>  finite element number                        C
C     H       <input>  element thickness (assumed constant)         C
C     ESTIFF  <output> element stiffness matrix                     C
C     MEDOF   <input>  maximum number of DOFs per FE                C
C     NU      <input>  Poisson's ratio                              C
C     X       <input>  nodal coordinates in the X-direction         C
C     Y       <input>  nodal coordinates in the Y-direction         C
C     Z       <input>  nodal coordinates in the Z-direction         C
C     NTTCO   <input>  number of attributes for laws of type-1      C
C     NTTLY   <input>  total number of composite layers             C
C     NCMPFR  <input>  number of frames for composite shells        C
C     CMPCO   <input>  constitutive coefficients of the attributes  C
C     IDLAY   <input>  identificators for the composite layers      C
C     MTLAY   <input>  material properties of the composite layers  C
C     CMPFR   <input>  storage of frames for composite shells       C
C     IATT    <input>  attribute number of the element              C
C     CTYP    <input>  composite attribute number                   C
C     CATT    <input>  starting address for the constitutive law    C
C     CFRM    <input>  frame number for definition of the shell     C
C     FLAG    <input>  integer specifying whether to return         C
C                      transformed element stiffness matrix or      C
C                      global element stiffness matrix              C
C                                                                   C
C                                                                   C
C     Computations =                                                C
C     --------------                                                C
C     This subroutine evaluates the stiffness matrix for the 18     C
C     degrees of freedom 3-node composite triangle obtained as a    C
C     combination of the Assumed Quadratic Rotations bending        C
C     triangle plus the membrane with driling degrees of freedom    C
C     developed by Militello, Felippa et al. For documentation,     C
C     see Carmelo Militello's doctoral dissertation, pp112-113.     C
C                                                                   C
C     The original version of the ANS shell element has been        C
C     generalized here to the case of a composite shell element     C
C     with complete bending-membrane coupling. Four types of        C
C     composite laws are available via the input:                   C
C        type-0: isotropic element                                  C
C        type-1: constitutive coefficients are given                C
C        type-2: properties of each layer are given and no coupling C
C                   is assumed between bending and membrane         C
C        type-3: properties of each layer are given and coupling    C
C                   between bending and membrane is assumed         C
C                                                                   C
C     The stiffness matrix [K] is assembled as the combination of   C
C     the basic stiffness and the higher order stiffness matrices.  C
C     There are two of such matrices for each physical effect.      C
C     In the most general case, the stiffness [K] is given as:      C
C        [K] = [K_basic       _Bending _Bending ]                   C
C            + [K_basic       _Membrane_Membrane]                   C
C            + [K_basic       _Bending _Membrane]                   C
C            + [K_basic       _Membrane_Bending ]                   C
C            + [K_higher_order_Bending _Bending ]                   C
C            + [K_higher_order_Membrane_Membrane]                   C
C            + [K_higher_order_Bending _Membrane]                   C
C            + [K_higher_order_Membrane_Bending ]                   C
C     In general, a particular stiffness matrix may be expressed as C
C     the product of the following quantities:                      C
C         [K] = sum_{i=1,2,3} [B1_i]^T [C] [B2_i]                   C
C     where [i] represents the numerical integration index; [B1_i]  C
C     and [B2_i] are either the moment-curvature or force-strain    C
C     matrices (for bending or membrane effect, respectively); and  C
C     [C] represents the constitutive law. See the following        C
C     routines for details:                                         C
C     "compbBB" assembly of [K_basic       _Bending _Bending ]      C
C     "compbMM" assembly of [K_basic       _Membrane_Membrane]      C
C     "compbBM" assembly of [K_basic       _Bending _Membrane]      C
C     "compbMB" assembly of [K_basic       _Membrane_Bending ]      C
C     "comphBB" assembly of [K_higher_order_Bending _Bending ]      C
C     "comphMM" assembly of [K_higher_order_Membrane_Membrane]      C
C     "comphBM" assembly of [K_higher_order_Bending _Membrane]      C
C     "comphMB" assembly of [K_higher_order_Membrane_Bending ]      C
C                                                                   C
C                                                                   C
C     Caution =                                                     C
C     ---------                                                     C
C     The finite element is assumed to have a constant thickness    C
C     so that no numerical interpolation is required.  The array    C
C     storing thicknesses at each nodal points, [H], has a size     C
C     equal to 3 here. Also, The outputed element stiffness matrix  C
C     is a 18 by 18 block stored in the left upper corner of the    C
C     output storage [ESTIFF]. Finally, the maximum number of       C
C     composite layers per finite element is assumed to be equal to C
C     one thousand (upgrade local parameter [MAXLAYER] if not large C
C     enough).                                                      C
C                                                                   C
C                                                                   C
C     Output =   no output.                                         C
C     --------                                                      C
C                                                                   C
C===================================================================C
C=Author  = Francois M. Hemez                                       C
C=Date    = June 9, 1995                                            C
C=Version = 2.0                                                     C
C=Comment =                                                         C
C===================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....GLOBAL VARIABLES
C
        integer medof, elm, idlay(5, nttly)
        integer nttco, nttly, ncmpfr, flag
        integer iatt, ctyp, catt, cfrm
C
        real*8 nu, h(3), e, estiff(medof, medof)
        real*8 x(3), y(3), z(3), cmpco(36, nttco)
        real*8 mtlay(9, nttly), cmpfr(9, ncmpfr)
C
C.....LOCAL DIMENSION OF THE STIFFNESS MATRIX
C.....(3 NODES AND 6 DOFS PER NODE = 18 DOFS TOTAL)
C
        integer ndof
        parameter (ndof = 18)
C
C.....LOCAL MAXIMUM NUMBER OF LAYERS PER COMPOSITE ELEMENT
C.....(SET TO 1,000 HERE)
C
        integer maxlayer
        parameter (maxlayer = 1000)
C
C.....LOCAL VARIABLES
C
        integer i, j, nlayer, ilayer
        integer rowb(9), colb(9), rowm(9), colm(9)
        integer idcmp23(5, maxlayer)
C
        real*8 rot(6, 6), xlp(3), ylp(3), zlp(3)
        real*8 zero, one, thick, cstcoef(36)
        real*8 cstbb(3, 3), cstmm(3, 3), onehalf
        real*8 cstbm(3, 3), cstmb(3, 3), point32
        real*8 kbbb(ndof, ndof), kbmm(ndof, ndof)
        real*8 kbbm(ndof, ndof), kbmb(ndof, ndof)
        real*8 khbb(ndof, ndof), khmm(ndof, ndof)
        real*8 khbm(ndof, ndof), khmb(ndof, ndof)
        real*8 lb(9, 3), pb(9, 3), mtcmp23(8, maxlayer)
        real*8 lh1(9, 3), lh2(9, 3), lh3(9, 3)
        real*8 ph1(9, 3), ph2(9, 3), ph3(9, 3)
        real*8 aframe(3, 3), eframe(3, 3)
        logical fastcal
C
C     ----
C     DATA
C     ----
C
       integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_estiff, ldg_h, ldg_cmpco, ldg_mtlay, ldg
     *_e, ldg_x, ldg_y, ldg_z
c
c manually inserted - begin	     
c
        parameter (g_p_=1,ldg_e=1,ldg_h=1,ldg_estiff=1,
     *             ldg_x=1,ldg_y=1,ldg_z=1,ldg_cmpco=1,ldg_mtlay=1)
c
c manually inserted - end      
c
        double precision d1_w, g_estiff(ldg_estiff, medof, medof), g_cst
     *bb(g_pmax_, 3, 3), g_cstmm(g_pmax_, 3, 3), g_cstbm(g_pmax_, 3, 3),
     * g_cstmb(g_pmax_, 3, 3), g_kbbb(g_pmax_, ndof, ndof), g_kbmm(g_pma
     *x_, ndof, ndof), g_kbbm(g_pmax_, ndof, ndof), g_kbmb(g_pmax_, ndof
     *, ndof)
        double precision g_khbb(g_pmax_, ndof, ndof), g_khmm(g_pmax_, nd
     *of, ndof), g_khbm(g_pmax_, ndof, ndof), g_khmb(g_pmax_, ndof, ndof
     *), g_cstcoef(g_pmax_, 36), g_xlp(g_pmax_, 3), g_ylp(g_pmax_, 3), g
     *_lb(g_pmax_, 9, 3), g_pb(g_pmax_, 9, 3), g_lh1(g_pmax_, 9, 3)
        double precision g_lh2(g_pmax_, 9, 3), g_lh3(g_pmax_, 9, 3), g_p
     *h1(g_pmax_, 9, 3), g_ph2(g_pmax_, 9, 3), g_ph3(g_pmax_, 9, 3), g_m
     *tcmp23(g_pmax_, 8, maxlayer), g_eframe(g_pmax_, 3, 3), g_thick(g_p
     *max_), g_h(ldg_h, 3), g_cmpco(ldg_cmpco, 36, nttco)
        double precision g_mtlay(ldg_mtlay, 9, nttly), g_rot(g_pmax_, 6,
     * 6), g_d1_w(g_pmax_), g_e(ldg_e), g_x(ldg_x, 3), g_y(ldg_y, 3), g_
     *z(ldg_z, 3)
        save g_ph1, g_ph2, g_ph3, g_mtcmp23, g_eframe, g_thick, g_rot, g
     *_d1_w
        save g_khbm, g_khmb, g_cstcoef, g_xlp, g_ylp, g_lb, g_pb, g_lh1,
     * g_lh2, g_lh3
        save g_cstbb, g_cstmm, g_cstbm, g_cstmb, g_kbbb, g_kbmm, g_kbbm,
     * g_kbmb, g_khbb, g_khmm
        external g_trirotation
        external g_comphmb
        external g_comphbm
        external g_comphmm
        external g_comphbb
        external g_compbmb
        external g_compbbm
        external g_compbmm
        external g_compbbb
        external g_compcst
        external g_compcrd
        data zero /0.000000d+00/
        data one /1.000000d+00/
        data onehalf /1.500000d+00/
        data point32 /0.320000d+00/
C
C.....INITIALIZE THE LOGICAL FOR ENFORCING FASTER COMPUTATIONS
C
        data fastcal /.true./
C
C     -----
C     LOGIC
C     -----
C
C.....CHECK DIMENSION OF THE ELEMENTAL STIFFNESS MATRIX
C
        integer g_ehfid
        data g_ehfid /0/
C

C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        if (medof .ne. ndof) then
          goto 100
        endif
C
C.....CLEAR THE OUTPUT ELEMENT STIFFNESS MATRIX
C
        do 99998 j = 1, medof
          do 99999 i = 1, medof
            do g_i_ = 1, g_p_
              g_estiff(g_i_, i, j) = 0.0d0
            enddo
            estiff(i, j) = zero
C--------
1002        continue
99999     continue
1001      continue
99998   continue
C
C.....CLEAR THE LOCAL CONSTITUTIVE MATRICES
C
        do 99993 j = 1, 3
          do 99997 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstbb(g_i_, i, j) = 0.0d0
            enddo
            cstbb(i, j) = zero
C--------
1004        continue
99997     continue
          do 99996 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstmm(g_i_, i, j) = 0.0d0
            enddo
            cstmm(i, j) = zero
C--------
1005        continue
99996     continue
          do 99995 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstbm(g_i_, i, j) = 0.0d0
            enddo
            cstbm(i, j) = zero
C--------
1006        continue
99995     continue
          do 99994 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstmb(g_i_, i, j) = 0.0d0
            enddo
            cstmb(i, j) = zero
C--------
1007        continue
99994     continue
1003      continue
99993   continue
C
C.....CLEAR THE LOCAL CONTRIBUTIONS TO THE STIFFNESS MATRIX
C
        do 99992 j = 1, ndof
C         do 1009 i=1,ndof
          do g_i_ = 1, g_p_
            g_kbbb(g_i_, i, j) = 0.0d0
          enddo
          kbbb(i, j) = zero
C--------
C 1009    continue
C         do 1010 i=1,ndof
          do g_i_ = 1, g_p_
            g_kbmm(g_i_, i, j) = 0.0d0
          enddo
          kbmm(i, j) = zero
C--------
C 1010    continue
C         do 1011 i=1,ndof
          do g_i_ = 1, g_p_
            g_kbbm(g_i_, i, j) = 0.0d0
          enddo
          kbbm(i, j) = zero
C--------
C 1011    continue
C         do 1012 i=1,ndof
          do g_i_ = 1, g_p_
            g_kbmb(g_i_, i, j) = 0.0d0
          enddo
          kbmb(i, j) = zero
C--------
C 1012    continue
C         do 1013 i=1,ndof
          do g_i_ = 1, g_p_
            g_khbb(g_i_, i, j) = 0.0d0
          enddo
          khbb(i, j) = zero
C--------
C 1013    continue
C         do 1014 i=1,ndof
          do g_i_ = 1, g_p_
            g_khmm(g_i_, i, j) = 0.0d0
          enddo
          khmm(i, j) = zero
C--------
C 1014    continue
C         do 1015 i=1,ndof
          do g_i_ = 1, g_p_
            g_khbm(g_i_, i, j) = 0.0d0
          enddo
          khbm(i, j) = zero
C--------
C 1015    continue
C         do 1016 i=1,ndof
          do g_i_ = 1, g_p_
            g_khmb(g_i_, i, j) = 0.0d0
          enddo
          khmb(i, j) = zero
C--------
C 1016    continue
1008      continue
99992   continue
C
C.....CLEAR THE LOCAL CONSTITUTIVE COEFFICIENT ARRAY
C
        do 99991 i = 1, 36
          do g_i_ = 1, g_p_
            g_cstcoef(g_i_, i) = 0.0d0
          enddo
          cstcoef(i) = zero
C--------
1017      continue
99991   continue
C
C.....CLEAR THE TRIANGULAR COORDINATES
C
        do 99990 i = 1, 3
          do g_i_ = 1, g_p_
            g_xlp(g_i_, i) = 0.0d0
          enddo
          xlp(i) = zero
C--------
1018      continue
99990   continue
C
        do 99989 i = 1, 3
          do g_i_ = 1, g_p_
            g_ylp(g_i_, i) = 0.0d0
          enddo
          ylp(i) = zero
C--------
1019      continue
99989   continue
C
        do 99988 i = 1, 3
          zlp(i) = zero
1020      continue
99988   continue
C
C.....CLEAR THE DEGREE OF FREEDOM POINTERS
C
C      do 1023 i=1,9
C         rowb(i) = 0
C 1023 continue
C
C      do 1024 i=1,9
C         colb(i) = 0
C 1024 continue
C
C      do 1025 i=1,9
C         rowm(i) = 0
C 1025 continue
C
C      do 1026 i=1,9
C         colm(i) = 0
C 1026 continue
C
C.....CLEAR THE LOCAL INTEGRATED STRAIN-TO-DISPLACEMENT AND
C.....CURVATURE-TO-DISPLACEMENT MATRICES
C
        do 99979 j = 1, 3
          do 99987 i = 1, 9
            do g_i_ = 1, g_p_
              g_lb(g_i_, i, j) = 0.0d0
            enddo
            lb(i, j) = zero
C--------
1028        continue
99987     continue
          do 99986 i = 1, 9
            do g_i_ = 1, g_p_
              g_pb(g_i_, i, j) = 0.0d0
            enddo
            pb(i, j) = zero
C--------
1029        continue
99986     continue
          do 99985 i = 1, 9
            do g_i_ = 1, g_p_
              g_lh1(g_i_, i, j) = 0.0d0
            enddo
            lh1(i, j) = zero
C--------
1030        continue
99985     continue
          do 99984 i = 1, 9
            do g_i_ = 1, g_p_
              g_lh2(g_i_, i, j) = 0.0d0
            enddo
            lh2(i, j) = zero
C--------
1031        continue
99984     continue
          do 99983 i = 1, 9
            do g_i_ = 1, g_p_
              g_lh3(g_i_, i, j) = 0.0d0
            enddo
            lh3(i, j) = zero
C--------
1032        continue
99983     continue
          do 99982 i = 1, 9
            do g_i_ = 1, g_p_
              g_ph1(g_i_, i, j) = 0.0d0
            enddo
            ph1(i, j) = zero
C--------
1033        continue
99982     continue
          do 99981 i = 1, 9
            do g_i_ = 1, g_p_
              g_ph2(g_i_, i, j) = 0.0d0
            enddo
            ph2(i, j) = zero
C--------
1034        continue
99981     continue
          do 99980 i = 1, 9
            do g_i_ = 1, g_p_
              g_ph3(g_i_, i, j) = 0.0d0
            enddo
            ph3(i, j) = zero
C--------
1035        continue
99980     continue
1027      continue
99979   continue
C
C.....CLEAR THE NUMBER OF LAYERS OF THE COMPOSITE ELEMENT
C
        nlayer = 0
C
C.....CLEAR THE LOCAL STORAGES FOR LAYER PROPERTIES
C
        do 99976 j = 1, maxlayer
          do 99978 i = 1, 5
            idcmp23(i, j) = 0
1037        continue
99978     continue
          do 99977 i = 1, 8
            do g_i_ = 1, g_p_
              g_mtcmp23(g_i_, i, j) = 0.0d0
            enddo
            mtcmp23(i, j) = zero
C--------
1038        continue
99977     continue
1036      continue
99976   continue
C
C.....CLEAR THE ARBITRARY FRAME OF THE COMPOSITE ELEMENT
C
        do 99974 j = 1, 3
          do 99975 i = 1, 3
            aframe(i, j) = zero
1040        continue
99975     continue
1039      continue
99974   continue
C
C.....CLEAR THE ELEMENT LEVEL FRAME
C
        do 99972 j = 1, 3
          do 99973 i = 1, 3
            do g_i_ = 1, g_p_
              g_eframe(g_i_, i, j) = 0.0d0
            enddo
            eframe(i, j) = zero
C--------
1042        continue
99973     continue
1041      continue
99972   continue
C
C.....INITIALIZE THE ELEMENT'S CONSTANT THICKNESS
C
        do g_i_ = 1, g_p_
          g_thick(g_i_) = g_h(g_i_, 1)
        enddo
        thick = h(1)
C--------
C
C.....CHECK THE TYPE OF CONSTITUTIVE LAW
C
        if ((ctyp .ne. 0) .and. (ctyp .ne. 1) .and. (ctyp .ne. 2) .and. 
     *(ctyp .ne. 3)) then
          goto 200
        endif
C
C.....CHECK THE ADDRESSING IN STORAGE [CMPCO]
C
        if (ctyp .eq. 1) then
          if ((catt .lt. 1) .or. (catt .gt. nttco)) then
            goto 300
          endif
        endif
C
C.....CHECK THE ADDRESSING IN ARRAYS [IDLAY] AND [MTLAY]
C
        if ((ctyp .eq. 2) .or. (ctyp .eq. 3)) then
          if ((catt .lt. 1) .or. (catt .gt. nttly)) then
            goto 400
          endif
        endif
C
C.....CHECK THE ADDRESSING IN ARRAY [CMPFR]
C
        if ((ctyp .eq. 1) .or. (ctyp .eq. 2) .or. (ctyp .eq. 3)) then
          if ((cfrm .lt. 0) .or. (cfrm .gt. ncmpfr)) then
            goto 500
          endif
        endif
C
C.....INITIALIZE THE CONSTITUTIVE COEFFICIENTS IN CASE OF A TYPE-1 LAW
C
        if (ctyp .eq. 1) then
          do 99971 i = 1, 36
            do g_i_ = 1, g_p_
              g_cstcoef(g_i_, i) = g_cmpco(g_i_, i, catt)
            enddo
            cstcoef(i) = cmpco(i, catt)
C--------
2001        continue
99971     continue
        endif
C
C.....INITIALIZE THE LAYER PROPERTIES FOR TYPE-2 AND TYPE-3 LAWS
C
        if ((ctyp .eq. 2) .or. (ctyp .eq. 3)) then
          nlayer = idlay(2, catt)
          if (nlayer .gt. maxlayer) then
            goto 600
          endif
          do 99970 i = catt, (catt + nlayer - 1)
            ilayer = idlay(3, i)
            if (ilayer .gt. nlayer) then
              goto 700
            endif
            idcmp23(1, ilayer) = idlay(1, i)
            idcmp23(2, ilayer) = idlay(2, i)
            idcmp23(3, ilayer) = idlay(3, i)
            idcmp23(4, ilayer) = idlay(4, i)
            idcmp23(5, ilayer) = idlay(5, i)
            do g_i_ = 1, g_p_
              g_mtcmp23(g_i_, 1, ilayer) = g_mtlay(g_i_, 1, i)
            enddo
            mtcmp23(1, ilayer) = mtlay(1, i)
C--------
            do g_i_ = 1, g_p_
              g_mtcmp23(g_i_, 2, ilayer) = g_mtlay(g_i_, 2, i)
            enddo
            mtcmp23(2, ilayer) = mtlay(2, i)
C--------
            do g_i_ = 1, g_p_
              g_mtcmp23(g_i_, 3, ilayer) = g_mtlay(g_i_, 3, i)
            enddo
            mtcmp23(3, ilayer) = mtlay(3, i)
C--------
            do g_i_ = 1, g_p_
              g_mtcmp23(g_i_, 4, ilayer) = g_mtlay(g_i_, 4, i)
            enddo
            mtcmp23(4, ilayer) = mtlay(4, i)
C--------
            do g_i_ = 1, g_p_
              g_mtcmp23(g_i_, 5, ilayer) = g_mtlay(g_i_, 5, i)
            enddo
            mtcmp23(5, ilayer) = mtlay(5, i)
C--------
            do g_i_ = 1, g_p_
              g_mtcmp23(g_i_, 6, ilayer) = g_mtlay(g_i_, 6, i)
            enddo
            mtcmp23(6, ilayer) = mtlay(6, i)
C--------
            do g_i_ = 1, g_p_
              g_mtcmp23(g_i_, 7, ilayer) = g_mtlay(g_i_, 8, i)
            enddo
            mtcmp23(7, ilayer) = mtlay(8, i)
C--------
            do g_i_ = 1, g_p_
              g_mtcmp23(g_i_, 8, ilayer) = g_mtlay(g_i_, 9, i)
            enddo
            mtcmp23(8, ilayer) = mtlay(9, i)
C--------
2002        continue
99970     continue
        endif
C
C.....CHECK CONSISTENCY OF BASIC PARAMETERS FOR TYPES 2 AND 3
C
        if ((ctyp .eq. 2) .or. (ctyp .eq. 3)) then
          do 99969 i = 1, nlayer
            if (idcmp23(1, i) .ne. idcmp23(1, 1)) then
              goto 800
            endif
            if (idcmp23(2, i) .ne. nlayer) then
              goto 800
            endif
2003        continue
99969     continue
          do 99968 i = (nlayer + 1), maxlayer
            if (idcmp23(1, i) .ne. 0) then
              goto 800
            endif
            if (idcmp23(2, i) .ne. 0) then
              goto 800
            endif
2004        continue
99968     continue
        endif
C
C.....GET THE ARBITRARY FRAME FOR DEFINITION OF THE CONSTITUTIVE LAW
C
        if (cfrm .eq. 0) then
C
C.....INITIALIZE WITH THE IDENTITY IF THE FRAME NUMBER IS ZERO
C
          aframe(1, 1) = one
          aframe(2, 1) = zero
          aframe(3, 1) = zero
          aframe(1, 2) = zero
          aframe(2, 2) = one
          aframe(3, 2) = zero
          aframe(1, 3) = zero
          aframe(2, 3) = zero
          aframe(3, 3) = one
C
        else
C
          aframe(1, 1) = cmpfr(1, cfrm)
          aframe(2, 1) = cmpfr(2, cfrm)
          aframe(3, 1) = cmpfr(3, cfrm)
          aframe(1, 2) = cmpfr(4, cfrm)
          aframe(2, 2) = cmpfr(5, cfrm)
          aframe(3, 2) = cmpfr(6, cfrm)
          aframe(1, 3) = cmpfr(7, cfrm)
          aframe(2, 3) = cmpfr(8, cfrm)
          aframe(3, 3) = cmpfr(9, cfrm)
C
        endif
C
C.....GET THE ELEMENT TRIANGULAR COORDINATES
C.....GET THE ROTATION MATRIX
C.....GET THE DEGREE OF FREEDOM POINTERS
C
        call g_compcrd(g_p_, elm, ctyp, x, g_x, ldg_x, y, g_y, ldg_y, z,
     * g_z, ldg_z, rot, g_rot, g_pmax_, xlp, g_xlp, g_pmax_, ylp, g_ylp,
     * g_pmax_, zlp, rowb, colb, rowm, colm)
C
C.....GET THE ELEMENT LEVEL FRAME FROM THE 6 x 6 ROTATION MATRIX
C
        do g_i_ = 1, g_p_
          g_eframe(g_i_, 1, 1) = g_rot(g_i_, 1, 1)
        enddo
        eframe(1, 1) = rot(1, 1)
C--------
        do g_i_ = 1, g_p_
          g_eframe(g_i_, 2, 1) = g_rot(g_i_, 2, 1)
        enddo
        eframe(2, 1) = rot(2, 1)
C--------
        do g_i_ = 1, g_p_
          g_eframe(g_i_, 3, 1) = g_rot(g_i_, 3, 1)
        enddo
        eframe(3, 1) = rot(3, 1)
C--------
        do g_i_ = 1, g_p_
          g_eframe(g_i_, 1, 2) = g_rot(g_i_, 1, 2)
        enddo
        eframe(1, 2) = rot(1, 2)
C--------
        do g_i_ = 1, g_p_
          g_eframe(g_i_, 2, 2) = g_rot(g_i_, 2, 2)
        enddo
        eframe(2, 2) = rot(2, 2)
C--------
        do g_i_ = 1, g_p_
          g_eframe(g_i_, 3, 2) = g_rot(g_i_, 3, 2)
        enddo
        eframe(3, 2) = rot(3, 2)
C--------
        do g_i_ = 1, g_p_
          g_eframe(g_i_, 1, 3) = g_rot(g_i_, 1, 3)
        enddo
        eframe(1, 3) = rot(1, 3)
C--------
        do g_i_ = 1, g_p_
          g_eframe(g_i_, 2, 3) = g_rot(g_i_, 2, 3)
        enddo
        eframe(2, 3) = rot(2, 3)
C--------
        do g_i_ = 1, g_p_
          g_eframe(g_i_, 3, 3) = g_rot(g_i_, 3, 3)
        enddo
        eframe(3, 3) = rot(3, 3)
C--------
C
C.....GET THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
        call g_compcst(g_p_, e, g_e, ldg_e, thick, g_thick, g_pmax_, nu,
     * cstcoef, g_cstcoef, g_pmax_, nlayer, idcmp23, mtcmp23, g_mtcmp23,
     * g_pmax_, xlp, g_xlp, g_pmax_, ylp, g_ylp, g_pmax_, zlp, cstbb, g_
     *cstbb, g_pmax_, ctyp, eframe, g_eframe, g_pmax_, aframe, 'BB')
C
C.....GET THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
        call g_compcst(g_p_, e, g_e, ldg_e, thick, g_thick, g_pmax_, nu,
     * cstcoef, g_cstcoef, g_pmax_, nlayer, idcmp23, mtcmp23, g_mtcmp23,
     * g_pmax_, xlp, g_xlp, g_pmax_, ylp, g_ylp, g_pmax_, zlp, cstmm, g_
     *cstmm, g_pmax_, ctyp, eframe, g_eframe, g_pmax_, aframe, 'MM')
C
C.....GET THE CONSTITUTIVE MATRIX FOR BENDING-MEMBRANE COUPLING
C
        call g_compcst(g_p_, e, g_e, ldg_e, thick, g_thick, g_pmax_, nu,
     * cstcoef, g_cstcoef, g_pmax_, nlayer, idcmp23, mtcmp23, g_mtcmp23,
     * g_pmax_, xlp, g_xlp, g_pmax_, ylp, g_ylp, g_pmax_, zlp, cstbm, g_
     *cstbm, g_pmax_, ctyp, eframe, g_eframe, g_pmax_, aframe, 'BM')
C
C.....GET THE CONSTITUTIVE MATRIX FOR MEMBRANE-BENDING COUPLING
C
        if (fastcal) then
C
          do g_i_ = 1, g_p_
            g_cstmb(g_i_, 1, 1) = g_cstbm(g_i_, 1, 1)
          enddo
          cstmb(1, 1) = cstbm(1, 1)
C--------
          do g_i_ = 1, g_p_
            g_cstmb(g_i_, 2, 1) = g_cstbm(g_i_, 1, 2)
          enddo
          cstmb(2, 1) = cstbm(1, 2)
C--------
          do g_i_ = 1, g_p_
            g_cstmb(g_i_, 3, 1) = g_cstbm(g_i_, 1, 3)
          enddo
          cstmb(3, 1) = cstbm(1, 3)
C--------
          do g_i_ = 1, g_p_
            g_cstmb(g_i_, 1, 2) = g_cstbm(g_i_, 2, 1)
          enddo
          cstmb(1, 2) = cstbm(2, 1)
C--------
          do g_i_ = 1, g_p_
            g_cstmb(g_i_, 2, 2) = g_cstbm(g_i_, 2, 2)
          enddo
          cstmb(2, 2) = cstbm(2, 2)
C--------
          do g_i_ = 1, g_p_
            g_cstmb(g_i_, 3, 2) = g_cstbm(g_i_, 2, 3)
          enddo
          cstmb(3, 2) = cstbm(2, 3)
C--------
          do g_i_ = 1, g_p_
            g_cstmb(g_i_, 1, 3) = g_cstbm(g_i_, 3, 1)
          enddo
          cstmb(1, 3) = cstbm(3, 1)
C--------
          do g_i_ = 1, g_p_
            g_cstmb(g_i_, 2, 3) = g_cstbm(g_i_, 3, 2)
          enddo
          cstmb(2, 3) = cstbm(3, 2)
C--------
          do g_i_ = 1, g_p_
            g_cstmb(g_i_, 3, 3) = g_cstbm(g_i_, 3, 3)
          enddo
          cstmb(3, 3) = cstbm(3, 3)
C--------
C
        else
C
          call g_compcst(g_p_, e, g_e, ldg_e, thick, g_thick, g_pmax_, n
     *u, cstcoef, g_cstcoef, g_pmax_, nlayer, idcmp23, mtcmp23, g_mtcmp2
     *3, g_pmax_, xlp, g_xlp, g_pmax_, ylp, g_ylp, g_pmax_, zlp, cstmb, 
     *g_cstmb, g_pmax_, ctyp, eframe, g_eframe, g_pmax_, aframe, 'MB')
C
        endif
C
C.....CHECK THE CONSTITUTIVE MATRIX
C.....(COMMENTED OUT HERE - USED FOR DEBUGGING ONLY)
C
C     call compchk( elm , cstbb , cstmm , cstbm , cstmb )
C
C.....FORM THE LOCAL BASIC STIFFNESS FOR PURE BENDING
C
        call g_compbbb(g_p_, elm, ctyp, xlp, g_xlp, g_pmax_, ylp, g_ylp,
     * g_pmax_, cstbb, g_cstbb, g_pmax_, one, zero, one, rowb, colb, rot
     *, g_rot, g_pmax_, lb, g_lb, g_pmax_, kbbb, g_kbbb, g_pmax_)
C
C.....FORM THE LOCAL BASIC STIFFNESS FOR PURE MEMBRANE
C
        call g_compbmm(g_p_, elm, ctyp, xlp, g_xlp, g_pmax_, ylp, g_ylp,
     * g_pmax_, cstmm, g_cstmm, g_pmax_, onehalf, one, rowm, colm, rot, 
     *g_rot, g_pmax_, pb, g_pb, g_pmax_, kbmm, g_kbmm, g_pmax_)
C
C.....FORM THE LOCAL BASIC STIFFNESS FOR BENDING-MEMBRANE COUPLING
C
        call g_compbbm(g_p_, elm, ctyp, xlp, g_xlp, g_pmax_, ylp, g_ylp,
     * g_pmax_, cstbm, g_cstbm, g_pmax_, one, zero, one, onehalf, one, r
     *owb, colm, rot, g_rot, g_pmax_, lb, g_lb, g_pmax_, pb, g_pb, g_pma
     *x_, fastcal, kbbm, g_kbbm, g_pmax_)
C
C.....FORM THE LOCAL BASIC STIFFNESS FOR MEMBRANE-BENDING COUPLING
C
        call g_compbmb(g_p_, elm, ctyp, xlp, g_xlp, g_pmax_, ylp, g_ylp,
     * g_pmax_, cstmb, g_cstmb, g_pmax_, one, zero, one, onehalf, one, r
     *owm, colb, rot, g_rot, g_pmax_, lb, g_lb, g_pmax_, pb, g_pb, g_pma
     *x_, fastcal, kbmb, g_kbmb, g_pmax_)
C
C.....FORM THE LOCAL HIGHER ORDER STIFFNESS FOR PURE BENDING
C
        call g_comphbb(g_p_, elm, ctyp, xlp, g_xlp, g_pmax_, ylp, g_ylp,
     * g_pmax_, cstbb, g_cstbb, g_pmax_, one, rowb, colb, rot, g_rot, g_
     *pmax_, lh1, g_lh1, g_pmax_, lh2, g_lh2, g_pmax_, lh3, g_lh3, g_pma
     *x_, khbb, g_khbb, g_pmax_)
C
C.....FORM THE LOCAL HIGHER ORDER STIFFNESS FOR PURE MEMBRANE
C
        call g_comphmm(g_p_, elm, ctyp, xlp, g_xlp, g_pmax_, ylp, g_ylp,
     * g_pmax_, cstmm, g_cstmm, g_pmax_, point32, rowm, colm, rot, g_rot
     *, g_pmax_, ph1, g_ph1, g_pmax_, ph2, g_ph2, g_pmax_, ph3, g_ph3, g
     *_pmax_, khmm, g_khmm, g_pmax_)
C
C.....FORM THE LOCAL HIGHER ORDER STIFFNESS FOR BENDING-MEMBRANE COUPLING
C
C
        call g_comphbm(g_p_, elm, ctyp, xlp, g_xlp, g_pmax_, ylp, g_ylp,
     * g_pmax_, cstbm, g_cstbm, g_pmax_, one, point32, rowb, colm, rot, 
     *g_rot, g_pmax_, lh1, g_lh1, g_pmax_, lh2, g_lh2, g_pmax_, lh3, g_l
     *h3, g_pmax_, ph1, g_ph1, g_pmax_, ph2, g_ph2, g_pmax_, ph3, g_ph3,
     * g_pmax_, fastcal, khbm, g_khbm, g_pmax_)
C
C.....FORM THE LOCAL HIGHER ORDER STIFFNESS FOR MEMBRANE-BENDING COUPLING
C
C
        call g_comphmb(g_p_, elm, ctyp, xlp, g_xlp, g_pmax_, ylp, g_ylp,
     * g_pmax_, cstmb, g_cstmb, g_pmax_, one, point32, rowm, colb, rot, 
     *g_rot, g_pmax_, lh1, g_lh1, g_pmax_, lh2, g_lh2, g_pmax_, lh3, g_l
     *h3, g_pmax_, ph1, g_ph1, g_pmax_, ph2, g_ph2, g_pmax_, ph3, g_ph3,
     * g_pmax_, fastcal, khmb, g_khmb, g_pmax_)
C
C.....ADD ALL LOCAL STIFFNESS MATRICES
C
        do 99966 j = 1, ndof
          do 99967 i = 1, ndof
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = g_khbb(g_i_, i, j) + g_kbmb(g_i_, i, j) + g
     *_kbbm(g_i_, i, j) + g_kbmm(g_i_, i, j) + g_kbbb(g_i_, i, j)
            enddo
            d1_w = kbbb(i, j) + kbmm(i, j) + kbbm(i, j) + kbmb(i, j) + k
     *hbb(i, j)
            do g_i_ = 1, g_p_
              g_estiff(g_i_, i, j) = g_khmb(g_i_, i, j) + g_khbm(g_i_, i
     *, j) + g_khmm(g_i_, i, j) + g_d1_w(g_i_)
            enddo
            estiff(i, j) = d1_w + khmm(i, j) + khbm(i, j) + khmb(i, j)
C--------
3002        continue
99967     continue
3001      continue
99966   continue
C
C.....ROTATE ELEMENT STIFFNESS TO GLOBAL COORDINATES
C
C      call compmrot(estiff, rot, rot, rot)
C     
C     if flag equals 1, transform from local to global frame
C     else return local element stiffness matrix.
C
        if (flag .eq. 1) then
          call g_trirotation(g_p_, estiff, g_estiff, ldg_estiff, eframe,
     * g_eframe, g_pmax_)
        endif
C
C.....CHECK THE POSITIVITY OF THE OUTPUT STIFFNESS MATRIX
C.....(COMMENTED OUT HERE - USED FOR DEBUGGING ONLY)
C
C     call compchk2( elm , estiff )
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
C.....ERROR-MESSAGE IF DIMENSIONS DO NOT AGREE
C
100     continue
        write (*, *) '*** FATAL ERROR in Routine COMPST     ***'
        write (*, *) '*** The Elemental Stiffness Matrix    ***'
        write (*, *) '*** Must Be a Square 18 by 18 Matrix: ***'
        write (*, *) '*** Dimensions Do Not Check Out...    ***'
        write (*, *) '*** Execution Terminated Here         ***'
        stop
C
C.....ERROR-MESSAGE IF UNKNOWN TYPE OF CONSTITUTIVE LAW
C
200     continue
        write (*, *) '*** FATAL ERROR in Routine COMPST ***'
        write (*, *) '*** The Type of Constitutive Law  ***'
        write (*, *) '*** Must Either be 0, 1, 2, or 3. ***'
        write (*, *) '*** Execution Terminated Here     ***'
        stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [CMPCO] IS NOT CORRECT
C
300     continue
        write (*, *) '*** FATAL ERROR in Routine COMPST ***'
        write (*, *) '*** The Address in Array [CMPCO]  ***'
        write (*, *) '*** is Out-of-Bounds.             ***'
        write (*, *) '*** Execution Terminated Here     ***'
        stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [IDLAY]/[MTLAY] IS NOT CORRECT
C
400     continue
        write (*, *) '*** FATAL ERROR in Routine COMPST ***'
        write (*, *) '*** The Address in Arrays [IDLAY] ***'
        write (*, *) '*** and [MTLAY] is Out-of-Bounds. ***'
        write (*, *) '*** Execution Terminated Here     ***'
        stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [CMPFR] IS NOT CORRECT
C
500     continue
        write (*, *) '*** FATAL ERROR in Routine COMPST ***'
        write (*, *) '*** The Address in Array [CMPFR]  ***'
        write (*, *) '*** is Out-of-Bounds.             ***'
        write (*, *) '*** Execution Terminated Here     ***'
        stop
C
C.....ERROR-MESSAGE IF THE MAXIMUM NUMBER OF LAYERS HAS BEEN EXCEEDED
C
600     continue
        write (*, *) '*** FATAL ERROR in Routine COMPST         ***'
        write (*, *) '*** The Maximum Number of Layers Has Been ***'
        write (*, *) '*** Exceeded: Boost Parameter [MAXLAYER]. ***'
        write (*, *) '*** Execution Terminated Here             ***'
        stop
C
C.....ERROR-MESSAGE IF THE TOTAL NUMBER OF LAYERS IS NOT CONSISTENT
C
700     continue
        write (*, *) '*** FATAL ERROR in Routine COMPST       ***'
        write (*, *) '*** A Layer Number Exceeds the Total    ***'
        write (*, *) '*** Number of Layers Stored in [IDLAY]. ***'
        write (*, *) '*** Execution Terminated Here           ***'
        stop
C
C.....ERROR-MESSAGE IF BASIC PARAMETERS FOR TYPES 2 & 3 ARE INCONSISTENT
C
800     continue
        write (*, *) '*** FATAL ERROR in Routine COMPST         ***'
        write (*, *) '*** The Element Numbers and Total Number  ***'
        write (*, *) '*** of Composite Layers Do Not Check Out. ***'
        write (*, *) '*** Execution Terminated Here             ***'
        stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPST"
C=======================C
