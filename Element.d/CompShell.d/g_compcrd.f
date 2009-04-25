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
      subroutine g_compcrd(g_p_, elm, type, x, g_x, ldg_x, y, g_y, ldg_y
     *, z, g_z, ldg_z, rot, g_rot, ldg_rot, xlp, g_xlp, ldg_xlp, ylp, g_
     *ylp, ldg_ylp, zlp, rowb, colb, rowm, colm)
C=====================================================================C
C                                                                     C
C     Perform =    This subroutine computes basic quantities needed   C
C     ---------    for the assembly of the basic and higher order     C
C                  composite stiffness matrices.                      C
C                                                                     C
C                                                                     C
C     Inputs/Outputs =                                                C
C     ----------------                                                C
C     ELM   <input>   finite element number                           C
C     TYPE  <input>   type of constitutive law                        C
C     X     <input>   nodal coordinates in the X-direction            C
C     Y     <input>   nodal coordinates in the Y-direction            C
C     Z     <input>   nodal coordinates in the Z-direction            C
C     ROT   <output>  rotation matrix obtained from nodal points      C
C     XLP   <output>  triangular coordinates in the X-direction       C
C     YLP   <output>  triangular coordinates in the Y-direction       C
C     ZLP   <output>  triangular coordinates in the Z-direction       C
C     ROWB  <output>  row pointer for bending DOFs                    C
C     COLB  <output>  column pointer for bending DOFs                 C
C     ROWM  <output>  row pointer for membrane DOFs                   C
C     COLM  <output>  column pointer for membrane DOFs                C
C                                                                     C
C                                                                     C
C     Computations =                                                  C
C     --------------                                                  C
C                                                                     C
C                                                                     C
C     Caution =                                                       C
C     ---------                                                       C
C                                                                     C
C                                                                     C
C     Outputs = no output.                                            C
C     ---------                                                       C
C                                                                     C
C=====================================================================C
C=Author  = Francois M. Hemez                                         C
C=Date    = June 9th, 1994                                            C
C=Version = 1.0                                                       C
C=Comment =                                                           C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
        integer elm, type, rowb(9), colb(9), rowm(9), colm(9)
        real*8 x(3), y(3), z(3)
        real*8 rot(6, 6), xlp(3), ylp(3), zlp(3)
C
C.....Local Variables
C
        integer i, j
        real*8 zero, side21length, side32length, projection
        real*8 x21, y21, z21, x13, y13, z13, x32, y32, z32
        real*8 signedarea, area, lengthy, lengthz, one
        real*8 xp(3), yp(3), zp(3), xcg, ycg, zcg
        real*8 xlcg, ylcg, zlcg, v1n(3), v2n(3), v3n(3)
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_xlp, ldg_ylp, ldg_x, ldg_y, ldg_z, ldg_r
     *ot
        double precision d2_w, d2_b, d3_b, d1_w, d1_p, d5_b, d4_b, d2_v,
     * d3_v, g_xlp(ldg_xlp, 3)
        double precision g_ylp(ldg_ylp, 3), g_x21(g_pmax_), g_x(ldg_x, 3
     *), g_y21(g_pmax_), g_y(ldg_y, 3), g_z21(g_pmax_), g_z(ldg_z, 3), g
     *_x32(g_pmax_), g_y32(g_pmax_), g_z32(g_pmax_)
        double precision g_d2_w(g_pmax_), g_d1_w(g_pmax_), g_side21lengt
     *h(g_pmax_), g_xp(g_pmax_, 3), g_zp(g_pmax_, 3), g_lengthz(g_pmax_)
     *, g_yp(g_pmax_, 3), g_lengthy(g_pmax_), g_xcg(g_pmax_), g_ycg(g_pm
     *ax_)
        double precision g_zcg(g_pmax_), g_xlcg(g_pmax_), g_ylcg(g_pmax_
     *), g_zlcg(g_pmax_), g_rot(ldg_rot, 6, 6)
        save g_zp, g_lengthz, g_yp, g_lengthy, g_xcg, g_ycg, g_zcg, g_xl
     *cg, g_ylcg, g_zlcg
        save g_x21, g_y21, g_z21, g_x32, g_y32, g_z32, g_d2_w, g_d1_w, g
     *_side21length, g_xp
        external g_compfrot
        data zero /0.000000d+00/
        data one /1.000000d+00/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE OUTPUT COORDINATES
C
        integer g_ehfid
        data g_ehfid /0/
C

C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do 99999 i = 1, 3
          do g_i_ = 1, g_p_
            g_xlp(g_i_, i) = 0.0d0
          enddo
          xlp(i) = zero
C--------
1001      continue
99999   continue
C
        do 99998 i = 1, 3
          do g_i_ = 1, g_p_
            g_ylp(g_i_, i) = 0.0d0
          enddo
          ylp(i) = zero
C--------
1002      continue
99998   continue
C
        do 99997 i = 1, 3
          zlp(i) = zero
1003      continue
99997   continue
C
C.....CLEAR THE DEGREE OF FREEDOM POINTERS
C
        do 99996 i = 1, 9
          rowb(i) = 0
1006      continue
99996   continue
C
        do 99995 i = 1, 9
          colb(i) = 0
1007      continue
99995   continue
C
        do 99994 i = 1, 9
          rowm(i) = 0
1008      continue
99994   continue
C
        do 99993 i = 1, 9
          colm(i) = 0
1009      continue
99993   continue
C
C.....COMPUTE THE NODAL COORDINATE DIFFERENCES
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
        x13 = x(1) - x(3)
        y13 = y(1) - y(3)
        z13 = z(1) - z(3)
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
C.....COMPUTE THE LENGTH OF SIDE 2-1
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
     +'g_compcrd.f',
     +228)
        endif
        do g_i_ = 1, g_p_
          g_side21length(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        side21length = d2_v
C--------
C
C.....CHECK IF LENGTH 2-1 IS DIFFERENT FROM ZERO
C
        if (side21length .eq. zero) then
          goto 100
        endif
C
C.....COMPUTE THE LENGTH OF SIDE 3-2
C
        side32length = sqrt(x32 * x32 + y32 * y32 + z32 * z32)
C
C.....CHECK IF LENGTH 3-2 IS DIFFERENT FROM ZERO
C
        if (side32length .eq. zero) then
          goto 200
        endif
C
C.....COMPUTE THE DISTANCE OF THE OPPOSING NODE 3 TO SIDE 2-1
C
        projection = abs(x21 * x32 + y21 * y32 + z21 * z32) / side21leng
     *th
C
C.....GET THE AREA OF THE TRIANGLE
C
        signedarea = side32length * side32length - projection * projecti
     *on
C
        if (signedarea .le. zero) then
          goto 300
        endif
C
        area = 0.50d+00 * side21length * sqrt(signedarea)
C
C.....COMPUTE THE DIRECTION COSINES OF THE LOCAL SYSTEM
C.....DIRECTION [X] IS DIRECTED PARALLEL TO THE SIDE 2-1
C.....DIRECTION [Z] IS THE EXTERNAL NORMAL (COUNTERCLOCKWISE)
C.....DIRECTION [Y] IS COMPUTED AS [Z] x [X] (TENSORIAL PRODUCT)
C
        d3_v = x21 / side21length
        d2_b = 1.0d0 / side21length
        d3_b = -d3_v / side21length
        do g_i_ = 1, g_p_
          g_xp(g_i_, 1) = d3_b * g_side21length(g_i_) + d2_b * g_x21(g_i
     *_)
        enddo
        xp(1) = d3_v
C--------
        d3_v = y21 / side21length
        d2_b = 1.0d0 / side21length
        d3_b = -d3_v / side21length
        do g_i_ = 1, g_p_
          g_xp(g_i_, 2) = d3_b * g_side21length(g_i_) + d2_b * g_y21(g_i
     *_)
        enddo
        xp(2) = d3_v
C--------
        d3_v = z21 / side21length
        d2_b = 1.0d0 / side21length
        d3_b = -d3_v / side21length
        do g_i_ = 1, g_p_
          g_xp(g_i_, 3) = d3_b * g_side21length(g_i_) + d2_b * g_z21(g_i
     *_)
        enddo
        xp(3) = d3_v
C--------
C
        do g_i_ = 1, g_p_
          g_zp(g_i_, 1) = -z21 * g_y32(g_i_) + (-y32) * g_z21(g_i_) + y2
     *1 * g_z32(g_i_) + z32 * g_y21(g_i_)
        enddo
        zp(1) = y21 * z32 - z21 * y32
C--------
        do g_i_ = 1, g_p_
          g_zp(g_i_, 2) = -x21 * g_z32(g_i_) + (-z32) * g_x21(g_i_) + z2
     *1 * g_x32(g_i_) + x32 * g_z21(g_i_)
        enddo
        zp(2) = z21 * x32 - x21 * z32
C--------
        do g_i_ = 1, g_p_
          g_zp(g_i_, 3) = -y21 * g_x32(g_i_) + (-x32) * g_y21(g_i_) + x2
     *1 * g_y32(g_i_) + y32 * g_x21(g_i_)
        enddo
        zp(3) = x21 * y32 - y21 * x32
C--------
C
        d4_b = zp(2) + zp(2)
        d5_b = zp(1) + zp(1)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_zp(g_i_, 2) + d5_b * g_zp(g_i_, 1)
        enddo
        d2_w = zp(1) * zp(1) + zp(2) * zp(2)
        d4_b = zp(3) + zp(3)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_zp(g_i_, 3) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + zp(3) * zp(3)
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compcrd.f',
     +337)
        endif
        do g_i_ = 1, g_p_
          g_lengthz(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        lengthz = d2_v
C--------
C
        if (lengthz .eq. zero) then
          goto 400
        endif
C
        d3_v = zp(1) / lengthz
        d2_b = 1.0d0 / lengthz
        d3_b = -d3_v / lengthz
        do g_i_ = 1, g_p_
          g_zp(g_i_, 1) = d3_b * g_lengthz(g_i_) + d2_b * g_zp(g_i_, 1)
        enddo
        zp(1) = d3_v
C--------
        d3_v = zp(2) / lengthz
        d2_b = 1.0d0 / lengthz
        d3_b = -d3_v / lengthz
        do g_i_ = 1, g_p_
          g_zp(g_i_, 2) = d3_b * g_lengthz(g_i_) + d2_b * g_zp(g_i_, 2)
        enddo
        zp(2) = d3_v
C--------
        d3_v = zp(3) / lengthz
        d2_b = 1.0d0 / lengthz
        d3_b = -d3_v / lengthz
        do g_i_ = 1, g_p_
          g_zp(g_i_, 3) = d3_b * g_lengthz(g_i_) + d2_b * g_zp(g_i_, 3)
        enddo
        zp(3) = d3_v
C--------
C
        do g_i_ = 1, g_p_
          g_yp(g_i_, 1) = -zp(3) * g_xp(g_i_, 2) + (-xp(2)) * g_zp(g_i_,
     * 3) + zp(2) * g_xp(g_i_, 3) + xp(3) * g_zp(g_i_, 2)
        enddo
        yp(1) = zp(2) * xp(3) - zp(3) * xp(2)
C--------
        do g_i_ = 1, g_p_
          g_yp(g_i_, 2) = -zp(1) * g_xp(g_i_, 3) + (-xp(3)) * g_zp(g_i_,
     * 1) + zp(3) * g_xp(g_i_, 1) + xp(1) * g_zp(g_i_, 3)
        enddo
        yp(2) = zp(3) * xp(1) - zp(1) * xp(3)
C--------
        do g_i_ = 1, g_p_
          g_yp(g_i_, 3) = -zp(2) * g_xp(g_i_, 1) + (-xp(1)) * g_zp(g_i_,
     * 2) + zp(1) * g_xp(g_i_, 2) + xp(2) * g_zp(g_i_, 1)
        enddo
        yp(3) = zp(1) * xp(2) - zp(2) * xp(1)
C--------
C
        d4_b = yp(2) + yp(2)
        d5_b = yp(1) + yp(1)
        do g_i_ = 1, g_p_
          g_d2_w(g_i_) = d4_b * g_yp(g_i_, 2) + d5_b * g_yp(g_i_, 1)
        enddo
        d2_w = yp(1) * yp(1) + yp(2) * yp(2)
        d4_b = yp(3) + yp(3)
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_yp(g_i_, 3) + g_d2_w(g_i_)
        enddo
        d1_w = d2_w + yp(3) * yp(3)
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compcrd.f',
     +410)
        endif
        do g_i_ = 1, g_p_
          g_lengthy(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        lengthy = d2_v
C--------
C
        if (lengthy .eq. zero) then
          goto 400
        endif
C
        d3_v = yp(1) / lengthy
        d2_b = 1.0d0 / lengthy
        d3_b = -d3_v / lengthy
        do g_i_ = 1, g_p_
          g_yp(g_i_, 1) = d3_b * g_lengthy(g_i_) + d2_b * g_yp(g_i_, 1)
        enddo
        yp(1) = d3_v
C--------
        d3_v = yp(2) / lengthy
        d2_b = 1.0d0 / lengthy
        d3_b = -d3_v / lengthy
        do g_i_ = 1, g_p_
          g_yp(g_i_, 2) = d3_b * g_lengthy(g_i_) + d2_b * g_yp(g_i_, 2)
        enddo
        yp(2) = d3_v
C--------
        d3_v = yp(3) / lengthy
        d2_b = 1.0d0 / lengthy
        d3_b = -d3_v / lengthy
        do g_i_ = 1, g_p_
          g_yp(g_i_, 3) = d3_b * g_lengthy(g_i_) + d2_b * g_yp(g_i_, 3)
        enddo
        yp(3) = d3_v
C--------
C
C.....COMPUTE THE COORDINATES FOR THE CENTER OF GRAVITY
C
        d2_b = 1.0d0 / 3.00d+00
        do g_i_ = 1, g_p_
          g_xcg(g_i_) = d2_b * g_x(g_i_, 3) + d2_b * g_x(g_i_, 2) + d2_b
     * * g_x(g_i_, 1)
        enddo
        xcg = (x(1) + x(2) + x(3)) / 3.00d+00
C--------
        d2_b = 1.0d0 / 3.00d+00
        do g_i_ = 1, g_p_
          g_ycg(g_i_) = d2_b * g_y(g_i_, 3) + d2_b * g_y(g_i_, 2) + d2_b
     * * g_y(g_i_, 1)
        enddo
        ycg = (y(1) + y(2) + y(3)) / 3.00d+00
C--------
        d2_b = 1.0d0 / 3.00d+00
        do g_i_ = 1, g_p_
          g_zcg(g_i_) = d2_b * g_z(g_i_, 3) + d2_b * g_z(g_i_, 2) + d2_b
     * * g_z(g_i_, 1)
        enddo
        zcg = (z(1) + z(2) + z(3)) / 3.00d+00
C--------
C
C.....COMPUTE THE LOCAL COORDINATES
C
        do 99992 i = 1, 3
          do g_i_ = 1, g_p_
            g_xlcg(g_i_) = -g_xcg(g_i_) + g_x(g_i_, i)
          enddo
          xlcg = x(i) - xcg
C--------
          do g_i_ = 1, g_p_
            g_ylcg(g_i_) = -g_ycg(g_i_) + g_y(g_i_, i)
          enddo
          ylcg = y(i) - ycg
C--------
          do g_i_ = 1, g_p_
            g_zlcg(g_i_) = -g_zcg(g_i_) + g_z(g_i_, i)
          enddo
          zlcg = z(i) - zcg
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = xp(2) * g_ylcg(g_i_) + ylcg * g_xp(g_i_, 2) +
     * xp(1) * g_xlcg(g_i_) + xlcg * g_xp(g_i_, 1)
          enddo
          d1_w = xp(1) * xlcg + xp(2) * ylcg
          do g_i_ = 1, g_p_
            g_xlp(g_i_, i) = xp(3) * g_zlcg(g_i_) + zlcg * g_xp(g_i_, 3)
     * + g_d1_w(g_i_)
          enddo
          xlp(i) = d1_w + xp(3) * zlcg
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = yp(2) * g_ylcg(g_i_) + ylcg * g_yp(g_i_, 2) +
     * yp(1) * g_xlcg(g_i_) + xlcg * g_yp(g_i_, 1)
          enddo
          d1_w = yp(1) * xlcg + yp(2) * ylcg
          do g_i_ = 1, g_p_
            g_ylp(g_i_, i) = yp(3) * g_zlcg(g_i_) + zlcg * g_yp(g_i_, 3)
     * + g_d1_w(g_i_)
          enddo
          ylp(i) = d1_w + yp(3) * zlcg
C--------
          zlp(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg
2001      continue
99992   continue
C
C.....COMPUTE THE NODAL ROTATION MATRIX
C
        v1n(1) = one
        v1n(2) = zero
        v1n(3) = zero
C
        v2n(1) = zero
        v2n(2) = one
        v2n(3) = zero
C
        v3n(1) = zero
        v3n(2) = zero
        v3n(3) = one
C
        call g_compfrot(g_p_, xp, g_xp, g_pmax_, yp, g_yp, g_pmax_, zp, 
     *g_zp, g_pmax_, v1n, v2n, v3n, rot, g_rot, ldg_rot)
C
C.....DEFINE ROW POINTER FOR BENDING CONTRIBUTIONS
C
        rowb(1) = 3
        rowb(2) = 4
        rowb(3) = 5
        rowb(4) = 9
        rowb(5) = 10
        rowb(6) = 11
        rowb(7) = 15
        rowb(8) = 16
        rowb(9) = 17
C
C.....DEFINE COLUMN POINTER FOR BENDING CONTRIBUTIONS
C
        colb(1) = 3
        colb(2) = 4
        colb(3) = 5
        colb(4) = 9
        colb(5) = 10
        colb(6) = 11
        colb(7) = 15
        colb(8) = 16
        colb(9) = 17
C
C.....DEFINE ROW POINTER FOR MEMBRANE CONTRIBUTION
C
        rowm(1) = 1
        rowm(2) = 2
        rowm(3) = 7
        rowm(4) = 8
        rowm(5) = 13
        rowm(6) = 14
        rowm(7) = 6
        rowm(8) = 12
        rowm(9) = 18
C
C.....DEFINE COLUMN POINTER FOR MEMBRANE CONTRIBUTION
C
        colm(1) = 1
        colm(2) = 2
        colm(3) = 7
        colm(4) = 8
        colm(5) = 13
        colm(6) = 14
        colm(7) = 6
        colm(8) = 12
        colm(9) = 18
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
C.....ERROR-MESSAGE IF SIDE 2-1 HAS ZERO-LENGTH
C
100     continue
        write (*, *) '*** FATAL ERROR in routine COMPCRD          ***'
        write (*, *) '*** Side Between Nodes 1 and 2 Has 0-Length ***'
        write (*, *) '*** Check Coordinates and FE Topology       ***'
        write (*, *) '*** EXECUTION TERMINATED HERE               ***'
        stop
C
C.....ERROR-MESSAGE IF SIDE 3-2 HAS ZERO-LENGTH
C
200     continue
        write (*, *) '*** FATAL ERROR in routine COMPCRD          ***'
        write (*, *) '*** Side Between Nodes 2 and 3 Has 0-Length ***'
        write (*, *) '*** Check Coordinates and FE Topology       ***'
        write (*, *) '*** EXECUTION TERMINATED HERE               ***'
        stop
C
C.....ERROR-MESSAGE IF THE AREA IS NEGATIVE OR ZERO
C
300     continue
        write (*, *) '*** FATAL ERROR in routine COMPCRD    ***'
        write (*, *) '*** The Area is Negative or Zero      ***'
        write (*, *) '*** Check Coordinates and FE Topology ***'
        write (*, *) '*** EXECUTION TERMINATED HERE         ***'
        stop
C
C.....ERROR-MESSAGE IF A LOCAL FRAME VECTOR HAS ZERO-LENGTH
C
400     continue
        write (*, *) '*** FATAL ERROR in routine COMPCRD     ***'
        write (*, *) '*** A Local Frame Vector has 0-Length  ***'
        write (*, *) '*** Check FE Topology and Computations ***'
        write (*, *) '*** EXECUTION TERMINATED HERE          ***'
        stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPCRD"
C========================C
