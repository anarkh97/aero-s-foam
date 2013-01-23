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
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine gxsands7(elm, a, g_a, e, g_e, efram
     *e, g_eframe, ix, g_ix, iy, g_iy, iz, g_iz, alphay, g_alphay, 
     *alphaz, g_alphaz, c1, g_c1, nu, g_nu, x, g_x, y, g_y, 
     *z, g_z, ug, g_ug, g_stress, 
     * numel, maxgus, maxstr, msize)
C=====================================================================C
C                                                                     C
C     This Routine Retreives the Stresses for the 3D Timoshenko Beam  C
C     Element.                                                        C
C                                                                     C
C     Francois M. Hemez - July 12th 1994 - Version 1.0                C
C                                                                     C
C=====================================================================C
C
C     ------------
C     DECLARATIONS
C     ------------
C
C.....GLOBAL VARIABLES
C
        integer elm, numel, maxstr, maxgus, msize
        real*8 a, e, nu, ix, iy, iz, eframe(*)
        real*8 ug(*), x(2), y(2), z(2), alphay, alphaz, c1
C
C.....LOCAL VARIABLES
C
        integer i, j, k, l
        real*8 t(9), zero, half
        real*8 jj, dx, dy, dz, length, g, one, two
        real*8 ue(12), estif(12, 12), fe(12), le(12, 12)
        real*8 c11, c22, c33, c12, c13, c23, eps
        real*8 b, gammay, gammaz, twelve, six, xke
        real*8 four, bendy, bendz, bendcy, bendcz
        real*8 locke(12, 12), s11, s22, s33, s26, s35
        real*8 s44, s55, s66, s17, s28, s39, s59, s68
        real*8 s77, s88, s99, s212, s311, s410, s511
        real*8 s612, s812, s911, s1010, s1111, s1212
        logical ortho, specialcasey, specialcasez
C
C.....LOCAL VARIABLES FOR FRAME COMPUTATION
C
C**   integer   index
C**   real*8    normT , Id(3,3) , xx , yy , zz , xnorm
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_x, ldg_y, ldg_z, ldg_eframe, ldg_ix, ldg
     *_e, ldg_nu, ldg_iy
        
	parameter (g_p_=1)
	parameter (ldg_x=1)
	parameter (ldg_y=1)
	parameter (ldg_z=1)
	parameter (ldg_eframe=1)
	parameter (ldg_ix=1)
	parameter (ldg_e=1)
	parameter (ldg_nu=1)
	parameter (ldg_iy=1)
     
        integer ldg_iz, ldg_a, ldg_c1, ldg_alphay, ldg_alphaz, ldg_ug, l
     *dg_stress
        
	parameter (ldg_iz=1)
	parameter (ldg_a=1)
	parameter (ldg_c1=1)
	parameter (ldg_alphay=1)
	parameter (ldg_alphaz=1)
	parameter (ldg_ug=1)
	parameter (ldg_stress=1)
     
        double precision d11_b, d11_v, d10_b, d9_b, d10_v, d9_v, d8_b, d
     *7_b, d2_v, d3_v
        double precision d6_b, d2_b, d3_b, d1_w, d2_w, d4_v, d5_v, d4_b,
     * d5_b, d1_p
        double precision d6_v, d7_v, d8_v, g_estif(g_pmax_, 12, 12), g_t
     *(g_pmax_, 9), g_locke(g_pmax_, 12, 12), g_le(g_pmax_, 12, 12), g_u
     *e(g_pmax_, 12), g_fe(g_pmax_, 12), g_dx(g_pmax_)
        double precision g_x(ldg_x, 2), g_dy(g_pmax_), g_y(ldg_y, 2), g_
     *dz(g_pmax_), g_z(ldg_z, 2), g_d2_w(g_pmax_), g_d1_w(g_pmax_), g_le
     *ngth(g_pmax_), g_eframe(ldg_eframe, *), g_jj(g_pmax_)
        double precision g_ix(ldg_ix), g_g(g_pmax_), g_e(ldg_e), g_nu(ld
     *g_nu), g_s11(g_pmax_), g_s22(g_pmax_), g_s33(g_pmax_), g_s26(g_pma
     *x_), g_s35(g_pmax_), g_s44(g_pmax_)
        double precision g_s55(g_pmax_), g_s66(g_pmax_), g_s17(g_pmax_),
     * g_s28(g_pmax_), g_s39(g_pmax_), g_s59(g_pmax_), g_s68(g_pmax_), g
     *_s77(g_pmax_), g_s88(g_pmax_), g_s99(g_pmax_)
        double precision g_s212(g_pmax_), g_s311(g_pmax_), g_s410(g_pmax
     *_), g_s511(g_pmax_), g_s612(g_pmax_), g_s812(g_pmax_), g_s911(g_pm
     *ax_), g_s1010(g_pmax_), g_s1111(g_pmax_), g_s1212(g_pmax_)
        double precision g_iy(ldg_iy), g_iz(ldg_iz), g_a(ldg_a), g_c1(ld
     *g_c1), g_b(g_pmax_), g_gammay(g_pmax_), g_bendy(g_pmax_), g_alphay
     *(ldg_alphay), g_bendcy(g_pmax_), g_gammaz(g_pmax_)
        double precision g_bendz(g_pmax_), g_alphaz(ldg_alphaz), g_bendc
     *z(g_pmax_), g_xke(g_pmax_), g_ug(ldg_ug, *), g_stress(ldg_stress, 
     *msize, maxstr, maxgus)
        save g_b, g_gammay, g_bendy, g_bendcy, g_gammaz, g_bendz, g_bend
     *cz, g_xke
        save g_s212, g_s311, g_s410, g_s511, g_s612, g_s812, g_s911, g_s
     *1010, g_s1111, g_s1212
        save g_s55, g_s66, g_s17, g_s28, g_s39, g_s59, g_s68, g_s77, g_s
     *88, g_s99
        save g_d1_w, g_length, g_jj, g_g, g_s11, g_s22, g_s33, g_s26, g_
     *s35, g_s44
        save g_estif, g_t, g_locke, g_le, g_ue, g_fe, g_dx, g_dy, g_dz, 
     *g_d2_w
        data zero /0.000000d+00/
        data half /0.500000d+00/
        data one /1.000000d+00/
        data two /2.000000d+00/
        data four /4.000000d+00/
        data six /6.000000d+00/
        data twelve /1.200000d+01/
C
C.....ZERO-MACHINE FOR FRAME ORTHONORMALITY CHECK
C
        data eps /1.000000d-10/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE LOCAL STIFFNESS MATRIX
C
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'sands7','g_sands7.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do 99998 j = 1, 12
          do 99999 i = 1, 12
            do g_i_ = 1, g_p_
              g_estif(g_i_, i, j) = 0.0d0
            enddo
            estif(i, j) = zero
C--------
1002        continue
99999     continue
1001      continue
99998   continue
C
C.....CLEAR THE LOCAL ARRAYS
C
        do 99997 i = 1, 9
          do g_i_ = 1, g_p_
            g_t(g_i_, i) = 0.0d0
          enddo
          t(i) = zero
C--------
1003      continue
99997   continue
C
        do 99995 j = 1, 12
          do 99996 i = 1, 12
            do g_i_ = 1, g_p_
              g_locke(g_i_, i, j) = 0.0d0
            enddo
            locke(i, j) = zero
C--------
1005        continue
99996     continue
1004      continue
99995   continue
C
        do 99993 j = 1, 12
          do 99994 i = 1, 12
            do g_i_ = 1, g_p_
              g_le(g_i_, i, j) = 0.0d0
            enddo
            le(i, j) = zero
C--------
1007        continue
99994     continue
1006      continue
99993   continue
C
        do 99992 i = 1, 12
          do g_i_ = 1, g_p_
            g_ue(g_i_, i) = 0.0d0
          enddo
          ue(i) = zero
C--------
1008      continue
99992   continue
C
        do 99991 i = 1, 12
          do g_i_ = 1, g_p_
            g_fe(g_i_, i) = 0.0d0
          enddo
          fe(i) = zero
C--------
1009      continue
99991   continue
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
     +237)
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
C.....CHECK IF MATERIAL PROPERTIES ARE POSITIVE OR ZERO
C
        if (e .le. zero) then
          goto 200
        endif
        if (a .lt. zero) then
          goto 200
        endif
        if (nu .lt. zero) then
          goto 200
        endif
        if (ix .lt. zero) then
          goto 200
        endif
        if (iy .lt. zero) then
          goto 200
        endif
        if (iz .lt. zero) then
          goto 200
        endif
        if (alphay .lt. zero) then
          goto 200
        endif
        if (alphaz .lt. zero) then
          goto 200
        endif
        if (c1 .lt. zero) then
          goto 200
        endif
C
C.....EXTRACT THE ROTATION MATRIX FROM [EFRAME]
C.....       [ x_X y_X z_X ]
C..... [T] = [ x_Y y_Y z_Y ]
C.....
C
        do 99990 i = 1, 9
          do g_i_ = 1, g_p_
            g_t(g_i_, i) = g_eframe(g_i_, 9 * (elm - 1) + i)
          enddo
          t(i) = eframe(9 * (elm - 1) + i)
C--------
2001      continue
99990   continue
C
C.....CALCULATE THE ROTATION MATRIX IF THE INPUT ONE IS NILL
C.....ARCHTUNG!!! ASSUMES THAT THE BEAM ELEMENT IS SYMMETRIC
C.....(THAT IS: [Iy] = [Iz]) WITH CIRCULAR CROSS-SECTION SO
C.....THAT THE LOCAL [y] AND [z] AXES MAY BE ORIENTED
C.....ARBITRARILY IN THE PLANE NORMAL TO THE BEAM'S AXIS [x]
C
C**   normT = zero
C
C**   do 2002 i=1,9
C**      normT = normT + T(i)*T(i)
C2002 continue
C
C**   if ( ( normT.eq.zero ).and.( Iy.eq.Iz ) ) then
C
C.....INITIALIZE THE IDENTITY MATRIX
C
C**   do 2003 i=1,3
C**      do 2004 j=1,3
C**         Id(i,j) = zero
C2004    continue
C**      Id(i,i) = one
C2003 continue
C
C.....COMPUTE DIRECTION [x] STORED IN FIRST 3 ENTRIES OF [T]
C
C**   T(1) = dx/length
C**   T(2) = dy/length
C**   T(3) = dz/length
C
C.....FIND WHICH AXIS HAS A MAXIMUM ANGLE WITH DIRECTION OF THE BEAM
C
C**   index = 1
C
C**   if ( abs(T(2)).lt.abs(T(1)).and.abs(T(2)).ne.zero ) then
C**      index = 2
C**   endif
C
C**   if ( abs(T(3)).lt.abs(T(1)) ) then
C**      if ( abs(T(3)).lt.abs(T(2)).and.abs(T(3)).ne.zero ) then
C**         index = 3
C**      endif
C**   endif
C
C.....COMPUTE THE CROSS-PRODUCT [axis^x]
C
C**   xx = Id(2,index)*T(3) - Id(3,index)*T(2)
C**   yy = Id(3,index)*T(1) - Id(1,index)*T(3)
C**   zz = Id(1,index)*T(2) - Id(2,index)*T(1)
C
C.....COMPUTE ITS NORM
C
C**   xnorm = dsqrt( (xx*xx) + (yy*yy) + (zz*zz) )
C
C.....TREATMENT FOR A ZERO-NORM
C
C**   if ( xnorm.eq.zero ) then
C
C.....THE BEAM IS IN THE [+/-X]-DIRECTION AND [y]=[Y] IS GIVEN
C
C**      if ( index.eq.1 ) then
C**         xx    = zero
C**         yy    = one
C**         zz    = zero
C**         xnorm = one
C**      endif
C
C.....THE BEAM IS IN THE [+/-Y]-DIRECTION AND [y]=[Z] IS GIVEN
C
C**      if ( index.eq.2 ) then
C**         xx    = zero
C**         yy    = zero
C**         zz    = one
C**         xnorm = one
C**      endif
C
C.....THE BEAM IS IN THE [+/-Z]-DIRECTION AND [y]=[X] IS GIVEN
C
C**      if ( index.eq.3 ) then
C**         xx    = one
C**         yy    = zero
C**         zz    = zero
C**         xnorm = one
C**      endif
C
C.....END OF EXCEPTION TREATMENT (IF [xnorm] IS ZERO)
C
C**   endif
C
C.....NORMALIZE [y] AND STORE IN ENTRIES 4 TO 6 OF [T]
C
C**   T(4) = xx/xnorm
C**   T(5) = yy/xnorm
C**   T(6) = zz/xnorm
C
C.....COMPUTE [z=x^y] AND STORE IN ENTRIES 7 TO 9 OF [T]
C
C**   T(7) = T(2)*T(6) - T(3)*T(5)
C**   T(8) = T(3)*T(4) - T(1)*T(6)
C**   T(9) = T(1)*T(5) - T(2)*T(4)
C
C.....CHECK ORTHOGONALITY OF ROTATION MATRIX [T]
C
        c11 = t(1) * t(1) + t(2) * t(2) + t(3) * t(3)
        c22 = t(4) * t(4) + t(5) * t(5) + t(6) * t(6)
        c33 = t(7) * t(7) + t(8) * t(8) + t(9) * t(9)
        c12 = t(1) * t(4) + t(2) * t(5) + t(3) * t(6)
        c13 = t(1) * t(7) + t(2) * t(8) + t(3) * t(9)
        c23 = t(4) * t(7) + t(5) * t(8) + t(6) * t(9)
C
C**   write(*,*) elm,(T(i),i=1,9)
C**   write(*,*) elm,c11,c22,c33,c12,c13,c23
C
        ortho = (abs(one - c11) .le. eps .and. abs(one - c22) .le. eps .
     *and. abs(one - c33) .le. eps .and. abs(c12) .le. eps .and. abs(c13
     *) .le. eps .and. abs(c23) .le. eps)
C
C**   if ( .not.ortho ) then
C**      write(*,*) elm,c11,c22,c33,c12,c13,c23
C**      go to 300
C**   endif
C
C.....END OF TREATMENT FOR ROTATION MATRIX OF A SYMMETRIC BEAM
C
C**   endif
C
C.....ERROR-MESSAGE IF THE FRAME IS NOT AVAILABLE FOR A
C.....NON-SYMMETRIC TIMOSHENKO BEAM ELEMENT
C
C**   if ( ( normT.eq.zero ).and.( Iy.ne.Iz ) ) go to 400
C
C.....ASSEMBLE THE TRANFORMATION MATRIX
C.....
C.....        [ [T]  0   0   0  ]
C.....        [  0  [T]  0   0  ]
C..... [Le] = [  0   0  [T]  0  ]
C.....        [  0   0   0  [T] ]
C
        do 99987 k = 1, 4
          do 99988 j = 1, 3
            do 99989 i = 1, 3
              do g_i_ = 1, g_p_
                g_le(g_i_, 3 * (k - 1) + i, 3 * (k - 1) + j) = g_t(g_i_,
     * 3 * (j - 1) + i)
              enddo
              le(3 * (k - 1) + i, 3 * (k - 1) + j) = t(3 * (j - 1) + i)
C--------
2007          continue
99989       continue
2006        continue
99988     continue
2005      continue
99987   continue
C
C.....INITIALIZE THE VARIABLE [JJ] FOR AN ARBITRARY CROSS SECTION
C
        do g_i_ = 1, g_p_
          g_jj(g_i_) = g_ix(g_i_)
        enddo
        jj = ix
C--------
C
C.....INITIALIZE THE VARIABLE [JJ] FOR A CIRCULAR CROSS SECTION
C
C**   JJ = half*A*A/acos(-one)
C
C.....INITIALIZE THE TRANSVERSE SHEAR MODULUS FOR ISOTROPIC MATERIAL
C
        d4_v = two * (one + nu)
        d5_v = e / d4_v
        d2_b = 1.0d0 / d4_v
        d5_b = -d5_v / d4_v * two
        do g_i_ = 1, g_p_
          g_g(g_i_) = d5_b * g_nu(g_i_) + d2_b * g_e(g_i_)
        enddo
        g = d5_v
C--------
C
C.....INITIALIZE THE ENTRIES OF THE LOCAL STIFFNESS MATRIX
C
        do g_i_ = 1, g_p_
          g_s11(g_i_) = 0.0d0
        enddo
        s11 = zero
C--------
        do g_i_ = 1, g_p_
          g_s22(g_i_) = 0.0d0
        enddo
        s22 = zero
C--------
        do g_i_ = 1, g_p_
          g_s33(g_i_) = 0.0d0
        enddo
        s33 = zero
C--------
        do g_i_ = 1, g_p_
          g_s26(g_i_) = 0.0d0
        enddo
        s26 = zero
C--------
        do g_i_ = 1, g_p_
          g_s35(g_i_) = 0.0d0
        enddo
        s35 = zero
C--------
        do g_i_ = 1, g_p_
          g_s44(g_i_) = 0.0d0
        enddo
        s44 = zero
C--------
        do g_i_ = 1, g_p_
          g_s55(g_i_) = 0.0d0
        enddo
        s55 = zero
C--------
        do g_i_ = 1, g_p_
          g_s66(g_i_) = 0.0d0
        enddo
        s66 = zero
C--------
        do g_i_ = 1, g_p_
          g_s17(g_i_) = 0.0d0
        enddo
        s17 = zero
C--------
        do g_i_ = 1, g_p_
          g_s28(g_i_) = 0.0d0
        enddo
        s28 = zero
C--------
        do g_i_ = 1, g_p_
          g_s39(g_i_) = 0.0d0
        enddo
        s39 = zero
C--------
        do g_i_ = 1, g_p_
          g_s59(g_i_) = 0.0d0
        enddo
        s59 = zero
C--------
        do g_i_ = 1, g_p_
          g_s68(g_i_) = 0.0d0
        enddo
        s68 = zero
C--------
        do g_i_ = 1, g_p_
          g_s77(g_i_) = 0.0d0
        enddo
        s77 = zero
C--------
        do g_i_ = 1, g_p_
          g_s88(g_i_) = 0.0d0
        enddo
        s88 = zero
C--------
        do g_i_ = 1, g_p_
          g_s99(g_i_) = 0.0d0
        enddo
        s99 = zero
C--------
        do g_i_ = 1, g_p_
          g_s212(g_i_) = 0.0d0
        enddo
        s212 = zero
C--------
        do g_i_ = 1, g_p_
          g_s311(g_i_) = 0.0d0
        enddo
        s311 = zero
C--------
        do g_i_ = 1, g_p_
          g_s410(g_i_) = 0.0d0
        enddo
        s410 = zero
C--------
        do g_i_ = 1, g_p_
          g_s511(g_i_) = 0.0d0
        enddo
        s511 = zero
C--------
        do g_i_ = 1, g_p_
          g_s612(g_i_) = 0.0d0
        enddo
        s612 = zero
C--------
        do g_i_ = 1, g_p_
          g_s812(g_i_) = 0.0d0
        enddo
        s812 = zero
C--------
        do g_i_ = 1, g_p_
          g_s911(g_i_) = 0.0d0
        enddo
        s911 = zero
C--------
        do g_i_ = 1, g_p_
          g_s1010(g_i_) = 0.0d0
        enddo
        s1010 = zero
C--------
        do g_i_ = 1, g_p_
          g_s1111(g_i_) = 0.0d0
        enddo
        s1111 = zero
C--------
        do g_i_ = 1, g_p_
          g_s1212(g_i_) = 0.0d0
        enddo
        s1212 = zero
C--------
C
C.....INITIALIZE LOGICALS FOR PARTICULAR CASES
C
        specialcasey = .false.
        specialcasez = .false.
C
C.....INITIALIZE IF [ALPHA_Y] IS ZERO
C
        if (alphay .eq. zero) then
          d2_v = twelve * e
          d6_v = length * length
          d7_v = d6_v * length
          d8_v = d2_v * iy / d7_v
          d2_b = 1.0d0 / d7_v
          d3_b = -d8_v / d7_v
          d4_b = d3_b * length
          d5_b = d3_b * d6_v + d4_b * length + d4_b * length
          d7_b = d2_b * d2_v
          d8_b = d2_b * iy * twelve
          do g_i_ = 1, g_p_
            g_s33(g_i_) = d5_b * g_length(g_i_) + d7_b * g_iy(g_i_) + d8
     *_b * g_e(g_i_)
          enddo
          s33 = d8_v
C--------
          d2_v = -six * e
          d6_v = length * length
          d7_v = d2_v * iy / d6_v
          d2_b = 1.0d0 / d6_v
          d3_b = -d7_v / d6_v
          d4_b = d3_b * length + d3_b * length
          d6_b = d2_b * d2_v
          d7_b = d2_b * iy * (-six)
          do g_i_ = 1, g_p_
            g_s35(g_i_) = d4_b * g_length(g_i_) + d6_b * g_iy(g_i_) + d7
     *_b * g_e(g_i_)
          enddo
          s35 = d7_v
C--------
          d2_v = four * e
          d6_v = d2_v * iy / length
          d2_b = 1.0d0 / length
          d3_b = -d6_v / length
          d5_b = d2_b * d2_v
          d6_b = d2_b * iy * four
          do g_i_ = 1, g_p_
            g_s55(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iy(g_i_) + d6
     *_b * g_e(g_i_)
          enddo
          s55 = d6_v
C--------
          d2_v = -twelve * e
          d6_v = length * length
          d7_v = d6_v * length
          d8_v = d2_v * iy / d7_v
          d2_b = 1.0d0 / d7_v
          d3_b = -d8_v / d7_v
          d4_b = d3_b * length
          d5_b = d3_b * d6_v + d4_b * length + d4_b * length
          d7_b = d2_b * d2_v
          d8_b = d2_b * iy * (-twelve)
          do g_i_ = 1, g_p_
            g_s39(g_i_) = d5_b * g_length(g_i_) + d7_b * g_iy(g_i_) + d8
     *_b * g_e(g_i_)
          enddo
          s39 = d8_v
C--------
          d2_v = six * e
          d6_v = length * length
          d7_v = d2_v * iy / d6_v
          d2_b = 1.0d0 / d6_v
          d3_b = -d7_v / d6_v
          d4_b = d3_b * length + d3_b * length
          d6_b = d2_b * d2_v
          d7_b = d2_b * iy * six
          do g_i_ = 1, g_p_
            g_s59(g_i_) = d4_b * g_length(g_i_) + d6_b * g_iy(g_i_) + d7
     *_b * g_e(g_i_)
          enddo
          s59 = d7_v
C--------
          d2_v = twelve * e
          d6_v = length * length
          d7_v = d6_v * length
          d8_v = d2_v * iy / d7_v
          d2_b = 1.0d0 / d7_v
          d3_b = -d8_v / d7_v
          d4_b = d3_b * length
          d5_b = d3_b * d6_v + d4_b * length + d4_b * length
          d7_b = d2_b * d2_v
          d8_b = d2_b * iy * twelve
          do g_i_ = 1, g_p_
            g_s99(g_i_) = d5_b * g_length(g_i_) + d7_b * g_iy(g_i_) + d8
     *_b * g_e(g_i_)
          enddo
          s99 = d8_v
C--------
          d2_v = -six * e
          d6_v = length * length
          d7_v = d2_v * iy / d6_v
          d2_b = 1.0d0 / d6_v
          d3_b = -d7_v / d6_v
          d4_b = d3_b * length + d3_b * length
          d6_b = d2_b * d2_v
          d7_b = d2_b * iy * (-six)
          do g_i_ = 1, g_p_
            g_s311(g_i_) = d4_b * g_length(g_i_) + d6_b * g_iy(g_i_) + d
     *7_b * g_e(g_i_)
          enddo
          s311 = d7_v
C--------
          d2_v = two * e
          d6_v = d2_v * iy / length
          d2_b = 1.0d0 / length
          d3_b = -d6_v / length
          d5_b = d2_b * d2_v
          d6_b = d2_b * iy * two
          do g_i_ = 1, g_p_
            g_s511(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iy(g_i_) + d
     *6_b * g_e(g_i_)
          enddo
          s511 = d6_v
C--------
          d2_v = six * e
          d6_v = length * length
          d7_v = d2_v * iy / d6_v
          d2_b = 1.0d0 / d6_v
          d3_b = -d7_v / d6_v
          d4_b = d3_b * length + d3_b * length
          d6_b = d2_b * d2_v
          d7_b = d2_b * iy * six
          do g_i_ = 1, g_p_
            g_s911(g_i_) = d4_b * g_length(g_i_) + d6_b * g_iy(g_i_) + d
     *7_b * g_e(g_i_)
          enddo
          s911 = d7_v
C--------
          d2_v = four * e
          d6_v = d2_v * iy / length
          d2_b = 1.0d0 / length
          d3_b = -d6_v / length
          d5_b = d2_b * d2_v
          d6_b = d2_b * iy * four
          do g_i_ = 1, g_p_
            g_s1111(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iy(g_i_) + 
     *d6_b * g_e(g_i_)
          enddo
          s1111 = d6_v
C--------
          specialcasey = .true.
        endif
C
C.....INITIALIZE IF [ALPHA_Z] IS ZERO
C
        if (alphaz .eq. zero) then
          d2_v = twelve * e
          d6_v = length * length
          d7_v = d6_v * length
          d8_v = d2_v * iz / d7_v
          d2_b = 1.0d0 / d7_v
          d3_b = -d8_v / d7_v
          d4_b = d3_b * length
          d5_b = d3_b * d6_v + d4_b * length + d4_b * length
          d7_b = d2_b * d2_v
          d8_b = d2_b * iz * twelve
          do g_i_ = 1, g_p_
            g_s22(g_i_) = d5_b * g_length(g_i_) + d7_b * g_iz(g_i_) + d8
     *_b * g_e(g_i_)
          enddo
          s22 = d8_v
C--------
          d2_v = six * e
          d6_v = length * length
          d7_v = d2_v * iz / d6_v
          d2_b = 1.0d0 / d6_v
          d3_b = -d7_v / d6_v
          d4_b = d3_b * length + d3_b * length
          d6_b = d2_b * d2_v
          d7_b = d2_b * iz * six
          do g_i_ = 1, g_p_
            g_s26(g_i_) = d4_b * g_length(g_i_) + d6_b * g_iz(g_i_) + d7
     *_b * g_e(g_i_)
          enddo
          s26 = d7_v
C--------
          d2_v = four * e
          d6_v = d2_v * iz / length
          d2_b = 1.0d0 / length
          d3_b = -d6_v / length
          d5_b = d2_b * d2_v
          d6_b = d2_b * iz * four
          do g_i_ = 1, g_p_
            g_s66(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iz(g_i_) + d6
     *_b * g_e(g_i_)
          enddo
          s66 = d6_v
C--------
          d2_v = -twelve * e
          d6_v = length * length
          d7_v = d6_v * length
          d8_v = d2_v * iz / d7_v
          d2_b = 1.0d0 / d7_v
          d3_b = -d8_v / d7_v
          d4_b = d3_b * length
          d5_b = d3_b * d6_v + d4_b * length + d4_b * length
          d7_b = d2_b * d2_v
          d8_b = d2_b * iz * (-twelve)
          do g_i_ = 1, g_p_
            g_s28(g_i_) = d5_b * g_length(g_i_) + d7_b * g_iz(g_i_) + d8
     *_b * g_e(g_i_)
          enddo
          s28 = d8_v
C--------
          d2_v = -six * e
          d6_v = length * length
          d7_v = d2_v * iz / d6_v
          d2_b = 1.0d0 / d6_v
          d3_b = -d7_v / d6_v
          d4_b = d3_b * length + d3_b * length
          d6_b = d2_b * d2_v
          d7_b = d2_b * iz * (-six)
          do g_i_ = 1, g_p_
            g_s68(g_i_) = d4_b * g_length(g_i_) + d6_b * g_iz(g_i_) + d7
     *_b * g_e(g_i_)
          enddo
          s68 = d7_v
C--------
          d2_v = twelve * e
          d6_v = length * length
          d7_v = d6_v * length
          d8_v = d2_v * iz / d7_v
          d2_b = 1.0d0 / d7_v
          d3_b = -d8_v / d7_v
          d4_b = d3_b * length
          d5_b = d3_b * d6_v + d4_b * length + d4_b * length
          d7_b = d2_b * d2_v
          d8_b = d2_b * iz * twelve
          do g_i_ = 1, g_p_
            g_s88(g_i_) = d5_b * g_length(g_i_) + d7_b * g_iz(g_i_) + d8
     *_b * g_e(g_i_)
          enddo
          s88 = d8_v
C--------
          d2_v = six * e
          d6_v = length * length
          d7_v = d2_v * iz / d6_v
          d2_b = 1.0d0 / d6_v
          d3_b = -d7_v / d6_v
          d4_b = d3_b * length + d3_b * length
          d6_b = d2_b * d2_v
          d7_b = d2_b * iz * six
          do g_i_ = 1, g_p_
            g_s212(g_i_) = d4_b * g_length(g_i_) + d6_b * g_iz(g_i_) + d
     *7_b * g_e(g_i_)
          enddo
          s212 = d7_v
C--------
          d2_v = two * e
          d6_v = d2_v * iz / length
          d2_b = 1.0d0 / length
          d3_b = -d6_v / length
          d5_b = d2_b * d2_v
          d6_b = d2_b * iz * two
          do g_i_ = 1, g_p_
            g_s612(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iz(g_i_) + d
     *6_b * g_e(g_i_)
          enddo
          s612 = d6_v
C--------
          d2_v = -six * e
          d6_v = length * length
          d7_v = d2_v * iz / d6_v
          d2_b = 1.0d0 / d6_v
          d3_b = -d7_v / d6_v
          d4_b = d3_b * length + d3_b * length
          d6_b = d2_b * d2_v
          d7_b = d2_b * iz * (-six)
          do g_i_ = 1, g_p_
            g_s812(g_i_) = d4_b * g_length(g_i_) + d6_b * g_iz(g_i_) + d
     *7_b * g_e(g_i_)
          enddo
          s812 = d7_v
C--------
          d2_v = four * e
          d6_v = d2_v * iz / length
          d2_b = 1.0d0 / length
          d3_b = -d6_v / length
          d5_b = d2_b * d2_v
          d6_b = d2_b * iz * four
          do g_i_ = 1, g_p_
            g_s1212(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iz(g_i_) + 
     *d6_b * g_e(g_i_)
          enddo
          s1212 = d6_v
C--------
          specialcasez = .true.
        endif
C
C.....INITIALIZE IF [ALPHA_Y] IS NONZERO BUT [A] IS ZERO
C
        if ((alphay .ne. zero) .and. (a .eq. zero)) then
          d5_v = e * iy / length
          d2_b = 1.0d0 / length
          d3_b = -d5_v / length
          d4_b = d2_b * iy
          d5_b = d2_b * e
          do g_i_ = 1, g_p_
            g_s55(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iy(g_i_) + d4
     *_b * g_e(g_i_)
          enddo
          s55 = d5_v
C--------
          d6_v = -e * iy / length
          d2_b = 1.0d0 / length
          d3_b = -d6_v / length
          d5_b = d2_b * (-e)
          d6_b = -(d2_b * iy)
          do g_i_ = 1, g_p_
            g_s511(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iy(g_i_) + d
     *6_b * g_e(g_i_)
          enddo
          s511 = d6_v
C--------
          d5_v = e * iy / length
          d2_b = 1.0d0 / length
          d3_b = -d5_v / length
          d4_b = d2_b * iy
          d5_b = d2_b * e
          do g_i_ = 1, g_p_
            g_s1111(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iy(g_i_) + 
     *d4_b * g_e(g_i_)
          enddo
          s1111 = d5_v
C--------
          specialcasey = .true.
        endif
C
C.....INITIALIZE IF [ALPHA_Z] IS NONZERO BUT [A] IS ZERO
C
        if ((alphaz .ne. zero) .and. (a .eq. zero)) then
          d5_v = e * iz / length
          d2_b = 1.0d0 / length
          d3_b = -d5_v / length
          d4_b = d2_b * iz
          d5_b = d2_b * e
          do g_i_ = 1, g_p_
            g_s66(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iz(g_i_) + d4
     *_b * g_e(g_i_)
          enddo
          s66 = d5_v
C--------
          d6_v = -e * iz / length
          d2_b = 1.0d0 / length
          d3_b = -d6_v / length
          d5_b = d2_b * (-e)
          d6_b = -(d2_b * iz)
          do g_i_ = 1, g_p_
            g_s612(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iz(g_i_) + d
     *6_b * g_e(g_i_)
          enddo
          s612 = d6_v
C--------
          d5_v = e * iz / length
          d2_b = 1.0d0 / length
          d3_b = -d5_v / length
          d4_b = d2_b * iz
          d5_b = d2_b * e
          do g_i_ = 1, g_p_
            g_s1212(g_i_) = d3_b * g_length(g_i_) + d5_b * g_iz(g_i_) + 
     *d4_b * g_e(g_i_)
          enddo
          s1212 = d5_v
C--------
          specialcasez = .true.
        endif
C
C.....INITIALIZE IF [I_Y] IS ZERO
C
        if (iy .eq. zero) then
          specialcasey = .true.
        endif
C
C.....INITIALIZE IF [I_Z] IS ZERO
C
        if (iz .eq. zero) then
          specialcasez = .true.
        endif
C
C.....INITIALIZE THE STIFFNESS ENTRIES IN ALL OTHER CASES
C
        d5_v = e * a / length
        d2_b = 1.0d0 / length
        d3_b = -d5_v / length
        d4_b = d2_b * a
        d5_b = d2_b * e
        do g_i_ = 1, g_p_
          g_s11(g_i_) = d3_b * g_length(g_i_) + d5_b * g_a(g_i_) + d4_b 
     ** g_e(g_i_)
        enddo
        s11 = d5_v
C--------
        d6_v = -e * a / length
        d2_b = 1.0d0 / length
        d3_b = -d6_v / length
        d5_b = d2_b * (-e)
        d6_b = -(d2_b * a)
        do g_i_ = 1, g_p_
          g_s17(g_i_) = d3_b * g_length(g_i_) + d5_b * g_a(g_i_) + d6_b 
     ** g_e(g_i_)
        enddo
        s17 = d6_v
C--------
        d5_v = e * a / length
        d2_b = 1.0d0 / length
        d3_b = -d5_v / length
        d4_b = d2_b * a
        d5_b = d2_b * e
        do g_i_ = 1, g_p_
          g_s77(g_i_) = d3_b * g_length(g_i_) + d5_b * g_a(g_i_) + d4_b 
     ** g_e(g_i_)
        enddo
        s77 = d5_v
C--------
C
        if (jj .eq. zero) then
          do g_i_ = 1, g_p_
            g_s44(g_i_) = 0.0d0
          enddo
          s44 = zero
C--------
          do g_i_ = 1, g_p_
            g_s410(g_i_) = 0.0d0
          enddo
          s410 = zero
C--------
          do g_i_ = 1, g_p_
            g_s1010(g_i_) = 0.0d0
          enddo
          s1010 = zero
C--------
        else
          if (c1 .eq. zero) then
            d5_v = g * jj / length
            d2_b = 1.0d0 / length
            d3_b = -d5_v / length
            d4_b = d2_b * jj
            d5_b = d2_b * g
            do g_i_ = 1, g_p_
              g_s44(g_i_) = d3_b * g_length(g_i_) + d5_b * g_jj(g_i_) + 
     *d4_b * g_g(g_i_)
            enddo
            s44 = d5_v
C--------
            d6_v = -g * jj / length
            d2_b = 1.0d0 / length
            d3_b = -d6_v / length
            d5_b = d2_b * (-g)
            d6_b = -(d2_b * jj)
            do g_i_ = 1, g_p_
              g_s410(g_i_) = d3_b * g_length(g_i_) + d5_b * g_jj(g_i_) +
     * d6_b * g_g(g_i_)
            enddo
            s410 = d6_v
C--------
            d5_v = g * jj / length
            d2_b = 1.0d0 / length
            d3_b = -d5_v / length
            d4_b = d2_b * jj
            d5_b = d2_b * g
            do g_i_ = 1, g_p_
              g_s1010(g_i_) = d3_b * g_length(g_i_) + d5_b * g_jj(g_i_) 
     *+ d4_b * g_g(g_i_)
            enddo
            s1010 = d5_v
C--------
          else
            d5_v = g * jj / c1
            d2_b = 1.0d0 / c1
            d3_b = -d5_v / c1
            d4_b = d2_b * jj
            d5_b = d2_b * g
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = d3_b * g_c1(g_i_) + d5_b * g_jj(g_i_) + d4_
     *b * g_g(g_i_)
            enddo
            d1_w = d5_v
            d2_v = half * length
            d4_v = sqrt(d1_w)
            if ( d1_w .gt. 0.0d0 ) then
               d1_p = 1.0d0 / (2.0d0 *  d4_v)
            else
               call ehufDO (9,d1_w, d4_v, d1_p,
     +g_ehfid,
     +1095)
            endif
            d4_b = d2_v * d1_p
            d5_b = d4_v * half
            do g_i_ = 1, g_p_
              g_b(g_i_) = d4_b * g_d1_w(g_i_) + d5_b * g_length(g_i_)
            enddo
            b = d2_v * d4_v
C--------
            d6_v = tanh (b)
            d1_p = 1.0d0 - ( d6_v *  d6_v)
            d7_v = d6_v / b
            d8_v = one - d7_v
            d9_v = length * d8_v
            d10_v = g * jj / d9_v
            d2_b = 1.0d0 / d9_v
            d3_b = -d10_v / d9_v
            d4_b = d3_b * d8_v
            d6_b = -(d3_b * length)
            d8_b = d6_b * (-d7_v / b) + d6_b * (1.0d0 / b) * d1_p
            d9_b = d2_b * jj
            d10_b = d2_b * g
            do g_i_ = 1, g_p_
              g_s44(g_i_) = d8_b * g_b(g_i_) + d4_b * g_length(g_i_) + d
     *10_b * g_jj(g_i_) + d9_b * g_g(g_i_)
            enddo
            s44 = d10_v
C--------
            d7_v = tanh (b)
            d1_p = 1.0d0 - ( d7_v *  d7_v)
            d8_v = d7_v / b
            d9_v = one - d8_v
            d10_v = length * d9_v
            d11_v = -(g * jj) / d10_v
            d3_b = -d11_v / d10_v
            d4_b = d3_b * d9_v
            d6_b = -(d3_b * length)
            d8_b = d6_b * (-d8_v / b) + d6_b * (1.0d0 / b) * d1_p
            d9_b = -(1.0d0 / d10_v)
            d10_b = d9_b * jj
            d11_b = d9_b * g
            do g_i_ = 1, g_p_
              g_s410(g_i_) = d8_b * g_b(g_i_) + d4_b * g_length(g_i_) + 
     *d11_b * g_jj(g_i_) + d10_b * g_g(g_i_)
            enddo
            s410 = d11_v
C--------
            d6_v = tanh (b)
            d1_p = 1.0d0 - ( d6_v *  d6_v)
            d7_v = d6_v / b
            d8_v = one - d7_v
            d9_v = length * d8_v
            d10_v = g * jj / d9_v
            d2_b = 1.0d0 / d9_v
            d3_b = -d10_v / d9_v
            d4_b = d3_b * d8_v
            d6_b = -(d3_b * length)
            d8_b = d6_b * (-d7_v / b) + d6_b * (1.0d0 / b) * d1_p
            d9_b = d2_b * jj
            d10_b = d2_b * g
            do g_i_ = 1, g_p_
              g_s1010(g_i_) = d8_b * g_b(g_i_) + d4_b * g_length(g_i_) +
     * d10_b * g_jj(g_i_) + d9_b * g_g(g_i_)
            enddo
            s1010 = d10_v
C--------
          endif
        endif
C
        if (.not. specialcasey) then
          d3_v = g * a
          d5_v = d3_v * length
          d4_b = length * length
          d3_b = d5_v + length * d3_v
          d5_b = d4_b * a
          d6_b = d4_b * g
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_length(g_i_) + d6_b * g_a(g_i_) + d5
     *_b * g_g(g_i_)
          enddo
          d1_w = d5_v * length
          d3_v = twelve * e
          d5_v = d3_v * iy
          d6_v = d1_w / d5_v
          d2_b = 1.0d0 / d5_v
          d3_b = -d6_v / d5_v
          d5_b = d3_b * d3_v
          d6_b = d3_b * iy * twelve
          do g_i_ = 1, g_p_
            g_gammay(g_i_) = d5_b * g_iy(g_i_) + d6_b * g_e(g_i_) + d2_b
     * * g_d1_w(g_i_)
          enddo
          gammay = d6_v
C--------
          d5_v = twelve * gammay
          d6_v = alphay + gammay
          d7_v = d5_v * d6_v
          d8_v = (alphay + four * gammay) / d7_v
          d2_b = 1.0d0 / d7_v
          d3_b = -d8_v / d7_v
          d5_b = d3_b * d5_v
          d6_b = d5_b + d2_b
          d7_b = d5_b + d3_b * d6_v * twelve + d2_b * four
          do g_i_ = 1, g_p_
            g_bendy(g_i_) = d7_b * g_gammay(g_i_) + d6_b * g_alphay(g_i_
     *)
          enddo
          bendy = d8_v
C--------
          d6_v = twelve * gammay
          d7_v = alphay + gammay
          d8_v = d6_v * d7_v
          d9_v = (-alphay + two * gammay) / d8_v
          d2_b = 1.0d0 / d8_v
          d3_b = -d9_v / d8_v
          d5_b = d3_b * d6_v
          d7_b = d5_b + d3_b * d7_v * twelve + d2_b * two
          d6_b = d5_b + (-d2_b)
          do g_i_ = 1, g_p_
            g_bendcy(g_i_) = d7_b * g_gammay(g_i_) + d6_b * g_alphay(g_i
     *_)
          enddo
          bendcy = d9_v
C--------
          d7_v = alphay + gammay
          d8_v = length * d7_v
          d9_v = g * a / d8_v
          d2_b = 1.0d0 / d8_v
          d3_b = -d9_v / d8_v
          d4_b = d3_b * d7_v
          d5_b = d3_b * length
          d8_b = d2_b * a
          d9_b = d2_b * g
          do g_i_ = 1, g_p_
            g_s33(g_i_) = d5_b * g_gammay(g_i_) + d5_b * g_alphay(g_i_) 
     *+ d4_b * g_length(g_i_) + d9_b * g_a(g_i_) + d8_b * g_g(g_i_)
          enddo
          s33 = d9_v
C--------
          d8_v = two * (alphay + gammay)
          d9_v = -g * a / d8_v
          d2_b = 1.0d0 / d8_v
          d4_b = -d9_v / d8_v * two
          d8_b = d2_b * (-g)
          d9_b = -(d2_b * a)
          do g_i_ = 1, g_p_
            g_s35(g_i_) = d4_b * g_gammay(g_i_) + d4_b * g_alphay(g_i_) 
     *+ d8_b * g_a(g_i_) + d9_b * g_g(g_i_)
          enddo
          s35 = d9_v
C--------
          d3_v = bendy * g
          d5_v = d3_v * a
          d4_b = length * a
          d5_b = length * d3_v
          d6_b = d4_b * g
          d7_b = d4_b * bendy
          do g_i_ = 1, g_p_
            g_s55(g_i_) = d5_v * g_length(g_i_) + d5_b * g_a(g_i_) + d7_
     *b * g_g(g_i_) + d6_b * g_bendy(g_i_)
          enddo
          s55 = d5_v * length
C--------
          d8_v = alphay + gammay
          d9_v = length * d8_v
          d10_v = -g * a / d9_v
          d2_b = 1.0d0 / d9_v
          d3_b = -d10_v / d9_v
          d4_b = d3_b * d8_v
          d5_b = d3_b * length
          d9_b = d2_b * (-g)
          d10_b = -(d2_b * a)
          do g_i_ = 1, g_p_
            g_s39(g_i_) = d5_b * g_gammay(g_i_) + d5_b * g_alphay(g_i_) 
     *+ d4_b * g_length(g_i_) + d9_b * g_a(g_i_) + d10_b * g_g(g_i_)
          enddo
          s39 = d10_v
C--------
          d7_v = two * (alphay + gammay)
          d8_v = g * a / d7_v
          d2_b = 1.0d0 / d7_v
          d4_b = -d8_v / d7_v * two
          d7_b = d2_b * a
          d8_b = d2_b * g
          do g_i_ = 1, g_p_
            g_s59(g_i_) = d4_b * g_gammay(g_i_) + d4_b * g_alphay(g_i_) 
     *+ d8_b * g_a(g_i_) + d7_b * g_g(g_i_)
          enddo
          s59 = d8_v
C--------
          d7_v = alphay + gammay
          d8_v = length * d7_v
          d9_v = g * a / d8_v
          d2_b = 1.0d0 / d8_v
          d3_b = -d9_v / d8_v
          d4_b = d3_b * d7_v
          d5_b = d3_b * length
          d8_b = d2_b * a
          d9_b = d2_b * g
          do g_i_ = 1, g_p_
            g_s99(g_i_) = d5_b * g_gammay(g_i_) + d5_b * g_alphay(g_i_) 
     *+ d4_b * g_length(g_i_) + d9_b * g_a(g_i_) + d8_b * g_g(g_i_)
          enddo
          s99 = d9_v
C--------
          d8_v = two * (alphay + gammay)
          d9_v = -g * a / d8_v
          d2_b = 1.0d0 / d8_v
          d4_b = -d9_v / d8_v * two
          d8_b = d2_b * (-g)
          d9_b = -(d2_b * a)
          do g_i_ = 1, g_p_
            g_s311(g_i_) = d4_b * g_gammay(g_i_) + d4_b * g_alphay(g_i_)
     * + d8_b * g_a(g_i_) + d9_b * g_g(g_i_)
          enddo
          s311 = d9_v
C--------
          d3_v = bendcy * g
          d5_v = d3_v * a
          d4_b = length * a
          d5_b = length * d3_v
          d6_b = d4_b * g
          d7_b = d4_b * bendcy
          do g_i_ = 1, g_p_
            g_s511(g_i_) = d5_v * g_length(g_i_) + d5_b * g_a(g_i_) + d7
     *_b * g_g(g_i_) + d6_b * g_bendcy(g_i_)
          enddo
          s511 = d5_v * length
C--------
          d7_v = two * (alphay + gammay)
          d8_v = g * a / d7_v
          d2_b = 1.0d0 / d7_v
          d4_b = -d8_v / d7_v * two
          d7_b = d2_b * a
          d8_b = d2_b * g
          do g_i_ = 1, g_p_
            g_s911(g_i_) = d4_b * g_gammay(g_i_) + d4_b * g_alphay(g_i_)
     * + d8_b * g_a(g_i_) + d7_b * g_g(g_i_)
          enddo
          s911 = d8_v
C--------
          d3_v = bendy * g
          d5_v = d3_v * a
          d4_b = length * a
          d5_b = length * d3_v
          d6_b = d4_b * g
          d7_b = d4_b * bendy
          do g_i_ = 1, g_p_
            g_s1111(g_i_) = d5_v * g_length(g_i_) + d5_b * g_a(g_i_) + d
     *7_b * g_g(g_i_) + d6_b * g_bendy(g_i_)
          enddo
          s1111 = d5_v * length
C--------
        endif
C
        if (.not. specialcasez) then
          d3_v = g * a
          d5_v = d3_v * length
          d4_b = length * length
          d3_b = d5_v + length * d3_v
          d5_b = d4_b * a
          d6_b = d4_b * g
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d3_b * g_length(g_i_) + d6_b * g_a(g_i_) + d5
     *_b * g_g(g_i_)
          enddo
          d1_w = d5_v * length
          d3_v = twelve * e
          d5_v = d3_v * iz
          d6_v = d1_w / d5_v
          d2_b = 1.0d0 / d5_v
          d3_b = -d6_v / d5_v
          d5_b = d3_b * d3_v
          d6_b = d3_b * iz * twelve
          do g_i_ = 1, g_p_
            g_gammaz(g_i_) = d5_b * g_iz(g_i_) + d6_b * g_e(g_i_) + d2_b
     * * g_d1_w(g_i_)
          enddo
          gammaz = d6_v
C--------
          d5_v = twelve * gammaz
          d6_v = alphaz + gammaz
          d7_v = d5_v * d6_v
          d8_v = (alphaz + four * gammaz) / d7_v
          d2_b = 1.0d0 / d7_v
          d3_b = -d8_v / d7_v
          d5_b = d3_b * d5_v
          d6_b = d5_b + d2_b
          d7_b = d5_b + d3_b * d6_v * twelve + d2_b * four
          do g_i_ = 1, g_p_
            g_bendz(g_i_) = d7_b * g_gammaz(g_i_) + d6_b * g_alphaz(g_i_
     *)
          enddo
          bendz = d8_v
C--------
          d6_v = twelve * gammaz
          d7_v = alphaz + gammaz
          d8_v = d6_v * d7_v
          d9_v = (-alphaz + two * gammaz) / d8_v
          d2_b = 1.0d0 / d8_v
          d3_b = -d9_v / d8_v
          d5_b = d3_b * d6_v
          d7_b = d5_b + d3_b * d7_v * twelve + d2_b * two
          d6_b = d5_b + (-d2_b)
          do g_i_ = 1, g_p_
            g_bendcz(g_i_) = d7_b * g_gammaz(g_i_) + d6_b * g_alphaz(g_i
     *_)
          enddo
          bendcz = d9_v
C--------
          d7_v = alphaz + gammaz
          d8_v = length * d7_v
          d9_v = g * a / d8_v
          d2_b = 1.0d0 / d8_v
          d3_b = -d9_v / d8_v
          d4_b = d3_b * d7_v
          d5_b = d3_b * length
          d8_b = d2_b * a
          d9_b = d2_b * g
          do g_i_ = 1, g_p_
            g_s22(g_i_) = d5_b * g_gammaz(g_i_) + d5_b * g_alphaz(g_i_) 
     *+ d4_b * g_length(g_i_) + d9_b * g_a(g_i_) + d8_b * g_g(g_i_)
          enddo
          s22 = d9_v
C--------
          d7_v = two * (alphaz + gammaz)
          d8_v = g * a / d7_v
          d2_b = 1.0d0 / d7_v
          d4_b = -d8_v / d7_v * two
          d7_b = d2_b * a
          d8_b = d2_b * g
          do g_i_ = 1, g_p_
            g_s26(g_i_) = d4_b * g_gammaz(g_i_) + d4_b * g_alphaz(g_i_) 
     *+ d8_b * g_a(g_i_) + d7_b * g_g(g_i_)
          enddo
          s26 = d8_v
C--------
          d3_v = bendz * g
          d5_v = d3_v * a
          d4_b = length * a
          d5_b = length * d3_v
          d6_b = d4_b * g
          d7_b = d4_b * bendz
          do g_i_ = 1, g_p_
            g_s66(g_i_) = d5_v * g_length(g_i_) + d5_b * g_a(g_i_) + d7_
     *b * g_g(g_i_) + d6_b * g_bendz(g_i_)
          enddo
          s66 = d5_v * length
C--------
          d8_v = alphaz + gammaz
          d9_v = length * d8_v
          d10_v = -g * a / d9_v
          d2_b = 1.0d0 / d9_v
          d3_b = -d10_v / d9_v
          d4_b = d3_b * d8_v
          d5_b = d3_b * length
          d9_b = d2_b * (-g)
          d10_b = -(d2_b * a)
          do g_i_ = 1, g_p_
            g_s28(g_i_) = d5_b * g_gammaz(g_i_) + d5_b * g_alphaz(g_i_) 
     *+ d4_b * g_length(g_i_) + d9_b * g_a(g_i_) + d10_b * g_g(g_i_)
          enddo
          s28 = d10_v
C--------
          d8_v = two * (alphaz + gammaz)
          d9_v = -g * a / d8_v
          d2_b = 1.0d0 / d8_v
          d4_b = -d9_v / d8_v * two
          d8_b = d2_b * (-g)
          d9_b = -(d2_b * a)
          do g_i_ = 1, g_p_
            g_s68(g_i_) = d4_b * g_gammaz(g_i_) + d4_b * g_alphaz(g_i_) 
     *+ d8_b * g_a(g_i_) + d9_b * g_g(g_i_)
          enddo
          s68 = d9_v
C--------
          d7_v = alphaz + gammaz
          d8_v = length * d7_v
          d9_v = g * a / d8_v
          d2_b = 1.0d0 / d8_v
          d3_b = -d9_v / d8_v
          d4_b = d3_b * d7_v
          d5_b = d3_b * length
          d8_b = d2_b * a
          d9_b = d2_b * g
          do g_i_ = 1, g_p_
            g_s88(g_i_) = d5_b * g_gammaz(g_i_) + d5_b * g_alphaz(g_i_) 
     *+ d4_b * g_length(g_i_) + d9_b * g_a(g_i_) + d8_b * g_g(g_i_)
          enddo
          s88 = d9_v
C--------
          d7_v = two * (alphaz + gammaz)
          d8_v = g * a / d7_v
          d2_b = 1.0d0 / d7_v
          d4_b = -d8_v / d7_v * two
          d7_b = d2_b * a
          d8_b = d2_b * g
          do g_i_ = 1, g_p_
            g_s212(g_i_) = d4_b * g_gammaz(g_i_) + d4_b * g_alphaz(g_i_)
     * + d8_b * g_a(g_i_) + d7_b * g_g(g_i_)
          enddo
          s212 = d8_v
C--------
          d3_v = bendcz * g
          d5_v = d3_v * a
          d4_b = length * a
          d5_b = length * d3_v
          d6_b = d4_b * g
          d7_b = d4_b * bendcz
          do g_i_ = 1, g_p_
            g_s612(g_i_) = d5_v * g_length(g_i_) + d5_b * g_a(g_i_) + d7
     *_b * g_g(g_i_) + d6_b * g_bendcz(g_i_)
          enddo
          s612 = d5_v * length
C--------
          d8_v = two * (alphaz + gammaz)
          d9_v = -g * a / d8_v
          d2_b = 1.0d0 / d8_v
          d4_b = -d9_v / d8_v * two
          d8_b = d2_b * (-g)
          d9_b = -(d2_b * a)
          do g_i_ = 1, g_p_
            g_s812(g_i_) = d4_b * g_gammaz(g_i_) + d4_b * g_alphaz(g_i_)
     * + d8_b * g_a(g_i_) + d9_b * g_g(g_i_)
          enddo
          s812 = d9_v
C--------
          d3_v = bendz * g
          d5_v = d3_v * a
          d4_b = length * a
          d5_b = length * d3_v
          d6_b = d4_b * g
          d7_b = d4_b * bendz
          do g_i_ = 1, g_p_
            g_s1212(g_i_) = d5_v * g_length(g_i_) + d5_b * g_a(g_i_) + d
     *7_b * g_g(g_i_) + d6_b * g_bendz(g_i_)
          enddo
          s1212 = d5_v * length
C--------
        endif
C
C.....INITIALIZE THE LOCAL STIFFNESS MATRIX
C
        do g_i_ = 1, g_p_
          g_locke(g_i_, 1, 1) = g_s11(g_i_)
        enddo
        locke(1, 1) = s11
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 2, 2) = g_s22(g_i_)
        enddo
        locke(2, 2) = s22
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 3, 3) = g_s33(g_i_)
        enddo
        locke(3, 3) = s33
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 2, 6) = g_s26(g_i_)
        enddo
        locke(2, 6) = s26
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 3, 5) = g_s35(g_i_)
        enddo
        locke(3, 5) = s35
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 4, 4) = g_s44(g_i_)
        enddo
        locke(4, 4) = s44
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 5, 5) = g_s55(g_i_)
        enddo
        locke(5, 5) = s55
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 6, 6) = g_s66(g_i_)
        enddo
        locke(6, 6) = s66
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 1, 7) = g_s17(g_i_)
        enddo
        locke(1, 7) = s17
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 2, 8) = g_s28(g_i_)
        enddo
        locke(2, 8) = s28
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 3, 9) = g_s39(g_i_)
        enddo
        locke(3, 9) = s39
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 5, 9) = g_s59(g_i_)
        enddo
        locke(5, 9) = s59
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 6, 8) = g_s68(g_i_)
        enddo
        locke(6, 8) = s68
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 7, 7) = g_s77(g_i_)
        enddo
        locke(7, 7) = s77
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 8, 8) = g_s88(g_i_)
        enddo
        locke(8, 8) = s88
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 9, 9) = g_s99(g_i_)
        enddo
        locke(9, 9) = s99
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 2, 12) = g_s212(g_i_)
        enddo
        locke(2, 12) = s212
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 3, 11) = g_s311(g_i_)
        enddo
        locke(3, 11) = s311
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 4, 10) = g_s410(g_i_)
        enddo
        locke(4, 10) = s410
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 5, 11) = g_s511(g_i_)
        enddo
        locke(5, 11) = s511
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 6, 12) = g_s612(g_i_)
        enddo
        locke(6, 12) = s612
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 8, 12) = g_s812(g_i_)
        enddo
        locke(8, 12) = s812
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 9, 11) = g_s911(g_i_)
        enddo
        locke(9, 11) = s911
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 10, 10) = g_s1010(g_i_)
        enddo
        locke(10, 10) = s1010
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 11, 11) = g_s1111(g_i_)
        enddo
        locke(11, 11) = s1111
C--------
        do g_i_ = 1, g_p_
          g_locke(g_i_, 12, 12) = g_s1212(g_i_)
        enddo
        locke(12, 12) = s1212
C--------
C
        do 99985 j = 1, 11
          do 99986 i = (j + 1), 12
            do g_i_ = 1, g_p_
              g_locke(g_i_, i, j) = g_locke(g_i_, j, i)
            enddo
            locke(i, j) = locke(j, i)
C--------
3002        continue
99986     continue
3001      continue
99985   continue
C
C.....ASSEMBLE THE ELEMENTAL STIFFNESS MATRIX
C.....
C..... [ESTIF] = [Le] * [LOCKE] * [Le]^T
C.....
C
C
C..... KHP: ask Michel about this 12x12x12x12 loop
C..... COULDN'T IT BE DONE MUCH FASTER?
C
        do 99981 l = 1, 12
          do 99982 k = 1, 12
            do g_i_ = 1, g_p_
              g_xke(g_i_) = g_locke(g_i_, k, l)
            enddo
            xke = locke(k, l)
C--------
            do 99983 j = 1, 12
              do 99984 i = 1, 12
                d4_v = le(i, k) * xke
                d6_b = le(j, l) * xke
                d7_b = le(j, l) * le(i, k)
                do g_i_ = 1, g_p_
                  g_estif(g_i_, i, j) = d4_v * g_le(g_i_, j, l) + d7_b *
     * g_xke(g_i_) + d6_b * g_le(g_i_, i, k) + g_estif(g_i_, i, j)
                enddo
                estif(i, j) = estif(i, j) + d4_v * le(j, l)
C--------
4004            continue
99984         continue
4003          continue
99983       continue
4002        continue
99982     continue
4001      continue
99981   continue
C
C.....ROTATE THE GLOBAL DISPLACEMENT VECTOR BACK TO LOCAL FRAME
C
C.....dgemv is a BLAS routine
C     dgemv( trans,m,n,alpha,a,lda,x,incx,beta,y,incy )
C
C      call dgemv('N',12,12,1.0d00,Le,12,Ug,1,0.0d00,Ue,1)
        do l = 1, 12
          do k = 1, 12
            do g_i_ = 1, g_p_
              g_ue(g_i_, l) = le(l, k) * g_ug(g_i_, k) + ug(k) * g_le(g_
     *i_, l, k) + g_ue(g_i_, l)
            enddo
            ue(l) = ue(l) + le(l, k) * ug(k)
C--------
          enddo
        enddo
C
C.....COMPUTE THE INTERNAL FORCE RESULTANTS
C
C      call dgemv('N',12,12,1.0d00,estif,12,Ue,1,0.0d00,Fe,1)
        do l = 1, 12
          do k = 1, 12
            do g_i_ = 1, g_p_
              g_fe(g_i_, l) = estif(l, k) * g_ue(g_i_, k) + ue(k) * g_es
     *tif(g_i_, l, k) + g_fe(g_i_, l)
            enddo
            fe(l) = fe(l) + estif(l, k) * ue(k)
C--------
          enddo
        enddo
C
C.....WRITE THE RESULTANTS INTO THE STRESS ARRAY WHERE:
C.....FORCE_X  = STRESSXX
C.....FORCE_Y  = STRESSYY
C.....FORCE_Z  = STRESSZZ
C.....MOMENT_X = STRESSXY
C.....MOMENT_Y = STRESSXZ
C.....MOMENT_Z = STRESSYZ
C
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 1, 1) = g_fe(g_i_, 1)
        enddo
c        stress(elm, 1, 1) = fe(1)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 2, 1) = g_fe(g_i_, 2)
        enddo
c        stress(elm, 2, 1) = fe(2)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 3, 1) = g_fe(g_i_, 3)
        enddo
c        stress(elm, 3, 1) = fe(3)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 4, 1) = g_fe(g_i_, 4)
        enddo
c        stress(elm, 4, 1) = fe(4)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 5, 1) = g_fe(g_i_, 5)
        enddo
c        stress(elm, 5, 1) = fe(5)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 6, 1) = g_fe(g_i_, 6)
        enddo
c        stress(elm, 6, 1) = fe(6)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 1, 2) = g_fe(g_i_, 7)
        enddo
c        stress(elm, 1, 2) = fe(7)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 2, 2) = g_fe(g_i_, 8)
        enddo
c        stress(elm, 2, 2) = fe(8)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 3, 2) = g_fe(g_i_, 9)
        enddo
c        stress(elm, 3, 2) = fe(9)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 4, 2) = g_fe(g_i_, 10)
        enddo
c        stress(elm, 4, 2) = fe(10)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 5, 2) = g_fe(g_i_, 11)
        enddo
c        stress(elm, 5, 2) = fe(11)
C--------
        do g_i_ = 1, g_p_
          g_stress(g_i_, elm, 6, 2) = g_fe(g_i_, 12)
        enddo
c        stress(elm, 6, 2) = fe(12)
C--------
C
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
        write (*, *) '*** FATAL ERROR in Routine SANDS7 ***'
        write (*, *) '*** The Timoschenko Beam Element  ***'
        write (*, *) '*** Has Zero Length!              ***'
        write (*, *) '*** ... All Treatments Terminated ***'
        stop
C
C.....ERROR-MESSAGE IF A MATERIAL/GEOMETRICAL PROPERTY IS WRONG
C
200     continue
        write (*, *) '*** FATAL ERROR in Routine SANDS7      ***'
        write (*, *) '*** A Material or Geometrical Property ***'
        write (*, *) '*** is Not Positive: Check Input Data  ***'
        write (*, *) '*** ... All Treatments Terminated      ***'
        stop
C
C.....ERROR-MESSAGE IF THE ROTATION MATRIX IS NOT ORTHOGONAL
C
300     continue
        write (*, *) '*** FATAL ERROR in Routine SANDS7      ***'
        write (*, *) '*** The Rotation Matrix Computed For a ***'
        write (*, *) '*** Timoshenko Beam is Not Orthogonal! ***'
        write (*, *) '*** ... All Treatments Terminated      ***'
        stop
C
C.....ERROR-MESSAGE IF THE ROTATION MATRIX IS NOT AVAILABLE
C.....FOR A BEAM WITH NON-SYMMETRIC CROSS-SECTIONAL AREA
C
400     continue
        write (*, *) '*** FATAL ERROR in Routine SANDS7    ***'
        write (*, *) '*** There is No Rotation Matrix      ***'
        write (*, *) '*** Available For a Timoshenko Beam  ***'
        write (*, *) '*** With Non-Symmetric Cross-section ***'
        write (*, *) '*** ... All Treatments Terminated    ***'
        stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "SANDS7"
C=======================C
