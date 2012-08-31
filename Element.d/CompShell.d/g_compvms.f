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
      subroutine gxcompvms(elm,type,numel,maxstr,maxgus,e,g_e,nu,h,g_h,
     *                     globalx,g_globalx,globaly,g_globaly,globalz,
     *                     g_globalz,globalu,g_globalu,stress,g_stress,
     *                     laysigm,nttco,nttly,ncmpfr,cmpco,g_cmpco,
     *                     idlay,mtlay,g_mtlay,cmpfr,maxgly,laysgid,
     *                     iatt,ctyp,catt,cfrm,nfely,msize,strainflg,
     *                     surface)
C=====================================================================C
C                                                                     C
C     -----------------                                               C
C     V A R I A B L E S                                               C
C     -----------------                                               C
C                                                                     C
C     elm      <input>   Finite Element Number                        C
C     type     <input>   Finite Element Type                          C
C     numel    <input>   Number of Finite Elements in the Mesh        C
C     maxstr   <input>   Maximum Number of Stresses                   C
C     maxgus   <input>   Maximum Number of Gauss Points for Stresses  C
C     E        <input>   Young Modulus (for an Isotropic Element)     C
C     nu       <input>   Poisson's Ratio (for an Isotropic Element)   C
C     h        <input>   Constant Thickness of the Element            C
C     globalX  <input>   X- Nodal Coordinates                         C
C     globalY  <input>   Y- Nodal Coordinates                         C
C     globalZ  <input>   Z- Nodal Coordinates                         C
C     globalU  <input>   Global Displacements at the Nodal Joints     C
C     stress   <output>  Stresses (Von Mises Stress) of the Element   C
C     laysigm  <output>  Stresses (Von Mises) per Layers of the Shell C
C     nttco    <input>   Number of Attributes of Type-1               C
C     nttly    <input>   Total Number of Layers for the Elements      C
C     ncmpfr   <input>   Total Number of Composite Frames             C
C     cmpco    <input>   Constitutive Coefficients                    C
C     idlay    <input>   Layer Identificators                         C
C     mtlay    <input>   Material Properties of the Layers            C
C     cmpfr    <input>   Composite Frames                             C
C     maxgly   <input>   Maximum Number of Gauss Points (Layers)      C
C     laysgid  <input>   Identificators for the Stresses per Layer    C
C     iatt     <input>   Attribute Number of the Element              C
C     ctyp     <input>   Type of Constitutive Law (0, 1, 2, or 3)     C
C     catt     <input>   Addressing for Composite Material Properties C
C     cfrm     <input>   Composite Frame Number                       C
C     nfely    <input>   Total Number of Composite Layers in the Mesh C
C     msize    <input>   Maximum Length of Stress Output Array        C
C                                                                     C
C=====================================================================C
C=Author   = Francois M. Hemez                                        C
C=Date     = June 10th, 1995                                          C
C=Version  = 2.0                                                      C
C=Modified = K. H. Pierson                                            C
C=Date     = April 11, 1997                                           C
C=Reason   = Added stress calculations for sigmaxx, sigmayy, sigmaxy  C
C            and von mises stress at top, median and bottom surfaces  C
C            Also added strain calculations for epsilonxx, epsilonyy, C
C            epsilonzz, epsilonxy and an equivalent strain at top,    C
C            median and bottom surfaces.                              C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....GLOBAL VARIABLES
C
        integer elm, type, numel, maxstr
        integer maxgus, maxgly, nttco, nttly, ncmpfr
        integer idlay(5, nttly), iatt, ctyp, catt, cfrm, surface
        integer laysgid(5, nfely), nfely, msize, strainflg
C
        real*8 ebar, epsxx, epsyy, epszz, epsxy, t2
        real*8 globalx(3), globaly(3), globalz(3), globalu(18)
        real*8 cmpco(36, nttco), mtlay(9, nttly), cmpfr(9, ncmpfr)
        real*8 stress(msize, maxstr, maxgus), e, nu, h(3)
        real*8 str(6), xp(3), yp(3), zp(3)
        real*8 xg(3), yg(3), zg(3)
        real*8 laysigm(nfely, maxstr, maxgly)
C
C.....SET THE MAXIMUM NUMBER OF LAYERS OF THE ELEMENT
C
        integer maxlayer
        parameter (maxlayer = 1000)
C
C.....LOCAL VARIABLES
C
        integer i, j, row, dimp, idcmp23(5, maxlayer)
        integer rowb(9), colb(9), rowm(9), colm(9), nlayer
        integer ilayer, layerpos
C
        real*8 area, twicearea, factor, alpha, clr, cqr
        real*8 llr(9, 3), lqr(9, 3), l(18, 3), p(18, 3)
        real*8 reducedl(9, 3), reducedp(9, 3)
        real*8 x(3), y(3), z(3), rot(6, 6), localu(18)
        real*8 x21, x32, x13, y21, y32, y13
        real*8 x12, x23, x31, y12, y23, y31
        real*8 x0, y0, dist12, dist23, dist31
        real*8 c12, c23, c31, s12, s23, s31
        real*8 cc12, cc23, cc31, ss12, ss23, ss31
        real*8 cs12, cs23, cs31, thick, x1, x2
        real*8 zero, one, two, three, six, xkj, xsj
        real*8 cstbb(3, 3), cstmm(3, 3), elem(3)
        real*8 cstbm(3, 3), cstmb(3, 3), elen(3)
        real*8 mtcmp23(8, maxlayer), cstcoef(36)
        real*8 elecrv(3), elestr(3), vonmises
        real*8 ups(3), mds(3), lws(3), upr, lwr, ups1, ups2
        real*8 mdr, mds1, mds2, mdvms
        real*8 lws1, lws2, upvms, lwvms, appxh2
        real*8 eframe(3, 3), aframe(3, 3)
C
        logical detecterror, fastcal
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_cmpco, ldg_mtlay, ldg_h, ldg_globalu, ld
     *g_stress, ldg_e, ldg_globalx, ldg_globaly
        integer ldg_globalz
c
c manually inserted - begin	     
c
        parameter (g_p_=1,ldg_e=1,ldg_h=1,ldg_globalx=1,ldg_globaly=1,
     *             ldg_globalz=1,ldg_globalu=1,ldg_stress=1,ldg_cmpco=1,
     *             ldg_mtlay=1)
c
c manually inserted - end      
c
        double precision d2_p, d9_b, d7_b, d2_v, d3_v, d4_v, d8_b, d2_b,
     * d3_b, d4_b
        double precision d1_p, d1_w, d6_v, d7_v, d5_b, d6_b, g_rot(g_pma
     *x_, 6, 6), g_llr(g_pmax_, 9, 3), g_lqr(g_pmax_, 9, 3), g_reducedl(
     *g_pmax_, 9, 3)
        double precision g_reducedp(g_pmax_, 9, 3), g_l(g_pmax_, 18, 3),
     * g_p(g_pmax_, 18, 3), g_x(g_pmax_, 3), g_y(g_pmax_, 3), g_cstbb(g_
     *pmax_, 3, 3), g_cstmm(g_pmax_, 3, 3), g_cstbm(g_pmax_, 3, 3), g_cs
     *tmb(g_pmax_, 3, 3), g_cstcoef(g_pmax_, 36)
        double precision g_mtcmp23(g_pmax_, 8, maxlayer), g_localu(g_pma
     *x_, 18), g_elecrv(g_pmax_, 3), g_elestr(g_pmax_, 3), g_elem(g_pmax
     *_, 3), g_elen(g_pmax_, 3), g_ups(g_pmax_, 3), g_lws(g_pmax_, 3), g
     *_eframe(g_pmax_, 3, 3), g_cmpco(ldg_cmpco, 36, nttco)
        double precision g_mtlay(ldg_mtlay, 9, nttly), g_thick(g_pmax_),
     * g_h(ldg_h, 3), g_appxh2(g_pmax_), g_x21(g_pmax_), g_x12(g_pmax_),
     * g_x32(g_pmax_), g_x23(g_pmax_), g_x13(g_pmax_), g_x31(g_pmax_)
        double precision g_y21(g_pmax_), g_y12(g_pmax_), g_y32(g_pmax_),
     * g_y23(g_pmax_), g_y13(g_pmax_), g_y31(g_pmax_), g_twicearea(g_pma
     *x_), g_area(g_pmax_), g_d1_w(g_pmax_), g_dist12(g_pmax_)
        double precision g_dist23(g_pmax_), g_dist31(g_pmax_), g_c12(g_p
     *max_), g_s12(g_pmax_), g_c23(g_pmax_), g_s23(g_pmax_), g_c31(g_pma
     *x_), g_s31(g_pmax_), g_cc12(g_pmax_), g_cc23(g_pmax_)
        double precision g_cc31(g_pmax_), g_ss12(g_pmax_), g_ss23(g_pmax
     *_), g_ss31(g_pmax_), g_cs12(g_pmax_), g_cs23(g_pmax_), g_cs31(g_pm
     *ax_), g_factor(g_pmax_), g_globalu(ldg_globalu, 18), g_t2(g_pmax_)
        double precision g_epsxx(g_pmax_), g_epsyy(g_pmax_), g_epsxy(g_p
     *max_), g_epszz(g_pmax_), g_str(g_pmax_, 6), g_stress(ldg_stress, m
     *size, maxstr, maxgus), g_ebar(g_pmax_), g_xkj(g_pmax_), g_xsj(g_pm
     *ax_), g_vonmises(g_pmax_)
        double precision g_x1(g_pmax_), g_x2(g_pmax_), g_upr(g_pmax_), g
     *_ups1(g_pmax_), g_ups2(g_pmax_), g_upvms(g_pmax_), g_lwr(g_pmax_),
     * g_lws1(g_pmax_), g_lws2(g_pmax_), g_lwvms(g_pmax_)
        double precision g_mds(g_pmax_, 3), g_mdr(g_pmax_), g_mds1(g_pma
     *x_), g_mds2(g_pmax_), g_mdvms(g_pmax_), g_e(ldg_e), g_globalx(ldg_
     *globalx, 3), g_globaly(ldg_globaly, 3), g_globalz(ldg_globalz, 3),
     * g_xp(g_pmax_, 3)
        double precision g_yp(g_pmax_, 3), g_zp(g_pmax_, 3)
        save g_mdr, g_mds1, g_mds2, g_mdvms, g_xp, g_yp, g_zp
        save g_x2, g_upr, g_ups1, g_ups2, g_upvms, g_lwr, g_lws1, g_lws2
     *, g_lwvms, g_mds
        save g_epsxx, g_epsyy, g_epsxy, g_epszz, g_str, g_ebar, g_xkj, g
     *_xsj, g_vonmises, g_x1
        save g_cc23, g_cc31, g_ss12, g_ss23, g_ss31, g_cs12, g_cs23, g_c
     *s31, g_factor, g_t2
        save g_dist12, g_dist23, g_dist31, g_c12, g_s12, g_c23, g_s23, g
     *_c31, g_s31, g_cc12
        save g_x31, g_y21, g_y12, g_y32, g_y23, g_y13, g_y31, g_twiceare
     *a, g_area, g_d1_w
        save g_ups, g_lws, g_eframe, g_thick, g_appxh2, g_x21, g_x12, g_
     *x32, g_x23, g_x13
        save g_cstmm, g_cstbm, g_cstmb, g_cstcoef, g_mtcmp23, g_localu, 
     *g_elecrv, g_elestr, g_elem, g_elen
        save g_rot, g_llr, g_lqr, g_reducedl, g_reducedp, g_l, g_p, g_x,
     * g_y, g_cstbb
        external g_complay
        external g_compcst
        external g_straineq
        external g_transform
        external g_compcrd2
        intrinsic dble
        data zero /0.000000d+00/
        data one /1.000000d+00/
        data two /2.000000d+00/
        data three /3.000000d+00/
        data six /6.000000d+00/
        data xg /1.0, 0.0, 0.0/
        data yg /0.0, 1.0, 0.0/
        data zg /0.0, 0.0, 1.0/
C
C.....INITIALIZE THE WEIGHTS FOR PURE BENDING CONTRIBUTION
C
        data clr /0.000000d+00/
        data cqr /1.000000d+00/
C
C.....INITIALIZE THE WEIGHT FOR PURE MEMBRANE CONTRIBUTION
C
        data alpha /1.500000d+00/
C
C.....INITIALIZE THE LOGICAL FOR ENFORCING FASTER COMPUTATIONS
C
        data fastcal /.true./
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE LOCAL MATRICES
C
C      do 1001 i=1,9
C         rowb(i) = 0
C         rowm(i) = 0
C         colb(i) = 0
C         colm(i) = 0
C 1001 continue
C
        integer g_ehfid
        data g_ehfid /0/
C

C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do 99998 j = 1, 6
          do 99999 i = 1, 6
            do g_i_ = 1, g_p_
              g_rot(g_i_, i, j) = 0.0d0
            enddo
            rot(i, j) = zero
C--------
1003        continue
99999     continue
1002      continue
99998   continue
C
        do 99991 j = 1, 3
          do 99997 i = 1, 9
            do g_i_ = 1, g_p_
              g_llr(g_i_, i, j) = 0.0d0
            enddo
            llr(i, j) = zero
C--------
1005        continue
99997     continue
          do 99996 i = 1, 9
            do g_i_ = 1, g_p_
              g_lqr(g_i_, i, j) = 0.0d0
            enddo
            lqr(i, j) = zero
C--------
1006        continue
99996     continue
          do 99995 i = 1, 9
            do g_i_ = 1, g_p_
              g_reducedl(g_i_, i, j) = 0.0d0
            enddo
            reducedl(i, j) = zero
C--------
1007        continue
99995     continue
          do 99994 i = 1, 9
            do g_i_ = 1, g_p_
              g_reducedp(g_i_, i, j) = 0.0d0
            enddo
            reducedp(i, j) = zero
C--------
1008        continue
99994     continue
          do 99993 i = 1, 18
            do g_i_ = 1, g_p_
              g_l(g_i_, i, j) = 0.0d0
            enddo
            l(i, j) = zero
C--------
1009        continue
99993     continue
          do 99992 i = 1, 18
            do g_i_ = 1, g_p_
              g_p(g_i_, i, j) = 0.0d0
            enddo
            p(i, j) = zero
C--------
1010        continue
99992     continue
1004      continue
99991   continue
C
        do 99990 i = 1, 3
          do g_i_ = 1, g_p_
            g_x(g_i_, i) = 0.0d0
          enddo
          x(i) = zero
C--------
          do g_i_ = 1, g_p_
            g_y(g_i_, i) = 0.0d0
          enddo
          y(i) = zero
C--------
          z(i) = zero
1011      continue
99990   continue
C
        do 99985 j = 1, 3
          do 99989 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstbb(g_i_, i, j) = 0.0d0
            enddo
            cstbb(i, j) = zero
C--------
1013        continue
99989     continue
          do 99988 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstmm(g_i_, i, j) = 0.0d0
            enddo
            cstmm(i, j) = zero
C--------
1014        continue
99988     continue
          do 99987 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstbm(g_i_, i, j) = 0.0d0
            enddo
            cstbm(i, j) = zero
C--------
1015        continue
99987     continue
          do 99986 i = 1, 3
            do g_i_ = 1, g_p_
              g_cstmb(g_i_, i, j) = 0.0d0
            enddo
            cstmb(i, j) = zero
C--------
1016        continue
99986     continue
1012      continue
99985   continue
C
        do 99984 i = 1, 36
          do g_i_ = 1, g_p_
            g_cstcoef(g_i_, i) = 0.0d0
          enddo
          cstcoef(i) = zero
C--------
1017      continue
99984   continue
C
        nlayer = 0
C
        do 99981 j = 1, maxlayer
          do 99983 i = 1, 5
            idcmp23(i, j) = 0
1019        continue
99983     continue
          do 99982 i = 1, 8
            do g_i_ = 1, g_p_
              g_mtcmp23(g_i_, i, j) = 0.0d0
            enddo
            mtcmp23(i, j) = zero
C--------
1020        continue
99982     continue
1018      continue
99981   continue
C
        do 99980 i = 1, 18
          do g_i_ = 1, g_p_
            g_localu(g_i_, i) = 0.0d0
          enddo
          localu(i) = zero
C--------
1021      continue
99980   continue
C
        do 99979 i = 1, 3
          do g_i_ = 1, g_p_
            g_elecrv(g_i_, i) = 0.0d0
          enddo
          elecrv(i) = zero
C--------
1022      continue
99979   continue
C
        do 99978 i = 1, 3
          do g_i_ = 1, g_p_
            g_elestr(g_i_, i) = 0.0d0
          enddo
          elestr(i) = zero
C--------
1023      continue
99978   continue
C
        do 99977 i = 1, 3
          do g_i_ = 1, g_p_
            g_elem(g_i_, i) = 0.0d0
          enddo
          elem(i) = zero
C--------
1024      continue
99977   continue
C
        do 99976 i = 1, 3
          do g_i_ = 1, g_p_
            g_elen(g_i_, i) = 0.0d0
          enddo
          elen(i) = zero
C--------
1025      continue
99976   continue
C
        do 99975 i = 1, 3
          do g_i_ = 1, g_p_
            g_ups(g_i_, i) = 0.0d0
          enddo
          ups(i) = zero
C--------
1026      continue
99975   continue
C
        do 99974 i = 1, 3
          do g_i_ = 1, g_p_
            g_lws(g_i_, i) = 0.0d0
          enddo
          lws(i) = zero
C--------
1027      continue
99974   continue
C
        do 99972 j = 1, 3
          do 99973 i = 1, 3
            aframe(i, j) = zero
1029        continue
99973     continue
1028      continue
99972   continue
C
        do 99970 j = 1, 3
          do 99971 i = 1, 3
            do g_i_ = 1, g_p_
              g_eframe(g_i_, i, j) = 0.0d0
            enddo
            eframe(i, j) = zero
C--------
1031        continue
99971     continue
1030      continue
99970   continue
C
C.....CHECK THE TYPE OF CONSTITUTIVE LAW
C
        if ((ctyp .ne. 0) .and. (ctyp .ne. 1) .and. (ctyp .ne. 2) .and. 
     *(ctyp .ne. 3)) then
          goto 100
        endif
C
C.....CHECK THE ADDRESSING IN STORAGE [CMPCO]
C
        if (ctyp .eq. 1) then
          if ((catt .lt. 1) .or. (catt .gt. nttco)) then
            goto 200
          endif
        endif
C
C.....CHECK THE ADDRESSING IN ARRAYS [IDLAY] AND [MTLAY]
C
        if ((ctyp .eq. 2) .or. (ctyp .eq. 3)) then
          if ((catt .lt. 1) .or. (catt .gt. nttly)) then
            goto 300
          endif
        endif
C
C.....CHECK THE ADDRESSING IN ARRAY [CMPFR]
C
        if ((ctyp .eq. 1) .or. (ctyp .eq. 2) .or. (ctyp .eq. 3)) then
          if ((cfrm .lt. 0) .or. (cfrm .gt. ncmpfr)) then
            goto 400
          endif
        endif
C
C.....INITIALIZE THE CONSTITUTIVE COEFFICIENTS IN CASE OF A TYPE-1 LAW
C
        if (ctyp .eq. 1) then
          do 99969 i = 1, 36
            do g_i_ = 1, g_p_
              g_cstcoef(g_i_, i) = g_cmpco(g_i_, i, catt)
            enddo
            cstcoef(i) = cmpco(i, catt)
C--------
2001        continue
99969     continue
        endif
C
C.....INITIALIZE THE LAYER PROPERTIES FOR TYPE-2 AND TYPE-3 LAWS
C
        if ((ctyp .eq. 2) .or. (ctyp .eq. 3)) then
          nlayer = idlay(2, catt)
          if (nlayer .gt. maxlayer) then
            goto 500
          endif
          do 99968 i = catt, (catt + nlayer - 1)
            ilayer = idlay(3, i)
            if (ilayer .gt. nlayer) then
              goto 600
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
99968     continue
        endif
C
C.....CHECK CONSISTENCY OF BASIC PARAMETERS FOR TYPES 2 AND 3
C
        if ((ctyp .eq. 2) .or. (ctyp .eq. 3)) then
          do 99967 i = 1, nlayer
            if (idcmp23(1, i) .ne. idcmp23(1, 1)) then
              goto 700
            endif
            if (idcmp23(2, i) .ne. nlayer) then
              goto 700
            endif
2003        continue
99967     continue
          do 99966 i = (nlayer + 1), maxlayer
            if (idcmp23(1, i) .ne. 0) then
              goto 700
            endif
            if (idcmp23(2, i) .ne. 0) then
              goto 700
            endif
2004        continue
99966     continue
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
C.....INITIALIZE THE THICKNESS FOR A TYPE-0 CONSTITUTIVE LAW
C.....IT IS ASSUMED CONSTANT HERE
C
        if (ctyp .eq. 0) then
          do g_i_ = 1, g_p_
            g_thick(g_i_) = g_h(g_i_, 1)
          enddo
          thick = h(1)
C--------
        endif
C
C.....INITIALIZE THE THICKNESS FOR A TYPE-1 CONSTITUTIVE LAW
C.....TAKE THE DEFAULT THICKNESS ASSUMED CONSTANT
C.....IF ZERO, ESTIMATE THE CONSTANT THICKNESS USING THE
C.....FIRST COEFFICIENTS OF EXTENTIONAL AND BENDING STIFFNESSES
C
        if (ctyp .eq. 1) then
C
          do g_i_ = 1, g_p_
            g_thick(g_i_) = g_h(g_i_, 1)
          enddo
          thick = h(1)
C--------
C
          if (thick .eq. zero) then
            if (cstcoef(1) .eq. zero) then
              goto 800
            endif
            d4_v = three * cstcoef(22) / cstcoef(1)
            d3_b = -d4_v / cstcoef(1)
            d4_b = 1.0d0 / cstcoef(1) * three
            do g_i_ = 1, g_p_
              g_appxh2(g_i_) = d3_b * g_cstcoef(g_i_, 1) + d4_b * g_cstc
     *oef(g_i_, 22)
            enddo
            appxh2 = d4_v
C--------
            if (appxh2 .le. zero) then
              goto 900
            endif
            d2_v = sqrt(appxh2)
            if ( appxh2 .gt. 0.0d0 ) then
               d1_p = 1.0d0 / (2.0d0 *  d2_v)
            else
               call ehufDV (9,appxh2, d2_v, d1_p,
     +'g_compvms.f',
     +670)
            endif
            do g_i_ = 1, g_p_
              g_thick(g_i_) = d1_p * g_appxh2(g_i_)
            enddo
            thick = d2_v
C--------
          endif
C
        endif
C
C.....INITIALIZE THE THICKNESS FOR TYPE-2 AND TYPE-3 CONSTITUTIVE LAWS
C.....IT IS ASSUMED CONSTANT AND EQUAL TO THE SUM OF EACH LAYER'S THICKNESS
C
        if ((ctyp .eq. 2) .or. (ctyp .eq. 3)) then
C
          do g_i_ = 1, g_p_
            g_thick(g_i_) = 0.0d0
          enddo
          thick = zero
C--------
C
          do 99965 ilayer = 1, nlayer
            do g_i_ = 1, g_p_
              g_thick(g_i_) = g_mtcmp23(g_i_, 7, ilayer) + g_thick(g_i_)
            enddo
            thick = thick + mtcmp23(7, ilayer)
C--------
2005        continue
99965     continue
C
        endif
C
C     ----------------------------------
C     STEP 1
C     COMPUTE THE TRIANGULAR COORDINATES
C     ----------------------------------
C
C.....GET THE ELEMENT TRIANGULAR COORDINATES
C.....GET THE ROTATION MATRIX
C.....GET THE DEGREE OF FREEDOM POINTERS
C
        call g_compcrd2(g_p_, elm, ctyp, globalx, g_globalx, ldg_globalx
     *, globaly, g_globaly, ldg_globaly, globalz, g_globalz, ldg_globalz
     *, rot, g_rot, g_pmax_, x, g_x, g_pmax_, y, g_y, g_pmax_, z, rowb, 
     *colb, rowm, colm, xp, g_xp, g_pmax_, yp, g_yp, g_pmax_, zp, g_zp, 
     *g_pmax_)
C
C.....GET THE ELEMENT LEVEL FRAME
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
C     ---------------------------------------------------
C     STEP 2
C     COMPUTE THE INTEGRATED CURVATURE-NODAL DISPLACEMENT
C     MATRIX FOR PURE BENDING (BASIC STIFFNESS MATRIX)
C     ---------------------------------------------------
C
C.....CHECK IF [CLR] AND [CQR] SATISFY THE CONSTRAINT [CLR]+[CQR]=1
C
        if ((clr + cqr) .ne. one) then
          goto 910
        endif
C
C.....GET THE DISTANCES BETWEEN NODAL POINT X- AND Y- COORDINATES
C
        do g_i_ = 1, g_p_
          g_x21(g_i_) = -g_x(g_i_, 1) + g_x(g_i_, 2)
        enddo
        x21 = x(2) - x(1)
C--------
        do g_i_ = 1, g_p_
          g_x12(g_i_) = -g_x21(g_i_)
        enddo
        x12 = -x21
C--------
        do g_i_ = 1, g_p_
          g_x32(g_i_) = -g_x(g_i_, 2) + g_x(g_i_, 3)
        enddo
        x32 = x(3) - x(2)
C--------
        do g_i_ = 1, g_p_
          g_x23(g_i_) = -g_x32(g_i_)
        enddo
        x23 = -x32
C--------
        do g_i_ = 1, g_p_
          g_x13(g_i_) = -g_x(g_i_, 3) + g_x(g_i_, 1)
        enddo
        x13 = x(1) - x(3)
C--------
        do g_i_ = 1, g_p_
          g_x31(g_i_) = -g_x13(g_i_)
        enddo
        x31 = -x13
C--------
        do g_i_ = 1, g_p_
          g_y21(g_i_) = -g_y(g_i_, 1) + g_y(g_i_, 2)
        enddo
        y21 = y(2) - y(1)
C--------
        do g_i_ = 1, g_p_
          g_y12(g_i_) = -g_y21(g_i_)
        enddo
        y12 = -y21
C--------
        do g_i_ = 1, g_p_
          g_y32(g_i_) = -g_y(g_i_, 2) + g_y(g_i_, 3)
        enddo
        y32 = y(3) - y(2)
C--------
        do g_i_ = 1, g_p_
          g_y23(g_i_) = -g_y32(g_i_)
        enddo
        y23 = -y32
C--------
        do g_i_ = 1, g_p_
          g_y13(g_i_) = -g_y(g_i_, 3) + g_y(g_i_, 1)
        enddo
        y13 = y(1) - y(3)
C--------
        do g_i_ = 1, g_p_
          g_y31(g_i_) = -g_y13(g_i_)
        enddo
        y31 = -y13
C--------
C
C.....CALCULATE TWICE THE AREA OF THE TRIANGLE
C
        do g_i_ = 1, g_p_
          g_twicearea(g_i_) = -x21 * g_y13(g_i_) + (-y13) * g_x21(g_i_) 
     *+ y21 * g_x13(g_i_) + x13 * g_y21(g_i_)
        enddo
        twicearea = y21 * x13 - x21 * y13
C--------
C
C.....CALCULATE THE AREA OF THE TRIANGLE
C
        d2_b = 1.0d0 / two
        do g_i_ = 1, g_p_
          g_area(g_i_) = d2_b * g_twicearea(g_i_)
        enddo
        area = twicearea / two
C--------
C
C.....CHECK THE AREA (ERROR IF NOT POSITIVE)
C
        if (twicearea .le. zero) then
          goto 920
        endif
C
C.....GET THE COORDINATES OF THE CENTROID OF THE TRIANGLE
C
        x0 = (x(1) + x(2) + x(3)) / three
        y0 = (y(1) + y(2) + y(3)) / three
C
C.....GET THE DISTANCES BETWEEN NODES 1-2, 2-3 AND 3-1
C
        d4_b = y12 + y12
        d5_b = x12 + x12
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y12(g_i_) + d5_b * g_x12(g_i_)
        enddo
        d1_w = x12 * x12 + y12 * y12
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +884)
        endif
        do g_i_ = 1, g_p_
          g_dist12(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        dist12 = d2_v
C--------
        d4_b = y23 + y23
        d5_b = x23 + x23
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y23(g_i_) + d5_b * g_x23(g_i_)
        enddo
        d1_w = x23 * x23 + y23 * y23
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +903)
        endif
        do g_i_ = 1, g_p_
          g_dist23(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        dist23 = d2_v
C--------
        d4_b = y31 + y31
        d5_b = x31 + x31
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d4_b * g_y31(g_i_) + d5_b * g_x31(g_i_)
        enddo
        d1_w = x31 * x31 + y31 * y31
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +922)
        endif
        do g_i_ = 1, g_p_
          g_dist31(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        dist31 = d2_v
C--------
C
C.....ASSEMBLE THE LOCAL MATRIX [LLR] W/ SHAPE FUNCTION DERIVATIVES
C
        if (clr .ne. zero) then
C
          d2_b = 1.0d0 / two
          do g_i_ = 1, g_p_
            g_llr(g_i_, 3, 1) = d2_b * g_y32(g_i_)
          enddo
          llr(3, 1) = y32 / two
C--------
          d2_b = 1.0d0 / two
          do g_i_ = 1, g_p_
            g_llr(g_i_, 6, 1) = d2_b * g_y13(g_i_)
          enddo
          llr(6, 1) = y13 / two
C--------
          d2_b = 1.0d0 / two
          do g_i_ = 1, g_p_
            g_llr(g_i_, 9, 1) = d2_b * g_y21(g_i_)
          enddo
          llr(9, 1) = y21 / two
C--------
C
          d2_b = 1.0d0 / two
          do g_i_ = 1, g_p_
            g_llr(g_i_, 2, 2) = d2_b * g_x32(g_i_)
          enddo
          llr(2, 2) = x32 / two
C--------
          d2_b = 1.0d0 / two
          do g_i_ = 1, g_p_
            g_llr(g_i_, 5, 2) = d2_b * g_x13(g_i_)
          enddo
          llr(5, 2) = x13 / two
C--------
          d2_b = 1.0d0 / two
          do g_i_ = 1, g_p_
            g_llr(g_i_, 8, 2) = d2_b * g_x21(g_i_)
          enddo
          llr(8, 2) = x21 / two
C--------
C
          d3_b = -(1.0d0 / two)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 2, 3) = d3_b * g_y32(g_i_)
          enddo
          llr(2, 3) = -y32 / two
C--------
          d3_b = -(1.0d0 / two)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 3, 3) = d3_b * g_x32(g_i_)
          enddo
          llr(3, 3) = -x32 / two
C--------
          d3_b = -(1.0d0 / two)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 5, 3) = d3_b * g_y13(g_i_)
          enddo
          llr(5, 3) = -y13 / two
C--------
          d3_b = -(1.0d0 / two)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 6, 3) = d3_b * g_x13(g_i_)
          enddo
          llr(6, 3) = -x13 / two
C--------
          d3_b = -(1.0d0 / two)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 8, 3) = d3_b * g_y21(g_i_)
          enddo
          llr(8, 3) = -y21 / two
C--------
          d3_b = -(1.0d0 / two)
          do g_i_ = 1, g_p_
            g_llr(g_i_, 9, 3) = d3_b * g_x21(g_i_)
          enddo
          llr(9, 3) = -x21 / two
C--------
C
        endif
C
C.....ASSEMBLE THE LOCAL MATRIX [LQR] W/ SHAPE FUNCTION DERIVATIVES
C
        if (cqr .ne. zero) then
C
          d3_v = y21 / dist12
          d2_b = 1.0d0 / dist12
          d3_b = -d3_v / dist12
          do g_i_ = 1, g_p_
            g_c12(g_i_) = d3_b * g_dist12(g_i_) + d2_b * g_y21(g_i_)
          enddo
          c12 = d3_v
C--------
          d3_v = x12 / dist12
          d2_b = 1.0d0 / dist12
          d3_b = -d3_v / dist12
          do g_i_ = 1, g_p_
            g_s12(g_i_) = d3_b * g_dist12(g_i_) + d2_b * g_x12(g_i_)
          enddo
          s12 = d3_v
C--------
          d3_v = y32 / dist23
          d2_b = 1.0d0 / dist23
          d3_b = -d3_v / dist23
          do g_i_ = 1, g_p_
            g_c23(g_i_) = d3_b * g_dist23(g_i_) + d2_b * g_y32(g_i_)
          enddo
          c23 = d3_v
C--------
          d3_v = x23 / dist23
          d2_b = 1.0d0 / dist23
          d3_b = -d3_v / dist23
          do g_i_ = 1, g_p_
            g_s23(g_i_) = d3_b * g_dist23(g_i_) + d2_b * g_x23(g_i_)
          enddo
          s23 = d3_v
C--------
          d3_v = y13 / dist31
          d2_b = 1.0d0 / dist31
          d3_b = -d3_v / dist31
          do g_i_ = 1, g_p_
            g_c31(g_i_) = d3_b * g_dist31(g_i_) + d2_b * g_y13(g_i_)
          enddo
          c31 = d3_v
C--------
          d3_v = x31 / dist31
          d2_b = 1.0d0 / dist31
          d3_b = -d3_v / dist31
          do g_i_ = 1, g_p_
            g_s31(g_i_) = d3_b * g_dist31(g_i_) + d2_b * g_x31(g_i_)
          enddo
          s31 = d3_v
C--------
C
          d2_b = c12 + c12
          do g_i_ = 1, g_p_
            g_cc12(g_i_) = d2_b * g_c12(g_i_)
          enddo
          cc12 = c12 * c12
C--------
          d2_b = c23 + c23
          do g_i_ = 1, g_p_
            g_cc23(g_i_) = d2_b * g_c23(g_i_)
          enddo
          cc23 = c23 * c23
C--------
          d2_b = c31 + c31
          do g_i_ = 1, g_p_
            g_cc31(g_i_) = d2_b * g_c31(g_i_)
          enddo
          cc31 = c31 * c31
C--------
          d2_b = s12 + s12
          do g_i_ = 1, g_p_
            g_ss12(g_i_) = d2_b * g_s12(g_i_)
          enddo
          ss12 = s12 * s12
C--------
          d2_b = s23 + s23
          do g_i_ = 1, g_p_
            g_ss23(g_i_) = d2_b * g_s23(g_i_)
          enddo
          ss23 = s23 * s23
C--------
          d2_b = s31 + s31
          do g_i_ = 1, g_p_
            g_ss31(g_i_) = d2_b * g_s31(g_i_)
          enddo
          ss31 = s31 * s31
C--------
          do g_i_ = 1, g_p_
            g_cs12(g_i_) = c12 * g_s12(g_i_) + s12 * g_c12(g_i_)
          enddo
          cs12 = c12 * s12
C--------
          do g_i_ = 1, g_p_
            g_cs23(g_i_) = c23 * g_s23(g_i_) + s23 * g_c23(g_i_)
          enddo
          cs23 = c23 * s23
C--------
          do g_i_ = 1, g_p_
            g_cs31(g_i_) = c31 * g_s31(g_i_) + s31 * g_c31(g_i_)
          enddo
          cs31 = c31 * s31
C--------
C
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 1, 1) = -g_cs31(g_i_) + g_cs12(g_i_)
          enddo
          lqr(1, 1) = cs12 - cs31
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 1, 2) = -g_lqr(g_i_, 1, 1)
          enddo
          lqr(1, 2) = -lqr(1, 1)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 1, 3) = g_ss12(g_i_) + (-g_cc12(g_i_)) + (-g_ss3
     *1(g_i_)) + g_cc31(g_i_)
          enddo
          lqr(1, 3) = cc31 - ss31 - (cc12 - ss12)
C--------
C
          d2_b = 1.0d0 / two
          d5_b = d2_b * x31
          d6_b = d2_b * cc31
          d7_b = d2_b * x12
          d8_b = d2_b * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 2, 1) = d6_b * g_x31(g_i_) + d5_b * g_cc31(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_cc12(g_i_)
          enddo
          lqr(2, 1) = (cc12 * x12 + cc31 * x31) / two
C--------
          d2_b = 1.0d0 / two
          d5_b = d2_b * x31
          d6_b = d2_b * ss31
          d7_b = d2_b * x12
          d8_b = d2_b * ss12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 2, 2) = d6_b * g_x31(g_i_) + d5_b * g_ss31(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_ss12(g_i_)
          enddo
          lqr(2, 2) = (ss12 * x12 + ss31 * x31) / two
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 2, 3) = ss31 * g_y13(g_i_) + y13 * g_ss31(g_i_) 
     *+ ss12 * g_y21(g_i_) + y21 * g_ss12(g_i_)
          enddo
          lqr(2, 3) = ss12 * y21 + ss31 * y13
C--------
C
          d3_b = -(1.0d0 / two)
          d6_b = d3_b * y13
          d7_b = d3_b * cc31
          d8_b = d3_b * y21
          d9_b = d3_b * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 3, 1) = d7_b * g_y13(g_i_) + d6_b * g_cc31(g_i_)
     * + d9_b * g_y21(g_i_) + d8_b * g_cc12(g_i_)
          enddo
          lqr(3, 1) = -(cc12 * y21 + cc31 * y13) / two
C--------
          d3_b = -(1.0d0 / two)
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 3, 2) = d3_b * g_lqr(g_i_, 2, 3)
          enddo
          lqr(3, 2) = -lqr(2, 3) / two
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 3, 3) = -two * g_lqr(g_i_, 2, 1)
          enddo
          lqr(3, 3) = -two * lqr(2, 1)
C--------
C
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 4, 1) = -g_cs12(g_i_) + g_cs23(g_i_)
          enddo
          lqr(4, 1) = cs23 - cs12
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 4, 2) = -g_lqr(g_i_, 4, 1)
          enddo
          lqr(4, 2) = -lqr(4, 1)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 4, 3) = g_ss23(g_i_) + (-g_cc23(g_i_)) + (-g_ss1
     *2(g_i_)) + g_cc12(g_i_)
          enddo
          lqr(4, 3) = cc12 - ss12 - (cc23 - ss23)
C--------
C
          d2_b = 1.0d0 / two
          d5_b = d2_b * x23
          d6_b = d2_b * cc23
          d7_b = d2_b * x12
          d8_b = d2_b * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 5, 1) = d6_b * g_x23(g_i_) + d5_b * g_cc23(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_cc12(g_i_)
          enddo
          lqr(5, 1) = (cc12 * x12 + cc23 * x23) / two
C--------
          d2_b = 1.0d0 / two
          d5_b = d2_b * x23
          d6_b = d2_b * ss23
          d7_b = d2_b * x12
          d8_b = d2_b * ss12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 5, 2) = d6_b * g_x23(g_i_) + d5_b * g_ss23(g_i_)
     * + d8_b * g_x12(g_i_) + d7_b * g_ss12(g_i_)
          enddo
          lqr(5, 2) = (ss12 * x12 + ss23 * x23) / two
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 5, 3) = ss23 * g_y32(g_i_) + y32 * g_ss23(g_i_) 
     *+ ss12 * g_y21(g_i_) + y21 * g_ss12(g_i_)
          enddo
          lqr(5, 3) = ss12 * y21 + ss23 * y32
C--------
C
          d3_b = -(1.0d0 / two)
          d6_b = d3_b * y32
          d7_b = d3_b * cc23
          d8_b = d3_b * y21
          d9_b = d3_b * cc12
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 6, 1) = d7_b * g_y32(g_i_) + d6_b * g_cc23(g_i_)
     * + d9_b * g_y21(g_i_) + d8_b * g_cc12(g_i_)
          enddo
          lqr(6, 1) = -(cc12 * y21 + cc23 * y32) / two
C--------
          d3_b = -(1.0d0 / two)
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 6, 2) = d3_b * g_lqr(g_i_, 5, 3)
          enddo
          lqr(6, 2) = -lqr(5, 3) / two
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 6, 3) = -two * g_lqr(g_i_, 5, 1)
          enddo
          lqr(6, 3) = -two * lqr(5, 1)
C--------
C
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 7, 1) = -g_cs23(g_i_) + g_cs31(g_i_)
          enddo
          lqr(7, 1) = cs31 - cs23
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 7, 2) = -g_lqr(g_i_, 7, 1)
          enddo
          lqr(7, 2) = -lqr(7, 1)
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 7, 3) = g_ss31(g_i_) + (-g_cc31(g_i_)) + (-g_ss2
     *3(g_i_)) + g_cc23(g_i_)
          enddo
          lqr(7, 3) = cc23 - ss23 - (cc31 - ss31)
C--------
C
          d2_b = 1.0d0 / two
          d5_b = d2_b * x31
          d6_b = d2_b * cc31
          d7_b = d2_b * x23
          d8_b = d2_b * cc23
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 8, 1) = d6_b * g_x31(g_i_) + d5_b * g_cc31(g_i_)
     * + d8_b * g_x23(g_i_) + d7_b * g_cc23(g_i_)
          enddo
          lqr(8, 1) = (cc23 * x23 + cc31 * x31) / two
C--------
          d2_b = 1.0d0 / two
          d5_b = d2_b * x31
          d6_b = d2_b * ss31
          d7_b = d2_b * x23
          d8_b = d2_b * ss23
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 8, 2) = d6_b * g_x31(g_i_) + d5_b * g_ss31(g_i_)
     * + d8_b * g_x23(g_i_) + d7_b * g_ss23(g_i_)
          enddo
          lqr(8, 2) = (ss23 * x23 + ss31 * x31) / two
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 8, 3) = ss31 * g_y13(g_i_) + y13 * g_ss31(g_i_) 
     *+ ss23 * g_y32(g_i_) + y32 * g_ss23(g_i_)
          enddo
          lqr(8, 3) = ss23 * y32 + ss31 * y13
C--------
C
          d3_b = -(1.0d0 / two)
          d6_b = d3_b * y13
          d7_b = d3_b * cc31
          d8_b = d3_b * y32
          d9_b = d3_b * cc23
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 9, 1) = d7_b * g_y13(g_i_) + d6_b * g_cc31(g_i_)
     * + d9_b * g_y32(g_i_) + d8_b * g_cc23(g_i_)
          enddo
          lqr(9, 1) = -(cc23 * y32 + cc31 * y13) / two
C--------
          d3_b = -(1.0d0 / two)
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 9, 2) = d3_b * g_lqr(g_i_, 8, 3)
          enddo
          lqr(9, 2) = -lqr(8, 3) / two
C--------
          do g_i_ = 1, g_p_
            g_lqr(g_i_, 9, 3) = -two * g_lqr(g_i_, 8, 1)
          enddo
          lqr(9, 3) = -two * lqr(8, 1)
C--------
C
        endif
C
C.....ASSEMBLE THE LOCAL MATRIX [L] AS [CLR]*[LLR] + [CQR]*[LQR]
C
        if (clr .eq. zero) then
          do 99963 j = 1, 3
            do 99964 i = 1, 9
              do g_i_ = 1, g_p_
                g_reducedl(g_i_, i, j) = g_lqr(g_i_, i, j)
              enddo
              reducedl(i, j) = lqr(i, j)
C--------
3002          continue
99964       continue
3001        continue
99963     continue
        endif
C
        if (cqr .eq. zero) then
          do 99961 j = 1, 3
            do 99962 i = 1, 9
              do g_i_ = 1, g_p_
                g_reducedl(g_i_, i, j) = g_llr(g_i_, i, j)
              enddo
              reducedl(i, j) = llr(i, j)
C--------
3004          continue
99962       continue
3003        continue
99961     continue
        endif
C
        if ((clr .ne. zero) .and. (cqr .ne. zero)) then
          do 99959 j = 1, 3
            do 99960 i = 1, 9
              do g_i_ = 1, g_p_
                g_reducedl(g_i_, i, j) = cqr * g_lqr(g_i_, i, j) + clr *
     * g_llr(g_i_, i, j)
              enddo
              reducedl(i, j) = clr * llr(i, j) + cqr * lqr(i, j)
C--------
3006          continue
99960       continue
3005        continue
99959     continue
        endif
C
C.....DIVIDE MATRIX [L] BY THE SQUARE ROOT OF THE AREA
C
        d2_v = one / area
        d2_b = -d2_v / area
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d2_b * g_area(g_i_)
        enddo
        d1_w = d2_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +1384)
        endif
        do g_i_ = 1, g_p_
          g_factor(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        factor = d2_v
C--------
C
        do 99957 j = 1, 3
          do 99958 i = 1, 9
            do g_i_ = 1, g_p_
              g_reducedl(g_i_, i, j) = factor * g_reducedl(g_i_, i, j) +
     * reducedl(i, j) * g_factor(g_i_)
            enddo
            reducedl(i, j) = factor * reducedl(i, j)
C--------
3008        continue
99958     continue
3007      continue
99957   continue
C
C.....TRANSFORM MATRIX [L] IN FEM-LIKE DOF NUMBERING
C
        do 99956 i = 1, 9
          row = rowb(i)
          do g_i_ = 1, g_p_
            g_l(g_i_, row, 1) = g_reducedl(g_i_, i, 1)
          enddo
          l(row, 1) = reducedl(i, 1)
C--------
          do g_i_ = 1, g_p_
            g_l(g_i_, row, 2) = g_reducedl(g_i_, i, 2)
          enddo
          l(row, 2) = reducedl(i, 2)
C--------
          do g_i_ = 1, g_p_
            g_l(g_i_, row, 3) = g_reducedl(g_i_, i, 3)
          enddo
          l(row, 3) = reducedl(i, 3)
C--------
3009      continue
99956   continue
C
C     -------------------------------------------------
C     STEP 3
C     COMPUTE THE INTEGRATED STRAIN-NODAL DISPLACEMENT
C     MATRIX FOR PURE MEMBRANE (BASIC STIFFNESS MATRIX)
C     -------------------------------------------------
C
C.....CHECK IF FACTOR [ALPHA] IS POSITIVE
C
        if (alpha .le. zero) then
          goto 930
        endif
C
C.....ASSEMBLE THE MATRIX [P] W/ SHAPE FUNCTION DERIVATIVES
C
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 1, 1) = g_y23(g_i_)
        enddo
        reducedp(1, 1) = y23
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 2, 1) = 0.0d0
        enddo
        reducedp(2, 1) = zero
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 3, 1) = g_y31(g_i_)
        enddo
        reducedp(3, 1) = y31
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 4, 1) = 0.0d0
        enddo
        reducedp(4, 1) = zero
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 5, 1) = g_y12(g_i_)
        enddo
        reducedp(5, 1) = y12
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 6, 1) = 0.0d0
        enddo
        reducedp(6, 1) = zero
C--------
C
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 1, 2) = 0.0d0
        enddo
        reducedp(1, 2) = zero
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 2, 2) = g_x32(g_i_)
        enddo
        reducedp(2, 2) = x32
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 3, 2) = 0.0d0
        enddo
        reducedp(3, 2) = zero
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 4, 2) = g_x13(g_i_)
        enddo
        reducedp(4, 2) = x13
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 5, 2) = 0.0d0
        enddo
        reducedp(5, 2) = zero
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 6, 2) = g_x21(g_i_)
        enddo
        reducedp(6, 2) = x21
C--------
C
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 1, 3) = g_x32(g_i_)
        enddo
        reducedp(1, 3) = x32
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 2, 3) = g_y23(g_i_)
        enddo
        reducedp(2, 3) = y23
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 3, 3) = g_x13(g_i_)
        enddo
        reducedp(3, 3) = x13
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 4, 3) = g_y31(g_i_)
        enddo
        reducedp(4, 3) = y31
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 5, 3) = g_x21(g_i_)
        enddo
        reducedp(5, 3) = x21
C--------
        do g_i_ = 1, g_p_
          g_reducedp(g_i_, 6, 3) = g_y12(g_i_)
        enddo
        reducedp(6, 3) = y12
C--------
C
        dimp = 6
C
        if (alpha .ne. zero) then
C
          d4_v = y13 - y21
          d3_b = 1.0d0 / six * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * y23
          do g_i_ = 1, g_p_
            g_reducedp(g_i_, 7, 1) = -d5_b * g_y21(g_i_) + d5_b * g_y13(
     *g_i_) + d4_b * g_y23(g_i_)
          enddo
          reducedp(7, 1) = y23 * d4_v * alpha / six
C--------
          d4_v = x31 - x12
          d3_b = 1.0d0 / six * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x32
          do g_i_ = 1, g_p_
            g_reducedp(g_i_, 7, 2) = -d5_b * g_x12(g_i_) + d5_b * g_x31(
     *g_i_) + d4_b * g_x32(g_i_)
          enddo
          reducedp(7, 2) = x32 * d4_v * alpha / six
C--------
          d3_b = 1.0d0 / three * alpha
          d6_b = -d3_b * y21
          d7_b = -d3_b * x12
          d8_b = d3_b * y13
          d9_b = d3_b * x31
          do g_i_ = 1, g_p_
            g_reducedp(g_i_, 7, 3) = d7_b * g_y21(g_i_) + d6_b * g_x12(g
     *_i_) + d9_b * g_y13(g_i_) + d8_b * g_x31(g_i_)
          enddo
          reducedp(7, 3) = (x31 * y13 - x12 * y21) * alpha / three
C--------
C
          d4_v = y21 - y32
          d3_b = 1.0d0 / six * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * y31
          do g_i_ = 1, g_p_
            g_reducedp(g_i_, 8, 1) = -d5_b * g_y32(g_i_) + d5_b * g_y21(
     *g_i_) + d4_b * g_y31(g_i_)
          enddo
          reducedp(8, 1) = y31 * d4_v * alpha / six
C--------
          d4_v = x12 - x23
          d3_b = 1.0d0 / six * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x13
          do g_i_ = 1, g_p_
            g_reducedp(g_i_, 8, 2) = -d5_b * g_x23(g_i_) + d5_b * g_x12(
     *g_i_) + d4_b * g_x13(g_i_)
          enddo
          reducedp(8, 2) = x13 * d4_v * alpha / six
C--------
          d3_b = 1.0d0 / three * alpha
          d6_b = -d3_b * y32
          d7_b = -d3_b * x23
          d8_b = d3_b * y21
          d9_b = d3_b * x12
          do g_i_ = 1, g_p_
            g_reducedp(g_i_, 8, 3) = d7_b * g_y32(g_i_) + d6_b * g_x23(g
     *_i_) + d9_b * g_y21(g_i_) + d8_b * g_x12(g_i_)
          enddo
          reducedp(8, 3) = (x12 * y21 - x23 * y32) * alpha / three
C--------
C
          d4_v = y32 - y13
          d3_b = 1.0d0 / six * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * y12
          do g_i_ = 1, g_p_
            g_reducedp(g_i_, 9, 1) = -d5_b * g_y13(g_i_) + d5_b * g_y32(
     *g_i_) + d4_b * g_y12(g_i_)
          enddo
          reducedp(9, 1) = y12 * d4_v * alpha / six
C--------
          d4_v = x23 - x31
          d3_b = 1.0d0 / six * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x21
          do g_i_ = 1, g_p_
            g_reducedp(g_i_, 9, 2) = -d5_b * g_x31(g_i_) + d5_b * g_x23(
     *g_i_) + d4_b * g_x21(g_i_)
          enddo
          reducedp(9, 2) = x21 * d4_v * alpha / six
C--------
          d3_b = 1.0d0 / three * alpha
          d6_b = -d3_b * y13
          d7_b = -d3_b * x31
          d8_b = d3_b * y32
          d9_b = d3_b * x23
          do g_i_ = 1, g_p_
            g_reducedp(g_i_, 9, 3) = d7_b * g_y13(g_i_) + d6_b * g_x31(g
     *_i_) + d9_b * g_y32(g_i_) + d8_b * g_x23(g_i_)
          enddo
          reducedp(9, 3) = (x23 * y32 - x31 * y13) * alpha / three
C--------
C
          dimp = 9
C
        endif
C
C.....DIVIDE MATRIX [P] BY THE SQUARE ROOT OF FOUR TIMES THE AREA
C
        d2_v = two * twicearea
        d3_v = one / d2_v
        d3_b = -d3_v / d2_v * two
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_twicearea(g_i_)
        enddo
        d1_w = d3_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +1653)
        endif
        do g_i_ = 1, g_p_
          g_factor(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        factor = d2_v
C--------
C
        do 99954 j = 1, 3
          do 99955 i = 1, dimp
            do g_i_ = 1, g_p_
              g_reducedp(g_i_, i, j) = factor * g_reducedp(g_i_, i, j) +
     * reducedp(i, j) * g_factor(g_i_)
            enddo
            reducedp(i, j) = factor * reducedp(i, j)
C--------
4002        continue
99955     continue
4001      continue
99954   continue
C
C.....TRANSFORM MATRIX [L] IN FEM-LIKE DOF NUMBERING
C
        do 99953 i = 1, 9
          row = rowm(i)
          do g_i_ = 1, g_p_
            g_p(g_i_, row, 1) = g_reducedp(g_i_, i, 1)
          enddo
          p(row, 1) = reducedp(i, 1)
C--------
          do g_i_ = 1, g_p_
            g_p(g_i_, row, 2) = g_reducedp(g_i_, i, 2)
          enddo
          p(row, 2) = reducedp(i, 2)
C--------
          do g_i_ = 1, g_p_
            g_p(g_i_, row, 3) = g_reducedp(g_i_, i, 3)
          enddo
          p(row, 3) = reducedp(i, 3)
C--------
4003      continue
99953   continue
C
C     -------------------------------------------
C     STEP 4
C     ROTATE THE NODAL DISPLACEMENTS TO THE LOCAL
C     FRAME SYSTEM (LOCAL TO THE SHELL ELEMENT)
C     -------------------------------------------
C
        do 99949 i = 1, 6
          do 99952 j = 1, 6
            do g_i_ = 1, g_p_
              g_localu(g_i_, i) = rot(j, i) * g_globalu(g_i_, j) + globa
     *lu(j) * g_rot(g_i_, j, i) + g_localu(g_i_, i)
            enddo
            localu(i) = localu(i) + rot(j, i) * globalu(j)
C--------
5002        continue
99952     continue
          do 99951 j = 1, 6
            do g_i_ = 1, g_p_
              g_localu(g_i_, i + 6) = rot(j, i) * g_globalu(g_i_, j + 6)
     * + globalu(j + 6) * g_rot(g_i_, j, i) + g_localu(g_i_, i + 6)
            enddo
            localu(i + 6) = localu(i + 6) + rot(j, i) * globalu(j + 6)
C--------
5003        continue
99951     continue
          do 99950 j = 1, 6
            do g_i_ = 1, g_p_
              g_localu(g_i_, i + 12) = rot(j, i) * g_globalu(g_i_, j + 1
     *2) + globalu(j + 12) * g_rot(g_i_, j, i) + g_localu(g_i_, i + 12)
            enddo
            localu(i + 12) = localu(i + 12) + rot(j, i) * globalu(j + 12
     *)
C--------
5004        continue
99950     continue
5001      continue
99949   continue
C
C     --------------------------------------------------
C     STEP 5
C     COMPUTE THE ELEMENTAL STRAIN AND CURVATURE VECTORS
C     --------------------------------------------------
C
C.....ELEMENTAL CURVATURE COMPUTATION (1/radius)
C
        do 99947 i = 1, 3
          do 99948 j = 1, 18
            do g_i_ = 1, g_p_
              g_elecrv(g_i_, i) = l(j, i) * g_localu(g_i_, j) + localu(j
     *) * g_l(g_i_, j, i) + g_elecrv(g_i_, i)
            enddo
            elecrv(i) = elecrv(i) + l(j, i) * localu(j)
C--------
6002        continue
99948     continue
6001      continue
99947   continue
C
C.....ELEMENTAL STRAIN COMPUTATION
C
        do 99945 i = 1, 3
          do 99946 j = 1, 18
            do g_i_ = 1, g_p_
              g_elestr(g_i_, i) = p(j, i) * g_localu(g_i_, j) + localu(j
     *) * g_p(g_i_, j, i) + g_elestr(g_i_, i)
            enddo
            elestr(i) = elestr(i) + p(j, i) * localu(j)
C--------
6004        continue
99946     continue
6003      continue
99945   continue
C
C COMPUTE EQUIVALENT STRAIN (VON MISES)
C
        if (strainflg .eq. 1) then
          d2_b = dble(0.5)
          do g_i_ = 1, g_p_
            g_t2(g_i_) = d2_b * g_thick(g_i_)
          enddo
          t2 = dble(0.5) * thick
C--------
C
          do g_i_ = 1, g_p_
            g_elecrv(g_i_, 1) = t2 * g_elecrv(g_i_, 1) + elecrv(1) * g_t
     *2(g_i_)
          enddo
          elecrv(1) = t2 * elecrv(1)
C--------
          do g_i_ = 1, g_p_
            g_elecrv(g_i_, 2) = t2 * g_elecrv(g_i_, 2) + elecrv(2) * g_t
     *2(g_i_)
          enddo
          elecrv(2) = t2 * elecrv(2)
C--------
          do g_i_ = 1, g_p_
            g_elecrv(g_i_, 3) = t2 * g_elecrv(g_i_, 3) + elecrv(3) * g_t
     *2(g_i_)
          enddo
          elecrv(3) = t2 * elecrv(3)
C--------
C
          if (surface .eq. 2) then
            do g_i_ = 1, g_p_
              g_epsxx(g_i_) = g_elestr(g_i_, 1)
            enddo
            epsxx = elestr(1)
C--------
            do g_i_ = 1, g_p_
              g_epsyy(g_i_) = g_elestr(g_i_, 2)
            enddo
            epsyy = elestr(2)
C--------
            d2_b = dble(0.5)
            do g_i_ = 1, g_p_
              g_epsxy(g_i_) = d2_b * g_elestr(g_i_, 3)
            enddo
            epsxy = dble(0.5) * elestr(3)
C--------
          else
            if (surface .eq. 3) then
              do g_i_ = 1, g_p_
                g_epsxx(g_i_) = -g_elecrv(g_i_, 1) + g_elestr(g_i_, 1)
              enddo
              epsxx = elestr(1) - elecrv(1)
C--------
              do g_i_ = 1, g_p_
                g_epsyy(g_i_) = -g_elecrv(g_i_, 2) + g_elestr(g_i_, 2)
              enddo
              epsyy = elestr(2) - elecrv(2)
C--------
              d2_b = dble(0.5)
              do g_i_ = 1, g_p_
                g_epsxy(g_i_) = -d2_b * g_elecrv(g_i_, 3) + d2_b * g_ele
     *str(g_i_, 3)
              enddo
              epsxy = dble(0.5) * (elestr(3) - elecrv(3))
C--------
            else
              do g_i_ = 1, g_p_
                g_epsxx(g_i_) = g_elecrv(g_i_, 1) + g_elestr(g_i_, 1)
              enddo
              epsxx = elestr(1) + elecrv(1)
C--------
              do g_i_ = 1, g_p_
                g_epsyy(g_i_) = g_elecrv(g_i_, 2) + g_elestr(g_i_, 2)
              enddo
              epsyy = elestr(2) + elecrv(2)
C--------
              d2_b = dble(0.5)
              do g_i_ = 1, g_p_
                g_epsxy(g_i_) = d2_b * g_elecrv(g_i_, 3) + d2_b * g_eles
     *tr(g_i_, 3)
              enddo
              epsxy = dble(0.5) * (elestr(3) + elecrv(3))
C--------
            endif
          endif
C
          d2_b = -nu / (1.0d0 - nu)
          do g_i_ = 1, g_p_
            g_epszz(g_i_) = d2_b * g_epsyy(g_i_) + d2_b * g_epsxx(g_i_)
          enddo
          epszz = -nu / (1.0d0 - nu) * (epsxx + epsyy)
C--------
C
          do g_i_ = 1, g_p_
            g_str(g_i_, 1) = g_epsxx(g_i_)
          enddo
          str(1) = epsxx
C--------
          do g_i_ = 1, g_p_
            g_str(g_i_, 2) = g_epsyy(g_i_)
          enddo
          str(2) = epsyy
C--------
          do g_i_ = 1, g_p_
            g_str(g_i_, 3) = g_epszz(g_i_)
          enddo
          str(3) = epszz
C--------
          do g_i_ = 1, g_p_
            g_str(g_i_, 4) = g_epsxy(g_i_)
          enddo
          str(4) = epsxy
C--------
          do g_i_ = 1, g_p_
            g_str(g_i_, 5) = 0.0d0
          enddo
          str(5) = 0.0d0
C--------
          do g_i_ = 1, g_p_
            g_str(g_i_, 6) = 0.0d0
          enddo
          str(6) = 0.0d0
C--------
C
          call g_transform(g_p_, xp, g_xp, g_pmax_, yp, g_yp, g_pmax_, z
     *p, g_zp, g_pmax_, xg, yg, zg, str, g_str, g_pmax_)
C
          do 99944 i = 1, 6
            do g_i_ = 1, g_p_
              g_stress(g_i_, elm, i, 1) = g_str(g_i_, i)
            enddo
            stress(elm, i, 1) = str(i)
C--------
            do g_i_ = 1, g_p_
              g_stress(g_i_, elm, i, 2) = g_str(g_i_, i)
            enddo
            stress(elm, i, 2) = str(i)
C--------
            do g_i_ = 1, g_p_
              g_stress(g_i_, elm, i, 3) = g_str(g_i_, i)
            enddo
            stress(elm, i, 3) = str(i)
C--------
101         continue
99944     continue
C
          call g_straineq(g_p_, elecrv(1), g_elecrv(1, 1), g_pmax_, elec
     *rv(2), g_elecrv(1, 2), g_pmax_, elecrv(3), g_elecrv(1, 3), g_pmax_
     *, elestr(1), g_elestr(1, 1), g_pmax_, elestr(2), g_elestr(1, 2), g
     *_pmax_, elestr(3), g_elestr(1, 3), g_pmax_, thick, g_thick, g_pmax
     *_, surface, ebar, g_ebar, g_pmax_)
C
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 7, 1) = g_ebar(g_i_)
          enddo
          stress(elm, 7, 1) = ebar
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 7, 2) = g_ebar(g_i_)
          enddo
          stress(elm, 7, 2) = ebar
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 7, 3) = g_ebar(g_i_)
          enddo
          stress(elm, 7, 3) = ebar
C--------
C
          return
        endif
C     ----------------------------------------------
C     STEP 6
C     COMPUTE THE ELEMENTAL MOMENT AND FORCE
C     PER UNIT LENGTH VECTORS AND COMPUTE THE
C     CENTROIDAL VON MISES STRESS RESULTANT
C     INTEGRATED OVER THE THICKNESS OF THE COMPOSITE
C     ----------------------------------------------
C
C.....GET THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
        call g_compcst(g_p_, e, g_e, ldg_e, thick, g_thick, g_pmax_, nu,
     * cstcoef, g_cstcoef, g_pmax_, nlayer, idcmp23, mtcmp23, g_mtcmp23,
     * g_pmax_, x, g_x, g_pmax_, y, g_y, g_pmax_, z, cstbb, g_cstbb, g_p
     *max_, ctyp, eframe, g_eframe, g_pmax_, aframe, 'BB')
C
C.....GET THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
        call g_compcst(g_p_, e, g_e, ldg_e, thick, g_thick, g_pmax_, nu,
     * cstcoef, g_cstcoef, g_pmax_, nlayer, idcmp23, mtcmp23, g_mtcmp23,
     * g_pmax_, x, g_x, g_pmax_, y, g_y, g_pmax_, z, cstmm, g_cstmm, g_p
     *max_, ctyp, eframe, g_eframe, g_pmax_, aframe, 'MM')
C
C.....GET THE CONSTITUTIVE MATRIX FOR BENDING-MEMBRANE COUPLING
C
        call g_compcst(g_p_, e, g_e, ldg_e, thick, g_thick, g_pmax_, nu,
     * cstcoef, g_cstcoef, g_pmax_, nlayer, idcmp23, mtcmp23, g_mtcmp23,
     * g_pmax_, x, g_x, g_pmax_, y, g_y, g_pmax_, z, cstbm, g_cstbm, g_p
     *max_, ctyp, eframe, g_eframe, g_pmax_, aframe, 'BM')
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
     *3, g_pmax_, x, g_x, g_pmax_, y, g_y, g_pmax_, z, cstmb, g_cstmb, g
     *_pmax_, ctyp, eframe, g_eframe, g_pmax_, aframe, 'MB')
C
        endif
C
C.....CHECK THE CONSTITUTIVE MATRIX
C.....(COMMENTED OUT HERE - USED FOR DEBUGGING ONLY)
C
C**   call compchk( elm , cstbb , cstmm , cstbm , cstmb )
C
C.....ESTIMATE THE ELEMENT'S MOMENTS PER UNIT LENGTH
C
        do 99942 j = 1, 3
          do g_i_ = 1, g_p_
            g_xkj(g_i_) = g_elecrv(g_i_, j)
          enddo
          xkj = elecrv(j)
C--------
          do g_i_ = 1, g_p_
            g_xsj(g_i_) = g_elestr(g_i_, j)
          enddo
          xsj = elestr(j)
C--------
          do 99943 i = 1, 3
            do g_i_ = 1, g_p_
              g_elem(g_i_, i) = cstbm(i, j) * g_xsj(g_i_) + xsj * g_cstb
     *m(g_i_, i, j) + cstbb(i, j) * g_xkj(g_i_) + xkj * g_cstbb(g_i_, i,
     * j) + g_elem(g_i_, i)
            enddo
            elem(i) = elem(i) + cstbb(i, j) * xkj + cstbm(i, j) * xsj
C--------
6102        continue
99943     continue
6101      continue
99942   continue
C
C.....ESTIMATE THE ELEMENT'S FORCES PER UNIT LENGTH
C
        do 99940 j = 1, 3
          do g_i_ = 1, g_p_
            g_xkj(g_i_) = g_elecrv(g_i_, j)
          enddo
          xkj = elecrv(j)
C--------
          do g_i_ = 1, g_p_
            g_xsj(g_i_) = g_elestr(g_i_, j)
          enddo
          xsj = elestr(j)
C--------
          do 99941 i = 1, 3
            do g_i_ = 1, g_p_
              g_elen(g_i_, i) = cstmm(i, j) * g_xsj(g_i_) + xsj * g_cstm
     *m(g_i_, i, j) + cstmb(i, j) * g_xkj(g_i_) + xkj * g_cstmb(g_i_, i,
     * j) + g_elen(g_i_, i)
            enddo
            elen(i) = elen(i) + cstmb(i, j) * xkj + cstmm(i, j) * xsj
C--------
6104        continue
99941     continue
6103      continue
99940   continue
C
C.....INITIALIZE VON MISES STRESS RESULTANT TO ZERO
C
        do g_i_ = 1, g_p_
          g_vonmises(g_i_) = 0.0d0
        enddo
        vonmises = zero
C--------
C
C.....ESTIMATE THE STRESSES ON THE UPPER SURFACE
C
        do 99939 i = 1, 3
          d3_v = elen(i) / thick
          d6_v = thick * thick
          d7_v = six * elem(i) / d6_v
          d5_b = -(-d7_v / d6_v)
          d7_b = -(1.0d0 / d6_v) * six
          d8_b = 1.0d0 / thick
          d6_b = d5_b * thick + d5_b * thick + (-d3_v) / thick
          do g_i_ = 1, g_p_
            g_ups(g_i_, i) = d7_b * g_elem(g_i_, i) + d6_b * g_thick(g_i
     *_) + d8_b * g_elen(g_i_, i)
          enddo
          ups(i) = d3_v - d7_v
C--------
6105      continue
99939   continue
C
C.....STORE SIGMAXX, SIGMAYY, SIGMAXY (UPPER)
C
        do 99938 i = 1, maxgus
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 1, i) = g_ups(g_i_, 1)
          enddo
          stress(elm, 1, i) = ups(1)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 2, i) = g_ups(g_i_, 2)
          enddo
          stress(elm, 2, i) = ups(2)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 4, i) = g_ups(g_i_, 3)
          enddo
          stress(elm, 4, i) = ups(3)
C--------
10        continue
99938   continue
C
C.....CALCULATE THE RADIUS OF MOHR CIRCLE ON THE UPPER SURFACE
C
        d4_v = (ups(1) - ups(2)) / two
        d6_v = (ups(1) - ups(2)) / two
        d4_b = d4_v * (1.0d0 / two)
        d7_b = d6_v * (1.0d0 / two)
        d5_b = d4_b + d7_b
        d6_b = -d4_b + (-d7_b)
        do g_i_ = 1, g_p_
          g_x1(g_i_) = d6_b * g_ups(g_i_, 2) + d5_b * g_ups(g_i_, 1)
        enddo
        x1 = d4_v * d6_v
C--------
        d2_b = ups(3) + ups(3)
        do g_i_ = 1, g_p_
          g_x2(g_i_) = d2_b * g_ups(g_i_, 3)
        enddo
        x2 = ups(3) * ups(3)
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = g_x2(g_i_) + g_x1(g_i_)
        enddo
        d1_w = x1 + x2
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +2161)
        endif
        do g_i_ = 1, g_p_
          g_upr(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        upr = d2_v
C--------
C
C.....CALCULATE THE PRINCIPAL STRESSES ON THE UPPER SURFACE
C
        d4_b = 1.0d0 / two
        do g_i_ = 1, g_p_
          g_ups1(g_i_) = g_upr(g_i_) + d4_b * g_ups(g_i_, 2) + d4_b * g_
     *ups(g_i_, 1)
        enddo
        ups1 = (ups(1) + ups(2)) / two + upr
C--------
        d4_b = 1.0d0 / two
        do g_i_ = 1, g_p_
          g_ups2(g_i_) = -g_upr(g_i_) + d4_b * g_ups(g_i_, 2) + d4_b * g
     *_ups(g_i_, 1)
        enddo
        ups2 = (ups(1) + ups(2)) / two - upr
C--------
C
C.....CALCULATE VON MISES STRESS OF THE UPPER LAYER
C
        d4_b = ups2 + ups2
        d5_b = ups1 + ups1
        do g_i_ = 1, g_p_
          g_x1(g_i_) = d4_b * g_ups2(g_i_) + d5_b * g_ups1(g_i_)
        enddo
        x1 = ups1 * ups1 + ups2 * ups2
C--------
        d3_v = ups1 - ups2
        d4_v = ups1 - ups2
        d4_b = d3_v + d4_v
        d5_b = -d3_v + (-d4_v)
        do g_i_ = 1, g_p_
          g_x2(g_i_) = d5_b * g_ups2(g_i_) + d4_b * g_ups1(g_i_)
        enddo
        x2 = d3_v * d4_v
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = g_x2(g_i_) + g_x1(g_i_)
        enddo
        d1_w = x1 + x2
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +2214)
        endif
        d3_b = 1.0d0 / sqrt(two) * d1_p
        do g_i_ = 1, g_p_
          g_upvms(g_i_) = d3_b * g_d1_w(g_i_)
        enddo
        upvms = d2_v / sqrt(two)
C--------
C
C.....ESTIMATE THE STRESSES ON THE LOWER SURFACE
C
        do 99937 i = 1, 3
          d3_v = elen(i) / thick
          d6_v = thick * thick
          d7_v = six * elem(i) / d6_v
          d5_b = -d7_v / d6_v
          d7_b = 1.0d0 / d6_v * six
          d8_b = 1.0d0 / thick
          d6_b = d5_b * thick + d5_b * thick + (-d3_v) / thick
          do g_i_ = 1, g_p_
            g_lws(g_i_, i) = d7_b * g_elem(g_i_, i) + d6_b * g_thick(g_i
     *_) + d8_b * g_elen(g_i_, i)
          enddo
          lws(i) = d3_v + d7_v
C--------
6106      continue
99937   continue
C
C.....CALCULATE THE RADIUS OF MOHR CIRCLE ON THE LOWER SURFACE
C
        d4_v = (lws(1) - lws(2)) / two
        d6_v = (lws(1) - lws(2)) / two
        d4_b = d4_v * (1.0d0 / two)
        d7_b = d6_v * (1.0d0 / two)
        d5_b = d4_b + d7_b
        d6_b = -d4_b + (-d7_b)
        do g_i_ = 1, g_p_
          g_x1(g_i_) = d6_b * g_lws(g_i_, 2) + d5_b * g_lws(g_i_, 1)
        enddo
        x1 = d4_v * d6_v
C--------
        d2_b = lws(3) + lws(3)
        do g_i_ = 1, g_p_
          g_x2(g_i_) = d2_b * g_lws(g_i_, 3)
        enddo
        x2 = lws(3) * lws(3)
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = g_x2(g_i_) + g_x1(g_i_)
        enddo
        d1_w = x1 + x2
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +2271)
        endif
        do g_i_ = 1, g_p_
          g_lwr(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        lwr = d2_v
C--------
C
C.....CALCULATE THE PRINCIPAL STRESSES ON THE LOWER SURFACE
C
        d4_b = 1.0d0 / two
        do g_i_ = 1, g_p_
          g_lws1(g_i_) = g_lwr(g_i_) + d4_b * g_lws(g_i_, 2) + d4_b * g_
     *lws(g_i_, 1)
        enddo
        lws1 = (lws(1) + lws(2)) / two + lwr
C--------
        d4_b = 1.0d0 / two
        do g_i_ = 1, g_p_
          g_lws2(g_i_) = -g_lwr(g_i_) + d4_b * g_lws(g_i_, 2) + d4_b * g
     *_lws(g_i_, 1)
        enddo
        lws2 = (lws(1) + lws(2)) / two - lwr
C--------
C
C.....CALCULATE VON MISES STRESS OF THE LOWER SURFACE
C
        d4_b = lws2 + lws2
        d5_b = lws1 + lws1
        do g_i_ = 1, g_p_
          g_x1(g_i_) = d4_b * g_lws2(g_i_) + d5_b * g_lws1(g_i_)
        enddo
        x1 = lws1 * lws1 + lws2 * lws2
C--------
        d3_v = lws1 - lws2
        d4_v = lws1 - lws2
        d4_b = d3_v + d4_v
        d5_b = -d3_v + (-d4_v)
        do g_i_ = 1, g_p_
          g_x2(g_i_) = d5_b * g_lws2(g_i_) + d4_b * g_lws1(g_i_)
        enddo
        x2 = d3_v * d4_v
C--------
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = g_x2(g_i_) + g_x1(g_i_)
        enddo
        d1_w = x1 + x2
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +2324)
        endif
        d3_b = 1.0d0 / sqrt(two) * d1_p
        do g_i_ = 1, g_p_
          g_lwvms(g_i_) = d3_b * g_d1_w(g_i_)
        enddo
        lwvms = d2_v / sqrt(two)
C--------
C
C.....THE VON MISES STRESS IS THE MAXIMUM OF THE TWO
C
        d3_v = max (upvms, lwvms)
        if (upvms .gt.  lwvms) then
           d1_p = 1.0d0
           d2_p = 0.0d0
        else if (upvms .lt.  lwvms) then
           d1_p = 0.0d0
           d2_p = 1.0d0
        else
           call ehbfDV (7,upvms, lwvms, d3_v, d1_p, d2_p,
     +'g_compvms.f',
     +2345)
           d2_p = 1.0d0 -  d1_p
        endif
        do g_i_ = 1, g_p_
          g_vonmises(g_i_) = d2_p * g_lwvms(g_i_) + d1_p * g_upvms(g_i_)
        enddo
        vonmises = d3_v
C--------
C
C.....IF UPPER SURFACE IS REQUESTED, RETURN upvms
C
        if (surface .eq. 1) then
          do g_i_ = 1, g_p_
            g_vonmises(g_i_) = g_upvms(g_i_)
          enddo
          vonmises = upvms
C--------
        endif
C
C.....IF lOWER SURFACE IS REQUESTED, RETURN lwvms
C
        if (surface .eq. 3) then
          do g_i_ = 1, g_p_
            g_vonmises(g_i_) = g_lwvms(g_i_)
          enddo
          vonmises = lwvms
C--------
          do 99936 i = 1, maxgus
            do g_i_ = 1, g_p_
              g_stress(g_i_, elm, 1, i) = g_lws(g_i_, 1)
            enddo
            stress(elm, 1, i) = lws(1)
C--------
            do g_i_ = 1, g_p_
              g_stress(g_i_, elm, 2, i) = g_lws(g_i_, 2)
            enddo
            stress(elm, 2, i) = lws(2)
C--------
            do g_i_ = 1, g_p_
              g_stress(g_i_, elm, 4, i) = g_lws(g_i_, 3)
            enddo
            stress(elm, 4, i) = lws(3)
C--------
11          continue
99936     continue
        endif
C
C
C.....IF MEDIAN SURFACE IS REQUESTED, RETURN mds
C
        if (surface .eq. 2) then
          do i = 1, 3
            d3_v = elen(i) / thick
            d2_b = 1.0d0 / thick
            d3_b = -d3_v / thick
            do g_i_ = 1, g_p_
              g_mds(g_i_, i) = d3_b * g_thick(g_i_) + d2_b * g_elen(g_i_
     *, i)
            enddo
            mds(i) = d3_v
C--------
          enddo
C
          do 99935 i = 1, maxgus
            do g_i_ = 1, g_p_
              g_stress(g_i_, elm, 1, i) = g_mds(g_i_, 1)
            enddo
            stress(elm, 1, i) = mds(1)
C--------
            do g_i_ = 1, g_p_
              g_stress(g_i_, elm, 2, i) = g_mds(g_i_, 2)
            enddo
            stress(elm, 2, i) = mds(2)
C--------
            do g_i_ = 1, g_p_
              g_stress(g_i_, elm, 4, i) = g_mds(g_i_, 3)
            enddo
            stress(elm, 4, i) = mds(3)
C--------
12          continue
99935     continue
C
C
C.....CALCULATE THE RADIUS OF MOHR CIRCLE ON THE MEDIAN SURFACE
C
          d4_v = (mds(1) - mds(2)) / two
          d6_v = (mds(1) - mds(2)) / two
          d4_b = d4_v * (1.0d0 / two)
          d7_b = d6_v * (1.0d0 / two)
          d5_b = d4_b + d7_b
          d6_b = -d4_b + (-d7_b)
          do g_i_ = 1, g_p_
            g_x1(g_i_) = d6_b * g_mds(g_i_, 2) + d5_b * g_mds(g_i_, 1)
          enddo
          x1 = d4_v * d6_v
C--------
          d2_b = mds(3) + mds(3)
          do g_i_ = 1, g_p_
            g_x2(g_i_) = d2_b * g_mds(g_i_, 3)
          enddo
          x2 = mds(3) * mds(3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = g_x2(g_i_) + g_x1(g_i_)
          enddo
          d1_w = x1 + x2
          d2_v = sqrt(d1_w)
          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d2_v)
          else
             call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +2457)
          endif
          do g_i_ = 1, g_p_
            g_mdr(g_i_) = d1_p * g_d1_w(g_i_)
          enddo
          mdr = d2_v
C--------
C
C.....CALCULATE THE PRINCIPAL STRESSES ON THE MEDIAN SURFACE
C
          d4_b = 1.0d0 / two
          do g_i_ = 1, g_p_
            g_mds1(g_i_) = g_mdr(g_i_) + d4_b * g_mds(g_i_, 2) + d4_b * 
     *g_mds(g_i_, 1)
          enddo
          mds1 = (mds(1) + mds(2)) / two + mdr
C--------
          d4_b = 1.0d0 / two
          do g_i_ = 1, g_p_
            g_mds2(g_i_) = -g_mdr(g_i_) + d4_b * g_mds(g_i_, 2) + d4_b *
     * g_mds(g_i_, 1)
          enddo
          mds2 = (mds(1) + mds(2)) / two - mdr
C--------
C
C.....CALCULATE VON MISES STRESS OF THE MEDIAN LAYER
C
          d4_b = mds2 + mds2
          d5_b = mds1 + mds1
          do g_i_ = 1, g_p_
            g_x1(g_i_) = d4_b * g_mds2(g_i_) + d5_b * g_mds1(g_i_)
          enddo
          x1 = mds1 * mds1 + mds2 * mds2
C--------
          d3_v = mds1 - mds2
          d4_v = mds1 - mds2
          d4_b = d3_v + d4_v
          d5_b = -d3_v + (-d4_v)
          do g_i_ = 1, g_p_
            g_x2(g_i_) = d5_b * g_mds2(g_i_) + d4_b * g_mds1(g_i_)
          enddo
          x2 = d3_v * d4_v
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = g_x2(g_i_) + g_x1(g_i_)
          enddo
          d1_w = x1 + x2
          d2_v = sqrt(d1_w)
          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d2_v)
          else
             call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +2510)
          endif
          d3_b = 1.0d0 / sqrt(two) * d1_p
          do g_i_ = 1, g_p_
            g_mdvms(g_i_) = d3_b * g_d1_w(g_i_)
          enddo
          mdvms = d2_v / sqrt(two)
C--------
C
          do g_i_ = 1, g_p_
            g_vonmises(g_i_) = g_mdvms(g_i_)
          enddo
          vonmises = mdvms
C--------
C
        endif
C
C
C.....ROTATE LOCAL STRESSES TO GLOBAL
C
        do g_i_ = 1, g_p_
          g_str(g_i_, 1) = g_stress(g_i_, elm, 1, 1)
        enddo
        str(1) = stress(elm, 1, 1)
C--------
        do g_i_ = 1, g_p_
          g_str(g_i_, 2) = g_stress(g_i_, elm, 2, 1)
        enddo
        str(2) = stress(elm, 2, 1)
C--------
        do g_i_ = 1, g_p_
          g_str(g_i_, 3) = 0.0d0
        enddo
        str(3) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_str(g_i_, 4) = g_stress(g_i_, elm, 4, 1)
        enddo
        str(4) = stress(elm, 4, 1)
C--------
        do g_i_ = 1, g_p_
          g_str(g_i_, 5) = 0.0d0
        enddo
        str(5) = 0.0d0
C--------
        do g_i_ = 1, g_p_
          g_str(g_i_, 6) = 0.0d0
        enddo
        str(6) = 0.0d0
C--------
C
        call g_transform(g_p_, xp, g_xp, g_pmax_, yp, g_yp, g_pmax_, zp,
     * g_zp, g_pmax_, xg, yg, zg, str, g_str, g_pmax_)
C
        do 99934 i = 1, 6
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, i, 1) = g_str(g_i_, i)
          enddo
          stress(elm, i, 1) = str(i)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, i, 2) = g_str(g_i_, i)
          enddo
          stress(elm, i, 2) = str(i)
C--------
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, i, 3) = g_str(g_i_, i)
          enddo
          stress(elm, i, 3) = str(i)
C--------
102       continue
99934   continue
C
C
C
C.....STORE VON MISES STRESS FOR EACH ONE OF THE GAUSS POINT
C
        do 99933 i = 1, maxgus
          do g_i_ = 1, g_p_
            g_stress(g_i_, elm, 7, i) = g_vonmises(g_i_)
          enddo
          stress(elm, 7, i) = vonmises
C--------
6107      continue
99933   continue
C
C     -------------------------------------------------
C     STEP 7
C     COMPUTE THE ELEMENTAL MOMENT AND FORCE
C     PER UNIT LENGTH VECTORS PER LAYER
C     COMPUTE THE CENTROIDAL VON MISES STRESS RESULTANT
C     LAYER PER LAYER IN THE COMPOSITE SHELL ELEMENT
C     -------------------------------------------------
C
        if ((ctyp .eq. 2) .or. (ctyp .eq. 3)) then
C
C.....LOOP ON THE LAYERS OF THE COMPOSITE SHELL ELEMENT
C
          do 99915 ilayer = 1, nlayer
C
C.....EXTRACT THE LAYER POSITION
C
            layerpos = 0
C
            detecterror = .false.
C
            do 99932 i = 1, nfely
              if (laysgid(1, i) .eq. elm) then
                if (laysgid(5, i) .eq. ilayer) then
                  if (detecterror) then
                    goto 940
                  endif
                  layerpos = i
                  detecterror = .true.
                endif
              endif
7002          continue
99932       continue
C
C.....INITIALIZE THE LAYER'S CONSTANT THICKNESS
C
            do g_i_ = 1, g_p_
              g_thick(g_i_) = g_mtcmp23(g_i_, 7, ilayer)
            enddo
            thick = mtcmp23(7, ilayer)
C--------
C
C.....GET THE CONSTITUTIVE LAW FOR THE LAYER
C
            do 99927 j = 1, 3
              do 99931 i = 1, 3
                do g_i_ = 1, g_p_
                  g_cstbb(g_i_, i, j) = 0.0d0
                enddo
                cstbb(i, j) = zero
C--------
7004            continue
99931         continue
              do 99930 i = 1, 3
                do g_i_ = 1, g_p_
                  g_cstmm(g_i_, i, j) = 0.0d0
                enddo
                cstmm(i, j) = zero
C--------
7005            continue
99930         continue
              do 99929 i = 1, 3
                do g_i_ = 1, g_p_
                  g_cstbm(g_i_, i, j) = 0.0d0
                enddo
                cstbm(i, j) = zero
C--------
7006            continue
99929         continue
              do 99928 i = 1, 3
                do g_i_ = 1, g_p_
                  g_cstmb(g_i_, i, j) = 0.0d0
                enddo
                cstmb(i, j) = zero
C--------
7007            continue
99928         continue
7003          continue
99927       continue
C
            call g_complay(g_p_, nlayer, idcmp23, mtcmp23, g_mtcmp23, g_
     *pmax_, x, g_x, g_pmax_, y, g_y, g_pmax_, z, ctyp, ilayer, cstbb, g
     *_cstbb, g_pmax_, cstmm, g_cstmm, g_pmax_, cstbm, g_cstbm, g_pmax_,
     * cstmb, g_cstmb, g_pmax_, eframe, g_eframe, g_pmax_, aframe)
C
C.....ESTIMATE THE LAYER'S MOMENTS PER UNIT LENGTH
C
            do 99926 i = 1, 3
              do g_i_ = 1, g_p_
                g_elem(g_i_, i) = 0.0d0
              enddo
              elem(i) = zero
C--------
7101          continue
99926       continue
C
            do 99924 j = 1, 3
              do g_i_ = 1, g_p_
                g_xkj(g_i_) = g_elecrv(g_i_, j)
              enddo
              xkj = elecrv(j)
C--------
              do g_i_ = 1, g_p_
                g_xsj(g_i_) = g_elestr(g_i_, j)
              enddo
              xsj = elestr(j)
C--------
              do 99925 i = 1, 3
                do g_i_ = 1, g_p_
                  g_elem(g_i_, i) = cstbm(i, j) * g_xsj(g_i_) + xsj * g_
     *cstbm(g_i_, i, j) + cstbb(i, j) * g_xkj(g_i_) + xkj * g_cstbb(g_i_
     *, i, j) + g_elem(g_i_, i)
                enddo
                elem(i) = elem(i) + cstbb(i, j) * xkj + cstbm(i, j) * xs
     *j
C--------
7103            continue
99925         continue
7102          continue
99924       continue
C
C.....ESTIMATE THE LAYER'S FORCES PER UNIT LENGTH
C
            do 99923 i = 1, 3
              do g_i_ = 1, g_p_
                g_elen(g_i_, i) = 0.0d0
              enddo
              elen(i) = zero
C--------
7201          continue
99923       continue
C
            do 99921 j = 1, 3
              do g_i_ = 1, g_p_
                g_xkj(g_i_) = g_elecrv(g_i_, j)
              enddo
              xkj = elecrv(j)
C--------
              do g_i_ = 1, g_p_
                g_xsj(g_i_) = g_elestr(g_i_, j)
              enddo
              xsj = elestr(j)
C--------
              do 99922 i = 1, 3
                do g_i_ = 1, g_p_
                  g_elen(g_i_, i) = cstmm(i, j) * g_xsj(g_i_) + xsj * g_
     *cstmm(g_i_, i, j) + cstmb(i, j) * g_xkj(g_i_) + xkj * g_cstmb(g_i_
     *, i, j) + g_elen(g_i_, i)
                enddo
                elen(i) = elen(i) + cstmb(i, j) * xkj + cstmm(i, j) * xs
     *j
C--------
7203            continue
99922         continue
7202          continue
99921       continue
C
C.....INITIALIZE VON MISES STRESS RESULTANT TO ZERO
C
            do g_i_ = 1, g_p_
              g_vonmises(g_i_) = 0.0d0
            enddo
            vonmises = zero
C--------
C
C.....ESTIMATE THE STRESSES ON THE LAYER'S UPPER SURFACE
C
            do 99920 i = 1, 3
              do g_i_ = 1, g_p_
                g_ups(g_i_, i) = 0.0d0
              enddo
              ups(i) = zero
C--------
7301          continue
99920       continue
C
            do 99919 i = 1, 3
              d3_v = elen(i) / thick
              d6_v = thick * thick
              d7_v = six * elem(i) / d6_v
              d5_b = -(-d7_v / d6_v)
              d7_b = -(1.0d0 / d6_v) * six
              d8_b = 1.0d0 / thick
              d6_b = d5_b * thick + d5_b * thick + (-d3_v) / thick
              do g_i_ = 1, g_p_
                g_ups(g_i_, i) = d7_b * g_elem(g_i_, i) + d6_b * g_thick
     *(g_i_) + d8_b * g_elen(g_i_, i)
              enddo
              ups(i) = d3_v - d7_v
C--------
7302          continue
99919       continue
C
C.....CALCULATE THE RADIUS OF MOHR CIRCLE ON THE LAYER'S UPPER SURFACE
C
            d4_v = (ups(1) - ups(2)) / two
            d6_v = (ups(1) - ups(2)) / two
            d4_b = d4_v * (1.0d0 / two)
            d7_b = d6_v * (1.0d0 / two)
            d5_b = d4_b + d7_b
            d6_b = -d4_b + (-d7_b)
            do g_i_ = 1, g_p_
              g_x1(g_i_) = d6_b * g_ups(g_i_, 2) + d5_b * g_ups(g_i_, 1)
            enddo
            x1 = d4_v * d6_v
C--------
            d2_b = ups(3) + ups(3)
            do g_i_ = 1, g_p_
              g_x2(g_i_) = d2_b * g_ups(g_i_, 3)
            enddo
            x2 = ups(3) * ups(3)
C--------
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = g_x2(g_i_) + g_x1(g_i_)
            enddo
            d1_w = x1 + x2
            d2_v = sqrt(d1_w)
            if ( d1_w .gt. 0.0d0 ) then
               d1_p = 1.0d0 / (2.0d0 *  d2_v)
            else
               call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +2817)
            endif
            do g_i_ = 1, g_p_
              g_upr(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            upr = d2_v
C--------
C
C.....CALCULATE THE PRINCIPAL STRESSES ON THE LAYER'S UPPER SURFACE
C
            d4_b = 1.0d0 / two
            do g_i_ = 1, g_p_
              g_ups1(g_i_) = g_upr(g_i_) + d4_b * g_ups(g_i_, 2) + d4_b 
     ** g_ups(g_i_, 1)
            enddo
            ups1 = (ups(1) + ups(2)) / two + upr
C--------
            d4_b = 1.0d0 / two
            do g_i_ = 1, g_p_
              g_ups2(g_i_) = -g_upr(g_i_) + d4_b * g_ups(g_i_, 2) + d4_b
     * * g_ups(g_i_, 1)
            enddo
            ups2 = (ups(1) + ups(2)) / two - upr
C--------
C
C.....CALCULATE VON MISES STRESS OF THE LAYER'S UPPER LAYER
C
            d4_b = ups2 + ups2
            d5_b = ups1 + ups1
            do g_i_ = 1, g_p_
              g_x1(g_i_) = d4_b * g_ups2(g_i_) + d5_b * g_ups1(g_i_)
            enddo
            x1 = ups1 * ups1 + ups2 * ups2
C--------
            d3_v = ups1 - ups2
            d4_v = ups1 - ups2
            d4_b = d3_v + d4_v
            d5_b = -d3_v + (-d4_v)
            do g_i_ = 1, g_p_
              g_x2(g_i_) = d5_b * g_ups2(g_i_) + d4_b * g_ups1(g_i_)
            enddo
            x2 = d3_v * d4_v
C--------
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = g_x2(g_i_) + g_x1(g_i_)
            enddo
            d1_w = x1 + x2
            d2_v = sqrt(d1_w)
            if ( d1_w .gt. 0.0d0 ) then
               d1_p = 1.0d0 / (2.0d0 *  d2_v)
            else
               call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +2870)
            endif
            d3_b = 1.0d0 / sqrt(two) * d1_p
            do g_i_ = 1, g_p_
              g_upvms(g_i_) = d3_b * g_d1_w(g_i_)
            enddo
            upvms = d2_v / sqrt(two)
C--------
C
C.....ESTIMATE THE STRESSES ON THE LAYER'S LOWER SURFACE
C
            do 99918 i = 1, 3
              do g_i_ = 1, g_p_
                g_lws(g_i_, i) = 0.0d0
              enddo
              lws(i) = zero
C--------
7401          continue
99918       continue
C
            do 99917 i = 1, 3
              d3_v = elen(i) / thick
              d6_v = thick * thick
              d7_v = six * elem(i) / d6_v
              d5_b = -d7_v / d6_v
              d7_b = 1.0d0 / d6_v * six
              d8_b = 1.0d0 / thick
              d6_b = d5_b * thick + d5_b * thick + (-d3_v) / thick
              do g_i_ = 1, g_p_
                g_lws(g_i_, i) = d7_b * g_elem(g_i_, i) + d6_b * g_thick
     *(g_i_) + d8_b * g_elen(g_i_, i)
              enddo
              lws(i) = d3_v + d7_v
C--------
7402          continue
99917       continue
C
C.....CALCULATE THE RADIUS OF MOHR CIRCLE ON THE LAYER'S LOWER SURFACE
C
            d4_v = (lws(1) - lws(2)) / two
            d6_v = (lws(1) - lws(2)) / two
            d4_b = d4_v * (1.0d0 / two)
            d7_b = d6_v * (1.0d0 / two)
            d5_b = d4_b + d7_b
            d6_b = -d4_b + (-d7_b)
            do g_i_ = 1, g_p_
              g_x1(g_i_) = d6_b * g_lws(g_i_, 2) + d5_b * g_lws(g_i_, 1)
            enddo
            x1 = d4_v * d6_v
C--------
            d2_b = lws(3) + lws(3)
            do g_i_ = 1, g_p_
              g_x2(g_i_) = d2_b * g_lws(g_i_, 3)
            enddo
            x2 = lws(3) * lws(3)
C--------
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = g_x2(g_i_) + g_x1(g_i_)
            enddo
            d1_w = x1 + x2
            d2_v = sqrt(d1_w)
            if ( d1_w .gt. 0.0d0 ) then
               d1_p = 1.0d0 / (2.0d0 *  d2_v)
            else
               call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +2936)
            endif
            do g_i_ = 1, g_p_
              g_lwr(g_i_) = d1_p * g_d1_w(g_i_)
            enddo
            lwr = d2_v
C--------
C
C.....CALCULATE THE PRINCIPAL STRESSES ON THE LAYER'S LOWER SURFACE
C
            d4_b = 1.0d0 / two
            do g_i_ = 1, g_p_
              g_lws1(g_i_) = g_lwr(g_i_) + d4_b * g_lws(g_i_, 2) + d4_b 
     ** g_lws(g_i_, 1)
            enddo
            lws1 = (lws(1) + lws(2)) / two + lwr
C--------
            d4_b = 1.0d0 / two
            do g_i_ = 1, g_p_
              g_lws2(g_i_) = -g_lwr(g_i_) + d4_b * g_lws(g_i_, 2) + d4_b
     * * g_lws(g_i_, 1)
            enddo
            lws2 = (lws(1) + lws(2)) / two - lwr
C--------
C
C.....CALCULATE VON MISES STRESS OF THE LAYER'S LOWER LAYER
C
            d4_b = lws2 + lws2
            d5_b = lws1 + lws1
            do g_i_ = 1, g_p_
              g_x1(g_i_) = d4_b * g_lws2(g_i_) + d5_b * g_lws1(g_i_)
            enddo
            x1 = lws1 * lws1 + lws2 * lws2
C--------
            d3_v = lws1 - lws2
            d4_v = lws1 - lws2
            d4_b = d3_v + d4_v
            d5_b = -d3_v + (-d4_v)
            do g_i_ = 1, g_p_
              g_x2(g_i_) = d5_b * g_lws2(g_i_) + d4_b * g_lws1(g_i_)
            enddo
            x2 = d3_v * d4_v
C--------
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = g_x2(g_i_) + g_x1(g_i_)
            enddo
            d1_w = x1 + x2
            d2_v = sqrt(d1_w)
            if ( d1_w .gt. 0.0d0 ) then
               d1_p = 1.0d0 / (2.0d0 *  d2_v)
            else
               call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compvms.f',
     +2989)
            endif
            d3_b = 1.0d0 / sqrt(two) * d1_p
            do g_i_ = 1, g_p_
              g_lwvms(g_i_) = d3_b * g_d1_w(g_i_)
            enddo
            lwvms = d2_v / sqrt(two)
C--------
C
C.....THE LAYER'S VON MISES STRESS IS THE MAXIMUM OF THE TWO
C
            d3_v = max (upvms, lwvms)
            if (upvms .gt.  lwvms) then
               d1_p = 1.0d0
               d2_p = 0.0d0
            else if (upvms .lt.  lwvms) then
               d1_p = 0.0d0
               d2_p = 1.0d0
            else
               call ehbfDV (7,upvms, lwvms, d3_v, d1_p, d2_p,
     +'g_compvms.f',
     +3010)
               d2_p = 1.0d0 -  d1_p
            endif
            do g_i_ = 1, g_p_
              g_vonmises(g_i_) = d2_p * g_lwvms(g_i_) + d1_p * g_upvms(g
     *_i_)
            enddo
            vonmises = d3_v
C--------
C
C.....STORE THE LAYER'S VON MISES STRESS IN 7TH POSITION
C.....(POSITIONS 1-6 ARE FOR Sx Sy Sz Sxy Sxz Syz) FOR EACH
C.....ONE OF THE [maxgly] GAUSS POINTS
C
            do 99916 i = 1, maxgly
              laysigm(layerpos, 7, i) = vonmises
C*ADIFOR* I/O statement contains active variables
              write (*, *) vonmises
7501          continue
99916       continue
C
C.....END OF LOOP ON THE LAYERS OF THE COMPOSITE SHELL ELEMENT
C
7001        continue
99915     continue
C
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
        write (*, *) '*** FATAL ERROR in Routine COMPVMS ***'
        write (*, *) '*** The Type of Constitutive Law   ***'
        write (*, *) '*** Must Either be 0, 1, 2, or 3.  ***'
        write (*, *) '*** Execution Terminated Here      ***'
        stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [CMPCO] IS NOT CORRECT
C
200     continue
        write (*, *) '*** FATAL ERROR in Routine COMPVMS ***'
        write (*, *) '*** The Address in Array [CMPCO]   ***'
        write (*, *) '*** is Out-of-Bounds.              ***'
        write (*, *) '*** Execution Terminated Here      ***'
        stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [IDLAY]/[MTLAY] IS NOT CORRECT
C
300     continue
        write (*, *) '*** FATAL ERROR in Routine COMPVMS ***'
        write (*, *) '*** The Address in Arrays [IDLAY]  ***'
        write (*, *) '*** and [MTLAY] is Out-of-Bounds.  ***'
        write (*, *) '*** Execution Terminated Here      ***'
        stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [CMPFR] IS NOT CORRECT
C
400     continue
        write (*, *) '*** FATAL ERROR in Routine COMPVMS ***'
        write (*, *) '*** The Address in Array [CMPFR]   ***'
        write (*, *) '*** is Out-of-Bounds.              ***'
        write (*, *) '*** Execution Terminated Here      ***'
        stop
C
C.....ERROR-MESSAGE IF THE MAXIMUM NUMBER OF LAYERS HAS BEEN EXCEEDED
C
500     continue
        write (*, *) '*** FATAL ERROR in Routine COMPVMS        ***'
        write (*, *) '*** The Maximum Number of Layers Has Been ***'
        write (*, *) '*** Exceeded: Boost Parameter [MAXLAYER]. ***'
        write (*, *) '*** Execution Terminated Here             ***'
        stop
C
C.....ERROR-MESSAGE IF THE TOTAL NUMBER OF LAYERS IS NOT CONSISTENT
C
600     continue
        write (*, *) '*** FATAL ERROR in Routine COMPVMS      ***'
        write (*, *) '*** A Layer Number Exceeds the Total    ***'
        write (*, *) '*** Number of Layers Stored in [IDLAY]. ***'
        write (*, *) '*** Execution Terminated Here           ***'
        stop
C
C.....ERROR-MESSAGE IF BASIC PARAMETERS FOR TYPES 2 & 3 ARE INCONSISTENT
C
700     continue
        write (*, *) '*** FATAL ERROR in Routine COMPVMS        ***'
        write (*, *) '*** The Element Numbers and Total Number  ***'
        write (*, *) '*** of Composite Layers Do Not Check Out. ***'
        write (*, *) '*** Execution Terminated Here             ***'
        stop
C
C.....ERROR-MESSAGE IF THE FIRST COEFFICIENT OF EXTENTION IS ZERO
C
800     continue
        write (*, *) '*** FATAL ERROR in routine COMPVMS       ***'
        write (*, *) '*** The First Coefficient of Extentional ***'
        write (*, *) '*** Stiffness Cbb(1,1) is Equal to Zero! ***'
        write (*, *) '*** Can Not Estimate the Thickness...    ***'
        write (*, *) '*** EXECUTION TERMINATED RIGHT HERE      ***'
        stop
C
C.....ERROR-MESSAGE IF THE RATIO BETWEEN FIRST COEFFICIENTS OF
C.....BENDING STIFFNESS OVER EXTENTIONAL STIFFNESS IS NEGATIVE
C
900     continue
        write (*, *) '*** FATAL ERROR in routine COMPVMS            ***'
        write (*, *) '*** The Ratio Between the First Coefficient   ***'
        write (*, *) '*** of Bending Stiffness and the First        ***'
        write (*, *) '*** Coefficient of Extentional Stiffness is   ***'
        write (*, *) '*** Negative or Zero: Can Not Take the Square ***'
        write (*, *) '*** Root and Estimate the Shell Thickness.    ***'
        write (*, *) '*** EXECUTION TERMINATED RIGHT HERE           ***'
        stop
C
C.....ERROR-MESSAGE IF [CLR]+[CQR] IS DIFFERENT FROM ONE
C
910     continue
        write (*, *) '*** FATAL ERROR in routine COMPVMS      ***'
        write (*, *) '*** The Factors [clr] and [cqr] Violate ***'
        write (*, *) '*** the Constraint [clr]+[cqr]=1:       ***'
        write (*, *) '*** Check the Calling Sequence.         ***'
        write (*, *) '*** EXECUTION TERMINATED RIGHT HERE     ***'
        stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
920     continue
        write (*, *) '*** FATAL ERROR in routine COMPVMS         ***'
        write (*, *) '*** The Triangle Area is Found Negative or ***'
        write (*, *) '*** Zero: Check the Nodal Point Numbering  ***'
        write (*, *) '*** ... Counterclockwise?                  ***'
        write (*, *) '*** EXECUTION TERMINATED RIGHT HERE        ***'
        stop
C
C.....ERROR-MESSAGE IF [ALPHA] IS NEGATIVE
C
930     continue
        write (*, *) '*** FATAL ERROR in routine COMPVMS   ***'
        write (*, *) '*** The Factor [alpha] is Negative   ***'
        write (*, *) '*** Check the Data Sequence: Factor  ***'
        write (*, *) '*** [alpha] Must be Positive or Zero ***'
        write (*, *) '*** EXECUTION TERMINATED RIGHT HERE  ***'
        stop
C
C.....ERROR-MESSAGE IF THE LAYER NUMBER IS NOT UNIQUE
C
940     continue
        write (*, *) '*** FATAL ERROR in routine COMPVMS       ***'
        write (*, *) '*** The Same Layer Number Is Encountered ***'
        write (*, *) '*** More than Once for a Given Element   ***'
        write (*, *) '*** Can Not Happen: Check Input File.    ***'
        write (*, *) '*** EXECUTION TERMINATED RIGHT HERE      ***'
        stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPVMS"
C========================C
