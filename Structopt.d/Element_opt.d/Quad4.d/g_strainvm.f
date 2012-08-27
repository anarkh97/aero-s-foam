C                           DISCLAIMER
C
C   This file was generated on 06/27/02 by the version of
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
      subroutine gxstrainvm(strain, g_strain, maxgus,
     *maxstr, msize, numnod)
C**************************************************************
C       THIS ROUTINE CALCULATES EQUIVLAENT VON MISES STRAIN   *
C**************************************************************
C                                                             *
C     AUTHOR : K.H. PIERSON          *
C     DATE : MARCH 1997         *
C     VERSION :       FEM-C++ 1.00            *
C
C       modified by:
C
C       AUTHOR  :       G. BROWN
C       DATE    :       DEC. 2000
C                                                             *
C**************************************************************
C                                                             *
C     STRAIN = 3-D ARRRAY HOLDING ELEMENTAL STRAINS         *
C     MAXGUS = THIRD DIMENSION OF STRAIN                    *
C     MAXSTR = SECOND DIMENSION OF STRAIN ARRAY             *
C      MSIZE = LEADING DIMENSION OF STRAIN ARRAY            *
C       NUMNOD = # OF NODES IN THE ELEMENT                    *
C                                                             *
C**************************************************************
C                                                             *
C     CALLED BY :  SANDS19.F                                *
C                                                             *
C**************************************************************
C
C.... INTEGER CONSTANTS
C
        integer maxgus, maxstr, msize, numnod
C
C.... REAL CONSTANTS
C
        double precision nno
C
C.... REAL ARRAYS
C
        double precision strain(msize, maxstr, maxgus)
C
C.... DECLARE ALL LOCAL VARIABLES
C
        integer n
        double precision exx, eyy, ezz, exy, exz, eyz
        double precision dexx, deyy, dezz, dexy, dexz, deyz
        double precision j2, comp, vms
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_strain

        parameter (g_p_= 1, ldg_strain= 1)

        double precision d1_p, d8_b, d7_b, d4_b, d2_b, d2_w, d1_w, d3_v,
     * d5_b, g_vms(g_pmax_)
        double precision g_exx(g_pmax_), g_strain(ldg_strain, msize, max
     *str, maxgus), g_eyy(g_pmax_), g_ezz(g_pmax_), g_exy(g_pmax_), g_ey
     *z(g_pmax_), g_exz(g_pmax_), g_comp(g_pmax_), g_dexx(g_pmax_), g_de
     *yy(g_pmax_)
        double precision g_dezz(g_pmax_), g_dexy(g_pmax_), g_deyz(g_pmax
     *_), g_dexz(g_pmax_), g_d1_w(g_pmax_), g_d2_w(g_pmax_), g_j2(g_pmax
     *_)
        save g_dezz, g_dexy, g_deyz, g_dexz, g_d1_w, g_d2_w, g_j2
        save g_vms, g_exx, g_eyy, g_ezz, g_exy, g_eyz, g_exz, g_comp, g_
     *dexx, g_deyy
        integer g_ehfid
        data g_ehfid /0/
C
c        call ehsfid(g_ehfid, 'strainvm','g_strainvm.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do g_i_ = 1, g_p_
          g_vms(g_i_) = 0.0d0
        enddo
        vms = 0.0d0
C--------
        nno = int(numnod)
C
C.... LOOP OVER ALL NODES
C
        do 99999 n = 1, nno
C
          do g_i_ = 1, g_p_
            g_exx(g_i_) = g_strain(g_i_, 1, 1, n)
          enddo
          exx = strain(1, 1, n)
C--------
          do g_i_ = 1, g_p_
            g_eyy(g_i_) = g_strain(g_i_, 1, 2, n)
          enddo
          eyy = strain(1, 2, n)
C--------
          do g_i_ = 1, g_p_
            g_ezz(g_i_) = g_strain(g_i_, 1, 3, n)
          enddo
          ezz = strain(1, 3, n)
C--------
C convert to strain tensor from engineering strain
          d2_b = 1.0d0 / 2.0d0
          do g_i_ = 1, g_p_
            g_exy(g_i_) = d2_b * g_strain(g_i_, 1, 4, n)
          enddo
          exy = strain(1, 4, n) / 2.0d0
C--------
          d2_b = 1.0d0 / 2.0d0
          do g_i_ = 1, g_p_
            g_eyz(g_i_) = d2_b * g_strain(g_i_, 1, 5, n)
          enddo
          eyz = strain(1, 5, n) / 2.0d0
C--------
          d2_b = 1.0d0 / 2.0d0
          do g_i_ = 1, g_p_
            g_exz(g_i_) = d2_b * g_strain(g_i_, 1, 6, n)
          enddo
          exz = strain(1, 6, n) / 2.0d0
C--------
C
C
C ... COMPUTE THE MEAN HYDROSTATIC STRAIN
C
          d2_b = 1.0d0 / 3.0d0
          do g_i_ = 1, g_p_
            g_comp(g_i_) = d2_b * g_ezz(g_i_) + d2_b * g_eyy(g_i_) + d2_
     *b * g_exx(g_i_)
          enddo
          comp = (exx + eyy + ezz) / 3.0d0
C--------
C
C.... COMPUTE THE FIRST DEVIATORIC STRAINS
C
          do g_i_ = 1, g_p_
            g_dexx(g_i_) = -g_comp(g_i_) + g_exx(g_i_)
          enddo
          dexx = exx - comp
C--------
          do g_i_ = 1, g_p_
            g_deyy(g_i_) = -g_comp(g_i_) + g_eyy(g_i_)
          enddo
          deyy = eyy - comp
C--------
          do g_i_ = 1, g_p_
            g_dezz(g_i_) = -g_comp(g_i_) + g_ezz(g_i_)
          enddo
          dezz = ezz - comp
C--------
          do g_i_ = 1, g_p_
            g_dexy(g_i_) = g_exy(g_i_)
          enddo
          dexy = exy
C--------
          do g_i_ = 1, g_p_
            g_deyz(g_i_) = g_eyz(g_i_)
          enddo
          deyz = eyz
C--------
          do g_i_ = 1, g_p_
            g_dexz(g_i_) = g_exz(g_i_)
          enddo
          dexz = exz
C--------
C
C.... COMPUTE THE SECOND DEVIATORIC STRAIN
C
          d4_b = deyy + deyy
          d5_b = dexx + dexx
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = d4_b * g_deyy(g_i_) + d5_b * g_dexx(g_i_)
          enddo
          d1_w = dexx * dexx + deyy * deyy
          d4_b = dexy + dexy
          d5_b = 1.0d0 / 2.0d0
          d8_b = d5_b * dezz + d5_b * dezz
          do g_i_ = 1, g_p_
            g_d2_w(g_i_) = d4_b * g_dexy(g_i_) + d8_b * g_dezz(g_i_) + d
     *5_b * g_d1_w(g_i_)
          enddo
          d2_w = (d1_w + dezz * dezz) / 2.0d0 + dexy * dexy
          d4_b = dexz + dexz
          d7_b = deyz + deyz
          do g_i_ = 1, g_p_
            g_j2(g_i_) = d4_b * g_dexz(g_i_) + d7_b * g_deyz(g_i_) + g_d
     *2_w(g_i_)
          enddo
          j2 = d2_w + deyz * deyz + dexz * dexz
C--------
C
C.... COMPUTE THE VON MISES STRAIN
C
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = 3.0d0 * g_j2(g_i_)
          enddo
          d1_w = 3.0d0 * j2
          d3_v = sqrt(d1_w)
          if ( d1_w .gt. 0.0d0 ) then
             d1_p = 1.0d0 / (2.0d0 *  d3_v)
          else
             call ehufDO (9,d1_w, d3_v, d1_p,
     +g_ehfid,
     +215)
          endif
          do g_i_ = 1, g_p_
            g_vms(g_i_) = d1_p * g_d1_w(g_i_) + g_vms(g_i_)
          enddo
          vms = vms + d3_v
C--------
C
10        continue
99999   continue
C
        d2_b = 1.0d0 / nno
        do g_i_ = 1, g_p_
          g_vms(g_i_) = d2_b * g_vms(g_i_)
        enddo
        vms = vms / nno
C--------
C
C.... DISTRIBUTE OUT TO THE NODES
C
        do 99998 n = 1, numnod
          do g_i_ = 1, g_p_
            g_strain(g_i_, 1, 7, n) = g_vms(g_i_)
          enddo
          strain(1, 7, n) = vms
C--------
20        continue
99998   continue
C
        return
      end
