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
      subroutine g_compbmm(g_p_, elm, type, x, g_x, ldg_x, y, g_y, ldg_y
     *, dm, g_dm, ldg_dm, alpha, f, rowm, colm, rot, g_rot, ldg_rot, p, 
     *g_p, ldg_p, kbmm, g_kbmm, ldg_kbmm)
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
        integer elm, type, rowm(9), colm(9)
        real*8 x(3), y(3), dm(3, 3), alpha, f
        real*8 kbmm(18, 18), rot(6, 6), p(9, 3)
C
C.....Local Variables
C
        integer i, j, row, col, dimp
        real*8 zero, three, six
        real*8 twicearea, factor
        real*8 x21, x32, x13, y21, y32, y13
        real*8 x12, x23, x31, y12, y23, y31
        real*8 dmpt1, dmpt2, dmpt3
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_kbmm, ldg_p, ldg_x, ldg_y, ldg_dm, ldg_r
     *ot
        double precision d1_p, d1_w, d9_b, d8_b, d3_b, d4_v, d7_b, d6_b,
     * d2_v, d3_v
        double precision d5_b, d4_b, g_kbmm(ldg_kbmm, 18, 18), g_p(ldg_p
     *, 9, 3), g_x21(g_pmax_), g_x(ldg_x, 3), g_x12(g_pmax_), g_x32(g_pm
     *ax_), g_x23(g_pmax_), g_x13(g_pmax_)
        double precision g_x31(g_pmax_), g_y21(g_pmax_), g_y(ldg_y, 3), 
     *g_y12(g_pmax_), g_y32(g_pmax_), g_y23(g_pmax_), g_y13(g_pmax_), g_
     *y31(g_pmax_), g_twicearea(g_pmax_), g_d1_w(g_pmax_)
        double precision g_factor(g_pmax_), g_dm(ldg_dm, 3, 3), g_dmpt1(
     *g_pmax_), g_dmpt2(g_pmax_), g_dmpt3(g_pmax_), g_rot(ldg_rot, 6, 6)
        save g_y13, g_y31, g_twicearea, g_d1_w, g_factor, g_dmpt1, g_dmp
     *t2, g_dmpt3
        save g_x21, g_x12, g_x32, g_x23, g_x13, g_x31, g_y21, g_y12, g_y
     *32, g_y23
        data zero /0.000000d+00/
        data three /3.000000d+00/
        data six /6.000000d+00/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE OUTPUT STIFFNESS MATRIX
C
        integer g_ehfid
        data g_ehfid /0/
C

C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do 99998 j = 1, 18
          do 99999 i = 1, 18
            do g_i_ = 1, g_p_
              g_kbmm(g_i_, i, j) = 0.0d0
            enddo
            kbmm(i, j) = zero
C--------
1002        continue
99999     continue
1001      continue
99998   continue
C
C.....CLEAR THE LOCAL MATRICES
C
        do 99996 j = 1, 3
          do 99997 i = 1, 9
            do g_i_ = 1, g_p_
              g_p(g_i_, i, j) = 0.0d0
            enddo
            p(i, j) = zero
C--------
1004        continue
99997     continue
1003      continue
99996   continue
C
C.....RETURN IF THE STIFFNESS FACTOR IS ZERO
C
        if (f .eq. zero) then
          return
        endif
C
C.....CHECK IF THE STIFFNESS FACTOR [F] IS POSITIVE
C
        if (f .lt. zero) then
          goto 100
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
C.....CHECK THE AREA (ERROR IF NOT POSITIVE)
C
        if (twicearea .le. zero) then
          goto 200
        endif
C
C.....ASSEMBLE THE MATRIX [P] W/ SHAPE FUNCTION DERIVATIVES
C
        do g_i_ = 1, g_p_
          g_p(g_i_, 1, 1) = g_y23(g_i_)
        enddo
        p(1, 1) = y23
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 2, 1) = 0.0d0
        enddo
        p(2, 1) = zero
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 3, 1) = g_y31(g_i_)
        enddo
        p(3, 1) = y31
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 4, 1) = 0.0d0
        enddo
        p(4, 1) = zero
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 5, 1) = g_y12(g_i_)
        enddo
        p(5, 1) = y12
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 6, 1) = 0.0d0
        enddo
        p(6, 1) = zero
C--------
C
        do g_i_ = 1, g_p_
          g_p(g_i_, 1, 2) = 0.0d0
        enddo
        p(1, 2) = zero
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 2, 2) = g_x32(g_i_)
        enddo
        p(2, 2) = x32
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 3, 2) = 0.0d0
        enddo
        p(3, 2) = zero
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 4, 2) = g_x13(g_i_)
        enddo
        p(4, 2) = x13
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 5, 2) = 0.0d0
        enddo
        p(5, 2) = zero
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 6, 2) = g_x21(g_i_)
        enddo
        p(6, 2) = x21
C--------
C
        do g_i_ = 1, g_p_
          g_p(g_i_, 1, 3) = g_x32(g_i_)
        enddo
        p(1, 3) = x32
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 2, 3) = g_y23(g_i_)
        enddo
        p(2, 3) = y23
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 3, 3) = g_x13(g_i_)
        enddo
        p(3, 3) = x13
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 4, 3) = g_y31(g_i_)
        enddo
        p(4, 3) = y31
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 5, 3) = g_x21(g_i_)
        enddo
        p(5, 3) = x21
C--------
        do g_i_ = 1, g_p_
          g_p(g_i_, 6, 3) = g_y12(g_i_)
        enddo
        p(6, 3) = y12
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
            g_p(g_i_, 7, 1) = -d5_b * g_y21(g_i_) + d5_b * g_y13(g_i_) +
     * d4_b * g_y23(g_i_)
          enddo
          p(7, 1) = y23 * d4_v * alpha / six
C--------
          d4_v = x31 - x12
          d3_b = 1.0d0 / six * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x32
          do g_i_ = 1, g_p_
            g_p(g_i_, 7, 2) = -d5_b * g_x12(g_i_) + d5_b * g_x31(g_i_) +
     * d4_b * g_x32(g_i_)
          enddo
          p(7, 2) = x32 * d4_v * alpha / six
C--------
          d3_b = 1.0d0 / three * alpha
          d6_b = -d3_b * y21
          d7_b = -d3_b * x12
          d8_b = d3_b * y13
          d9_b = d3_b * x31
          do g_i_ = 1, g_p_
            g_p(g_i_, 7, 3) = d7_b * g_y21(g_i_) + d6_b * g_x12(g_i_) + 
     *d9_b * g_y13(g_i_) + d8_b * g_x31(g_i_)
          enddo
          p(7, 3) = (x31 * y13 - x12 * y21) * alpha / three
C--------
C
          d4_v = y21 - y32
          d3_b = 1.0d0 / six * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * y31
          do g_i_ = 1, g_p_
            g_p(g_i_, 8, 1) = -d5_b * g_y32(g_i_) + d5_b * g_y21(g_i_) +
     * d4_b * g_y31(g_i_)
          enddo
          p(8, 1) = y31 * d4_v * alpha / six
C--------
          d4_v = x12 - x23
          d3_b = 1.0d0 / six * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x13
          do g_i_ = 1, g_p_
            g_p(g_i_, 8, 2) = -d5_b * g_x23(g_i_) + d5_b * g_x12(g_i_) +
     * d4_b * g_x13(g_i_)
          enddo
          p(8, 2) = x13 * d4_v * alpha / six
C--------
          d3_b = 1.0d0 / three * alpha
          d6_b = -d3_b * y32
          d7_b = -d3_b * x23
          d8_b = d3_b * y21
          d9_b = d3_b * x12
          do g_i_ = 1, g_p_
            g_p(g_i_, 8, 3) = d7_b * g_y32(g_i_) + d6_b * g_x23(g_i_) + 
     *d9_b * g_y21(g_i_) + d8_b * g_x12(g_i_)
          enddo
          p(8, 3) = (x12 * y21 - x23 * y32) * alpha / three
C--------
C
          d4_v = y32 - y13
          d3_b = 1.0d0 / six * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * y12
          do g_i_ = 1, g_p_
            g_p(g_i_, 9, 1) = -d5_b * g_y13(g_i_) + d5_b * g_y32(g_i_) +
     * d4_b * g_y12(g_i_)
          enddo
          p(9, 1) = y12 * d4_v * alpha / six
C--------
          d4_v = x23 - x31
          d3_b = 1.0d0 / six * alpha
          d4_b = d3_b * d4_v
          d5_b = d3_b * x21
          do g_i_ = 1, g_p_
            g_p(g_i_, 9, 2) = -d5_b * g_x31(g_i_) + d5_b * g_x23(g_i_) +
     * d4_b * g_x21(g_i_)
          enddo
          p(9, 2) = x21 * d4_v * alpha / six
C--------
          d3_b = 1.0d0 / three * alpha
          d6_b = -d3_b * y13
          d7_b = -d3_b * x31
          d8_b = d3_b * y32
          d9_b = d3_b * x23
          do g_i_ = 1, g_p_
            g_p(g_i_, 9, 3) = d7_b * g_y13(g_i_) + d6_b * g_x31(g_i_) + 
     *d9_b * g_y32(g_i_) + d8_b * g_x23(g_i_)
          enddo
          p(9, 3) = (x23 * y32 - x31 * y13) * alpha / three
C--------
C
          dimp = 9
C
        endif
C
C.....MULTIPLY MATRIX [P] BY THE SQUARE ROOT OF THE STIFFNESS FACTOR
C.....AND DIVIDE BY THE SQUARE ROOT OF FOUR TIMES THE AREA
C
        d2_v = 2.00d+00 * twicearea
        d3_v = f / d2_v
        d3_b = -d3_v / d2_v * 2.00d+00
        do g_i_ = 1, g_p_
          g_d1_w(g_i_) = d3_b * g_twicearea(g_i_)
        enddo
        d1_w = d3_v
        d2_v = sqrt(d1_w)
        if ( d1_w .gt. 0.0d0 ) then
           d1_p = 1.0d0 / (2.0d0 *  d2_v)
        else
           call ehufDV (9,d1_w, d2_v, d1_p,
     +'g_compbMM.f',
     +416)
        endif
        do g_i_ = 1, g_p_
          g_factor(g_i_) = d1_p * g_d1_w(g_i_)
        enddo
        factor = d2_v
C--------
C
        do 99994 j = 1, 3
          do 99995 i = 1, dimp
            do g_i_ = 1, g_p_
              g_p(g_i_, i, j) = factor * g_p(g_i_, i, j) + p(i, j) * g_f
     *actor(g_i_)
            enddo
            p(i, j) = factor * p(i, j)
C--------
2002        continue
99995     continue
2001      continue
99994   continue
C
C.....ASSEMBLE THE OUTPUT STIFFNESS SUCH THAT [kbMM] = [P]*[dm]*[P]^T
C
        do 99992 i = 1, dimp
C
          col = colm(i)
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dm(1, 2) * g_p(g_i_, i, 2) + p(i, 2) * g_dm(g
     *_i_, 1, 2) + dm(1, 1) * g_p(g_i_, i, 1) + p(i, 1) * g_dm(g_i_, 1, 
     *1)
          enddo
          d1_w = dm(1, 1) * p(i, 1) + dm(1, 2) * p(i, 2)
          do g_i_ = 1, g_p_
            g_dmpt1(g_i_) = dm(1, 3) * g_p(g_i_, i, 3) + p(i, 3) * g_dm(
     *g_i_, 1, 3) + g_d1_w(g_i_)
          enddo
          dmpt1 = d1_w + dm(1, 3) * p(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dm(2, 2) * g_p(g_i_, i, 2) + p(i, 2) * g_dm(g
     *_i_, 2, 2) + dm(2, 1) * g_p(g_i_, i, 1) + p(i, 1) * g_dm(g_i_, 2, 
     *1)
          enddo
          d1_w = dm(2, 1) * p(i, 1) + dm(2, 2) * p(i, 2)
          do g_i_ = 1, g_p_
            g_dmpt2(g_i_) = dm(2, 3) * g_p(g_i_, i, 3) + p(i, 3) * g_dm(
     *g_i_, 2, 3) + g_d1_w(g_i_)
          enddo
          dmpt2 = d1_w + dm(2, 3) * p(i, 3)
C--------
          do g_i_ = 1, g_p_
            g_d1_w(g_i_) = dm(3, 2) * g_p(g_i_, i, 2) + p(i, 2) * g_dm(g
     *_i_, 3, 2) + dm(3, 1) * g_p(g_i_, i, 1) + p(i, 1) * g_dm(g_i_, 3, 
     *1)
          enddo
          d1_w = dm(3, 1) * p(i, 1) + dm(3, 2) * p(i, 2)
          do g_i_ = 1, g_p_
            g_dmpt3(g_i_) = dm(3, 3) * g_p(g_i_, i, 3) + p(i, 3) * g_dm(
     *g_i_, 3, 3) + g_d1_w(g_i_)
          enddo
          dmpt3 = d1_w + dm(3, 3) * p(i, 3)
C--------
C
          do 99993 j = 1, i
            row = rowm(j)
            do g_i_ = 1, g_p_
              g_d1_w(g_i_) = dmpt2 * g_p(g_i_, j, 2) + p(j, 2) * g_dmpt2
     *(g_i_) + dmpt1 * g_p(g_i_, j, 1) + p(j, 1) * g_dmpt1(g_i_) + g_kbm
     *m(g_i_, row, col)
            enddo
            d1_w = kbmm(row, col) + dmpt1 * p(j, 1) + dmpt2 * p(j, 2)
            do g_i_ = 1, g_p_
              g_kbmm(g_i_, row, col) = dmpt3 * g_p(g_i_, j, 3) + p(j, 3)
     * * g_dmpt3(g_i_) + g_d1_w(g_i_)
            enddo
            kbmm(row, col) = d1_w + dmpt3 * p(j, 3)
C--------
            do g_i_ = 1, g_p_
              g_kbmm(g_i_, col, row) = g_kbmm(g_i_, row, col)
            enddo
            kbmm(col, row) = kbmm(row, col)
C--------
3002        continue
99993     continue
C
3001      continue
99992   continue
C
C.....OUTPUT THE MATRIX PRIOR TO ROTATION (FOR DEBUGING ONLY)
C
C     open(unit=90,file="kbMM.m")
C     write(90,*) "kbMM=["
C     do 991 i=1,18
C 991 write(90,9) (kbMM(i,j),j=1,18)
C     write(90,*) "     ];"
C     close(90)
C   9 format(18(1x,E16.9))
C
C.....ROTATE THE OUTPUT STIFFNESS MATRIX
C
C      call compmrot( kbMM , rot , rot , rot )
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
C.....ERROR-MESSAGE IF THE STIFFNESS FACTOR [F] IS NEGATIVE
C
100     continue
        write (*, *) '*** FATAL ERROR in routine COMPBMM       ***'
        write (*, *) '*** The Stiffness Factor [f] is Negative ***'
        write (*, *) '*** Check the Calling Sequence:          ***'
        write (*, *) '*** Factor [f] Must be Positive or Zero  ***'
        write (*, *) '*** EXECUTION TERNINATED RIGHT HERE      ***'
        stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
200     continue
        write (*, *) '*** FATAL ERROR in routine COMPBMM         ***'
        write (*, *) '*** The Triangle Area is Found Negative or ***'
        write (*, *) '*** Zero: Check the Nodal Point Numbering  ***'
        write (*, *) '*** ... Counterclockwise?                  ***'
        write (*, *) '*** EXECUTION TERNINATED RIGHT HERE        ***'
        stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPBMM"
C========================C
