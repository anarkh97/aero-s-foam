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
      subroutine g_compfrot(g_p_, v1o, g_v1o, ldg_v1o, v2o, g_v2o, ldg_v
     *2o, v3o, g_v3o, ldg_v3o, v1n, v2n, v3n, rot, g_rot, ldg_rot)
C=====================================================================C
C                                                                     C
C     This Subroutine Forms the Rotation Matrix From Local Elemental  C
C     Frames to Global (Computational) Axes such that the Rotation    C
C     Matrix is:                                                      C
C                                                                     C
C           [ [rot]   0     0   ]                                     C
C     [R] = [   0   [rot]   0   ]   and:   [Ke] = [R] * [ke] * [R]^T  C
C           [   0     0   [rot] ]                                     C
C                                                                     C
C     where [Ke] and [ke] Represent the (Same) Elemental Stiffness    C
C     Matrix Expressed in the Global and Local Frame Systems, Resp.   C
C     (See Routine "compmrot.f".)                                     C
C                                                                     C
C     The Rotation Provides the Transformation Matrix [rot] for:      C
C                                                                     C
C     [v] = [rot] * [x]                                               C
C                                                                     C
C     where:                                                          C
C                                                                     C
C     [x] are the Coordinates in the Old System                       C
C     [v] are the Coordinates in the New System                       C
C                                                                     C
C     The Input [v1o], [v2o] and [v3o] Represent the Old Basis and    C
C     the Input [v1n], [v2n] and [v3n] Represent the New Basis Both   C
C     in Cartesian Coordinates.                                       C
C                                                                     C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
        real*8 v1o(3), v2o(3), v3o(3)
        real*8 v1n(3), v2n(3), v3n(3)
        real*8 rot(6, 6)
C
C.....Local Variables
C
        integer i, j
        real*8 zero
C
C     ----
C     DATA
C     ----
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_rot, ldg_v1o, ldg_v2o, ldg_v3o
        double precision g_rot(ldg_rot, 6, 6), g_v1o(ldg_v1o, 3), g_v2o(
     *ldg_v2o, 3), g_v3o(ldg_v3o, 3)
        data zero /0.000000d+00/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE OUTPUT ROTATION MATRIX
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
1002        continue
99999     continue
1001      continue
99998   continue
C
C.....COMPUTE THE ROTATION MATRIX FOR TRANSLATIONAL DOFS
C
        do 99997 i = 1, 3
          do g_i_ = 1, g_p_
            g_rot(g_i_, 1, 1) = v1n(i) * g_v1o(g_i_, i) + g_rot(g_i_, 1,
     * 1)
          enddo
          rot(1, 1) = rot(1, 1) + v1n(i) * v1o(i)
C--------
          do g_i_ = 1, g_p_
            g_rot(g_i_, 1, 2) = v1n(i) * g_v2o(g_i_, i) + g_rot(g_i_, 1,
     * 2)
          enddo
          rot(1, 2) = rot(1, 2) + v1n(i) * v2o(i)
C--------
          do g_i_ = 1, g_p_
            g_rot(g_i_, 1, 3) = v1n(i) * g_v3o(g_i_, i) + g_rot(g_i_, 1,
     * 3)
          enddo
          rot(1, 3) = rot(1, 3) + v1n(i) * v3o(i)
C--------
          do g_i_ = 1, g_p_
            g_rot(g_i_, 2, 1) = v2n(i) * g_v1o(g_i_, i) + g_rot(g_i_, 2,
     * 1)
          enddo
          rot(2, 1) = rot(2, 1) + v2n(i) * v1o(i)
C--------
          do g_i_ = 1, g_p_
            g_rot(g_i_, 2, 2) = v2n(i) * g_v2o(g_i_, i) + g_rot(g_i_, 2,
     * 2)
          enddo
          rot(2, 2) = rot(2, 2) + v2n(i) * v2o(i)
C--------
          do g_i_ = 1, g_p_
            g_rot(g_i_, 2, 3) = v2n(i) * g_v3o(g_i_, i) + g_rot(g_i_, 2,
     * 3)
          enddo
          rot(2, 3) = rot(2, 3) + v2n(i) * v3o(i)
C--------
          do g_i_ = 1, g_p_
            g_rot(g_i_, 3, 1) = v3n(i) * g_v1o(g_i_, i) + g_rot(g_i_, 3,
     * 1)
          enddo
          rot(3, 1) = rot(3, 1) + v3n(i) * v1o(i)
C--------
          do g_i_ = 1, g_p_
            g_rot(g_i_, 3, 2) = v3n(i) * g_v2o(g_i_, i) + g_rot(g_i_, 3,
     * 2)
          enddo
          rot(3, 2) = rot(3, 2) + v3n(i) * v2o(i)
C--------
          do g_i_ = 1, g_p_
            g_rot(g_i_, 3, 3) = v3n(i) * g_v3o(g_i_, i) + g_rot(g_i_, 3,
     * 3)
          enddo
          rot(3, 3) = rot(3, 3) + v3n(i) * v3o(i)
C--------
2001      continue
99997   continue
C
C.....FILL OUT THE ROTATION MATRIX FOR ROTATIONAL DOFS
C
        do 99995 j = 1, 3
          do 99996 i = 1, 3
            do g_i_ = 1, g_p_
              g_rot(g_i_, i + 3, j + 3) = g_rot(g_i_, i, j)
            enddo
            rot(i + 3, j + 3) = rot(i, j)
C--------
3002        continue
99996     continue
3001      continue
99995   continue
C
C     ------
C     RETURN
C     ------
C
        return
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPFROT"
C=========================C
