C                           DISCLAIMER
C
C   This file was generated on 10/14/98 by the version of
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
C
C Authors: Michel Lesoinne and Kendall Pierson
C
C Modified routine: 3-31-98
C
C This routine provides the transformation
C matrix  R for:
C
C         vp= R x
C
C  x  are the coordinates in the old system
C  vp are the coordinates in the new system
C
C  Input:
C       v1o,v2o,v3o old basis described in cartesian
C                   components
C       v1n,v2n,v3n new basis
C
      subroutine g_rotation(v1o, g_v1o, v2o, g_v2o, 
     *v3o, g_v3o, v1n, v2n, v3n, r, g_r)
C
C Declarations
C
        real*8 v1o(3), v2o(3), v3o(3)
        real*8 v1n(3), v2n(3), v3n(3), r(3, 3)
C
C Local Declarations
C
        integer i, j
C
C Zero rotation matrix
C
        double precision g_r(3, 3), g_v1o(3), g_v2o(3)
     *, g_v3o(3)
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'rotation','g_rotation.f')
C
        do i = 1, 3
          do j = 1, 3
              g_r(i, j) = 0.0d0
            r(i, j) = 0.0d+00
C--------
          enddo
        enddo
C
C Compute rotation matrix
C
        do i = 1, 3
            g_r(1, 1) = v1n(i) * g_v1o(i) + g_r(1, 1)
          r(1, 1) = r(1, 1) + v1n(i) * v1o(i)
C--------
            g_r(1, 2) = v1n(i) * g_v2o(i) + g_r(1, 2)
          r(1, 2) = r(1, 2) + v1n(i) * v2o(i)
C--------
            g_r(1, 3) = v1n(i) * g_v3o(i) + g_r(1, 3)
          r(1, 3) = r(1, 3) + v1n(i) * v3o(i)
C--------
            g_r(2, 1) = v2n(i) * g_v1o(i) + g_r(2, 1)
          r(2, 1) = r(2, 1) + v2n(i) * v1o(i)
C--------
            g_r(2, 2) = v2n(i) * g_v2o(i) + g_r(2, 2)
          r(2, 2) = r(2, 2) + v2n(i) * v2o(i)
C--------
            g_r(2, 3) = v2n(i) * g_v3o(i) + g_r(2, 3)
          r(2, 3) = r(2, 3) + v2n(i) * v3o(i)
C--------
            g_r(3, 1) = v3n(i) * g_v1o(i) + g_r(3, 1)
          r(3, 1) = r(3, 1) + v3n(i) * v1o(i)
C--------
            g_r(3, 2) = v3n(i) * g_v2o(i) + g_r(3, 2)
          r(3, 2) = r(3, 2) + v3n(i) * v2o(i)
C--------
            g_r(3, 3) = v3n(i) * g_v3o(i) + g_r(3, 3)
          r(3, 3) = r(3, 3) + v3n(i) * v3o(i)
C--------
        enddo
C
        return
      end
