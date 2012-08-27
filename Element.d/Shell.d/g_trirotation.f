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
C Modifed: 3-31-98
C
C Authors: Michel Lesoinne and Kendall Pierson
C
C************************************************
C This subroutine rotates the local stiffness
C matrix to global (or computational) axis defined
C in each corner node
C************************************************
C
      subroutine g_trirotation(g_p_, sm, g_sm, ldg_sm, r1, g_r1, ldg_r1)
C
C Declarations
C
        real*8 r1(3, 3), sm(18, 18)
C
C Local Declarations
C
        real*8 prod1(18)
        integer i, j, k
C
C premultiplication by rotation matrix
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_r1, ldg_sm
        double precision g_prod1(g_pmax_, 18), g_r1(ldg_r1, 3, 3), g_sm(
     *ldg_sm, 18, 18)
        save g_prod1
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'trirotation','g_trirotation.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do 99996 j = 1, 18
          do 99998 k = 1, 3
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k) = 0.0d0
            enddo
            prod1(k) = 0.0d+00
C--------
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k + 3) = 0.0d0
            enddo
            prod1(k + 3) = 0.0d+00
C--------
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k + 6) = 0.0d0
            enddo
            prod1(k + 6) = 0.0d+00
C--------
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k + 9) = 0.0d0
            enddo
            prod1(k + 9) = 0.0d+00
C--------
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k + 12) = 0.0d0
            enddo
            prod1(k + 12) = 0.0d+00
C--------
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k + 15) = 0.0d0
            enddo
            prod1(k + 15) = 0.0d+00
C--------
            do 99999 i = 1, 3
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k) = r1(k, i) * g_sm(g_i_, i, j) + sm(i, j
     *) * g_r1(g_i_, k, i) + g_prod1(g_i_, k)
              enddo
              prod1(k) = prod1(k) + r1(k, i) * sm(i, j)
C--------
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k + 3) = r1(k, i) * g_sm(g_i_, i + 3, j) +
     * sm(i + 3, j) * g_r1(g_i_, k, i) + g_prod1(g_i_, k + 3)
              enddo
              prod1(k + 3) = prod1(k + 3) + r1(k, i) * sm(i + 3, j)
C--------
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k + 6) = r1(k, i) * g_sm(g_i_, i + 6, j) +
     * sm(i + 6, j) * g_r1(g_i_, k, i) + g_prod1(g_i_, k + 6)
              enddo
              prod1(k + 6) = prod1(k + 6) + r1(k, i) * sm(i + 6, j)
C--------
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k + 9) = r1(k, i) * g_sm(g_i_, i + 9, j) +
     * sm(i + 9, j) * g_r1(g_i_, k, i) + g_prod1(g_i_, k + 9)
              enddo
              prod1(k + 9) = prod1(k + 9) + r1(k, i) * sm(i + 9, j)
C--------
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k + 12) = r1(k, i) * g_sm(g_i_, i + 12, j)
     * + sm(i + 12, j) * g_r1(g_i_, k, i) + g_prod1(g_i_, k + 12)
              enddo
              prod1(k + 12) = prod1(k + 12) + r1(k, i) * sm(i + 12, j)
C--------
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k + 15) = r1(k, i) * g_sm(g_i_, i + 15, j)
     * + sm(i + 15, j) * g_r1(g_i_, k, i) + g_prod1(g_i_, k + 15)
              enddo
              prod1(k + 15) = prod1(k + 15) + r1(k, i) * sm(i + 15, j)
C--------
20            continue
99999       continue
99998     continue
C
          do 99997 k = 1, 18
            do g_i_ = 1, g_p_
              g_sm(g_i_, k, j) = g_prod1(g_i_, k)
            enddo
            sm(k, j) = prod1(k)
C--------
30          continue
99997     continue
10        continue
99996   continue
C
C postmultiplication by rotation matrix transposed
C
        do 99992 j = 1, 18
          do 99994 k = 1, 3
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k) = 0.0d0
            enddo
            prod1(k) = 0.0d+00
C--------
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k + 3) = 0.0d0
            enddo
            prod1(k + 3) = 0.0d+00
C--------
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k + 6) = 0.0d0
            enddo
            prod1(k + 6) = 0.0d+00
C--------
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k + 9) = 0.0d0
            enddo
            prod1(k + 9) = 0.0d+00
C--------
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k + 12) = 0.0d0
            enddo
            prod1(k + 12) = 0.0d+00
C--------
            do g_i_ = 1, g_p_
              g_prod1(g_i_, k + 15) = 0.0d0
            enddo
            prod1(k + 15) = 0.0d+00
C--------
            do 99995 i = 1, 3
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k) = sm(j, i) * g_r1(g_i_, k, i) + r1(k, i
     *) * g_sm(g_i_, j, i) + g_prod1(g_i_, k)
              enddo
              prod1(k) = prod1(k) + sm(j, i) * r1(k, i)
C--------
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k + 3) = sm(j, i + 3) * g_r1(g_i_, k, i) +
     * r1(k, i) * g_sm(g_i_, j, i + 3) + g_prod1(g_i_, k + 3)
              enddo
              prod1(k + 3) = prod1(k + 3) + sm(j, i + 3) * r1(k, i)
C--------
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k + 6) = sm(j, i + 6) * g_r1(g_i_, k, i) +
     * r1(k, i) * g_sm(g_i_, j, i + 6) + g_prod1(g_i_, k + 6)
              enddo
              prod1(k + 6) = prod1(k + 6) + sm(j, i + 6) * r1(k, i)
C--------
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k + 9) = sm(j, i + 9) * g_r1(g_i_, k, i) +
     * r1(k, i) * g_sm(g_i_, j, i + 9) + g_prod1(g_i_, k + 9)
              enddo
              prod1(k + 9) = prod1(k + 9) + sm(j, i + 9) * r1(k, i)
C--------
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k + 12) = sm(j, i + 12) * g_r1(g_i_, k, i)
     * + r1(k, i) * g_sm(g_i_, j, i + 12) + g_prod1(g_i_, k + 12)
              enddo
              prod1(k + 12) = prod1(k + 12) + sm(j, i + 12) * r1(k, i)
C--------
              do g_i_ = 1, g_p_
                g_prod1(g_i_, k + 15) = sm(j, i + 15) * g_r1(g_i_, k, i)
     * + r1(k, i) * g_sm(g_i_, j, i + 15) + g_prod1(g_i_, k + 15)
              enddo
              prod1(k + 15) = prod1(k + 15) + sm(j, i + 15) * r1(k, i)
C--------
50            continue
99995       continue
99994     continue
C
          do 99993 k = 1, 18
            do g_i_ = 1, g_p_
              g_sm(g_i_, j, k) = g_prod1(g_i_, k)
            enddo
            sm(j, k) = prod1(k)
C--------
60          continue
99993     continue
40        continue
99992   continue
C
        return
      end
