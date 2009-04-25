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
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C MODIFIED BY PAUL STERN MARCH 7 1990
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     QUAD2D forms the coupling  matrix of a
C     four-node quadrilateral. 
C
C
C     The calling sequence is
C
C       CALL   QUAD2D ( X, Y, C, P, SM, M )
C
C     where the input arguments are
C
C       X         (4 x 1) array of x coordinates of quadrilateral nodes
C       Y         (4 x 1) array of y coordinates of quadrilateral nodes
C       C         element coupling coefficient
C       P         Gauss quadrature rule (no. of points)
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        (8 x 4) computed element coupling  matrix.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine gxquad2d(x, g_x, y, g_y, c, g_c, p, cm, g_cm, m)
C
C                   A R G U M E N T S
C
        integer p, m
        real*8 x(*), y(*)
        real*8 c
        real*8 cm(8, *)
C
C                   L O C A L   V A R I A B L E S
C
        real*8 q(4), qx(4), qy(4)
        real*8 xi, eta, det, weight, w
        integer i, j, k, l
C
C                   L O G I C
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_cm, ldg_c, ldg_x, ldg_y
        parameter (g_p_= 1, ldg_cm= 1, ldg_c= 1, ldg_x= 1, ldg_y= 1)
        double precision d8_b, d7_b, d6_b, d5_v, g_cm(ldg_cm, m, *), g_w
     *(g_pmax_), g_det(g_pmax_), g_c(ldg_c), g_qx(g_pmax_, 4), g_qy(g_pm
     *ax_, 4)
        double precision g_x(ldg_x, 1), g_y(ldg_y, 1)
        save g_w, g_det, g_qx, g_qy
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
c        call ehsfid(g_ehfid, 'quad2d','g_quad2d.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        do 99998 j = 1, 4
          do 99999 i = 1, 8
            do g_i_ = 1, g_p_
              g_cm(g_i_, i, j) = 0.0d0
            enddo
            cm(i, j) = 0.0d0
C--------
200         continue
99999     continue
100       continue
99998   continue
C
C
C.... COMPUTE THE COUPLING ELEMENT MATRIX            
C
        do 99996 k = 1, p
          do 99997 l = 1, p
C
C.... COMPUTE THE SHAPE FUNCTIONS & DERIVATIVES
C.... NOTE THE THERMAL AND MECHANICS USE THE SAME SFS
C
            call qgauss(p, k, p, l, xi, eta, weight)
            call gxq4shpe(xi, eta, x, g_x, y, g_y, q,
     * qx, g_qx, qy, g_qy, det, g_det)
C
C.... CHECK THE DETERMINANT TO SEE IF THE QUAD IS TO DISTORTED
C
            if (det .le. 0.0) then
              write (6, *) 'Negative Jacobian determinant'
              if (det .eq. 0.0) then
                write (6, *) 'Zero Jacobian determinant'
              endif
              stop
            endif
            do g_i_ = 1, g_p_
              g_w(g_i_) = weight * g_det(g_i_)
            enddo
            w = weight * det
C--------
C
C.... ASSEMBLE THE ELEMENT COUPLING MATRIX
C
            d5_v = c * qx(1) * q(1)
            d6_b = w * q(1)
            d7_b = d6_b * qx(1)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 1, 1) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 1)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 1, 1)
            enddo
            cm(1, 1) = cm(1, 1) + d5_v * w
C--------
            d5_v = c * qx(1) * q(2)
            d6_b = w * q(2)
            d7_b = d6_b * qx(1)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 1, 2) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 1)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 1, 2)
            enddo
            cm(1, 2) = cm(1, 2) + d5_v * w
C--------
            d5_v = c * qx(1) * q(3)
            d6_b = w * q(3)
            d7_b = d6_b * qx(1)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 1, 3) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 1)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 1, 3)
            enddo
            cm(1, 3) = cm(1, 3) + d5_v * w
C--------
            d5_v = c * qx(1) * q(4)
            d6_b = w * q(4)
            d7_b = d6_b * qx(1)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 1, 4) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 1)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 1, 4)
            enddo
            cm(1, 4) = cm(1, 4) + d5_v * w
C--------
C
            d5_v = c * qy(1) * q(1)
            d6_b = w * q(1)
            d7_b = d6_b * qy(1)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 2, 1) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 1)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 2, 1)
            enddo
            cm(2, 1) = cm(2, 1) + d5_v * w
C--------
            d5_v = c * qy(1) * q(2)
            d6_b = w * q(2)
            d7_b = d6_b * qy(1)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 2, 2) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 1)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 2, 2)
            enddo
            cm(2, 2) = cm(2, 2) + d5_v * w
C--------
            d5_v = c * qy(1) * q(3)
            d6_b = w * q(3)
            d7_b = d6_b * qy(1)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 2, 3) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 1)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 2, 3)
            enddo
            cm(2, 3) = cm(2, 3) + d5_v * w
C--------
            d5_v = c * qy(1) * q(4)
            d6_b = w * q(4)
            d7_b = d6_b * qy(1)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 2, 4) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 1)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 2, 4)
            enddo
            cm(2, 4) = cm(2, 4) + d5_v * w
C--------
C
            d5_v = c * qx(2) * q(1)
            d6_b = w * q(1)
            d7_b = d6_b * qx(2)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 3, 1) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 2)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 3, 1)
            enddo
            cm(3, 1) = cm(3, 1) + d5_v * w
C--------
            d5_v = c * qx(2) * q(2)
            d6_b = w * q(2)
            d7_b = d6_b * qx(2)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 3, 2) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 2)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 3, 2)
            enddo
            cm(3, 2) = cm(3, 2) + d5_v * w
C--------
            d5_v = c * qx(2) * q(3)
            d6_b = w * q(3)
            d7_b = d6_b * qx(2)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 3, 3) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 2)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 3, 3)
            enddo
            cm(3, 3) = cm(3, 3) + d5_v * w
C--------
            d5_v = c * qx(2) * q(4)
            d6_b = w * q(4)
            d7_b = d6_b * qx(2)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 3, 4) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 2)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 3, 4)
            enddo
            cm(3, 4) = cm(3, 4) + d5_v * w
C--------
C
            d5_v = c * qy(2) * q(1)
            d6_b = w * q(1)
            d7_b = d6_b * qy(2)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 4, 1) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 2)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 4, 1)
            enddo
            cm(4, 1) = cm(4, 1) + d5_v * w
C--------
            d5_v = c * qy(2) * q(2)
            d6_b = w * q(2)
            d7_b = d6_b * qy(2)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 4, 2) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 2)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 4, 2)
            enddo
            cm(4, 2) = cm(4, 2) + d5_v * w
C--------
            d5_v = c * qy(2) * q(3)
            d6_b = w * q(3)
            d7_b = d6_b * qy(2)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 4, 3) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 2)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 4, 3)
            enddo
            cm(4, 3) = cm(4, 3) + d5_v * w
C--------
            d5_v = c * qy(2) * q(4)
            d6_b = w * q(4)
            d7_b = d6_b * qy(2)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 4, 4) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 2)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 4, 4)
            enddo
            cm(4, 4) = cm(4, 4) + d5_v * w
C--------
C
            d5_v = c * qx(3) * q(1)
            d6_b = w * q(1)
            d7_b = d6_b * qx(3)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 5, 1) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 3)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 5, 1)
            enddo
            cm(5, 1) = cm(5, 1) + d5_v * w
C--------
            d5_v = c * qx(3) * q(2)
            d6_b = w * q(2)
            d7_b = d6_b * qx(3)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 5, 2) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 3)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 5, 2)
            enddo
            cm(5, 2) = cm(5, 2) + d5_v * w
C--------
            d5_v = c * qx(3) * q(3)
            d6_b = w * q(3)
            d7_b = d6_b * qx(3)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 5, 3) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 3)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 5, 3)
            enddo
            cm(5, 3) = cm(5, 3) + d5_v * w
C--------
            d5_v = c * qx(3) * q(4)
            d6_b = w * q(4)
            d7_b = d6_b * qx(3)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 5, 4) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 3)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 5, 4)
            enddo
            cm(5, 4) = cm(5, 4) + d5_v * w
C--------
C
            d5_v = c * qy(3) * q(1)
            d6_b = w * q(1)
            d7_b = d6_b * qy(3)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 6, 1) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 3)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 6, 1)
            enddo
            cm(6, 1) = cm(6, 1) + d5_v * w
C--------
            d5_v = c * qy(3) * q(2)
            d6_b = w * q(2)
            d7_b = d6_b * qy(3)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 6, 2) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 3)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 6, 2)
            enddo
            cm(6, 2) = cm(6, 2) + d5_v * w
C--------
            d5_v = c * qy(3) * q(3)
            d6_b = w * q(3)
            d7_b = d6_b * qy(3)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 6, 3) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 3)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 6, 3)
            enddo
            cm(6, 3) = cm(6, 3) + d5_v * w
C--------
            d5_v = c * qy(3) * q(4)
            d6_b = w * q(4)
            d7_b = d6_b * qy(3)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 6, 4) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 3)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 6, 4)
            enddo
            cm(6, 4) = cm(6, 4) + d5_v * w
C--------
C
            d5_v = c * qx(4) * q(1)
            d6_b = w * q(1)
            d7_b = d6_b * qx(4)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 7, 1) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 4)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 7, 1)
            enddo
            cm(7, 1) = cm(7, 1) + d5_v * w
C--------
            d5_v = c * qx(4) * q(2)
            d6_b = w * q(2)
            d7_b = d6_b * qx(4)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 7, 2) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 4)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 7, 2)
            enddo
            cm(7, 2) = cm(7, 2) + d5_v * w
C--------
            d5_v = c * qx(4) * q(3)
            d6_b = w * q(3)
            d7_b = d6_b * qx(4)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 7, 3) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 4)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 7, 3)
            enddo
            cm(7, 3) = cm(7, 3) + d5_v * w
C--------
            d5_v = c * qx(4) * q(4)
            d6_b = w * q(4)
            d7_b = d6_b * qx(4)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 7, 4) = d5_v * g_w(g_i_) + d8_b * g_qx(g_i_, 4)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 7, 4)
            enddo
            cm(7, 4) = cm(7, 4) + d5_v * w
C--------
C
            d5_v = c * qy(4) * q(1)
            d6_b = w * q(1)
            d7_b = d6_b * qy(4)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 8, 1) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 4)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 8, 1)
            enddo
            cm(8, 1) = cm(8, 1) + d5_v * w
C--------
            d5_v = c * qy(4) * q(2)
            d6_b = w * q(2)
            d7_b = d6_b * qy(4)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 8, 2) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 4)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 8, 2)
            enddo
            cm(8, 2) = cm(8, 2) + d5_v * w
C--------
            d5_v = c * qy(4) * q(3)
            d6_b = w * q(3)
            d7_b = d6_b * qy(4)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 8, 3) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 4)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 8, 3)
            enddo
            cm(8, 3) = cm(8, 3) + d5_v * w
C--------
            d5_v = c * qy(4) * q(4)
            d6_b = w * q(4)
            d7_b = d6_b * qy(4)
            d8_b = d6_b * c
            do g_i_ = 1, g_p_
              g_cm(g_i_, 8, 4) = d5_v * g_w(g_i_) + d8_b * g_qy(g_i_, 4)
     * + d7_b * g_c(g_i_) + g_cm(g_i_, 8, 4)
            enddo
            cm(8, 4) = cm(8, 4) + d5_v * w
C--------
C
110         continue
99997     continue
210       continue
99996   continue
C
        return
      end
