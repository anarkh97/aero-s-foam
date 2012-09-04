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
C
C     QUAD4SHAPE computes the value of the shape functions for a
C     four-noded isoparametric quadrilateral and its
C     x-y derivatives, at a sample  point given by its quadrilateral
C     coordinates (xi,eta)
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       call    q4shpe (xi, eta, x, y, s, sx, sy, det)
C
C     Input arguments:
C
C       XI, ETA   Quadrilateral coordinates of given point
C       X         (4 x 1) array of x coordinates of quadrilateral corners
C       Y         (4 x 1) array of y coordinates of quadrilateral corners
C
C     Outputs arguments:
C
C       S         (4 x 1) array of shape function values
C       SX        (4 x 1) array of shape function x-derivatives
C       SY        (4 x 1) array of shape function y-derivatives
C       DET        Value of Jacobian determinant
C
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine gxq4shpe(xi, eta, x, g_x, y, g_y, s,
     * sx, g_sx, sy, g_sy, det, g_det)
C
C                   A R G U M E N T S
C
        double precision xi, eta, x(*), y(*), s(*), sx(*), sy(*), det
C
C                   L O C A L   V A R I A B L E S
C
        integer i
        double precision d1, d2, d3, d4, d1h, d2h, d3h, d4h, cdet
        double precision s1(4), s2(4), xd1, yd1, xd2, yd2
C
C                   L O G I C
C
C
C.... COMPUTE THE SHAPE FUNCTIONS
C
        integer g_pmax_
        parameter (g_pmax_ = 1)
        integer g_i_, g_p_, ldg_x, ldg_y, ldg_det, ldg_sx, ldg_sy
        parameter (g_p_= 1, ldg_x= 1, ldg_y= 1, ldg_det= 1, ldg_sx= 1,
     &             ldg_sy= 1)
        double precision d2_b, d8_b, d7_b, d6_b, d2_v, d7_v, d6_v, g_xd1
     *(g_pmax_), g_x(ldg_x, *), g_yd1(g_pmax_)
        double precision g_y(ldg_y, *), g_xd2(g_pmax_), g_yd2(g_pmax_), 
     *g_det(ldg_det), g_cdet(g_pmax_), g_sx(ldg_sx, *), g_sy(ldg_sy, *)
        save g_xd1, g_yd1, g_xd2, g_yd2, g_cdet
        intrinsic dble
        integer g_ehfid
        data g_ehfid /0/
C
        call ehsfid(g_ehfid, 'q4shpe','g_q4shpe.f')
C
        if (g_p_ .gt. g_pmax_) then
          print *, 'Parameter g_p_ is greater than g_pmax_'
          stop
        endif
        d1 = 0.5 * (1.0 + xi)
        d2 = 0.5 * (1.0 + eta)
        d3 = 1.0 - d1
        d4 = 1.0 - d2
        s(1) = d3 * d4
        s(2) = d4 * d1
        s(3) = d1 * d2
        s(4) = d2 * d3
C
C.... COMPUTE THE SHAPE FUNCTION DERIVATIVES
C
        d1h = 0.5 * d1
        d2h = 0.5 * d2
        d3h = 0.5 * d3
        d4h = 0.5 * d4
        s1(1) = -d4h
        s1(2) = d4h
        s1(3) = d2h
        s1(4) = -d2h
        s2(1) = -d3h
        s2(2) = -d1h
        s2(3) = d1h
        s2(4) = d3h
        do g_i_ = 1, g_p_
          g_xd1(g_i_) = d2h * g_x(g_i_, 3) + (-d2h) * g_x(g_i_, 4) + (-d
     *4h) * g_x(g_i_, 1) + d4h * g_x(g_i_, 2)
        enddo
        xd1 = (x(2) - x(1)) * d4h - (x(4) - x(3)) * d2h
C--------
        do g_i_ = 1, g_p_
          g_yd1(g_i_) = d2h * g_y(g_i_, 3) + (-d2h) * g_y(g_i_, 4) + (-d
     *4h) * g_y(g_i_, 1) + d4h * g_y(g_i_, 2)
        enddo
        yd1 = (y(2) - y(1)) * d4h - (y(4) - y(3)) * d2h
C--------
        do g_i_ = 1, g_p_
          g_xd2(g_i_) = d3h * g_x(g_i_, 4) + (-d3h) * g_x(g_i_, 1) + (-d
     *1h) * g_x(g_i_, 2) + d1h * g_x(g_i_, 3)
        enddo
        xd2 = (x(3) - x(2)) * d1h - (x(1) - x(4)) * d3h
C--------
        do g_i_ = 1, g_p_
          g_yd2(g_i_) = d3h * g_y(g_i_, 4) + (-d3h) * g_y(g_i_, 1) + (-d
     *1h) * g_y(g_i_, 2) + d1h * g_y(g_i_, 3)
        enddo
        yd2 = (y(3) - y(2)) * d1h - (y(1) - y(4)) * d3h
C--------
C
C.... COMPUTE THE DETERMINANT OF THE JACOBIAN
C
        do g_i_ = 1, g_p_
          g_det(g_i_) = -yd1 * g_xd2(g_i_) + (-xd2) * g_yd1(g_i_) + xd1 
     ** g_yd2(g_i_) + yd2 * g_xd1(g_i_)
        enddo
        det = xd1 * yd2 - yd1 * xd2
C--------
        if (det .eq. 0.0) then
          return
        endif
        d2_v = 1.0d0 / det
        d2_b = -d2_v / det
        do g_i_ = 1, g_p_
          g_cdet(g_i_) = d2_b * g_det(g_i_)
        enddo
        cdet = d2_v
C--------
        do 99999 i = 1, 4
          d6_v = yd2 * s1(i) - yd1 * s2(i)
          d6_b = -cdet * s2(i)
          d7_b = cdet * s1(i)
          do g_i_ = 1, g_p_
            g_sx(g_i_, i) = d6_b * g_yd1(g_i_) + d7_b * g_yd2(g_i_) + d6
     *_v * g_cdet(g_i_)
          enddo
          sx(i) = cdet * d6_v
C--------
          d7_v = -xd2 * s1(i) + xd1 * s2(i)
          d6_b = cdet * s2(i)
          d8_b = -(cdet * s1(i))
          do g_i_ = 1, g_p_
            g_sy(g_i_, i) = d6_b * g_xd1(g_i_) + d8_b * g_xd2(g_i_) + d7
     *_v * g_cdet(g_i_)
          enddo
          sy(i) = cdet * d7_v
C--------
2000      continue
99999   continue
        return
      end
C=END FORTRAN
