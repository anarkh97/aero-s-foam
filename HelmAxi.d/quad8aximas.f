C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     This routine forms the element consistent mass matrix of a
C     eight-node quadrilateral for 1 dof per node described
C     in cylindrical coordinates. So, we compute
C                       Ni(r,z) Nj(r,z) rdrdz
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     The calling sequence is
C
C       CALL   QUAD8AXIMAS (R, Z, P, MM, M)
C
C     where the input arguments are
C
C       R         (8 x 1) array of x coordinates of quadrilateral nodes
C       Z         (8 x 1) array of y coordinates of quadrilateral nodes
C       P         Gauss quadrature rule (no. of points)
C       M         First dimension of MM in calling program.
C
C     The outputs are:
C
C       MM        (8 x 8) computed element mass matrix.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine quad8aximas(r, z, p, mm, m )
C
C                   A R G U M E N T S
C
      integer           p, m
      double precision  r(*), z(*)
      double precision  mm(m,*)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  q(8), qr(8), qz(8)
      double precision  xi, eta, det, w, weight
      double precision  c1r
      double precision  rk
      integer           i, ir, j, jr, k, l
      integer           ls(8)
C
C                   D A T A
C
      data    ls /1,2,3,4,5,6,7,8/
C
C                   L O G I C
C
      do 1200  j = 1,8
        do 1100  i = 1,8
          mm(i,j) = 0.0
 1100     continue
 1200   continue
C
      do 3000  k = 1,p
        do 2500  l = 1,p
C
          call     QGAUSS (p, k, p, l, xi, eta, weight)
          call     QUAD8SHAPE (' ', xi, eta, r, z, q, qr, qz, det)
C
          if (det .le. 0.0)        then
            write(6,*) 'Negative Jacobian determinant'
            if (det .eq. 0.0)      then
              write(6,*) 'Zero Jacobian determinant'
            end if
            return
          end if
C
          w = weight * det * (q(1)+q(2)+q(3)+q(4)+q(5)+q(6)+q(7)+q(8))
C
          rk = 0.0
          do i = 1,8
            rk = rk + r(i)*q(i) 
          end do
C
          w  = w * rk
C
          do 2000  j = 1,8
            jr =    ls(j)
            c1r =  q(j) * w
            do 1500  i = j,8
              ir =     ls(i)
              mm(ir,jr) =  mm(ir,jr) + q(i)*c1r 
              mm(jr,ir) =  mm(ir,jr)
 1500         continue
 2000       continue
 2500     continue
 3000   continue
C
      return
      end
