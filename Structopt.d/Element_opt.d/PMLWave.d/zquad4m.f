C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     ZQUAD4MSTIF forms the element stiffness matrix of a
C     four-node quadrilateral in plane stress.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL   ZQUAD4M ( X, Y, H, C, P, SM, M)
C
C     where the input arguments are
C
C       X         (4 x 1) array of x coordinates of quadrilateral nodes
C       Y         (4 x 1) array of y coordinates of quadrilateral nodes
C       H         (4 x 1) array of thicknesses at quadrilateral nodes
C       C         (4 x 4) constitutive material matrix (not
C                         integrated thorugh the thickness)
C       P         Gauss quadrature rule (no. of points)
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        (8 x 8) computed element stiffness matrix.  The
C                 arrangement of rows and columns pertains to node
C                 displacements arranged in the order
C                  (vx1, vy1, vx2, ... vy4)
C                 This particular ordering is obtained by setting array
C                 LS to 1,3,5,7,2,4,6,8  (see DATA statement below)
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    zquad4m(x, y, h, c, p, sm, m)
      implicit none
C
C                   A R G U M E N T S
C      
      integer           p, m
      double precision  x(*), y(*), h(*)
      complex*16        c(4,*)
      complex*16        sm(m,*)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  q(4), qx(4), qy(4)
      double precision  xi, eta, det, weight
      complex*16        B(4,8), CB(4,8)
      complex*16        w
      integer           i, j, k, l
C
C                   L O G I C
C
      do 1200  j = 1,8
         do 1100  i = 1,8
            sm(i,j) = 0.0d0
 1100    continue
 1200 continue
      do 1300 j=1,8
         do 1400 i=1,4
            B(i,j) = 0.0d0
 1400    continue
 1300 continue
C      
      do 3000  k = 1,p
         do 2500  l = 1,p
            call     QGAUSS (p, k, p, l, xi, eta, weight)
            call     Q4SHPE (xi, eta, x, y, q, qx, qy, det)
C                 
            if (det .le. 0.0)        then
               write(6,*) 'Negative Jacobian determinant in 4 node quad'
               if (det .eq. 0.0)      then
                  write(6,*) 'Zero Jacobian determinant in 4 node quad'
               end if
               stop 
            end if            
            w =    weight * det *
     $           (h(1)*q(1)+h(2)*q(2)+h(3)*q(3)+h(4)*q(4))
C
            B(1,1) = qx(1)
            B(1,3) = qx(2)
            B(1,5) = qx(3)
            B(1,7) = qx(4)
            B(2,2) = qy(1)
            B(2,4) = qy(2)
            B(2,6) = qy(3)
            B(2,8) = qy(4)
            B(3,1) = qy(1)
            B(3,2) = qx(1)
            B(3,3) = qy(2)
            B(3,4) = qx(2)
            B(3,5) = qy(3)
            B(3,6) = qx(3)
            B(3,7) = qy(4)
            B(3,8) = qx(4)
            B(4,1) = qy(1)
            B(4,2) = -qx(1)
            B(4,3) = qy(2)
            B(4,4) = -qx(2)
            B(4,5) = qy(3)
            B(4,6) = -qx(3)
            B(4,7) = qy(4)
            B(4,8) = -qx(4)
            call zsymm('l', 'u', 4, 8, (1.0d0, 0.0d0), C, 4, B, 4,
     $           (0.0d0, 0.0d0), CB, 4)
            call zgemm('T', 'N', 8, 8, 4,              w, B, 4, CB, 4, 
     $           (1.0d0, 0.0d0), sm, m)
 2500    continue
 3000 continue
      
      return
      end
      
