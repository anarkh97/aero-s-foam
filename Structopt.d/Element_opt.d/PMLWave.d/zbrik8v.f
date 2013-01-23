C=PURPOSE Form stiffness of 8-node hexahedron in 3-D  stress
C
C     HEXA8STIF forms the element stiffness matrix of a
C     eight-node hexahedron in 3-D stress.
C
C
C     where the input arguments are
C
C       OPT       Option letter argument, presently ignored.
C       X         (8 x 1) array of x coordinates of hexahedron nodes
C       Y         (8 x 1) array of y coordinates of hexahedron nodes
C       Z         (8 x 1) array of z coordinates of hexahedron nodes
C       C         (6 x 6) constitutive material matrix 
C       P         Gauss quadrature rule (no. of points)
C       JREF      Jacobian det from initial stiffness computation of this elem
C
C     The outputs are:
C
C       SM        (24 x 24) computed element stiffness matrix.  The
C                 arrangement of rows and columns pertains to node
C                 displacements arranged in the order
C                  (vx1, vy1, vz1, vx2, ... vz8)
C                 This particular ordering is obtained by setting array
C                 LS to 1,4,7,10,13,16,19,22,2,5,8,11,14,
C                 17,20,23,3,6,9,12,15,18,21,24
C                 ( see DATA statement below)
C
C       STATUS    Status integer variable.  Zero if no error
C                 detected.
C
      subroutine    zbrik8v(x, y, z, c, p, sm, status)
C
C                   A R G U M E N T S
C
      implicit none
      integer p,status
      real*8  x(8), y(8), z(8)
      complex*16 c(9,9)
      complex*16 sm(24,24)
C
C                   L O C A L   V A R I A B L E S
C
      complex*16  w
      complex*16  B(9,24), CB(9,24)
      real*8  q(8), qx(8), qy(8),qz(8)
      real*8  xi, eta, emu, det, weight
      integer           i, j, k, l
      integer           jj
C
C                   L O G I C
C
C
C     Initialize no error
      status = 0
C
      do 1200  j = 1,24
         do 1100  i = 1,24
            sm(i,j) = 0.0d0
 1100    continue
 1200 continue
      do 1300 j=1,24
         do 1400 i=1,9
            B(i,j) = 0.0d0
 1400    continue
 1300 continue
C     
      do 3000  k = 1,p
         do 2500  l = 1,p
            do 2400  jj = 1,p
               call hxgaus(p,k,p,l,p,jj,xi,eta,emu,weight)
               call h8shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
C     
               if (det .le. 0.0) then
                  write(6,*) 'Negative Jacobian determinant in brik8v.f'
C     
C     det = -det 
C     
                  status = -1
                  return
               end if
C     
               w = weight * det
               B( 1, 1) = qx(1)
               B( 1, 4) = qx(2)
               B( 1, 7) = qx(3)
               B( 1,10) = qx(4)
               B( 1,13) = qx(5)
               B( 1,16) = qx(6)
               B( 1,19) = qx(7)
               B( 1,22) = qx(8)
               B( 2, 2) = qy(1)
               B( 2, 5) = qy(2)
               B( 2, 8) = qy(3)
               B( 2,11) = qy(4)
               B( 2,14) = qy(5)
               B( 2,17) = qy(6)
               B( 2,20) = qy(7)
               B( 2,23) = qy(8)
               B( 3, 3) = qz(1)
               B( 3, 6) = qz(2)
               B( 3, 9) = qz(3)
               B( 3,12) = qz(4)
               B( 3,15) = qz(5)
               B( 3,18) = qz(6)
               B( 3,21) = qz(7)
               B( 3,24) = qz(8)
               B( 4, 1) = qy(1)
               B( 4, 2) = qx(1)
               B( 4, 4) = qy(2)
               B( 4, 5) = qx(2)
               B( 4, 7) = qy(3)
               B( 4, 8) = qx(3)
               B( 4,10) = qy(4)
               B( 4,11) = qx(4)
               B( 4,13) = qy(5)
               B( 4,14) = qx(5)
               B( 4,16) = qy(6)
               B( 4,17) = qx(6)
               B( 4,19) = qy(7)
               B( 4,20) = qx(7)
               B( 4,22) = qy(8)
               B( 4,23) = qx(8)
               B( 5, 1) = qz(1)
               B( 5, 3) = qx(1)
               B( 5, 4) = qz(2)
               B( 5, 6) = qx(2)
               B( 5, 7) = qz(3)
               B( 5, 9) = qx(3)
               B( 5,10) = qz(4)
               B( 5,12) = qx(4)
               B( 5,13) = qz(5)
               B( 5,15) = qx(5)
               B( 5,16) = qz(6)
               B( 5,18) = qx(6)
               B( 5,19) = qz(7)
               B( 5,21) = qx(7)
               B( 5,22) = qz(8)
               B( 5,24) = qx(8)
               B( 6, 2) = qz(1)
               B( 6, 3) = qy(1)
               B( 6, 5) = qz(2)
               B( 6, 6) = qy(2)
               B( 6, 8) = qz(3)
               B( 6, 9) = qy(3)
               B( 6,11) = qz(4)
               B( 6,12) = qy(4)
               B( 6,14) = qz(5)
               B( 6,15) = qy(5)
               B( 6,17) = qz(6)
               B( 6,18) = qy(6)
               B( 6,20) = qz(7)
               B( 6,21) = qy(7)
               B( 6,23) = qz(8)
               B( 6,24) = qy(8)
               B( 7, 1) = qy(1)
               B( 7, 2) = -qx(1)
               B( 7, 4) = qy(2)
               B( 7, 5) = -qx(2)
               B( 7, 7) = qy(3)
               B( 7, 8) = -qx(3)
               B( 7,10) = qy(4)
               B( 7,11) = -qx(4)
               B( 7,13) = qy(5)
               B( 7,14) = -qx(5)
               B( 7,16) = qy(6)
               B( 7,17) = -qx(6)
               B( 7,19) = qy(7)
               B( 7,20) = -qx(7)
               B( 7,22) = qy(8)
               B( 7,23) = -qx(8)
               B( 8, 1) = qz(1)
               B( 8, 3) = -qx(1)
               B( 8, 4) = qz(2)
               B( 8, 6) = -qx(2)
               B( 8, 7) = qz(3)
               B( 8, 9) = -qx(3)
               B( 8,10) = qz(4)
               B( 8,12) = -qx(4)
               B( 8,13) = qz(5)
               B( 8,15) = -qx(5)
               B( 8,16) = qz(6)
               B( 8,18) = -qx(6)
               B( 8,19) = qz(7)
               B( 8,21) = -qx(7)
               B( 8,22) = qz(8)
               B( 8,24) = -qx(8)
               B( 9, 2) = qz(1)
               B( 9, 3) = -qy(1)
               B( 9, 5) = qz(2)
               B( 9, 6) = -qy(2)
               B( 9, 8) = qz(3)
               B( 9, 9) = -qy(3)
               B( 9,11) = qz(4)
               B( 9,12) = -qy(4)
               B( 9,14) = qz(5)
               B( 9,15) = -qy(5)
               B( 9,17) = qz(6)
               B( 9,18) = -qy(6)
               B( 9,20) = qz(7)
               B( 9,21) = -qy(7)
               B( 9,23) = qz(8)
               B( 9,24) = -qy(8)
               call zgemm('N', 'N',  9, 24, 9, (1.0d0, 0.0d0), C, 9, 
     $              B,  9,  (0.0d0, 0.0d0), CB, 9)
               call zgemm('T', 'N', 24, 24, 9,              w, B, 9, 
     $              CB, 9,  (1.0d0, 0.0d0), sm, 24)
 2400       continue
 2500    continue
 3000 continue
      
      return
      end
