C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     This routine forms the element stiffness matrix of a
C     four-node quadrilateral for 1 dof per node described
C     in cylindrical coordinates. So, we compute 
C               ( Ni,r * Nj,r + Ni,z * Nj,z ) rdrdz
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL   QUADAXISTIF1 ( R, Z, P, SM, M )
C
C     where the input arguments are
C
C       R         (4 x 1) array of r coordinates of quadrilateral nodes
C       Z         (4 x 1) array of z coordinates of quadrilateral nodes
C       P         Gauss quadrature rule (no. of points)
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        (4 x 4) computed element stiffness matrix.
C                 As there is only one dof per node, we set  
C                 LS to 1,2,3,4  (see DATA statement below)
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    quadaxistif1(r, z, p, sm, m )
C
C                   A R G U M E N T S
C
      integer           p, m
      double precision  r(*), z(*)
      double precision  sm(m,*)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  q(4), qr(4), qz(4)
      double precision  xi, eta, det, w, weight
      double precision  rk, alpha1, alpha2, alpha3, alpha4
      double precision  c1r, c3z
      integer           i, ir, j, jr, k, l
      integer           ls(4)
C
C                   D A T A
C
      data              ls /1,2,3,4/
C
C                   L O G I C
C
      do 1200  j = 1,4
        do 1100  i = 1,4
          sm(i,j) = 0.0
 1100     continue
 1200   continue
C
      alpha1 = (r(1)+r(2)+r(3)+r(4))*0.25
      alpha2 = (r(3)+r(4)-r(2)-r(1))*0.25
      alpha3 = (r(2)+r(3)-r(1)-r(4))*0.25
      alpha4 = (r(2)+r(4)-r(1)-r(3))*0.25
C
      do 3000  k = 1,p
        do 2500  l = 1,p
C
          call     QGAUSS (p, k, p, l, xi, eta, weight)
          call     Q4SHPE (xi, eta, r, z, q, qr, qz, det)
C
          if (det .le. 0.0)        then
            write(6,*)  'Negative Jacobian determinant'
            if (det .eq. 0.0)      then
              write(6,*) 'Zero Jacobian determinant'
            end if
            stop 
          end if
C
          w =    weight * det *
     $          (q(1)+q(2)+q(3)+q(4))
C
          rk = alpha1 + eta*alpha2 + xi*(alpha3 - eta*alpha4)
C
          w  = w * rk
C
          do 2000  j = 1,4
            jr =   ls(j)
            c1r =  qr(j) * w
            c3z =  qz(j) * w
            do 1500  i = j,4
              ir =     ls(i)
              sm(ir,jr) =  sm(ir,jr) + qr(i)*c1r + qz(i)*c3z
              sm(jr,ir) =  sm(ir,jr)
 1500         continue
 2000       continue
 2500     continue
 3000   continue
C
      return
      end
