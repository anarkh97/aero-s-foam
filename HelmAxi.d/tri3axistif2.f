C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     This routine forms the stiffness element matrix of a
C     three-node triangle for 1 dof per node coming from
C     the derivation in theta of a Fourier mode. So, we compute
C                      Ni(r,z) Nj(r,z) / r drdz
C     No test is made to see if r=0, because r>=0 and the Gauss
C     points are always inside the quadrilateral.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL   TRI3AXISTIF2 (R, Z, P, SM, M)
C
C     where the input arguments are
C
C       R         (3 x 1) array of r coordinates of triangle nodes
C       Z         (3 x 1) array of z coordinates of triangle nodes
C       P         Identifies quadrature rule by value and sign
C                 (see TRIGGAUSSQ)
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        (3 x 3) computed element matrix
C       STATUS    Status character variable.  Blank if no error detected.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine   tri3axistif2(r, z, p, sm, m)
C
C                   A R G U M E N T S
C
      integer           p, m
      double precision  r(*), z(*)
      double precision  sm(m,*)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  q(3), qr(3), qz(3)
      double precision  zeta1, zeta2, zeta3, det, w, weight
      double precision  c1r
      double precision  rk 
      integer           i, ir, j, jr, k
      integer           ls(3)
      data              ls /1,2,3/
C
C                   L O G I C
C
      do 1200  j = 1,3
        do 1100  i = 1,3
          sm(i,j) = 0.0
 1100     continue
 1200   continue
C
      do 3000  k = 1, abs(p)
C
        call TRIGGAUSSQ(p, k, zeta1, zeta2, zeta3, weight)
        call TRIG3SHAPE(' ', zeta1, zeta2, zeta3, r, z, q, qr, qz, det)
C
        if (det .le. 0.0)        then
          write(6,*) 'Negative Jacobian determinant'
          if (det .eq. 0.0)      then
            write(6,*) 'Zero Jacobian determinant'
          end if
          return
        end if
C
        w = weight * (0.5*det) * (q(1)+q(2)+q(3))
C
        rk = zeta1*r(1) + zeta2*r(2) + zeta3*r(3)
        w = w / rk
C
        do 2000  j = 1,3
          jr =    ls(j)
          c1r =  q(j) * w
          do 1500  i = j,3
            ir =     ls(i)
            sm(ir,jr) =  sm(ir,jr) + q(i)*c1r
            sm(jr,ir) =  sm(ir,jr)
 1500       continue
 2000     continue
 3000   continue
C
      return
      end
