      subroutine zbrkcmt(e, nu, epsx, epsy, epsz, c)
      implicit none
      real*8       e, nu, lambda, mu
      complex*16   epsx, epsy, epsz, c11
      complex*16   c(9,9)

C     Order:
C     (ux,x, uy,y, uz,z, ux,y, ux,z, uy,x, uy,z, uz,x, uz,y)
c     This is how it looks for the usual homogeneous material
      integer      i, j
      do 20  i=1,9
         do 10  j=1,9
            c(i,j) = 0.0d0
 10      continue
 20   continue
      lambda = e*nu/((1+nu)*(1-2*nu))
      mu = e/(2*(1+nu))
      c11 = lambda + 2.0d0*mu
      c(1,1) =  c11*epsy*epsz/epsx
      c(2,2) =  c11*epsz*epsx/epsy
      c(3,3) =  c11*epsx*epsy/epsz
      c(1,2) =  lambda*epsz
      c(2,1) =  c(1,2)
      c(1,3) =  lambda*epsy
      c(3,1) =  c(1,3)
      c(2,3) =  lambda*epsx
      c(3,2) =  c(2,3)
      c(4,4) =mu*(epsz/2.0d0+epsy*epsz/4.0d0/epsx+epsx*epsz/4.0d0/epsy)
      c(4,7) =  mu/4.0d0*(epsy*epsz/epsx-epsx*epsz/epsy)
      c(7,4) =  c(4,7)
      c(5,5) =mu*(epsy/2.0d0+epsy*epsz/4.0d0/epsx+epsx*epsy/4.0d0/epsz)
      c(5,8) =  mu/4.0d0*(epsy*epsz/epsx-epsx*epsy/epsz)
      c(8,5) =  c(5,8)
      c(6,6) =mu*(epsx/2.0d0+epsx*epsz/4.0d0/epsy+epsx*epsy/4.0d0/epsz)
      c(6,9) =  mu/4.0d0*(epsx*epsz/epsy-epsx*epsy/epsz)
      c(9,6) =  c(6,9)
      c(7,7) =mu*(-epsz/2.0d0+epsy*epsz/4.0d0/epsx+epsx*epsz/4.0d0/epsy)
      c(8,8) =mu*(-epsy/2.0d0+epsy*epsz/4.0d0/epsx+epsx*epsy/4.0d0/epsz)
      c(9,9) =mu*(-epsx/2.0d0+epsx*epsz/4.0d0/epsy+epsx*epsy/4.0d0/epsz)
      return
      end
C=END FORTRAN
