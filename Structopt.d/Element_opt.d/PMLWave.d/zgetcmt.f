      subroutine ZGETCMT (rip, e, nu, epsx, epsy, c)
C
C     rip = 0 : plane stress
C     rip = 1 : plane strain
C
      implicit none
      double precision  rip, e, nu
      complex*16 epsx, epsy
      complex*16 c(4,4)
      double precision c11
      double precision  lambda, mu, lambdabar
      lambda = E*nu/((1+nu)*(1-2*nu))
      mu     = E/(2*(1+nu))
      c11 = lambda + 2.0d0*mu
C
C     Order:
C     (ux,x, uy,y, ux,y, uy,x)
c     This is how it looks for the usual homogeneous material
c$$$      c(1,1) =  e / (1.0d0-nu*nu)
c$$$      c(1,2) =  c(1,1)*nu
c$$$      c(1,3) =  0.0d0
c$$$      c(1,4) =  0.0d0
c$$$      c(2,1) =  c(1,2)
c$$$      c(2,2) =  c(1,1)
c$$$      c(2,3) =  0.0d0
c$$$      c(2,4) =  0.0d0
c$$$      c(3,1) =  0.0d0
c$$$      c(3,2) =  0.0d0
c$$$      c(3,3) =  0.5d0*c(1,1)*(1.0d0-nu)
c$$$      c(3,4) =  c(3,3)
c$$$      c(4,1) =  0.0d0
c$$$      c(4,2) =  0.0d0
c$$$      c(4,3) =  c(3,3)
c$$$      c(4,4) =  c(3,3)
C
C     *old PMLs ("stretching" coefficient applied directly
C                to plane stress constitutive matrix)
c$$$c     Now the anisotropic PML
c$$$      c11 = e / (1.0d0-nu*nu)
c$$$      c33 = 0.5d0*c11*(1.0d0-nu)
c$$$      c(1,1) = c11*epsy/epsx
c$$$      c(1,2) = c11*nu
c$$$      c(1,3) = 0.0d0
c$$$      c(1,4) = 0.0d0
c$$$      c(2,1) = c11*nu
c$$$      c(2,2) = c11*epsx/epsy
c$$$      c(2,3) = 0.0d0
c$$$      c(2,4) = 0.0d0
c$$$      c(3,1) = 0.0d0
c$$$      c(3,2) = 0.0d0
c$$$      c(3,3) = c33
c$$$      c(3,4) = c33*epsy/epsx
c$$$      c(4,1) = 0.0d0
c$$$      c(4,2) = 0.0d0
c$$$      c(4,3) = c33*epsx/epsy
c$$$      c(4,4) = c33
C     plane strain
c$$$      c(1,1) = (lambda+2*mu)*epsy/epsx
c$$$      c(1,2) = lambda
c$$$      c(1,3) = 0.0d0
c$$$      c(1,4) = 0.0d0
c$$$      c(2,1) = lambda
c$$$      c(2,2) = (lambda+2*mu)*epsx/epsy
c$$$      c(2,3) = 0.0d0
c$$$      c(2,4) = 0.0d0
c$$$      c(3,1) = 0.0d0
c$$$      c(3,2) = 0.0d0
c$$$      c(3,3) = mu
c$$$      c(3,4) = mu*epsy/epsx
c$$$      c(4,1) = 0.0d0
c$$$      c(4,2) = 0.0d0
c$$$      c(4,3) = mu*epsx/epsy
c$$$      c(4,4) = mu
C
C     
      if (rip.eq.0) then
C     Now the symmetric PML (plane stress)
         lambdabar=2.0d0*lambda*mu/c11
         c(1,1) = (lambdabar+2.0d0*mu)*epsy/epsx
         c(1,2) = lambdabar
         c(1,3) = 0.0d0
         c(1,4) = 0.0d0
         c(2,1) = lambdabar
         c(2,2) = (lambdabar+2.0d0*mu)*epsx/epsy
         c(2,3) = 0.0d0
         c(2,4) = 0.0d0
         c(3,1) = 0.0d0
         c(3,2) = 0.0d0
         c(3,3) = mu/4.0d0*(epsx/epsy+epsy/epsx+2.0d0)
         c(3,4) = mu/4.0d0*(epsx/epsy-epsy/epsx)
         c(4,1) = 0.0d0
         c(4,2) = 0.0d0
         c(4,3) = mu/4.0d0*(epsx/epsy-epsy/epsx)
         c(4,4) = mu/4.0d0*(epsx/epsy+epsy/epsx-2.0d0)
         return
      endif
      if (rip.eq.1) then
C
c     Now the symmetric PML (plane strain)
         c(1,1) = c11*epsy/epsx
         c(1,2) = lambda
         c(1,3) = 0.0d0
         c(1,4) = 0.0d0
         c(2,1) = lambda
         c(2,2) = c11*epsx/epsy
         c(2,3) = 0.0d0
         c(2,4) = 0.0d0
         c(3,1) = 0.0d0
         c(3,2) = 0.0d0
         c(3,3) = mu/4.0d0*(epsx/epsy+epsy/epsx+2.0d0)
         c(3,4) = mu/4.0d0*(epsx/epsy-epsy/epsx)
         c(4,1) = 0.0d0
         c(4,2) = 0.0d0
         c(4,3) = mu/4.0d0*(epsx/epsy-epsy/epsx)
         c(4,4) = mu/4.0d0*(epsx/epsy+epsy/epsx-2.0d0)
         return
      endif
C
C     error output
      write(*,*) " *** ERROR: wrong option in zgetcmt (PMLs)"
      stop
C
      end
