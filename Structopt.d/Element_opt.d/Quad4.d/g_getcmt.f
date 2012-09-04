      subroutine gxgetcmt(rip, e, de, nu, dnu, c, dc)
c
c     rip = 0 : plane stress
c     rip = 1 : plane strain
c
      implicit none

      double precision  rip
      double precision  e, nu, c(3,3)
      double precision  de, dnu, dc(3,3)

      double precision  omn,om2n,opn
      double precision  domn,dom2n,dopn

      if (rip.eq.0) then
         c(1,1) =  e / (1.0d0-nu*nu)
         c(2,2) =  c(1,1)
         c(3,3) =  0.5d0*c(1,1)*(1.0d0-nu)
         c(1,2) =  c(1,1)*nu
         c(2,1) =  c(1,2)
         c(1,3) =  0.0d0
         c(2,3) =  0.0d0
         c(3,1) =  0.0d0
         c(3,2) =  0.0d0

         dc(1,1) = (de + 2.0d0*e*nu*dnu/(1.0d0-nu*nu))/(1.0d0-nu*nu)
         dc(2,2) =  dc(1,1)
         dc(3,3) =  0.5d0*dc(1,1)*(1.0d0-nu) - 0.5d0*c(1,1)*dnu
         dc(1,2) =  dc(1,1)*nu + c(1,1)*dnu
         dc(2,1) =  dc(1,2)
         dc(1,3) =  0.0d0
         dc(2,3) =  0.0d0
         dc(3,1) =  0.0d0
         dc(3,2) =  0.0d0
         return
      endif

      if (rip.eq.1.0) then
      
         omn     = 1.0d0-nu
         om2n    = 1.0d0-2.0d0*nu
         opn     = 1.0d0+nu

         domn    = -dnu
         dom2n   = -2.0d0*dnu
         dopn    = dnu

         c(1,1) =  e*omn/opn/om2n
         c(2,2) =  c(1,1)
         c(3,3) =  e/opn/2.0d0
         c(1,2) =  e*nu/opn/om2n
         c(2,1) =  c(1,2)
         c(1,3) =  0.0d0
         c(2,3) =  0.0d0
         c(3,1) =  0.0d0
         c(3,2) =  0.0d0

         dc(1,1) = (de*omn+e*domn-e*omn*(dopn*om2n+opn*dom2n)
     &              /opn/om2n)/opn/om2n
         dc(2,2) =  dc(1,1)
         dc(3,3) =  (de-e*dopn/opn)/opn/2.0d0
         dc(1,2) =  (de*nu+e*dnu-e*nu*(dopn*om2n+opn*dom2n)
     &              /opn/om2n)/opn/om2n
         dc(2,1) =  dc(1,2)
         dc(1,3) =  0.0d0
         dc(2,3) =  0.0d0
         dc(3,1) =  0.0d0
         dc(3,2) =  0.0d0
	 
         return
      endif

      write(*,*) " *** ERROR: wrong option in getcmt"
      stop

      end
