C=PURPOSE Compute the gradient of the integral of (Sigma_vm / Sigma_bar) of 8-node hexa
C
C     input arguments are
C
C       X         (8 x 1) array of x coordinates of hexahedron nodes
C       Y         (8 x 1) array of y coordinates of hexahedron nodes
C       Z         (8 x 1) array of z coordinates of hexahedron nodes
C       C         (6 x 6) constitutive material matrix 
C       P         Gauss quadrature rule (no. of points)
C       v         (24x1) array of element node displacement arranged in
C                 ux1,uy1,uz1,ux2,uy2,uz2, ............., uz8
C       sigbar    given stress value
C
C     The outputs are:
C
C       VMINT     integral of (Von Mises stresse / stress_bar) 
C       VOL       volume of element
C       STATUS    Status integer variable.  Zero if no error
C                 detected.
C
      subroutine  gxbrvmint(x,y,z,dx,dy,dz,v,dv,c,dc,p,vmint,dvmint,
     *                      vol,dvol,sigbar,fac,iarea,status)
C
C                   A R G U M E N T S
C
C      implicit none
      integer p, status
      real*8  x(8), y(8), z(8), c(6,6), v(*)
      real*8  dx(8), dy(8), dz(8), dc(6,6), dv(*)
      real*8  vmint, dvmint, sigbar, fac, vol, dvol
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(8), qx(8), qy(8),qz(8)
      real*8  dqx(8), dqy(8),dqz(8), ddet, dw
      real*8  xi, eta, emu, det, w, weight
      real*8  epsxx, epsyy,epszz, gamxy, gamyz, gamxz
      real*8  depsxx, depsyy,depszz, dgamxy, dgamyz, dgamxz
      real*8  stress(6),dsxx, dsyy, dszz, dsxy, dsxz, dsyz, j2, derj2
      real*8  dstress(6),ddsxx, ddsyy, ddszz, ddsxy, ddsxz, ddsyz, dj2
      real*8  dfac, dervmint, a, b, cvar, pow, comp, dcomp
      integer k, l, jj, i
C
C
C                   L O G I C
C
C                             inline function
C
      pow(a,b) = exp(b*log(a))

C     check if any variation
      cvar = 0.0
      
      do 5000 i = 1, 8
        cvar = cvar + dx(i)**2 + dy(i)**2 +dz(i)**2 
 5000 continue
C
C                             Initialize 
      status = 0
      vmint  = 0.0d0
      vol    = 0.0d0
      dvmint = 0.0d0
      dvol   = 0.0d0
C
c                                                          cvar.eq.0       
      if (cvar.eq.0) then

      do 6000  k = 1,p
        do 5500  l = 1,p
          do 5400  jj = 1,p
            call hxgaus(p,k,p,l,p,jj,xi,eta,emu,weight)
            call h8shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
C
            if (det .le. 0.0d0) then
              write(6,*) 'Negative Jacobian determinant in gxbrintvm.f'
              status = -1
              return
            end if
C
            w  = weight * det 
C
            epsxx = qx(1)*v(1)  + qx(2)*v(4)  + qx(3)*v(7)  +
     $       	    qx(4)*v(10) + qx(5)*v(13) + qx(6)*v(16) + 
     $      	    qx(7)*v(19) + qx(8)*v(22)

            epsyy = qy(1)*v(2)  + qy(2)*v(5)  + qy(3)*v(8)  +
     $       	    qy(4)*v(11) + qy(5)*v(14) + qy(6)*v(17) + 
     $      	    qy(7)*v(20) + qy(8)*v(23)

            epszz = qz(1)*v(3)  + qz(2)*v(6)  + qz(3)*v(9)  + 
     $       	    qz(4)*v(12) + qz(5)*v(15) + qz(6)*v(18) + 
     $      	    qz(7)*v(21) + qz(8)*v(24)
            gamxy = qy(1)*v(1)  + qy(2)*v(4)  + qy(3)*v(7)  + 
     $       	    qy(4)*v(10) + qy(5)*v(13) + qy(6)*v(16) + 
     $       	    qy(7)*v(19) + qy(8)*v(22) + qx(1)*v(2)  + 
     $      	    qx(2)*v(5)  + qx(3)*v(8)  + qx(4)*v(11) +
     $       	    qx(5)*v(14) + qx(6)*v(17) + qx(7)*v(20) + 
     $      	    qx(8)*v(23)
            gamyz = qy(1)*v(3)  + qy(2)*v(6)  + qy(3)*v(9)  + 
     $       	    qy(4)*v(12) + qy(5)*v(15) + qy(6)*v(18) +
     $       	    qy(7)*v(21) + qy(8)*v(24) + qz(1)*v(2)  + 
     $      	    qz(2)*v(5)  + qz(3)*v(8)  + qz(4)*v(11) +
     $       	    qz(5)*v(14) + qz(6)*v(17) + qz(7)*v(20) + 
     $      	    qz(8)*v(23)
            gamxz = qz(1)*v(1)  + qz(2)*v(4)  + qz(3)*v(7)  +
     $       	    qz(4)*v(10) + qz(5)*v(13) + qz(6)*v(16) +
     $       	    qz(7)*v(19) + qz(8)*v(22) + qx(1)*v(3)  + 
     $      	    qx(2)*v(6)  + qx(3)*v(9)  + qx(4)*v(12) +
     $       	    qx(5)*v(15) + qx(6)*v(18) + qx(7)*v(21) + 
     $      	    qx(8)*v(24)
C
            depsxx =  qx(1)*dv(1) + qx(2)*dv(4)+ qx(3)*dv(7)  +
     $       	          qx(4)*dv(10) + qx(5)*dv(13) + qx(6)*dv(16) + 
     $      	          qx(7)*dv(19) + qx(8)*dv(22)
            depsyy =  qy(1)*dv(2) + qy(2)*dv(5)+ qy(3)*dv(8)  +
     $       	          qy(4)*dv(11) + qy(5)*dv(14) + qy(6)*dv(17) + 
     $      	          qy(7)*dv(20) + qy(8)*dv(23)
            depszz =  qz(1)*dv(3) + qz(2)*dv(6)+ qz(3)*dv(9)  + 
     $       	          qz(4)*dv(12) + qz(5)*dv(15) + qz(6)*dv(18) + 
     $      	          qz(7)*dv(21) + qz(8)*dv(24)
            dgamxy =  qy(1)*dv(1) + qy(2)*dv(4)+ qy(3)*dv(7)  + 
     $       	          qy(4)*dv(10) + qy(5)*dv(13) + qy(6)*dv(16) + 
     $       	          qy(7)*dv(19) + qy(8)*dv(22) + qx(1)*dv(2)  + 
     $      	          qx(2)*dv(5)  + qx(3)*dv(8)  + qx(4)*dv(11) +
     $       	          qx(5)*dv(14) + qx(6)*dv(17) + qx(7)*dv(20) + 
     $      	          qx(8)*dv(23)
            dgamyz =  qy(1)*dv(3) + qy(2)*dv(6)+ qy(3)*dv(9)  + 
     $       	          qy(4)*dv(12) + qy(5)*dv(15) + qy(6)*dv(18) +
     $       	          qy(7)*dv(21) + qy(8)*dv(24) + qz(1)*dv(2)  + 
     $      	          qz(2)*dv(5)  + qz(3)*dv(8)  + qz(4)*dv(11) +
     $       	          qz(5)*dv(14) + qz(6)*dv(17) + qz(7)*dv(20) + 
     $      	          qz(8)*dv(23)
            dgamxz =  qz(1)*dv(1) + qz(2)*dv(4)+ qz(3)*dv(7)  +
     $       	          qz(4)*dv(10) + qz(5)*dv(13) + qz(6)*dv(16) +
     $       	          qz(7)*dv(19) + qz(8)*dv(22) + qx(1)*dv(3)  + 
     $      	          qx(2)*dv(6)  + qx(3)*dv(9)  + qx(4)*dv(12) +
     $       	          qx(5)*dv(15) + qx(6)*dv(18) + qx(7)*dv(21) + 
     $      	          qx(8)*dv(24)
C
C.... ENGINEERING STRESS COMPUTATION
C
            stress(1) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*epszz
     $       	      + c(1,4)*gamxy + c(1,5)*gamyz + c(1,6)*gamxz
            stress(2) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*epszz
     $       	      + c(2,4)*gamxy + c(2,5)*gamyz + c(2,6)*gamxz
            stress(3) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*epszz
     $       	      + c(3,4)*gamxy + c(3,5)*gamyz + c(3,6)*gamxz
            stress(4) = c(4,1)*epsxx + c(4,2)*epsyy + c(4,3)*epszz
     $       	      + c(4,4)*gamxy + c(4,5)*gamyz + c(4,6)*gamxz
            stress(5) = c(5,1)*epsxx + c(5,2)*epsyy + c(5,3)*epszz
     $       	      + c(5,4)*gamxy + c(5,5)*gamyz + c(5,6)*gamxz
            stress(6) = c(6,1)*epsxx + c(6,2)*epsyy + c(6,3)*epszz
     $       	      + c(6,4)*gamxy + c(6,5)*gamyz + c(6,6)*gamxz

            dstress(1) = dc(1,1)*epsxx + dc(1,2)*epsyy + dc(1,3)*epszz
     $       	       + dc(1,4)*gamxy + dc(1,5)*gamyz + dc(1,6)*gamxz
            dstress(2) = dc(2,1)*epsxx + dc(2,2)*epsyy + dc(2,3)*epszz
     $       	       + dc(2,4)*gamxy + dc(2,5)*gamyz + dc(2,6)*gamxz
            dstress(3) = dc(3,1)*epsxx + dc(3,2)*epsyy + dc(3,3)*epszz
     $       	       + dc(3,4)*gamxy + dc(3,5)*gamyz + dc(3,6)*gamxz
            dstress(4) = dc(4,1)*epsxx + dc(4,2)*epsyy + dc(4,3)*epszz
     $       	       + dc(4,4)*gamxy + dc(4,5)*gamyz + dc(4,6)*gamxz
            dstress(5) = dc(5,1)*epsxx + dc(5,2)*epsyy + dc(5,3)*epszz
     $       	       + dc(5,4)*gamxy + dc(5,5)*gamyz + dc(5,6)*gamxz
            dstress(6) = dc(6,1)*epsxx + dc(6,2)*epsyy + dc(6,3)*epszz
     $       	       + dc(6,4)*gamxy + dc(6,5)*gamyz + dc(6,6)*gamxz

            dstress(1) = dstress(1)   + c(1,1)*depsxx + c(1,2)*depsyy+
     $       	        c(1,3)*depszz + c(1,4)*dgamxy + c(1,5)*dgamyz+ 
     $                  c(1,6)*dgamxz
            dstress(2) = dstress(2)  + c(2,1)*depsxx + c(2,2)*depsyy+ 
     $       	       c(2,3)*depszz + c(2,4)*dgamxy + c(2,5)*dgamyz+ 
     $                 c(2,6)*dgamxz
            dstress(3) = dstress(3)  + c(3,1)*depsxx + c(3,2)*depsyy+ 
     $       	       c(3,3)*depszz + c(3,4)*dgamxy + c(3,5)*dgamyz+ 
     $                 c(3,6)*dgamxz
            dstress(4) = dstress(4)  + c(4,1)*depsxx + c(4,2)*depsyy+ 
     $       	       c(4,3)*depszz + c(4,4)*dgamxy + c(4,5)*dgamyz+ 
     $                 c(4,6)*dgamxz
            dstress(5) = dstress(5)  + c(5,1)*depsxx + c(5,2)*depsyy+ 
     $       	       c(5,3)*depszz + c(5,4)*dgamxy + c(5,5)*dgamyz+ 
     $                 c(5,6)*dgamxz
            dstress(6) = dstress(6)  + c(6,1)*depsxx + c(6,2)*depsyy+ 
     $       	       c(6,3)*depszz + c(6,4)*dgamxy + c(6,5)*dgamyz+ 
     $                 c(6,6)*dgamxz
C
C.... COMPUTE THE FIRST DEVEATORIC STRESSES
C
	    comp = (stress(1) + stress(2) + stress(3))/3.0d0
	    dsxx = stress(1) - comp
            dsyy = stress(2) - comp
            dszz = stress(3) - comp
            dsxy = stress(4)
            dsyz = stress(5)
            dsxz = stress(6)

	    dcomp = (dstress(1) + dstress(2) + dstress(3))/3.0d0
	    ddsxx = dstress(1) - dcomp
            ddsyy = dstress(2) - dcomp
            ddszz = dstress(3) - dcomp
            ddsxy = dstress(4)
            ddsyz = dstress(5)
            ddsxz = dstress(6)
C
C.... COMPUTE THE SECOND DEVEATORIC STRESS
C
	     j2    = ((dsxx*dsxx)+(dsyy*dsyy)+(dszz*dszz))/2.0d0+
     &                (dsxy*dsxy)+(dsyz*dsyz)+(dsxz*dsxz)
             j2    = dsqrt(3.0d0*j2)

	     derj2 = (ddsxx*dsxx)+(ddsyy*dsyy)+(ddszz*dszz)+
     &              ((ddsxy*dsxy)+(ddsyz*dsyz)+(ddsxz*dsxz))*2.0d0

             if (j2.eq.0.0) then
               dj2 = 0.0D0
             else
	       dj2   = 1.5d0*derj2/j2
             endif
C
C.... COMPUTE THE VON MISES STRESS AND VOLUME
C
             if (iarea.eq.0) w  = 1.0d0

	     vmint   = vmint + w * pow((j2/sigbar),fac)

	     dfac    = fac - 1.0d0
             dervmint= w * fac * (dj2/sigbar) * pow((j2/sigbar),dfac) 
             dvmint  = dvmint + dervmint

	     vol     = vol  + w
C
 5400       continue
 5500     continue
 6000   continue

C                                                          cvar.ne.0 
      else

      do 3000  k = 1,p
        do 2500  l = 1,p
          do 2400  jj = 1,p
            call hxgaus(p,k,p,l,p,jj,xi,eta,emu,weight)
            call gxh8shpe(xi,eta,emu,x,y,z,dx,dy,dz,q,qx,qy,qz,
     *                    dqx,dqy,dqz,det,ddet)
C
            if (det .le. 0.0d0) then
              write(6,*) 'Negative Jacobian determinant in gxbrintvm.f'
              status = -1
              return
            end if
C
            w  = weight * det 
            dw = weight * ddet
C
            epsxx = qx(1)*v(1)  + qx(2)*v(4)  + qx(3)*v(7)  +
     $       	    qx(4)*v(10) + qx(5)*v(13) + qx(6)*v(16) + 
     $      	    qx(7)*v(19) + qx(8)*v(22)

            epsyy = qy(1)*v(2)  + qy(2)*v(5)  + qy(3)*v(8)  +
     $       	    qy(4)*v(11) + qy(5)*v(14) + qy(6)*v(17) + 
     $      	    qy(7)*v(20) + qy(8)*v(23)

            epszz = qz(1)*v(3)  + qz(2)*v(6)  + qz(3)*v(9)  + 
     $       	    qz(4)*v(12) + qz(5)*v(15) + qz(6)*v(18) + 
     $      	    qz(7)*v(21) + qz(8)*v(24)
            gamxy = qy(1)*v(1)  + qy(2)*v(4)  + qy(3)*v(7)  + 
     $       	    qy(4)*v(10) + qy(5)*v(13) + qy(6)*v(16) + 
     $       	    qy(7)*v(19) + qy(8)*v(22) + qx(1)*v(2)  + 
     $      	    qx(2)*v(5)  + qx(3)*v(8)  + qx(4)*v(11) +
     $       	    qx(5)*v(14) + qx(6)*v(17) + qx(7)*v(20) + 
     $      	    qx(8)*v(23)
            gamyz = qy(1)*v(3)  + qy(2)*v(6)  + qy(3)*v(9)  + 
     $       	    qy(4)*v(12) + qy(5)*v(15) + qy(6)*v(18) +
     $       	    qy(7)*v(21) + qy(8)*v(24) + qz(1)*v(2)  + 
     $      	    qz(2)*v(5)  + qz(3)*v(8)  + qz(4)*v(11) +
     $       	    qz(5)*v(14) + qz(6)*v(17) + qz(7)*v(20) + 
     $      	    qz(8)*v(23)
            gamxz = qz(1)*v(1)  + qz(2)*v(4)  + qz(3)*v(7)  +
     $       	    qz(4)*v(10) + qz(5)*v(13) + qz(6)*v(16) +
     $       	    qz(7)*v(19) + qz(8)*v(22) + qx(1)*v(3)  + 
     $      	    qx(2)*v(6)  + qx(3)*v(9)  + qx(4)*v(12) +
     $       	    qx(5)*v(15) + qx(6)*v(18) + qx(7)*v(21) + 
     $      	    qx(8)*v(24)
C
            depsxx =dqx(1)*v(1)  + dqx(2)*v(4)  + dqx(3)*v(7)  +
     $       	    dqx(4)*v(10) + dqx(5)*v(13) + dqx(6)*v(16) + 
     $      	    dqx(7)*v(19) + dqx(8)*v(22)
            depsyy =dqy(1)*v(2)  + dqy(2)*v(5)  + dqy(3)*v(8)  +
     $       	    dqy(4)*v(11) + dqy(5)*v(14) + dqy(6)*v(17) + 
     $      	    dqy(7)*v(20) + dqy(8)*v(23)
            depszz =dqz(1)*v(3)  + dqz(2)*v(6)  + dqz(3)*v(9)  + 
     $       	    dqz(4)*v(12) + dqz(5)*v(15) + dqz(6)*v(18) + 
     $      	    dqz(7)*v(21) + dqz(8)*v(24)
            dgamxy =dqy(1)*v(1)  + dqy(2)*v(4)  + dqy(3)*v(7)  + 
     $       	    dqy(4)*v(10) + dqy(5)*v(13) + dqy(6)*v(16) + 
     $       	    dqy(7)*v(19) + dqy(8)*v(22) + dqx(1)*v(2)  + 
     $      	    dqx(2)*v(5)  + dqx(3)*v(8)  + dqx(4)*v(11) +
     $       	    dqx(5)*v(14) + dqx(6)*v(17) + dqx(7)*v(20) + 
     $      	    dqx(8)*v(23)
            dgamyz =dqy(1)*v(3)  + dqy(2)*v(6)  + dqy(3)*v(9)  + 
     $       	    dqy(4)*v(12) + dqy(5)*v(15) + dqy(6)*v(18) +
     $       	    dqy(7)*v(21) + dqy(8)*v(24) + dqz(1)*v(2)  + 
     $      	    dqz(2)*v(5)  + dqz(3)*v(8)  + dqz(4)*v(11) +
     $       	    dqz(5)*v(14) + dqz(6)*v(17) + dqz(7)*v(20) + 
     $      	    dqz(8)*v(23)
            dgamxz =dqz(1)*v(1)  + dqz(2)*v(4)  + dqz(3)*v(7)  +
     $       	    dqz(4)*v(10) + dqz(5)*v(13) + dqz(6)*v(16) +
     $       	    dqz(7)*v(19) + dqz(8)*v(22) + dqx(1)*v(3)  + 
     $      	    dqx(2)*v(6)  + dqx(3)*v(9)  + dqx(4)*v(12) +
     $       	    dqx(5)*v(15) + dqx(6)*v(18) + dqx(7)*v(21) + 
     $      	    dqx(8)*v(24)
C
C
            depsxx = depsxx+ qx(1)*dv(1) + qx(2)*dv(4)+ qx(3)*dv(7)  +
     $       	          qx(4)*dv(10) + qx(5)*dv(13) + qx(6)*dv(16) + 
     $      	          qx(7)*dv(19) + qx(8)*dv(22)
            depsyy = depsyy+ qy(1)*dv(2) + qy(2)*dv(5)+ qy(3)*dv(8)  +
     $       	          qy(4)*dv(11) + qy(5)*dv(14) + qy(6)*dv(17) + 
     $      	          qy(7)*dv(20) + qy(8)*dv(23)
            depszz = depszz+ qz(1)*dv(3) + qz(2)*dv(6)+ qz(3)*dv(9)  + 
     $       	          qz(4)*dv(12) + qz(5)*dv(15) + qz(6)*dv(18) + 
     $      	          qz(7)*dv(21) + qz(8)*dv(24)
            dgamxy = dgamxy+ qy(1)*dv(1) + qy(2)*dv(4)+ qy(3)*dv(7)  + 
     $       	          qy(4)*dv(10) + qy(5)*dv(13) + qy(6)*dv(16) + 
     $       	          qy(7)*dv(19) + qy(8)*dv(22) + qx(1)*dv(2)  + 
     $      	          qx(2)*dv(5)  + qx(3)*dv(8)  + qx(4)*dv(11) +
     $       	          qx(5)*dv(14) + qx(6)*dv(17) + qx(7)*dv(20) + 
     $      	          qx(8)*dv(23)
            dgamyz = dgamyz+ qy(1)*dv(3) + qy(2)*dv(6)+ qy(3)*dv(9)  + 
     $       	          qy(4)*dv(12) + qy(5)*dv(15) + qy(6)*dv(18) +
     $       	          qy(7)*dv(21) + qy(8)*dv(24) + qz(1)*dv(2)  + 
     $      	          qz(2)*dv(5)  + qz(3)*dv(8)  + qz(4)*dv(11) +
     $       	          qz(5)*dv(14) + qz(6)*dv(17) + qz(7)*dv(20) + 
     $      	          qz(8)*dv(23)
            dgamxz = dgamxz+ qz(1)*dv(1) + qz(2)*dv(4)+ qz(3)*dv(7)  +
     $       	          qz(4)*dv(10) + qz(5)*dv(13) + qz(6)*dv(16) +
     $       	          qz(7)*dv(19) + qz(8)*dv(22) + qx(1)*dv(3)  + 
     $      	          qx(2)*dv(6)  + qx(3)*dv(9)  + qx(4)*dv(12) +
     $       	          qx(5)*dv(15) + qx(6)*dv(18) + qx(7)*dv(21) + 
     $      	          qx(8)*dv(24)
C
C.... ENGINEERING STRESS COMPUTATION
C
            stress(1) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*epszz
     $       	      + c(1,4)*gamxy + c(1,5)*gamyz + c(1,6)*gamxz
            stress(2) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*epszz
     $       	      + c(2,4)*gamxy + c(2,5)*gamyz + c(2,6)*gamxz
            stress(3) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*epszz
     $       	      + c(3,4)*gamxy + c(3,5)*gamyz + c(3,6)*gamxz
            stress(4) = c(4,1)*epsxx + c(4,2)*epsyy + c(4,3)*epszz
     $       	      + c(4,4)*gamxy + c(4,5)*gamyz + c(4,6)*gamxz
            stress(5) = c(5,1)*epsxx + c(5,2)*epsyy + c(5,3)*epszz
     $       	      + c(5,4)*gamxy + c(5,5)*gamyz + c(5,6)*gamxz
            stress(6) = c(6,1)*epsxx + c(6,2)*epsyy + c(6,3)*epszz
     $       	      + c(6,4)*gamxy + c(6,5)*gamyz + c(6,6)*gamxz

            dstress(1) = dc(1,1)*epsxx + dc(1,2)*epsyy + dc(1,3)*epszz
     $       	       + dc(1,4)*gamxy + dc(1,5)*gamyz + dc(1,6)*gamxz
            dstress(2) = dc(2,1)*epsxx + dc(2,2)*epsyy + dc(2,3)*epszz
     $       	       + dc(2,4)*gamxy + dc(2,5)*gamyz + dc(2,6)*gamxz
            dstress(3) = dc(3,1)*epsxx + dc(3,2)*epsyy + dc(3,3)*epszz
     $       	       + dc(3,4)*gamxy + dc(3,5)*gamyz + dc(3,6)*gamxz
            dstress(4) = dc(4,1)*epsxx + dc(4,2)*epsyy + dc(4,3)*epszz
     $       	       + dc(4,4)*gamxy + dc(4,5)*gamyz + dc(4,6)*gamxz
            dstress(5) = dc(5,1)*epsxx + dc(5,2)*epsyy + dc(5,3)*epszz
     $       	       + dc(5,4)*gamxy + dc(5,5)*gamyz + dc(5,6)*gamxz
            dstress(6) = dc(6,1)*epsxx + dc(6,2)*epsyy + dc(6,3)*epszz
     $       	       + dc(6,4)*gamxy + dc(6,5)*gamyz + dc(6,6)*gamxz

            dstress(1) = dstress(1)   + c(1,1)*depsxx + c(1,2)*depsyy+
     $       	        c(1,3)*depszz + c(1,4)*dgamxy + c(1,5)*dgamyz+ 
     $                  c(1,6)*dgamxz
            dstress(2) = dstress(2)  + c(2,1)*depsxx + c(2,2)*depsyy+ 
     $       	       c(2,3)*depszz + c(2,4)*dgamxy + c(2,5)*dgamyz+ 
     $                 c(2,6)*dgamxz
            dstress(3) = dstress(3)  + c(3,1)*depsxx + c(3,2)*depsyy+ 
     $       	       c(3,3)*depszz + c(3,4)*dgamxy + c(3,5)*dgamyz+ 
     $                 c(3,6)*dgamxz
            dstress(4) = dstress(4)  + c(4,1)*depsxx + c(4,2)*depsyy+ 
     $       	       c(4,3)*depszz + c(4,4)*dgamxy + c(4,5)*dgamyz+ 
     $                 c(4,6)*dgamxz
            dstress(5) = dstress(5)  + c(5,1)*depsxx + c(5,2)*depsyy+ 
     $       	       c(5,3)*depszz + c(5,4)*dgamxy + c(5,5)*dgamyz+ 
     $                 c(5,6)*dgamxz
            dstress(6) = dstress(6)  + c(6,1)*depsxx + c(6,2)*depsyy+ 
     $       	       c(6,3)*depszz + c(6,4)*dgamxy + c(6,5)*dgamyz+ 
     $                 c(6,6)*dgamxz
C
C.... COMPUTE THE FIRST DEVEATORIC STRESSES
C
	    comp = (stress(1) + stress(2) + stress(3))/3.0d0
	    dsxx = stress(1) - comp
            dsyy = stress(2) - comp
            dszz = stress(3) - comp
            dsxy = stress(4)
            dsyz = stress(5)
            dsxz = stress(6)

	    dcomp = (dstress(1) + dstress(2) + dstress(3))/3.0d0
	    ddsxx = dstress(1) - dcomp
            ddsyy = dstress(2) - dcomp
            ddszz = dstress(3) - dcomp
            ddsxy = dstress(4)
            ddsyz = dstress(5)
            ddsxz = dstress(6)
C
C.... COMPUTE THE SECOND DEVEATORIC STRESS
C
	     j2    = ((dsxx*dsxx)+(dsyy*dsyy)+(dszz*dszz))/2.0d0 +
     &                (dsxy*dsxy)+(dsyz*dsyz)+(dsxz*dsxz)
             j2    = dsqrt(3.0d0*j2)

	     derj2 = (ddsxx*dsxx)+(ddsyy*dsyy)+(ddszz*dszz)+
     &              ((ddsxy*dsxy)+(ddsyz*dsyz)+(ddsxz*dsxz))*2.0d0

	     dj2   = 1.5d0*derj2/j2
C
C.... COMPUTE THE VON MISES STRESS AND VOLUME
C
             if (iarea.eq.0) then 
               w  = 1.0d0
               dw = 0.0d0
             endif

	     vmint   = vmint + w * pow((j2/sigbar),fac)

	     dfac    = fac - 1.0d0
             dervmint= dw * pow((j2/sigbar),fac) + 
     &                  w * fac * (dj2/sigbar) * pow((j2/sigbar),dfac) 
             dvmint  = dvmint + dervmint

	     vol     = vol  + w
	     dvol    = dvol + dw
C
 2400       continue
 2500     continue
 3000   continue
      
      endif

C
C.... COMPUTE MAXIMUM STRESS VALUE AT GAUSS POINT; RETURN MAXVAL^FAC
C
       
      if (iarea.eq.0) then 
        vmint  = vmint /vol
        dvmint = dvmint/vol
        vol    = 1.0d0
        dvol   = 0.0d0
      endif

      return
      end
