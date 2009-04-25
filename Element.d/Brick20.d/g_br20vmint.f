C=PURPOSE Compute integral of (Sigma_vm / Sigma_bar) of 8-node hexa
C
C     input arguments are
C
C       X         (20 x 1) array of x coordinates of hexahedron nodes
C       Y         (20 x 1) array of y coordinates of hexahedron nodes
C       Z         (20 x 1) array of z coordinates of hexahedron nodes
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
      subroutine  gxbr20vmint(x,y,z,dx,dy,dz,v,dv,c,dc,p,vmint,dvmint,
     *                        vol,dvol,sigbar,fac,iarea,status)
C
C                   A R G U M E N T S
C
C      implicit none
      integer p,status
      real*8  x(20), y(20), z(20), c(6,6), v(*)
      real*8  dx(20), dy(20), dz(20), dc(6,6), dv(*)
      real*8  vmint, dvmint, sigbar, fac, vol, dvol
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(20), qx(20), qy(20),qz(20)
      real*8  dqx(20), dqy(20), dqz(20)
      real*8  xi, eta, emu, det, w, weight, ddet, dw, dfac, dervmint
      real*8  epsxx, epsyy,epszz, gamxy, gamyz, gamxz
      real*8  depsxx, depsyy,depszz, dgamxy, dgamyz, dgamxz
      real*8  stress(6),dsxx, dsyy, dszz, dsxy, dsxz, dsyz, j2, derj2
      real*8  dstress(6),ddsxx, ddsyy, ddszz, ddsxy, ddsxz, ddsyz, dj2
      real*8  a,b,cvar, dcomp, pow, comp
      integer k, l, jj, i
C
C
C                   L O G I C
C
C                             inline function
C
      pow(a,b) = exp(b*log(a))

C     Initialize 
      status = 0
      vmint  = 0.0d0
      dvmint = 0.0d0
      vol    = 0.0d0
      dvol   = 0.0d0

C     check if any variation
C       
C       do 5000 i = 1, 20
C         cvar = cvar + dx(i)**2 + dy(i)**2 +dz(i)**2 
C  5000 continue

      cvar = 1.0d0     
c                                                          cvar.eq.0       
      if (cvar.eq.0) then
C
      do 3000  k = 1,p
        do 2500  l = 1,p
          do 2400  jj = 1,p
            call hxgaus20(p,k,p,l,p,jj,xi,eta,emu,weight)
            call h20shpe (xi,eta,emu,x,y,z,q,qx,qy,qz,det)
C
            if (det .le. 0.0d0) then
              write(6,*) 'Negative Jacobian determinant in brintvm.f'
              status = -1
              return
            end if
C
            w  = weight * det
C
            epsxx = 0.0d0
            epsyy = 0.0d0
            epszz = 0.0d0
            gamxy = 0.0d0
            gamyz = 0.0d0
            gamxz = 0.0d0
	
            depsxx = 0.0d0
            depsyy = 0.0d0
            depszz = 0.0d0
            dgamxy = 0.0d0
            dgamyz = 0.0d0
            dgamxz = 0.0d0
	
            do 1500 i = 1,20
	 
              epsxx =epsxx +qx(i)*v(3*i-2)
              epsyy =epsyy +qy(i)*v(3*i-1) 
              epszz =epszz +qz(i)*v(3*i  )  
              gamxy =gamxy +qy(i)*v(3*i-2) +qx(i)*v(3*i-1)
              gamyz =gamyz +qz(i)*v(3*i-1) +qy(i)*v(3*i  )
              gamxz =gamxz +qx(i)*v(3*i  ) +qz(i)*v(3*i-2)

              depsxx =depsxx +qx(i)*dv(3*i-2)
              depsyy =depsyy +qy(i)*dv(3*i-1) 
              depszz =depszz +qz(i)*dv(3*i  )  
              dgamxy =dgamxy +qy(i)*dv(3*i-2) +qx(i)*dv(3*i-1)
              dgamyz =dgamyz +qz(i)*dv(3*i-1) +qy(i)*dv(3*i  )
              dgamxz =dgamxz +qx(i)*dv(3*i  ) +qz(i)*dv(3*i-2)
 1500       continue

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

            dstress(1) = dstress(1) 
     $                 + c(1,1)*depsxx + c(1,2)*depsyy + c(1,3)*depszz
     $       	       + c(1,4)*dgamxy + c(1,5)*dgamyz + c(1,6)*dgamxz
            dstress(2) = dstress(2)
     $                 + c(2,1)*depsxx + c(2,2)*depsyy + c(2,3)*depszz
     $       	       + c(2,4)*dgamxy + c(2,5)*dgamyz + c(2,6)*dgamxz
            dstress(3) = dstress(3)
     $                 + c(3,1)*depsxx + c(3,2)*depsyy + c(3,3)*depszz
     $       	       + c(3,4)*dgamxy + c(3,5)*dgamyz + c(3,6)*dgamxz
            dstress(4) = dstress(4)
     $                 + c(4,1)*depsxx + c(4,2)*depsyy + c(4,3)*depszz
     $       	       + c(4,4)*dgamxy + c(4,5)*dgamyz + c(4,6)*dgamxz
            dstress(5) = dstress(5)
     $                 + c(5,1)*depsxx + c(5,2)*depsyy + c(5,3)*depszz
     $       	       + c(5,4)*dgamxy + c(5,5)*dgamyz + c(5,6)*dgamxz
            dstress(6) = dstress(6)
     $                 + c(6,1)*depsxx + c(6,2)*depsyy + c(6,3)*depszz
     $       	       + c(6,4)*dgamxy + c(6,5)*dgamyz + c(6,6)*dgamxz
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
	       dj2   = 3.0d0*derj2/(2.0d0*j2)
             endif
C
C.... COMPUTE THE VON MISES STRESS AND VOLUME
C
             if (iarea.eq.0) w  = 1.0d0

	     vmint   = vmint + w*pow((j2/sigbar),fac)

	     dfac    = fac - 1.0d0
             dervmint= w * fac * (dj2/sigbar) * pow((j2/sigbar),dfac) 
             dvmint  = dvmint + dervmint

	     vol     = vol  + w
C
 2400       continue
 2500     continue
 3000   continue
C                                                          cvar.ne.0 
      else
      
      do 6000  k = 1,p
        do 5500  l = 1,p
          do 5400  jj = 1,p
            call hxgaus20(p,k,p,l,p,jj,xi,eta,emu,weight)
            call gxh20shpe (xi,eta,emu,x,y,z,dx,dy,dz,q,qx,qy,qz,dqx,
     *                      dqy,dqz,det,ddet)
C
            if (det .le. 0.0d0) then
              write(6,*) 'Negative Jacobian determinant in brintvm.f'
              status = -1
              return
            end if
C
            w  = weight * det
            dw = weight * ddet
C
            epsxx = 0.0d0
            epsyy = 0.0d0
            epszz = 0.0d0
            gamxy = 0.0d0
            gamyz = 0.0d0
            gamxz = 0.0d0
	
            depsxx = 0.0d0
            depsyy = 0.0d0
            depszz = 0.0d0
            dgamxy = 0.0d0
            dgamyz = 0.0d0
            dgamxz = 0.0d0
	
            do 9500 i = 1,20
	 
              epsxx =epsxx +qx(i)*v(3*i-2)
              epsyy =epsyy +qy(i)*v(3*i-1) 
              epszz =epszz +qz(i)*v(3*i  )  
              gamxy =gamxy +qy(i)*v(3*i-2) +qx(i)*v(3*i-1)
              gamyz =gamyz +qz(i)*v(3*i-1) +qy(i)*v(3*i  )
              gamxz =gamxz +qx(i)*v(3*i  ) +qz(i)*v(3*i-2)

              depsxx =depsxx +dqx(i)*v(3*i-2)
              depsyy =depsyy +dqy(i)*v(3*i-1) 
              depszz =depszz +dqz(i)*v(3*i  )  
              dgamxy =dgamxy +dqy(i)*v(3*i-2) +dqx(i)*v(3*i-1)
              dgamyz =dgamyz +dqz(i)*v(3*i-1) +dqy(i)*v(3*i  )
              dgamxz =dgamxz +dqx(i)*v(3*i  ) +dqz(i)*v(3*i-2)

              depsxx =depsxx +qx(i)*dv(3*i-2)
              depsyy =depsyy +qy(i)*dv(3*i-1) 
              depszz =depszz +qz(i)*dv(3*i  )  
              dgamxy =dgamxy +qy(i)*dv(3*i-2) +qx(i)*dv(3*i-1)
              dgamyz =dgamyz +qz(i)*dv(3*i-1) +qy(i)*dv(3*i  )
              dgamxz =dgamxz +qx(i)*dv(3*i  ) +qz(i)*dv(3*i-2)
 9500       continue

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

            dstress(1) = dstress(1) 
     $                 + c(1,1)*depsxx + c(1,2)*depsyy + c(1,3)*depszz
     $       	       + c(1,4)*dgamxy + c(1,5)*dgamyz + c(1,6)*dgamxz
            dstress(2) = dstress(2)
     $                 + c(2,1)*depsxx + c(2,2)*depsyy + c(2,3)*depszz
     $       	       + c(2,4)*dgamxy + c(2,5)*dgamyz + c(2,6)*dgamxz
            dstress(3) = dstress(3)
     $                 + c(3,1)*depsxx + c(3,2)*depsyy + c(3,3)*depszz
     $       	       + c(3,4)*dgamxy + c(3,5)*dgamyz + c(3,6)*dgamxz
            dstress(4) = dstress(4)
     $                 + c(4,1)*depsxx + c(4,2)*depsyy + c(4,3)*depszz
     $       	       + c(4,4)*dgamxy + c(4,5)*dgamyz + c(4,6)*dgamxz
            dstress(5) = dstress(5)
     $                 + c(5,1)*depsxx + c(5,2)*depsyy + c(5,3)*depszz
     $       	       + c(5,4)*dgamxy + c(5,5)*dgamyz + c(5,6)*dgamxz
            dstress(6) = dstress(6)
     $                 + c(6,1)*depsxx + c(6,2)*depsyy + c(6,3)*depszz
     $       	       + c(6,4)*dgamxy + c(6,5)*dgamyz + c(6,6)*dgamxz
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

	     dj2   = 3.0d0*derj2/(2.0d0*j2)
C
C.... COMPUTE THE VON MISES STRESS AND VOLUME
C
             if (iarea.eq.0) then 
               w  = 1.0d0
               dw = 0.0d0
             endif

	     vmint   = vmint + w*pow((j2/sigbar),fac)

	     dfac    = fac - 1.0d0
             dervmint= dw * pow((j2/sigbar),fac) + 
     *                  w * fac * (dj2/sigbar) * pow((j2/sigbar),dfac) 
             dvmint  = dvmint + dervmint

	     vol     = vol  + w
	     dvol    = dvol + dw

C
 5400       continue
 5500     continue
 6000   continue

      endif  

C
C.... COMPUTE MAXIMUM STRESS VALUE AT GAUSS POINT; RETURN MAXVAL^FAC
C
       
      if (iarea.eq.0) then 
        vmint  = vmint/vol
        dvmint = dvmint/vol
        vol    = 1.0d0
        dvol   = 0.0d0
      endif

      return
      end
