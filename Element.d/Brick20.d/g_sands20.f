C--------------------------------------------------------------------
C    Given the node displacement,computes the derivatives of
C    corner stresses on a 20-node hexahedron.
C---------------------------------------------------------------------
          subroutine gxsands20(elm,x,y,z,dx,dy,dz,c,dc,u,du,stress,
     &                       dstress,maxgus,maxstr,maxsze,vmflg,
     &                       strainFlg,strind)
C
C                A R G U M E N T S
C
         double precision  x(20), y(20), z(20), c(6,6), u(60)
         double precision  dx(20), dy(20), dz(20), dc(6,6), du(60)
         real*8  stress(maxsze,maxstr,maxgus) 
         real*8  dstress(maxsze,maxstr,maxgus)
         integer elm,maxgus,maxstr,maxsze, strind
         logical vmflg,strainFlg
C
C                L O C A L  V A R I A B L E S
C
         double precision  q(20), qx(20), qy(20), qz(20)
         double precision  dqx(20), dqy(20), dqz(20), ddet
         double precision  xicorn(8), etacorn(8), mucorn(8)
         double precision  xi, eta, mu, det, epsxx, epsyy
         double precision  epszz, gamxy, gamyz, gamxz
         double precision  depsxx, depsyy
         double precision  depszz, dgamxy, dgamyz, dgamxz, cvar
         double precision  half 
C          
         integer           ipolm(12),ipoll(12),ipolr(12)
         integer           i, n, nm, nl, nr
C
C                D A T A
C
         data    xicorn  /-1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0, -1.0/
         data    etacorn /-1.0, -1.0,  1.0,  1.0, -1.0, -1.0, 1.0,  1.0/
         data    mucorn  /-1.0, -1.0, -1.0, -1.0,  1.0,  1.0, 1.0,  1.0/
c
         data    ipolm   / 9,10,11,12,13,14,15,16,17,18,19,20/
         data    ipoll   / 1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 8/
         data    ipolr   / 2, 3, 4, 1, 5, 6, 7, 8, 6, 7, 8, 1/
c
         half = 0.5d0
C
C     check if any variation
      cvar = 0.0d0
      
      do 5000 i = 1, 8
        cvar = cvar + dx(i)**2 + dy(i)**2 +dz(i)**2 
 5000 continue

c------------------------------------------------ stress calculation
       if(strInd .le. 6) then
c                                                          cvar.eq.0       
       if (cvar.eq.0) then

         do  3000  n = 1,8
            xi  = xicorn(n)
            eta = etacorn(n)
            mu = mucorn(n)
            call h20shpe (xi,eta,mu,x,y,z,q,qx,qy,qz,det)
            if ( det .le. 0.0d0) then
              write(6,*) 'Negative  Jacobian determinant'
            if ( det .eq. 0.0d0) then
              write(6,*) 'Zero Jacobian determinant'
            endif
          return
         endif
	 
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
         do 2500 i = 1,20
            epsxx =   epsxx + qx(i)*u(3*i-2)
            epsyy =   epsyy + qy(i)*u(3*i-1) 
            epszz =   epszz + qz(i)*u(3*i  )  
            gamxy =   gamxy + qy(i)*u(3*i-2) + qx(i)*u(3*i-1)
            gamyz =   gamyz + qz(i)*u(3*i-1) + qy(i)*u(3*i  )
            gamxz =   gamxz + qx(i)*u(3*i  ) + qz(i)*u(3*i-2)
	    
            depsxx =   depsxx + qx(i)*du(3*i-2)
            depsyy =   depsyy + qy(i)*du(3*i-1) 
            depszz =   depszz + qz(i)*du(3*i  )  
            dgamxy =   dgamxy + qy(i)*du(3*i-2) + qx(i)*du(3*i-1)
            dgamyz =   dgamyz + qz(i)*du(3*i-1) + qy(i)*du(3*i  )
            dgamxz =   dgamxz + qx(i)*du(3*i  ) + qz(i)*du(3*i-2)
 2500       continue

         stress(elm,1,n) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*epszz
     $                   + c(1,4)*gamxy + c(1,5)*gamyz + c(1,6)*gamxz
         stress(elm,2,n) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*epszz
     $                   + c(2,4)*gamxy + c(2,5)*gamyz + c(2,6)*gamxz
         stress(elm,3,n) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*epszz
     $                   + c(3,4)*gamxy + c(3,5)*gamyz + c(3,6)*gamxz
         stress(elm,4,n) = c(4,1)*epsxx + c(4,2)*epsyy + c(4,3)*epszz
     $                   + c(4,4)*gamxy + c(4,5)*gamyz + c(4,6)*gamxz
         stress(elm,5,n) = c(5,1)*epsxx + c(5,2)*epsyy + c(5,3)*epszz
     $                   + c(5,4)*gamxy + c(5,5)*gamyz + c(5,6)*gamxz
         stress(elm,6,n) = c(6,1)*epsxx + c(6,2)*epsyy + c(6,3)*epszz
     $                   + c(6,4)*gamxy + c(6,5)*gamyz + c(6,6)*gamxz

         dstress(elm,1,n) =dc(1,1)*epsxx +dc(1,2)*epsyy +dc(1,3)*epszz
     $                    +dc(1,4)*gamxy +dc(1,5)*gamyz +dc(1,6)*gamxz
         dstress(elm,2,n) =dc(2,1)*epsxx +dc(2,2)*epsyy +dc(2,3)*epszz
     $                    +dc(2,4)*gamxy +dc(2,5)*gamyz +dc(2,6)*gamxz
         dstress(elm,3,n) =dc(3,1)*epsxx +dc(3,2)*epsyy +dc(3,3)*epszz
     $                    +dc(3,4)*gamxy +dc(3,5)*gamyz +dc(3,6)*gamxz
         dstress(elm,4,n) =dc(4,1)*epsxx +dc(4,2)*epsyy +dc(4,3)*epszz
     $                    +dc(4,4)*gamxy +dc(4,5)*gamyz +dc(4,6)*gamxz
         dstress(elm,5,n) =dc(5,1)*epsxx +dc(5,2)*epsyy +dc(5,3)*epszz
     $                    +dc(5,4)*gamxy +dc(5,5)*gamyz +dc(5,6)*gamxz
         dstress(elm,6,n) =dc(6,1)*epsxx +dc(6,2)*epsyy +dc(6,3)*epszz
     $                    +dc(6,4)*gamxy +dc(6,5)*gamyz +dc(6,6)*gamxz
     
         dstress(elm,1,n) =dstress(elm,1,n) + c(1,1)*depsxx
     $                    +c(1,2)*depsyy +c(1,3)*depszz
     $                    +c(1,4)*dgamxy +c(1,5)*dgamyz +c(1,6)*dgamxz
         dstress(elm,2,n) =dstress(elm,2,n) + c(2,1)*depsxx 
     $                    +c(2,2)*depsyy +c(2,3)*depszz
     $                    +c(2,4)*dgamxy +c(2,5)*dgamyz +c(2,6)*dgamxz
         dstress(elm,3,n) =dstress(elm,3,n) + c(3,1)*depsxx 
     $                    +c(3,2)*depsyy +c(3,3)*depszz
     $                    +c(3,4)*dgamxy +c(3,5)*dgamyz +c(3,6)*dgamxz
         dstress(elm,4,n) =dstress(elm,4,n) + c(4,1)*depsxx 
     $                    +c(4,2)*depsyy +c(4,3)*depszz
     $                    +c(4,4)*dgamxy +c(4,5)*dgamyz +c(4,6)*dgamxz
         dstress(elm,5,n) =dstress(elm,5,n) + c(5,1)*depsxx 
     $                    +c(5,2)*depsyy +c(5,3)*depszz
     $                    +c(5,4)*dgamxy +c(5,5)*dgamyz +c(5,6)*dgamxz
         dstress(elm,6,n) =dstress(elm,6,n) + c(6,1)*depsxx 
     $                    +c(6,2)*depsyy +c(6,3)*depszz
     $                    +c(6,4)*dgamxy +c(6,5)*dgamyz +c(6,6)*dgamxz

 3000    continue

C                                                          cvar.ne.0 
      else

         do  3010  n = 1,8
            xi  = xicorn(n)
            eta = etacorn(n)
            mu = mucorn(n)
            call gxh20shpe (xi,eta,mu,x,y,z,dx,dy,dz,q,qx,qy,qz,
     *                      dqx,dqy,dqz,det,ddet)
            if ( det .le. 0.0d0) then
              write(6,*) 'Negative  Jacobian determinant'
            if ( det .eq. 0.0d0) then
              write(6,*) 'Zero Jacobian determinant'
            endif
          return
         endif
	 
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
         do 2510 i = 1,20
            epsxx =   epsxx + qx(i)*u(3*i-2)
            epsyy =   epsyy + qy(i)*u(3*i-1) 
            epszz =   epszz + qz(i)*u(3*i  )  
            gamxy =   gamxy + qy(i)*u(3*i-2) + qx(i)*u(3*i-1)
            gamyz =   gamyz + qz(i)*u(3*i-1) + qy(i)*u(3*i  )
            gamxz =   gamxz + qx(i)*u(3*i  ) + qz(i)*u(3*i-2)
	    
            depsxx =   depsxx + dqx(i)*u(3*i-2)
            depsyy =   depsyy + dqy(i)*u(3*i-1) 
            depszz =   depszz + dqz(i)*u(3*i  )  
            dgamxy =   dgamxy + dqy(i)*u(3*i-2) + dqx(i)*u(3*i-1)
            dgamyz =   dgamyz + dqz(i)*u(3*i-1) + dqy(i)*u(3*i  )
            dgamxz =   dgamxz + dqx(i)*u(3*i  ) + dqz(i)*u(3*i-2)
	    
            depsxx =   depsxx + qx(i)*du(3*i-2)
            depsyy =   depsyy + qy(i)*du(3*i-1) 
            depszz =   depszz + qz(i)*du(3*i  )  
            dgamxy =   dgamxy + qy(i)*du(3*i-2) + qx(i)*du(3*i-1)
            dgamyz =   dgamyz + qz(i)*du(3*i-1) + qy(i)*du(3*i  )
            dgamxz =   dgamxz + qx(i)*du(3*i  ) + qz(i)*du(3*i-2)
 2510       continue

         stress(elm,1,n) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*epszz
     $                   + c(1,4)*gamxy + c(1,5)*gamyz + c(1,6)*gamxz
         stress(elm,2,n) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*epszz
     $                   + c(2,4)*gamxy + c(2,5)*gamyz + c(2,6)*gamxz
         stress(elm,3,n) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*epszz
     $                   + c(3,4)*gamxy + c(3,5)*gamyz + c(3,6)*gamxz
         stress(elm,4,n) = c(4,1)*epsxx + c(4,2)*epsyy + c(4,3)*epszz
     $                   + c(4,4)*gamxy + c(4,5)*gamyz + c(4,6)*gamxz
         stress(elm,5,n) = c(5,1)*epsxx + c(5,2)*epsyy + c(5,3)*epszz
     $                   + c(5,4)*gamxy + c(5,5)*gamyz + c(5,6)*gamxz
         stress(elm,6,n) = c(6,1)*epsxx + c(6,2)*epsyy + c(6,3)*epszz
     $                   + c(6,4)*gamxy + c(6,5)*gamyz + c(6,6)*gamxz

         dstress(elm,1,n) =dc(1,1)*epsxx +dc(1,2)*epsyy +dc(1,3)*epszz
     $                    +dc(1,4)*gamxy +dc(1,5)*gamyz +dc(1,6)*gamxz
         dstress(elm,2,n) =dc(2,1)*epsxx +dc(2,2)*epsyy +dc(2,3)*epszz
     $                    +dc(2,4)*gamxy +dc(2,5)*gamyz +dc(2,6)*gamxz
         dstress(elm,3,n) =dc(3,1)*epsxx +dc(3,2)*epsyy +dc(3,3)*epszz
     $                    +dc(3,4)*gamxy +dc(3,5)*gamyz +dc(3,6)*gamxz
         dstress(elm,4,n) =dc(4,1)*epsxx +dc(4,2)*epsyy +dc(4,3)*epszz
     $                    +dc(4,4)*gamxy +dc(4,5)*gamyz +dc(4,6)*gamxz
         dstress(elm,5,n) =dc(5,1)*epsxx +dc(5,2)*epsyy +dc(5,3)*epszz
     $                    +dc(5,4)*gamxy +dc(5,5)*gamyz +dc(5,6)*gamxz
         dstress(elm,6,n) =dc(6,1)*epsxx +dc(6,2)*epsyy +dc(6,3)*epszz
     $                    +dc(6,4)*gamxy +dc(6,5)*gamyz +dc(6,6)*gamxz
     
         dstress(elm,1,n) =dstress(elm,1,n) + c(1,1)*depsxx
     $                    +c(1,2)*depsyy +c(1,3)*depszz
     $                    +c(1,4)*dgamxy +c(1,5)*dgamyz +c(1,6)*dgamxz
         dstress(elm,2,n) =dstress(elm,2,n) + c(2,1)*depsxx 
     $                    +c(2,2)*depsyy +c(2,3)*depszz
     $                    +c(2,4)*dgamxy +c(2,5)*dgamyz +c(2,6)*dgamxz
         dstress(elm,3,n) =dstress(elm,3,n) + c(3,1)*depsxx 
     $                    +c(3,2)*depsyy +c(3,3)*depszz
     $                    +c(3,4)*dgamxy +c(3,5)*dgamyz +c(3,6)*dgamxz
         dstress(elm,4,n) =dstress(elm,4,n) + c(4,1)*depsxx 
     $                    +c(4,2)*depsyy +c(4,3)*depszz
     $                    +c(4,4)*dgamxy +c(4,5)*dgamyz +c(4,6)*dgamxz
         dstress(elm,5,n) =dstress(elm,5,n) + c(5,1)*depsxx 
     $                    +c(5,2)*depsyy +c(5,3)*depszz
     $                    +c(5,4)*dgamxy +c(5,5)*dgamyz +c(5,6)*dgamxz
         dstress(elm,6,n) =dstress(elm,6,n) + c(6,1)*depsxx 
     $                    +c(6,2)*depsyy +c(6,3)*depszz
     $                    +c(6,4)*dgamxy +c(6,5)*dgamyz +c(6,6)*dgamxz

 3010    continue

       endif	
C
C.... COMPUTE THE VON MISES STRESS
C
       if (vmflg) then
          call gxvmelmv(stress,dstress,maxgus,maxstr,maxsze,elm,8)
       endif
C------------------------------------------------- strain calculation      
      else
C                                                          cvar.eq.0        
       if (cvar.eq.0) then

         do  2000  n = 1,8
            xi  = xicorn(n)
            eta = etacorn(n)
            mu = mucorn(n)
            call h20shpe (xi,eta,mu,x,y,z,q,qx,qy,qz,det)
            if ( det .le. 0.0d0) then
              write(6,*) 'Negative  Jacobian determinant'
            if ( det .eq. 0.0d0) then
              write(6,*) 'Zero Jacobian determinant'
            endif
          return
         endif
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
            epsxx =   epsxx + qx(i)*u(3*i-2)
            epsyy =   epsyy + qy(i)*u(3*i-1) 
            epszz =   epszz + qz(i)*u(3*i  )  
            gamxy =   gamxy + qy(i)*u(3*i-2) + qx(i)*u(3*i-1)
            gamyz =   gamyz + qz(i)*u(3*i-1) + qy(i)*u(3*i  )
            gamxz =   gamxz + qx(i)*u(3*i  ) + qz(i)*u(3*i-2)
	    
            depsxx =   depsxx + qx(i)*du(3*i-2)
            depsyy =   depsyy + qy(i)*du(3*i-1) 
            depszz =   depszz + qz(i)*du(3*i  )  
            dgamxy =   dgamxy + qy(i)*du(3*i-2) + qx(i)*du(3*i-1)
            dgamyz =   dgamyz + qz(i)*du(3*i-1) + qy(i)*du(3*i  )
            dgamxz =   dgamxz + qx(i)*du(3*i  ) + qz(i)*du(3*i-2)
 1500       continue
 
         stress(elm,1,n) = epsxx
         stress(elm,2,n) = epsyy 
         stress(elm,3,n) = epszz 
         stress(elm,4,n) = gamxy 
         stress(elm,5,n) = gamyz 
         stress(elm,6,n) = gamxz 
C
         dstress(elm,1,n) = depsxx
         dstress(elm,2,n) = depsyy 
         dstress(elm,3,n) = depszz 
         dstress(elm,4,n) = dgamxy 
         dstress(elm,5,n) = dgamyz 
         dstress(elm,6,n) = dgamxz 
 2000    continue

C                                                          cvar.ne.0       
      else

         do  2010  n = 1,8
            xi  = xicorn(n)
            eta = etacorn(n)
            mu = mucorn(n)
            call gxh20shpe (xi,eta,mu,x,y,z,dx,dy,dz,q,qx,qy,qz,
     *                      dqx,dqy,dqz,det,ddet)
            if ( det .le. 0.0d0) then
              write(6,*) 'Negative  Jacobian determinant'
            if ( det .eq. 0.0d0) then
              write(6,*) 'Zero Jacobian determinant'
            endif
          return
         endif
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
         do 1510 i = 1,20
            epsxx =   epsxx + qx(i)*u(3*i-2)
            epsyy =   epsyy + qy(i)*u(3*i-1) 
            epszz =   epszz + qz(i)*u(3*i  )  
            gamxy =   gamxy + qy(i)*u(3*i-2) + qx(i)*u(3*i-1)
            gamyz =   gamyz + qz(i)*u(3*i-1) + qy(i)*u(3*i  )
            gamxz =   gamxz + qx(i)*u(3*i  ) + qz(i)*u(3*i-2)
	    
            depsxx =   depsxx + dqx(i)*u(3*i-2)
            depsyy =   depsyy + dqy(i)*u(3*i-1) 
            depszz =   depszz + dqz(i)*u(3*i  )  
            dgamxy =   dgamxy + dqy(i)*u(3*i-2) + dqx(i)*u(3*i-1)
            dgamyz =   dgamyz + dqz(i)*u(3*i-1) + dqy(i)*u(3*i  )
            dgamxz =   dgamxz + dqx(i)*u(3*i  ) + dqz(i)*u(3*i-2)
	    
            depsxx =   depsxx + qx(i)*du(3*i-2)
            depsyy =   depsyy + qy(i)*du(3*i-1) 
            depszz =   depszz + qz(i)*du(3*i  )  
            dgamxy =   dgamxy + qy(i)*du(3*i-2) + qx(i)*du(3*i-1)
            dgamyz =   dgamyz + qz(i)*du(3*i-1) + qy(i)*du(3*i  )
            dgamxz =   dgamxz + qx(i)*du(3*i  ) + qz(i)*du(3*i-2)
 1510       continue
 
         stress(elm,1,n) = epsxx
         stress(elm,2,n) = epsyy 
         stress(elm,3,n) = epszz 
         stress(elm,4,n) = gamxy 
         stress(elm,5,n) = gamyz 
         stress(elm,6,n) = gamxz 
C
         dstress(elm,1,n) = depsxx
         dstress(elm,2,n) = depsyy 
         dstress(elm,3,n) = depszz 
         dstress(elm,4,n) = dgamxy 
         dstress(elm,5,n) = dgamyz 
         dstress(elm,6,n) = dgamxz 
 2010    continue

       endif

C.... COMPUTE EQUIVALENT STRAIN
C
	if(strainFlg) then
	  call gxstrainvm(stress,dstress,maxgus,maxstr,maxsze,8)
	endif
      	
      endif
C
C.... Interpolate stresses
C
      do i=1,12
c
        nm=ipolm(i)
        nl=ipoll(i)
        nr=ipolr(i)
c
        dstress(elm,1,nm) = half*(dstress(elm,1,nl)+dstress(elm,1,nr))
        dstress(elm,2,nm) = half*(dstress(elm,2,nl)+dstress(elm,2,nr))
        dstress(elm,3,nm) = half*(dstress(elm,3,nl)+dstress(elm,3,nr))
        dstress(elm,4,nm) = half*(dstress(elm,4,nl)+dstress(elm,4,nr))
        dstress(elm,5,nm) = half*(dstress(elm,5,nl)+dstress(elm,5,nr))      
        dstress(elm,6,nm) = half*(dstress(elm,6,nl)+dstress(elm,6,nr))
        dstress(elm,7,nm) = half*(dstress(elm,7,nl)+dstress(elm,7,nr))

      enddo
      
      return
      end
C
