         subroutine gxsands17(elm,x,y,z,dx,dy,dz,c,dc,v,dv,
     &                        stress,dstress,maxgus,maxstr,maxsze,
     &                        vmflg,strainFlg,strind)
C
C                A R G U M E N T S
C
c         implicit none
	 integer maxsze,maxstr,maxgus,elm,strind
         real*8  x(8), y(8), z(8), c(6,6), v(24)
         real*8  dx(8), dy(8), dz(8), dc(6,6), dv(24)
         real*8  dstress(maxsze,maxstr,maxgus)
         real*8  stress(maxsze,maxstr,maxgus)
         logical vmflg,strainFlg
C
C            L O C A L  V A R I A B L E S
C
         real*8  q(8), qx(8), qy(8), qz(8)
         real*8  dqx(8), dqy(8), dqz(8), ddet
         real*8  xinod(8), etanod(8), emunod(8)
         real*8  xi, eta, emu, det, epsxx, epsyy
         real*8  epszz, gamxy, gamyz, gamxz
         real*8  depsxx, depsyy
         real*8  depszz, dgamxy, dgamyz, dgamxz, cvar
c
         integer n, i
C
C               D   A   T   A
C
         data    xinod/-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0/
         data    etanod/-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0/
         data    emunod/-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0/
C
C               L   O   G   I    C
C
C     check if any variation
      cvar = 0.0
      
      do 5000 i = 1, 8
        cvar = cvar + dx(i)**2 + dy(i)**2 +dz(i)**2 
 5000 continue

c------------------------------------------------ stress calculation
       if(strInd .le. 6) then
c                                                          cvar.eq.0       
       if (cvar.eq.0) then
         do 1000  n = 1, 8
           xi  = xinod(n)
           eta = etanod(n)
           emu = emunod(n)
           call h8shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
C
           if (det .le. 0.0) then
             write(6,*)'Negative Jacobian determinant in g_sands17.f'
             if (det .eq. 0.0) then
               write(6,*)'Zero  jacobian determinant'
             endif
             stop  
           endif
C
         epsxx = qx(1)*v(1)  + qx(2)*v(4)  + qx(3)*v(7)  + qx(4)*v(10)
     $         + qx(5)*v(13) + qx(6)*v(16) + qx(7)*v(19) + qx(8)*v(22)
         epsyy = qy(1)*v(2)  + qy(2)*v(5)  + qy(3)*v(8)  + qy(4)*v(11)
     $         + qy(5)*v(14) + qy(6)*v(17) + qy(7)*v(20) + qy(8)*v(23)
         epszz = qz(1)*v(3)  + qz(2)*v(6)  + qz(3)*v(9)  + qz(4)*v(12)
     $         + qz(5)*v(15) + qz(6)*v(18) + qz(7)*v(21) + qz(8)*v(24)
         gamxy = qy(1)*v(1)  + qy(2)*v(4)  + qy(3)*v(7)  + qy(4)*v(10)
     $         + qy(5)*v(13) + qy(6)*v(16) + qy(7)*v(19) + qy(8)*v(22)
     $         + qx(1)*v(2)  + qx(2)*v(5)  + qx(3)*v(8)  + qx(4)*v(11)
     $         + qx(5)*v(14) + qx(6)*v(17) + qx(7)*v(20) + qx(8)*v(23)
         gamyz = qy(1)*v(3)  + qy(2)*v(6)  + qy(3)*v(9)  + qy(4)*v(12)
     $         + qy(5)*v(15) + qy(6)*v(18) + qy(7)*v(21) + qy(8)*v(24)
     $         + qz(1)*v(2)  + qz(2)*v(5)  + qz(3)*v(8)  + qz(4)*v(11)
     $         + qz(5)*v(14) + qz(6)*v(17) + qz(7)*v(20) + qz(8)*v(23)
         gamxz = qz(1)*v(1)  + qz(2)*v(4)  + qz(3)*v(7)  + qz(4)*v(10)
     $         + qz(5)*v(13) + qz(6)*v(16) + qz(7)*v(19) + qz(8)*v(22)
     $         + qx(1)*v(3)  + qx(2)*v(6)  + qx(3)*v(9)  + qx(4)*v(12)
     $         + qx(5)*v(15) + qx(6)*v(18) + qx(7)*v(21) + qx(8)*v(24)
     
         depsxx = 
     $          qx(1)*dv(1)  +qx(2)*dv(4)  +qx(3)*dv(7)  +qx(4)*dv(10)
     $         +qx(5)*dv(13) +qx(6)*dv(16) +qx(7)*dv(19) +qx(8)*dv(22)
         depsyy = 
     $          qy(1)*dv(2)  +qy(2)*dv(5)  +qy(3)*dv(8)  +qy(4)*dv(11)
     $         +qy(5)*dv(14) +qy(6)*dv(17) +qy(7)*dv(20) +qy(8)*dv(23)
         depszz = 
     $          qz(1)*dv(3)  +qz(2)*dv(6)  +qz(3)*dv(9)  +qz(4)*dv(12)
     $         +qz(5)*dv(15) +qz(6)*dv(18) +qz(7)*dv(21) +qz(8)*dv(24)
         dgamxy = 
     $          qy(1)*dv(1)  +qy(2)*dv(4)  +qy(3)*dv(7)  +qy(4)*dv(10)
     $         +qy(5)*dv(13) +qy(6)*dv(16) +qy(7)*dv(19) +qy(8)*dv(22)
     $         +qx(1)*dv(2)  +qx(2)*dv(5)  +qx(3)*dv(8)  +qx(4)*dv(11)
     $         +qx(5)*dv(14) +qx(6)*dv(17) +qx(7)*dv(20) +qx(8)*dv(23)
         dgamyz = 
     $          qy(1)*dv(3)  +qy(2)*dv(6)  +qy(3)*dv(9)  +qy(4)*dv(12)
     $         +qy(5)*dv(15) +qy(6)*dv(18) +qy(7)*dv(21) +qy(8)*dv(24)
     $         +qz(1)*dv(2)  +qz(2)*dv(5)  +qz(3)*dv(8)  +qz(4)*dv(11)
     $         +qz(5)*dv(14) +qz(6)*dv(17) +qz(7)*dv(20) +qz(8)*dv(23)
         dgamxz =
     $          qz(1)*dv(1)  +qz(2)*dv(4)  +qz(3)*dv(7)  +qz(4)*dv(10)
     $         +qz(5)*dv(13) +qz(6)*dv(16) +qz(7)*dv(19) +qz(8)*dv(22)
     $         +qx(1)*dv(3)  +qx(2)*dv(6)  +qx(3)*dv(9)  +qx(4)*dv(12)
     $         +qx(5)*dv(15) +qx(6)*dv(18) +qx(7)*dv(21) +qx(8)*dv(24)
     
C
C.... ENGINEERING STRESS COMPUTATION
C
         stress(elm,1,n) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*epszz
     $             + c(1,4)*gamxy + c(1,5)*gamyz + c(1,6)*gamxz
         stress(elm,2,n) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*epszz
     $             + c(2,4)*gamxy + c(2,5)*gamyz + c(2,6)*gamxz
         stress(elm,3,n) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*epszz
     $             + c(3,4)*gamxy + c(3,5)*gamyz + c(3,6)*gamxz
         stress(elm,4,n) = c(4,1)*epsxx + c(4,2)*epsyy + c(4,3)*epszz
     $             + c(4,4)*gamxy + c(4,5)*gamyz + c(4,6)*gamxz
         stress(elm,5,n) = c(5,1)*epsxx + c(5,2)*epsyy + c(5,3)*epszz
     $             + c(5,4)*gamxy + c(5,5)*gamyz + c(5,6)*gamxz
         stress(elm,6,n) = c(6,1)*epsxx + c(6,2)*epsyy + c(6,3)*epszz
     $             + c(6,4)*gamxy + c(6,5)*gamyz + c(6,6)*gamxz
C
         dstress(elm,1,n) = dc(1,1)*epsxx 
     $    	   + dc(1,2)*epsyy + dc(1,3)*epszz
     $             + dc(1,4)*gamxy + dc(1,5)*gamyz + dc(1,6)*gamxz
         dstress(elm,2,n) =   dc(2,1)*epsxx 
     $    	   + dc(2,2)*epsyy + dc(2,3)*epszz
     $             + dc(2,4)*gamxy + dc(2,5)*gamyz + dc(2,6)*gamxz
         dstress(elm,3,n) =  dc(3,1)*epsxx 
     $    	   + dc(3,2)*epsyy + dc(3,3)*epszz
     $             + dc(3,4)*gamxy + dc(3,5)*gamyz + dc(3,6)*gamxz
         dstress(elm,4,n) =  dc(4,1)*epsxx 
     $    	   + dc(4,2)*epsyy + dc(4,3)*epszz
     $             + dc(4,4)*gamxy + dc(4,5)*gamyz + dc(4,6)*gamxz
         dstress(elm,5,n) =   dc(5,1)*epsxx 
     $    	   + dc(5,2)*epsyy + dc(5,3)*epszz
     $             + dc(5,4)*gamxy + dc(5,5)*gamyz + dc(5,6)*gamxz
         dstress(elm,6,n) =   dc(6,1)*epsxx
     $    	   + dc(6,2)*epsyy + dc(6,3)*epszz
     $             + dc(6,4)*gamxy + dc(6,5)*gamyz + dc(6,6)*gamxz
C
         dstress(elm,1,n) =  dstress(elm,1,n) + c(1,1)*depsxx 
     $    	   + c(1,2)*depsyy + c(1,3)*depszz
     $             + c(1,4)*dgamxy + c(1,5)*dgamyz + c(1,6)*dgamxz
         dstress(elm,2,n) = dstress(elm,2,n) + c(2,1)*depsxx 
     $    	   + c(2,2)*depsyy + c(2,3)*depszz
     $             + c(2,4)*dgamxy + c(2,5)*dgamyz + c(2,6)*dgamxz
         dstress(elm,3,n) =  dstress(elm,3,n) + c(3,1)*depsxx 
     $    	   + c(3,2)*depsyy + c(3,3)*depszz
     $             + c(3,4)*dgamxy + c(3,5)*dgamyz + c(3,6)*dgamxz
         dstress(elm,4,n) =  dstress(elm,4,n) + c(4,1)*depsxx 
     $    	   + c(4,2)*depsyy + c(4,3)*depszz
     $             + c(4,4)*dgamxy + c(4,5)*dgamyz + c(4,6)*dgamxz
         dstress(elm,5,n) =  dstress(elm,5,n) + c(5,1)*depsxx 
     $    	   + c(5,2)*depsyy + c(5,3)*depszz
     $             + c(5,4)*dgamxy + c(5,5)*dgamyz + c(5,6)*dgamxz
         dstress(elm,6,n) = dstress(elm,6,n) +  c(6,1)*depsxx 
     $    	   + c(6,2)*depsyy + c(6,3)*depszz
     $             + c(6,4)*dgamxy + c(6,5)*dgamyz + c(6,6)*dgamxz

C
 1000    continue
C                                                          cvar.ne.0 
      else
      
         do 2000  n = 1, 8
           xi  = xinod(n)
           eta = etanod(n)
           emu = emunod(n)
           call gxh8shpe(xi,eta,emu,x,y,z,dx,dy,dz,q,qx,qy,qz,
     *                    dqx,dqy,dqz,det,ddet)
C
           if (det .le. 0.0) then
             write(6,*)'Negative Jacobian determinant in g_sands17.f'
             if (det .eq. 0.0) then
               write(6,*)'Zero  jacobian determinant'
             endif
             stop  
           endif
C
         epsxx = qx(1)*v(1)  + qx(2)*v(4)  + qx(3)*v(7)  + qx(4)*v(10)
     $         + qx(5)*v(13) + qx(6)*v(16) + qx(7)*v(19) + qx(8)*v(22)
         epsyy = qy(1)*v(2)  + qy(2)*v(5)  + qy(3)*v(8)  + qy(4)*v(11)
     $         + qy(5)*v(14) + qy(6)*v(17) + qy(7)*v(20) + qy(8)*v(23)
         epszz = qz(1)*v(3)  + qz(2)*v(6)  + qz(3)*v(9)  + qz(4)*v(12)
     $         + qz(5)*v(15) + qz(6)*v(18) + qz(7)*v(21) + qz(8)*v(24)
         gamxy = qy(1)*v(1)  + qy(2)*v(4)  + qy(3)*v(7)  + qy(4)*v(10)
     $         + qy(5)*v(13) + qy(6)*v(16) + qy(7)*v(19) + qy(8)*v(22)
     $         + qx(1)*v(2)  + qx(2)*v(5)  + qx(3)*v(8)  + qx(4)*v(11)
     $         + qx(5)*v(14) + qx(6)*v(17) + qx(7)*v(20) + qx(8)*v(23)
         gamyz = qy(1)*v(3)  + qy(2)*v(6)  + qy(3)*v(9)  + qy(4)*v(12)
     $         + qy(5)*v(15) + qy(6)*v(18) + qy(7)*v(21) + qy(8)*v(24)
     $         + qz(1)*v(2)  + qz(2)*v(5)  + qz(3)*v(8)  + qz(4)*v(11)
     $         + qz(5)*v(14) + qz(6)*v(17) + qz(7)*v(20) + qz(8)*v(23)
         gamxz = qz(1)*v(1)  + qz(2)*v(4)  + qz(3)*v(7)  + qz(4)*v(10)
     $         + qz(5)*v(13) + qz(6)*v(16) + qz(7)*v(19) + qz(8)*v(22)
     $         + qx(1)*v(3)  + qx(2)*v(6)  + qx(3)*v(9)  + qx(4)*v(12)
     $         + qx(5)*v(15) + qx(6)*v(18) + qx(7)*v(21) + qx(8)*v(24)
     
         depsxx=dqx(1)*v(1)  +dqx(2)*v(4)  +dqx(3)*v(7)  +dqx(4)*v(10)
     $         +dqx(5)*v(13) +dqx(6)*v(16) +dqx(7)*v(19) +dqx(8)*v(22)
         depsyy=dqy(1)*v(2)  +dqy(2)*v(5)  +dqy(3)*v(8)  +dqy(4)*v(11)
     $         +dqy(5)*v(14) +dqy(6)*v(17) +dqy(7)*v(20) +dqy(8)*v(23)
         depszz=dqz(1)*v(3)  +dqz(2)*v(6)  +dqz(3)*v(9)  +dqz(4)*v(12)
     $         +dqz(5)*v(15) +dqz(6)*v(18) +dqz(7)*v(21) +dqz(8)*v(24)
         dgamxy=dqy(1)*v(1)  +dqy(2)*v(4)  +dqy(3)*v(7)  +dqy(4)*v(10)
     $         +dqy(5)*v(13) +dqy(6)*v(16) +dqy(7)*v(19) +dqy(8)*v(22)
     $         +dqx(1)*v(2)  +dqx(2)*v(5)  +dqx(3)*v(8)  +dqx(4)*v(11)
     $         +dqx(5)*v(14) +dqx(6)*v(17) +dqx(7)*v(20) +dqx(8)*v(23)
         dgamyz=dqy(1)*v(3)  +dqy(2)*v(6)  +dqy(3)*v(9)  +dqy(4)*v(12)
     $         +dqy(5)*v(15) +dqy(6)*v(18) +dqy(7)*v(21) +dqy(8)*v(24)
     $         +dqz(1)*v(2)  +dqz(2)*v(5)  +dqz(3)*v(8)  +dqz(4)*v(11)
     $         +dqz(5)*v(14) +dqz(6)*v(17) +dqz(7)*v(20) +dqz(8)*v(23)
         dgamxz=dqz(1)*v(1)  +dqz(2)*v(4)  +dqz(3)*v(7)  +dqz(4)*v(10)
     $         +dqz(5)*v(13) +dqz(6)*v(16) +dqz(7)*v(19) +dqz(8)*v(22)
     $         +dqx(1)*v(3)  +dqx(2)*v(6)  +dqx(3)*v(9)  +dqx(4)*v(12)
     $         +dqx(5)*v(15) +dqx(6)*v(18) +dqx(7)*v(21) +dqx(8)*v(24)
     
         depsxx = depsxx +
     $          qx(1)*dv(1)  +qx(2)*dv(4)  +qx(3)*dv(7)  +qx(4)*dv(10)
     $         +qx(5)*dv(13) +qx(6)*dv(16) +qx(7)*dv(19) +qx(8)*dv(22)
         depsyy = depsyy +
     $          qy(1)*dv(2)  +qy(2)*dv(5)  +qy(3)*dv(8)  +qy(4)*dv(11)
     $         +qy(5)*dv(14) +qy(6)*dv(17) +qy(7)*dv(20) +qy(8)*dv(23)
         depszz = depszz +
     $          qz(1)*dv(3)  +qz(2)*dv(6)  +qz(3)*dv(9)  +qz(4)*dv(12)
     $         +qz(5)*dv(15) +qz(6)*dv(18) +qz(7)*dv(21) +qz(8)*dv(24)
         dgamxy = dgamxy +
     $          qy(1)*dv(1)  +qy(2)*dv(4)  +qy(3)*dv(7)  +qy(4)*dv(10)
     $         +qy(5)*dv(13) +qy(6)*dv(16) +qy(7)*dv(19) +qy(8)*dv(22)
     $         +qx(1)*dv(2)  +qx(2)*dv(5)  +qx(3)*dv(8)  +qx(4)*dv(11)
     $         +qx(5)*dv(14) +qx(6)*dv(17) +qx(7)*dv(20) +qx(8)*dv(23)
         dgamyz = dgamyz +
     $          qy(1)*dv(3)  +qy(2)*dv(6)  +qy(3)*dv(9)  +qy(4)*dv(12)
     $         +qy(5)*dv(15) +qy(6)*dv(18) +qy(7)*dv(21) +qy(8)*dv(24)
     $         +qz(1)*dv(2)  +qz(2)*dv(5)  +qz(3)*dv(8)  +qz(4)*dv(11)
     $         +qz(5)*dv(14) +qz(6)*dv(17) +qz(7)*dv(20) +qz(8)*dv(23)
         dgamxz = dgamxz +
     $          qz(1)*dv(1)  +qz(2)*dv(4)  +qz(3)*dv(7)  +qz(4)*dv(10)
     $         +qz(5)*dv(13) +qz(6)*dv(16) +qz(7)*dv(19) +qz(8)*dv(22)
     $         +qx(1)*dv(3)  +qx(2)*dv(6)  +qx(3)*dv(9)  +qx(4)*dv(12)
     $         +qx(5)*dv(15) +qx(6)*dv(18) +qx(7)*dv(21) +qx(8)*dv(24)
     
C
C.... ENGINEERING STRESS COMPUTATION
C
         stress(elm,1,n) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*epszz
     $             + c(1,4)*gamxy + c(1,5)*gamyz + c(1,6)*gamxz
         stress(elm,2,n) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*epszz
     $             + c(2,4)*gamxy + c(2,5)*gamyz + c(2,6)*gamxz
         stress(elm,3,n) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*epszz
     $             + c(3,4)*gamxy + c(3,5)*gamyz + c(3,6)*gamxz
         stress(elm,4,n) = c(4,1)*epsxx + c(4,2)*epsyy + c(4,3)*epszz
     $             + c(4,4)*gamxy + c(4,5)*gamyz + c(4,6)*gamxz
         stress(elm,5,n) = c(5,1)*epsxx + c(5,2)*epsyy + c(5,3)*epszz
     $             + c(5,4)*gamxy + c(5,5)*gamyz + c(5,6)*gamxz
         stress(elm,6,n) = c(6,1)*epsxx + c(6,2)*epsyy + c(6,3)*epszz
     $             + c(6,4)*gamxy + c(6,5)*gamyz + c(6,6)*gamxz
C
         dstress(elm,1,n) = dc(1,1)*epsxx 
     $    	   + dc(1,2)*epsyy + dc(1,3)*epszz
     $             + dc(1,4)*gamxy + dc(1,5)*gamyz + dc(1,6)*gamxz
         dstress(elm,2,n) =   dc(2,1)*epsxx 
     $    	   + dc(2,2)*epsyy + dc(2,3)*epszz
     $             + dc(2,4)*gamxy + dc(2,5)*gamyz + dc(2,6)*gamxz
         dstress(elm,3,n) =  dc(3,1)*epsxx 
     $    	   + dc(3,2)*epsyy + dc(3,3)*epszz
     $             + dc(3,4)*gamxy + dc(3,5)*gamyz + dc(3,6)*gamxz
         dstress(elm,4,n) =  dc(4,1)*epsxx 
     $    	   + dc(4,2)*epsyy + dc(4,3)*epszz
     $             + dc(4,4)*gamxy + dc(4,5)*gamyz + dc(4,6)*gamxz
         dstress(elm,5,n) =   dc(5,1)*epsxx 
     $    	   + dc(5,2)*epsyy + dc(5,3)*epszz
     $             + dc(5,4)*gamxy + dc(5,5)*gamyz + dc(5,6)*gamxz
         dstress(elm,6,n) =   dc(6,1)*epsxx
     $    	   + dc(6,2)*epsyy + dc(6,3)*epszz
     $             + dc(6,4)*gamxy + dc(6,5)*gamyz + dc(6,6)*gamxz
C
         dstress(elm,1,n) =  dstress(elm,1,n) + c(1,1)*depsxx 
     $    	   + c(1,2)*depsyy + c(1,3)*depszz
     $             + c(1,4)*dgamxy + c(1,5)*dgamyz + c(1,6)*dgamxz
         dstress(elm,2,n) = dstress(elm,2,n) + c(2,1)*depsxx 
     $    	   + c(2,2)*depsyy + c(2,3)*depszz
     $             + c(2,4)*dgamxy + c(2,5)*dgamyz + c(2,6)*dgamxz
         dstress(elm,3,n) =  dstress(elm,3,n) + c(3,1)*depsxx 
     $    	   + c(3,2)*depsyy + c(3,3)*depszz
     $             + c(3,4)*dgamxy + c(3,5)*dgamyz + c(3,6)*dgamxz
         dstress(elm,4,n) =  dstress(elm,4,n) + c(4,1)*depsxx 
     $    	   + c(4,2)*depsyy + c(4,3)*depszz
     $             + c(4,4)*dgamxy + c(4,5)*dgamyz + c(4,6)*dgamxz
         dstress(elm,5,n) =  dstress(elm,5,n) + c(5,1)*depsxx 
     $    	   + c(5,2)*depsyy + c(5,3)*depszz
     $             + c(5,4)*dgamxy + c(5,5)*dgamyz + c(5,6)*dgamxz
         dstress(elm,6,n) = dstress(elm,6,n) +  c(6,1)*depsxx 
     $    	   + c(6,2)*depsyy + c(6,3)*depszz
     $             + c(6,4)*dgamxy + c(6,5)*dgamyz + c(6,6)*dgamxz

C
 2000    continue
 
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
         do 1500  n = 1, 8
           xi  = xinod(n)
           eta = etanod(n)
           emu = emunod(n)
           call h8shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
C
           if (det .le. 0.0) then
             write(6,*)'Negative Jacobian determinant in g_sands17.f'
             if (det .eq. 0.0) then
               write(6,*)'Zero  jacobian determinant'
             endif
             stop  
           endif
C
         epsxx = qx(1)*v(1)  + qx(2)*v(4)  + qx(3)*v(7)  + qx(4)*v(10)
     $         + qx(5)*v(13) + qx(6)*v(16) + qx(7)*v(19) + qx(8)*v(22)
         epsyy = qy(1)*v(2)  + qy(2)*v(5)  + qy(3)*v(8)  + qy(4)*v(11)
     $         + qy(5)*v(14) + qy(6)*v(17) + qy(7)*v(20) + qy(8)*v(23)
         epszz = qz(1)*v(3)  + qz(2)*v(6)  + qz(3)*v(9)  + qz(4)*v(12)
     $         + qz(5)*v(15) + qz(6)*v(18) + qz(7)*v(21) + qz(8)*v(24)
         gamxy = qy(1)*v(1)  + qy(2)*v(4)  + qy(3)*v(7)  + qy(4)*v(10)
     $         + qy(5)*v(13) + qy(6)*v(16) + qy(7)*v(19) + qy(8)*v(22)
     $         + qx(1)*v(2)  + qx(2)*v(5)  + qx(3)*v(8)  + qx(4)*v(11)
     $         + qx(5)*v(14) + qx(6)*v(17) + qx(7)*v(20) + qx(8)*v(23)
         gamyz = qy(1)*v(3)  + qy(2)*v(6)  + qy(3)*v(9)  + qy(4)*v(12)
     $         + qy(5)*v(15) + qy(6)*v(18) + qy(7)*v(21) + qy(8)*v(24)
     $         + qz(1)*v(2)  + qz(2)*v(5)  + qz(3)*v(8)  + qz(4)*v(11)
     $         + qz(5)*v(14) + qz(6)*v(17) + qz(7)*v(20) + qz(8)*v(23)
         gamxz = qz(1)*v(1)  + qz(2)*v(4)  + qz(3)*v(7)  + qz(4)*v(10)
     $         + qz(5)*v(13) + qz(6)*v(16) + qz(7)*v(19) + qz(8)*v(22)
     $         + qx(1)*v(3)  + qx(2)*v(6)  + qx(3)*v(9)  + qx(4)*v(12)
     $         + qx(5)*v(15) + qx(6)*v(18) + qx(7)*v(21) + qx(8)*v(24)
     
         depsxx = 
     $          qx(1)*dv(1)  +qx(2)*dv(4)  +qx(3)*dv(7)  +qx(4)*dv(10)
     $         +qx(5)*dv(13) +qx(6)*dv(16) +qx(7)*dv(19) +qx(8)*dv(22)
         depsyy = 
     $          qy(1)*dv(2)  +qy(2)*dv(5)  +qy(3)*dv(8)  +qy(4)*dv(11)
     $         +qy(5)*dv(14) +qy(6)*dv(17) +qy(7)*dv(20) +qy(8)*dv(23)
         depszz = 
     $          qz(1)*dv(3)  +qz(2)*dv(6)  +qz(3)*dv(9)  +qz(4)*dv(12)
     $         +qz(5)*dv(15) +qz(6)*dv(18) +qz(7)*dv(21) +qz(8)*dv(24)
         dgamxy = 
     $          qy(1)*dv(1)  +qy(2)*dv(4)  +qy(3)*dv(7)  +qy(4)*dv(10)
     $         +qy(5)*dv(13) +qy(6)*dv(16) +qy(7)*dv(19) +qy(8)*dv(22)
     $         +qx(1)*dv(2)  +qx(2)*dv(5)  +qx(3)*dv(8)  +qx(4)*dv(11)
     $         +qx(5)*dv(14) +qx(6)*dv(17) +qx(7)*dv(20) +qx(8)*dv(23)
         dgamyz = 
     $          qy(1)*dv(3)  +qy(2)*dv(6)  +qy(3)*dv(9)  +qy(4)*dv(12)
     $         +qy(5)*dv(15) +qy(6)*dv(18) +qy(7)*dv(21) +qy(8)*dv(24)
     $         +qz(1)*dv(2)  +qz(2)*dv(5)  +qz(3)*dv(8)  +qz(4)*dv(11)
     $         +qz(5)*dv(14) +qz(6)*dv(17) +qz(7)*dv(20) +qz(8)*dv(23)
         dgamxz =
     $          qz(1)*dv(1)  +qz(2)*dv(4)  +qz(3)*dv(7)  +qz(4)*dv(10)
     $         +qz(5)*dv(13) +qz(6)*dv(16) +qz(7)*dv(19) +qz(8)*dv(22)
     $         +qx(1)*dv(3)  +qx(2)*dv(6)  +qx(3)*dv(9)  +qx(4)*dv(12)
     $         +qx(5)*dv(15) +qx(6)*dv(18) +qx(7)*dv(21) +qx(8)*dv(24)

C     ! strain field is only named stress field

         stress(elm,1,n) = epsxx
         stress(elm,2,n) = epsyy 
         stress(elm,3,n) = epszz 
         stress(elm,4,n) = gamxy 
         stress(elm,5,n) = gamyz 
         stress(elm,6,n) = gamxz 
	 
         dstress(elm,1,n) = depsxx
         dstress(elm,2,n) = depsyy
         dstress(elm,3,n) = depszz
         dstress(elm,4,n) = dgamxy
         dstress(elm,5,n) = dgamyz
         dstress(elm,6,n) = dgamxz
C
 1500    continue

C                                                          cvar.ne.0       
      else
      
         do 2500  n = 1, 8
           xi  = xinod(n)
           eta = etanod(n)
           emu = emunod(n)
           call gxh8shpe(xi,eta,emu,x,y,z,dx,dy,dz,q,qx,qy,qz,
     *                    dqx,dqy,dqz,det,ddet)
C
           if (det .le. 0.0) then
             write(6,*)'Negative Jacobian determinant in g_sands17.f'
             if (det .eq. 0.0) then
               write(6,*)'Zero  jacobian determinant'
             endif
             stop  
           endif
C
         epsxx = qx(1)*v(1)  + qx(2)*v(4)  + qx(3)*v(7)  + qx(4)*v(10)
     $         + qx(5)*v(13) + qx(6)*v(16) + qx(7)*v(19) + qx(8)*v(22)
         epsyy = qy(1)*v(2)  + qy(2)*v(5)  + qy(3)*v(8)  + qy(4)*v(11)
     $         + qy(5)*v(14) + qy(6)*v(17) + qy(7)*v(20) + qy(8)*v(23)
         epszz = qz(1)*v(3)  + qz(2)*v(6)  + qz(3)*v(9)  + qz(4)*v(12)
     $         + qz(5)*v(15) + qz(6)*v(18) + qz(7)*v(21) + qz(8)*v(24)
         gamxy = qy(1)*v(1)  + qy(2)*v(4)  + qy(3)*v(7)  + qy(4)*v(10)
     $         + qy(5)*v(13) + qy(6)*v(16) + qy(7)*v(19) + qy(8)*v(22)
     $         + qx(1)*v(2)  + qx(2)*v(5)  + qx(3)*v(8)  + qx(4)*v(11)
     $         + qx(5)*v(14) + qx(6)*v(17) + qx(7)*v(20) + qx(8)*v(23)
         gamyz = qy(1)*v(3)  + qy(2)*v(6)  + qy(3)*v(9)  + qy(4)*v(12)
     $         + qy(5)*v(15) + qy(6)*v(18) + qy(7)*v(21) + qy(8)*v(24)
     $         + qz(1)*v(2)  + qz(2)*v(5)  + qz(3)*v(8)  + qz(4)*v(11)
     $         + qz(5)*v(14) + qz(6)*v(17) + qz(7)*v(20) + qz(8)*v(23)
         gamxz = qz(1)*v(1)  + qz(2)*v(4)  + qz(3)*v(7)  + qz(4)*v(10)
     $         + qz(5)*v(13) + qz(6)*v(16) + qz(7)*v(19) + qz(8)*v(22)
     $         + qx(1)*v(3)  + qx(2)*v(6)  + qx(3)*v(9)  + qx(4)*v(12)
     $         + qx(5)*v(15) + qx(6)*v(18) + qx(7)*v(21) + qx(8)*v(24)
     
         depsxx=dqx(1)*v(1)  +dqx(2)*v(4)  +dqx(3)*v(7)  +dqx(4)*v(10)
     $         +dqx(5)*v(13) +dqx(6)*v(16) +dqx(7)*v(19) +dqx(8)*v(22)
         depsyy=dqy(1)*v(2)  +dqy(2)*v(5)  +dqy(3)*v(8)  +dqy(4)*v(11)
     $         +dqy(5)*v(14) +dqy(6)*v(17) +dqy(7)*v(20) +dqy(8)*v(23)
         depszz=dqz(1)*v(3)  +dqz(2)*v(6)  +dqz(3)*v(9)  +dqz(4)*v(12)
     $         +dqz(5)*v(15) +dqz(6)*v(18) +dqz(7)*v(21) +dqz(8)*v(24)
         dgamxy=dqy(1)*v(1)  +dqy(2)*v(4)  +dqy(3)*v(7)  +dqy(4)*v(10)
     $         +dqy(5)*v(13) +dqy(6)*v(16) +dqy(7)*v(19) +dqy(8)*v(22)
     $         +dqx(1)*v(2)  +dqx(2)*v(5)  +dqx(3)*v(8)  +dqx(4)*v(11)
     $         +dqx(5)*v(14) +dqx(6)*v(17) +dqx(7)*v(20) +dqx(8)*v(23)
         dgamyz=dqy(1)*v(3)  +dqy(2)*v(6)  +dqy(3)*v(9)  +dqy(4)*v(12)
     $         +dqy(5)*v(15) +dqy(6)*v(18) +dqy(7)*v(21) +dqy(8)*v(24)
     $         +dqz(1)*v(2)  +dqz(2)*v(5)  +dqz(3)*v(8)  +dqz(4)*v(11)
     $         +dqz(5)*v(14) +dqz(6)*v(17) +dqz(7)*v(20) +dqz(8)*v(23)
         dgamxz=dqz(1)*v(1)  +dqz(2)*v(4)  +dqz(3)*v(7)  +dqz(4)*v(10)
     $         +dqz(5)*v(13) +dqz(6)*v(16) +dqz(7)*v(19) +dqz(8)*v(22)
     $         +dqx(1)*v(3)  +dqx(2)*v(6)  +dqx(3)*v(9)  +dqx(4)*v(12)
     $         +dqx(5)*v(15) +dqx(6)*v(18) +dqx(7)*v(21) +dqx(8)*v(24)
     
         depsxx = depsxx +
     $          qx(1)*dv(1)  +qx(2)*dv(4)  +qx(3)*dv(7)  +qx(4)*dv(10)
     $         +qx(5)*dv(13) +qx(6)*dv(16) +qx(7)*dv(19) +qx(8)*dv(22)
         depsyy = depsyy +
     $          qy(1)*dv(2)  +qy(2)*dv(5)  +qy(3)*dv(8)  +qy(4)*dv(11)
     $         +qy(5)*dv(14) +qy(6)*dv(17) +qy(7)*dv(20) +qy(8)*dv(23)
         depszz = depszz +
     $          qz(1)*dv(3)  +qz(2)*dv(6)  +qz(3)*dv(9)  +qz(4)*dv(12)
     $         +qz(5)*dv(15) +qz(6)*dv(18) +qz(7)*dv(21) +qz(8)*dv(24)
         dgamxy = dgamxy +
     $          qy(1)*dv(1)  +qy(2)*dv(4)  +qy(3)*dv(7)  +qy(4)*dv(10)
     $         +qy(5)*dv(13) +qy(6)*dv(16) +qy(7)*dv(19) +qy(8)*dv(22)
     $         +qx(1)*dv(2)  +qx(2)*dv(5)  +qx(3)*dv(8)  +qx(4)*dv(11)
     $         +qx(5)*dv(14) +qx(6)*dv(17) +qx(7)*dv(20) +qx(8)*dv(23)
         dgamyz = dgamyz +
     $          qy(1)*dv(3)  +qy(2)*dv(6)  +qy(3)*dv(9)  +qy(4)*dv(12)
     $         +qy(5)*dv(15) +qy(6)*dv(18) +qy(7)*dv(21) +qy(8)*dv(24)
     $         +qz(1)*dv(2)  +qz(2)*dv(5)  +qz(3)*dv(8)  +qz(4)*dv(11)
     $         +qz(5)*dv(14) +qz(6)*dv(17) +qz(7)*dv(20) +qz(8)*dv(23)
         dgamxz = dgamxz +
     $          qz(1)*dv(1)  +qz(2)*dv(4)  +qz(3)*dv(7)  +qz(4)*dv(10)
     $         +qz(5)*dv(13) +qz(6)*dv(16) +qz(7)*dv(19) +qz(8)*dv(22)
     $         +qx(1)*dv(3)  +qx(2)*dv(6)  +qx(3)*dv(9)  +qx(4)*dv(12)
     $         +qx(5)*dv(15) +qx(6)*dv(18) +qx(7)*dv(21) +qx(8)*dv(24)

C     ! strain field is only named stress field

         stress(elm,1,n) = epsxx
         stress(elm,2,n) = epsyy 
         stress(elm,3,n) = epszz 
         stress(elm,4,n) = gamxy 
         stress(elm,5,n) = gamyz 
         stress(elm,6,n) = gamxz 
	 
         dstress(elm,1,n) = depsxx
         dstress(elm,2,n) = depsyy
         dstress(elm,3,n) = depszz
         dstress(elm,4,n) = dgamxy
         dstress(elm,5,n) = dgamyz
         dstress(elm,6,n) = dgamxz
C
 2500    continue
 
      endif
      
C.... COMPUTE EQUIVALENT STRAIN
C
	if(strainFlg) then
	  call gxstrainvm(stress,dstress,maxgus,maxstr,maxsze,8)
	endif
      
      endif
C
      return
      end
