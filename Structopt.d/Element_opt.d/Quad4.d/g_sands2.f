      subroutine gxsands2(escm, x, g_x, y, g_y,  c, g_c,  
     &                    v, g_v,  stress, g_stress,  strain, g_strain,
     &                    maxgus, maxstr, elm, msize, vmflg, strainflg,
     &                    emod, g_emod, nu, g_nu, alpha, g_alpha, tref, 
     &                    ndtemps, g_ndtemps, rip)
C
C
      implicit none
C
C                   A R G U M E N T S
C
      integer           elm,msize,maxgus,maxstr
      character*(*)     escm
C
      real*8  x(4), y(4), c(3,3), v(8)
      real*8  g_x(4), g_y(4), g_c(3,3), g_v(8)
C
      real*8  stress(msize,maxstr,maxgus),strain(msize,maxstr,maxgus) 
      real*8  g_stress(msize,maxstr,maxgus)
      real*8  g_strain(msize,maxstr,maxgus) 
C
      real*8  emod, nu, alpha, ndtemps(4), tref, rip
      real*8  g_emod, g_nu, g_alpha, g_ndtemps(4)
C
      logical vmflg,strainFlg
C
C                   L O C A L   V A R I A B L E S
C
        real*8 q(4), qx(4), qy(4), g_qx(4), g_qy(4)
        real*8 xinod(4), etanod(4), cext(4,4), xi, eta
	real*8 sigauss(4), g_sigauss(4)
        real*8 det, epsxx, epsyy, epszz, gamxy
        real*8 g_det, g_epsxx, g_epsyy, g_epszz, g_gamxy
        real*8 tl(4), tgp, tc, eptxo, eptyo
        real*8 g_tl(4), g_tgp, g_tc, g_eptxo, g_eptyo
        integer i, n, ids(3)
C
C                   D A T A
C
      data xinod  /-1.0, 1.0, 1.0,-1.0/
      data etanod /-1.0,-1.0, 1.0, 1.0/
      data cext / 1.866025404,        -0.5, 0.133974596,	-0.5,
     $     		 -0.5, 1.866025404,	   -0.5, 0.133974596,
     $     	  0.133974596,        -0.5, 1.866025404,        -0.5,
     $     	         -0.5, 0.133974596,        -0.5,  1.866025404 /
      data ids /1,2,4/
C
C
C
C.... COMPUTE THE THERMAL FIELD
C
        do i = 1, 4
            tl(i) =   ndtemps(i) - tref
          g_tl(i) = g_ndtemps(i)
        enddo
c
c...  compute stress correction due to thermal strains
c      
      if (rip.eq.0) then
          tc = emod*alpha/(1.0d0-nu)
        g_tc = (g_emod*alpha+emod*g_alpha+emod*alpha/(1.0d0-nu))
     &       / (1.0d0-nu)
      else
          tc = emod*alpha/(1.0d0-2.0d0*nu)
        g_tc = (g_emod*alpha+emod*g_alpha
     &       +  2.0d0*emod*alpha/(1.0d0-2.0d0*nu))
     &       / (1.0d0-2.0d0*nu)
      endif	
C
C     stress evaluation at Gauss points
C
      if (escm(1:1) .eq. 'D') then
C
C     loop over all Gauss points
C      
        do 2000  n = 1,4
C
          xi  = xinod (n)
          eta = etanod(n)
C
          call gxq4shpe(xi, eta, x, g_x, y, g_y, q,
     &                  qx, g_qx, qy, g_qy, det, g_det)
C
C.... COMPUTE THE THERMAL STRESS
C
            tgp = q(1)*tl(1)+q(2)*tl(2)+q(3)*tl(3)+q(4)*tl(4)
          g_tgp = q(1)*g_tl(1)+q(2)*g_tl(2)+q(3)*g_tl(3)+q(4)*g_tl(4)
C
            eptxo = tc*tgp
          g_eptxo = g_tc*tgp + tc*g_tgp
            eptyo = tc*tgp
          g_eptyo = g_tc*tgp + tc*g_tgp
C
C.... COMPUTE THE TOTAL STRAIN FROM THE COUPLED SOLVE FOR DISPLACEMENTS
C
            epsxx = qx(1)*v(1) + qx(2)*v(3) + qx(3)*v(5) + qx(4)*v(7)
          g_epsxx = g_qx(1)*v(1) + g_qx(2)*v(3) 
     &            + g_qx(3)*v(5) + g_qx(4)*v(7)
     &            + qx(1)*g_v(1) + qx(2)*g_v(3) 
     &            + qx(3)*g_v(5) + qx(4)*g_v(7)
C     
            epsyy = qy(1)*v(2) + qy(2)*v(4) + qy(3)*v(6) + qy(4)*v(8)
          g_epsyy = g_qy(1)*v(2) + g_qy(2)*v(4) 
     &            + g_qy(3)*v(6) + g_qy(4)*v(8)
     &            + qy(1)*g_v(2) + qy(2)*g_v(4) 
     &            + qy(3)*g_v(6) + qy(4)*g_v(8)
C
            gamxy = qy(1)*v(1) + qy(2)*v(3) + qy(3)*v(5) + qy(4)*v(7)
     &            + qx(1)*v(2) + qx(2)*v(4) + qx(3)*v(6) + qx(4)*v(8)
          g_gamxy = g_qy(1)*v(1) + g_qy(2)*v(3) 
     &            + g_qy(3)*v(5) + g_qy(4)*v(7)
     &            + g_qx(1)*v(2) + g_qx(2)*v(4) 
     &            + g_qx(3)*v(6) + g_qx(4)*v(8)
     &            + qy(1)*g_v(1) + qy(2)*g_v(3) 
     &            + qy(3)*g_v(5) + qy(4)*g_v(7)
     &            + qx(1)*g_v(2) + qx(2)*g_v(4) 
     &            + qx(3)*g_v(6) + qx(4)*g_v(8)
C
C.... COMPUTE THE TOTAL STRESS
C
            stress(elm,1,n) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*gamxy
     &                      - eptxo
          g_stress(elm,1,n) = g_c(1,1)*epsxx + g_c(1,2)*epsyy 
     &                      + g_c(1,3)*gamxy + c(1,1)*g_epsxx
     &                      + c(1,2)*g_epsyy + c(1,3)*g_gamxy
     &                      - g_eptxo
C
            stress(elm,2,n) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*gamxy
     &                      - eptyo
          g_stress(elm,2,n) = g_c(2,1)*epsxx + g_c(2,2)*epsyy 
     &                      + g_c(2,3)*gamxy + c(2,1)*g_epsxx
     &                      + c(2,2)*g_epsyy + c(2,3)*g_gamxy
     &                      - g_eptyo
C
            stress(elm,4,n) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*gamxy
          g_stress(elm,4,n) = g_c(3,1)*epsxx + g_c(3,2)*epsyy 
     &                      + g_c(3,3)*gamxy + c(3,1)*g_epsxx 
     &                      + c(3,2)*g_epsyy + c(3,3)*g_gamxy
C
          strain(elm,1,n) = epsxx 
          strain(elm,2,n) = epsyy 
          strain(elm,4,n) = gamxy 
C
          g_strain(elm,1,n) = g_epsxx 
          g_strain(elm,2,n) = g_epsyy 
          g_strain(elm,4,n) = g_gamxy 
C
          if (rip.eq.0) then
              stress(elm,3,n) = 0.0d0
            g_stress(elm,3,n) = 0.0d0
              strain(elm,3,n) = -nu/emod*(stress(elm,1,n)
     &                        + stress(elm,2,n)) + alpha*tgp 
            g_strain(elm,3,n) = (-g_nu*(stress(elm,1,n)+stress(elm,2,n)) 
     &	                      - nu*(g_stress(elm,1,n)+g_stress(elm,2,n)) 
     &                        + nu*(stress(elm,1,n)+stress(elm,2,n))
     &                        * g_emod/emod)/emod 
     &                        + g_alpha*tgp + alpha*g_tgp
          else
	      stress(elm,3,n) = c(1,2)*epsxx+c(1,2)*epsyy+c(1,3)*gamxy
     &                        - eptxo
	    g_stress(elm,3,n) = g_c(1,2)*epsxx + g_c(1,2)*epsyy 
     &                        + g_c(1,3)*gamxy + c(1,2)*g_epsxx
     &                        + c(1,2)*g_epsyy + c(1,3)*g_gamxy
     &                        - g_eptxo
	      strain(elm,3,n) = 0.0d0	    
	    g_strain(elm,3,n) = 0.0d0	    
          endif
	  
 2000   continue
C
C.... EXTRAPOLATE FROM THE GAUSS POINTS
c
      else

        do 2200  n = 1,4
          stress(elm,1,n) = 0.0d0
          stress(elm,2,n) = 0.0d0
          stress(elm,3,n) = 0.0d0
          stress(elm,4,n) = 0.0d0
C
          g_stress(elm,1,n) = 0.0d0
          g_stress(elm,2,n) = 0.0d0
          g_stress(elm,3,n) = 0.0d0
          g_stress(elm,4,n) = 0.0d0
C
          strain(elm,1,n) = 0.0d0    
          strain(elm,2,n) = 0.0d0    
          strain(elm,3,n) = 0.0d0    
          strain(elm,4,n) = 0.0d0    
C
          g_strain(elm,1,n) = 0.0d0    
          g_strain(elm,2,n) = 0.0d0    
          g_strain(elm,3,n) = 0.0d0    
          g_strain(elm,4,n) = 0.0d0    
 2200   continue
C
        do 3000  i = 1,4
          xi  =     xinod(i)*0.577350269
          eta =    etanod(i)*0.577350269
C
          call gxq4shpe(xi, eta, x, g_x, y, g_y, q,
     &                  qx, g_qx, qy, g_qy, det, g_det)
C
C.... COMPUTE THE THERMAL STRESS
C
            tgp = q(1)*tl(1)+q(2)*tl(2)+q(3)*tl(3)+q(4)*tl(4)
          g_tgp = q(1)*g_tl(1)+q(2)*g_tl(2)+q(3)*g_tl(3)+q(4)*g_tl(4)
C
            eptxo = tc*tgp
          g_eptxo = g_tc*tgp + tc*g_tgp
            eptyo = tc*tgp
          g_eptyo = g_tc*tgp + tc*g_tgp
C
C.... COMPUTE THE TOTAL STRAIN FROM THE COUPLED SOLVE FOR DISPLACEMENTS
C
            epsxx = qx(1)*v(1) + qx(2)*v(3) + qx(3)*v(5) + qx(4)*v(7)
          g_epsxx = g_qx(1)*v(1) + g_qx(2)*v(3) 
     &            + g_qx(3)*v(5) + g_qx(4)*v(7)
     &            + qx(1)*g_v(1) + qx(2)*g_v(3) 
     &            + qx(3)*g_v(5) + qx(4)*g_v(7)
C     
            epsyy = qy(1)*v(2) + qy(2)*v(4) + qy(3)*v(6) + qy(4)*v(8)
          g_epsyy = g_qy(1)*v(2) + g_qy(2)*v(4) 
     &            + g_qy(3)*v(6) + g_qy(4)*v(8)
     &            + qy(1)*g_v(2) + qy(2)*g_v(4) 
     &            + qy(3)*g_v(6) + qy(4)*g_v(8)
C
            gamxy = qy(1)*v(1) + qy(2)*v(3) + qy(3)*v(5) + qy(4)*v(7)
     &            + qx(1)*v(2) + qx(2)*v(4) + qx(3)*v(6) + qx(4)*v(8)
          g_gamxy = g_qy(1)*v(1) + g_qy(2)*v(3) 
     &            + g_qy(3)*v(5) + g_qy(4)*v(7)
     &            + g_qx(1)*v(2) + g_qx(2)*v(4) 
     &            + g_qx(3)*v(6) + g_qx(4)*v(8)
     &            + qy(1)*g_v(1) + qy(2)*g_v(3) 
     &            + qy(3)*g_v(5) + qy(4)*g_v(7)
     &            + qx(1)*g_v(2) + qx(2)*g_v(4) 
     &            + qx(3)*g_v(6) + qx(4)*g_v(8)
C
C.... COMPUTE THE TOTAL STRESS
C
            sigauss(1) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*gamxy
     &                 - eptxo
          g_sigauss(1) = g_c(1,1)*epsxx + g_c(1,2)*epsyy 
     &                 + g_c(1,3)*gamxy + c(1,1)*g_epsxx
     &                 + c(1,2)*g_epsyy + c(1,3)*g_gamxy
     &                 - g_eptxo
C
            sigauss(2) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*gamxy
     &                 - eptyo
          g_sigauss(2) = g_c(2,1)*epsxx + g_c(2,2)*epsyy 
     &                 + g_c(2,3)*gamxy + c(2,1)*g_epsxx
     &                 + c(2,2)*g_epsyy + c(2,3)*g_gamxy
     &                 - g_eptyo
C
            sigauss(4) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*gamxy
          g_sigauss(4) = g_c(3,1)*epsxx + g_c(3,2)*epsyy 
     &                 + g_c(3,3)*gamxy + c(3,1)*g_epsxx 
     &                 + c(3,2)*g_epsyy + c(3,3)*g_gamxy
C
          if (rip.eq.0) then
              sigauss(3) = 0.0d0
            g_sigauss(3) = 0.0d0
              epszz = -nu/emod*(stress(elm,1,n)
     &              + stress(elm,2,n)) + alpha*tgp 
            g_epszz = (-g_nu*(stress(elm,1,n)+stress(elm,2,n)) 
     &	            - nu*(g_stress(elm,1,n)+g_stress(elm,2,n)) 
     &              + nu*(stress(elm,1,n)+stress(elm,2,n))
     &              * g_emod/emod)/emod 
     &              + g_alpha*tgp + alpha*g_tgp
          else
	      sigauss(3) = c(1,2)*epsxx+c(1,2)*epsyy+c(1,3)*gamxy
     &                   - eptxo
	    g_sigauss(3) = g_c(1,2)*epsxx + g_c(1,2)*epsyy 
     &                   + g_c(1,3)*gamxy + c(1,2)*g_epsxx
     &                   + c(1,2)*g_epsyy + c(1,3)*g_gamxy
     &                   - g_eptxo
	      epszz = 0.0d0	    
	    g_epszz = 0.0d0	    
          endif
C
          do 2500  n = 1,4
              stress(elm,1,n)=  stress(elm,1,n)+cext(i,n)*  sigauss(1)
            g_stress(elm,1,n)=g_stress(elm,1,n)+cext(i,n)*g_sigauss(1)
              strain(elm,1,n)=  epsxx 
            g_strain(elm,1,n)=g_epsxx 
              stress(elm,2,n)=  stress(elm,2,n)+cext(i,n)*  sigauss(2)
            g_stress(elm,2,n)=g_stress(elm,2,n)+cext(i,n)*g_sigauss(2)
              strain(elm,2,n)=  epsyy 
            g_strain(elm,2,n)=g_epsyy 
              stress(elm,3,n)=  stress(elm,3,n)+cext(i,n)*  sigauss(3)
            g_stress(elm,3,n)=g_stress(elm,3,n)+cext(i,n)*g_sigauss(3)
              strain(elm,3,n)=  epszz 
            g_strain(elm,3,n)=g_epszz 
              stress(elm,4,n)=  stress(elm,4,n)+cext(i,n)*  sigauss(4)
            g_stress(elm,4,n)=g_stress(elm,4,n)+cext(i,n)*g_sigauss(4)
              strain(elm,4,n)=  gamxy 
            g_strain(elm,4,n)=g_gamxy 
 2500       continue
 3000     continue
      end if
C
C     set remain stress/strain components to zero
C
      do  n = 1,4
        do i = 5,maxstr
            stress(elm,i,n) = 0.0d0
          g_stress(elm,i,n) = 0.0d0
            strain(elm,i,n) = 0.0 d0    
          g_strain(elm,i,n) = 0.0d0    
        enddo
      enddo
C
C.... COMPUTE THE VON MISES STRESS IN THE PLATE
C
      if (vmflg) then
        call gxvmelmv(stress, g_stress, maxgus, maxstr, msize, elm, 4)
      endif
C
C.... COMPUTE THE VON MISES STRAIN IN THE PLATE
C
      if (strainFlg) then
        call gxstrainvm(strain, g_strain, maxgus, maxstr, msize, 4)
      endif

      return
      end
