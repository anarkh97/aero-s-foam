	subroutine gxbr8mas(p,rho,g_rho,g_elmass,x,y,z,g_x,g_y,g_z,
     .                      gamma,grvfor,g_grvfor,
     .                      grvflg,g_totmas,masflg,mratio)
***************************************************************
* THIS SUBROUTINE COMPUTE THE ELEMENT MASS MATRIX FOR THE 8-  *
* NODE BRICK.                                                 *
*                                                             *
***************************************************************
*                                                             *
*		VARIABLES 				      *
*                                                             *
*	     P = PxPxP GAUSS QUADRATURE RULE                  *
*	   RHO = ELEMENT DENSITY                              *
*	ELMASS = THE ELEMENT MASS MATRIX                      *
*	MXNSEQ = LEADING DIMENSION OF ELMASS                  *
*	     X = X COORDINATE ARRAY                           *
*	     Y = Y COORDINATE ARRAY                           *
*	     Z = Z COORDINATE ARRAY                           *
*       GAMMA = ACCELERATION DUE TO GRAVITY VECTOR            *
*      GRVFOR = FORCE DUE TO GRAVITY VECTOR                   *
*      GRVFLG = LOGICAL FLAG FOR COMPUTING GRAVITY FORCE      *
*
***************************************************************
	implicit none
	integer p

	real*8 rho, g_rho, g_totmas, mratio
	real*8 z(*),x(*),y(*)
	real*8 gamma(*),grvfor(*),g_grvfor(*)
	real*8 g_z(*),g_x(*),g_y(*),g_elmass(24,24)
	logical grvflg, masflg

	integer k,l,jj,i,j,ix,iy,iz,jx,jy,jz
	real*8 q(8),qx(8),qy(8),qz(8)
	real*8 ls(24), elmass(24,24)
	real*8 g_qx(8),g_qy(8),g_qz(8),det,g_det
	real*8 weight, mass, g_mass, w, g_w
	real*8 gam, g_gam, diag, g_diag
	real*8 xi,eta,emu

	data ls /1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,
     $           23,3,6,9,12,15,18,21,24/

C
C.... INITIALIZE MASS MATRIX
C
        do i=1,24
	  do j=1,24
   	      elmass(i,j) = 0.0d0
   	    g_elmass(i,j) = 0.0d0
          enddo
        enddo
C
C.... GAUSS QUADRATURE
C
	mass   = 0.0d0
	g_mass = 0.0d0
C	
        do 10 k = 1, p
	   do 20 l = 1, p
	      do 30 jj = 1, p
		 call hxgaus(p,k,p,l,p,jj,xi,eta,emu,weight)
		 call h8shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
		 call gxh8shpe(xi,eta,emu,x,y,z,g_x,g_y,g_z,
     $	                       q,qx,qy,qz,g_qx,g_qy,g_qz,det,g_det)
C
		 w = weight * det * rho *
     $               (q(1)+q(2)+q(3)+q(4)+q(5)+q(6)+q(7)+q(8))
C
		 g_w = weight * (g_det * rho + det * g_rho) *
     $               (q(1)+q(2)+q(3)+q(4)+q(5)+q(6)+q(7)+q(8))
C
		 mass   = mass + w
		 g_mass = g_mass + g_w
C
		 do j = 1,8
		    jx = ls(j)
		    jy = ls(j+8)
		    jz = ls(j+16)
		    do i = j,8
		       ix = ls(i)
		       iy = ls(i+8)
		       iz = ls(i+16)
		       elmass(jx,ix) =  elmass(jx,ix) + q(j)*q(i)*w
		       elmass(ix,jx) =  elmass(jx,ix)
		       elmass(jy,iy) =  elmass(jx,ix) 
		       elmass(iy,jy) =  elmass(jx,ix) 
		       elmass(jz,iz) =  elmass(jx,ix) 
		       elmass(iz,jz) =  elmass(jx,ix) 
C       
		       g_elmass(jx,ix) =  g_elmass(jx,ix) + 
     $		                          q(j)*q(i)*g_w 
		       g_elmass(ix,jx) =  g_elmass(jx,ix)
		       g_elmass(jy,iy) =  g_elmass(jx,ix) 
		       g_elmass(iy,jy) =  g_elmass(jx,ix) 
		       g_elmass(jz,iz) =  g_elmass(jx,ix) 
		       g_elmass(iz,jz) =  g_elmass(jx,ix) 
		    enddo
		 enddo
 30	      continue
 20	   continue
 10	continue
C
C.... CREATE LUMPED MATRIX
C
        if (mratio.gt.0) then
C
            gam  = 0.0d0
          g_gam  = 0.0d0
	    diag = 0.0d0
	  g_diag = 0.0d0
C	
          do i=1,24
	    diag   =   diag +   elmass(i,i)
	    g_diag = g_diag + g_elmass(i,i)
	    do j=1,24
	        gam =   gam +   elmass(i,j)
	      g_gam = g_gam + g_elmass(i,j)
            enddo
	  enddo
C	
C	
	  g_gam = (g_gam - gam*g_diag/diag)/diag;
	    gam = gam/diag;
C	
          do i=1,24
            g_diag = g_elmass(i,i) * gam + elmass(i,i) * g_gam
	    do j=i+1,24
	      g_elmass(i,j) = (1.0d0-mratio)*g_elmass(i,j)
	      g_elmass(j,i) = g_elmass(i,j)
            enddo
            g_elmass(i,i) = mratio*g_diag + (1.0d0-mratio)*g_elmass(i,i)
	  enddo
C
        endif
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
        if (grvflg) then
          g_grvfor(1) = g_mass*gamma(1)
          g_grvfor(2) = g_mass*gamma(2)
        endif
C
C.... COMPUTE THE SUBDOMAIN TOTAL MASS
C
	if (masflg) then
	  g_totmas = g_totmas + g_mass
	endif
C
	return
	end
