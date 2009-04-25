      subroutine gxq4dmas(p, rho, g_rho, g_elmass,
     * h, g_h, x, g_x, y, g_y, gamma, grvfor, g_grvfor,  
     * grvflg, totmas, g_totmas, masflg, mratio)
C**************************************************************
C                                                             *
C              Derivative of q4dmas routine                   *
C                                                             *
C**************************************************************
C
        implicit none
C
C..... DECLARE GLOBAL VARIABLES
C
C.... INTEGER CONSTANTS
C
        integer p, mxnseq
C
C.... REAL CONSTANTS
C
        real*8 rho, totmas, mratio
C
C.... REAL ARRAYS
C
        real*8 g_elmass(8,8), h(*), x(*), y(*)
        real*8 gamma(*), grvfor(*)
C
C.... LOGICAL CONSTANTS
C
        logical grvflg, masflg
C
C.... DECLARE LOCAL VARIABLES FOR Q4DMAS
C
        integer k,l,ls(8),i,j,ix,iy,jx,jy
        real*8  mass, c1, xi, eta, weight, det,diag,gam,w
        real*8  q(4), qx(4), qy(4)
C
        double precision d6_b, d5_b, d4_b, d3_b, d2_b, d3_v
	double precision g_det, g_c1, g_rho, g_totmas, g_mass 
        double precision g_h(*), elmass(8,8),g_grvfor(*) 
	double precision g_x(*), g_y(*), g_diag, g_gam, g_w
	double precision g_qx(4),g_qy(4)
C
        data ls /1,3,5,7,2,4,6,8/
C
C.... INITIALIZE MASS MATRIX
C
        do i=1,8
	  do j=1,8
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
            call qgauss(p, k, p, l, xi, eta, weight)
            call gxq4shpe(xi,eta,x,
     &                    g_x,y,g_y,q,qx,g_qx,qy,g_qy,det,g_det)
C
            w = weight * det * rho *
     $          (h(1)*q(1)+h(2)*q(2)+h(3)*q(3)+h(4)*q(4))
C
            g_w = weight * (g_det * rho + det * g_rho) *
     $            (h(1)*q(1)+h(2)*q(2)+h(3)*q(3)+h(4)*q(4))
     $          + weight * det * rho *
     $            (g_h(1)*q(1)+g_h(2)*q(2)+g_h(3)*q(3)+g_h(4)*q(4))
C
	      mass = mass + w
	    g_mass = g_mass + g_w
C
            do j = 1,4
              jx =    ls(j)
              jy =    ls(j+4)
              do i = j,4
                ix =     ls(i)
                iy =     ls(i+4)
                elmass(jx,ix) =  elmass(jx,ix) + q(j)*q(i)*w
                elmass(ix,jx) =  elmass(jx,ix)
                elmass(jy,iy) =  elmass(jx,ix) 
                elmass(iy,jy) =  elmass(jx,ix) 
C
                g_elmass(jx,ix) =  g_elmass(jx,ix) + q(j)*q(i)*g_w 
                g_elmass(ix,jx) =  g_elmass(jx,ix)
                g_elmass(jy,iy) =  g_elmass(jx,ix) 
                g_elmass(iy,jy) =  g_elmass(jx,ix) 
              enddo
	    enddo
20	  continue
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
          do i=1,8
	    diag   =   diag +   elmass(i,i)
	    g_diag = g_diag + g_elmass(i,i)
	    do j=1,8
	        gam =   gam +   elmass(i,j)
	      g_gam = g_gam + g_elmass(i,j)
            enddo
	  enddo
C	
C	
	  g_gam = (g_gam - gam*g_diag/diag)/diag;
	    gam = gam/diag;
C	
          do i=1,8
            g_diag = g_elmass(i,i) * gam + elmass(i,i) * g_gam
	    do j=i+1,8
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
