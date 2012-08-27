	subroutine  gxtrithmfr(x,dx,y,dy,z,dz,t,dt,a,da,e,de,nu,dnu,h,dh,
     *                         alpha,df,globalflag)
c
c--------------------------------------------------------------------------
c
c This routine caclulates the derivatives of the membrane themally induced 
c mechanical force for the 18 dof three node triangle shell element. 
c Coded by Joe Pajot on 2/24/03
c 
c input variables (and their respective derivative terms preceeded by d):
c	x  = x coordinates of nodes
c       y  = y coordinates of nodes
c	z  = z coordinates of nodes
c       t  = mean temperature difference in element 
c	a  = coefficent of thermal expansion
c	e  = Young's modulus
c	nu = Poisson's ratio
c       h  = shell thickness
c       globalflag = flag to return to global coordinates
c
c output variables:
c	df  = derivative of thermally induced mechanical force
c
c-------------------------------------------------------------------------
C
	real *8 t,a,e,nu,h,alpha
	real *8 dt,da,de,dnu,dh
	real*8 x(3),y(3),z(3),f(18)
	real*8 dx(3),dy(3),dz(3),df(18)
	real*8 xp(3),yp(3),zp(3),xlp(3),ylp(3),zlp(3)
	real*8 dxp(3),dyp(3),dzp(3),dxlp(3),dylp(3),dzlp(3)
	real*8 str(3),temp(3),ftemp(9)
	real*8 dstr(3),dtemp(3),dftemp(9)
	real*8 dm(3,3),p(9,3),rot(3,3)
	real*8 ddm(3,3),dp(9,3),drot(3,3)
	real*8 v1n(3),v2n(3),v3n(3)
	real*8 x21,y21,z21,x32,y32,z32,x13,y13,z13
	real*8 dx21,dy21,dz21,dx32,dy32,dz32,dx13,dy13,dz13
	real*8 x12,y12,z12,x23,y23,z23,x31,y31,z31
	real*8 dx12,dy12,dz12,dx23,dy23,dz23,dx31,dy31,dz31
	real*8 cb,rlr,rlb,bpr,area,area2,coef1,coef2
	real*8 dcb,drlr,drlb,dbpr,darea,darea2
	real*8 dxlcg,dylcg,dzlcg,dylr,dzlr,dxcg,dycg,dzcg
	data v1n/1.0,0.0,0.0/
        data v2n/0.0,1.0,0.0/
        data v3n/0.0,0.0,1.0/
	integer i,j,globalflag
c
c  dimension variables
c
	x21 = x(2) - x(1)
	y21 = y(2) - y(1)
        z21 = z(2) - z(1)
	x32 = x(3) - x(2)
	y32 = y(3) - y(2)
        z32 = z(3) - z(2)
	x13 = x(1) - x(3)
	y13 = y(1) - y(3)
        z13 = z(1) - z(3)
c
	dx21 = dx(2) - dx(1)
	dy21 = dy(2) - dy(1)
        dz21 = dz(2) - dz(1)
	dx32 = dx(3) - dx(2)
	dy32 = dy(3) - dy(2)
        dz32 = dz(3) - dz(2)
	dx13 = dx(1) - dx(3)
	dy13 = dy(1) - dy(3)
        dz13 = dz(1) - dz(3)
c
c triangle in space : we compute the length of one side and the distance of the
c opposing node to that side to compute the area
c
        rlr = dsqrt( x21*x21 + y21*y21 + z21*z21 )
c       rlb = dsqrt( x32*x32 + y32*y32 + z32*z32 )
c       bpr = dsqrt((x21 * x32 + y21 * y32 + z21 *z32 )**2)/rlr
c       area= rlr*(dsqrt(rlb**2-bpr**2))/2.0d+00
c
        drlr = ( dx21*x21 + dy21*y21 + dz21*z21 )/rlr
c        drlb = ( dx32*x32 + dy32*y32 + dz32*z32 )/rlb
c        dbpr = (x21 * x32 + y21 * y32 + z21 *z32 )*(x21*dx32 + dx21*x32 +
c     *	        y21*dy32 + dy21*y32 + z21*dz32 +dz21*z32)/
c     *          (rlr*dsqrt((x21 * x32 + y21 * y32 + z21 *z32 )**2)) -
c     *          drlr*bpr/rlr
c	darea= drlr*(dsqrt(rlb**2-bpr**2))/2.0d+00 + 
c     *         rlr*(rlb*drlb - bpr*dbpr)/(2*dsqrt(rlb**2-bpr**2))
c
c direction cosines of the local system . X' is directed parallel 
c to the 2-1 side
c Z' is the external normal (counterclockwise). Y' computed as Z' x X'
c
        xp(1) = x21/rlr
        xp(2) = y21/rlr
        xp(3) = z21/rlr
c
	dxp(1) = dx21/rlr - drlr*x21/(rlr*rlr)
        dxp(2) = dy21/rlr - drlr*y21/(rlr*rlr)
        dxp(3) = dz21/rlr - drlr*z21/(rlr*rlr)
c
        zp(1)  = y21 * z32 - z21 * y32
        zp(2)  = z21 * x32 - x21 * z32
        zp(3)  = x21 * y32 - y21 * x32
	dzp(1) = dy21 * z32 - dz21 * y32 + y21 * dz32 - z21 * dy32
        dzp(2) = dz21 * x32 - dx21 * z32 + z21 * dx32 - x21 * dz32 
        dzp(3)  = dx21 * y32 - dy21 * x32 + x21 * dy32 - y21 * dx32
        zlr    = dsqrt( zp(1)*zp(1) + zp(2)*zp(2)+ zp(3)*zp(3) )
	dzlr   = ( dzp(1)*zp(1) + dzp(2)*zp(2)+ dzp(3)*zp(3) )/zlr
        dzp(1)  = dzp(1)/zlr - dzlr*zp(1)/(zlr*zlr)
        dzp(2)  = dzp(2)/zlr - dzlr*zp(2)/(zlr*zlr)
        dzp(3)  = dzp(3)/zlr - dzlr*zp(3)/(zlr*zlr)
	zp(1)  = zp(1)/zlr
        zp(2)  = zp(2)/zlr
        zp(3)  = zp(3)/zlr
c
	yp(1) = zp(2) * xp(3) - zp(3) * xp(2)
        yp(2) = zp(3) * xp(1) - zp(1) * xp(3)
        yp(3) = zp(1) * xp(2) - zp(2) * xp(1)
	dyp(1) = dzp(2) * xp(3) - dzp(3) * xp(2) + 
     *	          zp(2) * dxp(3) - zp(3) * dxp(2)
        dyp(2) = dzp(3) * xp(1) - dzp(1) * xp(3) + 
     *	          zp(3) * dxp(1) - zp(1) * dxp(3)
        dyp(3) = dzp(1) * xp(2) - dzp(2) * xp(1) + 
     *	          zp(1) * dxp(2) - zp(2) * dxp(1)
	ylr   = dsqrt( yp(1)*yp(1) + yp(2)*yp(2) + yp(3)*yp(3) )
	dylr   = ( dyp(1)*yp(1) + dyp(2)*yp(2) + dyp(3)*yp(3) )/ylr
	dyp(1)  = dyp(1)/ylr - dylr*yp(1)/(ylr*ylr)
        dyp(2)  = dyp(2)/ylr - dylr*yp(2)/(ylr*ylr)
        dyp(3)  = dyp(3)/ylr - dylr*yp(3)/(ylr*ylr)
	yp(1)  = yp(1)/ylr
        yp(2)  = yp(2)/ylr
        yp(3)  = yp(3)/ylr
c
c compute center of gravity
c
        xcg = (x(1) + x(2) + x(3))/3.0d+00
        ycg = (y(1) + y(2) + y(3))/3.0d+00
        zcg = (z(1) + z(2) + z(3))/3.0d+00
	dxcg = (dx(1) + dx(2) + dx(3))/3.0d+00
        dycg = (dy(1) + dy(2) + dy(3))/3.0d+00
        dzcg = (dz(1) + dz(2) + dz(3))/3.0d+00
c
c compute local coordinates 
        do  i=1,3
          xlcg   = x(i) - xcg
          ylcg   = y(i) - ycg
          zlcg   = z(i) - zcg
	  dxlcg   = dx(i) - dxcg
          dylcg   = dy(i) - dycg
          dzlcg   = dz(i) - dzcg
          xlp(i) = xp(1) * xlcg + xp(2) * ylcg + xp(3) * zlcg
          ylp(i) = yp(1) * xlcg + yp(2) * ylcg + yp(3) * zlcg
          zlp(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg
	  dxlp(i) = dxp(1) * xlcg + dxp(2) * ylcg + dxp(3) * zlcg +
     *               xp(1) * dxlcg + xp(2) * dylcg + xp(3) * dzlcg
          dylp(i) = dyp(1) * xlcg + dyp(2) * ylcg + dyp(3) * zlcg +
     *   	     yp(1) * dxlcg + yp(2) * dylcg + yp(3) * dzlcg
          dzlp(i) = dzp(1) * xlcg + dzp(2) * ylcg + dzp(3) * zlcg +
     *               zp(1) * dxlcg + zp(2) * dylcg + zp(3) * dzlcg
  	end do	
c
c  dimension variables in local coordinates
c		
      x21 =      xlp(2) - xlp(1)
      x12 =     -x21
      x32 =      xlp(3) - xlp(2)
      x23 =     -x32
      x13 =      xlp(1) - xlp(3)
      x31 =     -x13
      y21 =      ylp(2) - ylp(1)
      y12 =     -y21
      y32 =      ylp(3) - ylp(2)
      y23 =     -y32
      y13 =      ylp(1) - ylp(3)
      y31 =     -y13
      area2 =    y21*x13 - x21*y13
c
      dx21 =      dxlp(2) - dxlp(1)
      dx12 =     -dx21
      dx32 =      dxlp(3) - dxlp(2)
      dx23 =     -dx32
      dx13 =      dxlp(1) - dxlp(3)
      dx31 =     -dx13
      dy21 =      dylp(2) - dylp(1)
      dy12 =     -dy21
      dy32 =      dylp(3) - dylp(2)
      dy23 =     -dy32
      dy13 =      dylp(1) - dylp(3)
      dy31 =     -dy13
      darea2 =    dy21*x13 - dx21*y13 + y21*dx13 - x21*dy13
c
c membrane elastic matrix
c
      cb=e*(h/2)/(1.0-(nu*nu))
      dm(1,1) =    cb
      dm(1,2) = nu*cb
      dm(1,3) = 0.0
      dm(2,1) = dm(1,2)
      dm(2,2) = cb
      dm(2,3) = 0.0
      dm(3,1) = 0.0
      dm(3,2) = 0.0
      dm(3,3) = ((1.0-nu)/2.0)*cb	
c
      dcb=(de*(h/2)+e*(dh/2))/(1.0-(nu*nu)) - 
     *      e*(h/2)*(-2*nu*dnu)/(1- 2*nu*nu+ nu**4)
      ddm(1,1) =    dcb
      ddm(1,2) = dnu*cb + nu*dcb
      ddm(1,3) = 0.0
      ddm(2,1) = ddm(1,2)
      ddm(2,2) = dcb
      ddm(2,3) = 0.0
      ddm(3,1) = 0.0
      ddm(3,2) = 0.0
      ddm(3,3) = ((1.0-nu)/2.0)*dcb - cb*dnu/2	
c
c  create strain vector
c
      str(1) = a*t 
      str(2) = a*t 
      str(3) = 0.0d+00	
      dstr(1) = da*t + a*dt 
      dstr(2) = da*t + a*dt 
      dstr(3) = 0.0d+00	
c
c  create strain-displacement matrix p 
c
      p(1,1) =   y23
      p(2,1) =   0.0
      p(3,1) =   y31
      p(4,1) =   0.0
      p(5,1) =   y12
      p(6,1) =   0.0
      p(1,2) =   0.0
      p(2,2) =   x32
      p(3,2) =   0.0
      p(4,2) =   x13
      p(5,2) =   0.0
      p(6,2) =   x21
      p(1,3) =   x32
      p(2,3) =   y23
      p(3,3) =   x13
      p(4,3) =   y31
      p(5,3) =   x21
      p(6,3) =   y12
      coef1  = alpha/6.0
      coef2  = alpha/3.0
      p(7,1) =  y23*(y13-y21)*coef1
      p(7,2) =  x32*(x31-x12)*coef1
      p(7,3) =  (x31*y13-x12*y21)*coef2
      p(8,1) =  y31*(y21-y32)*coef1
      p(8,2) =  x13*(x12-x23)*coef1
      p(8,3) =  (x12*y21-x23*y32)*coef2
      p(9,1) =  y12*(y32-y13)*coef1
      p(9,2) =  x21*(x23-x31)*coef1
      p(9,3) =  (x23*y32-x31*y13)*coef2
c
      dp(1,1) =   dy23
      dp(2,1) =   0.0
      dp(3,1) =   dy31
      dp(4,1) =   0.0
      dp(5,1) =   dy12
      dp(6,1) =   0.0
      dp(1,2) =   0.0
      dp(2,2) =   dx32
      dp(3,2) =   0.0
      dp(4,2) =   dx13
      dp(5,2) =   0.0
      dp(6,2) =   dx21
      dp(1,3) =   dx32
      dp(2,3) =   dy23
      dp(3,3) =   dx13
      dp(4,3) =   dy31
      dp(5,3) =   dx21
      dp(6,3) =   dy12
      dp(7,1) =  coef1*(dy23*(y13-y21) + y23*(dy13-dy21))
      dp(7,2) =  coef1*(dx32*(x31-x12) + x32*(dx31-dx12))
      dp(7,3) =  coef2*(dx31*y13-dx12*y21 + x31*dy13-x12*dy21)
      dp(8,1) =  coef1*(dy31*(y21-y32) + y31*(dy21-dy32))
      dp(8,2) =  coef1*(dx13*(x12-x23) + x13*(dx12-dx23))
      dp(8,3) =  coef2*(dx12*y21-dx23*y32 + x12*dy21-x23*dy32)
      dp(9,1) =  coef1*(dy12*(y32-y13) + y12*(dy32-dy13))
      dp(9,2) =  coef1*(dx21*(x23-x31) + x21*(dx23-dx31))
      dp(9,3) =  coef2*(dx23*y32-dx31*y13 + x23*dy32-x31*dy13)
c
c create vector = p'*D*str
c            
	do i=1,3
	   temp(i) = dm(i,1)*str(1) + dm(i,2)*str(2) + dm(i,3)*str(3)
	  dtemp(i) = ddm(i,1)*str(1) + ddm(i,2)*str(2) + ddm(i,3)*str(3) +
     *               dm(i,1)*dstr(1) + dm(i,2)*dstr(2) + dm(i,3)*dstr(3)
	end do
	do j=1,9
	   ftemp(j) = p(j,1)*temp(1) + p(j,2)*temp(2) + p(j,3)*temp(3)
	  dftemp(j) = dp(j,1)*temp(1) + dp(j,2)*temp(2) + dp(j,3)*temp(3) +
     *               p(j,1)*dtemp(1) + p(j,2)*dtemp(2) + p(j,3)*dtemp(3)
	end do
c
      if (globalflag .eq. 0) then
c transform to global coordinates	
c	first create the transformation matricies  rot and drot,
c	next carry out this multipication node by node which has 
c       forces [flx fly] 
        call g_rotation(xp,dxp,yp,dyp,zp,dzp, v1n,v2n,v3n,rot,drot)
	do i =1,3
	  df(i*6-5) =  drot(1,1)*ftemp(i*2-1) + drot(1,2)*ftemp(i*2) +
     *                 rot(1,1)*dftemp(i*2-1) + rot(1,2)*dftemp(i*2)
	  df(i*6-4) =  drot(2,1)*ftemp(i*2-1) + drot(2,2)*ftemp(i*2) +
     *                 rot(2,1)*dftemp(i*2-1) + rot(2,2)*dftemp(i*2)
	  df(i*6-3) =  drot(3,1)*ftemp(i*2-1) + drot(3,2)*ftemp(i*2) +
     *                 rot(3,1)*dftemp(i*2-1) + rot(3,2)*dftemp(i*2) 
	  df(i*6-2) =  drot(1,3)*ftemp(6+i) + rot(1,3)*dftemp(6+i)
	  df(i*6-1) =  drot(2,3)*ftemp(6+i) + rot(2,3)*dftemp(6+i)
	  df(i*6  ) =  drot(3,3)*ftemp(6+i) + rot(3,3)*dftemp(6+i)
	end do
      else
        do i=1,3
	   df(i*6-5) = dftemp(i*2-1)
	   df(i*6-4) = dftemp(i*2  )
	   df(i*6-3) = 0
	   df(i*6-2) = 0
	   df(i*6-1) = 0
	   df(i*6  ) = dftemp(6+i)
	  end do
      end if
      end

	
	
