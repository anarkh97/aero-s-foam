***************************************************************
* THIS SUBROUTINE COMPUTE THE ELEMENT MASS MATRIX FOR THE 20- *
* NODE BRICK IN A LUMPPED FORM.                               *
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
	subroutine gxbr20mas(p,rho,drho,delmass,x,y,z,dx,dy,dz,
     *                       gamma,grvfor,grvflg,dtotmas)
C
C..... DECLARE GLOBAL VARIABLES
C
C.... INTEGER CONSTANTS
C
C       implicit none
	integer p
C
C.... REAL CONSTANTS
C
	real*8 rho, drho, dtotmas, coef, cvar
C
C.... REAL ARRAYS
C
	real*8 delmass(60,60),z(*),x(*),y(*),dz(*),dx(*),dy(*)
        real*8 gamma(*),grvfor(*)
C
C.... LOGICAL CONSTANT
C
        logical grvflg
C
C.... DECLARE LOCAL VARIABLES FOR Q4DMAS
C
	integer k,l,jj, i
	real*8 v,dv,dc1,xi,eta,mu,weight,det
	real*8 q(20),qx(20),qy(20),qz(20)
	real*8 dqx(20),dqy(20),dqz(20),ddet
C
C.... DETERMINE THE VOLUME OF THE BRICK 
C
	v  = 0.0d0
	dv = 0.0d0

C     check if any variation
      cvar = 0.0
      
      do 5000 i = 1, 8
        cvar = cvar + dx(i)**2 + dy(i)**2 +dz(i)**2 
 5000 continue

      if (cvar.eq.0) then
C
	do 10 k=1, p
          do 20 l=1, p
            do 30 jj = 1, p
              call hxgaus20(p,k,p,l,p,jj,xi,eta,mu,weight)
              call h20shpe(xi,eta,mu,x,y,z,q,qx,qy,qz,det)
	      v  =  v + det
30	    continue
20	  continue
10	continue

C
C ... KHP
C ... Modified 5-6-98
C
	if(v .le. 0.0d0) then
          write(*,*) 'Negative/zero volume'
        endif
C
C.... COMPUTE THE TOTAL ELEMENT MASS
C
         dtotmas = drho*v 

      else
C
	do 110 k=1, p
          do 120 l=1, p
            do 130 jj = 1, p
              call hxgaus20(p,k,p,l,p,jj,xi,eta,mu,weight)
              call gxh20shpe(xi,eta,mu,x,y,z,dx,dy,dz,q,qx,
     *                       qy,qz,dqx,dqy,dqz,det,ddet)
	      v  =  v + det
	      dv = dv + ddet
130	    continue
120	  continue
110	continue

C
C ... KHP
C ... Modified 5-6-98
C
	if(v .le. 0.0d0) then
          write(*,*) 'Negative/zero volume'
        endif
C
C.... COMPUTE THE TOTAL ELEMENT MASS
C
         dtotmas = drho*v + rho*dv

      endif
C
C.... COMPUTE THE MATRIX COEFFICIENT
C
	dc1 = 0.125d0*dtotmas
C
C.... COMPUTE BRICK ELEMENT MASS MATRIX
C
	do 40 k = 1, 60
          do 50 l = 1, 60
            delmass(k,l) = 0.0d0
50	  continue
40	continue

	do 60 k=1, 60
	  delmass(k,k) = dc1
60	continue
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
C        if (grvflg) then
C	  coef = 8.0d0*c1
C          grvfor(1) = coef*gamma(1)
C          grvfor(2) = coef*gamma(2)
C          grvfor(3) = coef*gamma(3)
C        endif
C
	return
	end
