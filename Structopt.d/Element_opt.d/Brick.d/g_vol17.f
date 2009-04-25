      subroutine gxvol17(p,x,y,z,dx,dy,dz,volume,dvolume)
C
C      implicit none
      integer p
      real*8  volume,dvolume
      real*8  z(*),x(*),y(*)
      real*8  dz(*),dx(*),dy(*)

      integer k,l,jj,i
      real*8 v,dv,xi,eta,emu,weight,det,ddet
      real*8 q(8),qx(8),qy(8),qz(8)
      real*8 dqx(8),dqy(8),dqz(8), cvar
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
	do 110 k=1, p
          do 120 l=1, p
            do 130 jj = 1, p
              call hxgaus(p,k,p,l,p,jj,xi,eta,emu,weight)
              call h8shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
	      v = v + det
130	    continue
120	  continue
110	continue

	volume  =  v
 
      else
      
	do 10 k=1, p
          do 20 l=1, p
            do 30 jj = 1, p
              call hxgaus(p,k,p,l,p,jj,xi,eta,emu,weight)
              call gxh8shpe(xi,eta,emu,x,y,z,dx,dy,dz,q,qx,qy,qz,
     *                    dqx,dqy,dqz,det,ddet)
	      v = v + det
	      dv = dv + ddet
30	    continue
20	  continue
10	continue
        
	volume  =  v
        dvolume = dv
	
      endif

      return
      end
