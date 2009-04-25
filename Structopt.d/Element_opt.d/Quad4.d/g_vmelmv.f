	subroutine gxvmelmv(stress,dstress,maxgus,maxstr,msize,
     *                      elm,nno)
*
C.... DECLARE ALL GLOBAL VARIABLES
C
C.... INTEGER CONSTANTS
C
	integer maxgus,maxstr,msize,elm,nno
C
C.... REAL ARRAYS
C
	real*8  stress(msize,maxstr,maxgus)
	real*8 dstress(msize,maxstr,maxgus)

C
C.... DECLARE ALL LOCAL VARIABLES
C
	integer n
	real*8  sxx,syy,szz,sxy,sxz,syz
        real*8  dxsxx,dxsyy,dxszz,dxsxy,dxsxz,dxsyz
	real*8  dsxx,dsyy,dszz,dsxy,dsxz,dsyz
        real*8  dxdsxx,dxdsyy,dxdszz,dxdsxy,dxdsxz,dxdsyz
	real*8  j2,comp
	real*8  dj2,dcomp,dvms
C
C.... INITIALIZE CENTROIDAL STRESS TO ZERO
C
       dvms = 0.0d0
C
C.... COMPUTE THE CENTROIDAL STRESS FROM THE NODAL STRESSES
C
	do 10 n = 1, nno
	  sxx = stress(elm,1,n)
          syy = stress(elm,2,n)
          szz = stress(elm,3,n)
          sxy = stress(elm,4,n)
          syz = stress(elm,5,n)
          sxz = stress(elm,6,n)
	  
	  dxsxx = dstress(elm,1,n)
          dxsyy = dstress(elm,2,n)
          dxszz = dstress(elm,3,n)
          dxsxy = dstress(elm,4,n)
          dxsyz = dstress(elm,5,n)
          dxsxz = dstress(elm,6,n)
C
C.... COMPUTE THE FIRST DEVEATORIC STRESSES
C
 	  comp = (sxx + syy + szz)/3.0d0
	  dsxx = sxx - comp
          dsyy = syy - comp
          dszz = szz - comp
          dsxy = sxy
          dsyz = syz
          dsxz = sxz
	
	  dcomp = (dxsxx + dxsyy + dxszz)/3.0d0
	  dxdsxx = dxsxx - dcomp
          dxdsyy = dxsyy - dcomp
          dxdszz = dxszz - dcomp
          dxdsxy = dxsxy
          dxdsyz = dxsyz
          dxdsxz = dxsxz
C
C.... COMPUTE THE SECOND DEVEATORIC STRESS
C
	  j2 = ((dsxx*dsxx)+(dsyy*dsyy)+(dszz*dszz))/2.0d0+
     &          (dsxy*dsxy)+(dsyz*dsyz)+(dsxz*dsxz)
          j2 = dsqrt(3.0d0*j2)
     
	  dj2 = ((dsxx*dxdsxx)+(dsyy*dxdsyy)+(dszz*dxdszz)) +
     &          ((dsxy*dxdsxy)+(dsyz*dxdsyz)+(dsxz*dxdsxz))*2.0d0
C
C.... COMPUTE THE VON MISES STRESS
C
          if (j2 .le. 0.0) then
            write(6,*)'Negative / zero equivalent stress in g_vmelmv.f'
          else
   	    dvms = dvms + 1.5d0*dj2/j2
          endif
C
   10   continue
C
        dvms = dvms/nno
C
C.... DISTRIBUTE OUT TO THE NODES 
C
	do 20 n = 1, nno
          dstress(elm,7,n) = dvms
20	continue
C
	return
	end

