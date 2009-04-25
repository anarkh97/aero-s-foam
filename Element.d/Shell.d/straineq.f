        subroutine straineq(rmx,rmy,rmxy,rnx,rny,rnxy,t,surface,ebar)
***************************************************************
*       THIS ROUTINE CALCULATES THE EQUIVALENT VON MISES      *
*       STRAIN FOR THE 3 NODE SHELL ELEMENT                   *
***************************************************************
*                                                             *
*       AUTHOR  :       K.H. PIERSON                          *
*       DATE    :       MARCH 1997                            *
*       VERSION :       FEM-C++ 1.00                          *
*                                                             *
***************************************************************
*                                                             *
*	rmx  = centroidal bending curvature (kxxc)            *
*	rmy  = centroidal bending curvature (kyyc)            *
*	rmxy = centroidal bending curvature (kxyc)            *
*	rnx  = centroidal membrane strain   (exxc)            *
*	rny  = centroidal membrane strain   (eyyc)            *
*	rnxy = centroidal membrane strain   (exyc)            *
*	t    = element thickness                              *
*	ebar = equivalent strain                              *
*                                                             *
***************************************************************
*                                                             *
*       CALLED BY :  SANDS8.F                                 *
*                                                             *
***************************************************************

C ... ARGUMENTS
        integer surface
        double precision rmx,rmy,rmxy,rnx,rny,rnxy,t,ebar

C ... LOCAL VARIABLES
        double precision ex,ey,exy,etop,ebot,emid
c       double precision th

C ... RETURN IF MEDIAN SURFACE IS REQUESTED
        if(surface .eq. 2) then
c	  call equiv(rmx,rmy,rmxy,emid)
          call equiv(rnx,rny,0.5*rnxy,emid)
          ebar = emid
          return
        endif

C ... DIVIDE THICKNESS BY 2
c 	th = 0.5*t

C ... CALCULATE STRAINS AT TOP SURFACE
        ex  =  rnx +  rmx
        ey  =  rny +  rmy
C Strain Tensor Versus Engineering Strain
        exy = 0.5*(rnxy + rmxy)

C ... COMPUTE EQUIVALENT STRAIN AT TOP SURFACE
        call equiv(ex,ey,exy,etop)

C ... RETURN IF TOP SURFACE VALUE IS REQUESTED
        if(surface .eq. 1) then
          ebar = etop
          return
        endif

C ... CALCULATE STRAINS AT BOTTOM SURFACE
        ex   =  rnx -  rmx
        ey   =  rny -  rmy
        exy  = 0.5*(rnxy - rmxy)

C ... COMPUTE EQUIVALENT STRAIN AT BOTTOM SURFACE
        call equiv(ex,ey,exy,ebot)

C ... RETURN IF BOTTOM SURFACE VALUE IS REQUESTED
        if(surface .eq. 3) then
          ebar = ebot
          return
        endif


C ... RETURN THE MAXIMUM EQUIVALENT STRAIN
        ebar = max(etop,ebot)

        return
        end

C
C ... SUBROUTINE TO CALCULATE EQUIVALENT STRAIN
C
        subroutine equiv(ex,ey,exy,eq)

C ... ARGUMENTS
        double precision ex,ey,exy,eq

C ... LOCAL VARIABLES
        double precision e0,dex,dey,dez
        
C ... COMPUTE MEAN HYDROSTATIC STRAIN
        e0 = (ex + ey)/3.0d0

C ... COMPUTE DEVIATORIC STRAINS
        dex = ex - e0
        dey = ey - e0
        dez =    - e0

C ... COMPUTE EQUIVALENT STRAIN
C THIS COMPUTATION OF VON-MISES STRAIN MAKES NO SENSE
C       eq = dsqrt((2.0d0/3.0d0)*dex*dex + dey*dey + dez*dez)
C COMPUTE LIKE VM STRESS
        eq = ((dex*dex + dey*dey + dez*dez)/2.0d0) + (exy*exy)
        eq = dsqrt(3.0d0 * eq)

        return
        end
