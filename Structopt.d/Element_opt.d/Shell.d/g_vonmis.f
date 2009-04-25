      subroutine  g_vonmis(rmx,  grmx,  rmy,  grmy,
     *                     rmxy, grmxy, rnx,  grnx, 
     *                     rny,  grny,  rnxy, grnxy, 
     *                     t,    gt,    sv,   gsv,   surface)

C.... GLOBAL VARIABLES
      integer surface
      double precision rmx,rmy,rmxy,rnx,rny,rnxy,t,sv
      double precision grmx,grmy,grmxy,grnx,grny,grnxy,gt,gsv

C.... LOCAL VARIABLES
C st = von mises stress in top surface
C sm = von mises stress in median surface
C sb = von mises stress in bottom surface

      double precision sx,sy,sxy,st,sb,sm,t2,sq3
      double precision gsx,gsy,gsxy,gst,gsb,gsm,gt2
      double precision rnxt,rnyt,rnxyt,rmxt,rmyt,rmxyt
      double precision grnxt,grnyt,grnxyt,grmxt,grmyt,grmxyt
      double precision xt,gxt
c
      sq3  = dsqrt(3.0d0)
c
      xt   = dabs(t)
      t2   = xt*xt
c
      if (t .gt. 0.0d0) then
        gxt =  gt
      else if (t .lt. 0.0d0) then
        gxt = -gt
      else 
        gxt = 0.0d0
      endif
c
      t    = xt
c
      gt2  = 2.0d0 * xt * gxt
c
      rnxt  =  rnx/xt
      rnyt  =  rny/xt
      rnxyt = rnxy/xt
c
      grnxt  = (grnx  -  gxt*rnx/xt )/xt
      grnyt  = (grny  -  gxt*rny/xt )/xt
      grnxyt = (grnxy -  gxt*rnxy/xt)/xt
c
      rmxt  = 6.0d0* rmx/t2
      rmyt  = 6.0d0* rmy/t2
      rmxyt = 6.0d0*rmxy/t2
c
      grmxt  = 6.0d0* ( grmx - gt2*rmx/t2)  /t2 
      grmyt  = 6.0d0* ( grmy - gt2*rmy/t2)  /t2 
      grmxyt = 6.0d0* (grmxy - gt2*rmxy/t2) /t2

C ... COMPUTE VON MISES STRESS IN BOTTOM SURFACE

      sx  =  rnxt -  rmxt
      sy  =  rnyt -  rmyt
      sxy = rnxyt - rmxyt

      gsx  =  grnxt -  grmxt
      gsy  =  grnyt -  grmyt
      gsxy = grnxyt - grmxyt

      call g_compj2(sx,gsx,sy,gsy,sxy,gsxy,sb,gsb)

      sb  = sq3 * sb
      gsb = sq3 * gsb

      if(surface .eq. 3) then 
        sv  = sb
        gsv = gsb
	return
      endif

C ... COMPUTE VON MISES STRESS IN MEDIAN SURFACE

      if(surface .eq. 2) then

        sx  = rnxt
        sy  = rnyt
        sxy = rnxyt

        gsx  = grnxt
        gsy  = grnyt
        gsxy = grnxyt

        call g_compj2(sx,gsx,sy,gsy,sxy,gsxy,sm,gsm)

	sv  = sq3 * sm
	gsv = sq3 * gsm

	return

      endif

C ... COMPUTE VON MISES STRESS IN TOP SURFACE

      sx  =  rnxt +  rmxt
      sy  =  rnyt +  rmyt
      sxy = rnxyt + rmxyt

      gsx  =  grnxt +  grmxt
      gsy  =  grnyt +  grmyt
      gsxy = grnxyt + grmxyt

      call g_compj2(sx,gsx,sy,gsy,sxy,gsxy,st,gst)

      st  = sq3 *  st
      gst = sq3 * gst

      if(surface .eq. 1) then 
        sv =  st
        gsv = gst
	return
      endif

C ... VON MISES STRESS = max(sb,st)

      if (sb.gt.st) then
        sv   = sb
        gsv  = gsb
      elseif (sb.lt.st) then
        sv   = st
        gsv  = gst
      else
        sv   = sb
        gsv  = 0.5d0 * (gsb+gst)
        write(*,*) 'gsv set to average'
      endif

      return
      end
C
C ... SUBROUTINE TO CALCULATE J2
C
      subroutine g_compj2(sx,gsx,sy,gsy,sxy,gsxy,svm,gsvm)

C ... GLOBAL VARIABLES
      double precision sx,sy,sxy,svm
      double precision gsx,gsy,gsxy,gsvm

C ... LOCAL VARIABLES
      double precision sz,s0,dsx,dsy,dsz,j2
      double precision gsz,gs0,gdsx,gdsy,gdsz,gj2

C ... SET sz = 0 TO REMIND USER OF THIS ASSUMPTION
      sz  = 0.0d0
      gsz = 0.0d0

C ... COMPUTE AVERAGE HYDROSTATIC STRESS
      s0  = (sx + sy + sz)/3.0d0
      gs0 = (gsx + gsy + gsz)/3.0d0

C ... COMPUTE DEVIATORIC STRESSES
      dsx  = sx - s0
      dsy  = sy - s0
      dsz  = sz - s0 

      gdsx  = gsx - gs0
      gdsy  = gsy - gs0
      gdsz  = gsz - gs0 
C
C ... COMPUTE J2

      j2 = 0.5d0*((dsx*dsx) + (dsy*dsy) + (dsz*dsz)) + (sxy*sxy)

      gj2 = ( gdsx*dsx + gdsy*dsy + gdsz*dsz ) + 2.0d0*gsxy*sxy

C ... COMPUTE VON MISES STRESS

      svm  = dsqrt(j2)

      if (j2.ne.0.0d0) then
        gsvm = gj2/2.0d0/svm
      else
        gsvm = 0.0d0
      endif

      return
      end
