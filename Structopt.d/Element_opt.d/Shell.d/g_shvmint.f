       subroutine gxshvmint(xl,yl,zl,dxl,dyl,dzl,e,de,nu,h,dh,
     &                    v,dv,vmint,dvmint,vol,dvol,sigbar,fac,iarea,
     &                    alpha,beta)
*-------------------------------------------------------------------*
*  CALLED BY : ThreeNodeShell.C 
*
*-------------------------------------------------------------------* 
C
C t = thickness of triangle
c
c       implicit none
       double precision  e,de,nu,alpha,beta
       double precision  xl(*),yl(*),zl(*),h(*),v(*)
       double precision  dxl(*),dyl(*),dzl(*),dh(*),dv(*)
       double precision  dm(3,3),ddm(3,3),db(3,3),ddb(3,3)
       double precision  xp(3),yp(3),zp(3),xlp(3),ylp(3),zlp(3)
       double precision  yptmp(3),zptmp(3)
       double precision  r1(3,3),dll(18)
       double precision  dxp(3),dyp(3),dzp(3),dxlp(3),dylp(3),dzlp(3)
       double precision  dr1(3,3),ddll(18)
       double precision  xg(3),yg(3),zg(3)
       double precision  dxg(3),dyg(3),dzg(3)
       double precision  cb,dcb,x21,y21,z21,x32,y32,z32,x13,y13,z13
       double precision  dx21,dy21,dz21,dx32,dy32,dz32,dx13,dy13,dz13
       double precision  rx,ry,rz,bx,by,bz,rlr,rlb,bpr,area,ylr,zlr
       double precision  drx,dry,drz,dbxd,dby,dbz,drlr,drlb,dbpr
       double precision  darea,dylr,dzlr,sqrl,dsqrl
       double precision  ycg,xcg,zcg,xlcg,ylcg,zlcg,t,t2,t3,dt
       double precision  dycg,dxcg,dzcg,dxlcg,dylcg,dzlcg
       double precision  rmem(18,3),drmem(18,3),a,b
       double precision  rmom(18,3),drmom(18,3)
       double precision  rmx,rmy,rmxy,rnx,rny,rnxy
       double precision  rmmx,rmmy,rmmxy,rnnx,rnny,rnnxy,sbf
       double precision  drmx,drmy,drmxy,drnx,drny,drnxy
       double precision  drmmx,drmmy,drmmxy,drnnx,drnny,drnnxy
       double precision  vol, sigbar, vmint, dvmint, fac, dvol, dfac
       double precision  tsx,tsy,tsxy,msx,msy,msxy
       double precision  bsx,bsy,bsxy,dbsx,dbsy,dbsxy
       double precision  dtsx,dtsy,dtsxy,dmsx,dmsy,dmsxy
       double precision  tdsx,tdsy,tdsz,dtdsx,dtdsy,dtdsz
       double precision  mdsx,mdsy,mdsz,dmdsx,dmdsy,dmdsz
       double precision  bdsx,bdsy,bdsz,dbdsx,dbdsy,dbdsz
       double precision  ts0,ms0,bs0,dts0,dms0,dbs0
       double precision  tj2,mj2,bj2,tj2tmp,mj2tmp,bj2tmp,dintj2,intj2
       double precision  dtj2,dmj2,dbj2,dtj2tmp,dmj2tmp,dbj2tmp,cvar
       double precision  rnxt,rnyt,rnxyt,drnxt,drnyt,drnxyt
       double precision  rmxt,rmyt,rmxyt,drmxt,drmyt,drmxyt,dbx,dsq3
       character*10      status
       integer lb(9),le(9)
       integer i,j
       data lb/3,4,5,9,10,11,15,16,17/
       data le/1,2,7,8,13,14,6,12,18/
c
c	 PARAMETER(alpha = 1.50d+00)
c	 PARAMETER(beta  = 0.32d+00)
c       
c	PARAMETER(alpha = 0.0d+00)
c	PARAMETER(beta  = 0.0d+00)
c
C
       pow(a,b) = exp(b*log(a))
C
       vol    = 0.0d0
       dvol   = 0.0d0
       vmint  = 0.0d0
       dvmint = 0.0d0
C
       do 5 i = 1, 18
         rmom(i,1)  = 0.0d0
         rmom(i,2)  = 0.0d0
         rmom(i,3)  = 0.0d0
         drmom(i,1) = 0.0d0
         drmom(i,2) = 0.0d0
         drmom(i,3) = 0.0d0
         rmem(i,1)  = 0.0d0
         rmem(i,2)  = 0.0d0
         rmem(i,3)  = 0.0d0
         drmem(i,1) = 0.0d0
         drmem(i,2) = 0.0d0
         drmem(i,3) = 0.0d0
5      continue
C
	t    = h(1)
	dt   = dh(1)
C
C set the bending constitutive matrix
C
   	cb=e*(t*t*t)/12.0d0/(1.0d0-nu*nu)
	db(1,1) = cb
	db(1,2) = nu*cb
	db(1,3) = 0.0d0
	db(2,1) = db(1,2)
	db(2,2) = cb
	db(2,3) = 0.0d0
	db(3,1) = 0.0d0
	db(3,2) = 0.0d0
	db(3,3) = 0.5d0*(1.0d0-nu)*cb

   	dcb= de*t*t*t/12.0d0/(1.0d0-nu*nu)+ 
     *       e*3.0d0*t*t*dt/12.0d0/(1.0d0-nu*nu)
	ddb(1,1) = dcb
	ddb(1,2) = nu*dcb
	ddb(1,3) = 0.0d0
	ddb(2,1) = ddb(1,2)
	ddb(2,2) = dcb
	ddb(2,3) = 0.0d0
	ddb(3,1) = 0.0d0
	ddb(3,2) = 0.0d0
	ddb(3,3) = 0.5d0*(1.0d0-nu)*dcb
C
C set the membrane constitutive matrix
C
   	cb = e*t/(1.0d0-nu*nu)
	dm(1,1) = cb
	dm(1,2) = nu*cb
	dm(1,3) = 0.0D0
	dm(2,1) = nu*e*(t*t*t)/12.0d0/(1.0d0-nu*nu)
	dm(2,2) = cb
	dm(2,3) = 0.0d0
	dm(3,1) = 0.0d0
	dm(3,2) = 0.0d0
	dm(3,3) = 0.5d0*(1.0d0-nu)*cb

   	dcb = de*t/(1.0d0-nu*nu) + e*dt/(1.0d0-nu*nu)
	ddm(1,1) = dcb
	ddm(1,2) = nu*dcb
	ddm(1,3) = 0.0D0
	ddm(2,1) = nu*de*(t*t*t)/12.0d0/(1.0d0-nu*nu)
     *           + nu*e*3*(t*t)*dt/12.0d0/(1.0d0-nu*nu)
	ddm(2,2) = dcb
	ddm(2,3) = 0.0d0
	ddm(3,1) = 0.0d0
	ddm(3,2) = 0.0d0
	ddm(3,3) = 0.5d0*(1.0d0-nu)*dcb

C     check if any variation
      cvar = 0.0d0
      
      do 5000 i = 1,3
        cvar = cvar + dxl(i)**2 + dyl(i)**2 +dzl(i)**2 
 5000 continue

c                                                          cvar.eq.0       
      if (cvar.eq.0) then
C
C triangular dimension variables
C
	x21 = xl(2) - xl(1)
	y21 = yl(2) - yl(1)
        z21 = zl(2) - zl(1)
	x32 = xl(3) - xl(2)
	y32 = yl(3) - yl(2)
        z32 = zl(3) - zl(2)
	x13 = xl(1) - xl(3)
	y13 = yl(1) - yl(3)
        z13 = zl(1) - zl(3)

C  triangle in space : we compute the length of one side and the distance of the
C  opposing node to that side to compute the area
       rx   = x21
       ry   = y21
       rz   = z21
       bx   = x32
       by   = y32
       bz   = z32
       rlr  = dsqrt( rx*rx + ry*ry + rz*rz )
       rlb  = dsqrt( bx*bx + by*by + bz*bz )
       bpr  = abs(rx * bx + ry * by + rz *bz )/rlr
       area = 0.5d0*rlr*dsqrt(rlb*rlb-bpr*bpr)
       vol  = area*t

       darea = 0.0d0
       dvol  = area*dt
C
C Direction cosines of the local system . X' is directed parallel to the 
C 2-1 side Z' is the external normal (counterclockwise). 
C Y' computed as Z' x X'
C
       xp(1) = x21/rlr
       xp(2) = y21/rlr
       xp(3) = z21/rlr 
C
C Z' local axis
C
       zp(1) = y21 * z32 - z21 * y32
       zp(2) = z21 * x32 - x21 * z32
       zp(3) = x21 * y32 - y21 * x32
       zlr   = dsqrt( zp(1)*zp(1) + zp(2)*zp(2) + zp(3)*zp(3) )
       zp(1) = zp(1)/zlr
       zp(2) = zp(2)/zlr
       zp(3) = zp(3)/zlr
C		
C Y' local axis
C
       yp(1) = zp(2) * xp(3) - zp(3) * xp(2)
       yp(2) = zp(3) * xp(1) - zp(1) * xp(3)
       yp(3) = zp(1) * xp(2) - zp(2) * xp(1)
       ylr   = dsqrt( yp(1)*yp(1) + yp(2)*yp(2) + yp(3)*yp(3) )
       yp(1) = yp(1)/ylr
       yp(2) = yp(2)/ylr
       yp(3) = yp(3)/ylr

C
C center of gravity
C
       xcg = (xl(1) + xl(2) + xl(3))/3.0d+00
       ycg = (yl(1) + yl(2) + yl(3))/3.0d+00
       zcg = (zl(1) + zl(2) + zl(3))/3.0d+00
C
C computing local coordinates
C
       do 43 i=1,3
         xlcg   = xl(i) - xcg
         ylcg   = yl(i) - ycg
         zlcg   = zl(i) - zcg
         xlp(i) = xp(1) * xlcg + xp(2) * ylcg + xp(3) * zlcg
         ylp(i) = yp(1) * xlcg + yp(2) * ylcg + yp(3) * zlcg
         zlp(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg

 43    continue
C
C Set Global axes
C
       xg(1) = 1.0d0
       xg(2) = 0.0d0
       xg(3) = 0.0d0
       yg(1) = 0.0d0
       yg(2) = 1.0d0
       yg(3) = 0.0d0
       zg(1) = 0.0d0
       zg(2) = 0.0d0
       zg(3) = 1.0d0
C
C computing nodal rotation matrix
C
       call rotation(xp,yp,zp,xg,yg,zg,r1)
C
C compute the von mises stress
C rotate nodal displacements to local system
C
       do 270 i=1,18
          dll(i) = 0.0d+00
          ddll(i)= 0.0d+00
270    continue
       do 280 i=1,3
	 do 280 j=1,3
	   dll(i)    = dll(i)	 + r1(j,i)*v(j)
	   dll(i+3)  = dll(i+3)  + r1(j,i)*v(j+3)
	   dll(i+6)  = dll(i+6)  + r1(j,i)*v(j+6)
	   dll(i+9)  = dll(i+9)  + r1(j,i)*v(j+9)
	   dll(i+12) = dll(i+12) + r1(j,i)*v(j+12)
	   dll(i+15) = dll(i+15) + r1(j,i)*v(j+15)

	   ddll(i)    = ddll(i)    +r1(j,i)*dv(j)
	   ddll(i+3)  = ddll(i+3)  +r1(j,i)*dv(j+3)
	   ddll(i+6)  = ddll(i+6)  +r1(j,i)*dv(j+6)
	   ddll(i+9)  = ddll(i+9)  +r1(j,i)*dv(j+9)
	   ddll(i+12) = ddll(i+12) +r1(j,i)*dv(j+12)
	   ddll(i+15) = ddll(i+15) +r1(j,i)*dv(j+15)
280    continue

c      compute centroidal membrane strains
       call membra(xlp,ylp,alpha,le,rmem,status)

c      compute centroidal bending strains (curvatures (1/radius))
       call momen (xlp,ylp,lb,rmom,status)

       rnx  = 0.0d0
       rny  = 0.0d0
       rnxy = 0.0d0
       rmx  = 0.0d0
       rmy  = 0.0d0
       rmxy = 0.0d0

       drnx  = 0.0d0
       drny  = 0.0d0
       drnxy = 0.0d0
       drmx  = 0.0d0
       drmy  = 0.0d0
       drmxy = 0.0d0

       do 290 j=1,18
	  rnx  =  rnx + rmem(j,1)*dll(j)
	  rny  =  rny + rmem(j,2)*dll(j)
	  rnxy = rnxy + rmem(j,3)*dll(j)
          rmx  =  rmx + rmom(j,1)*dll(j)
          rmy  =  rmy + rmom(j,2)*dll(j)
          rmxy = rmxy + rmom(j,3)*dll(j)

	  drnx  =  drnx + rmem(j,1)*ddll(j)
	  drny  =  drny + rmem(j,2)*ddll(j)
	  drnxy = drnxy + rmem(j,3)*ddll(j)
          drmx  =  drmx + rmom(j,1)*ddll(j)
          drmy  =  drmy + rmom(j,2)*ddll(j)
          drmxy = drmxy + rmom(j,3)*ddll(j)
290    continue
C
C     compute centroidal stress resultants

C     membrane resultants
      rnnx  = dm(1,1)*rnx + dm(1,2)*rny + dm(1,3)*rnxy
      rnny  = dm(2,1)*rnx + dm(2,2)*rny + dm(2,3)*rnxy
      rnnxy = dm(3,1)*rnx + dm(3,2)*rny + dm(3,3)*rnxy 

      drnnx  = ddm(1,1)*rnx + ddm(1,2)*rny + ddm(1,3)*rnxy
      drnny  = ddm(2,1)*rnx + ddm(2,2)*rny + ddm(2,3)*rnxy
      drnnxy = ddm(3,1)*rnx + ddm(3,2)*rny + ddm(3,3)*rnxy 
      drnnx  = drnnx +  dm(1,1)*drnx + dm(1,2)*drny + dm(1,3)*drnxy
      drnny  = drnny +  dm(2,1)*drnx + dm(2,2)*drny + dm(2,3)*drnxy
      drnnxy = drnnxy+  dm(3,1)*drnx + dm(3,2)*drny + dm(3,3)*drnxy 

C     bending resultants
      rmmx  = db(1,1)*rmx + db(1,2)*rmy + db(1,3)*rmxy
      rmmy  = db(2,1)*rmx + db(2,2)*rmy + db(2,3)*rmxy
      rmmxy = db(3,1)*rmx + db(3,2)*rmy + db(3,3)*rmxy

      drmmx  = ddb(1,1)*rmx + ddb(1,2)*rmy + ddb(1,3)*rmxy
      drmmy  = ddb(2,1)*rmx + ddb(2,2)*rmy + ddb(2,3)*rmxy
      drmmxy = ddb(3,1)*rmx + ddb(3,2)*rmy + ddb(3,3)*rmxy
      drmmx  = drmmx  +db(1,1)*drmx + db(1,2)*drmy + db(1,3)*drmxy
      drmmy  = drmmy  +db(2,1)*drmx + db(2,2)*drmy + db(2,3)*drmxy
      drmmxy = drmmxy +db(3,1)*drmx + db(3,2)*drmy + db(3,3)*drmxy

      t   = abs(t)
      dt  = abs(dt)
      t2  = t*t
      t3  = t*t*t

      rnxt  =  rnnx/t
      rnyt  =  rnny/t
      rnxyt = rnnxy/t

      drnxt  =  drnnx/t -  rnnx/t2*dt
      drnyt  =  drnny/t -  rnny/t2*dt
      drnxyt = drnnxy/t - rnnxy/t2*dt

      rmxt  = 6.0* rmmx / t2
      rmyt  = 6.0* rmmy / t2
      rmxyt = 6.0*rmmxy / t2

      drmxt  = 6.0*( drmmx/t2 -  2.0d0*dt*rmmx / t3)
      drmyt  = 6.0*( drmmy/t2 -  2.0d0*dt*rmmy / t3)
      drmxyt = 6.0*(drmmxy/t2 - 2.0d0*dt*rmmxy / t3)


C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT TOP SURFACE (DEFAULT)
      tsx  =  rnxt +  rmxt
      tsy  =  rnyt +  rmyt
      tsxy = rnxyt + rmxyt

      dtsx  =  drnxt +  drmxt
      dtsy  =  drnyt +  drmyt
      dtsxy = drnxyt + drmxyt

C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT MEDIAN SURFACE
      msx  = rnxt
      msy  = rnyt
      msxy = rnxyt

      dmsx  = drnxt
      dmsy  = drnyt
      dmsxy = drnxyt

C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT BOTTOM SURFACE
      bsx  =  rnxt -  rmxt
      bsy  =  rnyt -  rmyt
      bsxy = rnxyt - rmxyt

      dbsx  =  drnxt -  drmxt
      dbsy  =  drnyt -  drmyt
      dbsxy = drnxyt - drmxyt

C ... SET sz = 0 TO REMIND USER OF THIS ASSUMPTION
C      sz  = 0.0d0
C      dersz = 0.0d0

C ... COMPUTE AVERAGE HYDROSTATIC STRESS
      ts0 = (tsx + tsy)/3.0d0
      ms0 = (msx + msy)/3.0d0
      bs0 = (bsx + bsy)/3.0d0

      dts0 = (dtsx + dtsy)/3.0d0
      dms0 = (dmsx + dmsy)/3.0d0
      dbs0 = (dbsx + dbsy)/3.0d0

C ... COMPUTE DEVIATORIC STRESSES
      tdsx  = tsx - ts0
      tdsy  = tsy - ts0
      tdsz  = -1.0d0 * ts0 
      dtdsx  = dtsx - dts0
      dtdsy  = dtsy - dts0
      dtdsz  = -1.0d0 * dts0 

      mdsx  = msx - ms0
      mdsy  = msy - ms0
      mdsz  = -1.0d0 * ms0 
      dmdsx  = dmsx - dms0
      dmdsy  = dmsy - dms0
      dmdsz  = -1.0d0 * dms0 

      bdsx  = bsx - bs0
      bdsy  = bsy - bs0
      bdsz  = -1.0d0 * bs0 
      dbdsx  = dbsx - dbs0
      dbdsy  = dbsy - dbs0
      dbdsz  = -1.0d0 * dbs0 

C
C ... RCFEM IGNORED THE DEVIATORIC Z COMPONENT: 
C     dsz  = 0.0

C ... COMPUTE J2

      dsq3   = dsqrt(3.0d0)

      tj2tmp = 0.5d0*((tdsx*tdsx)+(tdsy*tdsy)+(tdsz*tdsz))+(tsxy*tsxy)
      mj2tmp = 0.5d0*((mdsx*mdsx)+(mdsy*mdsy)+(mdsz*mdsz))+(msxy*msxy)
      bj2tmp = 0.5d0*((bdsx*bdsx)+(bdsy*bdsy)+(bdsz*bdsz))+(bsxy*bsxy)
      tj2    = dsq3 * dsqrt(tj2tmp)
      mj2    = dsq3 * dsqrt(mj2tmp)
      bj2    = dsq3 * dsqrt(bj2tmp)

      dtj2tmp =(dtdsx*tdsx)+(dtdsy*tdsy)+(dtdsz*tdsz)+2.0d0*(dtsxy*tsxy)
      dmj2tmp =(dmdsx*mdsx)+(dmdsy*mdsy)+(dmdsz*mdsz)+2.0d0*(dmsxy*msxy)
      dbj2tmp =(dbdsx*bdsx)+(dbdsy*bdsy)+(dbdsz*bdsz)+2.0d0*(dbsxy*bsxy)

      if (tj2.eq.0.0) then
        dtj2 = 0.0D0
      else
       dtj2  = 1.5d0*dtj2tmp/tj2
      endif

      if (mj2.eq.0.0) then
        dmj2 = 0.0D0
      else
        dmj2 = 1.5d0*dmj2tmp/mj2
      endif

      if (bj2.eq.0.0) then
        dbj2 = 0.0D0
      else
        dbj2 = 1.5d0*dbj2tmp/bj2
      endif

      intj2  = (tj2  + 4.0d0*mj2  + bj2 )/6.0d0
      dintj2 = (dtj2 + 4.0d0*dmj2 + dbj2)/6.0d0
      
      if (iarea.eq.0) then
        vol  = 1.0d0
        dvol = 0.0d0
      endif

      vmint = vol*pow((intj2/sigbar),fac)

      dfac = fac - 1.0d0
      dvmint = dvol*pow((intj2/sigbar),fac) + 
     *         vol*fac*(dintj2/sigbar)*pow((intj2/sigbar),dfac)

c                                                          cvar.ne.0       
       else
C
C triangular dimension variables
C
	x21 = xl(2) - xl(1)
	y21 = yl(2) - yl(1)
        z21 = zl(2) - zl(1)
	x32 = xl(3) - xl(2)
	y32 = yl(3) - yl(2)
        z32 = zl(3) - zl(2)
	x13 = xl(1) - xl(3)
	y13 = yl(1) - yl(3)
        z13 = zl(1) - zl(3)

	dx21 = dxl(2) - dxl(1)
	dy21 = dyl(2) - dyl(1)
        dz21 = dzl(2) - dzl(1)
	dx32 = dxl(3) - dxl(2)
	dy32 = dyl(3) - dyl(2)
        dz32 = dzl(3) - dzl(2)
	dx13 = dxl(1) - dxl(3)
	dy13 = dyl(1) - dyl(3)
        dz13 = dzl(1) - dzl(3)
C  triangle in space : we compute the length of one side and the distance of the
C  opposing node to that side to compute the area
       rx   = x21
       ry   = y21
       rz   = z21
       bx   = x32
       by   = y32
       bz   = z32
       rlr  = dsqrt( rx*rx + ry*ry + rz*rz )
       rlb  = dsqrt( bx*bx + by*by + bz*bz )
       bpr  = abs(rx * bx + ry * by + rz *bz )/rlr
       area = 0.5d0*rlr*dsqrt(rlb*rlb-bpr*bpr)
       vol  = area*t

       drx   = dx21
       dry   = dy21
       drz   = dz21
       dbx   = dx32
       dby   = dy32
       dbz   = dz32
       drlr  = (rx*drx + ry*dry + rz*drz)/rlr
       drlb  = (bx*dbx + by*dby + bz*dbz )/rlb
       dbpr  = (abs(drx * bx + dry * by + drz *bz + 
     *          rx * dbx + ry * dby + rz *dbz ))/rlr - 
     *         (abs(rx * bx + ry * by + rz *bz ))/(rlr*rlr)*drlr
       sqrl  = dsqrt(rlb*rlb-bpr*bpr)
       dsqrl = (rlb*drlb-bpr*dbpr)/sqrl
       darea = 0.5d0*drlr*sqrl + 0.5d0*rlr*dsqrl
       dvol  = darea*t + area*dt
C
C Direction cosines of the local system . X' is directed parallel to the 
C 2-1 side Z' is the external normal (counterclockwise). 
C Y' computed as Z' x X'
C
       xp(1) = x21/rlr
       xp(2) = y21/rlr
       xp(3) = z21/rlr 
C
       dxp(1) = dx21/rlr - x21/(rlr*rlr)*drlr
       dxp(2) = dy21/rlr - y21/(rlr*rlr)*drlr
       dxp(3) = dz21/rlr - z21/(rlr*rlr)*drlr
C
C Z' local axis
C
       zp(1) = y21 * z32 - z21 * y32
       zp(2) = z21 * x32 - x21 * z32
       zp(3) = x21 * y32 - y21 * x32
       zptmp(1) = zp(1)
       zptmp(2) = zp(2)
       zptmp(3) = zp(3)
       zlr   = dsqrt( zp(1)*zp(1) + zp(2)*zp(2) + zp(3)*zp(3) )
       zp(1) = zp(1)/zlr
       zp(2) = zp(2)/zlr
       zp(3) = zp(3)/zlr

       dzp(1) = dy21 * z32 - dz21 * y32 + y21 * dz32 - z21 * dy32
       dzp(2) = dz21 * x32 - dx21 * z32 + z21 * dx32 - x21 * dz32
       dzp(3) = dx21 * y32 - dy21 * x32 + x21 * dy32 - y21 * dx32
       dzlr   = (zptmp(1)*dzp(1)+zptmp(2)*dzp(2)+zptmp(3)*dzp(3))/zlr
       dzp(1) = dzp(1)/zlr - zptmp(1)/(zlr*zlr)*dzlr
       dzp(2) = dzp(2)/zlr - zptmp(2)/(zlr*zlr)*dzlr
       dzp(3) = dzp(3)/zlr - zptmp(3)/(zlr*zlr)*dzlr
C
C Y' local axis
C
       yp(1) = zp(2) * xp(3) - zp(3) * xp(2)
       yp(2) = zp(3) * xp(1) - zp(1) * xp(3)
       yp(3) = zp(1) * xp(2) - zp(2) * xp(1)
       yptmp(1) = yp(1)
       yptmp(2) = yp(2)
       yptmp(3) = yp(3)
       ylr   = dsqrt( yp(1)*yp(1) + yp(2)*yp(2) + yp(3)*yp(3) )
       yp(1) = yp(1)/ylr
       yp(2) = yp(2)/ylr
       yp(3) = yp(3)/ylr

       dyp(1) = dzp(2)*xp(3) -dzp(3)*xp(2) +zp(2)*dxp(3) -zp(3)*dxp(2)
       dyp(2) = dzp(3)*xp(1) -dzp(1)*xp(3) +zp(3)*dxp(1) -zp(1)*dxp(3)
       dyp(3) = dzp(1)*xp(2) -dzp(2)*xp(1) +zp(1)*dxp(2) -zp(2)*dxp(1)
       dylr   = (yptmp(1)*dyp(1) +yptmp(2)*dyp(2)+yptmp(3)*dyp(3))/ylr
       dyp(1) = dyp(1)/ylr - yptmp(1)/(ylr*ylr)*dylr
       dyp(2) = dyp(2)/ylr - yptmp(2)/(ylr*ylr)*dylr
       dyp(3) = dyp(3)/ylr - yptmp(3)/(ylr*ylr)*dylr
C
C center of gravity
C
       xcg = (xl(1) + xl(2) + xl(3))/3.0d+00
       ycg = (yl(1) + yl(2) + yl(3))/3.0d+00
       zcg = (zl(1) + zl(2) + zl(3))/3.0d+00

       dxcg = (dxl(1) + dxl(2) + dxl(3))/3.0d+00
       dycg = (dyl(1) + dyl(2) + dyl(3))/3.0d+00
       dzcg = (dzl(1) + dzl(2) + dzl(3))/3.0d+00
C
C computing local coordinates
C
       do 143 i=1,3
         xlcg   = xl(i) - xcg
         ylcg   = yl(i) - ycg
         zlcg   = zl(i) - zcg
         xlp(i) = xp(1) * xlcg + xp(2) * ylcg + xp(3) * zlcg
         ylp(i) = yp(1) * xlcg + yp(2) * ylcg + yp(3) * zlcg
         zlp(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg

         dxlcg   = dxl(i) - dxcg
         dylcg   = dyl(i) - dycg
         dzlcg   = dzl(i) - dzcg
         dxlp(i) = dxp(1)*xlcg + dxp(2)*ylcg + dxp(3)*zlcg
         dylp(i) = dyp(1)*xlcg + dyp(2)*ylcg + dyp(3)*zlcg
         dzlp(i) = dzp(1)*xlcg + dzp(2)*ylcg + dzp(3)*zlcg
         dxlp(i) = dxlp(i) + xp(1)*dxlcg + xp(2)*dylcg + xp(3)*dzlcg
         dylp(i) = dylp(i) + yp(1)*dxlcg + yp(2)*dylcg + yp(3)*dzlcg
         dzlp(i) = dzlp(i) + zp(1)*dxlcg + zp(2)*dylcg + zp(3)*dzlcg
 143    continue
C
C Set Global axes
C
       xg(1) = 1.0d0
       xg(2) = 0.0d0
       xg(3) = 0.0d0
       yg(1) = 0.0d0
       yg(2) = 1.0d0
       yg(3) = 0.0d0
       zg(1) = 0.0d0
       zg(2) = 0.0d0
       zg(3) = 1.0d0
C
C computing nodal rotation matrix
C
       call g_rotation(xp,dxp,yp,dyp,zp,dzp,xg,yg,zg,r1,dr1)
C
C compute the von mises stress
C rotate nodal displacements to local system
C
       do 1270 i=1,18
          dll(i) = 0.0d+00
          ddll(i)= 0.0d+00
1270    continue
       do 1280 i=1,3
	 do 1280 j=1,3
	   dll(i)    = dll(i)	 + r1(j,i)*v(j)
	   dll(i+3)  = dll(i+3)  + r1(j,i)*v(j+3)
	   dll(i+6)  = dll(i+6)  + r1(j,i)*v(j+6)
	   dll(i+9)  = dll(i+9)  + r1(j,i)*v(j+9)
	   dll(i+12) = dll(i+12) + r1(j,i)*v(j+12)
	   dll(i+15) = dll(i+15) + r1(j,i)*v(j+15)

	   ddll(i)    = ddll(i)    +dr1(j,i)*v(j)    +r1(j,i)*dv(j)
	   ddll(i+3)  = ddll(i+3)  +dr1(j,i)*v(j+3)  +r1(j,i)*dv(j+3)
	   ddll(i+6)  = ddll(i+6)  +dr1(j,i)*v(j+6)  +r1(j,i)*dv(j+6)
	   ddll(i+9)  = ddll(i+9)  +dr1(j,i)*v(j+9)  +r1(j,i)*dv(j+9)
	   ddll(i+12) = ddll(i+12) +dr1(j,i)*v(j+12) +r1(j,i)*dv(j+12)
	   ddll(i+15) = ddll(i+15) +dr1(j,i)*v(j+15) +r1(j,i)*dv(j+15)
1280    continue

c      compute centroidal membrane strains
       call g_membra(xlp,dxlp,ylp,dylp,alpha,le,rmem,drmem,status)

c      compute centroidal bending strains (curvatures (1/radius))
       call g_momen (xlp,dxlp,ylp,dylp,lb,rmom,drmom,status)

       rnx  = 0.0d0
       rny  = 0.0d0
       rnxy = 0.0d0
       rmx  = 0.0d0
       rmy  = 0.0d0
       rmxy = 0.0d0

       drnx  = 0.0d0
       drny  = 0.0d0
       drnxy = 0.0d0
       drmx  = 0.0d0
       drmy  = 0.0d0
       drmxy = 0.0d0

       do 1290 j=1,18
	  rnx  =  rnx + rmem(j,1)*dll(j)
	  rny  =  rny + rmem(j,2)*dll(j)
	  rnxy = rnxy + rmem(j,3)*dll(j)
          rmx  =  rmx + rmom(j,1)*dll(j)
          rmy  =  rmy + rmom(j,2)*dll(j)
          rmxy = rmxy + rmom(j,3)*dll(j)

	  drnx  =  drnx + drmem(j,1)*dll(j)+ rmem(j,1)*ddll(j)
	  drny  =  drny + drmem(j,2)*dll(j)+ rmem(j,2)*ddll(j)
	  drnxy = drnxy + drmem(j,3)*dll(j)+ rmem(j,3)*ddll(j)
          drmx  =  drmx + drmom(j,1)*dll(j)+ rmom(j,1)*ddll(j)
          drmy  =  drmy + drmom(j,2)*dll(j)+ rmom(j,2)*ddll(j)
          drmxy = drmxy + drmom(j,3)*dll(j)+ rmom(j,3)*ddll(j)
1290    continue
C
C     compute centroidal stress resultants

C     membrane resultants
      rnnx  = dm(1,1)*rnx + dm(1,2)*rny + dm(1,3)*rnxy
      rnny  = dm(2,1)*rnx + dm(2,2)*rny + dm(2,3)*rnxy
      rnnxy = dm(3,1)*rnx + dm(3,2)*rny + dm(3,3)*rnxy 

      drnnx  = ddm(1,1)*rnx + ddm(1,2)*rny + ddm(1,3)*rnxy
      drnny  = ddm(2,1)*rnx + ddm(2,2)*rny + ddm(2,3)*rnxy
      drnnxy = ddm(3,1)*rnx + ddm(3,2)*rny + ddm(3,3)*rnxy 
      drnnx  = drnnx +  dm(1,1)*drnx + dm(1,2)*drny + dm(1,3)*drnxy
      drnny  = drnny +  dm(2,1)*drnx + dm(2,2)*drny + dm(2,3)*drnxy
      drnnxy = drnnxy+  dm(3,1)*drnx + dm(3,2)*drny + dm(3,3)*drnxy 

C     bending resultants
      rmmx  = db(1,1)*rmx + db(1,2)*rmy + db(1,3)*rmxy
      rmmy  = db(2,1)*rmx + db(2,2)*rmy + db(2,3)*rmxy
      rmmxy = db(3,1)*rmx + db(3,2)*rmy + db(3,3)*rmxy

      drmmx  = ddb(1,1)*rmx + ddb(1,2)*rmy + ddb(1,3)*rmxy
      drmmy  = ddb(2,1)*rmx + ddb(2,2)*rmy + ddb(2,3)*rmxy
      drmmxy = ddb(3,1)*rmx + ddb(3,2)*rmy + ddb(3,3)*rmxy
      drmmx  = drmmx  +db(1,1)*drmx + db(1,2)*drmy + db(1,3)*drmxy
      drmmy  = drmmy  +db(2,1)*drmx + db(2,2)*drmy + db(2,3)*drmxy
      drmmxy = drmmxy +db(3,1)*drmx + db(3,2)*drmy + db(3,3)*drmxy

      t   = abs(t)
      dt  = abs(dt)
      t2  = t*t
      t3  = t*t*t

      rnxt  =  rnnx/t
      rnyt  =  rnny/t
      rnxyt = rnnxy/t

      drnxt  =  drnnx/t -  rnnx/t2*dt
      drnyt  =  drnny/t -  rnny/t2*dt
      drnxyt = drnnxy/t - rnnxy/t2*dt

      rmxt  = 6.0* rmmx / t2
      rmyt  = 6.0* rmmy / t2
      rmxyt = 6.0*rmmxy / t2

      drmxt  = 6.0*( drmmx/t2 -  2.0d0*dt*rmmx / t3)
      drmyt  = 6.0*( drmmy/t2 -  2.0d0*dt*rmmy / t3)
      drmxyt = 6.0*(drmmxy/t2 - 2.0d0*dt*rmmxy / t3)


C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT TOP SURFACE (DEFAULT)
      tsx  =  rnxt +  rmxt
      tsy  =  rnyt +  rmyt
      tsxy = rnxyt + rmxyt

      dtsx  =  drnxt +  drmxt
      dtsy  =  drnyt +  drmyt
      dtsxy = drnxyt + drmxyt

C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT MEDIAN SURFACE
      msx  = rnxt
      msy  = rnyt
      msxy = rnxyt

      dmsx  = drnxt
      dmsy  = drnyt
      dmsxy = drnxyt

C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT BOTTOM SURFACE
      bsx  =  rnxt -  rmxt
      bsy  =  rnyt -  rmyt
      bsxy = rnxyt - rmxyt

      dbsx  =  drnxt -  drmxt
      dbsy  =  drnyt -  drmyt
      dbsxy = drnxyt - drmxyt

C ... SET sz = 0 TO REMIND USER OF THIS ASSUMPTION
C      sz  = 0.0d0
C      dersz = 0.0d0

C ... COMPUTE AVERAGE HYDROSTATIC STRESS
      ts0 = (tsx + tsy )/3.0d0
      ms0 = (msx + msy )/3.0d0
      bs0 = (bsx + bsy )/3.0d0

      dts0 = (dtsx + dtsy)/3.0d0
      dms0 = (dmsx + dmsy)/3.0d0
      dbs0 = (dbsx + dbsy)/3.0d0

C ... COMPUTE DEVIATORIC STRESSES
      tdsx  = tsx - ts0
      tdsy  = tsy - ts0
      tdsz  = -1.0d0 * ts0 
      dtdsx  = dtsx - dts0
      dtdsy  = dtsy - dts0
      dtdsz  = -1.0d0 * dts0 

      mdsx  = msx - ms0
      mdsy  = msy - ms0
      mdsz  = -1.0d0 * ms0 
      dmdsx  = dmsx - dms0
      dmdsy  = dmsy - dms0
      dmdsz  = -1.0d0 * dms0 

      bdsx  = bsx - bs0
      bdsy  = bsy - bs0
      bdsz  = -1.0d0 * bs0 
      dbdsx  = dbsx - dbs0
      dbdsy  = dbsy - dbs0
      dbdsz  = -1.0d0 * dbs0 

C
C ... RCFEM IGNORED THE DEVIATORIC Z COMPONENT: 
C     dsz  = 0.0

C ... COMPUTE J2

      dsq3   = dsqrt(3.0d0)

      tj2tmp = 0.5d0*((tdsx*tdsx)+(tdsy*tdsy)+(tdsz*tdsz))+(tsxy*tsxy)
      mj2tmp = 0.5d0*((mdsx*mdsx)+(mdsy*mdsy)+(mdsz*mdsz))+(msxy*msxy)
      bj2tmp = 0.5d0*((bdsx*bdsx)+(bdsy*bdsy)+(bdsz*bdsz))+(bsxy*bsxy)
      tj2    = dsq3 * dsqrt(tj2tmp)
      mj2    = dsq3 * dsqrt(mj2tmp)
      bj2    = dsq3 * dsqrt(bj2tmp)

      dtj2tmp =(dtdsx*tdsx)+(dtdsy*tdsy)+(dtdsz*tdsz)+2.0d0*(dtsxy*tsxy)
      dmj2tmp =(dmdsx*mdsx)+(dmdsy*mdsy)+(dmdsz*mdsz)+2.0d0*(dmsxy*msxy)
      dbj2tmp =(dbdsx*bdsx)+(dbdsy*bdsy)+(dbdsz*bdsz)+2.0d0*(dbsxy*bsxy)

      if (tj2.eq.0.0) then
        dtj2 = 0.0D0
      else
       dtj2  = 1.5d0*dtj2tmp/tj2
      endif

      if (mj2.eq.0.0) then
        dmj2 = 0.0D0
      else
        dmj2 = 1.5d0*dmj2tmp/mj2
      endif

      if (bj2.eq.0.0) then
        dbj2 = 0.0D0
      else
        dbj2 = 1.5d0*dbj2tmp/bj2
      endif

      intj2  = (tj2  + 4.0d0*mj2  + bj2 )/6.0d0
      dintj2 = (dtj2 + 4.0d0*dmj2 + dbj2)/6.0d0

      if (iarea.eq.0) then
        vol  = 1.0d0
        dvol = 0.0d0
      endif
      
      vmint = vol*pow((intj2/sigbar),fac)

      dfac = fac - 1.0d0
      dvmint = dvol*pow((intj2/sigbar),fac) + 
     *         vol*fac*(dintj2/sigbar)*pow((intj2/sigbar),dfac)

      endif
      
      return
      end
