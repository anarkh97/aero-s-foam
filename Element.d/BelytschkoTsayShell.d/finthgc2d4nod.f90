! ==================================
! hourglass control force: 2d 4 node quad
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine          elefhgc2d4nod                 (optpty,prmhgc,ematpro,ecord,edisp,evelo, efhgc)
! 2.  subroutine          getgamma4nod                  (ecordloc, gamma)
!
! =========================================================================================================



subroutine elefhgc2d4nod(optpty,prmhgc,ematpro,ecord,edisp,evelo, efhgc)
  !=======================================================================
  !  elefhgc2d4nod = compute hourglass control force for 2d 4 node quad
  !
  !                  feature:
  !                  -------
  !                  frame invariant is considered
  !                  for transient analysis, corotational coordinate is used
  !
  !                  note:
  !                  ----
  !                  daniel and belytschko, IJNME, 2005, vol. , pp. -
  !                  supperession of spurious intermediate frequency modes in
  !                  under-integrated elements by combined stiffness/viscous stabilization             
  !
  !                  belytschko and bachrach, IJNME, 1986, vol. 54, pp. 279-301
  !                  efficient implimemtation of quadrilaterals with high coarse-mesh accuracy
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optpty : 2d problem type
  !           optpty=1 : plane-stress
  !           optpty=2 : plane-strain
  !
  !  prmhgc(*) : hourglass control parameter
  !
  !  ematpro(*) : material property
  !
  !  ecord(2,4) : element nodal coordinate data
  !                     
  !  edisp(2,4) : element nodal displacement data
  !
  !  evelo(2,4) : element nodal displacement 
  !
  !  output:
  !  ------
  !  efhgc(8,1) : hourglass mode control force
  !                            
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optpty
  real(8), dimension(*), intent(in) :: prmhgc
  real(8), dimension(*), intent(in) :: ematpro

  real(8), dimension(2,4), intent(in) :: ecord
  real(8), dimension(2,4), intent(in) :: edisp
  real(8), dimension(2,4), intent(in) :: evelo

  real(8), dimension(8,1), intent(out) :: efhgc
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: rk
  real(8) :: young, poiss, denst, thick

  real(8) :: e, nu, mu
  real(8), dimension(4) :: h

  real(8), dimension(2,4) :: ecurn
  
  real(8) :: bxmax, bymay, cxmdx, dymcy
  real(8) :: xlab, xlcd
  real(8) :: theta
  real(8), dimension(2,4) :: ecordloc

  real(8), dimension(4,1) :: gamma
  real(8), dimension(4,4) :: ggt

  real(8), dimension(2,4) :: gqpoin
  real(8), dimension(4) :: gqweigt

  real(8), dimension(4) :: shap
  real(8), dimension(2,4) :: deriv
  
  real(8), dimension(2,4) :: cartd, xjaci
  real(8) :: djacob

  real(8) :: psi,eta
  real(8) :: hx, hy, hxx, hyy, hxy
  real(8) :: hstar, c1, c2, c3
  real(8) :: cs, rc, area

  real(8), dimension(8,8) :: khgc
  real(8), dimension(8,8) :: trans
  real(8), dimension(8,8) :: glbkhgc
  real(8), dimension(8,8) :: temp1
  real(8), dimension(8,1) :: hgcdisp
  real(8), dimension(8,1) :: efhgc2

  ! loop index
  integer :: inode, jnode, igaus
  ! ====================================

  ! initialize
  efhgc(:,:)= 0.0d0


  ! -------------------------------------------
  ! stiffness propotional stabilization parameter
  rk= prmhgc(1)

  ! get material properties
  young= ematpro(1)
  poiss= ematpro(2)
  denst= ematpro(3)
  thick= ematpro(20) 
  ! -------------------------------------------

  ! compute material constant: e, mu, nu
  ! -------------------------
  select case(optpty)
  case(1) ! plane stress
     e= young
     nu= poiss
     mu= young/( 2.0d0*(1.0d0 + poiss) )

  case(2) ! plane strain
     e= young/(1.0d0 - poiss**2)
     nu= poiss/(1.0d0 - poiss)
     mu= young/( 2.0d0*(1.0d0 + poiss) )

  end select


  ! define hourglass vector: (Eq. 8)
  ! -----------------------
  data h /1.0d0, -1.0d0, 1.0d0, -1.0d0/


  ! compute local coordinate
  ! --------------------------
  ! compute current coordinate: x= u + X
  ecurn(1:2,1:4)= ecord(1:2,1:4) + edisp(1:2,1:4)

  ! compute rotation angle
  bxmax= ( ecurn(1,3) + ecurn(1,2) - ecurn(1,4) - ecurn(1,1) ) / 2.0d0
  bymay= ( ecurn(2,3) + ecurn(2,2) - ecurn(2,4) - ecurn(2,1) ) / 2.0d0
  cxmdx= ( ecurn(1,2) + ecurn(1,1) - ecurn(1,3) - ecurn(1,4) ) / 2.0d0
  dymcy= ( ecurn(2,3) + ecurn(2,4) - ecurn(2,2) - ecurn(2,1) ) / 2.0d0
  xlab= dsqrt(bxmax**2 + bymay**2)
  xlcd= dsqrt(cxmdx**2 + dymcy**2)
  theta= ( xlab*atan(bymay/bxmax) + xlcd*atan(cxmdx/dymcy) ) / ( xlab + xlcd )

  ! get element local coordinate
  do inode=1, 4
     ecordloc(1,inode)= ecurn(1,inode)*dcos(theta) + ecurn(2,inode)*dsin(theta)
     ecordloc(2,inode)= ecurn(1,inode)*(-dsin(theta)) + ecurn(2,inode)*dcos(theta)
  end do

  ! compute gamma projector
  ! -----------------------
  call getgamma4nod(ecordloc, gamma)
     ! input : ecordloc
     ! output : gamma

  ! compute gamma.gamma^t: ggt= gamma.gamma^t
  ! ---------------------
  call matprd(4,1,0, 4,1,1, 4,4, gamma,gamma, ggt)
     ! input : 4,1,0, 4,1,1, 4,4, gamma,gamma
     ! output : ggt


  ! ------------------------------------------------------------------------------------
  ! compute Hxx, Hyy, Hxy: using optimal gq rule
  ! ---------------------
  hxx= 0.0d0 ! initialize
  hyy= 0.0d0
  hxy= 0.0d0

  ! set optimal gq for 4 node quad
  ! ----------------------
  call getgqele(2,2,2,4, gqpoin,gqweigt)
     ! input : 2(optele=quad),2(ngqdim),2(ngqpt),4(mgaus)
     ! output : gqpoin,gqweigt


  do igaus=1, 4

     psi= gqpoin(1,igaus)
     eta= gqpoin(2,igaus)

     ! compute 2d shape function
     call getshape2d(2,4,psi,eta, shap,deriv)
        ! input : 2(optele=quad),4(nnode),psi,eta
        ! output : shap,deriv

     ! compute inverse of jacobian matrix: current coordinate
     call jacob2(2,4,deriv,ecordloc, djacob,xjaci,cartd)
        ! input : 2(ndime),4(nnode),deriv,ecordloc
        ! output : djacb,xjaci,cartd

! ### working
!if ( djacob <= 0.0d0 ) then
!  efhgc(:,:)= 0.0d0
!  return
!end if


     ! compute hx= a (psi*eta) / a x, hy= a (psi*eta) / a y
     hx= xjaci(1,1)*eta + xjaci(1,2)*psi ! hx= a (psi*eta) / a x= (a psi / a x)*eta + (a eta/ a x)*psi
     hy= xjaci(2,1)*eta + xjaci(2,2)*psi ! hy= a (psi*eta) / a y= (a psi / a y)*eta + (a eta/ a y)*psi

     ! compute hij=int( hi*hj det(J) ) : (Eq. 10)
     hxx= hxx + hx*hx *gqweigt(igaus) *djacob
     hyy= hyy + hy*hy *gqweigt(igaus) *djacob
     hxy= hxy + hx*hy *gqweigt(igaus) *djacob

  end do
  ! ------------------------------------------------------------------------------------

  ! compute coefficient: (Eq. 7) (QBI: see IJNME, 1986, vol. 54, pp. 279-301 )
  ! -------------------
! ### working
!if ( hxx*hyy - (nu*hxy)**2 == 0.0d0 ) then
!  efhgc(:,:)= 0.0d0
!  return
!end if

  hstar= (hxx*hyy) / ( hxx*hyy - (nu*hxy)**2 )
  c1= hxx*hstar
  c2= nu*hxy*hstar
  c3= hyy*hstar


  ! compute viscous damping parameter: (Eq. 18)
  ! ---------------------------------
  cs= dsqrt(mu/denst) ! shear wave speed
  rc= dsqrt( area(4,ecordloc) ) / cs * dsqrt( 3.0d0*rk / (2.0d0*(1.0d0+nu)) ) ! viscos propotional daming parameter


  ! construct hourglass control stiffness matrix: (Eq. 6)
  ! --------------------------------------------
  khgc(1:4,1:4)= c1 *ggt(1:4,1:4) ! row: 1-4, col: 1-4
  khgc(1:4,5:8)= c2 *ggt(1:4,1:4) ! row: 1-4, col: nnode+1-nnode+4
  khgc(5:8,1:4)= c2 *ggt(1:4,1:4) ! row: nnode+1-nnode+4, col: 1-4 
  khgc(5:8,5:8)= c3 *ggt(1:4,1:4) ! row: nnode+1-nnode+4, col: nnode+1-nnode+4 

  khgc(1:8,1:8)= 2.0d0*mu*(1.0d0+nu)*khgc(1:8,1:8)


  ! construct transformation matrix: corresponding to  [ k_xx k_xy ] shpae stiffness matrix
  ! -------------------------------                    [ k_xy k_yy ]
  do inode=1, 4
     trans(inode,inode)= dcos(theta) ! row: 1-4, col: 1-4
     trans(inode,4+inode)= dsin(theta) ! row: 1-4, col: nnode+1-nnode+4
     trans(4+inode,inode)= -dsin(theta) ! row: nnode+1-nnode+4, col: 1-4 
     trans(4+inode,4+inode)= dcos(theta) ! row: nnode+1-nnode+4, col: nnode+1-nnode+4 
  end do

  ! compute transformed stiffness matrix: gkstab= trans^t.kstab.trans
  ! ---------------------
  call matprd(8,8,1, 8,8,0, 8,8, trans,khgc, temp1) ! temp1= trans^t.khgc
     ! input : 8,8,1, 8,8,0, 8,8, trans,khgc
     ! output : temp1

  call matprd(8,8,0, 8,8,0, 8,8, temp1,trans, glbkhgc) ! glbkhgc= [temp1].trans
     ! input : 8,8,0, 8,8,0, 8,8, temp1,trans
     ! output : glbkhgc


  ! compute hourglass displacement vector: [u1x, u2x, u3x, u4x | u1y, u2y, u3y, u4y]
  ! -------------------------------------
  hgcdisp(1:4,1)= edisp(1,1:4) + rc*evelo(1,1:4)
  hgcdisp(5:8,1)= edisp(2,1:4) + rc*evelo(2,1:4)


  ! compute stabilization force: [f1x, f2x, f3x, f4x | f1y, f2y, f3y, f4y]
  ! ---------------------------
  call matprd(8,8,0, 8,1,0, 8,1, glbkhgc,hgcdisp, efhgc2)
     ! input : 8,8,0, 8,1,0, 8,1, glbkhgc,hgcdisp
     ! output : efhgc2


  ! ----------------------------------------------- 
  ! change component address of stabilization force: [f1x, f1y | f2x, f2y | f3x, f3y | f4x, f4y]
  ! -----------------------------------------------
  do inode=1, 4
     jnode=2*(inode-1)+1
     efhgc(jnode,1)= efhgc2(inode,1)
     efhgc(jnode+1,1)= efhgc2(4+inode,1)
  end do

  ! consider thickness
  efhgc(1:8,1)= efhgc(1:8,1) * thick



  return
end subroutine elefhgc2d4nod




subroutine getgamma4nod(ecordloc, gamma)
  !=======================================================================
  !  getgamma4nod = compute gamma projection operator of 4 node quad. element
  !
  !                 note:
  !                 ----
  !                 flanagan and belytschko, IJNME, 1981, vol. 17, pp. 679-706
  !                 a uniform strain hexahedron and quadrilateral
  !                 with orthogonal hourglass control
  !                 (see, eq (72))
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ecordloc(2,4) : local element nodal coordinate
  !
  !  output:
  !  ------
  !  gamma(4,1) : gamma projection operator
  !                            
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(2,4), intent(in) :: ecordloc

  real(8), dimension(4,1), intent(out) :: gamma
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: area, const
  real(8) :: x1,x2,x3,x4
  real(8) :: y1,y2,y3,y4

  ! ====================================

  ! initialize
  gamma(:,:)= 0.0d0


  ! define components
  x1= ecordloc(1,1)
  x2= ecordloc(1,2)
  x3= ecordloc(1,3)
  x4= ecordloc(1,4)

  y1= ecordloc(2,1)
  y2= ecordloc(2,2)
  y3= ecordloc(2,3)
  y4= ecordloc(2,4)

  ! compute area
  area= 0.50d0 * ( ( x3 - x1 ) * ( y4 - y2 ) + ( x2 - x4 ) * ( y3 - y1 ) )

  ! check current element configuration
  if ( area <= 0.0d0 ) then
     write(*,*) "current element has zero or negative area: getgamma4nod"
     write(nout5,*) "current element has zero or negative area: getgamma4nod"
     stop
  end if

  ! compute constant
  const= 1.0d0 / ( 4.0d0 * area )

  ! set gamm projection operator for one point integration
  gamma(1,1)= const*( x2*(y3-y4) + x3*(y4-y2) + x4*(y2-y3) )
  gamma(2,1)= const*( x3*(y1-y4) + x4*(y3-y1) + x1*(y4-y3) )
  gamma(3,1)= const*( x4*(y1-y2) + x1*(y2-y4) + x2*(y4-y1) )
  gamma(4,1)= const*( x1*(y3-y2) + x2*(y1-y3) + x3*(y2-y1) )



  return
end subroutine getgamma4nod