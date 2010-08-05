! ==================================
! internal force: bt shell / elasto viscoplastic 
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine        elefintevpbt1          (optcri,delt,ematpro,nndof,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo,eaccl, &
!                                               evar1,evoit2,evoit3, &
!                                               evar2,efint)
! 2.  subroutine        elefintevpbt0          (optcri,delt,ematpro,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo,eaccl, &
!                                               evar1,evoit2,evoit3, &
!                                               evar2,efint)
! 3.  subroutine        updevpstrsbt           (delt,ematpro,gqpoin,gqweigt,locbvec,ecord,edisp,evelo,eaccl, &
!                                               effpstrn,sigvoitloc,strnvoitloc, &
!                                               effstrs)
! 4.  subroutine        updsigevpbt1           (delt,ematpro,zeta,ecordloc,edisploc,eveloloc,eacclloc, &
!                                               effpstrn,sigvoitloc,strnvoitloc, &
!                                               effstrs)
!
! =========================================================================================================

subroutine elefintevpbt1(optcri,delt,ematpro,nndof,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo,eaccl, &
                         evar1,evoit2,evoit3, &
                         evar2,efint)
  !=======================================================================
  !  elefintevpbt1 = compute internal force matrix for bt shell
  !                   note:
  !                   ----
  !                   rotation projection for the drilling dof is considered
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optcri(*) : fracture criterion option handler
  !
  !  delt : time increment
  ! 
  !  ematpro(*) : material property
  !
  !  nndof : the total number of dof per node
  !
  !  mgaus3 : total number of gq through thickness
  !
  !  mgqpt1 : maximum number of gq point in regular element
  !
  !  gqpoin3(mgaus3), gqweigt3(mgaus3) : through thickness gq point and weight
  !
  !  ecord(3,4) : element nodal coordinate
  !
  !  edisp(nndof,4) : element nodal displacement data : d_x, d_y, d_z, r_x, r_y, [r_z]
  !
  !  evelo(nndof,4) : element nodal velocity
  !
  !  eaccl(nndof,4) : element nodal acceleration
  !
  !  inoutput:
  !  --------
  !  evar1(5,mgqpt1) : evar1(1,:) effective plastic strain
  !
  !  evoit2(6,mgqpt1) : cauchy stress( local )
  !
  !  evoit3(6,mgqpt1) : strain(local )
  !
  !  output:
  !  ------
  !  evar2(5,mgqpt1) : evar2(1,:) effective stress
  !
  !  efint(nndof*4,1) : element nodal internal force
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, dimension(*), intent(in) :: optcri
  real(8), intent(in) :: delt
  real(8), dimension(20), intent(in) :: ematpro
  integer, intent(in) :: nndof,mgaus3,mgqpt1
  real(8), dimension(mgaus3), intent(in) :: gqpoin3
  real(8), dimension(mgaus3), intent(in) :: gqweigt3
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(nndof,4), intent(in) :: edisp
  real(8), dimension(nndof,4), intent(in) :: evelo
  real(8), dimension(nndof,4), intent(in) :: eaccl

  ! ------------------------------------
  real(8), dimension(5,mgqpt1), intent(inout) :: evar1 ! update
  real(8), dimension(6,mgqpt1), intent(inout) :: evoit2
  real(8), dimension(6,mgqpt1), intent(inout) :: evoit3
  ! ------------------------------------
 
  real(8), dimension(5,mgqpt1), intent(out) :: evar2
  ! ------------------------------------
  real(8), dimension(nndof*4,1), intent(out) :: efint
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(5,4) :: edisp0
  real(8), dimension(5,4) :: evelo0
  real(8), dimension(5,4) :: eaccl0

  real(8), dimension(20,1) :: efint0
  ! ====================================

  ! initialize
  evar2(:,:)= 0.0d0
  efint(:,:)= 0.0d0


  ! rotation projection
  ! -------------------
  if ( nndof == 5 ) then ! 5 dof
     ! just use given edisp and evelo
     edisp0(1:5,1:4)= edisp(1:5,1:4)
     evelo0(1:5,1:4)= evelo(1:5,1:4)
     eaccl0(1:5,1:4)= eaccl(1:5,1:4)

  else if ( nndof == 6) then ! 6 dof
     ! project 6 dof edisp(6,4) and evelo(6,4) to 5 dof edisp(5,4) and evelo(5,4)
     call rotprojbt2(ecord,edisp,evelo,eaccl, edisp0,evelo0,eaccl0)
        ! input : ecord,edisp,evelo,eaccl
        ! output : edisp0,evelo0

  end if


  ! ----------------------------------------------------------------------------------------------
  ! 5 dof bt shell
  ! --------------
  call elefintevpbt0(optcri,delt,ematpro,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp0,evelo0,eaccl0, &
                     evar1,evoit2,evoit3, &
                     evar2,efint0)
        ! input : optcri,delt,ematpro,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp0,evelo0,eaccl0
        ! inoutput : evar1,evoit2(cauchy, local),evoit3(strain, local)
        ! output : evar2,efint0

  ! ----------------------------------------------------------------------------------------------


  ! rotation projection
  ! -------------------  
  if ( nndof == 5 ) then
     ! just use given edisp and evelo
     efint(1:20,1)= efint0(1:20,1)

  else if ( nndof == 6) then
     ! project 5 dof efint0(20,1) to 6 dof efint(24,1)
     call rotprojbt3(ecord,edisp,efint0, efint)
        ! input : ecord,edisp,efint0
        ! output : efint

  end if



  return
end subroutine elefintevpbt1




subroutine elefintevpbt0(optcri,delt,ematpro,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo,eaccl, &
                         evar1,evoit2,evoit3, &
                         evar2,efint)
  !=======================================================================
  !  elefintevpbt0 = compute element internal nodal force
  !                   of belytschko tsay shell element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optcri(*) : fracture criterion option handler
  !
  !  delt : time increment
  ! 
  !  ematpro(*) : material property
  !
  !  mgaus3 : total number of gq through thickness
  !
  !  mgqpt1 : maximum number of gq point in regular element
  !
  !  gqpoin3(mgaus3), gqweigt3(mgaus3) : through thickness gq point and weight
  !
  !  ecord(3,4) : element nodal coordinate
  !
  !  edisp(5,4) : element nodal displacement data
  !
  !  evelo(5,4) : element nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  inoutput:
  !  --------
  !  evar1(5,mgqpt1) : evar1(1,:) effective plastic strain
  !
  !  evoit2(6,mgqpt1) : cauchy stress( local )
  !
  !  evoit3(6,mgqpt1) : in plane strain( local )
  !
  !  output:
  !  ------
  !  evar2(5,mgqpt1) : evar2(1,:) effective stress
  !
  !  efint(20,1) : element nodal internal force
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, dimension(*), intent(in) :: optcri
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  integer, intent(in) :: mgaus3,mgqpt1
  real(8), dimension(mgaus3), intent(in) :: gqpoin3
  real(8), dimension(mgaus3), intent(in) :: gqweigt3
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(5,4), intent(in) :: edisp
  real(8), dimension(5,4), intent(in) :: evelo
  real(8), dimension(5,4), intent(in) :: eaccl

  ! -------------------------------------
  real(8), dimension(5,mgqpt1), intent(inout) :: evar1 ! update
  real(8), dimension(6,mgqpt1), intent(inout) :: evoit2
  real(8), dimension(6,mgqpt1), intent(inout) :: evoit3
  ! -------------------------------------

  real(8), dimension(5,mgqpt1), intent(out) :: evar2
  ! -------------------------------------
  real(8), dimension(20,1), intent(out) :: efint
  ! ====================================
  ! local variable
  ! ==============
  integer :: crityp

  real(8), dimension(3,3) :: locbvec
  real(8), dimension(3,4) :: ecurn

  real(8), dimension(6,1) :: sigvoitloc, strnvoitloc
  real(8) :: effpstrn, effstrs

  real(8), dimension(3,1) :: sigvoit2d, strnvoit2d
  real(8) :: crival, criang

  real(8), dimension(20,1) :: gqfint

  ! loop index
  integer :: igaus
  ! ====================================

  ! initialize
  evar2(:,:)= 0.0d0
  efint(:,:)= 0.0d0


  ! -------------------------------------------
  ! get fracture criterion parameters
  crityp= abs(optcri(1)) ! criterion type

  ! initialize
  crival= 0.0d0
  criang= 0.0d0
  ! -------------------------------------------

  ! get current nodal coordinate
  ecurn(1:3,1:4)= ecord(1:3,1:4) + edisp(1:3,1:4)

  ! compute co rotational local base vector: locbvec
  ! ---------------------------------------
  call getlocbvecbt(ecurn, locbvec)
     ! input : ecurn
     ! output : locbvec


  ! loop over gauss quadarture
  ! note: this loop is for through thickness integration
  !       in plane, we use 1 point rule
  do igaus=1, mgaus3 

     ! -------------------------------------
     ! get history variable
     ! --------------------
     ! effective plastic strain
     effpstrn= evar1(1,igaus)

     ! cauchy stress( local: co-rotational )
     sigvoitloc(1:6,1)= evoit2(1:6,igaus)

     ! strain(local)
     strnvoitloc(1:6,1)= evoit3(1:6,igaus)
     ! -------------------------------------


     ! update hypo stresses at gq
     ! --------------------------
     call updevpstrsbt(delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),locbvec,ecord,edisp,evelo,eaccl, &
                       effpstrn,sigvoitloc,strnvoitloc, &
                       effstrs)
        ! input : delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),locbvec,ecord,edisp,evelo,eaccl
        ! inoutput : effpstrn,sigvoitloc,strnvoitloc
        ! output : effstrs


     ! check fracture criterion
     ! ------------------------
     if ( crityp /= 0 ) then

        select case(crityp)
        case(2) ! critical effective plastic strain criterion
           ! set 2d in-plane stress
           sigvoit2d(1,1)= sigvoitloc(1,1) ! sig_xx
           sigvoit2d(2,1)= sigvoitloc(2,1) ! sig_yy
           sigvoit2d(3,1)= sigvoitloc(6,1) ! sig_xy

           call chkeffstrcri2d(0,effpstrn,sigvoit2d, crival,criang)
              ! input : 0(opt:stress),effpstrn,sigvoit2d
              ! output : crival,criang

        case default
           write(*,*) "not implemented fracture criterion: elefintevpbt0"
           write(nout5,*) "not implemented fracture criterion: elefintevpbt0"
           stop

        end select

     end if


     ! -------------------------------------
     ! set history variable
     ! --------------------
     ! effective plastic strain
     evar1(1,igaus)= effpstrn

     ! cauchy stress( local: co-rotational )
     evoit2(1:6,igaus)= sigvoitloc(1:6,1)

     ! strain(local)
     evoit3(1:6,igaus)= strnvoitloc(1:6,1)
     ! -------------------------------------

     ! set effective stress
     evar2(1,igaus)= effstrs

     ! set fracture criterion
     evar2(2,igaus)= crival

     ! set fracture angle
     evar2(3,igaus)= criang
     ! -------------------------------------


     ! compute local nodal internal forces at gq
     ! -----------------------------------------
     call gqfintbt(delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),locbvec,ecurn,sigvoitloc, gqfint)
        ! input : delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),locbvec,ecurn,sigvoitloc
        ! output : gqfint


     ! sum on internal force
     efint(1:20,1)= efint(1:20,1) + gqfint(1:20,1)

  end do


  return
end subroutine elefintevpbt0





subroutine updevpstrsbt(delt,ematpro,gqpoin,gqweigt,locbvec,ecord,edisp,evelo,eaccl, &
                        effpstrn,sigvoitloc,strnvoitloc, &
                        effstrs)
  !=======================================================================
  !  updevpstrsbt = update elasto visco plastic cauchy stress of belytschko tsay element 
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  delt : time increment
  ! 
  !  ematpro(*) : material property
  !
  !  gqpoin, gqweigt : through thickness gq point and weight
  !
  !  ecord(3,4) : global initial nodal coordinate
  !
  !  edisp(5,4) : global nodal displacement: d_x, d_y, d_z, r_x, r_y
  !
  !  evelo(5,4) : global nodal velocity
  !
  !  eaccl(5,4) : global nodal acceleration
  !
  !  inoutput:
  !  --------
  !  effpstrn : effective plastic strain
  !
  !  sigvoitloc(6,1) : co-rotational cauchy stress
  !
  !  strnvoitloc(6,1) : co-rotational in plane strain
  !
  !  output:
  !  ------
  !  effstrs : effective stress
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: gqpoin, gqweigt
  real(8), dimension(3,3), intent(in) :: locbvec
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(5,4), intent(in) :: edisp
  real(8), dimension(5,4), intent(in) :: evelo
  real(8), dimension(5,4), intent(in) :: eaccl

  ! ------------------------------------
  real(8), intent(inout) :: effpstrn
  real(8), dimension(6,1), intent(inout) :: sigvoitloc, strnvoitloc
  ! ------------------------------------

  real(8), intent(out) :: effstrs
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,4) :: ecordloc
  real(8), dimension(5,4) :: edisploc
  real(8), dimension(5,4) :: eveloloc
  real(8), dimension(5,4) :: eacclloc
  real(8), dimension(3,4) :: ecurnloc

  real(8) :: zeta

  real(8), dimension(3,1) :: ipstrnloc
  real(8), dimension(2,1) :: tsstrslocdot
  ! ====================================

  ! initialize
  effstrs= 0.0d0

  ! -------------------------------------------------------------------------------
  ! convert global vector to local vector
  ! -------------------------------------
  ! get local nodal coordinate
  call glb2locnodv(3,4,3,locbvec,ecord, ecordloc)
     ! input : 3(ndime),4(nnode),3(ntrndof),locbvec,ecord
     ! output : ecordloc

  ! get local nodal velocity
  call glb2locnodv(3,4,5,locbvec,edisp, edisploc)
     ! input : 3(ndime),4(nnode),5(ntrndof),locbvec,edisp
     ! output : edisploc

  ! get local nodal velocity
  call glb2locnodv(3,4,5,locbvec,evelo, eveloloc)
     ! input : 3(ndime),4(nnode),5(ntrndof),locbvec,evelo
     ! output : eveloloc

  ! get local nodal velocity
  call glb2locnodv(3,4,5,locbvec,eaccl, eacclloc)
     ! input : 3(ndime),4(nnode),5(ntrndof),locbvec,eaccl
     ! output : eacclloc

  ! get current nodal coordinate
  ecurnloc(1:3,1:4)= ecordloc(1:3,1:4) + edisploc(1:3,1:4)


  ! -------------------------------------------------------------------------------
  ! compute in plane stress rate and update
  ! ---------------------------------------
  ! pseudo-thickness parameter
  zeta= gqpoin

  ! extract in plane strain components
  ipstrnloc(1,1)= strnvoitloc(1,1)
  ipstrnloc(2,1)= strnvoitloc(2,1)
  ipstrnloc(3,1)= strnvoitloc(6,1)

  ! compute hypoelastic in-plane stress rate of belytschko tsay shell element
  call updsigevpbt1(delt,ematpro,zeta,ecordloc,edisploc,eveloloc,eacclloc, &
                    effpstrn,ipstrnloc,sigvoitloc, &
                    effstrs)
     ! input : delt,ematpro,zeta,ecordloc,edisploc,eveloloc,eacclloc
     ! inoutput : effpstrn,ipstrnloc,sigvoitloc
     ! output : effstrs

  ! set updated in plane strain components
  strnvoitloc(1,1)= ipstrnloc(1,1)
  strnvoitloc(2,1)= ipstrnloc(2,1)
  strnvoitloc(6,1)= ipstrnloc(3,1)

  ! -------------------------------------------------------------------------------
  ! compute tranverse shear stress rate and update
  ! ----------------------------------------------
  ! compute hypoelastic transverse shear stress rate of belytschko tsay shell element
  call getsighypobt2(ematpro,ecurnloc,eveloloc, tsstrslocdot)
     ! input : ematpro,ecurnloc,eveloloc
     ! output : tsstrslocdot

  ! update local cauchy stress
  sigvoitloc(3,1)= 0.0d0 ! sig_z (plane stress)
  sigvoitloc(4,1)= sigvoitloc(4,1) + tsstrslocdot(2,1) * delt ! sig_yz
  sigvoitloc(5,1)= sigvoitloc(5,1) + tsstrslocdot(1,1) * delt ! sig_xz
  ! -------------------------------------------------------------------------------


  return
end subroutine updevpstrsbt





subroutine updsigevpbt1(delt,ematpro,zeta,ecordloc,edisploc,eveloloc,eacclloc, &
                        effpstrn,ipstrn,sigvoitloc, &
                        effstrs)
  !=======================================================================
  !  updsigevpbt1 = compute vp stress rate of belytschko-tsay local element
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  delt : time integration step
  !
  !  ematpro(*) : material property
  !
  !  zeta : [-1, +1] : pseudo-thickness parameter
  !
  !  ecordloc(3,4) : local initial nodal coordinate
  !
  !  edisploc(5,4) : local nodal displacement: d_x, d_y, d_z, r_x, r_y
  !
  !  eveloloc(5,4) : local nodal velocity
  !
  !  eacclloc(5,4) : local nodal acceleration
  !
  !
  !  inoutput:
  !  --------
  !  effpstrn : effective plastic strain
  !
  !  ipstrn(3,1) : in plane strain
  !
  !  sigvoitloc(6,1) : voight form local cauchy stress
  !
  !  output:
  !  ------
  !  effstrs : effective stress
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: zeta
  real(8), dimension(3,4), intent(in) :: ecordloc
  real(8), dimension(5,4), intent(in) :: edisploc
  real(8), dimension(5,4), intent(in) :: eveloloc
  real(8), dimension(5,4), intent(in) :: eacclloc

  ! ------------------------------------
  real(8), intent(inout) :: effpstrn
  real(8), dimension(3,1), intent(inout) :: ipstrn
  real(8), dimension(6,1), intent(inout) :: sigvoitloc
  ! ------------------------------------

  real(8), intent(out) :: effstrs
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: thick

  real(8), dimension(3,4) :: ecurnloc
  real(8), dimension(5,4) :: mideveloloc
  real(8), dimension(3,4) :: midecurnloc

  real(8), dimension(2,4) :: midedisploc2d
  real(8), dimension(2,4) :: ecordloc2d

  real(8) :: psi, eta
  real(8), dimension(4) :: shap
  real(8), dimension(2,4) :: deriv
  real(8) :: djacob
  real(8), dimension(2,4) :: cartd
  real(8), dimension(2,2) :: ftens2d, midltens2d
  real(8), dimension(3,1) :: ipstrndot

  real(8), dimension(2,2) :: sigtensloc2d, kirtensloc2d
  real(8), dimension(2,2,2,2) :: ctantens
  ! ====================================

  ! initialize
  effstrs= 0.0d0

  ! ------------------------------------------------
  ! material properties
  thick= ematpro(20)
  ! ------------------------------------------------

  ! 3d
  ! --
  ! local current nodal coordinate: current step [n]
  ecurnloc(1:3,1:4)= ecordloc(1:3,1:4) + edisploc(1:3,1:4)

  ! linearly interpolated velocity field: mid step [n+1/2]
  mideveloloc(1:5,1:4)= eveloloc(1:5,1:4) + (0.50d0*delt)*eacclloc(1:5,1:4)
     
  ! linearly interpolated nodal coordinate: mid step [n+1/2]
  midecurnloc(1:3,1:4)= ecurnloc(1:3,1:4) + (0.50d0*delt)*mideveloloc(1:3,1:4)

  ! 2d
  ! --  
  ! linearly interpolated displacement field: mid step [n+1/2]
  midedisploc2d(1:2,1:4)= edisploc(1:2,1:4) + (0.50d0*delt)*mideveloloc(1:2,1:4)

  ! set 2d local initial nodal coordinate
  ecordloc2d(1:2,1:4)= ecordloc(1:2,1:4)


  ! ----------------------------------------------------------------
  ! set 2d cauchy stress tensor
  ! ---------------------------
  sigtensloc2d(1,1)= sigvoitloc(1,1) ! sig_x
  sigtensloc2d(2,2)= sigvoitloc(2,1) ! sig_y
  sigtensloc2d(1,2)= sigvoitloc(6,1) ! sig_xy
  sigtensloc2d(2,1)= sigvoitloc(6,1) ! sig_yx


  ! ----------------------------------------------------------------
  ! compute rate of deformation: mid step [n+1/2]
  ! ---------------------------
  ! ltens2d= [ vx,y  vx,y ]     ipstrndot= [d_x, d_y, 2d_xy]
  !          [ vy,x  vy,y ]
  call getstrndotbt1(thick,zeta,midecurnloc,mideveloloc, midltens2d,ipstrndot)
     ! input : thick,zeta,ecurnloc,eveloloc
     ! output : ltens2d,ipstrndot


  ! update strain
  ! -------------
  ipstrn(1:3,1)= ipstrn(1:3,1) + delt * ipstrndot(1:3,1)

  ! ----------------------------------------------------------------
  ! compute deformation gradient tensor
  ! -----------------------------------
  ! compute 2d shape function
  psi= 0.0d0
  eta= 0.0d0
  call getshape2d(3,4,psi,eta, shap,deriv)
     ! input : 3(optele),4(nnode),psi,eta
     ! output : shap,deriv

  ! compute determinant of jacobian
  call jacob1(2,4,deriv,ecordloc2d, djacob,cartd)
     ! input : 2(ndime),4(nnode),deriv,ecord
     ! output : djacb,cartd

  ! compute deformation gradient tensor
  call getftens(2,4,edisploc,cartd, ftens2d)
    ! input : 2(ndime),4,edisploc,cartd
    ! output : ftens2d


  ! ----------------------------------------------------------------
  ! update constitutive law
  ! -----------------------
  ! convert cauchy stress tensor to kirchoff stress tensor
  call caucy2kir(2,ftens2d,sigtensloc2d, kirtensloc2d)
     ! input : 2(nndex),ftens2d,sigtensloc2d
     ! output : kirtensloc2d


  ! update elasto-viscoplastic kirchoff stress 
  call updevp1(ematpro,2,midltens2d,delt, effpstrn,kirtensloc2d, ctantens,effstrs) 
     ! input : ematpro,2(nndex),midltens2d,delt
     ! inoutput : effpstrn,kirtens
     ! output : ctantens,effstrs


  ! convert kirchoff stress tensor to cauchy stress tensor
  call kir2caucy(2,ftens2d,kirtensloc2d, sigtensloc2d)
     ! input : 2(nndex),ftens2d,kirtens
     ! output: sigtensloc2d


  ! ----------------------------------------------------------------
  ! set results
  ! -----------
  sigvoitloc(1,1)= sigtensloc2d(1,1) ! sig_x
  sigvoitloc(2,1)= sigtensloc2d(2,2) ! sig_y
  sigvoitloc(6,1)= sigtensloc2d(1,2) ! sig_xy



  return
end subroutine updsigevpbt1