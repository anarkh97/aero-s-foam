! ==================================
! internal force: bt shell / j2 plasticity
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
!  1.  subroutine         elefintj2bt1           (optcri,delt,ematpro,nndof,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo, &
!                                                 evar1,evoit2,evoit3, &
!                                                 evar2,efint)
!  2.  subroutine         elefintj2bt0           (optcri,delt,ematpro,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo, &
!                                                 evar1,evoit2,evoit3, &
!                                                 evar2,efint)
!  3.  subroutine         updj2strsbt            (delt,ematpro,gqpoin,gqweigt,locbvec,ecurn,evelo, &
!                                                 effpstrn,hardvar,sigvoitloc,strnvoitloc, &
!                                                 effstrs)
!  4.  subroutine         getsigj2bt1            (delt,ematpro,zeta,ecurnloc,eveloloc, effpstrn,hardvar,ipstrn,ipstrs, effstrs)
!
! =========================================================================================================



subroutine elefintj2bt1(optcri,delt,ematpro,nndof,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo, &
                        evar1,evoit2,evoit3, &
                        evar2,efint)
  !=======================================================================
  !  elefintj2bt1 = compute internal force matrix for bt shell
  !
  !                 note:
  !                 ----
  !                 rotation projection for the drilling dof is considered
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optcri(*) : fracture criterion option
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
  !  edisp(nndof,4) : element nodal displacement data
  !
  !  evelo(nndof,4) : element nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  inoutput:
  !  --------
  !  evar1(5,mgqpt1) : evar1(1,mgqpt1) effective plastic strain
  !                    evar1(2,mgqpt1) hardening variavle
  !
  !  evoit2(6,mgqpt1) : cauchy stress( local )
  !
  !  evoit3(6,mgqpt1) : strain( local )
  !
  !  output:
  !  ------
  !  evar2(1,:) : effective stress
  !
  !  efint(nndof*4,1) : element nodal internal force
  !
  ! ======================================================================

  use preset
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

  ! ------------------------------------

  real(8), dimension(5,mgqpt1), intent(inout) :: evar1
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

  else if ( nndof == 6) then ! 6 dof
     ! project 6 dof edisp(6,4) and evelo(6,4) to 5 dof edisp(5,4) and evelo(5,4)
     call rotprojbt1(ecord,edisp,evelo, edisp0,evelo0)
        ! input : ecord,edisp,evelo
        ! output : edisp0,evelo0

  end if


  ! ----------------------------------------------------------------------------------------------
  ! 5 dof bt shell
  ! --------------
  call elefintj2bt0(optcri,delt,ematpro,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp0,evelo0, &
                    evar1,evoit2,evoit3, &
                    evar2,efint0)
        ! input : optcri,delt,ematpro,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp0,evelo0
        ! inoutput : evar1(effpstrn, hardening),evoit2(cauchy, local),evoit3(strain)
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
end subroutine elefintj2bt1





subroutine elefintj2bt0(optcri,delt,ematpro,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo, &
                        evar1,evoit2,evoit3, &
                        evar2,efint)
  !=======================================================================
  !  elefintj2bt0 = compute element internal nodal force
  !                 of belytschko tsay shell element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optcri(*) : fracture criterion option
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
  !  evar1(5,mgqpt1) : evar1(1,mgqpt1) effective plastic strain
  !                    evar1(2,mgqpt1) hardening variable
  !
  !  evoit2(6,mgqpt1) : cauchy stress( local )
  !
  !  evoit3(6,mgqpt1) : strain( local )
  !
  !  output:
  !  ------
  !  evar2(1,:) : effective stress
  !
  !  efint(20,1) : element nodal internal force
  !
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, dimension(*), intent(in) :: optcri
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  integer, intent(in) :: mgaus3, mgqpt1
  real(8), dimension(mgaus3), intent(in) :: gqpoin3
  real(8), dimension(mgaus3), intent(in) :: gqweigt3
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(5,4), intent(in) :: edisp
  real(8), dimension(5,4), intent(in) :: evelo
  ! -------------------------------------

  real(8), dimension(5,mgqpt1), intent(inout) :: evar1
  real(8), dimension(6,mgqpt1), intent(inout) :: evoit2
  real(8), dimension(6,mgqpt1), intent(inout) :: evoit3

  ! -------------------------------------
  real(8), dimension(5,mgqpt1), intent(out) :: evar2
  real(8), dimension(20,1), intent(out) :: efint
    ! ====================================
  ! local variable
  ! ==============
  integer :: crityp

  real(8), dimension(3,3) :: locbvec
  real(8), dimension(3,4) :: ecurn
  real(8) :: effpstrn, hardvar, effstrs
  real(8), dimension(6,1) :: sigvoitloc, strnvoitloc

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

     ! hardening variable
     hardvar= evar1(2,igaus)

     ! cauchy stress(local)
     sigvoitloc(1:6,1)= evoit2(1:6,igaus)

     ! strain(local)
     strnvoitloc(1:6,1)= evoit3(1:6,igaus)
     ! -------------------------------------

     ! update hypo stresses at gq
     ! --------------------------
     call updj2strsbt(delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),locbvec,ecurn,evelo, &
                      effpstrn,hardvar,sigvoitloc,strnvoitloc, &
                      effstrs)
        ! input : delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),locbvec,ecurn,evelo
        ! inoutput : effpstrn,hardvar,sigvoitloc,strnvoitloc
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
           write(*,*) "not implemented fracture criterion: elefintj2bt0"
           write(nout5,*) "not implemented fracture criterion: elefintj2bt0"
           stop

        end select

     end if


     ! -------------------------------------
     ! set history variable
     ! --------------------
     ! effective plastic strain
     evar1(1,igaus)= effpstrn

     ! hardening variable
     evar1(2,igaus)= hardvar

     ! cauchy stress( local )
     evoit2(1:6,igaus)= sigvoitloc(1:6,1)

     ! strain(local)
     evoit3(1:6,igaus)= strnvoitloc(1:6,1)

     ! --------------------
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
end subroutine elefintj2bt0





subroutine updj2strsbt(delt,ematpro,gqpoin,gqweigt,locbvec,ecurn,evelo, &
                       effpstrn,hardvar,sigvoitloc,strnvoitloc, &
                       effstrs)
  !=======================================================================
  !  updj2strsbt = update j2 cauchy stress of belytschko tsay element 
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
  !  ecurn(3,4) : global current nodal coordinate
  !
  !  evelo(5,4) : global nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  inoutput:
  !  --------
  !  effpstrn : effective plastic strain
  !
  !  hardvar : hardening variable
  !
  !  sigvoitloc(6,1) : co-rotational cauchy stress
  !
  !  strnvoitloc(6,1) : co-rotational strain
  !
  !  output:
  !  ------
  !  effstrs : effective stress
  !
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: gqpoin, gqweigt
  real(8), dimension(3,3), intent(in) :: locbvec
  real(8), dimension(3,4), intent(in) :: ecurn
  real(8), dimension(5,4), intent(in) :: evelo

  ! ------------------------------------
  real(8), intent(inout) :: effpstrn,hardvar
  real(8), dimension(6,1), intent(inout) :: sigvoitloc, strnvoitloc
  ! ------------------------------------

  real(8), intent(out) :: effstrs
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,4) :: ecurnloc
  real(8), dimension(5,4) :: eveloloc
  real(8) :: zeta

  real(8), dimension(3,1) :: ipstrnloc, ipstrsloc
  real(8), dimension(2,1) :: tsstrslocdot
  ! ====================================

  ! initialize
  effstrs= 0.0d0

  ! -------------------------------------------------------------------------------
  ! convert global vector to local vector
  ! -------------------------------------
  ! get local nodal coordinate
  call glb2locnodv(3,4,3,locbvec,ecurn, ecurnloc)
     ! input : 3(ndime),4(nnode),3(ntrndof),locbvec,ecurn
     ! output : ecurnloc

  ! get local nodal velocity
  call glb2locnodv(3,4,5,locbvec,evelo, eveloloc)
     ! input : 3(ndime),4(nnode),5(ntrndof),locbvec,evelo
     ! output : eveloloc

  ! -------------------------------------------------------------------------------
  ! compute in plane stress rate, strain and update
  ! -----------------------------------------------
  ! pseudo-thickness parameter
  zeta= gqpoin

  ! extract in plane strain components
  ipstrnloc(1,1)= strnvoitloc(1,1)
  ipstrnloc(2,1)= strnvoitloc(2,1)
  ipstrnloc(3,1)= strnvoitloc(6,1)

  ! extract in plane stress components
  ipstrsloc(1,1)= sigvoitloc(1,1)
  ipstrsloc(2,1)= sigvoitloc(2,1)
  ipstrsloc(3,1)= sigvoitloc(6,1)


  ! update j2 in-plane stress of belytschko tsay shell element
  call getsigj2bt1(delt,ematpro,zeta,ecurnloc,eveloloc, &
                   effpstrn,hardvar,ipstrnloc,ipstrsloc, &
                   effstrs)
     ! input : delt,ematpro,zeta,ecurnloc,eveloloc
     ! inoutput : effpstrn,hardvar,ipstrnloc,ipstrsloc
     ! output : effstrs

  ! set updated in plane strain components
  strnvoitloc(1,1)= ipstrnloc(1,1) ! eps_x
  strnvoitloc(2,1)= ipstrnloc(2,1) ! eps_y
  strnvoitloc(6,1)= ipstrnloc(3,1) ! eps_xy

  ! set updated in plane stress components
  sigvoitloc(1,1)= ipstrsloc(1,1) ! sig_x
  sigvoitloc(2,1)= ipstrsloc(2,1) ! sig_y
  sigvoitloc(6,1)= ipstrsloc(3,1) ! sig_xy


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
end subroutine updj2strsbt






subroutine getsigj2bt1(delt,ematpro,zeta,ecurnloc,eveloloc, effpstrn,hardvar,ipstrn,ipstrs, effstrs)
  !=======================================================================
  !  getsigj2bt1 = 
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  delt : integration time step
  !
  !  ematpro(*) : material property
  !
  !  zeta : [-1, +1] : pseudo-thickness parameter
  !
  !  ecurnloc(3,4) : local current nodal coordinate
  !
  !  eveloloc(5,4) : local nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  inoutput:
  !  --------
  !  effpstrn : effective plastic strain
  !
  !  hardvar : hardening variable
  !
  !  ipstrn(3,1) : in plane strain
  !                ipstrn(3,1)= strn_x
  !                ipstrn(3,2)= strn_y
  !                ipstrn(3,3)= strn_xy
  !
  !  ipstrs(3,1) : in plane stress
  !                ipstrs(3,1)= sig_x
  !                ipstrs(3,2)= sig_y
  !                ipstrs(3,3)= sig_xy
  !
  !  output:
  !  ------
  !  effstrs : effective stress
  !
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: zeta
  real(8), dimension(3,4), intent(in) :: ecurnloc
  real(8), dimension(5,4), intent(in) :: eveloloc

  ! ------------------------------------
  real(8), intent(inout) :: effpstrn,hardvar
  real(8), dimension(3,1), intent(inout) :: ipstrn, ipstrs
  ! ------------------------------------

  real(8), intent(out) :: effstrs
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: thick
  real(8), dimension(2,2) :: ltens2d
  real(8), dimension(3,1) :: ipstrndot, ipstrndel

  ! ====================================

  ! initialize
  effstrs= 0.0d0

  ! ------------------------------------------------
  ! get material properties
  thick= ematpro(20)
  ! ------------------------------------------------


  ! ----------------------------------------------------------------
  ! compute rate of deformation
  ! ---------------------------
  ! ltens2d= [ vx,y  vx,y ]     ipstrndot= [d_x, d_y, 2d_xy]
  !          [ vy,x  vy,y ]
  call getstrndotbt1(thick,zeta,ecurnloc,eveloloc, ltens2d,ipstrndot)
     ! input : thick,zeta,ecurnloc,eveloloc
     ! output : ltens2d,ipstrndot

  ! strain increment
  ! ----------------
  ipstrndel(1:3,1)= delt * ipstrndot(1:3,1)


  ! calculate j2 plasticity
  call updj2pexp(ematpro,ipstrndel, effpstrn,hardvar,ipstrn,ipstrs, effstrs) 
     ! input : ematpro,ipstrndel
     ! inoutput : effpstrn,hardvar,ipstrn,ipstrs
     ! output : effstrs



  return
end subroutine getsigj2bt1