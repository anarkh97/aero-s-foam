! ==================================
! internal force: bt shell / j2 plasticity explicit
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
!  1.  subroutine         elefintj2bt1           (optcri,delt,ematpro,nndof,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo, &
!                                                 evar1,evoit2,evoit3, &
!                                                 evar2,efint)
!  2.  subroutine         elefintj2bt0           (optcri,delt,ematpro,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo, &
!                                                 evar1,evoit2,evoit3, &
!                                                 evar2,efint)
!  3.  subroutine         updstrsbt              (delt,ematpro,gqpoin,gqweigt,locbvec,ecurn,evelo, &
!                                                 effpstrn,hardvar,sigvoitloc,strnvoitloc, &
!                                                 effstrs)
!  4.  subroutine         getsigj2bt1            (delt,ematpro,zeta,ecurnloc,eveloloc, effpstrn,hardvar,ipstrn,ipstrs, effstrs)
!
! =========================================================================================================



subroutine elefintbt1(optctv,optdmg,optcri,opttrc,optcor,prmhgc, &
                      delt,ematpro,nndof,mgaus3,mgqpt1,gqpoin3,gqweigt3, &
                      ecord,edisp,evelo,trac,tmftval, &
                      evar1,evoit1,evoit2,evoit3, &
                      evar2,efint)
  !=======================================================================
  !  elefintj2bt1 = compute internal force matrix for bt shell, including hourglass force
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

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optctv, optdmg
  integer, dimension(*), intent(in) :: optcri
  integer, intent(in) :: opttrc
  integer, dimension(2), intent(in) :: optcor
  real(8), dimension(*), intent(in) :: prmhgc
  real(8), intent(in) :: delt
  real(8), dimension(20), intent(in) :: ematpro
  integer, intent(in) :: nndof,mgaus3,mgqpt1
  real(8), dimension(mgaus3), intent(in) :: gqpoin3
  real(8), dimension(mgaus3), intent(in) :: gqweigt3
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(6,4), intent(in) :: edisp
  real(8), dimension(6,4), intent(in) :: evelo
  real(8), dimension(3,1), intent(in) :: trac
  real(8), intent(in) :: tmftval

  ! ------------------------------------

  real(8), dimension(5,mgqpt1), intent(inout) :: evar1
  real(8), dimension(6,1), intent(inout) :: evoit1
  real(8), dimension(6,mgqpt1), intent(inout) :: evoit2
  real(8), dimension(6,mgqpt1), intent(inout) :: evoit3

  ! ------------------------------------
 
  real(8), dimension(5,mgqpt1), intent(out) :: evar2
  ! ------------------------------------
  real(8), dimension(24,1), intent(out) :: efint
  ! ====================================
  ! local variable
  ! ==============
  integer :: crityp

  real(8), dimension(3,3) :: locbvec
  real(8), dimension(3,4) :: ecurn
  real(8) :: effpstrn, hardvar, effstrs

  real(8), dimension(3,1) :: sigvoit2d, strnvoit2d
  real(8) :: crival, criang
  real(8) :: etatri

  real(8), dimension(24,1) :: efintloc

  ! loop index
  integer :: igaus

  !pjsa
  real(8), dimension(3,4) :: ecurnloc
  real(8), dimension(5,4) :: eveloloc
  real(8), dimension(2,4) :: bmat1pt
  real(8), dimension(2,4) :: bcmat1pt
  real(8), dimension(2,3,4) :: bsmat1pt
  real(8) :: area
  real(8), dimension(4,1) :: gamma
  real(8) :: zgamma
  real(8), dimension(3,1) :: ipstrndot
  real(8), dimension(2,1) :: tsstrndot
  ! ====================================

  ! initialize
  efintloc(:,:)= 0.0d0

  ! -------------------------------------------
  ! get fracture criterion parameters
  crityp= abs(optcri(1)) ! criterion type

  ! get current nodal coordinate
  ecurn(1:3,1:4)= ecord(1:3,1:4) + edisp(1:3,1:4)

  ! compute co rotational local base vector: locbvec
  ! ---------------------------------------
  call getlocbvecbt(ecurn, locbvec)
     ! input : ecurn
     ! output : locbvec

  !pjsa get local nodal coordinates and velocity
  call dgemm('t','n',3,4,3,1.0d0,locbvec,3,ecurn,3,0.0d0,ecurnloc,3)
  call dgemm('t','n',3,4,3,1.0d0,locbvec,3,evelo(1,1),6,0.0d0,eveloloc(1,1),5)
  call dgemm('t','n',2,4,3,1.0d0,locbvec,3,evelo(4,1),6,0.0d0,eveloloc(4,1),5)

  ! compute area
  area= 0.50d0*( (ecurnloc(1,3)-ecurnloc(1,1))*(ecurnloc(2,4)-ecurnloc(2,2)) &
                +(ecurnloc(1,2)-ecurnloc(1,4))*(ecurnloc(2,3)-ecurnloc(2,1)) )

  ! check current element configuration
  if ( area <= 0.0d0 ) then
     write(*,*) "current element has negative or zero area: elefintj2bt1"
     write(nout5,*) "current element has negative or zero area: elefintj2bt1"
     stop

  end if

  !compute b matrix: b matrix, b^c matrix, and b^s matrix
  ! ----------------
  ! compute b matrix: one point integration
  call getbmat1pt(ecurnloc,area, bmat1pt)
    ! input : ecurnloc,area
    ! output : bmat1pt

  ! compute b^c matrix: warping correction
  call getgamma4nod(ecurnloc,area, gamma,zgamma)
  if(optcor(1) > 0) then
    call getbcmat1pt(ecurnloc,area,gamma,zgamma, bcmat1pt)
  end if
    ! input : ecurnloc,area,gamma,zgamma
    ! output : bcmat1pt

  ! compute b^s matrix: transverse shear projection
  if(optcor(2) > 0) then
    call getbsmat1pt(ecurnloc, bsmat1pt)
  end if
    ! input : ecurnloc
    ! output : bsmat1pt

  ! loop over gauss quadrature
  ! note: this loop is for through thickness integration
  !       in plane, we use 1 point rule
  do igaus=1, mgaus3 

     ! compute rates of deformation and update strain at gq
     ! ----------------------------------------------------
     call updstrnbt(optcor,delt,ematpro,gqpoin3(igaus),eveloloc,bmat1pt,bcmat1pt,bsmat1pt, &
                    evoit3(1,igaus), &
                    ipstrndot,tsstrndot)

     ! update hypo stresses at gq
     ! --------------------------
     call updstrsbt(optctv,optdmg,delt,ematpro,area,ipstrndot,tsstrndot, &
                    evar1(1,igaus),evar1(2,igaus),evoit2(1,igaus),evoit3(1,igaus), &
                    evar2(1,igaus))
        ! input : delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),locbvec,ecurn,evelo
        ! inoutput : effpstrn,hardvar,sigvoitloc,strnvoitloc
        ! output : effstrs

     ! check fracture criterion
     ! ------------------------
     if ( crityp /= 0 ) then
        ! set 2d in-plane stress
        sigvoit2d(1,1)= evoit2(1,igaus) ! sig_xx
        sigvoit2d(2,1)= evoit2(2,igaus) ! sig_yy
        sigvoit2d(3,1)= evoit2(6,igaus) ! sig_xy

        select case(crityp)
        case(2) ! critical effective plastic strain criterion
           call chkeffstrcri2d(0,evar1(1,igaus),sigvoit2d, evar2(2,igaus),evar2(3,igaus))
              ! input : 0(optstr:stress),effpstrn,sigvoit2d
              ! output : crival,criang

        case(5) ! xue-wierzbicki with mtps criterion for j2 plane stress
           call chkxwj2pstrscri2d(0,evar1(1,igaus),sigvoit2d, evar2(2,igaus),evar2(3,igaus),evar2(4,igaus))
              ! input : 0(optstr:stress),effpstrn,sigvoit2d
              ! output : crival,criang,etatri

        case default
           write(*,*) "not implemented fracture criterion: elefintj2bt0"
           write(nout5,*) "not implemented fracture criterion: elefintj2bt0"
           stop

        end select

     end if

     ! add local nodal internal forces at gq
     ! -----------------------------------------
     call gqfintbt(optcor,delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),area,evoit2(1,igaus), &
                   bmat1pt, bcmat1pt, bsmat1pt, efintloc)
        ! input : delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),locbvec,ecurn,sigvoitloc
        ! in/output : efintloc

  end do

  ! -------------------------------------
  ! update hourglass control stresses
  ! ---------------------------------
  call updhgcstrsbt(prmhgc,delt,ematpro,eveloloc,area,bmat1pt, &
                    gamma,zgamma, evoit1)
     ! input : prmhgc,delt,ematpro,eveloloc,area,bmat1pt,gamma,zgamma
     ! inoutput : hgcvoitloc

  ! add the hourglass control forces
  ! --------------------------------------
  call gqfhgcbt(gamma,zgamma,evoit1, efintloc)
     ! input : gamma,zgamma,evoit1
     ! inoutput : efintloc

  ! --------------------------------------------------------------
  ! subtract the local traction forces
  ! -------------------------------------
  if (opttrc >= 0) then
    call elefbc3dbrkshl2opt(area,trac,tmftval, efintloc)
       ! input : area,trac,tmftval
       ! inoutput : efintloc
  end if

  ! --------------------------------------------------------------
  ! convert local efintloc to global efint
  ! -------------------------------------
  call dgemm('n','n',3,8,3,1.0d0,locbvec,3,efintloc,3,0.0d0,efint,3)


  return
end subroutine elefintbt1




subroutine updstrsbt(optctv,optdmg,delt,ematpro,area,ipstrndot,tsstrndot, &
                     effpstrn,hardvar,sigvoitloc,strnvoitloc, &
                     effstrs)
  !=======================================================================
  !  updstrsbt = update j2 (or hypoelas) cauchy stress of belytschko tsay element 
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
  !  inoutput:
  !  --------
  !  effpstrn : effective plastic strain (j2) or effective strain (hypoelas)
  !
  !  hardvar : hardening variable (j2) or damage (hypoelas)
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

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optctv, optdmg
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: area
  real(8), dimension(3,1), intent(in) :: ipstrndot
  real(8), dimension(2,1), intent(in) :: tsstrndot

  ! ------------------------------------
  real(8), intent(inout) :: effpstrn,hardvar
  real(8), dimension(6,1), intent(inout) :: sigvoitloc, strnvoitloc
  ! ------------------------------------

  real(8), intent(out) :: effstrs
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,1) :: ipstrnloc, ipstrsloc
  real(8), dimension(2,1) :: tsstrslocdot
  real(8) :: ehleng
  ! ====================================

  ! -------------------------------------------------------------------------------
  ! update in plane stress 
  ! -----------------------------------------------

  ! extract in plane strain components
  ipstrnloc(1,1)= strnvoitloc(1,1)
  ipstrnloc(2,1)= strnvoitloc(2,1)
  ipstrnloc(3,1)= strnvoitloc(6,1)

  ! extract in plane stress components
  ipstrsloc(1,1)= sigvoitloc(1,1)
  ipstrsloc(2,1)= sigvoitloc(2,1)
  ipstrsloc(3,1)= sigvoitloc(6,1)

  select case(optctv)
    case(1)
    ! update hypoelastic in-plane stress of belytschko tsay shell element
    ehleng= dsqrt(area)
    call getsighypobt1(optdmg,delt,ematpro,ehleng,ipstrndot,ipstrnloc, &
                       hardvar,ipstrsloc)
     ! input : optdmg,delt,ematpro,ehleng,ipstrnloc
     ! inoutput : damage,ipstrsloc

    case(5)
    ! update j2 in-plane stress of belytschko tsay shell element
    call updj2pexp(1,2,ematpro,delt,ipstrndot, &
                   effpstrn,hardvar,ipstrsloc, &
                   effstrs)
       ! input : 1(optpty:p-strs),2(ndime),ematpro,delt,ipstrndot
       ! inoutput : effpstrn,hardvar,ipstrs
       ! output : effstrs

    case default
       write(*,*) "not implemented constitutive model: updstrsbt"
       write(nout5,*) "not implemented constitutive model: updstrsbt"
       stop

  end select

  ! set updated in plane stress components
  sigvoitloc(1,1)= ipstrsloc(1,1) ! sig_x
  sigvoitloc(2,1)= ipstrsloc(2,1) ! sig_y
  sigvoitloc(6,1)= ipstrsloc(3,1) ! sig_xy

  ! -------------------------------------------------------------------------------
  ! compute tranverse shear stress rate and update
  ! ----------------------------------------------
  ! compute hypoelastic transverse shear stress rate of belytschko tsay shell element
  call getsighypobt2(ematpro,tsstrndot, tsstrslocdot)
     ! input : ematpro,tsstrndot
     ! output : tsstrslocdot

  ! set updated transverse stress components
  sigvoitloc(3,1)= 0.0d0 ! sig_z (plane stress)
  sigvoitloc(4,1)= sigvoitloc(4,1) + tsstrslocdot(2,1) * delt ! sig_yz
  sigvoitloc(5,1)= sigvoitloc(5,1) + tsstrslocdot(1,1) * delt ! sig_xz
  ! -------------------------------------------------------------------------------
  
  return
end subroutine updstrsbt


subroutine updstrnbt(optcor,delt,ematpro,gqpoin,eveloloc,bmat1pt,bcmat1pt,bsmat1pt, &
                     strnvoitloc, &
                     ipstrndot,tsstrndot) 
  !=======================================================================
  !  updstrsbt = update strain of belytschko tsay element 
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
  !  eveloloc(5,4) : local nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  inoutput:
  !  --------
  !  strnvoitloc(6,1) : co-rotational strain
  !
  !  output:
  !  ------
  !  ipstrndot(3,1), tsstrndot(2,1) : in plane and transverse shear rate of deformation
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, dimension(2), intent(in) :: optcor
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: gqpoin
  real(8), dimension(6,4), intent(in) :: eveloloc
  real(8), dimension(2,4), intent(in) :: bmat1pt
  real(8), dimension(2,4), intent(in) :: bcmat1pt
  real(8), dimension(2,3,4), intent(in) :: bsmat1pt

  ! ------------------------------------
  real(8), dimension(6,1), intent(inout) :: strnvoitloc
  ! ------------------------------------

  real(8), dimension(3,1), intent(out) :: ipstrndot
  real(8), dimension(2,1), intent(out) :: tsstrndot
  ! ====================================

  ! ----------------------------------------------------------------
  ! compute in plane rate of deformation
  ! ---------------------------
  ! ipstrndot= [d_x, d_y, 2d_xy]
  call getstrndotbt1(optcor,ematpro(20),gqpoin,eveloloc,bmat1pt,bcmat1pt, ipstrndot)
     ! input : optcor,thick,zeta,eveloloc,bmat1pt,bcmat1pt
     ! output : ipstrndot

  ! -------------------------------------------------------------------------------
  ! update in plane strain
  ! -----------------------------------------------
  strnvoitloc(1,1)= strnvoitloc(1,1) + delt*ipstrndot(1,1) ! eps_x
  strnvoitloc(2,1)= strnvoitloc(2,1) + delt*ipstrndot(2,1) ! eps_y
  strnvoitloc(6,1)= strnvoitloc(6,1) + delt*ipstrndot(3,1) ! eps_xy

  ! ----------------------------------------------------------------
  ! compute transverse shear rate of deformation
  ! ---------------------------
  ! tsstrndot= [2d_xz, 2dyz]
  call getstrndotbt2(optcor,eveloloc,bmat1pt,bsmat1pt, tsstrndot)
     ! input : optcor,eveloloc,bmat1pt,bsmat1pt
     ! output : tsstrndot

  ! -------------------------------------------------------------------------------
  ! update transverse shear strain
  ! -----------------------------------------------
  ! TODO (currently this isn't used anywhere, except possibly for postprocessing
  
  return
end subroutine updstrnbt
