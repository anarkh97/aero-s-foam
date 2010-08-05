! ==================================
! internal force: bt shell / hypo-elastic
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine        elefinthypobt1         (optdmg,optcri, &
!                                               delt,ematpro,nndof,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo, &
!                                               evar1,evoit2,evoit3, &
!                                               evar2,efint)
! 2.  subroutine        updhypostrsbt          (optdmg, &
!                                               delt,ematpro,gqpoin,gqweigt,locbvec,ecurn,evelo, &
!                                               damage,sigvoitloc,strnvoitloc, &
!                                               effstrn,effstrs)
! 3.  subroutine        gqfintbt               (delt,ematpro,gqpoin,gqweigt,locbvec,ecurn,sigvoitloc, gqfint)
! 4.  subroutine        getsighypobt1          (optdmg,delt,ematpro,ehleng,zeta,ecurnloc,eveloloc, damage,ipstrn, ipstrsdot)
! 5.  subroutine        getsighypobt2          (ematpro,ecurnloc,eveloloc, tsstrsdot)
! 6.  subroutine        getstrndotbt1          (thick,zeta,ecurnloc,eveloloc, ltens2d,ipstrndot)
! 7.  subroutine        getstrndotbt2          (ecurnloc,eveloloc, tsstrndot)
!
! =========================================================================================================



subroutine elefinthypobt1(optdmg,optcri, &
                          delt,ematpro,nndof,mgaus3,mgqpt1,gqpoin3,gqweigt3,ecord,edisp,evelo, &
                          evar1,evoit2,evoit3, &
                          evar2,efint)
  !=======================================================================
  !  elefinthypobt1 = compute internal force matrix for bt shell
  !                   note:
  !                   ----
  !                   rotation projection for the drilling dof is considered
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optdmg : damage model option handler
  !
  !  optcri(*) : fracture criterion
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
  !  evar1(5,mgqpt1) : effective stress and damage
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
  integer, intent(in) :: optdmg
  integer, dimension(*), intent(in) :: optcri
  real(8), intent(in) :: delt
  real(8), dimension(20), intent(in) :: ematpro
  integer, intent(in) :: nndof,mgaus3,mgqpt1
  real(8), dimension(mgaus3), intent(in) :: gqpoin3
  real(8), dimension(mgaus3), intent(in) :: gqweigt3
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(6,4), intent(in) :: edisp
  real(8), dimension(6,4), intent(in) :: evelo

  ! ------------------------------------

  real(8), dimension(5,mgqpt1), intent(inout) :: evar1
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
  real(8) :: damage, effstrs, effstrn
  real(8), dimension(6,1) :: sigvoitloc, strnvoitloc

  real(8), dimension(3,1) :: sigvoit2d, strnvoit2d
  real(8) :: crival, criang

  real(8), dimension(24,1) :: gqfint

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
     ! damage
     damage= evar1(2,igaus)

     ! cauchy stress(local)
     sigvoitloc(1:6,1)= evoit2(1:6,igaus)

     ! strain(local)
     strnvoitloc(1:6,1)= evoit3(1:6,igaus)
     ! -------------------------------------

     ! update hypo stresses at gq
     ! --------------------------
     call updhypostrsbt(optdmg, &
                        delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),locbvec,ecurn,evelo, &
                        damage,sigvoitloc,strnvoitloc, &
                        effstrn,effstrs)
        ! input : optdmg
        !         delt,ematpro,gqpoin3(igaus),gqweigt3(igaus),locbvec,ecurn,evelo
        ! inoutput : damage,sigvoitloc,strnvoitloc
        ! output : effstrn,effstrs

     ! check fracture criterion
     ! ------------------------
     if ( crityp /= 0 ) then

        select case(crityp)
        case(2) ! critical effective plastic strain criterion
           ! set 2d in-plane stress
           sigvoit2d(1,1)= sigvoitloc(1,1) ! sig_xx
           sigvoit2d(2,1)= sigvoitloc(2,1) ! sig_yy
           sigvoit2d(3,1)= sigvoitloc(6,1) ! sig_xy

           call chkeffstrcri2d(0,effstrn,sigvoit2d, crival,criang)
              ! input : 0(opt:stress),effstrn,sigvoit2d
              ! output : crival,criang

        case default
           write(*,*) "not implemented fracture criterion: elefinthypobt0"
           write(nout5,*) "not implemented fracture criterion: elefinthypobt0"
           stop

        end select

     end if
     
     ! -------------------------------------
     ! set history variable
     ! --------------------
     ! damage
     evar1(2,igaus)= damage

     ! cauchy stress( local )
     evoit2(1:6,igaus)= sigvoitloc(1:6,1)

     ! strain(local)
     evoit3(1:6,igaus)= strnvoitloc(1:6,1)

     ! --------------------
     ! effective strain
     evar1(1,igaus)= effstrn

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
     efint(1:24,1)= efint(1:24,1) + gqfint(1:24,1)

  end do


  return
end subroutine elefinthypobt1








subroutine updhypostrsbt(optdmg, &
                         delt,ematpro,gqpoin,gqweigt,locbvec,ecurn,evelo, &
                         damage,sigvoitloc,strnvoitloc, &
                         effstrn,effstrs)
  !=======================================================================
  !  updhypostrsbt = update hypoelastic cauchy stress of belytschko tsay element 
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optdmg : damage model option handler
  !
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
  !  damage : material damage parameter
  !
  !  sigvoitloc(6,1) : co-rotational cauchy stress
  !
  !  sigvoitloc(6,1) : co-rotational strain
  !
  !  output:
  !  ------
  !  effstrn : effective strain
  !
  !  effstrs : effective stress
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optdmg
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: gqpoin, gqweigt
  real(8), dimension(3,3), intent(in) :: locbvec
  real(8), dimension(3,4), intent(in) :: ecurn
  real(8), dimension(6,4), intent(in) :: evelo

  ! ------------------------------------
  real(8), intent(inout) :: damage
  real(8), dimension(6,1), intent(inout) :: sigvoitloc, strnvoitloc
  ! ------------------------------------

  real(8), intent(out) :: effstrn, effstrs
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: ehleng
  real(8), dimension(3,4) :: ecurnloc
  real(8), dimension(6,4) :: eveloloc0
  real(8), dimension(5,4) :: eveloloc
  real(8) :: zeta

  real(8), dimension(3,1) :: ipstrslocdot, ipstrnloc
  real(8), dimension(2,1) :: tsstrslocdot

  real(8), dimension(2,2) :: devstrs2d, devstrn2d
  real(8), dimension(2,2) :: strs2dloc, strn2dloc
  ! ====================================

  ! initialize
  effstrn= 0.0d0
  effstrs= 0.0d0


  ! -------------------------------------------------------------------------------
  ! compute characteristic h length
  ! -------------------------------
  call elehleng0(3,3,ecurn, ehleng)
     ! input : 3(optele),3(ndime),ecurn
     ! output : ehleng

  ! -------------------------------------------------------------------------------
  ! convert global vector to local vector
  ! -------------------------------------
  ! get local nodal coordinate
  call glb2locnodv(3,4,3,locbvec,ecurn, ecurnloc)
     ! input : 3(ndime),4(nnode),3(ntrndof),locbvec,ecurn
     ! output : ecurnloc

  ! get local nodal velocity
  call glb2locnodv(3,4,6,locbvec,evelo, eveloloc0)
     ! input : 3(ndime),4(nnode),6(ntrndof),locbvec,evelo
     ! output : eveloloc0

  eveloloc(1:5,1:4)= eveloloc0(1:5,1:4)

  ! -------------------------------------------------------------------------------
  ! compute in plane stress rate, strain and update
  ! -----------------------------------------------
  ! pseudo-thickness parameter
  zeta= gqpoin

  ! extract in plane strain components
  ipstrnloc(1,1)= strnvoitloc(1,1)
  ipstrnloc(2,1)= strnvoitloc(2,1)
  ipstrnloc(3,1)= strnvoitloc(6,1)

  ! compute hypoelastic in-plane stress rate of belytschko tsay shell element
  call getsighypobt1(optdmg,delt,ematpro,ehleng,zeta,ecurnloc,eveloloc, &
                     damage,ipstrnloc, &
                     ipstrslocdot)
     ! input : optdmg,delt,ematpro,ehleng,zeta,ecurnloc,eveloloc
     ! inoutput : damage,ipstrnloc
     ! output : ipstrslocdot

  ! set updated in plane strain components
  strnvoitloc(1,1)= ipstrnloc(1,1)
  strnvoitloc(2,1)= ipstrnloc(2,1)
  strnvoitloc(6,1)= ipstrnloc(3,1)

  ! update local cauchy stress
  sigvoitloc(1,1)= sigvoitloc(1,1) + ipstrslocdot(1,1) * delt ! sig_x
  sigvoitloc(2,1)= sigvoitloc(2,1) + ipstrslocdot(2,1) * delt ! sig_y
  sigvoitloc(6,1)= sigvoitloc(6,1) + ipstrslocdot(3,1) * delt ! sig_xy


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
  ! compute effective strain
  ! ------------------------
  strn2dloc(1,1)= strnvoitloc(1,1)
  strn2dloc(1,2)= 0.50d0 * strnvoitloc(6,1) 
  strn2dloc(2,1)= 0.50d0 * strnvoitloc(6,1) 
  strn2dloc(2,2)= strnvoitloc(2,1) 

  ! compute deviatoric stress tensor
  call getdevtens(2,strn2dloc, devstrn2d)
     ! input : 2(nndex),strn2dloc
     ! output : devstrn2d

  ! compute effective strain  
  call geteffstrn(2,devstrn2d, effstrn)
     ! input : 2(nndex),devstrn2d
     ! output : effstrn

  ! -------------------------------------------------------------------------------
  ! compute effective stress
  ! ------------------------
  strs2dloc(1,1)= sigvoitloc(1,1)
  strs2dloc(1,2)= sigvoitloc(6,1)
  strs2dloc(2,1)= sigvoitloc(6,1)
  strs2dloc(2,2)= sigvoitloc(2,1)

  ! compute deviatoric stress tensor
  call getdevtens(2,strs2dloc, devstrs2d)
     ! input : 2(nndex),strs2dloc
     ! output : devstrs2d

  ! compute effective stress  
  call geteffstrs(2,devstrs2d, effstrs)
     ! input : 2(nndex),devstrs2d
     ! output : effstrs

  ! -------------------------------------------------------------------------------


  return
end subroutine updhypostrsbt





subroutine gqfintbt(delt,ematpro,gqpoin,gqweigt,locbvec,ecurn,sigvoitloc, gqfint)
  !=======================================================================
  !  gqfintbt = compute global internal force of belytschko tsay shell element
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
  !  locbvec(3,3) : local base vector
  !
  !  ecurn(3,4) : global current nodal coordinate
  !
  !  sigvoitloc(6,1) : voight form of local cauchy stress
  !
  !  output:
  !  ------
  !  gqfintloc(20,1) : element local nodal force
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
  real(8), dimension(3,4), intent(in) :: ecurn
  real(8), dimension(6,1), intent(in) :: sigvoitloc

  real(8), dimension(24,1), intent(out) :: gqfint
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: rk, thick

  real(8), dimension(3,4) :: ecurnloc
  real(8), dimension(3,1) :: ipstrsloc
  real(8), dimension(2,1) :: tsstrsloc

  real(8) :: psibar, dtdpsi

  real(8) :: fx, fy, fxy, fxz, fyz
  real(8) :: mx, my, mxy

  real(8), dimension(2,4) :: bmat1pt
  real(8), dimension(2,4) :: bcmat1pt
  real(8), dimension(2,3,4) :: bsmat1pt
  real(8) :: area

  real(8), dimension(24,1) :: gqfintloc
  real(8), dimension(6,1) :: locvec, glbvec

  integer :: iloc

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  gqfintloc(:,:)= 0.0d0


  ! ------------------------------------------------
  ! material properties
  rk= ematpro(19)     ! shear correction factor: 0.840d0 
  thick= ematpro(20)

  ! ------------------------------------------------

  ! get local nodal coordinate
  call glb2locnodv(3,4,3,locbvec,ecurn, ecurnloc)
     ! input : 3(ndime),4(nnode),3(ntrndof),locbvec,ecurn
     ! output : ecurnloc
  
  ! set rate form of cauchy local stress
  ipstrsloc(1,1)= sigvoitloc(1,1) ! sig_x
  ipstrsloc(2,1)= sigvoitloc(2,1) ! sig_y

  tsstrsloc(2,1)= sigvoitloc(4,1) ! sig_yz
  tsstrsloc(1,1)= sigvoitloc(5,1) ! sig_xz
  ipstrsloc(3,1)= sigvoitloc(6,1) ! sig_xy
  ! ------------------------------------------------

  ! compute b matrix: b matrix, b^c matrix, and b^s matrix
  ! ----------------
  ! compute b matrix: one point integration
  call getbmat1pt(ecurnloc, bmat1pt)
     ! input : ecurnloc
     ! output : bmat1pt

  ! compute b^c matrix: warpping correction
  call getbcmat1pt(ecurnloc, bcmat1pt)
     ! input : ecurnloc
     ! output : bcmat1pt

  ! compute b^s matrix: transverse shear projection
  call getbsmat1pt(ecurnloc, bsmat1pt)
     ! input : ecurnloc
     ! output : bsmat1pt

  ! compute area
  area= 0.50d0*( (ecurnloc(1,3)-ecurnloc(1,1))*(ecurnloc(2,4)-ecurnloc(2,2)) &
                +(ecurnloc(1,2)-ecurnloc(1,4))*(ecurnloc(2,3)-ecurnloc(2,1)) )


  ! -----------------------------------------------------------
  ! pseudo-thickness : eq (5)
  psibar= (gqpoin*thick)/2.0d0

  ! d thickness / d psi
  dtdpsi= thick / 2.0d0

  ! in plane stress
  fx=  ipstrsloc(1,1) *gqweigt *dtdpsi ! f_x
  fy=  ipstrsloc(2,1) *gqweigt *dtdpsi ! f_y
  fxy= ipstrsloc(3,1) *gqweigt *dtdpsi ! f_xy

  ! tranverse shear stress
  fxz= tsstrsloc(1,1) *gqweigt *dtdpsi ! f_xz
  fyz= tsstrsloc(2,1) *gqweigt *dtdpsi ! f_yz

  ! in plane moment
  mx=  psibar * fx ! m_x
  my=  psibar * fy ! m_y
  mxy= psibar * fxy ! m_xy

  ! -----------------------------------------------------------
  ! compute local internal nodal force
  ! ----------------------------------
  gqfintloc(:,:)= 0.0d0 ! initialize

  do inode=1, 4

     ! address
     iloc= 6*(inode-1)

     ! compute internal nodal force

! ### working
! --------------------------------------------------------------------------------------
! the original BLT element with warping and shear corrections
! -----------------------------------------------------------
!     gqfintloc(iloc+1,1)= area*( bmat1pt(1,inode)*fx + bmat1pt(2,inode)*fxy &
!                                + bcmat1pt(1,inode)*mx + bcmat1pt(2,inode)*mxy )  ! f_xi
!
!     gqfintloc(iloc+2,1)= area*( bmat1pt(2,inode)*fy + bmat1pt(1,inode)*fxy &
!                                + bcmat1pt(2,inode)*my + bcmat1pt(1,inode)*mxy )  ! f_yi
!
!     gqfintloc(iloc+3,1)= area*rk*( bsmat1pt(1,1,inode)*fxz + bsmat1pt(2,1,inode)*fyz ) ! f_zi
!
!     gqfintloc(iloc+4,1)= area*( rk*bsmat1pt(1,2,inode)*fxz + rk*bsmat1pt(2,2,inode)*fyz &
!                                - bmat1pt(2,inode)*my - bmat1pt(1,inode)*mxy)  ! m_xi
!
!     gqfintloc(iloc+5,1)= area*( rk*bsmat1pt(1,3,inode)*fxz + rk*bsmat1pt(2,3,inode)*fyz &
!                                + bmat1pt(1,inode)*mx + bmat1pt(2,inode)*mxy)  ! m_yi
!
!     gqfintloc(iloc+6,1)= 0.0d0  ! m_zi

! --------------------------------------------------------------------------------------
! the original BLT element with warping correction
! -----------------------------------------------------------
     gqfintloc(iloc+1,1)= area*( bmat1pt(1,inode)*fx + bmat1pt(2,inode)*fxy &
                                + bcmat1pt(1,inode)*mx + bcmat1pt(2,inode)*mxy )  ! f_xi

     gqfintloc(iloc+2,1)= area*( bmat1pt(2,inode)*fy + bmat1pt(1,inode)*fxy &
                                + bcmat1pt(2,inode)*my + bcmat1pt(1,inode)*mxy )  ! f_yi

     gqfintloc(iloc+3,1)= area*rk*( bmat1pt(1,inode)*fxz + bmat1pt(2,inode)*fyz ) ! f_zi

     gqfintloc(iloc+4,1)= area*( - bmat1pt(2,inode)*my - bmat1pt(1,inode)*mxy - 0.250d0*rk*fyz ) ! m_xi

     gqfintloc(iloc+5,1)= area*( bmat1pt(1,inode)*mx + bmat1pt(2,inode)*mxy + 0.250d0*rk*fxz )  ! m_yi

     gqfintloc(iloc+6,1)= 0.0d0  ! m_zi

! --------------------------------------------------------------------------------------
! the original BLT element
! ---------------------------
!     gqfintloc(iloc+1,1)= area*( bmat1pt(1,inode)*fx + bmat1pt(2,inode)*fxy )  ! f_xi
!
!     gqfintloc(iloc+2,1)= area*( bmat1pt(2,inode)*fy + bmat1pt(1,inode)*fxy )  ! f_yi
!
!     gqfintloc(iloc+3,1)= area*rk*( bmat1pt(1,inode)*fxz + bmat1pt(2,inode)*fyz ) ! f_zi
!
!     gqfintloc(iloc+4,1)= area*( - bmat1pt(2,inode)*my - bmat1pt(1,inode)*mxy - 0.250d0*rk*fyz ) ! m_xi
!
!     gqfintloc(iloc+5,1)= area*( bmat1pt(1,inode)*mx + bmat1pt(2,inode)*mxy + 0.250d0*rk*fxz )  ! m_yi
!
!     gqfintloc(iloc+6,1)= 0.0d0  ! m_zi

! --------------------------------------------------------------------------------------

  end do


  ! --------------------------------------------------------------
  ! convert local gqfint to global gqfint
  ! -------------------------------------
  do inode=1, 4

     ! address
     iloc= 6*(inode-1)

     ! get local force vector
     locvec(1:6,1)= gqfintloc(iloc+1:iloc+6,1)

     ! convert
     call loc2glbv(3,6,locbvec,locvec, glbvec)
        ! input : 3(ndime),6(ntrndof),locbvec,locvec
        ! output : glbvec
 
     ! set global force vector
     gqfint(iloc+1:iloc+6,1)= glbvec(1:6,1)

  end do
  ! --------------------------------------------------------------



  return
end subroutine gqfintbt





subroutine getsighypobt1(optdmg,delt,ematpro,ehleng,zeta,ecurnloc,eveloloc, &
                         damage,ipstrn, &
                         ipstrsdot)
  !=======================================================================
  !  getsighypobt1 = compute rate of belytschko-tsay local element cauchy stress 
  !
  !                  note:
  !                  ----
  !                  belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !                  advances in one point quadrature shell elements
  !                  (see, apendix a.3 )
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optdmg : damage model option handler
  !
  !  delt : integration time step
  !
  !  ematpro(*) : material property
  !
  !  ehleng : characteristic element length
  !
  !  zeta : [-1, +1] : pseudo-thickness parameter
  !
  !  ecurnloc(3,4) : local current nodal coordinate
  !
  !  eveloloc(5,4) : local nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  inoutput:
  !  --------
  !  damage : material damage parameter
  !
  !  ipstrn(3,1) : in plane strain
  !                ipstrn(3,1)= strn_x
  !                ipstrn(3,2)= strn_y
  !                ipstrn(3,3)= strn_xy
  !
  !  output:
  !  ------
  !  ipstrsdot(3,1) : in plane stress rate
  !                   ipstrsdot(1,1)= strs_x dot
  !                   ipstrsdot(2,1)= strs_y dot
  !                   ipstrsdot(3,1)= strs_xy dot
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optdmg
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: ehleng, zeta
  real(8), dimension(3,4), intent(in) :: ecurnloc
  real(8), dimension(5,4), intent(in) :: eveloloc

  ! ------------------------------------
  real(8), intent(inout) :: damage
  real(8), dimension(3,1), intent(inout) :: ipstrn
  ! ------------------------------------

  real(8), dimension(3,1), intent(out) :: ipstrsdot
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: young, poiss, thick

  real(8), dimension(2,2) :: ltens2d, etens2d
  real(8), dimension(3,1) :: ipstrndot

  real(8), dimension(3,3) :: elsvoit2d, csntvoit2d, ctanvoit2d
  real(8), dimension(2,2,2,2) :: ctantens2d
  ! ====================================

  ! initialize
  ipstrsdot(:,:)= 0.0d0


  ! ------------------------------------------------
  ! get material properties
  young= ematpro(1)
  poiss= ematpro(2)
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

  ! update strain
  ! -------------
  ipstrn(1:3,1)= ipstrn(1:3,1) + delt * ipstrndot(1:3,1)


  ! ----------------------------------------------------------------
  ! compute material damage
  ! -----------------------
  select case(optdmg)
  case(0)
     ! no damage
     damage= 0.0d0

     ! compute voight form elastic moduli matrix
     call getelsvoit2d(1,young,poiss, elsvoit2d)
       ! input : 1(optpty=plane-strs),young,poiss
       ! output : elsvoit2d

  case(1) ! lematire damage model
     etens2d(1,1)= ipstrn(1,1)
     etens2d(2,2)= ipstrn(2,1)
     etens2d(1,2)= 0.50d0*ipstrn(3,1)
     etens2d(2,1)= 0.50d0*ipstrn(3,1)

     ! compute tangent elastic moduli based on damage
     call getlemdmg2d0(1,ematpro,etens2d, damage, csntvoit2d,ctanvoit2d,ctantens2d)
        ! input : 1(optpty:p-strs),ematpro,etens2d
        ! inoutput : damage
        ! output : csntvoit2d,ctanvoit2d,ctantens2d

     ! compute voight form elastic moduli matrix
     elsvoit2d(1:3,1:3)= ctanvoit2d(1:3,1:3)

  case(2) ! linear softening with scaling
     etens2d(1,1)= ipstrn(1,1)
     etens2d(2,2)= ipstrn(2,1)
     etens2d(1,2)= 0.50d0*ipstrn(3,1)
     etens2d(2,1)= 0.50d0*ipstrn(3,1)

     ! compute tangent elastic moduli based on damage
     call getlineardmg2d0(1,ematpro,ehleng,etens2d, damage, csntvoit2d,ctanvoit2d,ctantens2d)
        ! input : 1(optpty:p-strs),ematpro,ehleng,etens2d
        ! inoutput : damage
        ! output : csntvoit2d,ctanvoit2d

     ! compute voight form elastic moduli matrix
     elsvoit2d(1:3,1:3)= ctanvoit2d(1:3,1:3)

  case default
     write(*,*) "not implemented damage model: gqfinthypo2d0"
     write(nout5,*) "not implemented damage model: gqfinthypo2d0"
     stop

  end select


  ! ----------------------------------------------------------------
  ! compute stress rate
  ! -------------------
  ! compute stress rate: plane stress law is used
  ! in-plane stress rate: [strs_x, strs_y. strs_xy]
  call matprd(3,3,0, 3,1,0, 3,1, elsvoit2d,ipstrndot, ipstrsdot)
     ! input : 3,3,0, 3,1,0, 3,1, elsvoit2d,ipstrndot
     ! output : ipstrsdot



  return
end subroutine getsighypobt1





subroutine getsighypobt2(ematpro,ecurnloc,eveloloc, tsstrsdot)
  !=======================================================================
  !  getsighypobt2 = compute rate of belytschko-tsay local element cauchy stress 
  !
  !                  note:
  !                  ----
  !                  belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !                  advances in one point quadrature shell elements
  !                  (see, apendix a.3 )
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : material property
  !
  !  ecurnloc(3,4) : local current nodal coordinate
  !
  !  eveloloc(5,4) : local nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  output:
  !  ------
  !  tsstrsdot(2,1) : transverse shear stress rate
  !                   tsstrsdot(1,1)= strs_xz dot
  !                   tsstrsdot(2,1)= strs_yz dot
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(3,4), intent(in) :: ecurnloc
  real(8), dimension(5,4), intent(in) :: eveloloc

  real(8), dimension(2,1), intent(out) :: tsstrsdot
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: young, poiss
  real(8) :: mu

  real(8), dimension(2,1) :: tsstrndot

  ! ====================================

  ! initialize
  tsstrsdot(:,:)= 0.0d0


  ! ------------------------------------------------
  ! get material properties
  young= ematpro(1)
  poiss= ematpro(2)

  ! ------------------------------------------------

  ! 2nd lame constant: shear modulus 
  mu= young / ( 2.0d0 * (1.0d0+poiss) )


  ! ----------------------------------------------------------------
  ! compute rate of deformation
  ! ---------------------------
  ! tsstrndot= [2d_xz, 2dyz]
  call getstrndotbt2(ecurnloc,eveloloc, tsstrndot)
     ! input : ecurnloc,eveloloc
     ! output : tsstrndot

  ! ----------------------------------------------------------------
  ! compute stress rate
  ! -------------------
  ! tranverse stress rate: [strs_xz, strs_yz]
  tsstrsdot(1,1)= mu * tsstrndot(1,1)
  tsstrsdot(2,1)= mu * tsstrndot(2,1)


  return
end subroutine getsighypobt2





subroutine getstrndotbt1(thick,zeta,ecurnloc,eveloloc, ltens2d,ipstrndot)
  !=======================================================================
  !  getstrndotbt1 = compute in plane local rate of deformation (velocity strain) of tb shell
  !                  with one point integration and local z method
  !                  for computational efficiency, explicit components form is used
  !
  !                 note:
  !                 ----
  !                 belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !                 advances in one point quadrature shell elements
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  thick : thickness of shell
  !
  !  zeta : [-1, +1] : pseudo-thickness parameter
  !
  !  ecurnloc(3,4) : local current nodal coordinate
  !
  !  eveloloc(5,4) : local nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  output:
  !  ------
  !  ltens2d(2,2) : velocity gradient tensor in 2d plane
  !  ipstrndot(3,1) : in plane strain rate [dx,dy,2dxy]
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: thick, zeta
  real(8), dimension(3,4), intent(in) :: ecurnloc
  real(8), dimension(5,4), intent(in) :: eveloloc

  real(8), dimension(2,2), intent(out) :: ltens2d
  real(8), dimension(3,1), intent(out) :: ipstrndot
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(2,4) :: bmat1pt
  real(8), dimension(2,4) :: bcmat1pt

  real(8) :: psibar
  real(8) :: dvxdx, dvxdy, dvydy, dvydx

  real(8) :: dx, dy, dxy

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  ltens2d(:,:)= 0.0d0
  ipstrndot(:,:)= 0.0d0


  ! compute b matrix: b matrix, b^c matrix, and b^s matrix
  ! ----------------
  ! compute b matrix: one point integration
  call getbmat1pt(ecurnloc, bmat1pt)
     ! input : ecurnloc
     ! output : bmat1pt

  ! compute b^c matrix: warpping correction
  call getbcmat1pt(ecurnloc, bcmat1pt)
     ! input : ecurnloc
     ! output : bcmat1pt

  ! ----------------------------------------------------------------
  ! compute in plane local velocity strain: eq (18) and see, research note
  ! --------------------------------------
  ! pseudo-thickness : eq (5)
  psibar= ( zeta * thick ) / 2.0d0

  dvxdx= 0.0d0 ! initialize
  dvydy= 0.0d0
  dvxdy =0.0d0
  dvydx =0.0d0

dx= 0.0d0
dy= 0.0d0
dxy= 0.0d0

  do inode=1, 4


! ### working
! --------------------------------------------------------------------------------------
! the original BLT element with warping correction
! ------------------------------------------------
     ! vx,x= b_xi v_xi + psibar*( b^c _xi v_xi + b_xi theta_yi )
     dvxdx= dvxdx + bmat1pt(1,inode)*eveloloc(1,inode) &
              + psibar*( bcmat1pt(1,inode)*eveloloc(1,inode) + bmat1pt(1,inode)*eveloloc(5,inode) )

     ! vy,y= b_yi v_yi + psibar*( b^c _yi v_yi - b_yi theta_xi )
     dvydy= dvydy + bmat1pt(2,inode)*eveloloc(2,inode) &
              + psibar*( bcmat1pt(2,inode)*eveloloc(2,inode) - bmat1pt(2,inode)*eveloloc(4,inode) )

     ! vx,y= b_yi v_xi + psibar*( b^c _yi v_xi + b_yi theta_yi)
     dvxdy= dvxdy + bmat1pt(2,inode)*eveloloc(1,inode) &
                  + psibar*( bcmat1pt(2,inode)*eveloloc(1,inode) + bmat1pt(2,inode)*eveloloc(5,inode) )

     ! vy,x= b_xi v_yi + psibar*( b^c _xi v_yi - b_xi theta_xi )
     dvydx= dvydx + bmat1pt(1,inode)*eveloloc(2,inode)  &
                  + psibar*( bcmat1pt(1,inode)*eveloloc(2,inode) - bmat1pt(1,inode)*eveloloc(4,inode) )


! --------------------------------------------------------------------------------------
! the original BLT element
! -------------------------
!     ! dx=
!     dx= dx + bmat1pt(1,inode)*eveloloc(1,inode) + psibar*( bmat1pt(1,inode)*eveloloc(5,inode) )
!
!     ! dy=
!     dy= dy + bmat1pt(2,inode)*eveloloc(2,inode) - psibar*( bmat1pt(2,inode)*eveloloc(4,inode) )
!
!     ! dxy=
!     dxy= dxy + 0.50d0 * ( bmat1pt(2,inode)*eveloloc(1,inode) + bmat1pt(1,inode)*eveloloc(2,inode) &
!                              + psibar*( bmat1pt(2,inode)*eveloloc(5,inode) - bmat1pt(1,inode)*eveloloc(4,inode) ) )

! --------------------------------------------------------------------------------------


  end do

! --------------------------------------------------------------------------------------
! the original BLT element with warping correction
! ------------------------------------------------

     ! set corotational velocity gradient tensir: ltens_ij = a v_i / a v_j 
     ltens2d(1,1)= dvxdx ! a v_x / a x  
     ltens2d(1,2)= dvxdy ! a v_x / a y

     ltens2d(2,1)= dvydx ! a v_y / a x
     ltens2d(2,2)= dvydy ! a v_y / a y

     ! set in plane velocity strain rate
     ipstrndot(1,1)= dvxdx ! d_x
     ipstrndot(2,1)= dvydy ! d_y
     ipstrndot(3,1)= dvxdy + dvydx ! 2.0 d_xy


! --------------------------------------------------------------------------------------
! the original BLT element
! ------------------------

!     ! set corotational velocity gradient tensir: ltens_ij = a v_i / a v_j 
!     ltens2d(1:2,1:2)= 0.0d0 ! temporary setting
!
!     ! set in plane velocity strain rate
!     ipstrndot(1,1)= dx
!     ipstrndot(2,1)= dy
!     ipstrndot(3,1)= 2.0d0 * dxy

! --------------------------------------------------------------------------------------




  return
end subroutine getstrndotbt1





subroutine getstrndotbt2(ecurnloc,eveloloc, tsstrndot)
  !=======================================================================
  !  getstrndotbt2 = compute transeverse shear local rate of deformation (velocity strain)
  !                  of tb shell with one point integration and local z method
  !                  for computational efficiency, explicit components form is used
  !
  !                  note:
  !                  ----
  !                  belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !                  advances in one point quadrature shell elements
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ecurnloc(3,4) : local current nodal coordinate
  !
  !  eveloloc(5,4) : local nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  output:
  !  ------
  !  tsstrndot(2,1) : tranverse shear strain rate [2dxz,2dyz]
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,4), intent(in) :: ecurnloc
  real(8), dimension(5,4), intent(in) :: eveloloc

  real(8), dimension(2,1), intent(out) :: tsstrndot
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(2,4) :: bmat1pt
  real(8), dimension(2,3,4) :: bsmat1pt

  real(8) :: dxz, dyz

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  tsstrndot(:,:)= 0.0d0


  ! compute b matrix: one point integration
  ! ----------------
  call getbmat1pt(ecurnloc, bmat1pt)
     ! input : ecurnloc
     ! output : bmat1pt

  ! compute b matrix: b^s matrix
  ! ----------------
  ! compute b^s matrix: transverse shear projection
  call getbsmat1pt(ecurnloc, bsmat1pt)
     ! input : ecurnloc
     ! output : bsmat1pt


  ! ----------------------------------------------------------------
  ! compute transverse shear local velocity strain: eq (34) and see, research note
  ! ----------------------------------------------
  dxz= 0.0d0 ! initialize
  dyz= 0.0d0

  do inode=1, 4

! --------------------------------------------------------------------------------------
! the original BLT element with shear correction
! ----------------------------------------------
!     ! dxz= 0.50d0 * ( b^s _x1i v_zi + b^s _x2i theta_xi + b^s _x3i theta_yi )
!     dxz= dxz + 0.50d0*( bsmat1pt(1,1,inode)*eveloloc(3,inode) &
!                       + bsmat1pt(1,2,inode)*eveloloc(4,inode) &
!                       + bsmat1pt(1,3,inode)*eveloloc(5,inode) )
!
!     ! dyz= 0.50d0 * ( b^s _y1i v_zi + b^s _y2i theta_xi + b^s _y3i theta_yi )
!     dyz= dyz + 0.50d0*( bsmat1pt(2,1,inode)*eveloloc(3,inode) &
!                       + bsmat1pt(2,2,inode)*eveloloc(4,inode) &
!                       + bsmat1pt(2,3,inode)*eveloloc(5,inode) )

! --------------------------------------------------------------------------------------
! the original BLT element
! ------------------------
     ! dxz=
     dxz= dxz + 0.50d0*( bmat1pt(1,inode)*eveloloc(3,inode) + 0.250d0*eveloloc(5,inode) )

     ! dyz=
     dyz= dyz + 0.50d0*( bmat1pt(2,inode)*eveloloc(3,inode) - 0.250d0*eveloloc(4,inode) )

! --------------------------------------------------------------------------------------

  end do

  ! ----------------------------------------------------------------
  ! set tranverse shear velocity strain
  ! -----------------------------------
  tsstrndot(1,1)= 2.0d0 * dxz
  tsstrndot(2,1)= 2.0d0 * dyz
  ! ----------------------------------------------------------------



  return
end subroutine getstrndotbt2