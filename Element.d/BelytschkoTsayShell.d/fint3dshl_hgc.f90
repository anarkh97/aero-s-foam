! ==================================
! hourglass control force: bt shell
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine         elefhgcbt1           (prmhgc,delt,ematpro,nndof,ecord,edisp,evelo, &
!                                              evoit1, &
!                                              efhgc)
! 2.  subroutine         updhgcstrsbt         (prmhgc,delt,ematpro,locbvec,ecurn,evelo, hgcvoitloc)
! 3.  subroutine         gethgcstrsdotbt      (prmhgc,ematpro,ecurnloc,eveloloc, hgcbstrsdot,hgcmstrsdot)
! 4.  subroutine         gethgcstrndotbt      (ecurnloc,eveloloc, hgcbstrndot,hgcmstrndot)
! 5.  subroutine         gethgconstbt         (prmhgc,ematpro,ecurnloc, c1,c2,c3)
! 6.  subroutine         gqfhgcbt             (locbvec,ecurn,evelo,hgcvoitloc, gqfhgc)
!
! =========================================================================================================



subroutine elefhgcbt1(prmhgc,delt,ematpro,nndof,ecord,edisp,evelo, &
                      evoit1, &
                      efhgc)
  !=======================================================================
  !  elefhgcbt1 = compute hourglass control nodal force for bt shell
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !
  !  inoutput:
  !  --------
  !  evoit1(6,*) : voight form hourglass control stress
  !
  !  output:
  !  ------
  !  efhgc(nndof*4,1) : hourglass control nodal force
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: prmhgc
  real(8), intent(in) :: delt
  real(8), dimension(20), intent(in) :: ematpro
  integer, intent(in) :: nndof
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(6,4), intent(in) :: edisp
  real(8), dimension(6,4), intent(in) :: evelo
  ! -------------------------------------

  real(8), dimension(6,*), intent(inout) :: evoit1

  ! -------------------------------------

  real(8), dimension(24,1), intent(out) :: efhgc
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,4) :: ecurn
  real(8), dimension(3,3) :: locbvec
  real(8), dimension(6,1) :: hgcvoitloc
  ! ====================================

  ! initialize
  efhgc(:,:)= 0.0d0


  ! get current nodal coordinate: mid surface
  ecurn(1:3,1:4)= ecord(1:3,1:4) + edisp(1:3,1:4)


  ! compute co rotational local base vector: locbvec
  ! ---------------------------------------
  call getlocbvecbt(ecurn, locbvec)
     ! input : ecurn
     ! output : locbvec

  ! -------------------------------------
  ! get history variable
  ! --------------------
  hgcvoitloc(1:6,1)= evoit1(1:6,1)
  ! -------------------------------------


  ! update hourglass control stresses
  ! ---------------------------------
  call updhgcstrsbt(prmhgc,delt,ematpro,locbvec,ecurn,evelo, hgcvoitloc)
     ! input : prmhgc,delt,ematpro,locbvec,ecurn,evelo
     ! output : hgcvoitloc


  ! -------------------------------------
  ! set history variable
  ! --------------------
  evoit1(1:6,1)= hgcvoitloc(1:6,1)
  ! -------------------------------------


  ! compute local hourglass control nodal forces
  ! --------------------------------------------
  ! local hourglass control nodal forces 
  call gqfhgcbt(locbvec,ecurn,evelo,hgcvoitloc, efhgc)
     ! input : locbvec,ecurn,evelo,hgcvoitloc
     ! output : efhgc



  return
end subroutine elefhgcbt1



subroutine updhgcstrsbt(prmhgc,delt,ematpro,locbvec,ecurn,evelo, hgcvoitloc)
  !=======================================================================
  !  updhgcstrsbt = update hourglass control stress
  !
  !               note:
  !               ----
  !               chiang, m.s. thesis, northwestern univ., 1992
  !               advances in one point quadrature shell element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  inoutput:
  !  --------
  !  hgcvoitloc(6,1) : voightform hourglass control stress
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: prmhgc
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(3,3), intent(in) :: locbvec
  real(8), dimension(3,4), intent(in) :: ecurn
  real(8), dimension(6,4), intent(in) :: evelo
  
  ! ------------------------------------

  real(8), dimension(6,1), intent(inout) :: hgcvoitloc
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,4) :: ecurnloc
  real(8), dimension(6,4) :: eveloloc0
  real(8), dimension(5,4) :: eveloloc

  real(8), dimension(3,1) :: hgcbstrsdot
  real(8), dimension(2,1) :: hgcmstrsdot
  ! ====================================

  ! initialize: do not initialize

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
  ! compute hourglass stress rate
  ! -----------------------------
  call gethgcstrsdotbt(prmhgc,ematpro,ecurnloc,eveloloc, hgcbstrsdot,hgcmstrsdot)
     ! input : prmhgc,ematpro,ecurnloc,eveloloc
     ! output : hgcbstrsdot,hgcmstrsdot


  ! -------------------------------------------------------------------------------
  ! update hourglass control stress
  ! -------------------------------
  hgcvoitloc(1,1)= hgcvoitloc(1,1) + hgcbstrsdot(1,1) * delt ! sig_x
  hgcvoitloc(2,1)= hgcvoitloc(2,1) + hgcbstrsdot(2,1) * delt ! sig_y
  hgcvoitloc(3,1)= 0.0d0 ! sig_z (plane stress)

  hgcvoitloc(4,1)= hgcvoitloc(4,1) + hgcmstrsdot(2,1) * delt ! sig_yz
  hgcvoitloc(5,1)= hgcvoitloc(5,1) + hgcmstrsdot(1,1) * delt ! sig_xz
  hgcvoitloc(6,1)= hgcvoitloc(6,1) + hgcbstrsdot(3,1) * delt ! sig_xy
  ! -------------------------------------------------------------



  return
end subroutine updhgcstrsbt





subroutine gethgcstrsdotbt(prmhgc,ematpro,ecurnloc,eveloloc, hgcbstrsdot,hgcmstrsdot)
  !=======================================================================
  !  gethgcstrsdotbt = compute hourglass stress rate
  !
  !                 note:
  !                 ----
  !                 chiang, m.s. thesis, northwestern univ., 1992
  !                 advances in one point quadrature shell element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !  hgcbstrsdot(3,1) : bending hourglass control stress rate
  !
  !  hgcmstrsdot(2,1) : membrane hourglass control stress rate
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: prmhgc
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(3,4), intent(in) :: ecurnloc
  real(8), dimension(5,4), intent(in) :: eveloloc

  real(8), dimension(3,1), intent(out) :: hgcbstrsdot
  real(8), dimension(2,1), intent(out) :: hgcmstrsdot
  ! ====================================
  ! local variable
  ! ==============
  real(8):: c1, c2, c3

  real(8), dimension(3,1) :: hgcbstrndot
  real(8), dimension(2,1) :: hgcmstrndot
  ! ====================================

  ! initialize
  hgcbstrsdot(:,:)= 0.0d0
  hgcmstrsdot(:,:)= 0.0d0

  ! compute hourglass mode control constant
  call gethgconstbt(prmhgc,ematpro,ecurnloc, c1,c2,c3)
     ! input : prmhgc,ematpro,ecurnloc
     ! output : c1,c2,c3


  ! compute hourglass strain rate
  call gethgcstrndotbt(ecurnloc,eveloloc, hgcbstrndot,hgcmstrndot)
     ! input : ecurnloc,eveloloc
     ! output : hgbstrndot,hgmstrndot 


  ! -------------------------------------------------------------
  ! compute hourglass stress rate
  ! bending
  hgcbstrsdot(1,1)= c1*hgcbstrndot(1,1)
  hgcbstrsdot(2,1)= c1*hgcbstrndot(2,1)
  hgcbstrsdot(3,1)= c2*hgcbstrndot(3,1)

  ! membrane
  hgcmstrsdot(1,1)= c3*hgcmstrndot(1,1)
  hgcmstrsdot(2,1)= c3*hgcmstrndot(2,1)
  ! -------------------------------------------------------------


  return
end subroutine gethgcstrsdotbt





subroutine gethgcstrndotbt(ecurnloc,eveloloc, hgcbstrndot,hgcmstrndot)
  !=======================================================================
  !  gethgcstrndotbt = compute hourglass strain rate for memebrane and bending
  !                  in membrane strain rate, warping correction is considered
  !
  !                  note:
  !                  ----
  !                  chiang, m.s. thesis, northwestern univ., 1992
  !                  advances in one point quadrature shell element
  !                  (see, eq (47), (48))
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !  hgcbstrndot(3,1) : bending hourglass strain rate
  !
  !  hgcmstrndot(2,1) : membrane hourglass strain rate
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,4), intent(in) :: ecurnloc
  real(8), dimension(5,4), intent(in) :: eveloloc

  real(8), dimension(3,1), intent(out) :: hgcbstrndot
  real(8), dimension(2,1), intent(out) :: hgcmstrndot
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(2,4) :: ecurn2dloc
  real(8), dimension(4,1) :: gamma
  real(8), dimension(4,1) :: s
  real(8) :: zgamma
  real(8), dimension(3,1) :: gammatv
  real(8), dimension(2,1) :: gammatrot
  real(8), dimension(2,1) :: strot

  ! loop index
  integer :: inode, idime
  ! ====================================

  ! initialize
  hgcbstrndot(:,:)= 0.0d0
  hgcmstrndot(:,:)= 0.0d0


  ! set 2d coordinate only
  ecurn2dloc(1:2,1:4)= ecurnloc(1:2,1:4)

  ! compute gamma projection oprator
  call getgamma4nod(ecurn2dloc, gamma)
     ! input : ecurn2dloc
     ! output : gamma

  ! define rigid body motion vector
  s(1:4,1)= 1.0d0

  ! compute z_i gamma_i
  zgamma= 0.0d0 ! initialize
  do inode=1, 4
     zgamma= zgamma + ecurnloc(3,inode)*gamma(inode,1)
  end do

 
  ! compute gamma^t v
  gammatv(:,:)= 0.0d0 ! initialize
  do inode=1, 4
     do idime=1, 3
        gammatv(idime,1)= gammatv(idime,1) + gamma(inode,1)*eveloloc(idime,inode)
     end do
  end do


  ! compute gamma^t theta
  gammatrot(:,:)= 0.0d0 ! initialize
  do inode=1, 4
     do idime=1, 2 ! only have x and y dof in rotation
        gammatrot(idime,1)= gammatrot(idime,1) + gamma(inode,1)*eveloloc(3+idime,inode)
     end do
  end do


  ! compute s ^t theta: s= [1 1 1 1] for rigid body motion
  strot(:,:)= 0.0d0 ! initialize
  do inode=1, 4
     do idime=1, 2 ! only have x and y dof in rotation
        strot(idime,1)= strot(idime,1) + s(inode,1)*eveloloc(3+idime,inode)
     end do
  end do


  ! -------------------------------------------------------------
  ! compute strain rate
  ! -------------------
  ! bending hourglass strain rate
  hgcbstrndot(1,1)= gammatrot(1,1) 
  hgcbstrndot(2,1)= gammatrot(2,1)
  hgcbstrndot(3,1)= gammatv(3,1)


  ! membrane hourglass strain rate
  hgcmstrndot(1,1)= gammatv(1,1) - (1.0d0/4.0d0)*zgamma*strot(2,1)
  hgcmstrndot(2,1)= gammatv(2,1) + (1.0d0/4.0d0)*zgamma*strot(1,1)
  ! -------------------------------------------------------------


  return
end subroutine gethgcstrndotbt





subroutine gethgconstbt(prmhgc,ematpro,ecurnloc, c1,c2,c3)
  !=======================================================================
  !  gethgconstbt = compute constant for hourglass mode control
  !                 warping correction is considered
  !
  !                 ref.:
  !                 ----
  !                 chiang, m.s. thesis, northwestern univ., 1992
  !                 advances in one point quadrature shell element
  !                 (see, eq (49))
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !  c1, c2, c3 : hourglass mode control constant
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: prmhgc
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(3,4), intent(in) :: ecurnloc

  real(8), intent(out) :: c1,c2,c3
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: rt, rw, rm, rk

  real(8) :: young, poiss, thick

  real(8) :: mu

  real(8) :: area

  real(8), dimension(2,4) :: bmat1pt

  real(8) :: btb

  real(8) :: const1, const2, const3

  ! loop index
  integer :: inode, idime
  ! ====================================

  ! initialize
  c1= 0.0d0
  c2= 0.0d0
  c3= 0.0d0

  ! ------------------------------------------------
  ! hourglass mode control parameter
  rt= prmhgc(1) ! 0.010d0 
  rw= prmhgc(2) ! 0.010d0 
  rm= prmhgc(3) ! 0.010d0
  ! ------------------------------------------------
  ! get material properties
  young= ematpro(1)
  poiss= ematpro(2)
  thick= ematpro(20)

  ! shear correction factor
  rk= ematpro(19) ! 0.840d0 
  ! ------------------------------------------------

  ! compute shear modulus
  mu= young /( 2.0d0*(1.0d0+poiss) )

  ! compute area
  area= 0.50d0*( (ecurnloc(1,3)-ecurnloc(1,1))*(ecurnloc(2,4)-ecurnloc(2,2)) &
                +(ecurnloc(1,2)-ecurnloc(1,4))*(ecurnloc(2,3)-ecurnloc(2,1)) )

  ! compute one point integrated b matrix
  call getbmat1pt(ecurnloc, bmat1pt)
     ! input : ecurnloc
     ! output : bmat1pt

  ! compute b_xi b_xi + b_yi b_yi
  btb= 0.0d0 ! initialize
  do inode=1, 4
    do idime=1, 2
        btb= btb + bmat1pt(idime,inode) * bmat1pt(idime,inode)
    end do
  end do


  ! -------------------------------------------------------------
  ! compute c1
  const1= rt/192.0d0
  const2= young * thick**3 * area 
  const3= 1.0d0 + (2.0d0*rk*area)/(3.0d0*thick**2)
  c1= const1 * const2 * const3 * btb

  ! compute c2
  const1= rw/12.0d0
  const2= rk * mu * thick**3
  c2= const1 * const2 * btb

  ! compute c3
  const1= rm/8.0d0
  const2= young * thick * area
  c3= const1 * const2 * btb
  ! -------------------------------------------------------------



  return
end subroutine gethgconstbt





subroutine gqfhgcbt(locbvec,ecurn,evelo,hgcvoitloc, gqfhgc) 
  !=======================================================================
  !  elefhgclocbt = compute hourglass control force for bt shell
  !
  !                 note:
  !                 ----
  !                 chiang, m.s. thesis, northwestern univ., 1992
  !                 advances in one point quadrature shell element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !  efhgcloc(20,1) : element hourglass control force in local coordinate
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,3), intent(in) :: locbvec
  real(8), dimension(3,4), intent(in) :: ecurn
  real(8), dimension(6,4), intent(in) :: evelo

  ! -----------------------------------

  real(8), dimension(6,1), intent(inout) :: hgcvoitloc

  ! -----------------------------------

  real(8), dimension(24,1), intent(out) :: gqfhgc
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,4) :: ecurnloc
  real(8), dimension(6,4) :: eveloloc0
  real(8), dimension(5,4) :: eveloloc

  real(8), dimension(3,1) :: hgcbstrsloc
  real(8), dimension(2,1) :: hgcmstrsloc

  real(8), dimension(2,4) :: ecurn2dloc
  real(8), dimension(4,1) :: gamma
  real(8), dimension(4,1) :: s
  real(8) :: zgamma

  real(8), dimension(24,1) :: gqfhgcloc
  real(8), dimension(6,1) :: locvec, glbvec
  integer :: iloc

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  gqfhgc(:,:)= 0.0d0

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
  ! set local hourglass control stress
  ! ----------------------------------
  hgcbstrsloc(1,1)= hgcvoitloc(1,1) ! sig_x
  hgcbstrsloc(2,1)= hgcvoitloc(2,1) ! sig_y

  hgcmstrsloc(2,1)= hgcvoitloc(4,1) ! sig_yz
  hgcmstrsloc(1,1)= hgcvoitloc(5,1) ! sig_xz
  hgcbstrsloc(3,1)= hgcvoitloc(6,1) ! sig_xy

  ! -------------------------------------------------------------------------------
  ! compute gamma projection oprator
  ! --------------------------------
  ! set 2d coordinate only
  ecurn2dloc(1:2,1:4)= ecurnloc(1:2,1:4)

  ! compute gamma projection
  call getgamma4nod(ecurn2dloc, gamma)
     ! input : ecurn2dloc
     ! output : gamma

  ! define rigid body motion vector
  s(1:4,1)= 1.0d0

  ! compute z_i gamma_i
  zgamma= 0.0d0 ! initialize
  do inode=1, 4
     zgamma= zgamma + ecurnloc(3,inode)*gamma(inode,1)
  end do


  ! -------------------------------------------------------------------------------
  ! set local hourglass stabilization force
  ! ---------------------------------------
  do inode=1, 4

     ! address
     iloc= 6*(inode-1)

     ! compute local hourglass control force
     gqfhgcloc(iloc+1,1)= gamma(inode,1)*hgcmstrsloc(1,1)  ! f^hgc_x
     gqfhgcloc(iloc+2,1)= gamma(inode,1)*hgcmstrsloc(2,1)  ! f^hgc_y
     gqfhgcloc(iloc+3,1)= gamma(inode,1)*hgcbstrsloc(3,1)  ! f^hgc_z

     ! compute local hourglass moment
     gqfhgcloc(iloc+4,1)= gamma(inode,1)*hgcbstrsloc(1,1) + (1.0d0/4.0d0)*zgamma*s(inode,1)*hgcmstrsloc(2,1)  ! m^hgc_x
     gqfhgcloc(iloc+5,1)= gamma(inode,1)*hgcbstrsloc(2,1) - (1.0d0/4.0d0)*zgamma*s(inode,1)*hgcmstrsloc(1,1)  ! m^hgc_y
     gqfhgcloc(iloc+6,1)= 0.0d0 ! m^hgc_z

  end do

  ! -------------------------------------------------------------------------------
  ! convert local to global
  ! -----------------------
  do inode=1, 4

     ! address
     iloc= 6*(inode-1)

     ! get local nodal mass
     locvec(1:6,1)= gqfhgcloc(iloc+1:iloc+6,1)

     ! convert
     call loc2glbv(3,6,locbvec,locvec, glbvec)
        ! input : 3(ndime),6(ntrndof),locbvec,locvec
        ! output : glbvec
 
     ! set global nodal mass
     gqfhgc(iloc+1:iloc+6,1)= glbvec(1:6,1)

  end do
  ! -------------------------------------------------------------------------------



  return
end subroutine gqfhgcbt