! ==========================
! apply boundary condition2 : element surface
! ==========================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine          getfbc3dbrkshl2      (optbc,opttrc,optele,msize, &
!                                               time,tmax,prmtft,adrs,coord3d,disp,bcnod,bcdof,bcval,bctmf, &
!                                               fext)
! 2.  subroutine          elefbc3dbrkshl2      (opttrc,optele,ecord3d,edisp3d,trac, efbc)
! 3.  subroutine          ele3dsurf2d4nod0     (ecord3d,edisp3d, area,nvec)
!
! =========================================================================================================



subroutine getfbc3dbrkshl2(optbc,opttrc,optele,msize, &
                           time,tmax,prmtft,adrs,coord3d,disp,bcnod,bcdof,bcval,bctmf, &
                           fext)
  !=======================================================================
  !  getfbc3dbrkshl2 = compute force bc for 3d bt shell and brick elements
  !
  !             note:
  !             -----
  !             element surface pressure
  !             element surface stress(traction)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optbc : boundary condition type option
  !          optbc=0 : initial bc
  !          optbc=1 : transient bc
  !
  !  opttrc : traction type option handler
  !           opttrc= 0 : pressure
  !           opttrc= 1 : traction
  !
  !  inoutput:
  !  --------
  !  fext(msize,1) : external force matrix
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optbc, opttrc
  integer, intent(in) :: optele
  integer, intent(in) :: msize
  real(8), intent(in) :: time, tmax
  real(8), dimension(*), intent(in) :: prmtft
  integer, dimension(*), intent(in) :: adrs
  real(8), dimension(3,*), intent(in) :: coord3d
  real(8), dimension(msize,1), intent(in) :: disp
  integer, dimension(8), intent(in) :: bcnod
  integer, dimension(6), intent(in) :: bcdof
  real(8), dimension(6), intent(in) :: bcval
  integer, dimension(6), intent(in) :: bctmf

  real(8), dimension(msize,1), intent(inout) :: fext
  ! ====================================
  ! local variable
  ! ==============
  integer, dimension(4) :: econc
  real(8), dimension(3,4) :: ecord3d
  integer, dimension(4) :: eadrs
  real(8), dimension(3,4) :: edisp3d
  real(8), dimension(3,1) :: trac

  real(8), dimension(12,1) :: efbc

  integer, dimension(12) :: sctr

  real(8) :: xdist
  real(8) :: tmftval
  integer :: nflagrmv

  integer :: iloc

  ! loop index
  integer :: inode, idime, idof
  ! ====================================

  ! initialize: do not initialize


  ! ---------------------------------------------------------------------------
  ! get element data
  ! ----------------
  ! get element surface nodal connectivity
  econc(1:4)= bcnod(1:4)
        
  ! get element nodal coordinate
  do inode=1, 4
     do idime=1, 3
        ecord3d(idime, inode)= coord3d(idime, econc(inode))
     end do
  end do

  ! get global element nodal address
  call geteadrs(4,adrs,econc, eadrs)
     ! input : 4(nnode),adrs,econc
     ! output : eadrs

  ! get element nodal displacement: translation only
  call geteleval1hybrid(msize,4,3,eadrs,disp, edisp3d)
     ! input : msize,4(nnode),3(nndof),eadrs,disp
     ! output : edisp3d

! ### sheperd's problem: for remote detonation
! ---------------------
!xdist= 0.0d0 ! initialzie
!do inode=1, 4
!   xdist= xdist + ecord3d(1,inode)
!end do
!xdist= xdist / 4.0d0

  ! set traction
  trac(1:3,1)= bcval(1:3)
  
  ! ---------------------------------------------------------------------------
  ! compute element nodal force
  ! ---------------------------
  call elefbc3dbrkshl2(opttrc,optele,ecord3d,edisp3d,trac, efbc)
     ! input : opttrc,optele,ecord3d,edisp3d,trac
     ! output : efbc

  ! ---------------------------------------------------------------------------
  ! set result: translational dof only
  ! ----------
  ! get scatter matrix
  call getsctrhybrid(4,3,eadrs, sctr)
     ! input : 4(nnode),3(nndof),eadrs
     ! output : sctr

  do idof=1, 3

     ! skip for unprescribed traction dof
     if( opttrc==1 .and. bcdof(idof)==0 ) cycle

     ! ------------------------------------------------------------------
     ! time function factor
     ! --------------------
     if( optbc == 0 ) then ! initial bc
        tmftval= 1.0d0
        nflagrmv= 0

     else if ( optbc == 1 .and. opttrc == 0 ) then ! transient pressure
        ! get time function
        call gettmftval(bctmf(1),prmtft,time,tmax,xdist, tmftval,nflagrmv)
           ! input : bctmf(1),prmtft,time,tmax,0.0d0(dist)
           ! output : tmftval,nflagrmv

     else if ( optbc == 1 .and. opttrc==1 ) then ! transient traction
        ! get time function
        call gettmftval(bctmf(idof),prmtft,time,tmax,xdist, tmftval,nflagrmv)
           ! input : bctmf(idof),prmtft,time,tmax,0.0d0(dist)
           ! output : tmftval,nflagrmv

     else
        write(*,*) "wrong optbc or opttrc: getfbc3dbrkshl2"
        write(nout5,*) "wrong optbc or opttrc: getfbc3dbrkshl2"
        stop
     
     end if

     ! ------------------------------------------------------------------
     ! add to external force vector
     ! ----------------------------
     do inode=1, 4

        ! get address
        iloc= 3*(inode-1) + idof

        ! add element fbc to fext
        fext(sctr(iloc),1)= fext(sctr(iloc),1) + efbc(iloc,1) * tmftval

     end do

     ! ------------------------------------------------------------------

  end do ! end do idof

  ! ---------------------------------------------------------------------------



  return
end subroutine getfbc3dbrkshl2





subroutine elefbc3dbrkshl2(opttrc,optele,ecord3d,edisp3d,trac, efbc)
  !=======================================================================
  !  elefbc3dbrkshl2 = compute 3d bt shell or birck element force due to
  !                    element surface stress or pressure
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opttrc : traction type option handler
  !           opttrc= 0 : pressure
  !           opttrc= 1 : traction
  !
  !  output:
  !  ------
  !  efbc(12,1) : element nodal force vector ( translational dof only)
  !                          
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opttrc, optele
  real(8), dimension(3,4), intent(in) :: ecord3d
  real(8), dimension(3,4), intent(in) :: edisp3d
  real(8), dimension(3,1), intent(in) :: trac

  real(8), dimension(12,1), intent(out) :: efbc
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: area
  real(8), dimension(3,1) :: nvec
  real(8), dimension(3,1) :: tracvec

  ! ====================================

  ! initialize
  efbc(:,:)= 0.0d0


  ! ---------------------------------------------------------------------
  ! compute surface area and normal vector
  ! --------------------------------------
  call ele3dsurf2d4nod0(ecord3d,edisp3d, area,nvec)
     ! input : ecord3d,edisp3d
     ! output : area,nvec

  ! ---------------------------------------------------------------------
  ! get global traction vector
  ! --------------------------
  select case(opttrc)
  case(0) ! pressure
     ! set pressure vector
     tracvec(1:3,1)= trac(1,1) * nvec(1:3,1)

  case(1) ! traction
     tracvec(1:3,1)= trac(1:3,1)

  case default
     write(*,*) "wrong opttrc: elefbc3dbrkshl2"
     write(nout5,*) "wrong opttrc: elefbc3dbrkshl2"
     stop

  end select


  ! ---------------------------------------------------------------------
  ! set element nodal force
  ! -----------------------
  efbc(1:3,1)= area *  tracvec(1:3,1) / 4.0d0    ! f_x, f_y and f_z at node 1
  efbc(4:6,1)= area *  tracvec(1:3,1) / 4.0d0    ! f_x, f_y and f_z at node 2
  efbc(7:9,1)= area *  tracvec(1:3,1) / 4.0d0    ! f_x, f_y and f_z at node 3
  efbc(10:12,1)= area *  tracvec(1:3,1) / 4.0d0  ! f_x, f_y and f_z at node 4

  ! ---------------------------------------------------------------------



  return
end subroutine elefbc3dbrkshl2





subroutine ele3dsurf2d4nod0(ecord3d,edisp3d, area,nvec)
  !=======================================================================
  !  ele3dsurf2d4nod0 = compute 3d bt shell or birck element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !                          
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,4), intent(in) :: ecord3d
  real(8), dimension(3,4), intent(in) :: edisp3d

  real(8), intent(out) :: area
  real(8), dimension(3,1), intent(out) :: nvec
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,4) :: ecurn
  real(8), dimension(3,3) :: locbvec
  real(8), dimension(3,4) :: ecurnloc

  real(8), dimension(3,1) :: norlvecloc

  ! loop index
  ! ====================================

  ! initialize
  area= 0.0d0
  nvec(:,:)= 0.0d0


  ! -------------------------------------------------------------------------------
  ! construct corotational coordinate
  ! ---------------------------------
  ! set element current coordinate
  ecurn(1:3,1:4)= ecord3d(1:3,1:4) + edisp3d(1:3,1:4)

  ! compute co rotational local base vector: locbvec
  call getlocbvecbt(ecurn, locbvec)
     ! input : ecurn
     ! output : locbvec

  ! -------------------------------------------------------------------------------
  ! compute surface area
  ! --------------------
  ! get local nodal coordinate: glb -> loc
  call glb2locnodv(3,4,3,locbvec,ecurn, ecurnloc)
     ! input : 3(ndime),4(nnode),3(ntrndof),locbvec,ecurn
     ! output : ecurnloc

  ! compute area
  area= 0.50d0*( (ecurnloc(1,3)-ecurnloc(1,1))*(ecurnloc(2,4)-ecurnloc(2,2)) &
                +(ecurnloc(1,2)-ecurnloc(1,4))*(ecurnloc(2,3)-ecurnloc(2,1)) )

  ! -------------------------------------------------------------------------------
  ! compute surface unit normal vector
  ! ---------------------------------
  ! set normal to corotational element surface
  norlvecloc(1,1)= 0.0d0  ! x^
  norlvecloc(2,1)= 0.0d0  ! y^
  norlvecloc(3,1)= -1.0d0 ! z^

  ! convert local to global
  call loc2glbv(3,3,locbvec,norlvecloc, nvec)
     ! input : 3(ndime),3(ntrndof),locbvec,norlvecloc
     ! output : nvec



  return
end subroutine ele3dsurf2d4nod0