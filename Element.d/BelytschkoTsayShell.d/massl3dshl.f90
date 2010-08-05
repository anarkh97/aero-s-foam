! ==================================
! lumped mass 3d shell
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine           getmaslbt1           (ielem,optele,msize,mnode,ndime,nnode,nndof,matpro,matnum,adrs,elestatus,conec,coord,disp, &
!                                                mass)
! 2.  subroutine           elemaslbt            (nndof,ematpro,ecord,edisp, emaslbt)
!
! =========================================================================================================



subroutine getmaslbt1(ielem,optele,msize,mnode,ndime,nnode,nndof,matpro,matnum,adrs,elestatus,conec,coord,disp, &
                      mass)
  !=======================================================================
  !  getmaslbt1 = compute lumped element mass matrix for bt shell
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  inoutput:
  !  --------
  !  mass(msize,1) : global lumped mass matrix
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ielem,optele,msize,mnode,ndime,nnode,nndof
  real(8), dimension(20,*), intent(in) :: matpro
  integer, dimension(*), intent(in) :: matnum
  integer, dimension(*), intent(in) :: adrs
  integer, dimension(*), intent(in) :: elestatus
  integer, dimension(mnode,*), intent(in) :: conec
  real(8), dimension(ndime,*), intent(in) :: coord
  real(8), dimension(msize,1), intent(in) :: disp
  ! ------------------------------------

  real(8), dimension(msize,1), intent(inout) :: mass
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(20) :: ematpro
  
  integer, dimension(nnode) :: eadrs
  integer, dimension(nnode) :: econc
  real(8), dimension(ndime,nnode) :: ecord
  real(8), dimension(nndof,nnode) :: edisp

  real(8), dimension(nndof*nnode,1) :: emasl
  integer, dimension(nndof*nnode) :: sctr

  ! loop index
  integer :: idof
  ! ====================================

  ! initialize: do not initialize


  ! ---------------------------------------------------------------
  ! get material properties
  ! -----------------------
  call getmatprop0(ielem,matpro,matnum, ematpro) 
     ! input : ielem,matpro,matnum
     ! output : ematpro

  ! get element nodal coordinates and connectivity
  ! ----------------------------------------------
  call geteledata1hybrid(ielem,mnode,ndime,nnode,conec,coord, econc,ecord) 
     ! input : ielem,mnode,ndime,nnode,conec,coord
     ! output : econc,ecord

  ! get global element nodal address
  call geteadrs(nnode,adrs,econc, eadrs)
     ! input : nnode,adrs,econc
     ! output : eadrs

  ! get element nodal displacement
  ! ------------------------------
  call geteleval1hybrid(msize,nnode,nndof,eadrs,disp, edisp)
     ! input : msize,nnode,nndof,eadrs,disp
     ! output : edisp


  ! get bt shell element lumped mass
  ! --------------------------------
  call elemaslbt(nndof,ematpro,ecord,edisp, emasl)
     ! input : nndof,ematpro,ecord,edisp
     ! output : emasl

  ! ---------------------------------------------------------------
  ! get scatter matrix
  ! ------------------
  call getsctrhybrid(nnode,nndof,eadrs, sctr)
     ! input : nnode,nndof,eadrs
     ! output : sctr

  ! set global lumped mass matrix
  ! -----------------------------
  do idof=1, nndof*nnode
     mass(sctr(idof),1)= mass(sctr(idof),1) + emasl(idof,1)
  end do
  ! ---------------------------------------------------------------



  return
end subroutine getmaslbt1





subroutine elemaslbt(nndof,ematpro,ecord,edisp, emaslbt)
  !=======================================================================
  !  elemaslbt = compute argumented rotation lumped mass matrix
  !              for 5 or 6 dof bt shell element
  !
  !              note:
  !              ---
  !              hughes, cohen and haroun, NED, 1978, vol. 46, pp. 203-222
  !              reduced and selective integration techniques in
  !              the finite element analysis of plates
  !              (see, ch 4.1)
  !
  !              kennedy, belytschko and lin, NED, 1986, vol. 97, pp. 1-24
  !              recent developments in explicit finite element techniques and
  !              their application to reactor structures
  !              (see, eq (16))
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndof : the total number of dof per node
  !
  !  ematpro(*) : element material properties
  !
  !  ecord(3,4) : element nodal coordinate
  !
  !  edisp(nndof,4) : element nodal displacement data
  !
  !  output:
  !  ------
  !  emaslbt(nndof*4,1) : lumped element mass matrix
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndof
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(nndof,4), intent(in) :: edisp

  real(8), dimension(nndof*4,1), intent(out) :: emaslbt
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: denst, thick

  real(8), dimension(3,4) :: ecurn
  real(8), dimension(3,4) :: ecordloc
  real(8), dimension(3,3) :: locbvec

  real(8), dimension(nndof,nndof) :: loctens, glbtens
  real(8) :: area, mpnd, alpha

  real(8), dimension(nndof*4,1) :: emassloc

  integer :: iloc

  ! loop index
  integer :: inode, idof
  ! ====================================

  ! initialize
  emaslbt(:,:)= 0.0d0

  ! -------------------------------------------
  ! get material properties
  denst= ematpro(3)
  thick= ematpro(20)
  ! -------------------------------------------

  ! get current nodal coordinate
  ecurn(1:3,1:4)= ecord(1:3,1:4) + edisp(1:3,1:4)


  ! compute co rotational local base vector: locbvec
  ! ---------------------------------------
  call getlocbvecbt(ecurn, locbvec)
     ! input : ecurn
     ! output : locbvec

  ! --------------------------------------------------------------
  ! convert global to local
  ! -----------------------
  ! get local nodal coordinate
  call glb2locnodv(3,4,3,locbvec,ecurn, ecordloc)
     ! input : 3(ndime),4(nnode),3(ntrndof),locbvec,ecurn
     ! output : ecordloc
  ! --------------------------------------------------------------

  ! compute area
  area= 0.50d0 * ( (ecordloc(1,3)-ecordloc(1,1))*(ecordloc(2,4)-ecordloc(2,2)) &
                  +(ecordloc(1,2)-ecordloc(1,4))*(ecordloc(2,3)-ecordloc(2,1)) )

  ! compute total element mass per node
  mpnd= area * denst * thick / 4.0d0

  ! compute scaling factor alpha
  alpha= area / 8.0d0 ! (thick**2 + area)/12.0d0

  ! set local element mass
  ! ----------------------
  do inode=1, 4

     ! address
     iloc= nndof*(inode-1)

     emassloc(iloc+1,1)= mpnd * 1.0d0 ! mass_trn.x
     emassloc(iloc+2,1)= mpnd * 1.0d0 ! mass_trn.y
     emassloc(iloc+3,1)= mpnd * 1.0d0 ! mass_trn.z

     emassloc(iloc+4,1)= mpnd * alpha ! mass_rot.x
     emassloc(iloc+5,1)= mpnd * alpha ! mass_rot.y

     if ( nndof == 6 ) then
        emassloc(iloc+6,1)= mpnd * 2.0d0* alpha ! mass_rot.z
     end if

  end do

  ! --------------------------------------------------------------
  ! convert local to global
  ! -----------------------
  do inode=1, 4

     ! address
     iloc= nndof*(inode-1)

     ! get local nodal mass
     loctens(:,:)= 0.0d0
     do idof=1, nndof
        loctens(idof,idof)= emassloc(iloc+idof,1)
     end do

     ! convert
     call loc2glbtens(3,nndof,locbvec,loctens, glbtens)
        ! input : 3(ndime),nndof,locbvec,loctens
        ! output : glbtens
 
     ! set global nodal mass
     do idof=1, nndof
        emaslbt(iloc+idof,1)= glbtens(idof,idof)
     end do

  end do
  ! --------------------------------------------------------------



  return
end subroutine elemaslbt