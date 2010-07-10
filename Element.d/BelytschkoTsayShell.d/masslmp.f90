! ==================================
! lumped mass matrix 
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine           elemasl1d          (nnode,nndof,ematpro,ecord, emasl1d)
! 2.  subroutine           elemasl2d          (nnode,nndof,ematpro,ecord, emasl2d)
! 3.  subroutine           elemaslbt          (nndof,ematpro,ecord,edisp, emaslbt)
! 4.  subroutine           elemasl3d          (nnode,nndof,ematpro,ecord, emasl3d)
!
! =========================================================================================================


subroutine elemasl1d(nnode,nndof,ematpro,ecord, emasl1d)
  !=======================================================================
  !  elemaslmp1 = compute lumped element mass matrix
  !
  !               1d element, total lagrangian
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nnode : node per element
  !
  !  nndof : total dof per node
  !
  !  ematpro(*) : element material properties
  !
  !  ecord(1,*) : element nodal coordinate
  !
  !  output:
  !  ------
  !  emasl1d(nndof*nnode,1) : lumped 1d element mass matrix
  !                            
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nnode, nndof
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(1,*), intent(in) :: ecord

  real(8), dimension(nndof*nnode,1), intent(out) :: emasl1d
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: denst, thick, xarea

  ! loop index
  integer :: irow, icol
  ! ====================================

  ! initialize
  emasl1d(:,:)= 0.0d0

  ! -------------------------------------------
  ! get material properties
  denst= ematpro(3)
  thick= ematpro(20)
  xarea= thick**2
  ! -------------------------------------------


  ! compute lumped mass matrix
  if( nnode==2 ) then ! 1d linear: 2 node
     emasl1d(1:nndof*nnode,1)= dsqrt( ( ecord(1,1)-ecord(1,2) )**2 ) * xarea * denst / 2.0d0

  else
     write(*,*) "not implemented yet: elemasl1d"
     write(nout5,*) "not implemented yet: elemasl1d"
     stop

  end if



  return
end subroutine elemasl1d





subroutine elemasl2d(nnode,nndof,ematpro,ecord, emasl2d)
  !=======================================================================
  !  elemasl2d = compute lumped element mass matrix
  !
  !              2d element, total lagrangian
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nnode : node per element
  !
  !  nndof : total dof per node
  !
  !  ematpro(*) : element material properties
  !
  !  ecord(2,*) : element nodal coordinate
  !
  !  output:
  !  ------
  !  emasl2d(nndof*nnode,1) : lumped 2d element mass matrix
  !                            
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nnode, nndof
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(2,*), intent(in) :: ecord

  real(8), dimension(nndof*nnode,1), intent(out) :: emasl2d
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: denst, thick
  real(8) :: area

  ! loop index
  integer :: irow, icol
  ! ====================================

  ! initialize
  emasl2d(:,:)= 0.0d0

  ! -------------------------------------------
  ! get material properties
  denst= ematpro(3)
  thick= ematpro(20)
  ! -------------------------------------------

  ! compute lumped mass matrix
  if( nnode==3 ) then ! 3 node tri.
     emasl2d(1:nndof*nnode,1)= area(nnode,ecord) * thick * denst / 3.0d0
 
  else if ( nnode==4 ) then ! 4 node quad.
     emasl2d(1:nndof*nnode,1)= area(nnode,ecord) * thick * denst / 4.0d0

  else
     write(*,*) "not available: elemasl2d"
     write(nout5,*) "not available: elemasl2d"
     stop

  end if



  return
end subroutine elemasl2d





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

  use preset
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





subroutine elemasl3d(nnode,nndof,ematpro,ecord, emasl3d)
  !=======================================================================
  !  elemasl3d = compute lumped element mass matrix
  !
  !              3d element, total lagrangian
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nnode : node per element
  !
  !  nndof : total dof per node
  !
  !  ematpro(*) : element material properties
  !
  !  ecord(3,*) : element nodal coordinate
  !
  !  output:
  !  ------
  !  emasl3d(nndof*nnode,1) : lumped 2d element mass matrix
  !                            
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nnode, nndof
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(3,*), intent(in) :: ecord

  real(8), dimension(nndof*nnode,1), intent(out) :: emasl3d
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: denst
  real(8) :: vol3d8nod

  ! loop index
  integer :: irow, icol
  ! ====================================

  ! initialize
  emasl3d(:,:)= 0.0d0


  ! -------------------------------------------
  ! get material properties
  denst= ematpro(3)
  ! -------------------------------------------

  ! compute lumped mass matrix
  if( nnode==8 ) then ! 8 node brick
     emasl3d(1:nndof*nnode,1)= vol3d8nod(ecord) * denst / 8.0d0

  else
     write(*,*) "not available: elemasl3d"
     write(nout5,*) "not available: elemasl3d"
     stop

  end if



  return
end subroutine elemasl3d