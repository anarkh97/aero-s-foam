! ==================================
! fem data manipulation  
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine          getsctr                    (nnode,nndof,econc, sctr)
! 2.  subroutine          geteledata0                (ielem,nnode,conec, econc) 
! 3.  subroutine          geteledata1                (ielem,ndime,nnode,conec,coord, econc,ecord) 
! 4.  subroutine          geteleval1                 (msize,nnode,nndof,econc,disp, edisp) 
! 5.  subroutine          geteleval2                 (msize,nnode,nndof,econc,disp,velo, edisp,evelo) 
! 6.  subroutine          geteleval3                 (msize,nnode,nndof,econc,disp,velo,accl, edisp,evelo,eaccl) 
! 7.  subroutine          getmatprop0                (ielem,matpro,matnum, ematpro) 
! 8.  subroutine          getmatprop1                (ielem,matpro,matnum, young,poiss,denst,thick) 
! 9.  integer function    ismbr1di                   (nndx,array1d,inumber)
! 10. integer function    ismbr2di                   (nrndx,ncndx,array2d,inumber)
! 11. integer function    ismbr1d1d0                 (nndx0,nndx1,array1d0,array1d1)
! 12. subroutine          geteleval1vec              (msize,nnode,nndof,econc,disp, edisp)
!
! =========================================================================================================



subroutine getsctr(nnode,nndof,econc, sctr)
  !=======================================================================
  !  getsctr= get scatter matrix
  !           : convert local node number to system equation address
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nnode : total number of node per element
  !  
  !  nndof : total number of dof per node
  !
  !  econc(*) : element connectivity matrix
  !
  !  output:
  !  ------
  !  sctr(nnode*nndof) : scatter matrix
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nnode, nndof 
  integer, dimension(*), intent(in) :: econc

  integer, dimension(nnode*nndof), intent(out) :: sctr
   ! ====================================
  ! local variable
  ! ==============
  integer :: indx
  
  ! loop index
  integer :: inode, idof
  ! ====================================

  ! initialize
  sctr(:)=0
  
  indx=0 ! initialize
  do inode=1, nnode
     do idof=1, nndof

	    indx= indx+1 ! count
		 sctr(indx)= nndof*(econc(inode)-1) + idof

     end do
  end do

  return
end subroutine getsctr





subroutine geteledata0(ielem,nnode,conec, econc) 
  !=======================================================================
  !  geteledata0= get element data: connectivity
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ielem : current element number
  !
  !  nnode : node per element
  !
  !  conec(nnode,*) : global element connectivity data
  !
  !  output:
  !  ------
  !  econc : element connectivity data
  !          econc(inode) : global node number of current i th node
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ielem
  integer, intent(in) :: nnode  
  integer, dimension(nnode,*), intent(in) :: conec
  
  integer, dimension(nnode), intent(out) :: econc
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  econc(:)= 0

  ! get element connectivity data
  do inode=1, nnode
     econc(inode)= conec(inode, ielem)
  end do



  return
end subroutine geteledata0





subroutine geteledata1(ielem,ndime,nnode,conec,coord, econc,ecord) 
  !=======================================================================
  !  geteledata= get element data: connectivity, coordinate
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ielem : current element number
  !
  !  ndime : dimension
  !
  !  nnode : node per element
  !
  !  conec(nnode,*) : global element connectivity data
  !
  !  coord(ndime,*) : global nodal coordinate data
  ! 
  !  output:
  !  ------
  !  econc : element connectivity data
  !          econc(inode) : global node number of current i th node
  !
  !  ecord : element nodal coordinate data
  !          ecord(idime, inode) : i direction coordinate of i th node
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ielem
  integer, intent(in) :: ndime, nnode  
  integer, dimension(nnode,*), intent(in) :: conec
  real(8), dimension(ndime,*), intent(in) :: coord
  
  integer, dimension(nnode), intent(out) :: econc
  real(8), dimension(ndime,nnode), intent(out) :: ecord
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: inode, idime
  ! ====================================

  ! initialize
  econc(:)= 0
  ecord(:,:)= 0.d0

  ! get element connectivity data
  do inode=1, nnode
     econc(inode)= conec(inode, ielem)
  end do

  ! get element nodal coordinate
  do inode=1, nnode
     do idime=1, ndime
        ecord(idime, inode)= coord(idime, econc(inode))
     end do
  end do



  return
end subroutine geteledata1





subroutine geteleval1(msize,nnode,nndof,econc,disp, edisp) 
  !=======================================================================
  !  geteleval1= get element nodal displacement value
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  msize :
  !
  !  ndime,nnode : basic dimensions of problem
  !
  !  econc(nnode) : element connectivity data
  !
  !  disp(msize,1) : global displacement vector
  !
  !  output:
  !  ------
  !  edisp(nndof,nnode) : element displacement
  !                       edisp(idof,inode) : i dof displacement of i th node
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: msize
  integer, intent(in) :: nnode, nndof
  integer, dimension(nnode), intent(in) :: econc
  real(8), dimension(msize,1), intent(in) :: disp

  real(8), dimension(nndof,nnode), intent(out) :: edisp
  ! ====================================
  ! local variable
  ! ==============
  integer :: jnode, iloc

  integer :: inode, idof
  ! ====================================

  ! initialize
  edisp(:,:)= 0.0d0


  do inode=1, nnode
     jnode= econc(inode)

     do idof=1, nndof

	     iloc=(jnode-1)*nndof + idof
        edisp(idof,inode)= disp(iloc,1) ! nodal displacement

      end do
  end do



  return
end subroutine geteleval1





subroutine geteleval2(msize,nnode,nndof,econc,disp,velo, edisp,evelo) 
  !=======================================================================
  !  geteleval2= get element nodal displacement, velocity and acceleration
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  msize :
  !
  !  ndime,nnode : basic dimensions of problem
  !
  !  econc(nnode) : element connectivity data
  !
  !  disp(msize,1) : global displacement vector
  !  
  !  velo(msize,1) : gloval velocity vector
  !
  !  output:
  !  ------
  !  edisp(nndof,nnode) : element displacement
  !                       edisp(idof,inode) : i dof displacement of i th node
  !
  !  evelo(nndof,nnode) : element velocity
  !                      evel(idof,inode) : i dof velocity of i th node
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: msize
  integer, intent(in) :: nnode, nndof
  integer, dimension(nnode), intent(in) :: econc
  real(8), dimension(msize,1), intent(in) :: disp
  real(8), dimension(msize,1), intent(in) :: velo

  real(8), dimension(nndof,nnode), intent(out) :: edisp
  real(8), dimension(nndof,nnode), intent(out) :: evelo
  ! ====================================
  ! local variable
  ! ==============
  integer :: jnode, iloc

  integer :: inode, idof
  ! ====================================

  ! initialize
  edisp(:,:)= 0.0d0
  evelo(:,:)= 0.0d0


  do inode=1, nnode
     jnode= econc(inode)

     do idof=1, nndof

	    iloc=(jnode-1)*nndof + idof
        edisp(idof,inode)= disp(iloc,1) ! nodal displacement
        evelo(idof,inode)= velo(iloc,1) ! nodal velocity

      end do
  end do



  return
end subroutine geteleval2





subroutine geteleval3(msize,nnode,nndof,econc,disp,velo,accl, edisp,evelo,eaccl) 
  !=======================================================================
  !  geteleval3= get element nodal displacement, velocity and acceleration
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  msize :
  !
  !  ndime,nnode : basic dimensions of problem
  !
  !  econc(nnode) : element connectivity data
  !
  !  disp(msize,1) : global displacement vector
  !  
  !  velo(msize,1) : gloval velocity vector
  !
  !  accl(msize,1) : global accerelation vector
  !
  !  output:
  !  ------
  !  edisp(nndof,nnode) : element displacement
  !                       edisp(idof,inode) : i dof displacement of i th node
  !
  !  evelo(nndof,nnode) : element velocity
  !                      evel(idof,inode) : i dof velocity of i th node
  !
  !  eaccl(nndof,nnode) : element accleration
  !                      eaccl(idof,inode) : i dof accleration of i th node
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: msize
  integer, intent(in) :: nnode, nndof
  integer, dimension(nnode), intent(in) :: econc
  real(8), dimension(msize,1), intent(in) :: disp
  real(8), dimension(msize,1), intent(in) :: velo
  real(8), dimension(msize,1), intent(in) :: accl

  real(8), dimension(nndof,nnode), intent(out) :: edisp
  real(8), dimension(nndof,nnode), intent(out) :: evelo
  real(8), dimension(nndof,nnode), intent(out) :: eaccl
  ! ====================================
  ! local variable
  ! ==============
  integer :: jnode, iloc

  integer :: inode, idof
  ! ====================================

  ! initialize
  edisp(:,:)= 0.0d0
  evelo(:,:)= 0.0d0
  eaccl(:,:)= 0.0d0


  do inode=1, nnode
     jnode= econc(inode)

     do idof=1, nndof

	    iloc=(jnode-1)*nndof + idof
        edisp(idof,inode)= disp(iloc,1) ! nodal displacement
        evelo(idof,inode)= velo(iloc,1) ! nodal velocity
        eaccl(idof,inode)= accl(iloc,1) ! nodal accleration

      end do

  end do



  return
end subroutine geteleval3





subroutine getmatprop0(ielem,matpro,matnum, ematpro) 
  !=======================================================================
  !  getmatprop0 = get material property
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ielem : current element number
  !
  !  matpro(20,*) : material property data
  !
  !  matnum(*) : assigned material number at each element
  !
  !  output:
  !  ------
  !  ematpro(20,1) : material properties
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ielem
  real(8), dimension(20,*), intent(in) :: matpro
  integer, dimension(*), intent(in) :: matnum

  real(8), dimension(20), intent(out) :: ematpro
  ! ====================================

  ! initialize
  ematpro(:)= 0.0d0


  ! set material property
  ematpro(1:20)= matpro(1:20,matnum(ielem))



  return
end subroutine getmatprop0





subroutine getmatprop1(ielem,matpro,matnum, young,poiss,denst,thick) 
  !=======================================================================
  !  getmatprop1 = get material property: E, nu, rho + thickness
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ielem : current element number
  !
  !  matpro(20,*) : material property data
  !
  !  matnum(*) : assigned material number at each element
  !
  !  output:
  !  ------
  !  young, poiss, denst : material properties
  !
  !  thick : thickness
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ielem
  real(8), dimension(20,*), intent(in) :: matpro
  integer, dimension(*), intent(in) :: matnum

  real(8), intent(out) :: young, poiss, denst, thick
  ! ====================================

  ! initialize
  young= 0.0d0
  poiss= 0.0d0
  denst= 0.0d0
  thick= 0.0d0


  ! set material property
  young= matpro(1,matnum(ielem))
  poiss= matpro(2,matnum(ielem))
  denst= matpro(3,matnum(ielem))

  thick= matpro(20,matnum(ielem))



  return
end subroutine getmatprop1





integer function ismbr1di(nndx,array1d,inumber)
  !=======================================================================
  !  ismember1di = check given inumber is member of array1d
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndx : size of array
  !
  !  array1d(nndx) : 1d integer array
  !
  !  inumber : given number
  !
  !  output:
  !  ------
  !  ismbr1di : result flag
  !             ismbr1di= 0 : inumber is not a member of array1d
  !             ismbr1di= 1 : inumber is a member of array1d
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndx
  integer, dimension(nndx), intent(in) :: array1d
  integer, intent(in) :: inumber

  ! ====================================
  ! local variable
  ! ==============
  
  ! loop indx
  integer :: indx
  ! ====================================

  ! initialize
  ismbr1di= 0


  ! check cracking condition
  do indx=1, nndx

     if( inumber == array1d(indx) ) then
        ismbr1di= 1
        exit
     end if

  end do



  return
end function ismbr1di





integer function ismbr2di(nrndx,ncndx,array2d,inumber)
  !=======================================================================
  !  ismbr2di = check given inumber is member of array2d
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nrndx, ncndx : size of array
  !
  !  array2d(nrndx,ncndx) : 2d integer array
  !
  !  inumber : given number
  !
  !  output:
  !  ------
  !  ismbr2di : result flag
  !             ismbr2di= 0 : inumber is not a member of array2d
  !             ismbr2di= 1 : inumber is a member of array2d
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nrndx, ncndx
  integer, dimension(nrndx,ncndx), intent(in) :: array2d
  integer, intent(in) :: inumber
  ! ====================================
  ! local variable
  ! ==============
  
  ! loop indx
  integer :: irndx, icndx
  ! ====================================

  ! initialize
  ismbr2di= 0

  ! check cracking condition
  do irndx=1, nrndx
     do icndx=1, ncndx  

        if( inumber == array2d(irndx,icndx) ) then
           ismbr2di= 1
           return
        end if

     end do
  end do



  return
end function ismbr2di





integer function ismbr1d1d0(nndx0,nndx1,array1d0,array1d1)
  !=======================================================================
  !  ismbr1d1d0 = check any component of given array1d1 is member of array1d
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndx0,nndx1 : size of array
  !
  !  array1d0(nndx) : 1d integer array
  !
  !  array1d1(nndx) : 1d integer array
  !
  !  output:
  !  ------
  !  ismbr1d1d : result flag
  !              ismbr1d1d= 0 : array1d1 is not a member of array1d
  !              ismbr1d1d= 1 : array1d1 is a member of array1d
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndx0, nndx1
  integer, dimension(nndx0), intent(in) :: array1d0
  integer, dimension(nndx1), intent(in) :: array1d1
  ! ====================================
  ! local variable
  ! ==============
  
  ! loop indx
  integer :: indx, jndx
  ! ====================================

  ! initialize
  ismbr1d1d0= 0

  ! check cracking condition
  do indx=1, nndx0

     do jndx=1, nndx1

        if( array1d0(indx) == array1d1(jndx) ) then
           ismbr1d1d0= 1
           exit
        end if

     end do

  end do



  return
end function ismbr1d1d0





subroutine geteleval1vec(msize,nnode,nndof,econc,disp, edisp) 
  !=======================================================================
  !  geteleval1vec= get element nodal displacement value as a column vector
  !
  !                 note:
  !                 ----
  !                 the returned edip is a column vector
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  msize :
  !
  !  ndime,nnode : basic dimensions of problem
  !
  !  econc(nnode) : element connectivity data
  !
  !  disp(msize,1) : global displacement vector
  !
  !  output:
  !  ------
  !  edisp(nndof*nnode,1) : element displacement
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: msize
  integer, intent(in) :: nnode, nndof
  integer, dimension(nnode), intent(in) :: econc
  real(8), dimension(msize,1), intent(in) :: disp

  real(8), dimension(nndof*nnode,1), intent(out) :: edisp
  ! ====================================
  ! local variable
  ! ==============
  integer :: jnode, iloc, jloc

  integer :: inode, idof
  ! ====================================

  ! initialize
  edisp(:,:)= 0.0d0


  jloc= 0 ! initialize
  do inode=1, nnode
     jnode= econc(inode)

     do idof=1, nndof

	     iloc=(jnode-1)*nndof + idof
        jloc= jloc + 1

        edisp(jloc,1)= disp(iloc,1) ! nodal displacement

      end do
  end do



  return
end subroutine geteleval1vec



subroutine cnvtevalvec2mat(nnode,nndof,evalvec, evalmat) 
  !=======================================================================
  !  cnvtevalvec2mat= convert element wise vector to element wise matrix
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime,nnode : basic dimensions of problem
  !
  !  econc(nnode) : element connectivity data
  !
  !  disp(msize,1) : global displacement vector
  !
  !  evalvec(nndof*nnode,1) : element displacement
  !
  !  output:
  !  ------
  !  evalmat(nndof*nnode,1) : element displacement
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nnode, nndof
  real(8), dimension(nndof*nnode,1), intent(in) :: evalvec

  real(8), dimension(nndof,nnode), intent(out) :: evalmat
  ! ====================================
  ! local variable
  ! ==============
  integer :: iloc

  integer :: inode, idof
  ! ====================================

  ! initialize
  evalmat(:,:)= 0.0d0


  iloc= 0 ! initialize
  do inode=1, nnode
     do idof=1, nndof

	     iloc= iloc + 1
        evalmat(idof,inode)= evalvec(iloc,1) ! nodal displacement

      end do
  end do



  return
end subroutine cnvtevalvec2mat