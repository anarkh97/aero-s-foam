! ==================================
! fem data manipulation  
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine          getdofs                    (optele,npoin,nelem,mnode,eletyp,conec, dofs)
! 2.  subroutine          getedofs                   (ielem,mnode,nnode,dofs,conec, edofs)
! 3.  subroutine          getadrs                    (npoin,dofs, size,adrs)
! 4.  subroutine          geteadrs                   (nnode,adrs,econc, eadrs)
! 5.  subroutine          getsctrhybrid              (nnode,nndof,eadrs, sctr)
! 6.  subroutine          geteledata0hybrid          (ielem,mnode,nnode,conec, econc)
! 7.  subroutine          geteledata1hybrid          (ielem,mnode,ndime,nnode,conec,coord, econc,ecord)
! 8.  subroutine          geteleval1hybrid           (msize,nnode,nndof,eadrs,disp, edisp)
! 9.  subroutine          geteleval2hybrid           (msize,nnode,nndof,eadrs,disp,velo, edisp,evelo)
! 10. subroutine          geteleval3hybrid           (msize,nnode,nndof,eadrs,disp,velo,accl, edisp,evelo,eaccl)
! 11. subroutine          getctvmodel                (ielem,matnum,ctvtyp, optctv,optpty,optdmg)
!
! =========================================================================================================



subroutine getdofs(optele,npoin,nelem,mnode,eletyp,conec, dofs)
  !=======================================================================
  !  getdofs =get dofs
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !                           
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele
  integer, intent(in) :: npoin,nelem,mnode
  integer, dimension(2,*), intent(in) :: eletyp
  integer, dimension(mnode,*), intent(in) :: conec

  integer, dimension(2,npoin), intent(out) :: dofs
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: ielem, inode
  ! ====================================
  
  ! initialize
  dofs(:,:)= 0

  
  select case(optele)
  case(1,2) ! 2d tri/quad
     dofs(1,1:npoin)= 2 ! ntdof
     dofs(2,1:npoin)= 0 ! nrdof

  case(3) ! bt shell
     dofs(1,1:npoin)= 3 ! ntdof
     dofs(2,1:npoin)= 3 ! nrdof
  
  case(4) ! 3d brick
     dofs(1,1:npoin)= 3 ! ntdof
     dofs(2,1:npoin)= 0 ! nrdof

  case(5) ! bt shell + 3d brick
     do ielem=1, nelem

        if ( eletyp(1,ielem) == 3 ) then ! bt shell
           do inode=1, 4
              dofs(1,conec(inode,ielem))= 3
              dofs(2,conec(inode,ielem))= 3
           end do
        
        else if  ( eletyp(1,ielem) == 4 ) then ! 3d brick
           do inode=1, 8
              dofs(1,conec(inode,ielem))= 3
              dofs(2,conec(inode,ielem))= 0
           end do
     
        else
           write(*,*) "not implemented eletyp: getdofs"
           write(nout5,*) "not implemented eletyp: getdofs"
        
        end if

     end do

  case(6) ! t3 + q4
     dofs(1,1:npoin)= 2 ! ntdof
     dofs(2,1:npoin)= 0 ! nrdof

  case default
     write(*,*) "not implemented optele: getdofs"
     write(nout5,*) "not implemented optele: getdofs"
  
  end select



  return
end subroutine getdofs





subroutine getedofs(ielem,mnode,nnode,dofs,conec, edofs)
  !=======================================================================
  !  getdofs =get dofs
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !                           
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ielem
  integer, intent(in) :: mnode, nnode
  integer, dimension(2,*), intent(in) :: dofs
  integer, dimension(mnode,*), intent(in) :: conec

  integer, dimension(2), intent(out) :: edofs
  ! ====================================
  ! local variable
  ! ==============
  integer, dimension(nnode) :: econc

  ! loop index
  integer :: inode
  ! ====================================
  
  ! initialize
  edofs(:)= 0

  
  ! get element nodal connectivity
  call geteledata0hybrid(ielem,mnode,nnode,conec, econc)
     ! input : ielem,mnode,nnode,conec
     ! output : econc
  
  do inode=1, nnode

     ! translational dof
     edofs(1)= edofs(1) + dofs(1,econc(inode))

     ! rotational dof
     edofs(2)= edofs(2) + dofs(2,econc(inode))
  
  end do

  ! translational dof
  edofs(1)= edofs(1) / nnode

  ! rotational dof
  edofs(2)= edofs(2) / nnode



  return
end subroutine getedofs





subroutine getadrs(npoin,dofs, size,adrs)
  !=======================================================================
  !  getadrs =get global nodal addresses
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !                           
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: npoin
  integer, dimension(2,npoin), intent(in) :: dofs

  integer, intent(out) :: size
  integer, dimension(npoin), intent(out) :: adrs
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: ipoin
  ! ====================================
  
  ! initialize
  adrs(:)= 0
  
  
  adrs(1)= 1 ! initialize
  size= dofs(1,1) + dofs(2,1)

  do ipoin=2, npoin

     adrs(ipoin)= adrs(ipoin-1) +  dofs(1,ipoin-1) + dofs(2,ipoin-1)

     size= size + dofs(1,ipoin) + dofs(2,ipoin)

  end do



  return
end subroutine getadrs





subroutine geteadrs(nnode,adrs,econc, eadrs)
  !=======================================================================
  !  getadrs =get global nodal addresses
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !                           
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nnode
  integer, dimension(*), intent(in) :: adrs
  integer, dimension(nnode) :: econc

  integer, dimension(nnode), intent(out) :: eadrs
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: inode
  ! ====================================
  
  ! initialize
  eadrs(:)= 0
  

  do inode=1, nnode
     eadrs(inode)= adrs(econc(inode))
  end do



  return
end subroutine geteadrs





subroutine getsctrhybrid(nnode,nndof,eadrs, sctr)
  !=======================================================================
  !  getsctrhybrid= get scatter matrix for hybrid mesh
  !                 : convert local node number to system equation address
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
  integer, dimension(*), intent(in) :: eadrs

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
	    
	    if( eadrs(inode) /= 0 ) then

          sctr(indx)= eadrs(inode) + idof - 1

       end if

     end do
  end do

  return
end subroutine getsctrhybrid





subroutine geteledata0hybrid(ielem,mnode,nnode,conec, econc)
  !=======================================================================
  !  geteledata0hybrid= get element data: connectivity for hybrid element mesh
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ielem : current element number
  !
  !  mnode : maximum node per element in hybrid element mesh
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
  integer, intent(in) :: mnode
  integer, intent(in) :: nnode
  integer, dimension(mnode,*), intent(in) :: conec
  
  integer, dimension(nnode), intent(out) :: econc
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: inode, idime
  ! ====================================

  ! initialize
  econc(:)= 0


  ! get element connectivity data
  do inode=1, nnode
     econc(inode)= conec(inode,ielem)
  end do



  return
end subroutine geteledata0hybrid





subroutine geteledata1hybrid(ielem,mnode,ndime,nnode,conec,coord, econc,ecord)
  !=======================================================================
  !  geteledata1hybrid= get element data: connectivity, coordinate for hybrid element mesh
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ielem : current element number
  !
  !  mnode : maximum node per element in hybrid element mesh
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
  integer, intent(in) :: mnode
  integer, intent(in) :: ndime, nnode  
  integer, dimension(mnode,*), intent(in) :: conec
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
end subroutine geteledata1hybrid





subroutine geteleval1hybrid(msize,nnode,nndof,eadrs,disp, edisp) 
  !=======================================================================
  !  geteleval1hybrid= get element nodal displacement value
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  msize :
  !
  !  ndime,nnode : basic dimensions of problem
  !
  !  eadrs(nnode) : element nodal global address data
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
  integer, dimension(nnode), intent(in) :: eadrs
  real(8), dimension(msize,1), intent(in) :: disp

  real(8), dimension(nndof,nnode), intent(out) :: edisp
  ! ====================================
  ! local variable
  ! ==============
  integer :: iloc

  integer :: inode, idof
  ! ====================================

  ! initialize
  edisp(:,:)= 0.0d0


  do inode=1, nnode

     if ( eadrs(inode) /= 0 ) then

        do idof=1, nndof

           iloc= eadrs(inode) + idof -1

           edisp(idof,inode)= disp(iloc,1) ! nodal displacement

        end do ! do idof=1

     end if ! if ( eadrs(inode) /= 0 )

  end do ! do inode=1



  return
end subroutine geteleval1hybrid





subroutine geteleval2hybrid(msize,nnode,nndof,eadrs,disp,velo, edisp,evelo) 
  !=======================================================================
  !  geteleval2hybrid= get element nodal displacement, velocity and acceleration
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  msize :
  !
  !  ndime,nnode : basic dimensions of problem
  !
  !  eadrs(nnode) : element nodal global address data
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
  integer, dimension(nnode), intent(in) :: eadrs
  real(8), dimension(msize,1), intent(in) :: disp
  real(8), dimension(msize,1), intent(in) :: velo

  real(8), dimension(nndof,nnode), intent(out) :: edisp
  real(8), dimension(nndof,nnode), intent(out) :: evelo
  ! ====================================
  ! local variable
  ! ==============
  integer :: iloc

  integer :: inode, idof
  ! ====================================

  ! initialize
  edisp(:,:)= 0.0d0
  evelo(:,:)= 0.0d0


  do inode=1, nnode

     if ( eadrs(inode) /= 0 ) then

        do idof=1, nndof

           iloc= eadrs(inode) + idof - 1
          
           edisp(idof,inode)= disp(iloc,1) ! nodal displacement
           evelo(idof,inode)= velo(iloc,1) ! nodal velocity

        end do ! do idof=1

     end if ! if ( eadrs(inode) /= 0 )

  end do ! do inode=1



  return
end subroutine geteleval2hybrid





subroutine geteleval3hybrid(msize,nnode,nndof,eadrs,disp,velo,accl, edisp,evelo,eaccl) 
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
  !  eadrs(nnode) : element nodal global address data
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
  integer, dimension(nnode), intent(in) :: eadrs
  real(8), dimension(msize,1), intent(in) :: disp
  real(8), dimension(msize,1), intent(in) :: velo
  real(8), dimension(msize,1), intent(in) :: accl

  real(8), dimension(nndof,nnode), intent(out) :: edisp
  real(8), dimension(nndof,nnode), intent(out) :: evelo
  real(8), dimension(nndof,nnode), intent(out) :: eaccl
  ! ====================================
  ! local variable
  ! ==============
  integer :: iloc

  integer :: inode, idof
  ! ====================================

  ! initialize
  edisp(:,:)= 0.0d0
  evelo(:,:)= 0.0d0
  eaccl(:,:)= 0.0d0


  do inode=1, nnode

     if ( eadrs(inode) /= 0 ) then

        do idof=1, nndof

           iloc= eadrs(inode) + idof - 1

           edisp(idof,inode)= disp(iloc,1) ! nodal displacement
           evelo(idof,inode)= velo(iloc,1) ! nodal velocity
           eaccl(idof,inode)= accl(iloc,1) ! nodal accleration

        end do ! do idof=1

     end if ! if ( eadrs(inode) /= 0 )

  end do ! do inode=1



  return
end subroutine geteleval3hybrid





subroutine getctvmodel(ielem,matnum,ctvtyp, optctv,optpty,optdmg) 
  !=======================================================================
  !  getctvmodel = get constitutive model data
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ielem : current element number
  !
  !  matnum(*) : assigned material number at each element
  !
  !  ctvtyp(10,*) : constitutive model data
  !
  !  output:
  !  ------
  !  optctv : constitutive model type
  !
  !  optpty : 2d condition
  !
  !  optdmg : material damage model
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ielem
  integer, dimension(*), intent(in) :: matnum
  integer, dimension(10,*), intent(in) :: ctvtyp

  integer, intent(out) :: optctv, optpty, optdmg
  ! ====================================

  ! initialize
  optpty= 0
  optctv= 0
  optdmg= 0


  ! set constitutive model flags

  ! constitutive model
  ! ------------------
  ! 1 = hypoelstic
  ! 2 = tevp1: thermo-elasto-visco-plastic
  ! 3 = evp1: elasto=visco-plastic
  ! 4 = gtn: gurson-tvergaard-needleman
  ! 5 = j2exp
  ! 6 = hypoelstic_w/_damage1 (lematire)
  ! 7 = hypoelstic_w/_damage2 (linear)
  ! 8 = visco hyperelastic1 (ogden)
  ! 9 = thermo-johnson-cook

  optctv= ctvtyp(1,matnum(ielem))

  ! 2d condition
  ! ------------
  ! 1 = plane-stress
  ! 2 = plane-strain
  optpty= ctvtyp(2,matnum(ielem))

  ! material damage option
  ! ----------------------
  optdmg= 0 ! initialize
  select case(optctv)
  case(6) ! hypoelstic_w/_damage1 (lematire)
     optctv= 1 ! switch to hypoelasticity
     optdmg= 1 ! with lemaitre damage
  
  case(7) ! hypoelstic_w/_damage1 (linear)
     optctv= 1 ! switch to hypoelasticity
     optdmg= 2 ! with linear damage
  
  end select



  return
end subroutine getctvmodel