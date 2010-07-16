! ==========================
! apply boundary condition2 : element surface
! ==========================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine          applybc2             (optbc,optdom,opttrc,optele,msize,ndime,nnode,nndof,ntdof,ngqpt4,ncond, &
!                                               time,tmax,prmtft,elestatus,conec,coord,disp,cond, &
!                                               fext)
! 2.  subroutine          elefbcbt2            (optdom,opttrc,ngaus,ecord,edisp,trac, efbc)
! 3.  subroutine          gqfbc2               (optele,ndime,nnode,gqpsi,gqeta,gqweigt,ecord2d,tracvec, gqfbc)
!
! =========================================================================================================



subroutine applybc2(optbc,optdom,opttrc,optele,msize,ndime,nnode,nndof,ntdof,ngqpt4,ncond, &
                    time,tmax,prmtft,elestatus,conec,coord,disp,cond, &
                    fext)
  !=======================================================================
  !  applybc2 = compute force bc
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
  !  optdom : integration domain option handler
  !           optdom= 0 : integrate over reference domain
  !           optdom= 1 : integrate over current domain
  !
  !  opttrc : traction type option handler
  !           opttrc= 0 : pressure
  !           opttrc= 1 : traction
  !
  !  optele : element option hadler
  !
  !  msize : size of result array
  !
  !  ndime,nndof,ntdof : problem demension data
  !
  !  ngqpt4 : integration rule
  !
  !  ncond : the number of transient bc nodes 
  !
  !  time, tmax : current and maximum simulation time
  !
  !  prmtft(*) : time function parameters
  !
  !  elestatus(*) : element status
  !
  !  conec(nnode,*) : connectivity matrix
  !
  !  coord(ndime,*) : nodal coordinate
  !
  !  disp(msize,1) : displacement
  !
  !  cond(structure)
  !  ---------------
  !     %ele : element number
  !     %nod(2) : boundary condition node
  !     %dof(6) : applied dof
  !     %val(6) : prescribed coundary condition value
  !     %tmf(6) : time function type
  !
  !  inoutput:
  !  --------
  !  fext(msize,1) : external force matrix
  !
  ! ======================================================================

  use typcond
  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optbc,optdom, opttrc
  integer, intent(in) :: optele
  integer, intent(in) :: msize,ndime,nnode,nndof,ntdof,ngqpt4,ncond
  real(8), intent(in) :: time, tmax
  real(8), dimension(*), intent(in) :: prmtft
  integer, dimension(*), intent(in) :: elestatus
  integer, dimension(nnode,*), intent(in) :: conec
  real(8), dimension(ndime,*), intent(in) :: coord
  real(8), dimension(msize,1), intent(in) :: disp
  type(cond_bc), dimension(*), intent(in) :: cond

  real(8), dimension(msize,1), intent(inout) :: fext
  ! ====================================
  ! local variable
  ! ==============
  integer :: ielem
  real(8), dimension(ndime,1) :: trac

  integer, dimension(nnode) :: econc
  real(8), dimension(ndime,nnode) :: ecord
  real(8), dimension(nndof,nnode) :: edisp

  real(8), dimension(nndof*nnode,1) :: efbc
  integer, dimension(nndof*nnode) :: sctr

  integer :: itmft, iloc, nflagrmv
  real(8) :: xdist
  real(8) :: tmftval

  ! loop index
  integer :: icond, idime, idof, inode
  ! ====================================

  ! initialize: do not initialize


  ! skip 2d element: tri. and quad.
  if ( optele == 1 .or. optele == 2 ) return

  do icond=1, ncond

     ! get element number
     ielem= cond(icond)%ele

     ! loop over: uncracked element, no cracking element
     if ( elestatus(ielem) == 0 .or. elestatus(ielem) == 1 ) then

        ! --------------------------------------------------------------------------------------------

        ! get element connectivity and coordinate
        call geteledata1(ielem,ndime,nnode,conec,coord, econc,ecord)
           ! input : ielem,ndime,nnode,conec,coord
           ! output : econc,ecord

        ! get element nodal displacement
        call geteleval1(msize,nnode,nndof,econc,disp, edisp)
           ! input : msize,nnode,nndof,econc,disp
           ! output : edisp

! ### working
xdist= 0.0d0 ! initialzie
do inode=1, nnode
xdist= xdist + ecord(1,inode)
end do
xdist= xdist / dble(nnode)

        ! get traction
        select case(opttrc)
        case(0) ! pressure
           trac(:,:)= 0.0d0
           trac(1,1)= cond(icond)%val(1)

        case(1) ! traction
           do idime=1, ndime
              trac(idime,1)= cond(icond)%val(idime)
           end do

        end select

        ! --------------------------------------------------------------------------------------------
        ! compute element external force
        ! ------------------------------
        select case(optele)
        case(3) ! bt shell
           call elefbcbt2(optdom,opttrc,nndof,ngqpt4,ecord,edisp,trac, &
                          efbc)
              ! input : optdom,opttrc,nndof,ngqpt4,ecord,edisp,trac
              ! output : efbc

        end select

        ! get scatter matrix
        call getsctr(nnode,nndof,econc, sctr)
           ! input : nnode,nndof,econc
           ! output : sctr

        ! --------------------------------------------------------------------------------------------
        ! set result
        ! ----------
        do idof=1, ntdof

           ! skip
           if( opttrc==1 .and. cond(icond)%dof(idof)==0 ) cycle

           ! time function factor
           if ( optbc== 1 ) then ! transient bc

              if ( opttrc == 0 ) then  ! pressure
                 itmft= cond(icond)%tmf(1)

              else if ( opttrc==1 ) then ! traction
                 itmft= cond(icond)%tmf(idof)

              end if

              ! get time function
              call gettmftval(itmft,prmtft,time,tmax,xdist, tmftval,nflagrmv)
                 ! input : itmft,prmtft,time,tmax,xdist
                 ! output : tmftval,nflagrmv

           else ! initial bc
              tmftval= 1.0d0

           end if

           ! set result
           do inode=1, nnode

             ! location
             iloc= nndof*(inode-1) + idof

             ! add element fbc to fext
             fext(sctr(iloc),1)= fext(sctr(iloc),1) + efbc(iloc,1)  * tmftval

          end do

        end do ! end do idof
        ! --------------------------------------------------------------------------------------------

     end if ! if ( elestatus(ielem) == 0 .or. elestatus(ielem) == 1 )

  end do ! end do icond



  return
end subroutine applybc2





subroutine elefbcbt2(optdom,opttrc,nndof,ngaus,ecord,edisp,trac, efbc)
  !=======================================================================
  !  elefbcbt2 = compute surface force for belytschko-tsay shell element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optdom : integration domain option handler
  !           optdom= 0 : integrate over reference domain
  !           optdom= 1 : integrate over current domain
  !
  !  opttrc : traction type option handler
  !           opttrc= 0 : pressure
  !           opttrc= 1 : traction
  !
  !  nndof : the total number of dof per node
  !
  !  ngaus : integration rule
  !
  !  ecord(3,4) : element nodal coordinate
  !
  !  edisp(nndof,4) : nodal displacement of original element
  !
  !  trac(3,1) : traction vector
  !              for pressure case, specify trac(1,1) with hydro-static pressure
  !  output:
  !  ------
  !  efbc(nndof*4,1) : bounday force
  !                          
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optdom, opttrc
  integer, intent(in) :: nndof, ngaus
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(nndof,4), intent(in) :: edisp
  real(8), dimension(3,1), intent(in) :: trac

  real(8), dimension(nndof*4,1), intent(out) :: efbc
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,4) :: ecurn
  real(8), dimension(3,3) :: locbvec
  real(8), dimension(3,4) :: ecordloc
  real(8), dimension(2,4) :: ecord2dloc

  real(8), dimension(3,1) :: norlvecloc, norlvec
  real(8), dimension(3,1) :: tracvec

  real(8), dimension(2,ngaus**2) :: gqpoin 
  real(8), dimension(ngaus**2) :: gqweigt 

  real(8), dimension(12,1) :: efbc0, gqfbc

  ! loop index
  integer :: igaus
  ! ====================================

  ! initialize
  efbc(:,:)= 0.0d0

  ! -------------------------------------------------------------------------------

  ! integration domain
  select case(optdom)
  case(0) ! reference domain
     ecurn(1:3,1:4)= ecord(1:3,1:4)

  case(1) ! current domain
     ecurn(1:3,1:4)= ecord(1:3,1:4) + edisp(1:3,1:4)

  end select

  ! compute co rotational local base vector: locbvec
  call getlocbvecbt(ecurn, locbvec)
     ! input : ecurn
     ! output : locbvec

  ! get local nodal coordinate: glb -> loc
  call glb2locnodv(3,4,3,locbvec,ecurn, ecordloc)
     ! input : 3(ndime),4(nnode),3(ntrndof),locbvec,ecurn
     ! output : ecordloc

  ecord2dloc(1:2,1:4)= ecordloc(1:2,1:4)

  ! -------------------------------------------------------------------------

  ! traction vector
  select case(opttrc)
  case(0) ! pressure

     ! set normal to element surface vector
     norlvecloc(1:2,1)= 0.0d0
     norlvecloc(3,1)= -1.0d0

     ! convert local to global
     call loc2glbv(3,3,locbvec,norlvecloc, norlvec)
        ! input : 3(ndime),3(ntrndof),zvecloc,norlvecloc
        ! output : norlvec

     ! set pressure vector
     tracvec(1:3,1)= trac(1,1) * norlvec(1:3,1)

  case(1) ! traction

     tracvec(1:3,1)= trac(1:3,1)

  end select

  ! -------------------------------------------------------------------------

  ! get 2d regular gq rule
  call getgqele(3,2,ngaus,ngaus**2, gqpoin,gqweigt)
     ! input : 3(optele),2(ngqdim),ngaus,ngaus**2
     ! output : gqpoin,gqweigt

  ! -------------------------------------------------------------------------

  efbc0(:,:)= 0.0d0 ! initialize
  do igaus=1, ngaus*ngaus

     ! compute cohesive force at gq
     call gqfbc2(3,3,4,gqpoin(1,igaus),gqpoin(2,igaus),gqweigt(igaus),ecord2dloc,tracvec, &
                 gqfbc)
        ! input : 3(optele),3(ndime),4(nnode),gqpsi,gqeta,gqweigt,ecord2dloc,tracvec
        ! output : gqfbc

     ! sum on cohesive force
     efbc0(1:12,1)= efbc0(1:12,1) + gqfbc(1:12,1)

  end do

  ! -------------------------------------------------------------------------
  ! set results
  ! -----------
  if ( nndof == 5 ) then
     efbc(1:3,1)= efbc0(1:3,1)     ! f_x, f_y and f_z at node 1
     efbc(6:8,1)= efbc0(4:6,1)     ! f_x, f_y and f_z at node 2
     efbc(11:13,1)= efbc0(7:9,1)   ! f_x, f_y and f_z at node 3
     efbc(16:18,1)= efbc0(10:12,1) ! f_x, f_y and f_z at node 4

  else if ( nndof == 6 ) then
     efbc(1:3,1)= efbc0(1:3,1)     ! f_x, f_y and f_z at node 1
     efbc(7:9,1)= efbc0(4:6,1)     ! f_x, f_y and f_z at node 2
     efbc(13:15,1)= efbc0(7:9,1)   ! f_x, f_y and f_z at node 3
     efbc(19:21,1)= efbc0(10:12,1) ! f_x, f_y and f_z at node 4

  end if



  return
end subroutine elefbcbt2





subroutine gqfbc2(optele,ndime,nnode,gqpsi,gqeta,gqweigt,ecord2d,tracvec, gqfbc)
  !=======================================================================
  !  gqfbc2 = compute external force matrix due to stress or pressure
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optele : option hadler
  !
  !  ndime : dimension
  !
  !  nnode : the total number of node per element
  !
  !  gqpsi, gqeta : integration point in parent domain
  !
  !  gqweigt : gq weight function value
  !
  !  ecord2d(2,nnode) : 2d nodal coordinate
  !
  !  tracvec(ndime,1) : traction vector
  !
  !  output:
  !  ------
  !  gqfbc(ndime*nnode,1) : external force at current gq
  !                      
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele
  integer, intent(in) :: ndime, nnode
  real(8), intent(in) :: gqpsi, gqeta, gqweigt
  real(8), dimension(2,nnode) :: ecord2d
  real(8), dimension(ndime,1), intent(in) :: tracvec

  real(8), dimension(ndime*nnode,1), intent(out) :: gqfbc
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nnode) :: shap
  real(8), dimension(2,nnode) :: deriv
  real(8) :: djacob

  real(8), dimension(ndime,ndime*nnode) :: nmat

  real(8), dimension(ndime*nnode,1) :: temp
  ! ====================================

  ! initialize
  gqfbc(:,:)= 0.0d0

  ! compute 2d shape function
  call getshape2d(optele,nnode,gqpsi,gqeta, shap,deriv)
     ! input : optele,nnode,gqpsi,gqeta
     ! output : shap,deriv

  ! compute determinant of jacobian
  call jacob0(2,nnode,deriv,ecord2d, djacob)
     ! input : 2(ndime),nnode,deriv,ecord2d
     ! output : djacb

  ! get n matrix
  call getnmat(nnode,ndime,shap, nmat)
     ! input : nnode,ndime,shap
     ! output : nmat

  ! ---------------------------------------------------
  ! compute nmat^t.tracvec
  call matprd(ndime,ndime*nnode,1, ndime,1,0, ndime*nnode,1, nmat,tracvec, temp)
     ! input : ndime,ndime*nnode,1, ndime,1,0, ndime*nnode,1, nmat,normlvec
     ! output : temp

  ! compute force vector
  gqfbc(1:ndime*nnode,1)= temp(1:ndime*nnode,1) * gqweigt * djacob



  return
end subroutine gqfbc2