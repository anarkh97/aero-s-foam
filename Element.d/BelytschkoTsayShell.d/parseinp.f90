! ==================================
! interprete reading input data 
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine          getsimdat                 (prbsiz,prbopt,prbcon, &
!                                                    npoin,nelem,nmatp,nxele,nxtip,ndime,nnode,ntdof,nrdof,ngqpt,ngqpt4,nibcd,ntbcd,ncrcd, &
!                                                    optmhd,optsim,optele,optpty,optctv,optdmg,optcoh,opthgc,optdmp,optcri,optprn)
! 2.  subroutine          getsimprm                 (prbprm, tmax,tmfct,prmdmg,prmcoh,prmhgc,prmcri,prmtft)
! 3.  subroutine          getinistatus              (optmhd,nelem,ncrcd,crctrl, elestatus)
! 4.  subroutine          getgqsize                 (optmhd,optele,ngqpt, mgaus,mgqpt)
! ------------------------------------------------------------------------------------- hybrid mesh
! 5.  subroutine          chkeledat                 (optele,nelem, eletyp)
! 6.  subroutine          getsimdat1                (prbsiz,prbopt,prbcon, &
!                                                    npoin,nelem,nmatp,nxele,nxtip,ndime,nnode,ntdof,nrdof,ngqpt,ngqpt4,nibcd,ntbcd,ncrcd, &
!                                                    optmhd,optsim,optele,optcnt,optcoh,opthgc,optdmp,optcri,optprn)
! 7.  subroutine          getsimprm1                (prbprm, tmax,tmfct,prmcnt,prmcoh,prmhgc,prmcri,prmtft,prmdmp)
! 8.  subroutine          getgqsize1                (optmhd,nelem,eletyp,ngqpt, mgaus,mgqpt)
! 9.  subroutine          parsctvinp                (line, optctv) 
!
! =========================================================================================================



subroutine getsimdat(prbsiz,prbopt,prbcon, &
                     npoin,nelem,nmatp,nxele,nxtip,ndime,nnode,ntdof,nrdof,ngqpt,ngqpt4,nibcd,ntbcd,ncrcd, &
                     optmhd,optsim,optele,optpty,optctv,optdmg,optcoh,opthgc,optdmp,optcri,optprn)
  !=======================================================================
  !  getprobinf = get problem information
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  prbsiz(10) : problem size data
  !
  !  prbopt(10,10) : problem option data
  !
  !  prbcon(3,10) : problem condition data
  !
  !  output:
  !  ------
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, dimension(10), intent(in) :: prbsiz
  integer, dimension(10,10), intent(in) :: prbopt
  integer, dimension(3,10), intent(in) :: prbcon

  integer, intent(out) :: npoin, nelem, nmatp, nxele, nxtip
  integer, intent(out) :: ndime, nnode, ntdof, nrdof
  integer, dimension(3), intent(out) :: ngqpt
  integer, intent(out) :: ngqpt4
  integer, dimension(10), intent(out) :: nibcd, ntbcd, ncrcd
  integer, intent(out) :: optmhd, optsim, optele, optpty, optctv, optdmg, optcoh, opthgc, optdmp
  integer, dimension(10), intent(out) :: optcri, optprn
  ! ====================================

  ! initialize
  ngqpt(:)= 0
  nibcd(:)= 0
  ntbcd(:)= 0
  ncrcd(:)= 0
  optcri(:)= 0
  optprn(:)= 0


  ! problem size
  npoin= prbsiz(1) ! node
  nelem= prbsiz(2) ! element
  nmatp= prbsiz(3) ! assigned material properties set
  nxele= prbsiz(4) ! max extra elements 
  nxtip= prbsiz(5) ! max extra tips

  ! 1. general
  optmhd= prbopt(1,1) ! numerical method
  ndime= prbopt(1,2) ! dimension
  optsim= prbopt(1,3) ! simulation type

  ! 2. element
  optele= prbopt(2,1) ! element
  nnode= prbopt(2,2) ! node per element
  ntdof= prbopt(2,3) ! trans. dof per node
  nrdof= prbopt(2,4) ! rot. dof per node

  ngqpt(1)= prbopt(2,5) ! gq rule for regular domain
  ngqpt(2)= prbopt(2,6) ! gq rule for sub domain
  ngqpt(3)= prbopt(2,7) ! gq rule for through thickness

  ngqpt4= prbopt(2,8) ! gq rule for bc or cohesive force integeation


  ! 3. constitutive law
  optpty= prbopt(3,1) ! problem type
  optctv= prbopt(3,2) ! constitutive

  ! 4. damage model
  optdmg= prbopt(4,1) ! damage model type

  ! 5. cohesive model
  optcoh= prbopt(5,1) ! cohesive model

  ! 6. hg control
  opthgc= prbopt(6,1) ! hourglass control

  ! 7. crack criterion
  optcri(1)= prbopt(7,1) ! criterion type
  optcri(2)= prbopt(7,2) ! criterion average

  ! 8. time function
  
  ! 9. output 
  optprn(1)= prbopt(9,1) ! result print interval
  optprn(2)= prbopt(9,2) ! crack path file
  optprn(3)= prbopt(9,3) ! log echo message file
  optprn(4)= prbopt(9,4) ! echo monitor message
  optprn(5)= prbopt(9,5) ! nodal displacement
  optprn(6)= prbopt(9,6) ! nodal velocity
  optprn(7)= prbopt(9,7) ! stress
  optprn(8)= prbopt(9,8) ! load-deflection curve

  ! 10. damping
  optdmp= prbopt(10,1) ! damping type

  ! condition
  nibcd(1)= prbcon(1,1) ! ic vel
  nibcd(2)= prbcon(1,2) ! ic frc
  nibcd(3)= prbcon(1,3) !
  nibcd(4)= prbcon(1,4) !
  nibcd(5)= prbcon(1,5) !
  nibcd(6)= prbcon(1,6) !
  ! ---------------
  ntbcd(1)= prbcon(2,1) ! bc vel
  ntbcd(2)= prbcon(2,2) ! bc frc
  ntbcd(3)= prbcon(2,3) ! bc prs
  ntbcd(4)= prbcon(2,4) !
  ntbcd(5)= prbcon(2,5) !
  ntbcd(6)= prbcon(2,6) !
  ! --------------

  ncrcd(1)= prbcon(3,1) ! total number of cracking control
  ncrcd(2)= prbcon(3,2) ! total number of load deflection curve point



  return
end subroutine getsimdat





subroutine getsimprm(prbprm, tmax,tmfct,prmdmg,prmcoh,prmhgc,prmcri,prmtft,prmdmp)
  !=======================================================================
  !  getprobprm1 = get general parameter
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  prbprm(10,10) : problem paramters
  !
  !  output:
  !  ------
  !  tmax : maximum simulation time
  !
  !  tmfct : time step control factor
  !
  !  prmdmg(10) : damage model parameters
  !
  !  prmcoh(10) : cohesive crack model parameters
  !
  !  prmhgc(10) : hourglass control parameters
  !
  !  prmcri(10) : crack criterion parameters
  !
  !  prmtft(10) : loading time function parameters
  !
  !  prmdmp(10) : damping parameter
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(10,10), intent(in) :: prbprm

  real(8), intent(out) :: tmax,tmfct
  real(8), dimension(10), intent(out) :: prmdmg
  real(8), dimension(10), intent(out) :: prmcoh
  real(8), dimension(10), intent(out) :: prmhgc
  real(8), dimension(10), intent(out) :: prmcri
  real(8), dimension(10), intent(out) :: prmtft
  real(8), dimension(10), intent(out) :: prmdmp
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: indx
  ! ====================================

  ! simulation parameter
  ! --------------------
  tmax= prbprm(1,1) ! simulatio time
  tmfct= prbprm(1,2)
  ! ------------------
  toler(1)= prbprm(1,3) ! tolerance
  toler(2)= prbprm(1,4) ! <- global variables 
  toler(3)= prbprm(1,5) 
  ! -----------------

  do indx=1, 10
     prmdmg(indx)= prbprm(4,indx)
     prmcoh(indx)= prbprm(5,indx)
     prmhgc(indx)= prbprm(6,indx)
     prmcri(indx)= prbprm(7,indx)
     prmtft(indx)= prbprm(8,indx)
     prmdmp(indx)= prbprm(10,indx)
  end do



  return
end subroutine getsimprm





subroutine getinistatus(opt,optmhd,nelem,ncrcd,crctrl, status)
  !=======================================================================
  !  getinistatus = get initial element status data
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt : option hadler
  !        opt=0 : check element cracking status
  !        opt=1 : check element dimensional coupling status
  !
  !  optmhd : simulation option
  !
  !  nelem : the total number of element
  !
  !  ncrcd(*) : the total number of cracking controled element
  !
  !  crctrl(3,*) : cracking control conditions
  !
  !  output:
  !  ------
  !  status(nelem) : element status flag
  !                  opt=0
  !                  status= -3 : initial crack tip element
  !                  status= -2 : initial crack face element
  !                  status= -1 : killed element
  !                  status=  0 : regular element
  !                  status=  1 : no cracking element
  !                  status=  2 : phantom crack face element
  !                  status=  3 : phantom crack tip element
  !
  !                  opt=1
  !                  status= 11 : dimensional adaptive coupling
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opt
  integer, intent(in) :: optmhd
  integer, intent(in) :: nelem
  integer, dimension(*), intent(in) :: ncrcd
  integer, dimension(3,*), intent(in) :: crctrl

  integer, dimension(nelem), intent(out) :: status
  ! ====================================
  ! local variable
  ! ==============
  integer :: jelem, jcond

  ! loop index
  integer :: icond
  ! ====================================

  ! initialize
  status(1:nelem)= 0


  if( opt==0 ) then ! check element status regarding cracking

     ! loop over whole cracking controled element
     do icond= 1, ncrcd(1)

        jelem= crctrl(1,icond) ! element number
        jcond= crctrl(2,icond) ! condition

        if ( jcond == -1 .or. jcond == 1 ) then
           status(jelem)= jcond

        else if ( jcond == 2 .or. jcond == 3) then
           status(jelem)= -jcond
           ! note: (-) is used to represent initial crack: see, coheisve force
        end if

     end do

  else if ( opt==1 ) then ! check element status regarding dimensional coupling

     ! loop over whole cracking controled element
     do icond= 1, ncrcd(1)

        jelem= crctrl(1,icond) ! element number
        jcond= crctrl(2,icond) ! condition

        select case(optmhd)
        case(7) ! dimensional coupling
           if ( jcond == 11 ) status(jelem)= jcond

        end select

     end do

  end if



  return
end subroutine getinistatus





subroutine getgqsize(optmhd,optele,ngqpt, mgaus,mgqpt)
  !=======================================================================
  !  getgqsize = get general parameter
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optmhd : numerical method option
  !
  !  optele : element type option
  !
  !  ngqpt(3) : gq rule
  !             ngqpt(1) : gq rule for regular element 
  !             ngqpt(2) : gq rule for enriched element 
  !             ngqpt(3) : gq rule through thickness
  !
  !  output:
  !  ------
  !  mgaus(3) : the number of quadrature point
  !             mgaus(1)= regular element
  !             mgaus(2)= enriched element
  !             mgaus(3)= through thickness
  !
  !  mgqpt(2) : the total number of quadrature point
  !             mgqpt(1)= mgaus(1) * mgaus(3) : regular element
  !             mgqpt(2)= mgaus(2) * mgaus(3) : enriched element
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optmhd, optele
  integer, dimension(3), intent(in) :: ngqpt
  
  integer, dimension(3), intent(out) :: mgaus
  integer, dimension(2), intent(out) :: mgqpt
  ! ====================================

  ! initialize
  mgaus(:)= 0
  mgqpt(:)= 0


  ! --------------------------------------------------------------------
  ! set mgaus(1:3)
  ! --------------

  ! in plane regular element: ngqpt(1)
  ! ------------------------
  if( optele==0 ) then ! 1d line
     mgaus(1)= ngqpt(1)

  else if( optele==1 ) then ! 2d tri
     mgaus(1)= ngqpt(1)

  else if( optele==2 ) then ! 2d quad
     mgaus(1)= ngqpt(1)**2
  
  else if( optele==3 ) then ! 3d tb shell
     mgaus(1)= ngqpt(1)**2

  else if( optele==4 ) then ! 3d brick
     mgaus(1)= ngqpt(1)**3

  end if
 
  ! in plane enriched element: ngqpt(2)
  ! -------------------------
  select case(optmhd)
  case(0:1) ! conventional fem or element delection
     mgaus(2)= 0

  case(2) ! discrete phantom nodes crack method
     mgaus(2)= 1

  case(3:4) ! phantom nodes method
     if( optele==1 ) then ! 2d tri
        mgaus(2)= ngqpt(1)

     else if( optele==2 ) then ! 2d quad
        mgaus(2)= ngqpt(1)**2

     else if( optele==3 ) then ! 3d tb shell
        mgaus(2)= ngqpt(1)**2

     else if( optele==4 ) then ! 3d brick
        mgaus(2)= ngqpt(1)**3

     end if

  case(5:6) ! xfem method
     if( optele==1 ) then ! 2d tri
        mgaus(2)= ngqpt(2)

     else if( optele==2 ) then ! 2d quad
        mgaus(2)= ngqpt(2)**2

     else if( optele==3 ) then ! 3d tb shell
        mgaus(2)= ngqpt(2)**2

     else if( optele==4 ) then ! 3d brick
        mgaus(2)= ngqpt(2)**3

     end if

  end select

  ! through thickness: ngqpt(3)
  ! -----------------
  if( optele==0 ) then ! 1d line
     mgaus(3)= 1

  else if( optele==1 ) then ! 2d tri
     mgaus(3)= 1

  else if( optele==2 ) then ! 2d quad
     mgaus(3)= 1

  else if( optele==3 ) then ! 3d tb shell
     mgaus(3)= ngqpt(3)

  else if( optele==4 ) then ! 3d brick
     mgaus(3)= 1

  end if

  ! --------------------------------------------------------------------
  ! set mgqpt(1:2)
  ! --------------

  ! total gauss point: total= inplane
  ! -----------------
  mgqpt(1)= mgaus(1) * mgaus(3) ! regular element
  mgqpt(2)= mgaus(2) * mgaus(3) ! enriched element

  ! --------------------------------------------------------------------



  return
end subroutine getgqsize





subroutine chkeledat(optele,nelem, eletyp)
  !=======================================================================
  !  chkeledat= check element type data from gid flag
  !
  !             note:
  !             -----
  !             - gid                       - code_xed: optele
  !               1 = linear                  0 = linear
  !               2 = triangle                1 = triangle
  !               3 = quadrilateral           2 = 2d quadrilateral
  !               4 = tetrahedra              3 = 3d bt shell
  !               5 = hexahedra               4 = 3d brick
  !               6 = prism                   5 = 3d hybrid: bt shell + brick
  !               7 = only points
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
  integer, intent(in) :: nelem

  integer, dimension(2,*), intent(inout) :: eletyp
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: ielem
  ! ====================================
  
  ! initialize: do not initialize


  select case(optele)
  case(0) ! 1d line
     eletyp(1,1:nelem)= 0

  case(1) ! 2d tri
     eletyp(1,1:nelem)= 1

  case(2) ! 2d quad
     eletyp(1,1:nelem)= 2

  case(3) ! bt shell
     eletyp(1,1:nelem)= 3 ! xed: bt shell
  
  case(4) ! 3d brick
     eletyp(1,1:nelem)= 4 ! xed: 3d brick

  ! --------------------------------------------------------------
  ! hybrid element mesh
  ! -------------------
  case(5) ! 3d hybrid: bt shell + brick
     do ielem=1, nelem
     
       if ( eletyp(1,ielem) == 3 ) then ! gid: quad.
          eletyp(1,ielem)= 3            ! xed: bt shell

       else if  ( eletyp(1,ielem) == 5 ) then ! gid: hexahedra
          eletyp(1,ielem)= 4                  ! xed: 3d brick

       else
          write(*,*) "wrong eletyp from gid: chkeledat"
          write(nout5,*) "wrong eletyp from gid: chkeledat"
       
       end if
  
     end do
  ! --------------------------------------------------------------

  case default
     write(*,*) "wrong optele: chkeledat"
     write(nout5,*) "wrong opteled: chkeledat"

  end select
  
  

  return
end subroutine chkeledat





subroutine getsimdat1(prbsiz,prbopt,prbcon, &
                      npoin,nelem,nmatp,nxele,nxtip,ndime,nnode,ntdof,nrdof,ngqpt,ngqpt4,nibcd,ntbcd,ncrcd, &
                      optmhd,optsim,optele,optcnt,optcoh,opthgc,optdmp,optcri,optprn)
  !=======================================================================
  !  getsimdat1 = get problem information
  !
  !               note:
  !               ----
  !               changed for code_xed3d
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  prbsiz(10) : problem size data
  !
  !  prbopt(10,10) : problem option data
  !
  !  prbcon(3,10) : problem condition data
  !
  !  output:
  !  ------
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, dimension(10), intent(in) :: prbsiz
  integer, dimension(10,10), intent(in) :: prbopt
  integer, dimension(3,10), intent(in) :: prbcon

  integer, intent(out) :: npoin, nelem, nmatp, nxele, nxtip
  integer, intent(out) :: ndime, nnode, ntdof, nrdof
  integer, dimension(3), intent(out) :: ngqpt
  integer, intent(out) :: ngqpt4
  integer, dimension(10), intent(out) :: nibcd, ntbcd, ncrcd
  integer, intent(out) :: optmhd, optsim, optele
  integer, dimension(10), intent(out) :: optcnt
  integer, intent(out) :: optcoh, opthgc, optdmp
  integer, dimension(10), intent(out) :: optcri, optprn
  ! ====================================

  ! initialize
  ngqpt(:)= 0
  nibcd(:)= 0
  ntbcd(:)= 0
  ncrcd(:)= 0
  optcnt(:)= 0
  optcri(:)= 0
  optprn(:)= 0


  ! ---------------------------------------------------------------
  ! parsing read size data
  ! ------------------------
  ! problem size
  npoin= prbsiz(1) ! node
  nelem= prbsiz(2) ! element
  nmatp= prbsiz(3) ! assigned material properties set
  nxele= prbsiz(4) ! max extra elements 
  nxtip= prbsiz(5) ! max extra tips

  ! ---------------------------------------------------------------
  ! parsing read option data
  ! ------------------------
  ! 1. general
  optmhd= prbopt(1,1) ! numerical method
  ndime= prbopt(1,2) ! dimension
  optsim= prbopt(1,3) ! simulation type

  ! 2. element
  optele= prbopt(2,1) ! element
  nnode= prbopt(2,2) ! max. node per element
  ntdof= prbopt(2,3) ! max. trans. dof per node
  nrdof= prbopt(2,4) ! max. rot. dof per node

  ngqpt(1)= prbopt(2,5) ! gq rule for regular domain
  ngqpt(2)= prbopt(2,6) ! gq rule for sub domain
  ngqpt(3)= prbopt(2,7) ! gq rule for through thickness
  ngqpt4= prbopt(2,8) ! gq rule for bc or cohesive force integeation


  ! 3. contact
  optcnt(1)= prbopt(3,1) ! contact search algorithm
  optcnt(2)= prbopt(3,2) ! contact body sets
  optcnt(3)= prbopt(3,3) ! x direction
  optcnt(4)= prbopt(3,4) ! y direction
  optcnt(5)= prbopt(3,5) ! z direction
 
  ! 4. cohesive model
  optcoh= prbopt(4,1) ! cohesive model

  ! 5. hg control
  opthgc= prbopt(5,1) ! hourglass control

  ! 6. crack criterion
  optcri(1)= prbopt(6,1) ! criterion type
  optcri(2)= prbopt(6,2) ! criterion average

  ! 7. time function
  
  ! 8. output 
  optprn(1)= prbopt(8,1) ! result print interval
  optprn(2)= prbopt(8,2) ! crack path file
  optprn(3)= prbopt(8,3) ! log echo message file
  optprn(4)= prbopt(8,4) ! echo monitor message
  optprn(5)= prbopt(8,5) ! nodal displacement
  optprn(6)= prbopt(8,6) ! nodal velocity
  optprn(7)= prbopt(8,7) ! stress
  optprn(8)= prbopt(8,8) ! load-deflection curve

  ! 9. damping
  optdmp= prbopt(9,1) ! damping type

  ! ---------------------------------------------------------------
  ! parsing read condition data
  ! ---------------------------
  ! condition
  nibcd(1)= prbcon(1,1) ! ic vel
  nibcd(2)= prbcon(1,2) ! ic frc
  nibcd(3)= prbcon(1,3) !
  nibcd(4)= prbcon(1,4) !
  nibcd(5)= prbcon(1,5) !
  nibcd(6)= prbcon(1,6) !
  ! ---------------
  ntbcd(1)= prbcon(2,1) ! bc vel
  ntbcd(2)= prbcon(2,2) ! bc frc
  ntbcd(3)= prbcon(2,3) ! bc prs
  ntbcd(4)= prbcon(2,4) !
  ntbcd(5)= prbcon(2,5) !
  ntbcd(6)= prbcon(2,6) !
  ! --------------

  ncrcd(1)= prbcon(3,1) ! total number of cracking control
  ncrcd(2)= prbcon(3,2) ! total number of load deflection curve point
  ! ---------------------------------------------------------------



  return
end subroutine getsimdat1





subroutine getsimprm1(prbprm, tmax,tmfct,prmcnt,prmcoh,prmhgc,prmcri,prmtft,prmdmp)
  !=======================================================================
  !  getprobprm1 = get general parameter
  ! 
  !                note:
  !                ----
  !                changed for code_xed3d
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  prbprm(10,10) : problem paramters
  !
  !  output:
  !  ------
  !  tmax : maximum simulation time
  !
  !  tmfct : time step control factor
  !
  !  prmcnt(10) : contact parameters
  !
  !  prmcoh(10) : cohesive crack model parameters
  !
  !  prmhgc(10) : hourglass control parameters
  !
  !  prmcri(10) : crack criterion parameters
  !
  !  prmtft(10) : loading time function parameters
  !
  !  prmdmp(10) : damping parameter
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(10,10), intent(in) :: prbprm

  real(8), intent(out) :: tmax,tmfct
  real(8), dimension(10), intent(out) :: prmcnt
  real(8), dimension(10), intent(out) :: prmcoh
  real(8), dimension(10), intent(out) :: prmhgc
  real(8), dimension(10), intent(out) :: prmcri
  real(8), dimension(10), intent(out) :: prmtft
  real(8), dimension(10), intent(out) :: prmdmp
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: indx
  ! ====================================

  ! ---------------------------------------------------------------
  ! parsing read parameter data
  ! ---------------------------
  ! 1. general
  tmax= prbprm(1,1)  ! simulation time / iteration
  tmfct= prbprm(1,2) ! courant number
  ! ------------------
  toler(1)= prbprm(1,3) ! tolerance : global variables
  toler(2)= prbprm(1,4)
  toler(3)= prbprm(1,5) 

  do indx=1, 10
     ! 3. contact
     prmcnt(indx)= prbprm(3,indx)

     ! 4. cohesive model
     prmcoh(indx)= prbprm(4,indx)

     ! 5. hg control
     prmhgc(indx)= prbprm(5,indx)

     ! 6. crack criterion
     prmcri(indx)= prbprm(6,indx)

     ! 7. time function
     prmtft(indx)= prbprm(7,indx)

     ! 9. damping
     prmdmp(indx)= prbprm(9,indx)
  end do

  ! note: no parameters for
  ! -----
  ! 2. element
  ! 8. output 



  return
end subroutine getsimprm1





subroutine getgqsize1(optmhd,nelem,eletyp,ngqpt, mgaus,mgqpt)
  !=======================================================================
  !  getgqsize1 = get general parameter
  !
  !               note:
  !               ----
  !               this subroutine gives elementwise gq data
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optmhd : numerical method option
  !
  !  nelem : total number of elements
  !
  !  eletyp(2,nelem) : element type flag
  !
  !  ngqpt(3) : gq rule
  !             ngqpt(1) : gq rule for regular element 
  !             ngqpt(2) : gq rule for enriched element 
  !             ngqpt(3) : gq rule through thickness
  !
  !  output:
  !  ------
  !  mgaus(3,nelem) : the number of quadrature point
  !                   mgaus(1)= regular element
  !                   mgaus(2)= enriched element
  !                   mgaus(3)= through thickness
  !
  !  mgqpt(2,nelem) : the total number of quadrature point
  !                   mgqpt(1)= mgaus(1) * mgaus(3) : regular element
  !                   mgqpt(2)= mgaus(2) * mgaus(3) : enriched element
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optmhd, nelem
  integer, dimension(2,*), intent(in) :: eletyp
  integer, dimension(3), intent(in) :: ngqpt
  
  integer, dimension(3,nelem), intent(out) :: mgaus
  integer, dimension(2,nelem), intent(out) :: mgqpt
  ! ====================================
  ! local variable
  ! ==============
  integer :: optele

  ! loop index
  integer :: ielem
  ! ====================================

  ! initialize
  mgaus(:,:)= 0
  mgqpt(:,:)= 0


  do ielem=1, nelem

     ! get element type
     ! ----------------
     optele= eletyp(1,ielem)

     ! --------------------------------------------------------------------
     ! set mgaus(1:3)
     ! --------------
     ! in plane regular element: ngqpt(1)
     ! ------------------------
     if( optele==0 ) then ! 1d line
        mgaus(1,ielem)= ngqpt(1)

     else if( optele==1 ) then ! 2d tri
        mgaus(1,ielem)= ngqpt(1)

     else if( optele==2 ) then ! 2d quad
        mgaus(1,ielem)= ngqpt(1)**2
     
     else if( optele==3 ) then ! 3d tb shell
        mgaus(1,ielem)= ngqpt(1)**2

     else if( optele==4 ) then ! 3d brick
        mgaus(1,ielem)= ngqpt(1)**3

     end if
    
     ! in plane enriched element: ngqpt(2)
     ! -------------------------
     select case(optmhd)
     case(0:1) ! conventional fem or element delection
        mgaus(2,ielem)= 0

     case(2) ! discrete phantom nodes crack method
        mgaus(2,ielem)= 1

     case(3:4) ! phantom nodes method
        if( optele==1 ) then ! 2d tri
           mgaus(2,ielem)= ngqpt(1)

        else if( optele==2 ) then ! 2d quad
           mgaus(2,ielem)= ngqpt(1)**2

        else if( optele==3 ) then ! 3d tb shell
           mgaus(2,ielem)= ngqpt(1)**2

        else if( optele==4 ) then ! 3d brick
           mgaus(2,ielem)= ngqpt(1)**3

        end if

     case(5:6) ! xfem method
        if( optele==1 ) then ! 2d tri
           mgaus(2,ielem)= ngqpt(2)

        else if( optele==2 ) then ! 2d quad
           mgaus(2,ielem)= ngqpt(2)**2

        else if( optele==3 ) then ! 3d tb shell
           mgaus(2,ielem)= ngqpt(2)**2

        else if( optele==4 ) then ! 3d brick
           mgaus(2,ielem)= ngqpt(2)**3

        end if

     end select

     ! through thickness: ngqpt(3)
     ! -----------------
     if( optele==0 ) then ! 1d line
        mgaus(3,ielem)= 1

     else if( optele==1 ) then ! 2d tri
        mgaus(3,ielem)= 1

     else if( optele==2 ) then ! 2d quad
        mgaus(3,ielem)= 1

     else if( optele==3 ) then ! 3d tb shell
        mgaus(3,ielem)= ngqpt(3)

     else if( optele==4 ) then ! 3d brick
        mgaus(3,ielem)= 1

     end if

     ! --------------------------------------------------------------------
     ! set mgqpt(1:2)
     ! --------------

     ! total gauss point: total= inplane
     ! -----------------
     mgqpt(1,ielem)= mgaus(1,ielem) * mgaus(3,ielem) ! regular element
     mgqpt(2,ielem)= mgaus(2,ielem) * mgaus(3,ielem) ! enriched element

     ! --------------------------------------------------------------------

  end do ! do ielem=1, nelem



  return
end subroutine getgqsize1





subroutine parsctvinp(line, optctv) 
  !=======================================================================
  !  parsctvinp = parse constitutive model name input
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  line(100) : character label of constitutive model name
  !
  !  output:
  !  ------
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  character(len=100), intent(in) :: line

  integer, intent(out) :: optctv
  ! ====================================
  character(len=100) :: lab

  ! loop index
  integer :: indx
  ! ====================================

  ! initialize
  optctv= 0


  ! constitutive model option 
  ! -------------------------
  !   label         optctv        description
  !   -----         ------        ------------
  !   hypoelas         1          hypo-elasticity
  !   tevp1            2          thermo-elasto-visco-plasticity
  !   evp1             3          elasto=visco-plasticity
  !   gtn              4          gurson-tvergaard-needleman
  !   j2exp            5          j2-plasticity with explicit integration
  !   dmg1             6          hypoelstic_w/_damage1 (lematire)
  !   dmg2             7          hypoelstic_w/_damage2 (linear)
  !   ve1              8          viscous hyperelastic1 (ogden)
  !   tjc1             9          thermo-johnson-cook

  ! loop over read line
  do indx=1, 100

     ! read until "/" appear
     if( line(indx:indx)=='/' ) then

        lab= line(1:indx-1)
        exit

     end if
  
  end do

  ! set constitutive model option number
  if( lab == 'hypoelas' ) then
     optctv= 1

  else if( lab == 'tevp1' ) then 
     optctv= 2

  else if( lab == 'evp1' ) then 
     optctv= 3

  else if( lab == 'gtn' ) then 
     optctv= 4

  else if( lab == 'j2exp' ) then 
     optctv= 5

  else if( lab == 'dmg1' ) then 
     optctv= 6

  else if( lab == 'dmg2' ) then 
     optctv= 7

  else if( lab == 've1' ) then 
     optctv= 8

  else if( lab == 'tjc1' ) then 
     optctv= 9

  else
     write(*,*) "wrong constitutive label: parsctvinp"
     write(nout5,*) "wrong constitutive label: parsctvinp"
  
  end if
  

  return
end subroutine parsctvinp