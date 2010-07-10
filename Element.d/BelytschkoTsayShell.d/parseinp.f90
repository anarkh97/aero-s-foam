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

  use preset
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

  use preset
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

  use preset
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

! ### working
if ( jcond == -1 .or. jcond == 1 ) then
   status(jelem)= jcond

else if ( jcond == 2 .or. jcond == 3) then
   status(jelem)= -jcond
   ! note: (-) is used to represent initial crack: see, coheisve force

end if



!        select case(optmhd)
!        case(0) ! conventional fem
!           if ( jcond == -1 ) status(jelem)= jcond
     
!        case(1) ! element deletion
!           if ( jcond == -1 .or. jcond == 1 ) status(jelem)= jcond

!        case(2:7) ! phantom nodes, xfem and dimensional coupling
!           if ( jcond == -1 .or. jcond == 1 ) status(jelem)= jcond
!           if ( jcond == 2 .or. jcond == 3) status(jelem)= -jcond
           ! note: (-) is used to represent initial crack: see, coheisve force

!        end select

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

  use preset
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