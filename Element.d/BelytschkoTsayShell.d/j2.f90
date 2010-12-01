! ================================
! j2 plasticity model model: explicit integration
! ================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine        getpropj2            (ematpro, e,nu,rho,sigy,h)
! 2.  subroutine        getj2dlambdaexp      (nvoit,h,j2dirvoit,elsvoid,dstrnvoit, dlambda)
! 3.  subroutine        updj2pexp            (optpty,ndime,ematpro,dstrnvoit, effpstrn,hardvar,strnvoit,strsvoit, effstrs) 
!
! =========================================================================================================



subroutine getpropj2(ematpro, e,nu,rho,sigy,h)
  !=======================================================================
  !  getpropj2 = get j2 model material properties
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : evp material properties
  !
  !  output:
  !  ------
  !  e,nu,rho,eps0dot,m,sig0,eps0,n : constitutive model properties
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro

  real(8), intent(out) :: e,nu,rho,sigy,h
  ! ====================================
 
  ! get material properties
  ! -----------------------
  e= ematpro(1)      ! young's modulus
  nu= ematpro(2)     ! poisson's ratio
  rho= ematpro(3)    ! mass density

  sigy= ematpro(4)   ! yield stress
  h= ematpro(5)      ! hardening modulus



  return
end subroutine getpropj2





subroutine getj2dlambdaexp(nvoit,h,j2dirvoit,elsvoit,dstrnvoit, dlambda)
  !=======================================================================
  !  getj2dlambdaexp = explicit calculation of dlambda (plastic flow)
  !                    in j2 plasticity
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
  integer, intent(in) :: nvoit
  real(8), intent(in) :: h
  real(8), dimension(nvoit,1), intent(in) :: j2dirvoit
  real(8), dimension(nvoit,nvoit), intent(in) :: elsvoit
  real(8), dimension(nvoit,1), intent(in) :: dstrnvoit

  real(8), intent(out) :: dlambda
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nvoit,1) :: temp1
  real(8) :: temp2
  real(8), dimension(nvoit,1) :: temp3
  real(8) :: temp4

  ! loop index
  integer :: ivoit
  ! ====================================

  ! initialize
  dlambda= 0.0d0


  ! elasvoit2d[nvoit x nvoit]*delstrnvoit[nvoit,1]
  call matprd(nvoit,nvoit,0, nvoit,1,0, nvoit,1, elsvoit,dstrnvoit, temp1)
     ! input : nvoit,nvoit,0, nvoit,1,0, nvoit,1, elsvoit,dstrnvoit
     ! output : temp1

  ! j2dirvoit[nvoit x nvoit].temp1[nvoit x 1]
  temp2= 0.0d0 ! initialize
  do ivoit=1, nvoit
     temp2= temp2 + j2dirvoit(ivoit,1)*temp1(ivoit,1)
  end do

  ! elasvoit[nvoit x nvoit].j2dirvoit[nvoit,1]
  call matprd(nvoit,nvoit,0, nvoit,1,0, nvoit,1, elsvoit,j2dirvoit, temp3)
     ! input : nvoit,nvoit,0, nvoit,1,0, nvoit,1, elsvoit,j2dirvoit
     ! output : temp3

  ! j2dirvoit[nvoit x nvoit].temp3[nvoit x 1]
  temp4= 0.0d0 ! initialize
  do ivoit=1, nvoit
     temp4= temp4 + j2dirvoit(ivoit,1)*temp3(ivoit,1)
  end do


  ! compute dlambda
  dlambda= dabs( temp2 / ( temp4 + h ) )



  return
end subroutine getj2dlambdaexp





subroutine updj2pexp(optpty,ndime,ematpro,delt,ipstrndot, effpstrn,hardvar,strsvoit, effstrs) 
  !=======================================================================
  !  updj2pexp = explicitly update j2 plasticity
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optpty : 2d problem type
  !
  !  ndime : dimension
  !
  !  ematpro(*) : material properties
  !
  !  delt : time step
  !
  !  ipstrndot(nvoit,1) : strain rate
  !
  !  inoutput:
  !  --------
  !  effpstrn : effective plastic strain
  !
  !  hardvar : hardening variable
  !
  !  strsvoit(nvoit,1) : voight form stress
  !
  !  output:
  !  ------
  !  effstrs : effective stress
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optpty
  integer, intent(in) :: ndime
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: delt
  real(8), dimension(ndime*(ndime+1)/2,1), intent(in) :: ipstrndot
  ! -------------------------------------

  real(8), intent(inout) :: effpstrn
  real(8), intent(inout) :: hardvar
  real(8), dimension(ndime*(ndime+1)/2,1), intent(inout) :: strsvoit

  ! -------------------------------------
  real(8), intent(out) :: effstrs
  ! ====================================
  ! local variable
  ! ==============
  integer :: nvoit

  real(8) :: e, nu, rho, sigy, h

  ! total strain increment in voigt form
  real(8), dimension(ndime*(ndime+1)/2,1) :: dstrnvoit

  ! voight form
  real(8), dimension(ndime*(ndime+1)/2,ndime*(ndime+1)/2) :: elsvoit

  ! tensor
  real(8), dimension(ndime,ndime) :: strstens
  real(8), dimension(ndime,ndime) :: devstrstens
  real(8), dimension(ndime,ndime) :: j2dirtens

  ! voight form
  real(8), dimension(ndime*(ndime+1)/2,1) :: j2dirvoit

  real(8) :: yield
  real(8) :: dlambda

  ! voight form
  real(8), dimension(ndime*(ndime+1)/2,1) :: dpstrnvoit
  real(8), dimension(ndime*(ndime+1)/2,1) :: destrnvoit
  real(8), dimension(ndime*(ndime+1)/2,1) :: destrsvoit
  ! ====================================

  ! compute number of voight form components
  nvoit= ndime * ( ndime + 1 ) / 2

  ! strain increment
  ! ----------------
  dstrnvoit(1:3,1)= delt * ipstrndot(1:3,1)

  ! --------------------------------------------------------------
  ! material properties
  ! -------------------
  ! get material properties
  call getpropj2(ematpro, e,nu,rho,sigy,h)
     ! input : ematpro
     ! output : e,nu,rho,sigy,h

  ! get virgin elastic moduli
  if ( ndime==2 ) then
     call getelsvoit2d(optpty,e,nu, elsvoit)
        ! input : optpty,e,nu
        ! output : elsvoit

  else if (ndime == 3 ) then
     call getelsvoit3d(e,nu, elsvoit)
        ! input : e,nu
        ! output : elsvoit
  end if


  ! TODO better to compute the effective stress after updating
  ! so we can use it for post processing. Then the effective
  ! stress used in the yield condition check below should be
  ! taken from evar1
  ! --------------------------------------------------------------
  ! compute effective stress
  ! ------------------------
  ! convert strsvoit to strstens
  call voit2ind(0,ndime,strsvoit, strstens)
     ! input : 0(ntype:kinetic),ndime,strsvoit
     ! output : strstens

  ! compute deviatoric stress
  call getdevtens(ndime,strstens, devstrstens)
     ! input : ndime,strstens
     ! output : devstrstens

  ! compute effective stress
  call geteffstrs(ndime,devstrstens, effstrs)
     ! input : ndime,devstrstens
     ! output : effstrn


  ! --------------------------------------------------------------
  ! check yield condition
  ! ---------------------
  yield= effstrs - hardvar - sigy

  if ( yield > 0.0d0 ) then ! yield criterion is satisfied

     ! calculate plastic flow direction
     call getj2dirtens(ndime,devstrstens,effstrs, j2dirtens)
        ! input : ndime,devstrstens,effstrs
        ! output : j2dirtens

     ! get voight form j2 plastic flow direction
     call ind2voit(0,ndime,j2dirtens, j2dirvoit)
        ! input : 0(ntype:kinetic),ndime,j2dirtens
        ! output : j2dirvoit  
  
     ! compute dlambda: plastic multiplier
     call getj2dlambdaexp(nvoit,h,j2dirvoit,elsvoit,dstrnvoit, dlambda)
        ! input : nvoit,h,j2dirvoit,elsvoit,dstrnvoit
        ! output : dlambda

  else
     dlambda= 0.0d0
     j2dirvoit(:,:)= 0.0d0
  end if


  ! --------------------------------------------------------------
  ! compute elastic stress increment
  ! --------------------------------
  ! plastic strain increment
  dpstrnvoit(1:nvoit,1)= dlambda * j2dirvoit(1:nvoit,1)

  ! elastic strain increment
  destrnvoit(1:nvoit,1)= dstrnvoit(1:nvoit,1) - dpstrnvoit(1:nvoit,1)

  ! elastic stress increment
  call dgemv('n',nvoit,nvoit,1.0d0,elsvoit,nvoit,destrnvoit,1,0.0d0,destrsvoit,1)
     ! input : nvoit,nvoit,0, nvoit,1,0, nvoit,1, elsvoit,destrnvoit
     ! output : destrsvoit

  ! --------------------------------------------------------------
  ! explicit update
  ! ---------------
  ! effective plastic strain
  effpstrn= effpstrn + dlambda

  ! hardening variable
  hardvar= hardvar + h*dlambda

  ! update stress
  strsvoit(1:nvoit,1)= strsvoit(1:nvoit,1) + destrsvoit(1:nvoit,1)
  ! --------------------------------------------------------------

  return
end subroutine updj2pexp
