! ================================
! j2 plasticity model model
! ================================
!      type                  name                              arguement
!      ----                  ----                              ---------
!  1. subroutine        getpropj21            (ematpro, e,nu,rho,sigy,h)
!  2. subroutine        getj2dlambda1         (h,j2dirvoit,elsvoid2d,ipstrndel, dlambda)
!  3. subroutine        updj2pexp             (ematpro,ipstrndel, effpstrn,hardvar,ipstrn,ipstrs, effstrs) 
!
! =========================================================================================================



subroutine getpropj21(ematpro, e,nu,rho,sigy,h)
  !=======================================================================
  !  getpropj21 = get j21 model material properties
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

  use preset
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
end subroutine getpropj21





subroutine getj2dlambda1(h,j2dirvoit,elsvoid2d,ipstrndel, dlambda)
  !=======================================================================
  !  getj2dlambda1 = explicit calculation of dlambda(plastic flow)
  !                 in j2 plasticity
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : material properties
  !
  !  ipstrndel(3,1) : increment of in plane strain tensor
  !
  !  inoutput:
  !  --------
  !  effpstrn : effective plastic strain
  !
  !  hardvar : hardening variable
  !
  !  ipstrn(3,1) : voight form in plane strain
  !
  !  ipstrs(3,1) : voight form in plane stress
  !
  !  output:
  !  ------
  !  effstrs : effective stress
  !
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: h
  real(8), dimension(3,1), intent(in) :: j2dirvoit
  real(8), dimension(3,3), intent(in) :: elsvoid2d
  real(8), dimension(3,1), intent(in) :: ipstrndel

  real(8), intent(out) :: dlambda
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,1) :: temp1
  real(8) :: temp2
  real(8), dimension(3,1) :: temp3
  real(8) :: temp4

  ! loop index
  ! ====================================

  ! initialize
  dlambda= 0.0d0

  ! elasvoit2d[3x3]*ipstrndel[3,1]
  call matprd(3,3,0, 3,1,0, 3,1, elsvoid2d,ipstrndel, temp1)
     ! input : 3,3,0, 3,1,0, 3,1, elsvoid2d,ipstrndel
     ! output : temp1

  ! n.elasvoit2d[3x3].ipstrndel[3,1]
  temp2= j2dirvoit(1,1)*temp1(1,1) + j2dirvoit(2,1)*temp1(2,1) + j2dirvoit(3,1)*temp1(3,1)

  ! elasvoit2d[3x3]*j2dirvoit[3,1]
  call matprd(3,3,0, 3,1,0, 3,1, elsvoid2d,j2dirvoit, temp3)
     ! input : 3,3,0, 3,1,0, 3,1, elsvoid2d,j2dirvoit
     ! output : temp3

  ! n.elasvoit2d[3x3].n
  temp4= j2dirvoit(1,1)*temp3(1,1) + j2dirvoit(2,1)*temp3(2,1) + j2dirvoit(3,1)*temp3(3,1)

  ! compute dlambda
  dlambda= dabs(temp2 / ( temp4 + h) )



  return
end subroutine getj2dlambda1





subroutine updj2pexp(ematpro,ipstrndel, effpstrn,hardvar,ipstrn,ipstrs, effstrs) 
  !=======================================================================
  !  updj2pexp = explicitly update j2 plasticity
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : material properties
  !
  !  ipstrndel(3,1) : increment of in plane strain tensor
  !
  !  inoutput:
  !  --------
  !  effpstrn : effective plastic strain
  !
  !  hardvar : hardening variable
  !
  !  ipstrn(3,1) : voight form in-plane strain
  !
  !  ipstrs(3,1) : voight form in-plane stress
  !
  !  output:
  !  ------
  !  effstrs : effective stress
  !
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(3,1), intent(in) :: ipstrndel

  real(8), intent(inout) :: effpstrn, hardvar
  real(8), dimension(3,1), intent(inout) :: ipstrn, ipstrs

  real(8), intent(out) :: effstrs
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: e, nu, rho, sigy, h

  real(8), dimension(3,3) :: elsvoid2d
  real(8), dimension(2,2) :: ipstrstens
  real(8), dimension(2,2) :: devstrstens
  real(8) :: yield

  real(8), dimension(2,2) :: j2dirtens
  real(8), dimension(3,1) :: j2dirvoit

  real(8) :: dlambda

  real(8), dimension(3,1) :: pstrndel, estrndel, estrsdel

  ! loop index
  ! ====================================

  ! initialize
  effstrs= 0.0d0


  ! get material properties
  call getpropj21(ematpro, e,nu,rho,sigy,h)
     ! input : ematpro
     ! output : e,nu,rho,sigy,h


  ! get virgin elastic moduli
  call getelsvoit2d(1,e,nu, elsvoid2d)
     ! input : 1(optpty:p strs),e,nu
     ! output : elsvoit2d


  ! change voight form to tensor
  ipstrstens(1,1)= ipstrs(1,1) 
  ipstrstens(2,2)= ipstrs(2,1) 
  ipstrstens(1,2)= ipstrs(3,1) 
  ipstrstens(2,1)= ipstrs(3,1) 


  ! compute deviatoric stress
  call getdevtens(2,ipstrstens, devstrstens)
     ! input : 2(nndex),ipstrstens
     ! output : devstrstens

  ! compute effective stress
  call geteffstrs(2,devstrstens, effstrs)
     ! input : 2(nndex),devstrstens
     ! output : effstrn

  ! check yield condition
  yield= effstrs - hardvar - sigy


  ! if yield criterion is satisfied
  if ( yield > 0.0d0 ) then

     ! calculate plastic flow direction
     call getj2dirtens(2,devstrstens,effstrs, j2dirtens)
        ! input : 2(nndex),devstrstens,effstrs
        ! output : j2dirtens

     ! voight form
     j2dirvoit(1,1)= j2dirtens(1,1)
     j2dirvoit(2,1)= j2dirtens(2,2)
     j2dirvoit(3,1)= j2dirtens(1,2)

     ! compute dlambda: plastic multiplier
     call getj2dlambda1(h,j2dirvoit,elsvoid2d,ipstrndel, dlambda)
        ! input : h,j2dirvoit,elsvoid2d,ipstrndel
        ! output : dlambda

  else
     dlambda= 0.0d0

  end if

  ! plastic strain increment
  pstrndel(1:3,1)= dlambda * j2dirvoit(1:3,1)

  ! elastic strain increment
  estrndel(1:3,1)= ipstrndel(1:3,1) - pstrndel(1:3,1)

  ! elastic stress increment
  call matprd(3,3,0, 3,1,0, 3,1, elsvoid2d,estrndel, estrsdel)
     ! input : 3,3,0, 3,1,0, 3,1, elsvoid2d,estrndel
     ! output : estrsdel

  ! explicit update
  effpstrn= effpstrn + dlambda ! effective plastic strain

  hardvar= hardvar + h*dlambda ! hardening variable

  ipstrn(1:3,1)= ipstrn(1:3,1) + ipstrndel(1:3,1) ! estrndel(1:3,1) ! update strain

  ipstrs(1:3,1)= ipstrs(1:3,1) + estrsdel(1:3,1) ! update stress



  return
end subroutine updj2pexp