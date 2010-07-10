! ================================
! elasto-viscoplastic model
! ================================
!      type                  name                              arguement
!      ----                  ----                              ---------
!  1. subroutine        getpropevp1           (ematpro, e,nu,rho,eps0dot,m,sig0,eps0,n)
!  2. subroutine        getetaevp1            (ematpro,effstrs,effpstrn, eta)
!  3. subroutine        gethnpsievp1          (opt,ematpro,delt,prveffstrs,mideffstrs,effpstrn, h,psi)
!  4. subroutine        gettantensevp1        (ematpro,nndex,lptens,h,psi,effstrs, tantens)
!  5. subroutine        getkirjevp1           (ematpro,nndex,deftens,pdirtens,delt,effstrs,effpstrn, tantens,kirjtens)
!  6. subroutine        updkirtevp1           (ematpro,nndex,middeftens,delt,effpstrn, kirtens, ctan,midkirtens,effstrs) 
!  7. subroutine        updpstrnevp1          (ematpro,nndex,middeftens,prvkirtens,midkirtens,delt, effpstrn) 
!  8. subroutine        updevp1               (ematpro,nndex,midltens,delt, effpstrn,kirtens, ctan,effstrs)
!
! =========================================================================================================



subroutine getpropevp1(ematpro, e,nu,rho,eps0dot,m,sig0,eps0,n)
  !=======================================================================
  !  getpropevp1 = get evp1 model material properties
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

  real(8), intent(out) :: e,nu,rho,eps0dot,m,sig0,eps0,n
  ! ====================================
 
  ! get material properties
  ! -----------------------
  e= ematpro(1)      ! young's modulus
  nu= ematpro(2)     ! poisson's ratio
  rho= ematpro(3)    ! mass density

  eps0dot= ematpro(4) ! reference strain rate
  m= ematpro(5)       ! rate senstivity parameter
  sig0= ematpro(6)    ! yield stress
  eps0= ematpro(7)    ! yield strain
  n= ematpro(8)       ! strain hardening exponent



  return
end subroutine getpropevp1



subroutine getetaevp1(ematpro,effstrs,effpstrn, eta)
  !=======================================================================
  !  getetaevp1 = compute effective plastic strain flow rate: eta
  !               for elasto-viscoplastic model
  !
  !               note:
  !               ----
  !               computer & structure, vol 18, 1984, pp. 875-887 (Eq. 4.1)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : tevp material properties
  !
  !  effstrs : effective stress
  !
  !  effpstrn : effective plastic strain
  !
  !  output:
  !  ------
  !  eta : effective plastic strain flow rate
  !
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: effstrs,effpstrn

  real(8), intent(out) :: eta
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: e,nu,rho,eps0dot,m,sig0,eps0,n

  real(8) :: gft
  ! ====================================

  ! initialize
  eta= 0.0d0


  ! get material properties
  ! -----------------------
  call getpropevp1(ematpro, e,nu,rho,eps0dot,m,sig0,eps0,n)
     ! input : ematpro
     ! output : e,nu,rho,eps0dot,m,sig0,eps0,n


  ! compute g function value: g(effpstrn,tmp)   (Eq. 4.1)
  ! ------------------------
  gft= sig0 * ( (1.0d0 + effpstrn/eps0)**n )


  ! compute etabardot: eta(effstrs,g(effpstrn,tmp))   (Eq. 3.10) 
  ! -----------------
  eta= min(1.0e10, eps0dot * ( effstrs / gft )**m ) ! effective plastic strain flow rate



  return
end subroutine getetaevp1



subroutine gethnpsievp1(opt,ematpro,delt,prveffstrs,mideffstrs,effpstrn, h,psi)
  !=======================================================================
  !  gethnpsievp1 = compute coefficient h and psi
  !
  !             note:
  !             ----
  !             opt=0 : inconsistent approximation: effpstrs at previous step [n] is used
  !             opt=1 : consistent approximation  : effpstrs at mid step [n+1/2] is used
  !
  !             computer & structure, vol 18, 1984, pp. 875-887 (Eq. 3.11)
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt : option handler
  !        opt=0 : inconsistent approximation: effpstrs at previous step [n] is used
  !        opt=1 : consistent approximation  : effpstrs at mid step [n+1/2] is used
  !
  !  ematpro(*) : tevp material properties
  !
  !  delt : time increment
  !
  !  prveffstrs : effective stress: previous step [n]
  !
  !  mideffstrs : effective stress: mid step [n+1/2]
  !
  !  effpstrn : effective plastic strain: previous step [n]
  !
  !  output:
  !  ------
  !  h, psi : compute coefficient h and psi (Eq. 3.11)
  !
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: delt
  real(8), intent(in) :: prveffstrs,mideffstrs, effpstrn

  real(8), intent(out) :: h, psi
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: e,nu,rho,eps0dot,m,sig0,eps0,n

  real(8) :: mu,lamd

  real(8) :: eta ! eta function value
  real(8) :: deffstrs ! d eta / d effstrs
  real(8) :: deffpstrn ! d eta / d effpstrn
  real(8) :: effstrs ! effective stress
  ! ====================================

  ! initialize
  h= 0.0d0
  psi= 0.0d0


  ! get material properties
  ! -----------------------
  call getpropevp1(ematpro, e,nu,rho,eps0dot,m,sig0,eps0,n)
     ! input : ematpro
     ! output : e,nu,rho,eps0dot,m,sig0,eps0,n
 
  ! get lame constant
  call getlameconst(e,nu, lamd,mu)
    ! input : e,nu
    ! output : lamd,mu

  ! -------------------------------------------------------------------
  ! compute eta function value for convenience: previous step [n]
  call getetaevp1(ematpro,effstrs,effpstrn, eta)
     ! input : ematpro,prveffstrs,effpstrn
     ! output : eta


  ! compute derivative of eta
  if ( prveffstrs == 0.0d0 ) then
     deffstrs= 0.0d0
     deffpstrn= 0.0d0

  else
     ! compute d eta/ d effstrs: previous step [n]
     deffstrs= ( eta * m ) / prveffstrs

     ! compute d eta/ d effpstrn: previous step [n]
     deffpstrn= ( -1.0d0 * eta * m * n ) / ( eps0 * (1.0d0 + effpstrn/eps0) )

  end if
  ! -------------------------------------------------------------------


  ! -----------------
  ! compute h and psi: (Eq. 3.7)
  ! -----------------
  select case(opt)
  case(0) ! inconsistent approximation
     effstrs= prveffstrs
  case(1) ! consistent approximation
     effstrs= mideffstrs
  end select

  if ( deffstrs == 0.0d0 ) then
     h= 0.0d0
     psi= 0.0d0
  else
     h= ( 3.0d0 * mu ) - (deffpstrn/deffstrs)
     psi= 0.50d0 * delt * deffstrs * h
  end if


  return
end subroutine gethnpsievp1



subroutine gettantensevp1(ematpro,nndex,lptens,h,psi,effstrs, tantens)
  !=======================================================================
  !  gettantensevp1 = compute adiabatic tangent stiffness c^tan tensor
  !
  !              note:
  !              ----
  !              computer & structure, vol 18, 1984, pp. 875-887 (Eq. 3.9)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : constitutive model parameter
  !
  !  nndex : dimension of tensor
  !
  !  lptens(nndex,nndex) : P tensor (Eq. 51)
  !
  !  h, psi : coefficient (Eq. 50, 52)
  !
  !  effstrs : effective stress
  !
  !  output:
  !  ------
  !  tantens(nndex,nndex,nndex,nndex) : tangent stiffness tensor
  !
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: lptens
  real(8), intent(in) :: h,psi
  real(8), intent(in) :: effstrs

  real(8), dimension(nndex,nndex,nndex,nndex), intent(out) :: tantens
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: e,nu,rho,eps0dot,m,sig0,eps0,n
  real(8) ::  mu, lamd

  real(8), dimension(nndex,nndex,nndex,nndex) :: elstens

  real(8) :: const1

  real(8) :: dij, krodelt

  ! loop index
  integer :: indx, jndx, kndx, lndx
  ! ====================================

  ! initialize
  tantens(:,:,:,:)= 0.0d0

  ! -------------------------------------------------------------------
  ! get material properties
  ! -----------------------
  call getpropevp1(ematpro, e,nu,rho,eps0dot,m,sig0,eps0,n)
     ! input : ematpro
     ! output : e,nu,rho,eps0dot,m,sig0,eps0,n

  ! get lame constant
  call getlameconst(e,nu, lamd,mu)
     ! input : e,nu
     ! output : lamd,mu

  ! compute elastic moduli tensor
  call getelstens(nndex,lamd,mu, elstens)
     ! input : nndex,lamd,mu
     ! output : elstens
  ! -------------------------------------------------------------------

  ! compute coefficient
  if ( h == 0.0d0 ) then
     const1= 0.0d0
  else
     const1= psi/(h*(1.0d0 + psi))
  end if

  ! --------------------------------
  ! compute tangent stiffness tensor:   (Eq. 54)
  ! --------------------------------
  do indx=1, nndex
     do jndx=1, nndex
	    do kndx=1, nndex
		   do lndx=1, nndex

              dij= krodelt(indx,jndx) ! get d_ij

              ! set tangent stiffness tensor component
		      tantens(indx,jndx,kndx,lndx)= elstens(indx,jndx,kndx,lndx) - const1* lptens(indx,jndx)*lptens(kndx,lndx) 

           end do
        end do
     end do
  end do



  return
end subroutine gettantensevp1



subroutine getkirjevp1(ematpro,nndex,deftens,pdirtens,delt,effstrs,effpstrn, tantens,kirjtens)
  !=======================================================================
  !  getkirjevp1 = compute jaumann rate of kirchoff stress at mid step
  !            
  !                 note:
  !                 ----
  !                 computer & structure, vol 18, 1984, pp. 875-887 (Eq. 3.9)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : constitutive model properties
  !
  !  nndex : dimension of tensor
  !
  !  deftens(nndex,nndex) : rate of deformation tensor
  !
  !  pdirtens(nndex,nndex) : plastic flow direction tensor
  !
  !  delt : time increment
  !
  !  effstrs : effective stress
  !
  !  effpstrn : effective plastic strain
  !
  !  output:
  !  ------
  !  tantens(nndex,nndex,nndex,nndex) : tangent modulus
  !
  !  kirjtens(nndex,nndex) : the jaumann rate of kirchoff stress (Eq. 53)
  !
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: deftens
  real(8), dimension(nndex,nndex), intent(in) :: pdirtens
  real(8), intent(in) :: delt
  real(8), intent(in) :: effstrs, effpstrn

  real(8), dimension(nndex,nndex,nndex,nndex), intent(out) :: tantens
  real(8), dimension(nndex,nndex), intent(out) :: kirjtens
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: e,nu,rho,eps0dot,m,sig0,eps0,n
  real(8) :: mu, lamd
  real(8), dimension(nndex,nndex,nndex,nndex) :: elstens

  real(8) :: h, psi
  real(8) :: eta

  real(8), dimension(nndex,nndex) :: lptens

  real(8) :: const1
  real(8), dimension(nndex,nndex) :: cdtens

  real(8) :: dij, krodelt

  ! loop index
  integer :: indx, jndx, kndx, lndx 
  ! ====================================

  ! initialize
  tantens(:,:,:,:)= 0.0d0
  kirjtens(:,:)= 0.0d0

  ! -------------------------------------------------------------------
  ! get material properties
  ! -----------------------
  call getpropevp1(ematpro, e,nu,rho,eps0dot,m,sig0,eps0,n)
     ! input : ematpro
     ! output : e,nu,rho,eps0dot,m,sig0,eps0,n

  ! get lame constant
  call getlameconst(e,nu, lamd,mu)
    ! input : e,nu
    ! output : lamd,mu

  ! compute elastic moduli tensor
  call getelstens(nndex,lamd,mu, elstens)
     ! input : nndex,lamd,mu
     ! output : elstens
  ! -------------------------------------------------------------------
  ! compute effective plastic strain rate: same as the value of eta function
  call getetaevp1(ematpro,effstrs,effpstrn, eta)
     ! input : ematpro,effstrs,effpstrn
     ! output : eta

  ! -------------------------------------------------------------------

  ! compute coefficient h and psi:   (Eq. 50, 52)
  call gethnpsievp1(0,ematpro,delt,effstrs,0.0d0,effpstrn, h,psi)
     ! input : 0(opt=inconsistent),ematpro,delt,effstrs,0.0d0,effpstrn
     ! output : h,psi
  ! -------------------------------------------------------------------

  ! compute "large" p tensor: p_ij= elas_ijkl * pdir_kl   (Eq. 2.9)
  call sprdtt(0,nndex,elstens,pdirtens, lptens)
     ! input : 0(opt=:),nndex,elstens,pdirtens
     ! output : lptens
  ! -------------------------------------------------------------------

  ! compute tangent stiffness tensor:   (Eq. 3.9)
  call gettantensevp1(ematpro,nndex,lptens,h,psi,effstrs, tantens)
     ! input : ematpro,nndex,lptens,h,psi,effstrs
     ! output : tantens
  ! -------------------------------------------------------------------

  ! calculate constant 
  const1= eta / (1.0d0+psi)

  ! compute scalar product: c^tan : D
  call sprdtt(0,nndex,tantens,deftens, cdtens)
     ! input : 0(opt),nndex,tantens,deftens
     ! output : cdtens
  ! -------------------------------------------------------------------


  ! ----------------------------------------
  ! compute jaumann rate of kirchhoff stress:   (Eq. 53)
  ! ----------------------------------------
  do indx=1, nndex
     do jndx=1, nndex
        dij= krodelt(indx,jndx) ! get d_ij
        kirjtens(indx,jndx)= cdtens(indx,jndx) - const1*lptens(indx,jndx)
     end do
  end do

  return
end subroutine getkirjevp1




subroutine updkirevp1(ematpro,nndex,middeftens,delt,effpstrn, kirtens, ctan,midkirtens,effstrs) 
  !=======================================================================
  !  updkirtevp1 = driver for updating kirchoff stress of evp model
  !
  !                kir_n+1= kir_n + (a kir_n+1/2 / a t)*delt
  !                where a kir_n+1/2 / a t = kir_n+1/2 ^jauman + w.kir_n + kir_n.w^t
  !
  !                note:
  !                ----
  !                D. Peirce, comp. & str., vol. 18, 1984, pp. 875-887
  !                S. Li, IJSS, vol. 39, 2002, pp. 1213-1240
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : constitutive model properties
  !
  !  nndex : dimension of tensor
  !
  !  middeftens(nndex,nndex) : rate of deformation tensor: mid step [n+1/2]
  !
  !  delt : time increment
  !
  !  effpstrn : effective plastic strain: previous step [n]
  !
  !
  !  inoutput:
  !  --------
  !  kirtens(nndex,nndex) : kirchoff stress tensor: update [n] -> [n+1]
  !
  !  output:
  !  --------
  !  ctan(nndex,nndex,nndex,nndex) : tangent modulus
  !
  !  midkirtens(nndex,nndex) : kirchoff stresstensor: mid step [n+1/2]
  !
  !  effstrs : effective stress
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: middeftens
  real(8), intent(in) :: delt
  real(8), intent(in) :: effpstrn

  real(8), dimension(nndex,nndex), intent(inout) :: kirtens

  real(8), dimension(nndex,nndex,nndex,nndex), intent(out) :: ctan
  real(8), dimension(nndex,nndex), intent(out) :: midkirtens
  real(8), intent(out) :: effstrs
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nndex,nndex) :: devstrstens
  real(8), dimension(nndex,nndex) :: pdirtens

  real(8), dimension(nndex,nndex) :: dkirtens, kirjtens
  ! ====================================
  ! initialize: do not initialize kirtens
  ctan(:,:,:,:)= 0.0d0
  midkirtens(:,:)= 0.0d0
  effstrs= 0.0d0

  ! -------------------------------------------------------------------
  ! compute deviatoric stress tensor: previous step [n]
  call getdevtens(nndex,kirtens, devstrstens)
    ! input : nndex,kirtens
    ! output : devstrstens

  ! compute effective stress: previous step [n]
  call geteffstrs(nndex,devstrstens, effstrs)
    ! input : nndex,devstrstens
    ! output : effstrs

  ! compute plastic strain flow direction tensor: previous step [n]
  call getj2dirtens(nndex,devstrstens,effstrs, pdirtens)
    ! input : nndex,devstrstens,effstrs
    ! output : pdirtens
  ! -------------------------------------------------------------------

  ! compute jauman rate of kirchoff stress: mid step [n+1/2]
  call getkirjevp1(ematpro,nndex,middeftens,pdirtens,delt,effstrs,effpstrn, ctan,kirjtens)
    ! input : ematpro,nndex,middeftens,pdirtens,delt,effstrs,effpstrn
    ! output : ctan,kirjtens
  ! -------------------------------------------------------------------

  ! compute increment of kirchoff stress: d_kir= delt * kir^jauman
  dkirtens(1:nndex,1:nndex)= delt* kirjtens(1:nndex,1:nndex)
  ! -------------------------------------------------------------------

  ! -----------------------
  ! compute kirchoff stress: mid step [n+1/2]
  ! -----------------------
  midkirtens(1:nndex,1:nndex)= kirtens(1:nndex,1:nndex) + 0.50d0 * dkirtens(1:nndex,1:nndex) ! at mid point [n+1/2]

  ! ----------------------------
  ! update kichoff stress tensor: current step [n+1]
  ! ----------------------------
  kirtens(1:nndex,1:nndex)= kirtens(1:nndex,1:nndex) + dkirtens(1:nndex,1:nndex) ! at current step [n+1]



  return
end subroutine updkirevp1





subroutine updpstrnevp1(ematpro,nndex,middeftens,prvkirtens,midkirtens,delt, effpstrn) 
  !=======================================================================
  !  updpstrnevp1 = update plastic strain
  !
  !                 tmp_n+1= tmp_n + (a tmp_n+1/2 / a t)*delt
  !                 where a tmp_n+1/2 / a t= (kai/rho*cp) * (effstrs_n+1/2)*(a effpstrn_n+1/2 / a t)
  !
  !                 note:
  !                 ----
  !                 D. Peirce, comp. & str., vol. 18, 1984, pp. 875-887
  !                 S. Li, IJSS, vol. 39, 2002, pp. 1213-1240
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : constitutive model properties
  !
  !  nndex : dimension of tensor
  !
  !  middeftens(nndex,nndex) : rate of deformation tensor: mid step [n+1/2]
  !
  !  prvkirtens(nndex,nndex) : kirchoff stresstensor: previous step [n]
  !
  !  midkirtens(nndex,nndex) : kirchoff stresstensor: mid step [n+1/2]
  !
  !  delt : time increment
  !
  !  inoutput:
  !  --------
  !  effpstrn : effective plastic strain: update [n] -> [n+1]
  !
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: middeftens
  real(8), dimension(nndex,nndex), intent(in) :: prvkirtens, midkirtens
  real(8), intent(in) :: delt

  real(8), intent(inout) :: effpstrn
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: e,nu,rho,eps0dot,m,sig0,eps0,n
  real(8) :: lamd,mu
  real(8), dimension(nndex,nndex,nndex,nndex) :: elstens

  real(8), dimension(nndex,nndex) :: devstrstens
  real(8) :: prveffstrs, mideffstrs

  real(8), dimension(nndex,nndex) :: pdirtens

  real(8), dimension(nndex,nndex) :: lptens
  real(8) :: lpd, sprdts

  real(8) :: eta, h,psi

  real(8) :: pstrndot
  ! ====================================
  ! initialize: do not initialize effpstrn


  ! -------------------------------------------------------------------
  ! get material properties
  ! -----------------------
  call getpropevp1(ematpro, e,nu,rho,eps0dot,m,sig0,eps0,n)
     ! input : ematpro
     ! output : e,nu,rho,eps0dot,m,sig0,eps0,n

  ! get lame constant
  call getlameconst(e,nu, lamd,mu)
    ! input : e,nu
    ! output : lamd,mu

  ! compute elastic moduli tensor
  call getelstens(nndex,lamd,mu, elstens)
     ! input : nndex,lamd,mu
     ! output : elstens
  ! -------------------------------------------------------------------

  ! compute deviatoric stress tensor: previous step [n]
  call getdevtens(nndex,prvkirtens, devstrstens)
     ! input : nndex,prvkirtens
     ! output : devstrstens

  ! compute effective stress: previous step [n]
  call geteffstrs(nndex,devstrstens, prveffstrs)
     ! input : nndex,devstrstens
     ! output : prveffstrs

  ! compute effective plastic strain rate: previous step [n]
  call getetaevp1(ematpro,prveffstrs,effpstrn, eta)
    ! input : ematpro,prveffstrs,effpstrn,
    ! output : eta
  ! -------------------------------------------------------------------

  ! compute deviatoric stress tensor: mid step [n+1/2]
  call getdevtens(nndex,midkirtens, devstrstens)
    ! input : nndex,midkirtens
    ! output : devstrstens

  ! compute effective stress: mid step [n+1/2]
  call geteffstrs(nndex,devstrstens, mideffstrs)
    ! input : nndex,devstrstens
    ! output : mideffstrs

  ! ---------------------------------
  ! compute plastic strain flow direction tensor: previous step [n+1/2]
  call getj2dirtens(nndex,devstrstens,mideffstrs, pdirtens)
    ! input : nndex,devstrstens,mideffstrs
    ! output : pdirtens

  ! compute "large" p tensor: p_ij= elas_ijkl * pdir_kl   (Eq. 51)
  call sprdtt(0,nndex,elstens,pdirtens, lptens)
    ! input : 0(opt=:),nndex,elstens,pdirtens
    ! output : lptens

  ! compute p:d
  lpd=sprdts(0,nndex,lptens,middeftens)
  ! ---------------------------------

  ! compute coefficient h and psi: mid step [n+1/2]   (Eq. 50, 52)
  call gethnpsievp1(1,ematpro,delt,prveffstrs,mideffstrs,effpstrn, h,psi)
    ! input : 1(opt=consistent),ematpro,delt,prveffstrs,mideffstrs,effpstrn
    ! output : h,psi

  ! compute effective plastic strain rate: mid step [n+1/2]   (Eq. 49)
  if ( h == 0.0d0 .or. psi==0.0d0 ) then
     pstrndot= 0.0d0
  else
     pstrndot= eta/(1.0d0+psi) + psi/( h*(1.0d0 + psi) )*lpd
  end if
  ! -------------------------------------------------------------------

  ! -------------------------------
  ! update effective plastic strain: [n] -> [n+1]
  ! -------------------------------
  effpstrn= effpstrn + dsqrt(pstrndot**2)*delt



  return
end subroutine updpstrnevp1



subroutine updevp1(ematpro,nndex,midltens,delt, effpstrn,kirtens, ctan,effstrs) 
  !=======================================================================
  !  updevp1 = main driver routine
  !            update elasto-visco plastic constitutive model and historic variables
  !
  !            note:
  !            ----
  !            D. Peirce, comp. & str., vol. 18, 1984, pp. 875-887
  !            S. Li, IJSS, vol. 39, 2002, pp. 1213-1240
  !
  !            note:
  !            ----
  !            kirchoff stress :  main update
  !            effpstrn : evolving internal variavles
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : constitutive model properties
  !
  !  nndex : dimension of tensor
  !
  !  midltens(nndex,nndex) : velocity gradient tensor at [n+1/2] : predicted value
  !
  !  delt : time increment
  !
  !  inoutput:
  !  --------
  !  effpstrn : effective plastic strain
  !
  !  kirtens(nndex,nndex) : kirchoff stress tensor
  !
  !  output:
  !  ------
  !  ctan(nndex,nndex,nndex,nndex) : tangent modulus
  !
  !  effstrs : effective stress
  !
  ! ======================================================================

  use preset
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: midltens
  real(8), intent(in) :: delt
  
  real(8), intent(inout) :: effpstrn
  real(8), dimension(nndex,nndex), intent(inout) :: kirtens

  real(8), dimension(nndex,nndex,nndex,nndex), intent(out) :: ctan
  real(8), intent(out) :: effstrs
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nndex,nndex) :: devstrstens
  real(8) :: eta

  real(8), dimension(nndex,nndex) :: middeftens
  real(8), dimension(nndex,nndex) :: prvkirtens ! kirchoff stress: previous step [n]
  real(8), dimension(nndex,nndex) :: midkirtens ! kirchoff stress: mid step [n+1/2]
  real(8) :: pstrndot ! consistently approximated effective plastic strain: mid step [n+1/2]  
  ! ====================================
  ! initialize : do not initialize effpstrn, kirtens and tmp
  ctan(:,:,:,:)= 0.0d0

  ! store kirchoff stress tensor: previous step [n]
  prvkirtens(:,:)= kirtens(:,:)

  ! compute rate of deformation tensor: mid step [n+1/2]
  call symtens(nndex,midltens, middeftens)
     ! input : nndex,midltens
     ! output : middeftens

  ! ----------------------
  ! update kirchoff stress: [n] -> [n+1]
  ! ----------------------
  call updkirevp1(ematpro,nndex,middeftens,delt,effpstrn, kirtens, ctan,midkirtens,effstrs) 
     ! input : ematpro,nndex,middeftens,delt,effpstrn
     ! inoutput: kirtens
     ! output : ctan,midkirtens,effstrs

  ! ---------------------
  ! update plastic strain: [n] -> [n+1]
  ! ---------------------
  call updpstrnevp1(ematpro,nndex,middeftens,prvkirtens,midkirtens,delt, effpstrn) 
     ! input : ematpro,nndex,middeftens,prvkirtens,midkirtens,delt
     ! inoutput : effpstrn



  return
end subroutine updevp1