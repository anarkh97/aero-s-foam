! ================================
! crack growth criterion
! ================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine           maxcri2d1             (prmcri,mgqpt1,crival,criang, nflagcr,crvalue,crangle)
! 2.  subroutine           chkmptcri2d           (optstr,svoit2d, crival,criang)
! 3.  subroutine           chkmptscri2d          (optstr,svoit2d, crival,criang)
! 4.  subroutine           getptstr2d            (optstr,svoit2d, pval1,pval2,pang1,pang2)
! 5.  subroutine           getpsstr2d            (optstr,svoit2d, pval1,pang1)
!
! =========================================================================================================



subroutine maxcri2d1(prmcri,mgqpt1,crival,criang, nflagcr,crvalue,crangle)
  !=======================================================================
  !  maxcri2d1= 
  !
  !            note:
  !            ----
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  prmcri(*) : criterion parameter
  !
  !  mgqpt1 : number of gq point
  !
  !  crival(mgqpt1) : pre computed crack criterion value
  !
  !  criang(mgqpt1) : pre computed crack propagation angle
  !
  !  output:
  !  ------
  !  nflagcr : crack initiation flag
  !            nflagcr=0 : crack criterion is not satisfied
  !            nflagcr=1 : crack criterion is satisfied
  !
  !  crvalue : crack propagation value
  !
  !  crangle : crack propagation angle
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: prmcri
  integer, intent(in) :: mgqpt1
  real(8), dimension(mgqpt1), intent(in) :: crival, criang

  integer, intent(out) :: nflagcr
  real(8), intent(out) :: crvalue
  real(8), intent(out) :: crangle
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: valcr

  ! loop index
  integer :: igaus
  ! ====================================

  ! initialize
  nflagcr= 0
  crvalue= 0.0d0
  crangle= 0.0d0


  ! -------------------------------------------
  ! get criterion parameters
  valcr= prmcri(5) ! critical value
  ! -------------------------------------------

  ! loop over gauss point
  do igaus=1, mgqpt1

     if ( crival(igaus) >= valcr ) then

        ! increase counter
        nflagcr= nflagcr + 1

        ! crack propagation
        crvalue= crvalue + crival(igaus)
        
        ! crack propagation angle
        crangle= crangle + criang(igaus)

     end if

  end do ! do igaus

  ! ------------------------------------------------------------------
  ! set flag and average angle
  if ( nflagcr /= 0 ) then

     crvalue= crvalue / real(nflagcr)

     crangle= crangle / real(nflagcr)

     nflagcr= 1

  end if
  ! ------------------------------------------------------------------



  return
end subroutine maxcri2d1





subroutine chkmptcri2d(optstr,svoit2d, crival,criang)
  !=======================================================================
  !  mptcri2d= 2d maximum principal tensile stress/strain criterion
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optstr : stress/strain option handler
  !           opt=0 : stress
  !           opt=1 : strain
  !
  !  svoit2d(3,1) : 2d voight form stress or strain 
  !
  !  output:
  !  ------
  !  crival : criterion value
  !
  !  criang : crack kink angle form criterion
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optstr
  real(8), dimension(3,1), intent(in) :: svoit2d

  real(8), intent(out) :: crival
  real(8), intent(out) :: criang
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: pval1, pval2
  real(8) :: pang1, pang2
  ! ====================================

  ! initialize
  crival= 0.0d0
  criang= 0.0d0


  ! compute principal tensile stress/strain
  call getptstr2d(optstr,svoit2d, pval1,pval2,pang1,pang2)
     ! input : optstr,svoit2d
     ! output : pval1,pval2,pang1,pang2

  ! set results
  ! maximum principal value
  crival= pval1

  ! perpendicular to maximum principal direction
  criang= pang2



  return
end subroutine chkmptcri2d





subroutine chkmptscri2d(optstr,svoit2d, crival,criang)
  !=======================================================================
  !  mptcri2d= 2d maximum principal tensile-shear stress/strain combined criterion
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optstr : stress/strain option handler
  !           opt=0 : stress
  !           opt=1 : strain
  !
  !  svoit2d(3,1) : 2d voight form stress or strain 
  !
  !  output:
  !  ------
  !  crival : maximum tensile stress / maximum shear stress
  !
  !  criang : crack kink angle : principal tensile stress criterion
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optstr
  real(8), dimension(3,1), intent(in) :: svoit2d

  real(8), intent(out) :: crival
  real(8), intent(out) :: criang
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: ptval1, ptval2, ptang1, ptang2
  real(8) :: psval1, psang1

  ! ====================================

  ! initialize
  crival= 0.0d0
  criang= 0.0d0


  ! compute principal tensile stress/strain
  call getptstr2d(optstr,svoit2d, ptval1,ptval2,ptang1,ptang2)
     ! input : optstr,svoit2d
     ! output : ptval1,ptval2,ptang1,ptang2

  ! compute principal shear stress/strain
  call getpsstr2d(optstr,svoit2d, psval1,psang1)
     ! input : optstr,svoit2d
     ! output : psval1,psang1


  ! ratio between principal maximum tensile stress and shear stress
  if ( psval1 /= 0.0d0 ) then
     crival= ptval1 / psval1
  else
     crival= 0.0d0

  end if

  ! set crack growth angle of 
  criang= ptang2



  return
end subroutine chkmptscri2d





subroutine getptstr2d(optstr,svoit2d, pval1,pval2,pang1,pang2)
  !=======================================================================
  !  getptstr2d = compute principal tensile stress/strain
  !
  !               note:
  !               ----
  !               minimum principal direction angle: -90 < pang(2) < 90
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optstr : stress/strain option handler
  !           opt=0 : stress
  !           opt=1 : strain
  !
  !  svoit2d(3,1) : 2d voight form stress or strain 
  ! 
  !  output:
  !  ------
  !  pval1= maximum principle tensile stress
  !
  !  pval2= minimum principle tensile stress
  !
  !  pang1= maximum principle tensile stress angle
  !
  !  pang2= minimum principle tensile stress angle
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optstr
  real(8), dimension(3,1), intent(in) :: svoit2d

  real(8), intent(out) :: pval1, pval2
  real(8), intent(out) :: pang1, pang2
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: s11, s12, s22
  ! ====================================

  ! initialize
  pval1= 0.0d0
  pval2= 0.0d0
  pang1= 0.0d0
  pang2= 0.0d0


  ! set component
  s11= svoit2d(1,1)
  s22= svoit2d(2,1)

  select case(optstr)
  case(0) ! stress
     s12= svoit2d(3,1)

  case(1) ! strain
     s12= 0.50d0 * svoit2d(3,1)

  case default
     write(*,*) "wrong opt: getptstr2d"
     write(nout5,*) "wrong opt: getptstr2d"
     stop

  end select

  ! compute principal stress/strain
  pval1= ( s11 + s22 ) / 2.0d0 + dsqrt( ( ( s11 - s22 ) / 2.0d0 )**2 + s12**2 ) ! maximum
  pval2= ( s11 + s22 ) / 2.0d0 - dsqrt( ( ( s11 - s22 ) / 2.0d0 )**2 + s12**2 ) ! minimum

  ! -----------------------------------
  ! note: description of generic ATAN2
  ! ----
  !
  ! The result type is the same as x. The value lies in the range -pi < ATAN2 (y, x) <= pi.
  ! If x /= zero, the result is approximately equal to the value of arctan (y/x). 
  ! -----------------------------------
  ! If y > zero, the result is positive. 
  ! If y < zero, the result is negative. 
  ! If y = zero, the result is zero (if x > zero) or pi (if x < zero). 
  ! If x = zero, the absolute value of the result is pi/2
  ! -----------------------------------

  ! maximum principal direction
  pang1= 0.50d0 * datan2( 2.0d0*s12,s11-s22 )

  ! minimum principal direction: -90 < pang2 <= 90
  ! pang1= 0 or (+/-)180
  if ( pang1 == 0.0d0 .or. pang1 >= 180.0d0*d2r .or. pang1 <= -180.0d0*d2r ) then
      pang2= 90.0d0 * d2r

  ! pang1=(+/-)90
  else if ( pang1 == 90.0d0*d2r .or. pang1 == -90.0d0*d2r ) then
      pang2= 0.0d0

  ! 0 < pang1 < 90
  else if ( 0.0d0 < pang1 .and. pang1 < 90.0d0*d2r ) then
      pang2= pang1 - 90.0d0 * d2r

  ! 90 < pang1 < 180
  else if ( 90.0d0*d2r < pang1 .and. pang1 < 180.0d0*d2r ) then
      pang2= pang1 - 90.0d0 * d2r

  ! -90 < pang1 < 0
  else if ( -90.0d0*d2r < pang1 .and. pang1< 0.0d0 ) then
      pang2= pang1 + 90.0d0 * d2r

  ! -180 < pang1 < -90
  else if ( -180.0d0*d2r < pang1 .and. pang1< -90.0d0*d2r ) then
      pang2= pang1 + 90.0d0 * d2r

  end if


  
  return
end subroutine getptstr2d





subroutine getpsstr2d(optstr,svoit2d, pval1,pang1)
  !=======================================================================
  !  getpsstr2d = compute principal shear stress/strain
  !
  !               note:
  !               ----
  !               minimum principal direction angle: -90 < pang(2) < 90
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optstr : stress/strain option handler
  !           opt=0 : stress
  !           opt=1 : strain
  !
  !  svoit2d(3,1) : 2d voight form stress or strain 
  ! 
  !  output:
  !  ------
  !  pval1= maximum principle tensile stress
  !
  !  pang1= maximum principle tensile stress angle
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optstr
  real(8), dimension(3,1), intent(in) :: svoit2d

  real(8), intent(out) :: pval1
  real(8), intent(out) :: pang1
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: s11, s12, s22
  ! ====================================

  ! initialize
  pval1= 0.0d0
  pang1= 0.0d0


  ! set component
  s11= svoit2d(1,1)
  s22= svoit2d(2,1)

  select case(optstr)
  case(0) ! stress
     s12= svoit2d(3,1)

  case(1) ! strain
     s12= 0.50d0 * svoit2d(3,1)

  case default
     write(*,*) "wrong opt: getptstr2d"
     write(nout5,*) "wrong opt: getptstr2d"
     stop

  end select

  ! compute principal stress/strain
  pval1= dsqrt( ( ( s11 - s22 ) / 2.0d0 )**2 + s12**2 ) ! maximum

  ! -----------------------------------
  ! note: description of generic ATAN2
  ! ----
  !
  ! The result type is the same as x. The value lies in the range -pi < ATAN2 (y, x) <= pi.
  ! If x /= zero, the result is approximately equal to the value of arctan (y/x). 
  ! -----------------------------------
  ! If y > zero, the result is positive. 
  ! If y < zero, the result is negative. 
  ! If y = zero, the result is zero (if x > zero) or pi (if x < zero). 
  ! If x = zero, the absolute value of the result is pi/2
  ! -----------------------------------

  ! maximum principal direction
  pang1= 0.50d0 * datan2( s11-s22, -2.0d0*s12 )


  
  return
end subroutine getpsstr2d