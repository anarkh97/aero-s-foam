! ==================================
! general geometry
! ==================================
!      type                  name                              arguement
!      ----                  ----                              ---------
! 1.  subroutine         parline                          (ndime,line,tpar, val)       ! parametric line eqn
! 2.  subroutine         getcrospt3dll                    (line1,line2, crospt)        ! line-line cross point
! 3.  subroutine         getcrospt3dpl                    (ppoin,pnvec,line2, crospt)  ! plane-line cross point
! 4.  subroutine         getline3dpar                     (line, c1,dc1,c2,dc2,c3,dc3) ! 3d line equation
! 5.  subroutine         getplane3dpn                     (ppoin,pnvec, c1,c2,c3,c4)   ! 3d plane equation
! ---------------------------------------------------------------------------------------------
! 6.  real(8) function   normv                            (opt,ndime, avec)
! 7.  subroutine         l2lmapping                       (ndime,x0,x0pt,x1, x)
! 8.  logical function   ptonlseg                         (ndime,segment,point)
! 9.  subroutine         elecntr                          (ndime,nnode,ecord, ecntrpt)
! 10. subroutine         divdln                           (optdvd,ndime,rm,rn,points, ptdivd)
! 11. subroutine         getsubtri2d0                     (mtri,npt,ptcord, ntri,triconc,triarea)  ! construct sub tri.
! ---------------------------------------------------------------------------------------------
! 12. real(8) function   distpt2ele                       (opt,ndime,nnode,pt,ecord)        ! point-element distance
! 13. real(8) function   distele                          (opt,ndime,nnode,ecord1,ecord2)   ! element-element distance 
! 14. real(8) function   distpt                           (opt,ndime, points)               ! point-point distance 
! 15. real(8) function   distpts                          (opt,ndime,point1,point2)         ! point-point distance 
! 16. subroutine         getdistptl2d0                    (lpts,pt, dist)                   ! line-point distance
! 17. subroutine         getdistptl2d1                    (lpt,langle,pt, dist)             ! line-point distance
! ---------------------------------------------------------------------------------------------
! 18. subroutine         getadjhavg0                      (ipoin,mnode,npoin,eletyp,conec,coord,nod2ele, havg)
! 19. subroutine         geoimpcyl                        (impsize,impmode,ndime,npoin, coord) ! geometric imperfection
!
! =========================================================================================================



subroutine parline(ndime,line,tpar, val)
  !=======================================================================
  !  parline= compute parametric line equation and corresponding point
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of given points
  !
  !  line(ndime,2) : start and end point coordinates of line
  !
  !  tpar : parametric value t
  !
  !  output:
  !  ------
  !  val(ndime,1) : coordinate corresponding to given parametric value t
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,2), intent(in) :: line
  real(8), intent(in) ::tpar

  real(8), dimension(ndime,1), intent(out) :: val
  ! ====================================
  ! local variable
  ! ==============

  ! ====================================

  ! initialize
  val(:,:)= 0.0d0


  ! compute point value
  val(1:ndime,1)= line(1:ndime,1) + ( line(1:ndime,2) - line(1:ndime,1) ) * tpar



  return
end subroutine parline





subroutine getcrospt3dll(line1,line2, crospt)
  !=======================================================================
  !  getcrospt3dll= compute intersection point between line and line
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  line1(3,2) : arbitarary points on 3d line1
  !
  !  line2(3,2) : arbitarary points on 3d line2
  !
  !  output:
  !  ------
  !  crospt(3,1) : cross point
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,2), intent(in) :: line1, line2

  real(8), dimension(3,1), intent(out) :: crospt
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,1) :: avec, bvec, cvec
  real(8), dimension(3,1) :: abvec, cbvec
  real(8) :: const1, dotprdt
  real(8) :: const2, normv

  ! loop index
  ! ====================================

  ! initialize
  crospt(:,:)= 0.0d0

  avec(1:3,1)= line1(1:3,2)-line1(1:3,1)
  bvec(1:3,1)= line2(1:3,2)-line2(1:3,1)
  cvec(1:3,1)= line2(1:3,1)-line1(1:3,1)

  ! compute cross product: cvec x bvec
  call crsprdt3d(0,cvec,bvec, cbvec)
     ! input : 0(opt),cvec,bvec
     ! output : cbvec

  ! compute cross product: avec x bvec
  call crsprdt3d(0,avec,bvec, abvec)
     ! input : 0(opt),avec,bvec
     ! output : abvec

  ! (cvec x bvec) . (avec x bvec)
  const1= dotprdt(3,cbvec,abvec)

  ! |avec x bvec|^2
  const2= normv(1,3, abvec)

  ! compute cross point
  crospt(1:3,1)= line1(1:3,1) + avec(1:3,1)*(const1/const2)


  return
end subroutine getcrospt3dll





subroutine getcrospt3dpl(ppoin,pnvec,line2, crospt)
  !=======================================================================
  !  getcrospt3dpl= compute intersection point between plane and line
  !
  !                 note:
  !                 ----
  !                 A. Bowyer and J. Woodwark, a programmer's geometry pp. 111
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ppoin(3,1) : arbitrary point on plane
  !
  !  pnvec(3,1) : plane normal vector
  !
  !  line2(3,2) : arbitarary points on 3d line2
  !
  !  output:
  !  ------
  !  crospt(3,1) : cross point
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,1), intent(in) :: ppoin
  real(8), dimension(3,1), intent(in) :: pnvec
  real(8), dimension(3,2), intent(in) :: line2

  real(8), dimension(3,1), intent(out) :: crospt
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: a, b, c, d
  real(8) :: x0, f, y0, g, z0, h

  real(8) :: denom, t, crosx, crosy, crosz
  ! ====================================

  ! initialize
  crospt(:,:)= 0.0d0

  ! get 3d plane equation coefficient
   call getplane3dpn(ppoin,pnvec, a,b,c,d)
     ! input : ppoin,pnvec
     ! output : a,b,c,d

  ! get 3d line equation coefficient in parameteric form
  call getline3dpar(line2, x0,f,y0,g,z0,h)
     ! input : line2
     ! output : x0,f,y0,g,z0,h

  denom= a*f + b*g + c*h

  if( abs(denom) <=  toler(1) ) then
     write(*,*) "plane and line are parallel: getcrospt3dpl"
     write(nout5,*) "plane and line are parallel: getcrospt3dpl"

  else
     t= -(a*x0 + b*y0 + c*z0 + d) / denom
     crosx= x0 + f*t
     crosy= y0 + g*t
     crosz= z0 + h*t

  end if

  ! set result
  crospt(1,1)= crosx 
  crospt(2,1)= crosy 
  crospt(3,1)= crosz 



  return
end subroutine getcrospt3dpl





subroutine getline3dpar(line, c1,dc1,c2,dc2,c3,dc3)
  !=======================================================================
  !  getplane3dpn= get 3d line equation coefficient in paramter form
  !
  !                note:
  !                ----
  !                x= c1 x + dc1 t
  !                y= c2 y + dc2 t
  !                z= c3 z + dc3 t
  ! 
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  line(3,2) : start and end points of line segment
  !
  !  output:
  !  ------
  !  c1,dc1,c2,dc2,c3,dc3 : 3d line equation coefficient
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,2), intent(in) :: line

  real(8), intent(out) :: c1,dc1,c2,dc2,c3,dc3
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,1) :: lvec

  ! loop index
  ! ====================================

  ! initialize
  c1= 0.0d0
  dc1= 0.0d0
  c2= 0.0d0
  dc2= 0.0d0
  c3= 0.0d0
  dc3= 0.0d0

  !                 lvec
  ! --------o------------------>o----------
  !        poin1               poin2
  !         t=0                 t=1

  ! get line vector
  lvec(1:3,1)= line(1:3,2) - line(1:3,1)

  ! compute line equation coefficient
  c1= line(1,1) 
  dc1= lvec(1,1)

  c2= line(2,1)
  dc2= lvec(2,1) 

  c3= line(3,1)
  dc3= lvec(3,1) 


  return
end subroutine getline3dpar




subroutine getplane3dpn(ppoin,pnvec, c1,c2,c3,c4)
  !=======================================================================
  !  getplane3dpn= get 3d plane equation coefficient from origin and normal vector
  !
  !                note:
  !                ----
  !                plane equation: c1 x + c2 y + c3 z + c4 = 0
  ! 
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ppoin(3,1) : arbitrary point on plane
  !
  !  pnvec(3,1) : plane normal vector
  !
  !  output:
  !  ------
  !  c1,c2,c3,c4 : 3d plane equation coefficient
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,1), intent(in) :: ppoin
  real(8), dimension(3,1), intent(in) :: pnvec

  real(8), intent(out) :: c1,c2,c3,c4
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,1) :: pnvec0


  ! loop index
  ! ====================================

  ! initialize
  c1= 0.0d0
  c2= 0.0d0
  c3= 0.0d0
  c4= 0.0d0

  ! normalize plane normal vector
  call unitvec2(3,pnvec, pnvec0)
     ! input : 3(ndime),pnvec
     ! output : pnvec0

  ! compute plane equation coefficient
  c1= pnvec0(1,1) 
  c2= pnvec0(2,1) 
  c3= pnvec0(3,1) 
  c4= -( pnvec0(1,1)*ppoin(1,1) +  pnvec0(2,1)*ppoin(2,1) +  pnvec0(3,1)*ppoin(3,1) ) 



  return
end subroutine getplane3dpn





real(8) function normv(opt,ndime, avec)
  !=======================================================================
  !  distance = compute distance between two points
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt= option handler
  !       opt=0 : ||avec||= sqrt( avec . avec )
  !       opt=1 : ||avec||^2= avec . avec
  !
  !  ndime : dimension
  !
  !  line(ndime,2) : start and end point of line 
  !
  !  output:
  !  ------
  !  length : length of line
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opt, ndime
  real(8), dimension(ndime,1), intent(in) :: avec
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: normv0

  ! loop index
  integer :: idime
  ! ====================================

  ! initialize
  normv= 0.0d0

  ! compute norm
  normv0= 0.0d0 ! initialize
  do idime=1, ndime
     normv0= normv0 + avec(idime,1)**2
  end do

  select case(opt)
  case (0)
     normv= dsqrt(normv0)

  case(1)
     normv= normv0

  case default
    write(*,*) "option is not available: normv"
    write(nout5,*) "option is not available: normv"
    stop

  end select



  return
end function normv




subroutine l2lmapping(ndime,x0,x0pt,x1, x)
  !=======================================================================
  !  l2lmapping=line to line linear mapping
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension
  !
  !  x0(ndime,2) : line coordinate
  !
  !  x0pt(ndime,1) : point
  !
  !  x1(ndime,2) : line coordinate
  !
  !  output: 
  !  ------
  !  x(ndime,1) : mapped point  
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,2), intent(in) :: x0
  real(8), dimension(ndime,1), intent(in) :: x0pt
  real(8), dimension(ndime,2), intent(in) :: x1

  real(8), dimension(ndime,1), intent(out) :: x  
  ! ====================================
  ! local variable
  ! ==============
  logical :: ptonlseg   
  real(8) :: t0, t1, t

  ! loop index
  integer :: idime
  ! ====================================

  ! initialize
  x(:,:)= 0.0d0

  !                x0pt
  !  o--------------|--------------o
  ! x01                           x02
  !  |--------------|
  !         t1
  !  |-----------------------------|
  !                t0

  !                 x
  !  o--------------|--------------o
  ! x1                            x2

  if ( ptonlseg(ndime,x0,x0pt) ) then

     ! compute parameter
     t0= 0.0d0 ! initialize
     t1= 0.0d0
     do idime=1, ndime
        t0= t0 + (x0(idime,1) - x0(idime,2))**2
        t1= t1 + (x0(idime,1) - x0pt(idime,1))**2
     end do
     t= dsqrt(t1/t0)

     ! get x point by parametric line equation
     do idime=1, ndime
        x(idime,1)= x1(idime,1) + ( x1(idime,2) - x1(idime,1) ) * t
     end do

  else
    write(*,*) "x0pt is not on x0 line: l2lmapping"
    write(nout5,*) "x0pt is not on x0 line: l2lmapping"
    stop

  end if



  return
end subroutine l2lmapping





logical function ptonlseg(ndime,segment,point)
  !=======================================================================
  !  ptonlseg = check wether point is on the line segment or not
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  line(ndime,2) : strar point of finite length line
  !  point(ndime,1): end point of finite length line
  ! 
  !  output:
  !  ------
  !  chkinline2d : logical value
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,2), intent(in) :: segment
  real(8), dimension(ndime,1), intent(in) :: point
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,2) :: segment1, segment2
  real(8) :: distpt, dist12, dist1, dist2

  ! loop index
  ! ==================================== 

  ! initialize
  ptonlseg=.false.

  ! o---------------X-----------------o
  ! node1          point             node2
  !
  ! |----------------|----------------|
  !       dist1            dist2
  ! |---------------------------------|
  !               dist12
 
  ! note:
  ! ----
  ! for computational accuracy, we used length^2 for comparision

  ! compute distance12
  dist12= distpt(0,ndime, segment)

  ! compute distance 1
  segment1(1:ndime,1)= segment(1:ndime,1) 
  segment1(1:ndime,2)= point(1:ndime,1)
  dist1= distpt(0,ndime, segment1)

  ! compute distance 2
  segment2(1:ndime,1)= segment(1:ndime,2) 
  segment2(1:ndime,2)= point(1:ndime,1)
  dist2= distpt(0,ndime, segment2)

  !-------------------------------------------

  ! check result for both of node
  if ( abs(1.0d0 - (dist1+dist2)/dist12 ) <= toler(1) ) then
     ptonlseg= .true. ! point is in segment
  else
     ptonlseg= .false.
  end if



  return
end function ptonlseg





subroutine elecntr(ndime,nnode,ecord, ecntrpt)
  !=======================================================================
  !  elecntr = compute cneter of element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension
  !
  !  nnode : number of node
  !
  !  ecord(ndime,nnode) : element nodal coordinate
  !
  !  output:
  !  ------
  !  ecntrpt(ndime,1) : center point of element
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, nnode
  real(8), dimension(ndime,nnode), intent(in) :: ecord

  real(8), dimension(ndime,1), intent(out) :: ecntrpt
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: inode, idime
  ! ====================================

  ! initialize
  ecntrpt(:,:)= 0.0d0


  ! get center of elements
  do inode=1, nnode
     do idime=1, ndime
        ecntrpt(idime,1)= ecntrpt(idime,1) + ecord(idime,inode)/ real(nnode)
     end do
  end do



  return
end subroutine elecntr





subroutine divdln(optdvd,ndime,rm,rn,points, ptdivd)
  !=======================================================================
  !  divdln = divide given line according to given ratio
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optdvd : divide option
  !           optdvd= 0 : internally divide given line
  !           optdvd= 1 : externally divide given line
  !
  !  ndime : dimension of given line
  !
  !  rm, rn : divide ratio (m:n)
  !
  !  points(ndime,2) : given line
  !
  !  output:
  !  ------
  !  ptdivd(ndime,2) : result
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optdvd
  integer, intent(in) :: ndime
  real(8), intent(in) :: rm, rn
  real(8), dimension(ndime,2), intent(in) :: points
  
  real(8), dimension(ndime,1), intent(out) :: ptdivd
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: idime 
  ! ====================================

  ! initialize
  ptdivd(:,:)= 0.0d0


!  internally divide
!
!            m          :      n
!  o--------------------x-----------o
!  x1                   x           x2
!
!  x= ( m*x2 + n*x1 ) / ( m + n )



  select case(optdvd)
  case(0) ! internalli divide
     do idime=1, ndime
        ptdivd(idime,1)= ( rm * points(idime,2) + rn * points(idime,1) ) / (rm + rn)
     end do

  case(1) ! externally divide
     do idime=1, ndime
        ptdivd(idime,1)= ( rm * points(idime,2) - rn * points(idime,1) ) / (rm - rn)
     end do

     write(*,*) "not implemented yet: divdln"
     write(nout5,*) "not implemented yet: divdln"
     stop

  case default
     write(*,*) "not implemented yet: divdln"
     write(nout5,*) "not implemented yet: divdln"
     stop

  end select


  return
end subroutine divdln





subroutine getsubtri2d0(mtri,npt,ptcord, ntri,triconc,triarea)
  !=======================================================================
  !  getsubtri2d0 = construct sub triangle
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  mtri : maximum number of triangle
  !
  !  npt : the total number of point
  !
  !  ptcord(2,npt) : coordinate
  !
  !  output:
  !  ------
  !  ntri : the total number of sub triangle
  !
  !  triconc(3:mtri) : connectivity of sub triangle
  !
  !  triarea(mtri) : area of sub triangle
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: mtri
  integer, intent(in) :: npt
  real(8), dimension(2,npt), intent(in) :: ptcord

  integer, intent(out) :: ntri
  integer, dimension(3,mtri), intent(out) :: triconc
  real(8), dimension(mtri), intent(out) :: triarea
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(npt) :: x, y

  integer, parameter :: nrow = 6 ! 6 : neighbour element, 9: constraint curve
  real(8), dimension(npt) :: ds ! work space
  integer :: ier ! error handler

  integer :: nt
  integer :: lnew
  integer, dimension(npt) :: lend
  integer, dimension(6*npt) :: list, lptr
  integer, dimension(nrow,2*npt) :: ltri
  integer, dimension(2*npt) :: nodes

  real(8) :: x1,y1,x2,y2,x3, y3
  logical :: ratio
  real(8) :: cx,cy,cr,sa,ar

  real(8) :: sumarea

  ! dummy array for constrain curve
  integer :: ncc = 0
  integer, dimension(1) :: lcc = (/0/)
  integer :: lct(1)

  ! loop index
  integer :: ipt, it
  ! ====================================

  ! initialize
  ntri= 0
  triconc(:,:)= 0
  triarea(:)= 0.0d0


  ! copy coordinate
  do ipt=1, npt
     x(ipt)= ptcord(1,ipt)
     y(ipt)= ptcord(2,ipt)
  end do

  
  ! delaunay triangulation
  call trmesh(npt,x,y, list,lptr,lend,lnew,nodes(1:npt),nodes(npt+1:2*npt),ds,ier)
     ! input : npt,x,y
     ! output : list,lptr,lend,lnew,nodes(1:n),nodes(n+1:2*n),ds,ier

  ! check error
  if ( ier /= 0 ) then
     write(*,*) "error in trmesh: getsubtri2d0"
     write(nout5,*) "error in trmesh: getsubtri2d0"
     stop
  end if
  
  ! converts a triangulation data structure to a triangle list
  call trlist(ncc,lcc,npt,list,lptr,lend,nrow, nt,ltri,lct,ier)
     ! input : ncc,lcc,npt,list,lptr,lend,nrow
     ! output : nt,ltri,lct,ier

  ! check error
  if ( ier /= 0 ) then
     write(*,*) "error in trlist: getsubtri2d0"
     write(nout5,*) "error in trlist: getsubtri2d0"
     stop
  end if

  ! check size of output array
  if ( nt > mtri ) then
     write(*,*) "increase mtri: getsubtri2d0"
     write(nout5,*) "increase mtri: getsubtri2d0"
     stop
  end if


  ! set results
  ! -----------
  ntri= 0 ! initialize
  sumarea= 0.0d0
  do it=1, nt

     ! get coordinate
     x1= x(ltri(1,it))
     y1= y(ltri(1,it))
     x2= x(ltri(2,it))
     y2= y(ltri(2,it))
     x3= x(ltri(3,it))
     y3= y(ltri(3,it))

     ! check triangle
     call circum(x1,y1,x2,y2,x3,y3, ratio,cx,cy,cr,sa,ar)
        ! input : x1,y1,x2,y2,x3,y3
        ! output : ratio,cx,cy,cr,sa,ar

     ! skip boundary null triangle
     if ( abs(sa) <= toler(3) ) cycle
     
     ! increase counter
     ntri= ntri + 1

     ! store connectivity
     triconc(1:3,ntri)= ltri(1:3,it)

     ! store area
     triarea(ntri)= sa
     sumarea= sumarea + sa

   end do

   ! check error
!   if ( abs( sumarea - 0.50d0 ) / 0.50d0 > toler(3) ) then
!      write(*,*) "error in summation of the sub triangle area: getsubtri2d0"
!      write(nout5,*) "error in summation of the sub triangle area: getsubtri2d0"
!      stop
!   end if



  return
end subroutine getsubtri2d0





real(8) function distpt2ele(opt,ndime,nnode,pt,ecord)
  !=======================================================================
  !  distpt2ele = compute distance between point and element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt= option handler
  !       opt=0 : return length
  !       opt=1 : return length^2
  !
  !  ndime : dimension
  !
  !  nnode : number of node
  !
  !  pt(ndime,1) : point coordinate
  !
  !  ecord(ndime,nnode) : element nodal coordinate
  !
  !  output:
  !  ------
  !  distpt2ele : distance between center of two elements
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opt, ndime, nnode
  real(8), dimension(ndime,1), intent(in) :: pt
  real(8), dimension(ndime,nnode), intent(in) :: ecord
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,2) :: points
  real(8) :: distpt

  ! loop index
  integer :: inode, idime
  ! ====================================

  ! initialize
  distpt2ele= 0.0d0

  ! set given point coordinate
  points(1:ndime,1)= pt(1:ndime,1)

  ! get center of elements
  points(:,:)= 0.0d0 ! initialize
  do inode=1, nnode
     do idime=1, ndime
        points(idime,2)= points(idime,2) + ecord(idime,inode)/real(nnode)
     end do
  end do


  ! get distacne
  distpt2ele= distpt(opt,ndime, points)



  return
end function distpt2ele





real(8) function distele(opt,ndime,nnode,ecord1,ecord2)
  !=======================================================================
  !  distele = compute distance between two elements
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt= option handler
  !       opt=0 : return length
  !       opt=1 : return length^2
  !
  !  ndime : dimension
  !
  !  nnode : number of node
  !
  !  ecord1(ndime,nnode), ecord2(ndime,nnode) : element nodal coordinate
  !
  !  output:
  !  ------
  !  distele : distance between center of two elements
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opt, ndime, nnode
  real(8), dimension(ndime,nnode), intent(in) :: ecord1, ecord2
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,2) :: points
  real(8) :: distpt

  ! loop index
  integer :: inode, idime
  ! ====================================

  ! initialize
  distele= 0.0d0


  ! get center of elements
  points(:,:)= 0.0d0 ! initialize
  do inode=1, nnode
     do idime=1, ndime
        points(idime,1)= points(idime,1) + ecord1(idime,inode)/real(nnode)
        points(idime,2)= points(idime,2) + ecord2(idime,inode)/real(nnode)
     end do
  end do


  ! get distacne
  distele= distpt(opt,ndime, points)



  return
end function distele





real(8) function distpt(opt,ndime, points)
  !=======================================================================
  !  distpt = compute distance between two points
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt= option handler
  !       opt=0 : return length
  !       opt=1 : return length^2
  !
  !  ndime : dimension
  !
  !  points(ndime,2) : start and end point of line 
  !
  !  output:
  !  ------
  !  distpt : distance between two points
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opt, ndime
  real(8), dimension(ndime,2), intent(in) :: points
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: distance0

  ! loop index
  integer :: idime
  ! ====================================

  ! initialize
  distpt= 0.0d0

  ! compute line length
  distance0= 0.0d0
  do idime=1, ndime
     distance0= distance0 + ( points(idime,2) - points(idime,1) )**2
  end do

  select case(opt)
  case (0)
     distpt= dsqrt(distance0)

  case(1)
     distpt= distance0

  end select



  return
end function distpt





real(8) function distpts(opt,ndime,point1,point2)
  !=======================================================================
  !  distpts = compute distance between two points
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt= option handler
  !       opt=0 : return length
  !       opt=1 : return length^2
  !
  !  ndime : dimension
  !
  !  point1(ndime,1), point2(ndime,1) : start and end point of line 
  !
  !  output:
  !  ------
  !  distpts : distance between two points
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opt, ndime
  real(8), dimension(ndime,1), intent(in) :: point1
  real(8), dimension(ndime,1), intent(in) :: point2
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: distance0

  ! loop index
  integer :: idime
  ! ====================================

  ! initialize
  distpts= 0.0d0

  ! compute line length
  distance0= 0.0d0
  do idime=1, ndime
     distance0= distance0 + ( point1(idime,1) - point2(idime,1) )**2
  end do

  select case(opt)
  case (0)
     distpts= dsqrt(distance0)

  case(1)
     distpts= distance0

  end select



  return
end function distpts





subroutine getdistptl2d0(lpts,pt, dist)
  !=======================================================================
  !  getdistptl2d0 = compute distance from the given point to line
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
  real(8), dimension(2,2), intent(in) :: lpts
  real(8), dimension(2,1), intent(in) :: pt

  real(8), intent(out) :: dist
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: x1, y1
  real(8) :: x2, y2
  real(8) :: x0, y0
  
  real(8) :: term1, term2

  ! ====================================

  ! initialize
  dist= 0.0d0
  
  !                (x2,y2)    
  !                -
  !              -
  !            -
  !          *
  !        -    *   dist
  !      -         *
  !    -              * 
  !  -                   o (x0,y0)
  ! (x1,y1)
  
  
  ! set points on line 
  x1= lpts(1,1)
  y1= lpts(2,1)
  x2= lpts(1,2)
  y2= lpts(2,2)
  
  ! set point
  x0= pt(1,1)
  y0= pt(2,1)

  
  ! compute term  
  
  term1= abs( (x2-x1)*(y1-y0) - (x1-x0)*(y2-y1) )
  term2= sqrt( (x2-x1)**2 + (y2-y1)**2 )

  ! compute distance  
  if ( term2 /= 0.0d0 ) then

     dist= term1 / term2
 
  else
     write(*,*) "lpts is wrong: getdistptl2d0"
     write(nout5,*) "lpts is wrong: getdistptl2d0"
     stop
  
  end if



  return
end subroutine getdistptl2d0





subroutine getdistptl2d1(lpt,langle,pt, dist)
  !=======================================================================
  !  getdistptl2d1 = compute distance from the given point to line
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
  real(8), dimension(2,1), intent(in) :: lpt
  real(8), intent(in) :: langle
  real(8), dimension(2,1), intent(in) :: pt

  real(8), intent(out) :: dist
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: x1, y1
  real(8) :: x2, y2
  real(8) :: x0, y0
  
  real(8), dimension(2,2) :: lpts

  ! ====================================

  ! initialize
  dist= 0.0d0
  
  !                (x2,y2)    
  !                -
  !              -
  !            -
  !          *
  !        -    *   dist
  !      -         *
  !    -              * 
  !  -                   o (x0,y0)
  ! (x1,y1)
  
  
  ! set points on line 
  lpts(1,1)= lpt(1,1)
  lpts(2,1)= lpt(2,1)
  
  lpts(1,2)= lpt(1,1) + dcos(langle)
  lpts(2,2)= lpt(2,1) + dsin(langle)
  
  ! compute distance
  call getdistptl2d0(lpts,pt, dist)
     ! input : lpts,pt
     ! output : dist

 

  return
end subroutine getdistptl2d1





subroutine getadjhavg0(ipoin,mnode,npoin,eletyp,conec,coord,nod2ele, havg)
  !=======================================================================
  !  getadjhavg0 = compute average element size
  !
  !                note:
  !                ----
  !                average adjacent elements
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  -----
  !                          
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ipoin
  integer, intent(in) :: mnode,npoin
  integer, dimension(2,*), intent(in) :: eletyp
  integer, dimension(mnode,*), intent(in) :: conec
  real(8), dimension(2,*), intent(in) :: coord
  integer, dimension(npoin,*), intent(in) :: nod2ele

  real(8), intent(out) :: havg
  ! ====================================
  ! local variable
  ! ==============
  integer :: adjele

  integer :: optele
  integer :: nnode
  
  real(8) :: ehleng

  ! loop index
  integer :: iadjele, inode
  ! ====================================

  ! initialize
  havg= 0.0d0


  ! loop over adjacent elements
  do iadjele=1, nod2ele(ipoin,10)
     
     ! get adjcent element
     adjele= nod2ele(ipoin,iadjele)

     ! get element type and total number of nodes
     optele= eletyp(1,adjele)
     nnode= eletyp(2,adjele)

     ! compute element characteristic length
     call elehleng1a(adjele,optele,mnode,2,nnode,conec,coord, ehleng)
        ! input : adjele,optele,mnode,2(ndime),nnode,conec,coord
        ! output : ehleng

     ! sum on
     havg= havg + ehleng

  end do


  ! take average
  havg= havg / real(nod2ele(ipoin,10))



  return
end subroutine getadjhavg0





subroutine geoimpcyl(impsize,impmode,ndime,npoin, coord)
  !=======================================================================
  !  geoimpcyl = geometiric imperfection in cylinder
  !
  !              note:
  !              ----
  !              x-direction should be aligned with cylinder axis
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  inoutput:
  !  --------
  !  coord(ndime,npoin) : modified nodal coordinates
  !
  ! ======================================================================
  
  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: impsize, impmode
  integer, intent(in) :: ndime, npoin
 
   real(8), dimension(ndime,npoin), intent(inout) :: coord
  ! ====================================
  ! local variables
  ! ===============
  real(8) :: x0, y0, z0
  real(8) :: theta
  real(8) :: alpha

  ! loop index
  integer :: ipoin
  ! ====================================

  ! initialize: do not initialize
  
  
  do ipoin= 1, npoin
  
     x0= coord(1,ipoin)
     y0= coord(2,ipoin)
     z0= coord(3,ipoin)

     theta= atan2(z0,y0)

     ! radial imperfection
     alpha= 1.0d0 - impsize * cos(impmode * theta)

     ! modified nodal position
     coord(2,ipoin)= y0 * alpha
     coord(3,ipoin)= z0 * alpha

  end do


  return
end subroutine geoimpcyl