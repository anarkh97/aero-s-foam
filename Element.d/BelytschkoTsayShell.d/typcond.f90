module typcond
   type cond_bc
      integer :: ele
      integer, dimension(2) :: nod
      integer, dimension(6) :: dof
      real(8), dimension(6) :: val
      integer, dimension(6) :: tmf
  end type cond_bc
end module typcond

