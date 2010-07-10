! =====================================
! predefined parameters
! =====================
module preset
   implicit none

   ! I/O device number
   ! --------------------------
   integer, parameter :: ninp1= 8  ! main input
   integer, parameter :: nout1= 10 ! main output file device number
   integer, parameter :: nout2= 11 ! crack path file
   integer, parameter :: nout3= 7  ! simulation log file
   integer, parameter :: nout4= 6  ! monitor
   integer, parameter :: nout5= 9  ! error and temporary message log file
   integer, parameter :: nout6= 15 ! mesh file
   integer, parameter :: nout7= 14 ! load-dflection curve

   integer, parameter :: ntmp1= 13 ! temporary file
   integer, parameter :: ntmp2= 12 ! temporary

   ! pre-defined constant
   ! -----------------------
   real(8), parameter :: pi= 3.141592653589790d0
   real(8), parameter :: d2r= 0.017453292519940d0  ! degree to radian (3.14/180.0)
   real(8), parameter :: r2d= 57.295779513082320d0 ! radian to degree (180.0/3.14)
   real(8), parameter :: s2ms= 1.0e6 ! sec to micro sec

   ! global variable
   ! ---------------
   real(8), dimension(3) :: toler ! tolerence
   common toler
end module preset

