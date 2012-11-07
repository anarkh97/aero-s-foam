# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH
    /home/avery/Codes/eigen
#    /usr/include/trilinos
    /home/avery/Codes/trilinos/trilinos-11.0.3-Source/packages/sacado/src
    /home/avery/Codes/trilinos/trilinos-11.0.3-Obj_cmake/include
    /home/avery/Codes/trilinos/trilinos-11.0.3-Source/packages/zoltan/src
#    /home/avery/Codes/stxxl-trunk/include
)
SET(CMAKE_LIBRARY_PATH
    /home/avery/Codes/ARPACK
    /home/avery/Codes/trilinos/trilinos-11.0.3-Obj_cmake/lib
#    /home/avery/Codes/stxxl-trunk/lib
)
#SET(EXTRALIB_MPI /usr/lib/libmpif77.so
#                 CACHE STRING "Extra MPI link parameters")
#SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
SET(BLAS_blas_LIBRARY "/home/avery/Codes/eigen-build/blas/libeigen_blas.so" CACHE FILEPATH "Path to a library.")
add_definitions(-D_AEROS_ASYCHRONOUS_IO)
