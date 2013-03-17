# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH
    /home/avery/Codes/eigen
    /home/avery/Codes/MUMPS_4.10.0/include
    /home/avery/Codes/trilinos/trilinos-10.10.2-Source/packages/sacado/src
    /home/avery/Codes/trilinos/trilinos-10.10.2-Obj_cmake/include
    /home/avery/Codes/trilinos/trilinos-10.10.2-Source/packages/zoltan/src)
SET(CMAKE_LIBRARY_PATH
    /home/avery/Codes/ARPACK
    /home/avery/Codes/MUMPS_4.10.0/lib
    /home/avery/Codes/trilinos/trilinos-10.10.2-Obj_cmake/lib)
SET(EXTRALIB_MPI /usr/lib/libmpif77.so
                 CACHE STRING "Extra MPI link parameters")
SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
SET(BLAS_blas_LIBRARY "/home/avery/Codes/eigen-build/blas/libeigen_blas.so" CACHE FILEPATH "Path to a library.")
add_definitions(-DEIGEN_SPARSELU_SUPPORT -DEIGEN_SPARSEQR_SUPPORT)
