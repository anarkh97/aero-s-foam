# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
#SET(CMAKE_INCLUDE_PATH /home/avery/Codes/eigen;/home/avery/Codes/Zoltan/include)
#SET(CMAKE_LIBRARY_PATH /home/avery/Codes/Zoltan/Obj_generic)
SET(CMAKE_INCLUDE_PATH /home/avery/Codes/eigen;
    /home/avery/Codes/trilinos/trilinos-10.8.3-Obj_cmake/include
    /home/avery/Codes/trilinos/trilinos-10.8.3-Source/packages/zoltan/src)
SET(CMAKE_LIBRARY_PATH /home/avery/Codes/trilinos/trilinos-10.8.3-Obj_cmake/lib)
SET(EXTRALIB_MPI /usr/lib/libmpif77.so
                 CACHE STRING "Extra MPI link parameters")
SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
