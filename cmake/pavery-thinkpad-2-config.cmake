# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH /home/avery/Codes/eigen;/home/avery/Codes/acme-2.5e/search;/home/avery/Codes/acme-2.5e/enforcement;/home/avery/Codes/Zoltan/include;/home/avery/Codes/trilinos/trilinos-10.10.2-Source/packages/sacado/src)
SET(CMAKE_LIBRARY_PATH /home/avery/Codes/acme-2.5e/lib;/home/avery/Codes/Zoltan/Obj_generic)
SET(EXTRALIB_MPI /usr/lib/libmpif77.so
                 CACHE STRING "Extra MPI link parameters")
SET(ACME_VERSION "2_5")
SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
