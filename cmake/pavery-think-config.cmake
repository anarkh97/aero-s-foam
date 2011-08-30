# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH /home/pavery/Codes/eigen;/home/pavery/Codes/Acme/acme-2.9/search;/home/pavery/Codes/Acme/acme-2.9/enforcement;/home/pavery/Codes/Zoltan/include)
SET(CMAKE_LIBRARY_PATH /home/pavery/Codes/Acme/acme-2.9/lib;/home/pavery/Codes/Zoltan/Obj_generic;/home/pavery/Codes/GotoBLAS2/lib-single-threaded;/home/pavery/Codes/ARPACK)
#SET(CMAKE_INCLUDE_PATH /home/pavery/Codes/Acme/acme-2.5e/search;/home/pavery/Codes/Acme/acme-2.5e/enforcement;/home/pavery/Codes/Zoltan/include)
#SET(CMAKE_LIBRARY_PATH /home/pavery/Codes/Acme/acme-2.5e/lib;/home/pavery/Codes/Zoltan/Obj_generic;/home/pavery/Codes/GotoBLAS2/lib-single-threaded;/home/pavery/Codes/ARPACK)
#SET(CMAKE_INCLUDE_PATH /home/pavery/Codes/Acme/myacme-1.3a/search;/home/pavery/Codes/Acme/myacme-1.3a/enforcement;/home/pavery/Codes/Acme/myacme-1.3a/drivers)
#SET(CMAKE_LIBRARY_PATH /home/pavery/Codes/Acme/myacme-1.3a/lib)
SET(EXTRALIB_MPI /usr/lib/libmpif77.so
                 CACHE STRING "Extra MPI link parameters")
SET(ACME_VERSION "2_9")
