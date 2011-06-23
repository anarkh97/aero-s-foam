# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH /home/pavery/lib/eigen
                       /home/pavery/lib/SPOOLES
                       /home/pavery/lib/MUMPS/include
                       /home/pavery/lib/ACME/search /home/pavery/lib/ACME/enforcement)
SET(CMAKE_LIBRARY_PATH /home/pavery/lib/ARPACK
                       /home/pavery/lib/SPOOLES
                       /home/pavery/lib/SPOOLES/MT/src
                       /home/pavery/lib/MUMPS/lib
                       /home/pavery/lib/ACME/lib)
SET(ACME_VERSION "2_5")
# add anything missed by cmake auto-detection
SET(EXTRALIB /usr/lib/gcc/x86_64-redhat-linux5E/4.1.2/libgfortran.a -lpthread
             CACHE STRING "Extra link parameters")

