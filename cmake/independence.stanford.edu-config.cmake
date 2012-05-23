# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
#SET(CMAKE_INCLUDE_PATH /home/avery/Codes/eigen
#                       /share/apps/lib/arpack/src
#                       /share/apps/lib/spooles-2.2/src
#                       /home/avery/Codes/mumps/MUMPS_4.9.2/include
#                       /share/apps/lib/ZOLTAN/src
#                       /share/apps/lib/ZOLTAN/src/include)
#SET(CMAKE_LIBRARY_PATH /usr/lib64
#                       /opt/intel/Compiler/11.1/064/mkl/lib/em64t
#                       /share/apps/lib/arpack/gcc_openmpi143
#                       /share/apps/lib/spooles-2.2/gcc_openmpi143
#                       /share/apps/lib/scalapack-1.8.0/gcc_openmpi143
##                       /share/apps/lib/MUMPS_4.9.2/gcc_openmpi143
#                       /home/avery/Codes/mumps/MUMPS_4.9.2/lib
#                       /share/apps/lib/metis-4.0.1/gcc
#                       /share/apps/lib/ZOLTAN/gcc_openmpi143)
## generic blas and lapack
#SET(BLAS_LIBRARIES /home/avery/Codes/lapack-3.3.0/BLAS/SRC/libblas.a)
#SET(LAPACK_LIBRARIES /home/avery/Codes/lapack-3.3.0/SRC/liblapack.a ${BLAS_LIBRARIES})
#SET(LAPACK_FOUND true)


SET(CMAKE_INCLUDE_PATH /home/avery/Intel/eigen
                       /home/avery/Intel/sacado/src
                       /home/avery/Intel/SPOOLES
                       /home/avery/Intel/MUMPS_4.10.0/include
                       /home/avery/Intel/Zoltan/Zoltan_v3.5/include
                       /home/avery/Intel/Zoltan/Zoltan_v3.5/src)
SET(CMAKE_LIBRARY_PATH /home/avery/Intel/ARPACK
                       /home/avery/Intel/SPOOLES
                       /home/avery/Intel/SPOOLES/MT/src
                       /home/avery/Intel/MUMPS_4.10.0/lib
                       /home/avery/Intel/MUMPS_4.10.0_seq/lib
                       /home/avery/Intel/Zoltan/Zoltan_v3.5/lib)
## blas and lapack without scalapack and blacs (using intel math kernel library)
#SET(LAPACK_LIBRARIES -Wl,--start-group /opt/intel/compilerpro-12.0.2.137/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/compilerpro-12.0.2.137/mkl/lib/intel64/libmkl_sequential.a /opt/intel/compilerpro-12.0.2.137/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread CACHE STRING "Path to a library.")
#SET(LAPACK_FOUND true)

## blas and lapack with scalapack and blacs (using intel math kernel library)
SET(LAPACK_LIBRARIES -Wl,--start-group /opt/intel/compilerpro-12.0.2.137/mkl/lib/intel64/libmkl_scalapack_lp64.a /opt/intel/compilerpro-12.0.2.137/mkl/lib/intel64/libmkl_blacs_openmpi_lp64.a  /opt/intel/compilerpro-12.0.2.137/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/compilerpro-12.0.2.137/mkl/lib/intel64/libmkl_sequential.a /opt/intel/compilerpro-12.0.2.137/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread CACHE STRING "Path to a library.")
SET(LAPACK_FOUND true)
SET(BLACS_LIBRARIES "")
SET(BLACS_FOUND TRUE)
SET(SCALAPACK_LIBRARY "")
SET(SCALAPACK_FOUND TRUE)

#
#SET(EXTRALIB /usr/lib64/libg2c.so.0
#              CACHE STRING "Extra link parameters")
#SET(EXTRALIB_MPI /usr/mpi/gcc/openmpi-1.4.3-1/lib64/libmpi_f77.so
#                 CACHE STRING "Extra MPI link parameters")
#
#SET(BLA_VENDOR Intel10_64lp)
