# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH /home/avery/Codes/eigen
                       /share/apps/lib/arpack/src
                       /share/apps/lib/spooles-2.2/src
#                       /share/apps/lib/MUMPS_4.9.2/src/include
                       /home/avery/Codes/mumps/MUMPS_4.9.2/include
                       /share/apps/lib/acme-2.9/src/search /share/apps/lib/acme-2.9/src/enforcement
                       /share/apps/lib/ZOLTAN/src/include)
SET(CMAKE_LIBRARY_PATH /usr/lib64
                       /opt/intel/Compiler/11.1/064/mkl/lib/em64t
                       /share/apps/lib/arpack/gcc_openmpi143
                       /share/apps/lib/spooles-2.2/gcc_openmpi143
                       /share/apps/lib/scalapack-1.8.0/gcc_openmpi143
#                       /share/apps/lib/MUMPS_4.9.2/gcc_openmpi143
                       /home/avery/Codes/mumps/MUMPS_4.9.2/lib
                       /share/apps/lib/metis-4.0.1/gcc
                       /share/apps/lib/acme-2.9/gcc_openmpi143
                       /share/apps/lib/ZOLTAN/gcc_openmpi143)
# generic blas and lapack
SET(BLAS_LIBRARIES /home/avery/Codes/lapack-3.3.0/BLAS/SRC/libblas.a)
SET(LAPACK_LIBRARIES /home/avery/Codes/lapack-3.3.0/SRC/liblapack.a ${BLAS_LIBRARIES})
SET(LAPACK_FOUND true)
# blas and lapack and scalapack and blacs (using intel math kernel library)
#SET(LAPACK_LIBRARIES -Wl,--start-group /opt/intel/composerxe-2011.2.137/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composerxe-2011.2.137/mkl/lib/intel64/libmkl_sequential.a /opt/intel/composerxe-2011.2.137/mkl/lib/intel64/libmkl_core.a -Wl,--end-group CACHE STRING "Path to a library.")
#SET(LAPACK_FOUND true)
#SET(BLACS_LIBRARIES /opt/intel/composerxe-2011.2.137/mkl/lib/intel64/libmkl_blacs_openmpi_lp64.a CACHE STRING "Path to a library.")
#SET(BLACS_FOUND TRUE)
#SET(SCALAPACK_LIBRARY /opt/intel/composerxe-2011.2.137/mkl/lib/intel64/libmkl_scalapack_lp64.a CACHE STRING "Path to a library.")
#SET(SCALAPACK_FOUND TRUE)

SET(EXTRALIB /usr/lib64/libg2c.so.0
              CACHE STRING "Extra link parameters")
SET(EXTRALIB_MPI /usr/mpi/gcc/openmpi-1.4.3-1/lib64/libmpi_f77.so
                 CACHE STRING "Extra MPI link parameters")

#SET(BLA_VENDOR Intel10_64lp)
