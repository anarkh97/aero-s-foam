# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH /home/avery/Codes
                       /share/apps/lib/arpack/src
                       /share/apps/lib/spooles-2.2/src
                       /share/apps/lib/MUMPS_4.9.2/src/include
                       /share/apps/lib/acme-2.9/src/search /share/apps/lib/acme-2.9/src/enforcement
                       /share/apps/lib/ZOLTAN/src/include)
SET(CMAKE_LIBRARY_PATH /usr/lib64
                       /opt/intel/Compiler/11.1/064/mkl/lib/em64t
                       /share/apps/lib/arpack/gcc_openmpi143
                       /share/apps/lib/spooles-2.2/gcc_openmpi143
                       /share/apps/lib/scalapack-1.8.0/gcc_openmpi143
                       /share/apps/lib/MUMPS_4.9.2/gcc_openmpi143
                       /share/apps/lib/metis-4.0.1/gcc
                       /share/apps/lib/acme-2.9/gcc_openmpi143
                       /share/apps/lib/ZOLTAN/gcc_openmpi143)
# blas and lapack (using intel math kernel library)
SET(LAPACK_LIBRARIES -Wl,--start-group /opt/intel/Compiler/11.1/064/mkl/lib/em64t/libmkl_intel_lp64.so /opt/intel/Compiler/11.1/064/mkl/lib/em64t/libmkl_sequential.so /opt/intel/Compiler/11.1/064/mkl/lib/em64t/libmkl_lapack.so /opt/intel/Compiler/11.1/064/mkl/lib/em64t/libmkl_core.so -Wl,--end-group CACHE STRING "Path to a library.")
SET(LAPACK_FOUND true)
# add anything missed by cmake auto-detection
SET(EXTRALIB /usr/lib64/libg2c.so.0
              CACHE STRING "Extra link parameters")
SET(EXTRALIB_MPI /usr/mpi/gcc/openmpi-1.4.3-1/lib64/libmpi_f77.so
                 CACHE STRING "Extra MPI link parameters")

#SET(BLA_VENDOR Intel10_64lp)
