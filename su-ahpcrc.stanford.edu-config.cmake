# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH /lustre/home/avery/code
                       /lustre/home/avery/code/spooles
                       /lustre/home/avery/code/MUMPS_4.9.1/include
                       /lustre/home/avery/code/acme-2.5e/search /lustre/home/avery/code/acme-2.5e/enforcement
                       /lustre/home/avery/code/Zoltan/include)
SET(CMAKE_LIBRARY_PATH /lustre/home/avery/code/arpack
                       /lustre/home/avery/code/spooles
                       /lustre/home/avery/code/spooles/MT/src
                       /lustre/home/avery/code/BLACS/LIB
                       /lustre/home/avery/code/scalapack-1.8.0
                       /lustre/home/avery/code/MUMPS_4.9.1/lib
                       /lustre/home/avery/code/metis-4.0.1
                       /lustre/home/avery/code/acme-2.5e/lib
                       /lustre/home/avery/code/Zoltan/Obj_generic)
# add anything missed by cmake auto-detection
SET(EXTRALIB /opt/intel/fce/10.1.015/lib/libifcore.so 
             /opt/intel/fce/10.1.015/lib/libsvml.so
             CACHE STRING "Extra link parameters")
SET(EXTRALIB_MPI /usr/mpi/gcc/openmpi-1.2.6/lib64/libmpi_f77.so
                 CACHE STRING "Extra MPI link parameters")

#SET(BLA_VENDOR Intel10_64lp)
