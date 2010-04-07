# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
SET(CMAKE_INCLUDE_PATH /lustre/home/avery/code/spooles
                       /lustre/home/avery/code/MUMPS_4.9.1/include)
SET(CMAKE_LIBRARY_PATH /lustre/home/avery/code/arpack
                       /lustre/home/avery/code/spooles
                       /lustre/home/avery/code/spooles/MT/src
                       /lustre/home/avery/code/BLACS/LIB
                       /lustre/home/avery/code/scalapack-1.8.0
                       /lustre/home/avery/code/MUMPS_4.9.1/lib
                       /lustre/home/avery/code/metis-4.0.1)
# add anything missed by cmake auto-detection
SET(EXTRALIB /opt/intel/fce/10.1.015/lib/libifcore.so 
             /opt/intel/fce/10.1.015/lib/libsvml.so
             CACHE STRING "Extra link parameters")
if(MPI_FOUND)
 SET(EXTRALIB ${EXTRALIB} /usr/mpi/gcc/openmpi-1.2.6/lib64/libmpi_f77.so)
endif(MPI_FOUND)

#SET(BLA_VENDOR Intel10_64lp)
