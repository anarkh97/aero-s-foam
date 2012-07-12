# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
# note: for openmpi 1.2.6 + gcc 4.1 use libs built in ~avery/code.copy
SET(CMAKE_INCLUDE_PATH 
    /lustre/home/avery/code
    /lustre/home/avery/sacado/src
    /lustre/home/avery/code/spooles
    /lustre/home/avery/code/MUMPS_4.10.0/include
    /lustre/home/avery/code/trilinos/trilinos-10.8.4-Obj_cmake/include
    /lustre/home/avery/code/trilinos/trilinos-10.8.4-Source/packages/zoltan/src
)
SET(CMAKE_LIBRARY_PATH 
    /lustre/home/avery/code/arpack
    /lustre/home/avery/code/spooles
    /lustre/home/avery/code/spooles/MT/src
    /lustre/home/avery/code/BLACS/LIB
    /lustre/home/avery/code/scalapack-1.8.0
    /lustre/home/avery/code/MUMPS_4.10.0/lib
    /lustre/home/avery/code/trilinos/trilinos-10.8.4-Obj_cmake/lib
)

