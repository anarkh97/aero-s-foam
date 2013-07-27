# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan
# note: for openmpi 1.2.6 + gcc 4.1 use libs built in ~pavery/code.copy
SET(CMAKE_INCLUDE_PATH 
    /lustre/home/pavery/code
    /lustre/home/pavery/Intel12/spooles
    /lustre/home/pavery/Intel12/MUMPS_4.10.0/include
    /lustre/home/pavery/Intel12/trilinos/trilinos-11.0.3-Obj_cmake/include
    /lustre/home/pavery/Intel12/trilinos/trilinos-11.0.3-Source/packages/zoltan/src
    /lustre/home/pavery/Intel12/trilinos/trilinos-11.0.3-Source/packages/sacado/src
    /lustre/home/mpotts/boost/include
)
SET(CMAKE_LIBRARY_PATH 
    /lustre/home/pavery/Intel12/ARPACK
    /lustre/home/pavery/Intel12/spooles
    /lustre/home/pavery/Intel12/spooles/MT/src
    /lustre/home/pavery/Intel12/MUMPS_4.10.0/lib
    /lustre/home/pavery/Intel12/trilinos/trilinos-11.0.3-Obj_cmake/lib
    /lustre/home/mpotts/boost/lib
)

