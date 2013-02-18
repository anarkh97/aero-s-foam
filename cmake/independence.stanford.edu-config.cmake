# use to tell cmake where to search for the include files and libraries of:
# acme, arpack, blacs, metis, mumps, scalapack, spooles, zoltan

SET(CMAKE_INCLUDE_PATH 
    /home/pavery/Intel/eigen
    /home/pavery/Intel/SPOOLES
    /home/pavery/Intel/MUMPS_4.10.0/include
    /home/pavery/Intel/trilinos-11.0.3-Source/packages/sacado/src
    /home/pavery/Intel/trilinos-11.0.3-Obj_cmake/include
    /home/pavery/Intel/trilinos-11.0.3-Source/packages/zoltan/src
    /home/pavery/Intel/stxxl-trunk/include
)
SET(CMAKE_LIBRARY_PATH 
    /home/pavery/Intel/ARPACK
    /home/pavery/Intel/SPOOLES
    /home/pavery/Intel/SPOOLES/MT/src
    /home/pavery/Intel/MUMPS_4.10.0/lib
    /home/pavery/Intel/MUMPS_4.10.0_seq/lib
    /home/pavery/Intel/trilinos-11.0.3-Obj_cmake/lib
    /home/pavery/Intel/stxxl-trunk/lib
)
## blas and lapack with scalapack and blacs (using intel math kernel library)
SET(LAPACK_LIBRARIES -Wl,--start-group /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64/libmkl_scalapack_lp64.a /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64/libmkl_blacs_openmpi_lp64.a  /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64/libmkl_sequential.a /opt/intel/composer_xe_2013.0.079/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread CACHE STRING "Path to a library.")
SET(LAPACK_FOUND true)
SET(BLACS_LIBRARIES "")
SET(BLACS_FOUND TRUE)
SET(SCALAPACK_LIBRARY "")
SET(SCALAPACK_FOUND TRUE)
#add_definitions(-DEIGEN_USE_MKL_ALL)
#include_directories("/opt/intel/composer_xe_2013.0.079/mkl/include")
