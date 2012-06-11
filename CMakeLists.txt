PROJECT (FEM)

#incomplete: parallel acme, zoltan, metis, spooles threads, acme version, apple, windows

cmake_minimum_required(VERSION 2.6  FATAL_ERROR)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# pjsa trying to build statically linked executable
#set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
#set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
#set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
#set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
#set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
#set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
#SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)
#SET(CMAKE_EXE_LINKER_FLAGS -static)

SITE_NAME(HOSTNAME)
FIND_PACKAGE(LOCAL NAMES ${HOSTNAME} PATHS . QUIET)
MARK_AS_ADVANCED(HOSTNAME LOCAL_DIR)

FIND_PACKAGE(OpenMP)
if(OpenMP_FLAG_DETECTED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OpenMP_FLAG_DETECTED)

SET(BLA_VENDOR Generic)
if(NOT LAPACK_FOUND)
 FIND_PACKAGE(LAPACK)
endif(NOT LAPACK_FOUND)

SET(INC_DIR .)
SET(BINDIR ${CMAKE_CURRENT_BINARY_DIR}/bin)
SET(LIBDIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
SET(EXTRALIB ${EXTRALIB} CACHE STRING "Extra link parameters")
SET(EXEC_DIR ${CMAKE_CURRENT_BINARY_DIR}/Regression.d)
SET(REG_SRC_DIR ${PROJECT_SOURCE_DIR}/Regression.d)
SET(BASELINE $EXEC_DIR/baseline)

if(AERO)
   MARK_AS_ADVANCED(CLEAR MPI_COMPILER)
   MARK_AS_ADVANCED(CLEAR MPI_LIBRARY)
   FIND_PACKAGE(MPI REQUIRED)
endif(AERO)

IF(UNIX)
  IF(${CMAKE_SYSTEM_NAME} STREQUAL Linux)
    add_definitions(-D__LINUX)

    SET(CMAKE_FIND_LIBRARY_PREFIXES "lib" "")

# zoltan
    IF(MPI_FOUND AND NOT WITH_ZOLTAN EQUAL 0)
      FIND_PATH(ZOLTAN_INCLUDE_PATH zoltan.h
                PATHS /usr/include /usr/local/include)
      FIND_PATH(ZOLTAN_SRC_PATH zz/zz_const.h
                PATHS /usr/include /usr/local/include)
      FIND_LIBRARY(ZOLTAN_LIBRARY NAMES zoltan
                   PATHS /usr/lib /usr/local/lib)
      IF(ZOLTAN_INCLUDE_PATH AND ZOLTAN_SRC_PATH AND ZOLTAN_LIBRARY)
        MESSAGE("-- A library with Zoltan API found.")
        SET(ZOLTAN_FOUND TRUE)
        SET(EXTRALIBRARIES ${ZOLTAN_LIBRARY} ${EXTRALIBRARIES})
        SET(INC_DIR ${INC_DIR} ${ZOLTAN_INCLUDE_PATH}
            ${ZOLTAN_SRC_PATH}/zz ${ZOLTAN_SRC_PATH}/lb ${ZOLTAN_SRC_PATH}/Utilities/shared
            ${ZOLTAN_SRC_PATH}/Utilities/Timer ${ZOLTAN_SRC_PATH}/par ${ZOLTAN_SRC_PATH}/params
            ${ZOLTAN_SRC_PATH}/order ${ZOLTAN_SRC_PATH}/tpls)
      ENDIF(ZOLTAN_INCLUDE_PATH AND ZOLTAN_SRC_PATH AND ZOLTAN_LIBRARY)
    ENDIF(MPI_FOUND AND NOT WITH_ZOLTAN EQUAL 0)

# acme
    IF(NOT MPI_FOUND)
      SET(ACME_FOUND TRUE)
      add_definitions(-DUSE_ACME -DCONTACT_NO_MPI -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0)
      SET(INC_DIR ${INC_DIR} "Acme.d/search" "Acme.d/enforcement")
    ELSE(NOT MPI_FOUND)
      IF(ZOLTAN_FOUND)
        SET(ACME_FOUND TRUE)
        add_definitions(-DUSE_ACME -DCONTACT_NO_EXODUS_OUTPUT -DCONTACT_DEBUG_PRINT_LEVEL=0)
        SET(INC_DIR ${INC_DIR} "Acme.d/search" "Acme.d/enforcement")
      ENDIF(ZOLTAN_FOUND)
    ENDIF(NOT MPI_FOUND)

# arpack
   IF(NOT WITH_ARPACK EQUAL 0)
    FIND_LIBRARY(ARPACK_arpack_LIBRARY NAMES arpack_LINUX arpack
                 PATHS /usr/lib /usr/local/lib)
    IF(NOT MPI_FOUND)
      MARK_AS_ADVANCED(FORCE ARPACK_parpack_LIBRARY)
      IF(ARPACK_arpack_LIBRARY)
        MESSAGE("-- A library with ARPACK API found.")
        SET(ARPACK_FOUND TRUE)
        add_definitions(-DUSE_ARPACK)
        SET(EXTRALIBRARIES ${ARPACK_arpack_LIBRARY} ${EXTRALIBRARIES})
      ENDIF(ARPACK_arpack_LIBRARY)
    ELSE(NOT MPI_FOUND)
      MARK_AS_ADVANCED(CLEAR ARPACK_parpack_LIBRARY)
      FIND_LIBRARY(ARPACK_parpack_LIBRARY NAMES parpack_MPI-LINUX parpack
                   PATHS /usr/lib /usr/local/lib)
      IF(ARPACK_arpack_LIBRARY AND ARPACK_parpack_LIBRARY)
        MESSAGE("-- A library with ARPACK API found.")
        SET(ARPACK_FOUND TRUE)
        add_definitions(-DUSE_ARPACK)
        SET(EXTRALIBRARIES ${ARPACK_parpack_LIBRARY} ${ARPACK_arpack_LIBRARY} ${EXTRALIBRARIES})
      ENDIF(ARPACK_arpack_LIBRARY AND ARPACK_parpack_LIBRARY)
    ENDIF(NOT MPI_FOUND)
   ELSE(NOT WITH_ARPACK EQUAL 0)
    SET(ARPACK_arpack_LIBRARY "ARPACK_arpack_LIBRARY-NOTFOUND")
    SET(ARPACK_parpack_LIBRARY "ARPACK_parpack_LIBRARY-NOTFOUND")
   ENDIF(NOT WITH_ARPACK EQUAL 0)

# spooles
    FIND_PATH(SPOOLES_INCLUDE_PATH FrontMtx.h
              PATHS /usr/include /usr/local/include /usr/include/spooles)
    #SET(CMAKE_FIND_LIBRARY_PREFIXES "lib" "")
    FIND_LIBRARY(SPOOLES_spooles_LIBRARY NAMES spooles
                 PATHS /usr/lib /usr/local/lib)
    FIND_LIBRARY(SPOOLES_spoolesMT_LIBRARY NAMES spoolesMT
                 PATHS /usr/lib /usr/local/lib)
    if(SPOOLES_INCLUDE_PATH AND SPOOLES_spooles_LIBRARY)
      MESSAGE("-- A library with SPOOLES API found.")
      SET(SPOOLES_FOUND TRUE)
    endif(SPOOLES_INCLUDE_PATH AND SPOOLES_spooles_LIBRARY)
    IF(SPOOLES_FOUND)
      add_definitions(-DUSE_SPOOLES)
      if(SPOOLES_spoolesMT_LIBRARY)
        SET(EXTRALIBRARIES ${SPOOLES_spoolesMT_LIBRARY} ${SPOOLES_spooles_LIBRARY} ${EXTRALIBRARIES})
      else(SPOOLES_spoolesMT_LIBRARY)
        SET(EXTRALIBRARIES ${SPOOLES_spooles_LIBRARY} ${EXTRALIBRARIES})
      endif(SPOOLES_spoolesMT_LIBRARY)
      SET(INC_DIR ${INC_DIR} ${SPOOLES_INCLUDE_PATH})
    ENDIF(SPOOLES_FOUND)

# eigen
   find_package(Eigen3 3.0.93)
   IF(EIGEN3_FOUND)
     include_directories(${EIGEN3_INCLUDE_DIR})
     add_definitions(-DUSE_EIGEN3)
     set(CMAKE_MODULE_PATH ${EIGEN3_INCLUDE_DIR}/cmake ${CMAKE_MODULE_PATH})
     find_package(Cholmod)
     if (CHOLMOD_INCLUDES AND CHOLMOD_LIBRARIES)
       include_directories(${CHOLMOD_INCLUDES})
       add_definitions(-DEIGEN_CHOLMOD_SUPPORT)
       SET(EXTRALIBRARIES ${CHOLMOD_LIBRARIES} ${EXTRALIBRARIES})
     endif(CHOLMOD_INCLUDES AND CHOLMOD_LIBRARIES)
     find_package(SuperLU)
     if (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)
       include_directories(${SUPERLU_INCLUDES})
       add_definitions(-DEIGEN_SUPERLU_SUPPORT)
       SET(EXTRALIBRARIES ${SUPERLU_LIBRARIES} ${EXTRALIBRARIES})
     endif(SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)
   ENDIF(EIGEN3_FOUND)

#sacado
    FIND_PATH(SACADO_INCLUDE_PATH NAMES Sacado.hpp )
    IF(SACADO_INCLUDE_PATH)
      MESSAGE("-- A library with SACADO API found. ")
      add_definitions(-DTRILINOS_NO_CONFIG_H -DUSE_SACADO)
      include_directories(${SACADO_INCLUDE_PATH} ${SACADO_INCLUDE_PATH}/mpl)
    ENDIF(SACADO_INCLUDE_PATH)

# metis
    FIND_LIBRARY(METIS_LIBRARY NAMES metis
                 PATHS /usr/lib /usr/local/lib)
    IF(METIS_LIBRARY)
      MESSAGE("-- A library with METIS API found.")
      add_definitions(-DUSE_METIS)
      SET(EXTRALIBRARIES ${METIS_LIBRARY} ${EXTRALIBRARIES})
    ENDIF(METIS_LIBRARY)

# mumps
    FIND_PATH(MUMPS_INCLUDE_PATH dmumps_c.h
              PATHS /usr/include usr/local/include)
    if (NOT MPI_FOUND)
      MARK_AS_ADVANCED(CLEAR MUMPS_dmumps_seq_LIBRARY MUMPS_zmumps_seq_LIBRARY
                       MUMPS_common_seq_LIBRARY MUMPS_pord_seq_LIBRARY
                       MUMPS_mpiseq_LIBRARY)
      MARK_AS_ADVANCED(FORCE MUMPS_dmumps_LIBRARY MUMPS_zmumps_LIBRARY
                       MUMPS_common_LIBRARY MUMPS_pord_LIBRARY)
      FIND_LIBRARY(MUMPS_dmumps_seq_LIBRARY NAMES dmumps_seq SEQ/libdmumps.a
                   PATHS /usr/lib /usr/local/lib)
      FIND_LIBRARY(MUMPS_zmumps_seq_LIBRARY NAMES zmumps_seq SEQ/libzmumps.a
                   PATHS /usr/lib /usr/local/lib)
      FIND_LIBRARY(MUMPS_common_seq_LIBRARY NAMES mumps_common_seq SEQ/libmumps_common.a
                   PATHS /usr/lib /usr/local/lib)
      FIND_LIBRARY(MUMPS_pord_seq_LIBRARY NAMES pord_seq SEQ/libpord.a
                   PATHS /usr/lib /usr/local/lib)
      FIND_LIBRARY(MUMPS_mpiseq_LIBRARY NAMES mpiseq mpiseq_seq SEQ/libmpiseq.a
                   PATHS /usr/lib /usr/local/lib)
      if(MUMPS_INCLUDE_PATH AND MUMPS_dmumps_seq_LIBRARY AND MUMPS_zmumps_seq_LIBRARY AND MUMPS_common_seq_LIBRARY AND MUMPS_pord_seq_LIBRARY AND MUMPS_mpiseq_LIBRARY)
        MESSAGE("-- A library with MUMPS API found.")
        SET(MUMPS_LIBRARIES ${MUMPS_dmumps_seq_LIBRARY} ${MUMPS_zmumps_seq_LIBRARY} ${MUMPS_common_seq_LIBRARY} ${MUMPS_pord_seq_LIBRARY} ${MUMPS_mpiseq_LIBRARY})
        SET(MUMPS_FOUND TRUE)
      endif(MUMPS_INCLUDE_PATH AND MUMPS_dmumps_seq_LIBRARY AND MUMPS_zmumps_seq_LIBRARY AND MUMPS_common_seq_LIBRARY AND MUMPS_pord_seq_LIBRARY AND MUMPS_mpiseq_LIBRARY)
      if(MUMPS_FOUND)
        add_definitions(-DUSE_MUMPS)
        SET(INC_DIR ${INC_DIR} ${MUMPS_INCLUDE_PATH})
        SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${MUMPS_LIBRARIES})
      endif(MUMPS_FOUND)
    else(NOT MPI_FOUND)
      MARK_AS_ADVANCED(FORCE MUMPS_dmumps_seq_LIBRARY MUMPS_zmumps_seq_LIBRARY
                       MUMPS_common_seq_LIBRARY MUMPS_pord_seq_LIBRARY
                       MUMPS_mpiseq_LIBRARY)
      MARK_AS_ADVANCED(CLEAR MUMPS_dmumps_LIBRARY MUMPS_zmumps_LIBRARY
                       MUMPS_common_LIBRARY MUMPS_pord_LIBRARY)
      FIND_LIBRARY(MUMPS_dmumps_LIBRARY NAMES dmumps
                   PATHS /usr/lib /usr/local/lib)
      FIND_LIBRARY(MUMPS_zmumps_LIBRARY NAMES zmumps        
                   PATHS /usr/lib /usr/local/lib)
      FIND_LIBRARY(MUMPS_common_LIBRARY NAMES mumps_common
                   PATHS /usr/lib /usr/local/lib)
      FIND_LIBRARY(MUMPS_pord_LIBRARY NAMES pord
                   PATHS /usr/lib /usr/local/lib)
      if(MUMPS_INCLUDE_PATH AND MUMPS_dmumps_LIBRARY AND MUMPS_zmumps_LIBRARY AND MUMPS_common_LIBRARY AND MUMPS_pord_LIBRARY)
        SET(MUMPS_LIBRARIES ${MUMPS_dmumps_LIBRARY} ${MUMPS_zmumps_LIBRARY} ${MUMPS_common_LIBRARY} ${MUMPS_pord_LIBRARY})
        MESSAGE("-- A library with MUMPS API found.")
        SET(MUMPS_FOUND TRUE)
      endif(MUMPS_INCLUDE_PATH AND MUMPS_dmumps_LIBRARY AND MUMPS_zmumps_LIBRARY AND MUMPS_common_LIBRARY AND MUMPS_pord_LIBRARY)
      # blacs is required for parallel mumps
      FIND_LIBRARY(BLACS_blacs_LIBRARY NAMES blacs blacs-openmpi blacs_MPI-LINUX-0 
                   PATHS /usr/lib /usr/local/lib)
      FIND_LIBRARY(BLACS_Cinit_LIBRARY NAMES blacsCinit blacsCinit-openmpi blacsCinit_MPI-LINUX-0 blacsC
                   PATHS /usr/lib /usr/local/lib)
      if(BLACS_blacs_LIBRARY AND BLACS_Cinit_LIBRARY)
        MESSAGE("-- A library with BLACS API found.")
        SET(BLACS_LIBRARIES ${BLACS_blacs_LIBRARY} ${BLACS_Cinit_LIBRARY} ${BLACS_blacs_LIBRARY})  
        SET(BLACS_FOUND TRUE)
      endif(BLACS_blacs_LIBRARY AND BLACS_Cinit_LIBRARY)
      # scalapack is required for parallel mumps
      FIND_LIBRARY(SCALAPACK_LIBRARY NAMES scalapack scalapack-openmpi
                   PATHS /usr/lib /usr/local/lib)
      if(SCALAPACK_LIBRARY)
        MESSAGE("-- A library with SCALAPACK API found.")
        SET(SCALAPACK_FOUND TRUE)
      endif(SCALAPACK_LIBRARY)

      if(BLACS_FOUND AND SCALAPACK_FOUND)
        if(MUMPS_FOUND)
          add_definitions(-DUSE_MUMPS)
          SET(INC_DIR ${INC_DIR} ${MUMPS_INCLUDE_PATH})
          SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${MUMPS_LIBRARIES})
        endif(MUMPS_FOUND)
        add_definitions(-DUSE_SCALAPACK)
        SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${SCALAPACK_LIBRARY} ${BLACS_LIBRARIES})
      endif(BLACS_FOUND AND SCALAPACK_FOUND)
    endif(NOT MPI_FOUND)

# metis
#    FIND_LIBRARY(METIS_LIBRARY NAMES metis
#                 PATHS /usr/lib /usr/local/lib)
#    IF(METIS_LIBRARY)
#      MESSAGE("-- A library with METIS API found.")
#      add_definitions(-DUSE_METIS)
#      SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${METIS_LIBRARY})
#    ENDIF(METIS_LIBRARY)

  FIND_PACKAGE(Boost 1.35)
  if(Boost_FOUND)
    add_definitions(-DUSE_BOOST)
    include_directories(${Boost_INCLUDE_DIRS})
  else(Boost_FOUND)
    unset(Boost_INCLUDE_DIR CACHE)
    unset(Boost_LIBRARY_DIRS CACHE)
  endif(Boost_FOUND)

  ENDIF(${CMAKE_SYSTEM_NAME} STREQUAL Linux)
ENDIF(UNIX)

IF(APPLE)
  add_definitions(-DNO_MALLOC_DOT_H)
ENDIF(APPLE)

if(MPI_FOUND)
   add_definitions(-DUSE_MPI -DDISTRIBUTED -DMPICH_IGNORE_CXX_SEEK)
#   add_definitions(-DUSE_MPI -DDISTRIBUTED -DAEROELASTIC -DMPICH_IGNORE_CXX_SEEK)
   SET(EXTRALIBDYN ${EXTRALIBRARIES})
   SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${MPI_LIBRARIES})
   SET(INC_DIR ${INC_DIR} ${MPI_INCLUDE_PATH})
   IF(NOT APPLE)
     ADD_LIBRARY(femlib main.C HgIdentity.C)
     target_link_libraries(femlib Parser
       Driver Elem Feti HelmAxi solver Corotational Mortar Math Driver
       Dec Threads linpa Sfem paral Problems Hetero Timers Comm Material Rom Util getopt
       ${EXTRALIBDYN} ${CMAKE_DL_LIBS} ${LAPACK_LIBRARIES})
     SET_TARGET_PROPERTIES(femlib
                           PROPERTIES 
                           ARCHIVE_OUTPUT_DIRECTORY ${LIBDIR}
                           RELEASE_OUTPUT_NAME aeros
                           DEBUG_OUTPUT_NAME aeros.debug
                           PREFIX "")
   ENDIF(NOT APPLE)
endif(MPI_FOUND)

FIND_PROGRAM(MERCURIAL_EXECUTABLE hg2xx)

IF(MERCURIAL_EXECUTABLE)
  #MESSAGE("Rule 1")
  MESSAGE("Mercurial is: " ${MERCURIAL_EXECUTABLE})
  add_custom_command(OUTPUT HgIdentity.C 
                     COMMAND echo "const char *THE_VERSION = \\\"" `${MERCURIAL_EXECUTABLE} identify -in` "\\\"" "\;" > HgIdentity.C
                     DEPENDS main.C
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
ELSE (MERCURIAL_EXECUTABLE)
  #MESSAGE("Rule 2")
  add_custom_command(OUTPUT HgIdentity.C 
                     COMMAND echo "const char *THE_VERSION = \\\"\\\"\;" > HgIdentity.C
                     DEPENDS main.C
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
ENDIF (MERCURIAL_EXECUTABLE)

add_definitions(-D_TEMPLATE_FIX_ -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS  -DMAP_MIN_MEMORY)

INCLUDE_DIRECTORIES( ${INC_DIR} )

add_subdirectory(Element.d)
add_subdirectory(Feti.d)
add_subdirectory(Driver.d)
add_subdirectory(Comm.d)
add_subdirectory(Corotational.d)
add_subdirectory(Dec.d)
add_subdirectory(HelmAxi.d)
add_subdirectory(Solvers.d)
add_subdirectory(Utils.d)
add_subdirectory(Parser.d)
add_subdirectory(Timers.d)
add_subdirectory(Threads.d)
add_subdirectory(Mortar.d)
add_subdirectory(Math.d)
add_subdirectory(Linpack.d)
add_subdirectory(Sfem.d)
add_subdirectory(Paral.d)
add_subdirectory(Problems.d)
add_subdirectory(GNU-getopt.d)
add_subdirectory(Hetero.d)
add_subdirectory(Material.d)
add_subdirectory(Rom.d)
add_subdirectory(Regression.d)

SET(ALLLIBS Parser
     Driver Elem Feti HelmAxi solver Corotational Mortar Math Driver
     Dec Threads linpa Sfem paral Problems Hetero Timers Comm Material Rom Util getopt)

IF(ACME_FOUND AND (ZOLTAN_FOUND OR NOT MPI_FOUND))
  add_subdirectory(Acme.d)
  SET(ALLLIBS ${ALLLIBS} Acme)
ENDIF(ACME_FOUND AND (ZOLTAN_FOUND OR NOT MPI_FOUND))

IF(MPI_FOUND)
  add_subdirectory(Pita.d)
  SET(ALLLIBS Pita ${ALLLIBS} )
ENDIF(MPI_FOUND)

SET_TARGET_PROPERTIES(${ALLLIBS}
                PROPERTIES COMPILE_FLAGS -fPIC)

add_executable(femexecutable  main.C HgIdentity.C)

SET_TARGET_PROPERTIES(femexecutable 
                      PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
                      RELEASE_OUTPUT_NAME aeros
                      DEBUG_OUTPUT_NAME aeros.debug
#                      LINK_SEARCH_END_STATIC 1
                      )

SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${EXTRALIB})
if(MPI_FOUND)
  SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${EXTRALIB_MPI})
endif(MPI_FOUND)

target_link_libraries(femexecutable ${ALLLIBS}
       ${EXTRALIBRARIES} ${CMAKE_DL_LIBS} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_executable(robcodec Rom.d/RobCodec.C)
set_target_properties(robcodec
                      PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
                      RELEASE_OUTPUT_NAME rob
                      DEBUG_OUTPUT_NAME rob.debug)

target_link_libraries(robcodec Rom Util)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING 
      "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(AllTargets femexecutable robcodec)
IF(MPI_FOUND AND NOT APPLE)
  SET(AllTargets ${AllTargets} femlib)
ENDIF(MPI_FOUND AND NOT APPLE)

INSTALL(TARGETS ${AllTargets}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)


MESSAGE("")
MESSAGE("=================================================")
MESSAGE("           Summary of build options")
MESSAGE("-------------------------------------------------")
IF(MPI_FOUND)
  MESSAGE("Distributed FETI:         YES")
  MESSAGE("Aeroelastic:              YES")
  IF(MUMPS_FOUND AND SCALAPACK_FOUND AND BLACS_FOUND)
    MESSAGE("Mumps:                    YES")
  ELSE(MUMPS_FOUND AND SCALAPACK_FOUND AND BLACS_FOUND)
    MESSAGE("Mumps:                    NO ")
  ENDIF(MUMPS_FOUND AND SCALAPACK_FOUND AND BLACS_FOUND)
  IF(SCALAPACK_FOUND AND BLACS_FOUND)
    MESSAGE("Scalapack:                YES")
  ELSE(SCALAPACK_FOUND AND BLACS_FOUND)
    MESSAGE("Scalapack:                NO ")
  ENDIF(SCALAPACK_FOUND AND BLACS_FOUND)
ELSE(MPI_FOUND)
  MESSAGE("Distributed FETI:         NO ")
  MESSAGE("Aeroelastic:              NO ")
  IF(MUMPS_FOUND)
    MESSAGE("Mumps:                    YES")
  ELSE(MUMPS_FOUND)
    MESSAGE("Mumps:                    NO ") 
  ENDIF(MUMPS_FOUND)
  MESSAGE("Scalapack:                NO ")
ENDIF(MPI_FOUND)
IF(ARPACK_FOUND)
  MESSAGE("Arpack:                   YES")
ELSE(ARPACK_FOUND)
  MESSAGE("Arpack:                   NO ")
ENDIF(ARPACK_FOUND)
IF(SPOOLES_FOUND)
  MESSAGE("Spooles:                  YES")
ELSE(SPOOLES_FOUND)
  MESSAGE("Spooles:                  NO ")
ENDIF(SPOOLES_FOUND)
IF(ACME_FOUND)
  MESSAGE("Acme:                     YES")
ELSE(ACME_FOUND)
  MESSAGE("Acme:                     NO ")
ENDIF(ACME_FOUND)
IF(METIS_LIBRARY)
  MESSAGE("Metis:                    YES")
ELSE(METIS_LIBRARY)
  MESSAGE("Metis:                    NO ")
ENDIF(METIS_LIBRARY)
if(OpenMP_FLAG_DETECTED)
  MESSAGE("OpenMP:                   YES")
else(OpenMP_FLAG_DETECTED)
  MESSAGE("OpenMP:                   NO ")
endif(OpenMP_FLAG_DETECTED)
MESSAGE("Build type:               " ${CMAKE_BUILD_TYPE})
MESSAGE("Extra libraries:          " ${EXTRALIB})
MESSAGE("=================================================")
MESSAGE("")

FIND_PATH(BASELINE REQUIRED)

MARK_AS_ADVANCED(FORCE BISON_EXECUTABLE FLEX_EXECUTABLE MERCURIAL_EXECUTABLE CMAKE_INSTALL_PREFIX)
MARK_AS_ADVANCED(CLEAR CMAKE_CXX_COMPILER CMAKE_Fortran_COMPILER)
ENABLE_TESTING()
SET (COPY_FILES ${REG_SRC_DIR}/copy_files.pl ${REG_SRC_DIR} ${EXEC_DIR} )
SET (TEST_SHORT ${EXEC_DIR}/dataComp.py -l short)
SET (BUILD_SHORT ${EXEC_DIR}/buildInputs.py ALL)

SET (TEST_DURATION_SHORT $ENV{TEST_DURATION_SHORT})  

ADD_CUSTOM_TARGET(test-short 
COMMAND ${COPY_FILES}
COMMAND ${BUILD_SHORT}
COMMAND ${TEST_SHORT}
)

if(TEST_DURATION_SHORT)
  ADD_TEST(copyfiles ${REG_SRC_DIR}/copy_files.pl ${REG_SRC_DIR} ${EXEC_DIR} )
  ADD_TEST(buildInputs ${EXEC_DIR}/buildInputs.py ALL)
  ADD_TEST(testShort ${EXEC_DIR}/dataComp.py -l -n short)
else(TEST_DURATION_SHORT)
  ADD_TEST(copyfiles ${REG_SRC_DIR}/copy_files.pl ${REG_SRC_DIR} ${EXEC_DIR} )
  ADD_TEST(buildInputs ${EXEC_DIR}/buildInputs.py ALL)
  ADD_TEST(testStatics ${EXEC_DIR}/dataComp.py -l -n statics)
  ADD_TEST(testNLstatics ${EXEC_DIR}/dataComp.py -l nlstatics)
  ADD_TEST(testEigen ${EXEC_DIR}/dataComp.py -l eigen)
  ADD_TEST(testDynamics ${EXEC_DIR}/dataComp.py -l dynamics)
  ADD_TEST(testNLDynamics ${EXEC_DIR}/dataComp.py -l nldynamics)
  ADD_TEST(testIMPE ${EXEC_DIR}/dataComp.py -l impe)
  ADD_TEST(testFreqSweep ${EXEC_DIR}/dataComp.py -l freqsweep)
  ADD_TEST(testTempStatics ${EXEC_DIR}/dataComp.py -l tempstatics)
  ADD_TEST(testTempNLStatics ${EXEC_DIR}/dataComp.py -l tempnlstatics)
  ADD_TEST(testTempDynamics ${EXEC_DIR}/dataComp.py -l tempdynamics)
  ADD_TEST(testTempNLDynamics ${EXEC_DIR}/dataComp.py -l tempnldynamics)
  ADD_TEST(testRtest1 ${EXEC_DIR}/dataComp.py -l dsvm1)
  ADD_TEST(testRtest2 ${EXEC_DIR}/dataComp.py -r dsvm2)
  ADD_TEST(testRtest11 ${EXEC_DIR}/dataComp.py -r dsvm11)
  ADD_TEST(testRtest13 ${EXEC_DIR}/dataComp.py -l dsvm13)
  #ADD_TEST(testRtest15 ${EXEC_DIR}/dataComp.py -r dsvm15)
  ADD_TEST(testRtest19 ${EXEC_DIR}/dataComp.py -l dsvm19)
  ADD_TEST(testRtest20 ${EXEC_DIR}/dataComp.py -l dsvm20)
  ADD_TEST(testRtest21 ${EXEC_DIR}/dataComp.py -l dsvm21)
  ADD_TEST(testRtest22 ${EXEC_DIR}/dataComp.py -l dsvm22)
  ADD_TEST(testRtest23 ${EXEC_DIR}/dataComp.py -l dsvm23)
  ADD_TEST(testRtest24 ${EXEC_DIR}/dataComp.py -r dsvm24)
  ADD_TEST(testRtest25 ${EXEC_DIR}/dataComp.py -l dsvm25)
  ADD_TEST(testRtest27a ${EXEC_DIR}/dataComp.py -l dsvm27a)
  ADD_TEST(testRtest27b ${EXEC_DIR}/dataComp.py -l dsvm27b)
  ADD_TEST(testRtest29 ${EXEC_DIR}/dataComp.py -r dsvm29)
  ADD_TEST(testRtest30 ${EXEC_DIR}/dataComp.py -l dsvm30)
  ADD_TEST(testRtest31 ${EXEC_DIR}/dataComp.py -l dsvm31)
  ADD_TEST(testRtest32 ${EXEC_DIR}/dataComp.py -l dsvm32)
  ADD_TEST(testRtest34 ${EXEC_DIR}/dataComp.py -l dsvm34)
  ADD_TEST(testRtest35a ${EXEC_DIR}/dataComp.py -l dsvm35a)
  ADD_TEST(testRtest35b ${EXEC_DIR}/dataComp.py -l dsvm35b)
  ADD_TEST(testRtest37 ${EXEC_DIR}/dataComp.py -r dsvm37)
#  ADD_TEST(testRtest38 ${EXEC_DIR}/dataComp.py -r dsvm38)
  ADD_TEST(testRtest39 ${EXEC_DIR}/dataComp.py -l dsvm39)
  # use the -s flag to send results via email
  ADD_TEST(testRtest40 ${EXEC_DIR}/dataComp.py -l -s dsvm40)
  ADD_TEST(sendResults ${EXEC_DIR}/sendResults.py )
endif(TEST_DURATION_SHORT)
  
