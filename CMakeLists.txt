PROJECT (FEM)

cmake_minimum_required(VERSION 2.6  FATAL_ERROR)
enable_language(Fortran)

IF(NOT EXTRALIB)
SET(EXTRALIB "")
ENDIF(NOT EXTRALIB)

SET(EXTRALIB ${EXTRALIB} CACHE STRING "Extra link parameters")

FIND_PACKAGE(BLAS)
SET(EXTRALIBRARIES ${BLAS_LIBRARIES})
FIND_PACKAGE(LAPACK)
FIND_PACKAGE(MPI)

MESSAGE("Found Blas library at " ${BLAS_LIBRARIES})

SET(INC_DIR .)
SET(BINDIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

if(AERO)
   set(USE_MPI 1)
   add_definitions(-DAEROELASTIC)
endif(AERO)

IF(NOT SPOOLESDIR)
  SET(SPOOLESDIR "")
ENDIF(NOT SPOOLESDIR)

SET(SPOOLESDIR ${SPOOLESDIR} CACHE STRING "Location of the SPOOLES library and includes")

IF(UNIX)
    MESSAGE(${CMAKE_SYSTEM_NAME})
  IF(${CMAKE_SYSTEM_NAME} STREQUAL Linux)
    FIND_LIBRARY(ARPACKLIB NAMES arpack arpack_Linux64
    PATHS  . /usr/lib /usr/local/lib)
    MESSAGE("Found ARPACK in " ${ARPACKLIB})
    IF(ARPACKLIB AND LAPACK_FOUND)
      add_definitions(-DUSE_ARPACK)
      SET(EXTRALIBRARIES ${ARPACKLIB} ${EXTRALIBRARIES})
    ENDIF(ARPACKLIB AND LAPACK_FOUND)
    SET(CMAKE_FIND_LIBRARY_PREFIXES "lib" "")
    FIND_LIBRARY(SPOOLESLIB NAMES spooles
                 PATHS  . /usr/lib /usr/local/lib /home/hg/lib/SPOOLES ${SPOOLESDIR})
    FIND_LIBRARY(SPOOLESMTLIB NAMES spoolesMT
                 PATHS  . /usr/lib /usr/local/lib /home/hg/lib/SPOOLES/MT/src
                       ${SPOOLESDIR}/MT/src)
    MESSAGE("Found SPOOLES in " ${SPOOLESLIB})
    MESSAGE("Found SPOOLESMT in " ${SPOOLESMTLIB})
    FIND_PATH(SPOOLESINC FrontMtx.h 
              PATHS . /usr/lib /usr/local/lib /home/hg/lib/SPOOLES ${SPOOLESDIR})
    MESSAGE("SPOOLES PATH in " ${SPOOLESINC})
    FIND_LIBRARY(PTHREAD pthread)
    MESSAGE("Found pthread in " ${PTHREAD})
    IF(SPOOLESLIB AND SPOOLESMTLIB AND SPOOLESINC AND PTHREAD)
      add_definitions(-DUSE_SPOOLES)
      SET(EXTRALIBRARIES ${SPOOLESLIB} ${SPOOLESMTLIB} ${EXTRALIBRARIES} ${PTHREAD})
      SET(INC_DIR ${INC_DIR} ${SPOOLESINC})
    ENDIF(SPOOLESLIB AND SPOOLESMTLIB AND SPOOLESINC AND PTHREAD)
  ENDIF(${CMAKE_SYSTEM_NAME} STREQUAL Linux)
ENDIF(UNIX)

IF(APPLE)
  add_definitions(-DNO_MALLOC_DOT_H)
ENDIF(APPLE)

if(USE_MPI)
   MESSAGE("Will build using MPI")
   add_definitions(-DUSE_MPI -DDISTRIBUTED -DMPICH_IGNORE_CXX_SEEK)
   SET(EXTRALIBDYN ${EXTRALIBRARIES})
   SET(EXTRALIBRARIES ${MPI_LIBRARIES} ${EXTRALIBRARIES})
   SET(INC_DIR ${INC_DIR} ${MPI_INCLUDE_PATH})
   IF(NOT APPLE)
     ADD_LIBRARY(femlib SHARED main.C HgIdentity.C)
     target_link_libraries(femlib Parser
       Driver Elem Feti HelmAxi solver Corrotational Mortar Math Util Driver
       Dec Threads linpa Sfem paral Problems Hetero Timers Comm getopt
       ${EXTRALIBDYN} ${CMAKE_DL_LIBS} ${LAPACK_LIBRARIES})
     SET_TARGET_PROPERTIES(femlib
                        PROPERTIES 
		        LIBRARY_OUTPUT_DIRECTORY ${BINDIR}
			RELEASE_OUTPUT_NAME aeros
			DEBUG_OUTPUT_NAME aeros.debug
			PREFIX ""
		     )
   ENDIF(NOT APPLE)
endif(USE_MPI)

FIND_PROGRAM(MERCURIAL_EXECUTABLE hg2xx)

MESSAGE("Mercurial is: " ${MERCURIAL_EXECUTABLE})

IF(MERCURIAL_EXECUTABLE)
  MESSAGE("Rule 1")
  add_custom_command(OUTPUT HgIdentity.C 
                   COMMAND echo "const char *THE_VERSION = \\\"" `${MERCURIAL_EXECUTABLE} identify -in` "\\\"" "\;" > HgIdentity.C
                 DEPENDS main.C
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ELSE (MERCURIAL_EXECUTABLE)
  MESSAGE("Rule 2")
  add_custom_command(OUTPUT HgIdentity.C 
                   COMMAND echo "const char *THE_VERSION = \\\"\\\"\;" > HgIdentity.C
                 DEPENDS main.C
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
ENDIF (MERCURIAL_EXECUTABLE)

add_definitions(-D_TEMPLATE_FIX_ -D__LINUX -DF_NEEDS_UNDSC -DMORTAR_LOCALNUMBERING -DSOWER_SURFS  -DMAP_MIN_MEMORY)

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

SET(ALLLIBS Parser
     Driver Elem Feti HelmAxi solver Corrotational Mortar Math Util Driver
     Dec Threads linpa Sfem paral Problems Hetero Timers Comm getopt)

IF(USE_MPI)
  add_subdirectory(Pita.d)
  add_subdirectory(OOPita.d)
  SET(ALLLIBS pita ${ALLLIBS} pitadriver )
ENDIF(USE_MPI)

SET_TARGET_PROPERTIES(${ALLLIBS}
                PROPERTIES COMPILE_FLAGS -fPIC)

add_executable(femexecutable  main.C HgIdentity.C)

SET_TARGET_PROPERTIES(femexecutable 
                      PROPERTIES
		         RUNTIME_OUTPUT_DIRECTORY ${BINDIR}
			 RELEASE_OUTPUT_NAME aeros
			 DEBUG_OUTPUT_NAME aeros.debug
		      )

MESSAGE(Lapack is ${LAPACK_LIBRARIES})

SET(EXTRALIBRARIES ${EXTRALIBRARIES} ${EXTRALIB})
target_link_libraries(femexecutable ${ALLLIBS}
       ${EXTRALIBRARIES} ${CMAKE_DL_LIBS} ${LAPACK_LIBRARIES})

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release)
ENDIF(NOT CMAKE_BUILD_TYPE)


SET(AllTargets femexecutable)
IF(USE_MPI AND NOT APPLE)
  SET(AllTargets ${AllTargets} femlib)
ENDIF(USE_MPI AND NOT APPLE)

INSTALL(TARGETS ${AllTargets}
        RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib)

MESSAGE("")
MESSAGE("=================================================")
MESSAGE("           Summary of build options")
MESSAGE("-------------------------------------------------")
IF(USE_MPI)
 MESSAGE("Distributed FETI:         YES")
 MESSAGE("Aeroelastic:              YES")
ELSE(USE_MPI)
 MESSAGE("Distributed FETI:         NO *******")
 MESSAGE("Aeroelastic:              NO *******")
ENDIF(USE_MPI)
IF(ARPACKLIB AND LAPACK_FOUND)
 MESSAGE("Arpack:                   YES")
ELSE(ARPACKLIB AND LAPACK_FOUND)
 MESSAGE("Arpack:                   NO *******")
ENDIF(ARPACKLIB AND LAPACK_FOUND)
IF(SPOOLESLIB AND SPOOLESMTLIB AND SPOOLESINC AND PTHREAD)
 MESSAGE("Spooles:                  YES")
ELSE(SPOOLESLIB AND SPOOLESMTLIB)
 MESSAGE("Spooles:                  NO *******")
ENDIF(SPOOLESLIB AND SPOOLESMTLIB AND SPOOLESINC AND PTHREAD)
MESSAGE("Build type:               " ${CMAKE_BUILD_TYPE})
MESSAGE("Extra libraries:          " ${EXTRALIB})
MESSAGE("=================================================")
MESSAGE("")



