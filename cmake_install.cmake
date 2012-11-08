# Install script for directory: /lustre/home/tac688/source_code/FEM3

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aeros" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aeros")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aeros"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/lustre/home/tac688/source_code/FEM3/bin/aeros")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aeros" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aeros")
    FILE(RPATH_REMOVE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aeros")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/aeros")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rob" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rob")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rob"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/lustre/home/tac688/source_code/FEM3/bin/rob")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rob" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rob")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/rob")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/relerr" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/relerr")
    FILE(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/relerr"
         RPATH "")
  ENDIF()
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/lustre/home/tac688/source_code/FEM3/bin/relerr")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/relerr" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/relerr")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/relerr")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/lustre/home/tac688/source_code/FEM3/lib/aeros.a")
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Element.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Feti.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Driver.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Comm.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Corotational.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Dec.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/HelmAxi.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Solvers.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Utils.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Parser.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Timers.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Threads.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Mortar.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Math.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Linpack.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Sfem.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Paral.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Problems.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/GNU-getopt.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Hetero.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Material.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Rom.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Regression.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Acme.d/cmake_install.cmake")
  INCLUDE("/lustre/home/tac688/source_code/FEM3/Pita.d/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

IF(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
ELSE(CMAKE_INSTALL_COMPONENT)
  SET(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
ENDIF(CMAKE_INSTALL_COMPONENT)

FILE(WRITE "/lustre/home/tac688/source_code/FEM3/${CMAKE_INSTALL_MANIFEST}" "")
FOREACH(file ${CMAKE_INSTALL_MANIFEST_FILES})
  FILE(APPEND "/lustre/home/tac688/source_code/FEM3/${CMAKE_INSTALL_MANIFEST}" "${file}\n")
ENDFOREACH(file)
