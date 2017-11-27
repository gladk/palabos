# ENABLE_MPI: enable MPI-modules (ON by default)
# ENABLE_POSIX: use POSIX (ON by default)
# ENABLE_SMP_PARALLEL: use SMP_PARALLEL (ON by default)
# VERSION: version number (1.4.1 by default)

PROJECT(PALABOS CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET(CMAKE_VERBOSE_MAKEFILE 1)
INCLUDE(GNUInstallDirs)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cMake")

#=======================================

SET(PALABOS_MAJOR_VERSION 2)
SET(PALABOS_MINOR_VERSION 0)
SET(PALABOS_PATCH_VERSION 0)

#=======================================

FIND_PACKAGE(TinyXML REQUIRED)
FIND_PACKAGE(Eigen3 REQUIRED)

#=======================================

IF (NOT CMAKE_CXX_FLAGS)
  #If flags are not set, add default flags
  SET(CMAKE_BUILD_TYPE Release)
  SET(CMAKE_CXX_FLAGS  " -Wall -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security")
ENDIF (NOT CMAKE_CXX_FLAGS)

#=======================================

SET(DEFAULT ON CACHE INTERNAL "Default value for enabled by default options")
OPTION(ENABLE_MPI "Enable MPI" ${DEFAULT})

IF(ENABLE_MPI)
  INCLUDE(FindMPI)
  IF(MPI_C_FOUND AND MPI_CXX_FOUND)
    SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
    SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
    ADD_DEFINITIONS("-DPLB_MPI_PARALLEL")
  ELSE(MPI_C_FOUND AND MPI_CXX_FOUND)
    #Check old FindMPI version
    IF(MPI_COMPILER)
      SET(CMAKE_C_COMPILER ${MPI_COMPILER})
      SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
      ADD_DEFINITIONS("-DPLB_MPI_PARALLEL")
    ELSE(MPI_COMPILER)
      MESSAGE(FATAL_ERROR "MPI-COMPILER NOT found!")
    ENDIF(MPI_COMPILER)
  ENDIF(MPI_C_FOUND AND MPI_CXX_FOUND)
ENDIF(ENABLE_MPI)

#=======================================

OPTION(ENABLE_POSIX "Enable POSIX" ${DEFAULT})

IF(ENABLE_POSIX)
  ADD_DEFINITIONS("-DPLB_USE_POSIX")
ENDIF(ENABLE_POSIX)

#=======================================

OPTION(ENABLE_SMP_PARALLEL "Enable SMP_PARALLEL" ${DEFAULT})

IF(ENABLE_SMP_PARALLEL)
  ADD_DEFINITIONS("-DPLB_SMP_PARALLEL")
ENDIF(ENABLE_SMP_PARALLEL)

#=======================================

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(${TINYXML_INCLUDE_DIR})

ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")

INSTALL(FILES cMake/FindPalabos.cmake cMake/FindTinyXML.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share/cmake-3.0/Modules
  )
