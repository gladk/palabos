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



if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)

    set(default_build_type "Release")
    message(WARNING "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})


#adding default options for all build types
LIST(APPEND WARNINGS_DEFAULT
    "-Wall"
    "-Wformat"
    "-Wformat-security"
    "-Werror=format-security"
    )

foreach(warning ${WARNINGS_DEFAULT})
    add_compile_options(${warning})
endforeach()

add_compile_options(
    "$<$<CONFIG:Debug>:-O0>"
    "$<$<CONFIG:Debug>:-fstack-protector>"
    "$<$<CONFIG:Debug>:--param=ssp-buffer-size=4>"

    "$<$<CONFIG:Release>:-O3>"

    "$<$<CONFIG:RelWithDebInfo>:-g>"
    "$<$<CONFIG:RelWithDebInfo>:-O3>"
)

#MESSAGE( STATUS "CMAKE_C_FLAGS: " ${CMAKE_C_FLAGS} )
#MESSAGE( STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS} )
#MESSAGE( STATUS "COMPILE_OPTIONS: " ${COMPILE_OPTIONS} )
#MESSAGE( STATUS "COMPILE_DEFINITIONS: " ${COMPILE_DEFINITIONS} )



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

SET(SITMO_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/externalLibraries")

ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")

INSTALL(FILES cMake/FindPalabos.cmake cMake/FindTinyXML.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX}/share/cmake-3.0/Modules
  )

