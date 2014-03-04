# Try to find the Palabos librairies
#  PALABOS_FOUND - system has GMP lib
#  PALABOS_INCLUDE_DIR - the GMP include directory
#  PALABOS_LIBRARIES - Libraries needed to use GMP

# Copyright (c) 2014, Anton Gladky <gladk@debian.org> 
#
# Redistribution and use is allowed according to the terms of the GPL license.

IF (PALABOS_INCLUDE_DIR AND PALABOS_LIBRARIES)
  SET(PALABOS_FIND_QUIETLY TRUE)
ENDIF (PALABOS_INCLUDE_DIR AND PALABOS_LIBRARIES)

FIND_PATH(PALABS_INCLUDE_DIR NAMES palabos3d.h palabos2d.h )
FIND_PATH(PALABOS_LIBRARIES NAMES palabos )
MESSAGE(STATUS "Palabos libs: " ${PALABOS_LIBRARIES} " )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Palabos DEFAULT_MSG PALABOS_INCLUDE_DIR PALABOS_LIBRARIES)
