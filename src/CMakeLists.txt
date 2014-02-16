FILE(GLOB_RECURSE SRC_FILES_SRC "*.cpp")

ADD_LIBRARY(libplb SHARED ${SRC_FILES_SRC})
SET_TARGET_PROPERTIES(libplb PROPERTIES 
  OUTPUT_NAME "libplb"
  VERSION ${PALABOS_MAJOR_VERSION}.${PALABOS_MINOR_VERSION}.${PALABOS_PATCH_VERSION}
  SOVERSION ${PALABOS_MAJOR_VERSION})
TARGET_LINK_LIBRARIES(libplb ${TINYXML_LIBRARIES})
INSTALL(TARGETS libplb DESTINATION "${CMAKE_INSTALL_LIBDIR}/")

#=======================================

FILE(GLOB SRC_FILES_HEAD "*.h")
FILE(GLOB SRC_FILES_HHEAD "*.hh")
INSTALL(FILES ${SRC_FILES_HEAD} ${SRC_FILES_HHEAD} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/palabos/")

SET(SUBDIRS  algorithm atomicBlock basicDynamics boundaryCondition complexDynamics 
  coProcessors core dataProcessors finiteDifference io latticeBoltzmann libraryInterfaces
  multiBlock multiGrid multiPhysics offLattice parallelism particles)

FOREACH(subdir ${SUBDIRS})
  FILE(GLOB SRC_FILES_HEAD_${subdir} "${subdir}/*.h")
  FILE(GLOB SRC_FILES_HHEAD_${subdir} "${subdir}/*.hh")
  INSTALL(FILES ${SRC_FILES_HEAD_${subdir}} ${SRC_FILES_HHEAD_${subdir}} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/palabos/${subdir}/")
ENDFOREACH()

#=======================================

