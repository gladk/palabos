FILE(GLOB_RECURSE SRC_FILES_SRC "*.cpp")

INCLUDE_DIRECTORIES(${SITMO_INCLUDE_DIR})

ADD_LIBRARY(plb SHARED ${SRC_FILES_SRC})
SET_TARGET_PROPERTIES(plb PROPERTIES 
  VERSION ${PALABOS_MAJOR_VERSION}.${PALABOS_MINOR_VERSION}.${PALABOS_PATCH_VERSION}
  SOVERSION ${PALABOS_MAJOR_VERSION})
TARGET_LINK_LIBRARIES(plb ${TINYXML_LIBRARIES})
INSTALL(TARGETS plb DESTINATION "${CMAKE_INSTALL_LIBDIR}/")

#=======================================

FILE(GLOB SRC_FILES_HEAD "*.h")
FILE(GLOB SRC_FILES_HHEAD "*.hh")
INSTALL(FILES ${SRC_FILES_HEAD} ${SRC_FILES_HHEAD} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/palabos/")

SET(SUBDIRS  algorithm atomicBlock basicDynamics boundaryCondition complexDynamics 
  coProcessors core dataProcessors finiteDifference gridRefinement io latticeBoltzmann libraryInterfaces
  multiBlock multiGrid multiPhysics offLattice parallelism particles)

FOREACH(subdir ${SUBDIRS})
  FILE(GLOB SRC_FILES_HEAD_${subdir} "${subdir}/*.h")
  FILE(GLOB SRC_FILES_HHEAD_${subdir} "${subdir}/*.hh")
  INSTALL(FILES ${SRC_FILES_HEAD_${subdir}} ${SRC_FILES_HHEAD_${subdir}} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/palabos/${subdir}/")
ENDFOREACH()
#=======================================

set(SITMO_FOLDER "${SITMO_INCLUDE_DIR}/sitmo/")
FILE(GLOB SITMO_HEADER "${SITMO_FOLDER}/*.hpp")
INSTALL(FILES ${SITMO_HEADER} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/palabos/sitmo/")
