######################################################################
## CMakeLists.txt --- gsUnstructuredSplines
## This file is part of the G+Smo library.
##
######################################################################

## gsUnstructuredSplines module
project(gsUnstructuredSplines)

# Apply G+Smo config
include(gsConfig)

## Collect files
aux_header_directory     (${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_H  )
aux_cpp_directory        (${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_CPP)
aux_tmpl_header_directory(${CMAKE_CURRENT_SOURCE_DIR}/src ${PROJECT_NAME}_HPP)

if( (NOT GISMO_BUILD_LIB) )
  aux_instance_directory (${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_NAME}_INS)
  if(${PROJECT_NAME}_INS)
    LIST( REMOVE_ITEM ${PROJECT_NAME}_CPP ${${PROJECT_NAME}_INS})
  endif()
endif()

# Add object library
add_library(${PROJECT_NAME} OBJECT
  ${${PROJECT_NAME}_H}
  ${${PROJECT_NAME}_HPP}
  ${${PROJECT_NAME}_CPP} )

set_target_properties(${PROJECT_NAME} PROPERTIES
  COMPILE_DEFINITIONS gismo_EXPORTS
  POSITION_INDEPENDENT_CODE ON
  LINKER_LANGUAGE CXX
  FOLDER "G+Smo modules" )

set(gismo_MODULES ${gismo_MODULES} $<TARGET_OBJECTS:${PROJECT_NAME}>
  CACHE INTERNAL "G+Smo modules" )

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/src"
  DESTINATION include/gismo
  FILES_MATCHING PATTERN "*.h" )

# add filedata folder
add_definitions(-DELAST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/filedata/")

# add example files
if(GISMO_BUILD_EXAMPLES)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tutorials)
else()
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/examples EXCLUDE_FROM_ALL)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tutorials EXCLUDE_FROM_ALL)
endif(GISMO_BUILD_EXAMPLES)

# add unittests
if(GISMO_BUILD_UNITTESTS)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unittests)
else()
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/unittests EXCLUDE_FROM_ALL)
endif(GISMO_BUILD_UNITTESTS)
