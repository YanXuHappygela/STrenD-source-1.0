cmake_minimum_required(VERSION 2.8.10)

project( STrenDAnalysis )

if("${STrenDAnalysis_BINARY_DIR}" MATCHES "^${STrenDAnalysis_SOURCE_DIR}$")
  message(FATAL_ERROR "STrenDAnalysis must be built with separate source and binary directories.")
endif("${STrenDAnalysis_BINARY_DIR}" MATCHES "^${STrenDAnalysis_SOURCE_DIR}$")

include_directories( ${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR} )

set( EXE_DIR ${STrenDAnalysis_BINARY_DIR}/exe CACHE PATH
  "Single Directory for all executables")
set( LIB_DIR ${STrenDAnalysis_BINARY_DIR}/libs CACHE PATH
  "Single Directory for all static Libraries")

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${EXE_DIR} )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB_DIR} )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR} )


#Find ITK
find_package(ITK REQUIRED)
include( ${ITK_USE_FILE} )

#Find VTK
find_package(VTK REQUIRED)
include( ${VTK_USE_FILE} )

#Find VXL
find_package(VXL REQUIRED)
set(VXL_PROVIDE_STANDARD_optionS 1)
include(${VXL_CMAKE_DIR}/UseVXL.cmake)
# UseITK.cmake utilizes include_directories(BEFORE..., so we do the same here
# because we want the VXL vxl_config.h to be found first.
include_directories(BEFORE ${VXL_VCL_INCLUDE_DIR} ${VXL_CORE_INCLUDE_DIR}
			${VXL_RPL_INCLUDE_DIR} ${VXL_MUL_INCLUDE_DIR})
# Find Qt
find_package(Qt4)
if(QT4_FOUND)
    include( ${QT_USE_FILE} )
else(QT4_FOUND)
    message(FATAL_ERROR "Cannot build without Qt4")
endif(QT4_FOUND)

# Find Boost
find_package(Boost)
include_directories ( ${Boost_INCLUDE_DIRS})

add_subdirectory(STrenDEditor)
add_subdirectory(ClusClus)
add_subdirectory(STrenD)

set( common_libraries
  ftkCommon
  ftkGUI
  ftkGraphs
  ftkImage
  ftkFeatures
  tinyxml
  )
  
foreach( _lib ${common_libraries} )
  include_directories( ${PROJECT_SOURCE_DIR}/${_lib} )
endforeach()
foreach( _lib ${common_libraries} )
  add_subdirectory( ${_lib} )
endforeach()
