# CMakeLists.txt for event package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(event)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package( ROOT REQUIRED COMPONENTS RIO Net Tree Hist )
#find_package(ROOT_INCLUDE_DIRS)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#ROOT_GENERATE_DICTIONARY(G__Event Event.h LINKDEF EventLinkDef.h)
#to take into account also the objects inside the src folder
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

#---Create a shared library with geneated dictionary
#add_library(Event SHARED ${SOURCES})
#target_link_libraries(Event ${ROOT_LIBRARIES})



#---Create  a main program using the library
add_executable(analysis_ester main.cc ${sources})
add_executable(CalibrateTracker CalibrateTracker.cc ${sources})
add_executable(analysis_follega analysis_follega.cc ${sources})

target_link_libraries(analysis_ester ${ROOT_LIBRARIES})
target_link_libraries(CalibrateTracker ${ROOT_LIBRARIES})
target_link_libraries(analysis_follega ${ROOT_LIBRARIES})
#target_link_libraries(Main Event)
