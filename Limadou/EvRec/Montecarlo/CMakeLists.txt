#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(HEPDSW)


#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# - Prepend local CMake Modules to the search path

set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake/Modules
    ${CMAKE_MODULE_PATH})



#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

set(hepdsw_include_dirs 
	    ${CMAKE_CURRENT_SOURCE_DIR}/source/Simulation/Core/include 
	    ${CMAKE_CURRENT_SOURCE_DIR}/source/Simulation/MCTruth/include 
	    ${CMAKE_CURRENT_SOURCE_DIR}/source/Simulation/Tracker/include 
	    ${CMAKE_CURRENT_SOURCE_DIR}/source/Simulation/Calorimeter/include 
            ${CMAKE_CURRENT_SOURCE_DIR}/source/Geometry/Common/include 
            ${CMAKE_CURRENT_SOURCE_DIR}/source/Geometry/Tracker/include 
            ${CMAKE_CURRENT_SOURCE_DIR}/source/Geometry/Calorimeter/include 
            ${CMAKE_CURRENT_SOURCE_DIR}/source/Geometry/Satellite/include 
            ${CMAKE_CURRENT_SOURCE_DIR}/source/Geometry/HEPDBox/include 
            ${CMAKE_CURRENT_SOURCE_DIR}/source/Dataformats/include 
	    ${CMAKE_CURRENT_SOURCE_DIR}/source/Producer/include 		   
	    ${Geant4_INCLUDE_DIR})   

include_directories(${hepdsw_include_dirs})
#----------------------------------------------------------------------------
# Source subdirectory
#
set(headers "")
set(sources "")
add_subdirectory(source)
get_directory_property(headers DIRECTORY source DEFINITION headers)
get_directory_property(sources DIRECTORY source DEFINITION sources)

include(FindROOT)
if(ROOT_FOUND)
	MESSAGE( STATUS "HEPDSW will be compiled with ROOT.")
	include_directories(${ROOT_INCLUDE_DIR})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_ANALYSIS_ROOT")
	set(EXTRA_LIBRARIES "${EXTRA_LIBRARIES} ${ROOT_LIBRARIES}")
	set (HEPDSW_BUILTWITH_ROOT "yes")
else()
	MESSAGE(FATAL_ERROR "ROOT installation not found!Installation interrupted!")
endif()

#---------------------------------------------------------
#set directories where library and executable will be set

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

set(DataFormats_LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/source/Dataformats/include/LinkDef.h )
set(DataFormats_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/DataFormatDict.cxx)
set(DataFormatsHeaders "")
get_directory_property(DataFormatsHeaders DIRECTORY source/Dataformats DEFINITION rootdfheaders)
ROOT_GENERATE_DICTIONARY("${DataFormatsHeaders}" "${DataFormats_LINKDEF}" "${DataFormats_DICTIONARY}" "${hepdsw_include_dirs}")

# add the dictionary to the list of source files
set(sources ${sources} ${DataFormats_DICTIONARY}) 

if(Geant4_gdml_FOUND)
    MESSAGE( STATUS "HEPDSW will be compiled with GDML!")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHEPDSW_USE_G4GDML") 
    find_package(XercesC)
    include_directories(${XERCESC_INCLUDE_DIR})
endif(Geant4_gdml_FOUND)


#----------------------------------------------------------------------------
# HEPDSW library
option(STATIC_BUILD "Build HEPDSW static library" ON)
STRING(TOUPPER ${STATIC_BUILD} staticbuild)
STRING(COMPARE EQUAL ${staticbuild} "ON" test_res)

if (test_res)
    add_library(hepdsw_lib STATIC ${sources} ${headers})

else()
    add_library(hepdsw_lib SHARED ${sources} ${headers})
    
endif()
set_target_properties(hepdsw_lib PROPERTIES OUTPUT_NAME "hepdsw")
target_link_libraries(hepdsw_lib  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES}) 


#----------------------------------------------------------------------------
# HEPDSW executable
add_executable(hepdsw_exe ${PROJECT_SOURCE_DIR}/source/Executable/HEPDSW.cc ${headers}) 
add_dependencies(hepdsw_exe hepdsw_lib)

set_target_properties(hepdsw_exe PROPERTIES OUTPUT_NAME "hepdsw")
target_link_libraries(hepdsw_exe hepdsw_lib)


#----------------------------------------------------------------------------
#set HEPDSW_INSTALL_PREFIX
#
set(HEPDSW_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE PATH "Prefix for the HEPDSW installation")

#--------------------------------------
#Building of hepdsw-env.sh 

set (HEPDSW_BASH_BINPATH_SETUP "export PATH=${HEPDSW_INSTALL_PREFIX}/bin:$PATH")
if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(libpathname DYLD_LIBRARY_PATH)
else()
    set(libpathname LD_LIBRARY_PATH)
endif()
set (HEPDSW_BASH_LIBPATH_SETUP "export ${libpathname}=${HEPDSW_INSTALL_PREFIX}/lib:$${libpathname}")   
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/hepdsw-env.sh.in
    ${PROJECT_BINARY_DIR}/config/hepdsw-env.sh
    @ONLY)

#--------------------------------------
#Building of hepdsw-env.csh 

set (HEPDSW_CSH_BINPATH_SETUP "setenv  PATH ${HEPDSW_INSTALL_PREFIX}/bin:$PATH")     
if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(libpathname DYLD_LIBRARY_PATH)
else()
    set(libpathname LD_LIBRARY_PATH)
endif()
set (HEPDSW_CSH_LIBPATH_SETUP "setenv ${libpathname} ${HEPDSW_INSTALL_PREFIX}/lib:$${libpathname}")   
configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/hepdsw-env.csh.in
    ${PROJECT_BINARY_DIR}/config/hepdsw-env.csh
    @ONLY) 


#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS hepdsw_exe DESTINATION ${HEPDSW_INSTALL_PREFIX}/bin)
install(TARGETS hepdsw_lib DESTINATION ${HEPDSW_INSTALL_PREFIX}/lib)
install(FILES ${PROJECT_BINARY_DIR}/config/hepdsw-env.sh
    DESTINATION ${HEPDSW_INSTALL_PREFIX}/bin
    PERMISSIONS
      OWNER_READ OWNER_WRITE OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE
      WORLD_READ WORLD_EXECUTE
    COMPONENT Development
    )
install(FILES ${PROJECT_BINARY_DIR}/config/hepdsw-env.csh
    DESTINATION ${HEPDSW_INSTALL_PREFIX}/bin
    PERMISSIONS
      OWNER_READ OWNER_WRITE OWNER_EXECUTE
      GROUP_READ GROUP_EXECUTE
      WORLD_READ WORLD_EXECUTE
    COMPONENT Development
    )
