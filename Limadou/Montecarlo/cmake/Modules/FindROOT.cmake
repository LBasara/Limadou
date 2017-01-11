# - Finds ROOT installation, copied from geant4/cmake/Modules
# This module sets up ROOT information 
# It defines:
# ROOT_FOUND          If the ROOT is found
# ROOT_INCLUDE_DIR    PATH to the include directory
# ROOT_LIBRARIES      Most common libraries
# ROOT_LIBRARY_DIR    PATH to the library directory 


find_program(ROOT_CONFIG_EXECUTABLE root-config
  PATHS $ENV{ROOTSYS}/bin)

if(NOT ROOT_CONFIG_EXECUTABLE)
  set(ROOT_FOUND FALSE)
else()    
  set(ROOT_FOUND TRUE)

  execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --prefix 
    OUTPUT_VARIABLE ROOTSYS 
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --version 
    OUTPUT_VARIABLE ROOT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --incdir
    OUTPUT_VARIABLE ROOT_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs --noldflags --noauxlibs --glibs
    OUTPUT_VARIABLE ROOT_LIBRARIES
    OUTPUT_STRIP_TRAILING_WHITESPACE)	

  execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --bindir 
    OUTPUT_VARIABLE ROOT_BINARY_DIR 
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(
    COMMAND ${ROOT_CONFIG_EXECUTABLE} --libdir
    OUTPUT_VARIABLE ROOT_LIBRARY_DIR_TMP
    OUTPUT_STRIP_TRAILING_WHITESPACE)	
    
  if(EXISTS "${ROOT_LIBRARY_DIR_TMP}")
    set(ROOT_LIBRARY_DIR ${ROOT_LIBRARY_DIR_TMP})
  else(EXISTS "${ROOT_LIBRARY_DIR_TMP}")
    MESSAGE("Warning: ROOT_CONFIG_EXECUTABLE reported ${ROOT_LIBRARY_DIR_TMP} as library path,")
    MESSAGE("Warning: but ${ROOT_LIBRARY_DIR_TMP} does NOT exist, ROOT must NOT be installed correctly.")
    set(ROOT_LIBRARY_DIR ${ROOTSYS}/lib)
  ENDIF(EXISTS "${ROOT_LIBRARY_DIR_TMP}")

  # Make variables ch angeble to the advanced user
  mark_as_advanced(ROOT_CONFIG_EXECUTABLE)

  FIND_PROGRAM(ROOT_CINT_EXECUTABLE
    NAMES rootcint
    PATHS ${ROOT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  if(NOT ROOT_FIND_QUIETLY)
    message(STATUS "Found ROOT ${ROOT_VERSION} in ${ROOTSYS}!!")
  endif()
endif()

  ###########################################
  #
  #       Macros for building ROOT dictionary
  #
  ###########################################

MACRO (ROOT_GENERATE_DICTIONARY INFILES LINKDEF_FILE OUTFILE INCLUDE_DIRS_IN)
 
  set(INCLUDE_DIRS)
  foreach (_current_FILE ${INCLUDE_DIRS_IN})
    set(INCLUDE_DIRS ${INCLUDE_DIRS} -I${_current_FILE})   
  endforeach (_current_FILE ${INCLUDE_DIRS_IN})
  
  STRING(REGEX REPLACE "^(.*)\\.(.*)$" "\\1.h" bla "${OUTFILE}")

  SET (OUTFILES ${OUTFILE} ${bla})


  if (CMAKE_SYSTEM_NAME MATCHES Linux)
    ADD_CUSTOM_COMMAND(OUTPUT ${OUTFILES}
       COMMAND LD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS} ${ROOT_CINT_EXECUTABLE}
       ARGS -f ${OUTFILE} -c -DHAVE_CONFIG_H ${INCLUDE_DIRS} ${INFILES} ${LINKDEF_FILE} DEPENDS ${INFILES} ${LINKDEF_FILE})
  else (CMAKE_SYSTEM_NAME MATCHES Linux)
    if (CMAKE_SYSTEM_NAME MATCHES Darwin)
      ADD_CUSTOM_COMMAND(OUTPUT ${OUTFILES}
       COMMAND DYLD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS} ${ROOT_CINT_EXECUTABLE}
       ARGS -f ${OUTFILE} -c -DHAVE_CONFIG_H ${INCLUDE_DIRS} ${INFILES} ${LINKDEF_FILE} DEPENDS ${INFILES} ${LINKDEF_FILE})
    endif (CMAKE_SYSTEM_NAME MATCHES Darwin)
  endif (CMAKE_SYSTEM_NAME MATCHES Linux)

ENDMACRO (ROOT_GENERATE_DICTIONARY)

