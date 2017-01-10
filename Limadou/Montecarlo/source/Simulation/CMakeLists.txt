#----------------------------------------------------------------------------
# process subdirectories
#
set(headers "")
set(sources "")

foreach(subdir Core MCTruth Calorimeter Tracker)
  message(STATUS "process subdir :" ${subdir})
  add_subdirectory(${subdir})
  get_directory_property(sub_headers DIRECTORY ${subdir} DEFINITION headers)
  get_directory_property(sub_sources DIRECTORY ${subdir} DEFINITION sources)
  set(headers ${headers} ${sub_headers})
  set(sources ${sources} ${sub_sources})
endforeach()



