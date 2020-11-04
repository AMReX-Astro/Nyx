# include(AMReXBuildInfo)

# #
# # Function to setup the Nyx executables
# #
# function (setup_nyx_executable _srcs _inputs)

#    cmake_parse_arguments( "" "HAS_FORTRAN_MODULES"
#       "BASE_NAME;RUNTIME_SUBDIR;EXTRA_DEFINITIONS;MAIN" "" ${ARGN} )

#    if (_BASE_NAME)
#       set(_base_name ${_BASE_NAME})
#    else ()
#       string(REGEX REPLACE ".*Exec/" "" _base_name ${CMAKE_CURRENT_LIST_DIR})
#       string(REPLACE "/" "_" _base_name ${_base_name})
#    endif ()

#    # Prepend "nyx_" to base name
#    set(_base_name "nyx_${_base_name}")

#    if (_RUNTIME_SUBDIR)
#       set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR}/${_RUNTIME_SUBDIR})
#    else ()
#       set(_exe_dir ${CMAKE_CURRENT_BINARY_DIR})
#    endif ()

#    set( _exe_name  ${_base_name} )

#    add_executable( ${_exe_name} $<TARGET_OBJECTS:nyx> )

#    target_sources( ${_exe_name} PRIVATE ${${_srcs}} )

#    if (_MAIN)
#       target_sources( ${_exe_name} PRIVATE ${_MAIN} )
#    else ()
#       target_sources( ${_exe_name} PRIVATE ${NYX_DEFAULT_MAIN} )
#    endif ()

#    set_target_properties( ${_exe_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${_exe_dir} )

#    if (_EXTRA_DEFINITIONS)
#       target_compile_definitions(${_exe_name} PRIVATE ${_EXTRA_DEFINITIONS})
#    endif ()

#    # Find out which include directory is needed
#    set(_includes ${${_srcs}})
#    list(FILTER _includes INCLUDE REGEX "\\.H$")
#    foreach(_item IN LISTS _includes)
#       get_filename_component( _include_dir ${_item} DIRECTORY)
#       if (NOT _include_dir)
#          set(_include_dir ${CMAKE_CURRENT_LIST_DIR})
#       endif ()
#       target_include_directories( ${_exe_name} PRIVATE  ${_include_dir} )
#    endforeach()

#    generate_buildinfo(${_exe_name} ${CMAKE_SOURCE_DIR} REQUIRED)

#    # if (_HAS_FORTRAN_MODULES)
#    #    target_include_directories(${_exe_name}
#    #       PRIVATE
#    #       ${CMAKE_CURRENT_BINARY_DIR}/mod_files)
#    #    set_target_properties( ${_exe_name}
#    #       PROPERTIES
#    #       Fortran_MODULE_DIRECTORY
#    #       ${CMAKE_CURRENT_BINARY_DIR}/mod_files )
#    # endif ()

#    # target_link_libraries( ${_exe_name}  )

#    if (AMReX_CUDA)
#       setup_target_for_cuda_compilation( ${_exe_name} )
#    endif ()

#    if (${_inputs})
#       file( COPY ${${_inputs}} DESTINATION ${_exe_dir} )
#    endif ()

# endfunction ()

include(NyxSetupExecutable)
add_subdirectory(AMR-density)
