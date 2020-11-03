cmake_minimum_required(VERSION 3.14)

project(Nyx
   # VERSION <version>
   DESCRIPTION  "A cosmological hydrodynamics/N-body code"
   HOMEPAGE_URL "https://amrex-astro.github.io/Nyx/"
   LANGUAGES    C CXX
   )

# Set the search path for cmake modules
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake )

#
# This sets CMAKE_CXX_FLAGS_<CONFIG> to a default value
# if the variable is empty
#
if( NOT CMAKE_CXX_FLAGS_DEBUG )
   set(CMAKE_CXX_FLAGS_DEBUG "-g")
endif()

if( NOT CMAKE_CXX_FLAGS_RELEASE )
   set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
endif()

#
# Set default build type to Release
#
if ( NOT CMAKE_BUILD_TYPE )
   message(STATUS "Setting build type to Release as none was specified.")
   set( CMAKE_BUILD_TYPE Release )
else ()
   message(STATUS "Build type set by user to '${CMAKE_BUILD_TYPE}'.")
endif()

##############################################################################
#                                                                            #
#                         Nyx options                                        #
#                                                                            #
##############################################################################
option( Nyx_Fortran    "Enable Fortran" OFF)   # This will go
option( Nyx_MPI        "Enable MPI"     YES)
option( Nyx_OMP        "Enable OpenMP"  NO)
option( Nyx_CUDA       "Enable CUDA"    NO)
option( Nyx_PARTICLES  "Enable particles support" OFF)

option( Nyx_REEBER             "" OFF)
option( Nyx_HYDRO              "" ON)
option( Nyx_SDC                "" OFF)
option( Nyx_GIMLET             "" OFF)
option( Nyx_AGN                "" OFF)
option( Nyx_CGS                "" OFF)
option( Nyx_HEATCOOL           "" OFF)
option( Nyx_FORCING            "" OFF)
option( Nyx_GRAVITY            "" OFF)
option( Nyx_CVODE              "" OFF)
option( Nyx_NEUTRINO_PARTICLES "" OFF)

# option( Nyx_FPE    "Enable Floating Point Exceptions checks" NO )

if (Nyx_CUDA)
    enable_language(CUDA)
endif ()

##############################################################################
#                                                                            #
#                         Nyx third party libraries                          #
#                                                                            #
##############################################################################
include(NyxSetupAMReX)


##############################################################################
#                                                                            #
#                                Nyx target                                  #
#                                                                            #
##############################################################################
add_executable(nyx)
add_subdirectory(Source)
target_link_libraries(nyx PRIVATE AMReX::amrex)


##############################################################################
#                                                                            #
#                                Nyx defines                                 #
#                                                                            #
##############################################################################
if (NOT Nyx_HYDRO)
   target_compile_definitions(nyx PUBLIC $<BUILD_INTERFACE:NO_HYDRO>)
endif()

if (NOT Nyx_HEATCOOL)
   target_compile_definitions(nyx PUBLIC $<BUILD_INTERFACE:HEATCOOL>)
endif()

# if (NOT Nyx_CONST_SPECIES)
#    target_compile_definitions(nyx PUBLIC $<BUILD_INTERFACE:CONST_SPECIES>)
# endif()

if (NOT Nyx_NEUTRINO_PARTICLES)
   target_compile_definitions(nyx PUBLIC $<BUILD_INTERFACE:NEUTRINO_PARTICLES>)
endif()

# if (NOT Nyx_NEUTRINO_DARK_PARTICLES)
#    target_compile_definitions(nyx PUBLIC $<BUILD_INTERFACE:NEUTRINO_DARK_PARTICLES>)
# endif()

if (NOT Nyx_SDC)
   target_compile_definitions(nyx PUBLIC $<BUILD_INTERFACE:SDC>)
endif()



# get_git_info( )

# # Set some cache variable for the whole project
# set( MAIN ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp CACHE INTERNAL "Main source file")

# if(ENABLE_CSG)
#   add_subdirectory(subprojects/csg-eb)
# endif()

# #
# # Core MFIX library
# #
# add_subdirectory(src)

# # Print out config summary
# print_mfix_configuration_summary(mfixcore)

# add_executable(mfix ${MAIN})
# target_link_libraries( mfix mfixcore )

# # We should check what this line does
# file( GLOB USR_OVERRIDES ${CMAKE_CURRENT_BINARY_DIR}/*.cpp)
# add_udfs_to_target(mfix ${USR_OVERRIDES})

# if (ENABLE_CUDA)
#    setup_target_for_cuda_compilation(mfix)
# endif ()

# # Tools
# set(MFIX_TOOLS_DIR ${PROJECT_BINARY_DIR}/tools/ CACHE INTERNAL "Path to tools directory")
# add_subdirectory(tools)

# # Tests
# enable_testing()
# add_subdirectory( tests )

# # Tags
# include(MFIX_Tags)
# add_tags_targets()

# # Add a target to generate API documentation with Doxygen
# find_package(Doxygen)
# if (DOXYGEN_FOUND)
#    set(DOC_DIR ${CMAKE_SOURCE_DIR}/doc)
#    configure_file(${DOC_DIR}/Doxyfile.in ${DOC_DIR}/Doxyfile)
#    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${DOC_DIR}/Doxyfile
#       WORKING_DIRECTORY ${DOC_DIR}
#       COMMENT "Generating API documentation with Doxygen" VERBATIM
#       )
# endif (DOXYGEN_FOUND)

# # Copy compile_commands.json from the build directory to the project root directory
# if ( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
#     execute_process( COMMAND ${CMAKE_COMMAND} -E copy_if_different
#         ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
#         ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
#   )
# endif ()

# if(ENABLE_CSG)
#   target_link_libraries(mfixcore PRIVATE csg)
# endif()
