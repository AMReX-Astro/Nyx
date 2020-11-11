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
include(CMakeDependentOption)

option( Nyx_Fortran    "Enable Fortran" NO)   # This will go away
option( Nyx_MPI        "Enable MPI"     YES)

cmake_dependent_option(Nyx_MPI_THREAD_MULTIPLE
   "Concurrent MPI calls from multiple threads"
   NO "Nyx_MPI" NO)

option( Nyx_OMP        "Enable OpenMP"  YES)
option( Nyx_CUDA       "Enable CUDA"    NO)
option( Nyx_PARTICLES  "Enable particles support" YES)

cmake_dependent_option(Nyx_SINGLE_PRECISION_PARTICLES
   "Use single precision particle data structures"
   YES "Nyx_PARTICLES" NO)
message(" ================================ Nyx_SINGLE_PRECISION_PARTICLES = ${Nyx_SINGLE_PRECISION_PARTICLES}")
option( Nyx_HEATCOOL   "" NO)

cmake_dependent_option(Nyx_HYDRO "" YES "NOT Nyx_HEATCOOL" YES)

# Define normal variables to overwrite values of some options define below
# if necessary
if (NOT Nyx_HYDRO)
   message(STATUS "Nyx_HYDRO set to ${Nyx_HYDRO}:"
      " disabling Nyx_SDC, Nyx_SUNDIALS, and Nyx_FUSED")
   set(Nyx_SDC NO)
   set(Nyx_SUNDIALS NO)
   set(Nyx_FUSED NO)
endif ()

if (Nyx_HEATCOOL)
   message(STATUS "Nyx_HEATCOOL set to ${Nyx_HEATCOOL}:"
      " enabling Nyx_SDC and Nyx_SUNDIALS")
   set(Nyx_SDC YES)
   set(Nyx_SUNDIALS YES)
   ###set(Nyx_FUSED NO)
endif ()

# option( Nyx_AGN                "" NO)
option( Nyx_SDC                "" NO)
message(" ================================ Nyx_SDC = ${Nyx_SDC}")
option( Nyx_SUNDIALS           "" NO)
option( Nyx_FUSED              "" NO)

option( Nyx_CONST_SPECIES      "" ON)
option( Nyx_REEBER             "" NO)
option( Nyx_GIMLET             "" NO)
option( Nyx_CGS                "" NO)
option( Nyx_FORCING            "" NO)
option( Nyx_GRAVITY            "" YES)
option( Nyx_CVODE              "" NO)
# option( Nyx_NEUTRINO_PARTICLES "" NO)
# option( Nyx_NEUTRINO_DARK_PARTICLES "" OFF)

if (Nyx_CUDA)
    enable_language(CUDA)
endif ()

##############################################################################
#                                                                            #
#                         Nyx third party libraries                          #
#                                                                            #
##############################################################################
if (Nyx_SUNDIALS)
   include(NyxSetupSUNDIALS)
endif ()
include(NyxSetupAMReX)


##############################################################################
#                                                                            #
#                                Nyx executables                             #
#                                                                            #
##############################################################################
add_subdirectory(Exec)

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
