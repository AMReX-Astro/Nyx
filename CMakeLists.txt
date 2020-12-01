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
option( Nyx_OMP        "Enable OpenMP"  YES)
option( Nyx_MPI        "Enable MPI"     YES)
cmake_dependent_option(Nyx_MPI_THREAD_MULTIPLE
   "Concurrent MPI calls from multiple threads" NO "Nyx_MPI" NO)


#
# GPU backends    =============================================================
#
set(Nyx_GPU_BACKEND_VALUES NONE SYCL CUDA HIP)
set(Nyx_GPU_BACKEND NONE CACHE STRING "On-node, accelerated GPU backend: <NONE,SYCL,CUDA,HIP>")
set_property(CACHE Nyx_GPU_BACKEND PROPERTY STRINGS ${Nyx_GPU_BACKEND_VALUES})
if (NOT Nyx_GPU_BACKEND IN_LIST Nyx_GPU_BACKEND_VALUES)
   message(FATAL_ERROR "Nyx_GPU_BACKEND=${Nyx_GPU_BACKEND} is not allowed."
      " Must be one of ${Nyx_GPU_BACKEND_VALUES}")
endif ()

if (NOT Nyx_GPU_BACKEND STREQUAL NONE)
   message( STATUS "   Nyx_GPU_BACKEND = ${Nyx_GPU_BACKEND}")
endif ()


#
# These used to be options but now they are set to true and should not
# be set by the user. We keep them around for legacy for the time being
#
set(Nyx_GRAVITY   YES CACHE INTERNAL "Enable gravity (should always be on)")
set(Nyx_PARTICLES YES CACHE INTERNAL "Enable particles (should always be on)")


#
#
#
option(Nyx_SINGLE_PRECISION_PARTICLES
   "Use single precision particle data structures" YES )


option( Nyx_HYDRO NO "")
cmake_dependent_option(Nyx_HEATCOOL "" NO
   "Nyx_HYDRO; NOT Nyx_GPU_BACKEND STREQUAL SYCL; NOT Nyx_GPU_BACKEND STREQUAL HIP" NO)

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

option( BUILD_SHARED_LIBS "Build AMReX and SUNDIALS as shared libraries" OFF )

# option( Nyx_AGN                "" NO)
option( Nyx_SDC                "" NO)

# Sundials options
cmake_dependent_option(Nyx_SUNDIALS "Use SUNDIALS CVODE libs" NO
   "NOT Nyx_GPU_BACKEND STREQUAL SYCL; NOT Nyx_GPU_BACKEND STREQUAL HIP" NO)
cmake_dependent_option(Nyx_SUNDIALS_FUSED "Use SUNDIALS FUSED libs" NO
   "Nyx_SUNDIALS" NO)

cmake_dependent_option(Nyx_GIMLET "" NO "Nyx_MPI;Nyx_OMP" NO)


option( Nyx_CONST_SPECIES      "" ON)
option( Nyx_REEBER             "" NO)
option( Nyx_CGS                "" NO)
option( Nyx_FORCING            "" NO)
# option( Nyx_NEUTRINO_PARTICLES "" NO)
# option( Nyx_NEUTRINO_DARK_PARTICLES "" OFF)


if (Nyx_GPU_BACKEND STREQUAL CUDA)
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
