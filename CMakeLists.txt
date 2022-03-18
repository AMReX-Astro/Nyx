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
# Other options   =============================================================
#
cmake_dependent_option( Nyx_OMP  "Enable OpenMP"  ON
   "Nyx_GPU_BACKEND STREQUAL NONE" OFF)

cmake_dependent_option( Nyx_MPI  "Enable MPI"  ON
   "NOT Nyx_GPU_BACKEND STREQUAL SYCL" OFF)

cmake_dependent_option(Nyx_MPI_THREAD_MULTIPLE
   "Concurrent MPI calls from multiple threads" NO "Nyx_MPI" NO)

option(Nyx_SINGLE_PRECISION_PARTICLES
   "Use single precision particle data structures" YES )

option( Nyx_HYDRO "Run with baryon hydro solve and fields" YES)
option(Nyx_HEATCOOL "Run with H and He heating-cooling effects" NO)

option( BUILD_SHARED_LIBS "Build AMReX and SUNDIALS as shared libraries" OFF )
cmake_dependent_option(Nyx_GIMLET "" NO "Nyx_MPI;Nyx_OMP" NO)


option( Nyx_CONST_SPECIES "Don't evolve H and He, treat them as constant" ON)
option( Nyx_CGS   "Evolve quantities in CGS units instead of code units" NO)
option( Nyx_REEBER             "" NO)

#  These two are only really used by LyA_Neutrinos
option( Nyx_NEUTRINO_PARTICLES "" NO)
option( Nyx_NEUTRINO_DARK_PARTICLES "" NO)


#
# Internal options   ==========================================================
#
#
if (NOT Nyx_HEATCOOL )
   set(Nyx_SDC             NO CACHE INTERNAL "")
   set(Nyx_SUNDIALS        NO CACHE INTERNAL "")
   set(Nyx_SUNDIALS_FUSED  NO CACHE INTERNAL "")
else ()
   set(Nyx_SDC      YES CACHE INTERNAL "")
   set(Nyx_SUNDIALS YES CACHE INTERNAL "")
   if (Nyx_GPU_BACKEND STREQUAL "CUDA")
      set(Nyx_SUNDIALS_FUSED    YES CACHE INTERNAL "")
   else ()
      set(Nyx_SUNDIALS_FUSED    NO  CACHE INTERNAL "")
   endif ()
   if (Nyx_GPU_BACKEND STREQUAL "HIP")
      set(Nyx_SUNDIALS_FUSED    YES CACHE INTERNAL "")
   else ()
      set(Nyx_SUNDIALS_FUSED    NO  CACHE INTERNAL "")
   endif ()
endif ()

set(Nyx_AGN  NO  CACHE INTERNAL "AGN is broken at the moment")

#
# These used to be options but now they are set to true and should not
# be set by the user. We keep them around for legacy for the time being
#
set(Nyx_GRAVITY   YES CACHE INTERNAL "Enable gravity (should always be on)")
set(Nyx_PARTICLES YES CACHE INTERNAL "Enable particles (should always be on)")


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
#                                Nyx core library                            #
#                                                                            #
##############################################################################
add_subdirectory(Source)


##############################################################################
#                                                                            #
#                                Nyx executables                             #
#                                                                            #
##############################################################################
add_subdirectory(Exec)
