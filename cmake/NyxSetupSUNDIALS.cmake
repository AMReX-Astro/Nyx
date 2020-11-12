#
# This module
#

set(SUNDIALS_MINIMUM_VERSION 5.5.0 CACHE INTERNAL "Minimum required SUNDIALS version")

# We first check if we can find an AMReX installation.
# If so, we proceed with STANDALONE mode
# If not, we proceed with SUPERBUILD MODE
# find_package( AMReX ${AMREX_MINIMUM_VERSION} CONFIG QUIET )

# #
# # ~~~~~~~~~~~~~~~~~~~~~~~ STANDALONE MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #
# if (AMReX_FOUND)

#    # We are in this branch because an AMReX installation has been found.
#    # In this scenario, we don't know how AMReX has been built. We only know
#    # that some AMReX components, namely 3D, DOUBLE, PARTICLES
#    # MUST be present in the installation.
#    set(AMREX_REQUIRED_COMPONENTS 3D DOUBLE PARTICLES)

#    if (Nyx_Fortran)
#       list(APPEND AMREX_REQUIRED_COMPONENTS FORTRAN)
#    endif ()
#    if (Nyx_MPI)
#       list(APPEND AMREX_REQUIRED_COMPONENTS MPI)
#       if (Nyx_MPI_THREAD_MULTIPLE)
#          list(APPEND AMREX_REQUIRED_COMPONENTS MPI_THREAD_MULTIPLE)
#       endif ()
#    endif ()
#    if (Nyx_OMP)
#       list(APPEND AMREX_REQUIRED_COMPONENTS OMP)
#    endif ()
#    if (Nyx_CUDA)
#       list(APPEND AMREX_REQUIRED_COMPONENTS CUDA)
#    endif ()
#    if (Nyx_GRAVITY)
#       list(APPEND AMREX_REQUIRED_COMPONENTS LSOLVERS)
#    endif ()
#    if (Nyx_SUNDIALS)
#       list(APPEND AMREX_REQUIRED_COMPONENTS SUNDIALS)
#    endif ()
#    if (Nyx_SINGLE_PRECISION_PARTICLES)
#       list(APPEND AMREX_REQUIRED_COMPONENTS PSINGLE)
#    else ()
#       list(APPEND AMREX_REQUIRED_COMPONENTS PDOUBLE)
#    endif ()

#    # We now check again for the AMReX package.
#    # This time we mark AMReX + its required components as REQUIRED.
#    # If the AMReX installation does not contain all the required components,
#    # the configuration step stops with an error message.
#    find_package(AMReX ${AMREX_MINIMUM_VERSION} CONFIG
#       REQUIRED ${AMREX_REQUIRED_COMPONENTS}
#       )

#    message(STATUS "AMReX found: configuration file located at ${AMReX_DIR}")

#    # IS it worth checking this???
#    if ( NOT ( "${CMAKE_BUILD_TYPE}" STREQUAL "${AMReX_BUILD_TYPE}" ) )
#       message (WARNING "Nyx build type (${CMAKE_BUILD_TYPE}) type does not match AMReX build type (${AMReX_BUILD_TYPE})")
#    endif ()

#    # We load this here so we have the CUDA helper functions
#    # available everywhere we need it
#    if (ENABLE_CUDA)
#       include(AMReXTargetHelpers)
#    endif ()
# else ()
# #
# # ~~~~~~~~~~~~~~~~~~~~~~~ SUPERBUILD MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #

   message(STATUS "No SUNDIALS installation found: cloning SUNDIALS repo")

   if (NOT EXISTS  "${PROJECT_SOURCE_DIR}/.git")
      message(FATAL_ERROR
         "${PROJECT_SOURCE_DIR} is not a Git repo: missing .git")
   endif ()

   set(SUNDIALS_SRC_DIR "${PROJECT_SOURCE_DIR}/subprojects/sundials"
     CACHE INTERNAL "Path to SUNDIALS source (submodule)")

   if (NOT EXISTS "${SUNDIALS_SRC_DIR}/.git")
      message(STATUS "Initializing git submodule for SUNDIALS")

      find_package(Git REQUIRED)

      execute_process(
         COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive subprojects/sundials
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         RESULT_VARIABLE GIT_SUBMOD_RESULT
         )

      if ( NOT GIT_SUBMOD_RESULT EQUAL "0")
         message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
      endif()

      unset(GIT_SUBMOD_RESULT)

   endif ()
   message("Nyx_ARKODE=${Nyx_ARKODE}")

   # Set build options for subproject
   set(EXAMPLES_ENABLE_C            OFF                        CACHE INTERNAL "" )
   set(EXAMPLES_ENABLE_CXX          OFF                        CACHE INTERNAL "" )
   set(ENABLE_MPI                   ${Nyx_MPI}                 CACHE INTERNAL "" )
   set(ENABLE_OPENMP                ${Nyx_OMP}                 CACHE INTERNAL "" )
   set(ENABLE_CUDA                  ${Nyx_CUDA}                CACHE INTERNAL "" )
   set(BUILD_ARKODE                 ${Nyx_ARKODE}              CACHE INTERNAL "" )
   set(BUILD_KINSOL                 OFF                        CACHE INTERNAL "" )
   set(BUILD_IDA                    OFF                        CACHE INTERNAL "" )
   set(BUILD_IDAS                   OFF                        CACHE INTERNAL "" )
   set(BUILD_CVODES                 OFF                        CACHE INTERNAL "" )
   set(BUILD_TESTING OFF)



   add_subdirectory(${SUNDIALS_SRC_DIR})

   add_library(SUNDIALS::cvode ALIAS sundials_cvode_static)

   if (Nyx_ARKODE)
      add_library(SUNDIALS::arkode ALIAS sundials_arkode_static)
   endif ()
# endif ()
