set(AMREX_MINIMUM_VERSION 20.11 CACHE INTERNAL "Minimum required AMReX version")

# We first check if we can find an AMReX installation.
# If so, we proceed with STANDALONE mode
# If not, we proceed with SUPERBUILD MODE
find_package( AMReX ${AMREX_MINIMUM_VERSION} CONFIG QUIET )

#
# ~~~~~~~~~~~~~~~~~~~~~~~ STANDALONE MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
if (AMReX_FOUND)

   # We are in this branch because an AMReX installation has been found.
   # In this scenario, we don't know how AMReX has been built. We only know
   # that some AMReX components, namely 3D, DOUBLE, PARTICLES
   # MUST be present in the installation.
   set(AMREX_REQUIRED_COMPONENTS 3D DOUBLE PARTICLES)

   if (Nyx_MPI)
      list(APPEND AMREX_REQUIRED_COMPONENTS MPI)
      if (Nyx_MPI_THREAD_MULTIPLE)
         list(APPEND AMREX_REQUIRED_COMPONENTS MPI_THREAD_MULTIPLE)
      endif ()
   endif ()
   if (Nyx_OMP)
      list(APPEND AMREX_REQUIRED_COMPONENTS OMP)
   endif ()
   if (Nyx_GPU_BACKEND STREQUAL CUDA)
      list(APPEND AMREX_REQUIRED_COMPONENTS CUDA)
   elseif (Nyx_GPU_BACKEND STREQUAL HIP)
      list(APPEND AMREX_REQUIRED_COMPONENTS HIP)
   elseif (Nyx_GPU_BACKEND STREQUAL SYCL)
      list(APPEND AMREX_REQUIRED_COMPONENTS SYCL)
   endif ()
   if (Nyx_GRAVITY)
      list(APPEND AMREX_REQUIRED_COMPONENTS LSOLVERS)
   endif ()
   if (Nyx_SINGLE_PRECISION_PARTICLES)
      list(APPEND AMREX_REQUIRED_COMPONENTS PSINGLE)
   else ()
      list(APPEND AMREX_REQUIRED_COMPONENTS PDOUBLE)
   endif ()

   # We now check again for the AMReX package.
   # This time we mark AMReX + its required components as REQUIRED.
   # If the AMReX installation does not contain all the required components,
   # the configuration step stops with an error message.
   find_package(AMReX ${AMREX_MINIMUM_VERSION} CONFIG
      REQUIRED ${AMREX_REQUIRED_COMPONENTS}
      )

   message(STATUS "AMReX found: configuration file located at ${AMReX_DIR}")

   # IS it worth checking this???
   if ( NOT ( "${CMAKE_BUILD_TYPE}" STREQUAL "${AMReX_BUILD_TYPE}" ) )
      message (WARNING "Nyx build type (${CMAKE_BUILD_TYPE}) type does not match AMReX build type (${AMReX_BUILD_TYPE})")
   endif ()

   # We load this here so we have the CUDA helper functions
   # available everywhere we need it
   if (ENABLE_CUDA)
      include(AMReXTargetHelpers)
   endif ()
else ()
#
# ~~~~~~~~~~~~~~~~~~~~~~~ SUPERBUILD MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

   message(STATUS "No AMReX installation found: cloning AMReX repo")

   if (NOT EXISTS  "${PROJECT_SOURCE_DIR}/.git")
      message(FATAL_ERROR
         "${PROJECT_SOURCE_DIR} is not a Git repo: missing .git")
   endif ()

   set(AMREX_SRC_DIR "${PROJECT_SOURCE_DIR}/subprojects/amrex"
     CACHE INTERNAL "Path to AMReX source (submodule)")

   if (NOT EXISTS "${AMREX_SRC_DIR}/.git")
      message(STATUS "Initializing git submodule for AMReX")

      find_package(Git REQUIRED)

      execute_process(
         COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive subprojects/amrex
         WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
         RESULT_VARIABLE GIT_SUBMOD_RESULT
         )

      if ( NOT GIT_SUBMOD_RESULT EQUAL "0")
         message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
      endif()

      unset(GIT_SUBMOD_RESULT)

   endif ()

   # Set build options for subproject
   set(AMReX_SPACEDIM             3                          CACHE INTERNAL "" )
   set(AMReX_PRECISION            DOUBLE                     CACHE INTERNAL "" )
   set(AMReX_FORTRAN              OFF                        CACHE INTERNAL "" )
   set(AMReX_PROBINIT             OFF                        CACHE INTERNAL "" )
   set(AMReX_MPI                  ${Nyx_MPI}                 CACHE INTERNAL "" )
   set(AMReX_MPI_THREAD_MULTIPLE  ${Nyx_MPI_THREAD_MULTIPLE} CACHE INTERNAL "" )
   set(AMReX_OMP                  ${Nyx_OMP}                 CACHE INTERNAL "" )
   set(AMReX_GPU_BACKEND          ${Nyx_GPU_BACKEND}         CACHE INTERNAL "" )
   set(AMReX_LINEAR_SOLVERS       ${Nyx_GRAVITY}             CACHE INTERNAL "" )
   set(AMReX_PARTICLES            ON                         CACHE INTERNAL "" )
   set(AMReX_SUNDIALS             ${Nyx_SUNDIALS}            CACHE INTERNAL "" )
   set(AMReX_BUILD_TUTORIALS      OFF                        CACHE INTERNAL "" )
   set(AMReX_INSTALL              OFF                        CACHE INTERNAL "" )
   if (Nyx_SINGLE_PRECISION_PARTICLES)
      set(AMReX_PARTICLES_PRECISION SINGLE CACHE INTERNAL "" )
   else ()
      set(AMReX_PARTICLES_PRECISION DOUBLE CACHE INTERNAL "" )
   endif ()

   list(APPEND CMAKE_MODULE_PATH ${AMREX_SRC_DIR}/Tools/CMake)

   # If CUDA is required, enable the language BEFORE adding the AMReX directory
   # Since AMReX_SetupCUDA has an include guard, it will be included only once here.
   # The reason for enabling CUDA before adding the AMReX subdirectory is that
   # the top-most directory needs to setup the CUDA language before a CUDA-enabled target
   # from a sub-project is included via add_subdirectory.
   # IMPORTANT: if you don't do this, AMReX will perform this step in a sub-scope and therefore
   # it will not setup CUDA here!
   if(Nyx_GPU_BACKEND STREQUAL CUDA)
      include(AMReX_SetupCUDA)
   endif ()

   add_subdirectory(${AMREX_SRC_DIR})
#   include(AMReXGenerateConfigHeader)  # provides add_amrex_defines
#   add_amrex_define(AMREX_USE_SUNDIALS NO_LEGACY IF AMReX_SUNDIALS)

endif ()
