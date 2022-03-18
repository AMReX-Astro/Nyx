set(SUNDIALS_MINIMUM_VERSION 6.0.0 CACHE INTERNAL "Minimum required SUNDIALS version")

# We first check if we can find an AMReX installation.
# If so, we proceed with STANDALONE mode
# If not, we proceed with SUPERBUILD MODE
find_package( SUNDIALS ${SUNDIALS_MINIMUM_VERSION} CONFIG QUIET )

#
# ~~~~~~~~~~~~~~~~~~~~~~~ STANDALONE MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
if (SUNDIALS_FOUND)

   # Unfortunately SUNDIALS config file doesn't seem to allow to check for
   # components, so we assume that the SUNDIALS installation, if found,
   # has all the necessary components
   message(STATUS "SUNDIALS found: configuration file located at ${SUNDIALS_DIR}")

else ()
#
# ~~~~~~~~~~~~~~~~~~~~~~~ SUPERBUILD MODE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

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

   # Set build options for subproject
   set(EXAMPLES_ENABLE_C            OFF                        CACHE INTERNAL "" )
   set(EXAMPLES_ENABLE_CXX          OFF                        CACHE INTERNAL "" )
   set(ENABLE_MPI                   OFF                        CACHE INTERNAL "" )
   set(ENABLE_OPENMP                ${Nyx_OMP}                 CACHE INTERNAL "" )
   if (Nyx_GPU_BACKEND STREQUAL CUDA)
     set(ENABLE_CUDA                  ON                      CACHE INTERNAL "" )
     set(SUNDIALS_INDEX_SIZE          32                      CACHE INTERNAL "" )
     set(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS     ON          CACHE INTERNAL "" )
    else ()
      set(ENABLE_CUDA                  OFF                     CACHE INTERNAL "" )
      if (Nyx_GPU_BACKEND STREQUAL HIP)
	set(ENABLE_HIP                   ON                      CACHE INTERNAL "" )
        set(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS     ON          CACHE INTERNAL "" )
      endif ()
      if (Nyx_GPU_BACKEND STREQUAL SYCL)
	set(ENABLE_SYCL                  ON                      CACHE INTERNAL "" )
	set(CMAKE_CXX_STANDARD           17                      CACHE INTERNAL "" )
      endif ()
   endif ()
   set(BUILD_ARKODE                 OFF                        CACHE INTERNAL "" )
   set(BUILD_KINSOL                 OFF                        CACHE INTERNAL "" )
   set(BUILD_IDA                    OFF                        CACHE INTERNAL "" )
   set(BUILD_IDAS                   OFF                        CACHE INTERNAL "" )
   set(BUILD_CVODES                 OFF                        CACHE INTERNAL "" )
   set(BUILD_TESTING OFF)

   #  Add  SUNDIALS sources to the build
   add_subdirectory(${SUNDIALS_SRC_DIR})

   # This is to use the same target name uses by the sundials exported targets
   if (BUILD_SHARED_LIBS)
   add_library(SUNDIALS::cvode      ALIAS sundials_cvode_shared)
   if (Nyx_OMP)
      add_library(SUNDIALS::nvecopenmp ALIAS sundials_nvecopenmp_shared)
   endif ()
   if (Nyx_GPU_BACKEND STREQUAL CUDA)
      add_library(SUNDIALS::nveccuda ALIAS sundials_nveccuda_shared)
      add_library(SUNDIALS::cvode_fused_cuda ALIAS sundials_cvode_fused_cuda_shared)
      get_target_property(SUNDIALS_INCLUDES sundials_nveccuda_shared SOURCES)
      list(FILTER SUNDIALS_INCLUDES INCLUDE REGEX "\\.h")
      set_target_properties(
      sundials_nveccuda_shared PROPERTIES PUBLIC_HEADER "${SUNDIALS_INCLUDES}")
   endif ()
   if (Nyx_GPU_BACKEND STREQUAL HIP)
     add_library(SUNDIALS::nvechip ALIAS sundials_nvechip_shared)
     add_library(SUNDIALS::cvode_fused_hip ALIAS sundials_cvode_fused_hip_shared)
     get_target_property(SUNDIALS_INCLUDES sundials_nvechip_shared SOURCES)
     list(FILTER SUNDIALS_INCLUDES INCLUDE REGEX "\\.h")
     set_target_properties(
     sundials_nvechip_shared PROPERTIES PUBLIC_HEADER "${SUNDIALS_INCLUDES}")
   endif ()
   if (Nyx_GPU_BACKEND STREQUAL SYCL)
      add_library(SUNDIALS::nvecsycl ALIAS sundials_nvecsycl_shared)
   endif ()
   else ()
      add_library(SUNDIALS::cvode      ALIAS sundials_cvode_static)
   if (Nyx_OMP)
      add_library(SUNDIALS::nvecopenmp ALIAS sundials_nvecopenmp_static)
   endif ()
   if (Nyx_GPU_BACKEND STREQUAL CUDA)
      add_library(SUNDIALS::nveccuda ALIAS sundials_nveccuda_static)
      add_library(SUNDIALS::cvode_fused_cuda ALIAS sundials_cvode_fused_cuda_static)
      get_target_property(SUNDIALS_INCLUDES sundials_nveccuda_static SOURCES)
      list(FILTER SUNDIALS_INCLUDES INCLUDE REGEX "\\.h")
      set_target_properties(
      sundials_nveccuda_static PROPERTIES PUBLIC_HEADER "${SUNDIALS_INCLUDES}")
   endif ()
   if (Nyx_GPU_BACKEND STREQUAL HIP)
     add_library(SUNDIALS::nvechip ALIAS sundials_nvechip_static)
     add_library(SUNDIALS::cvode_fused_hip ALIAS sundials_cvode_fused_hip_static)
     get_target_property(SUNDIALS_INCLUDES sundials_nvechip_static SOURCES)
     list(FILTER SUNDIALS_INCLUDES INCLUDE REGEX "\\.h")
     set_target_properties(
     sundials_nvechip_static PROPERTIES PUBLIC_HEADER "${SUNDIALS_INCLUDES}")
   endif ()
   if (Nyx_GPU_BACKEND STREQUAL SYCL)
      add_library(SUNDIALS::nvecsycl ALIAS sundials_nvecsycl_static)
   endif ()
   endif ()

   set(SUNDIALS_FOUND TRUE)


endif ()
