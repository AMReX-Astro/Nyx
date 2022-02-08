##############################################################################
#                                                                            #
#                                Nyx core library                            #
#                                                                            #
##############################################################################
add_library(nyxcore OBJECT)


##############################################################################
#                                                                            #
#                                Nyx defines                                 #
#                                                                            #
##############################################################################
target_compile_definitions( nyxcore PUBLIC CRSEGRNDOMP) # Coarse grid OMP approach

if (Nyx_HENSON)
   target_compile_definitions( nyxcore PUBLIC HENSON)
endif ()

if (Nyx_GRAVITY)
   target_compile_definitions( nyxcore PUBLIC GRAVITY)
endif ()

if (NOT Nyx_HYDRO)
   target_compile_definitions( nyxcore PUBLIC NO_HYDRO)
endif()

if (Nyx_HEATCOOL)
   target_compile_definitions( nyxcore PUBLIC HEATCOOL)
endif()

if (Nyx_CONST_SPECIES)
   target_compile_definitions(nyxcore PUBLIC CONST_SPECIES)
endif()

if (Nyx_NEUTRINO_PARTICLES)
   target_compile_definitions( nyxcore PUBLIC NEUTRINO_PARTICLES)
endif ()

if (Nyx_NEUTRINO_DARK_PARTICLES)
   target_compile_definitions( nyxcore PUBLIC NEUTRINO_DARK_PARTICLES)
endif ()

if (Nyx_SDC)
   target_compile_definitions( nyxcore PUBLIC SDC)
endif()

##############################################################################
#                                                                            #
#                          Nyx  compile flags                                #
#                                                                            #
##############################################################################
if (NOT CMAKE_CXX_FLAGS)
   target_link_libraries(nyxcore PUBLIC AMReX::Flags_CXX)
endif ()

##############################################################################
#                                                                            #
#                          Nyx third party libraries                         #
#                                                                            #
##############################################################################
target_link_libraries(nyxcore PUBLIC AMReX::amrex)

if (Nyx_SUNDIALS)
   target_link_libraries(nyxcore PUBLIC SUNDIALS::cvode)
   if (Nyx_OMP)
      target_link_libraries(nyxcore PUBLIC SUNDIALS::nvecopenmp)
   endif ()
   if (Nyx_GPU_BACKEND STREQUAL CUDA)
     target_link_libraries(nyxcore PUBLIC SUNDIALS::nveccuda)
     target_link_libraries(nyxcore PUBLIC SUNDIALS::cvode_fused_cuda)
    endif ()
   if (Nyx_GPU_BACKEND STREQUAL HIP)
      target_link_libraries(nyxcore PUBLIC SUNDIALS::nvechip)
      target_link_libraries(nyxcore PUBLIC SUNDIALS::cvode_fused_hip)
   endif ()
   if (Nyx_GPU_BACKEND STREQUAL SYCL)
      target_link_libraries(nyxcore PUBLIC SUNDIALS::nvecsycl)
   endif ()
endif ()


##############################################################################
#                                                                            #
#                               Nyx sources                                  #
#                                                                            #
##############################################################################
if (Nyx_AGN)
   add_subdirectory(AGN)
endif ()

add_subdirectory(DerivedQuantities)

add_subdirectory(Driver)

add_subdirectory(Forcing)

if (Nyx_GRAVITY)
   add_subdirectory(Gravity)
endif ()

if (Nyx_HYDRO)
   add_subdirectory(HeatCool)
   add_subdirectory(Hydro)
   add_subdirectory(EOS)
endif ()

add_subdirectory(Initialization)

add_subdirectory(IO)

if (Nyx_MPI)
   add_subdirectory(Monitors)
endif ()

add_subdirectory(Particle)

add_subdirectory(TimeStep)


##############################################################################
#                                                                            #
#                               Nyx CUDA and buildInfo                       #
#                                                                            #
##############################################################################
if (Nyx_GPU_BACKEND STREQUAL "CUDA")
   setup_target_for_cuda_compilation( nyxcore )
endif ()

include(AMReXBuildInfo)
generate_buildinfo(nyxcore ${CMAKE_SOURCE_DIR} REQUIRED)
