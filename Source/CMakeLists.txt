##############################################################################
#                                                                            #
#                                Nyx core library                            #
#                                                                            #
##############################################################################
add_library(nyxcore OBJECT)

# Variables to store genex for conditional compilation based on
# consuming target properties.
# Properties that the consuming target can set are:
#
# DEFAULT_MAIN
# AGN
# NEUTRINO_PARTICLES
# NEUTRINO_DARK_PARTICLES
#
set(_use_default_main  $<BOOL:$<TARGET_PROPERTY:DEFAULT_MAIN>>)
set(_has_agn           $<BOOL:$<TARGET_PROPERTY:AGN>>)
set(_has_neutrino      $<BOOL:$<TARGET_PROPERTY:NEUTRINO_PARTICLES>>)
set(_has_neutrino_dark $<BOOL:$<TARGET_PROPERTY:NEUTRINO_DARK_PARTICLES>>)


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

target_compile_definitions( nyxcore PUBLIC $<${_has_neutrino}:NEUTRINO_PARTICLES>)
target_compile_definitions( nyxcore PUBLIC $<${_has_neutrino_dark}:NEUTRINO_DARK_PARTICLES>)

if (Nyx_SDC)
   target_compile_definitions( nyxcore PUBLIC SDC)
endif()

if (Nyx_FORCING)
   target_compile_definitions( nyxcore PUBLIC FORCING)
endif ()

if (Nyx_SUNDIALS)
   target_compile_definitions( nyxcore PUBLIC AMREX_USE_SUNDIALS)
endif ()

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
   endif ()
endif ()


##############################################################################
#                                                                            #
#                               Nyx sources                                  #
#                                                                            #
##############################################################################
add_subdirectory(AGN)

add_subdirectory(DerivedQuantities)

add_subdirectory(Driver)

if (Nyx_FORCING)
   add_subdirectory(Forcing)
endif ()

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
