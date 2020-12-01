##############################################################################
#                                                                            #
#                                Nyx defines                                 #
#                                                                            #
##############################################################################
target_compile_definitions( ${EXE_NAME} PUBLIC $<BUILD_INTERFACE:CRSEGRNDOMP>) # Coarse grid OMP approach

if (Nyx_HENSON)
   target_compile_definitions( ${EXE_NAME} PUBLIC $<BUILD_INTERFACE:HENSON>)
endif ()

if (Nyx_GRAVITY)
   target_compile_definitions( ${EXE_NAME} PUBLIC $<BUILD_INTERFACE:GRAVITY>)
endif ()

if (NOT Nyx_HYDRO)
   target_compile_definitions( ${EXE_NAME} PUBLIC $<BUILD_INTERFACE:NO_HYDRO>)
endif()

if (Nyx_HEATCOOL)
   target_compile_definitions( ${EXE_NAME} PUBLIC $<BUILD_INTERFACE:HEATCOOL>)
endif()

if (Nyx_CONST_SPECIES)
   target_compile_definitions(${EXE_NAME} PUBLIC $<BUILD_INTERFACE:CONST_SPECIES>)
endif()

if (Nyx_NEUTRINO_PARTICLES)
   target_compile_definitions( ${EXE_NAME} PUBLIC $<BUILD_INTERFACE:NEUTRINO_PARTICLES>)
endif()

if (Nyx_NEUTRINO_DARK_PARTICLES)
   target_compile_definitions( ${EXE_NAME} PUBLIC $<BUILD_INTERFACE:NEUTRINO_DARK_PARTICLES>)
endif()

if (Nyx_SDC)
   target_compile_definitions( ${EXE_NAME} PUBLIC $<BUILD_INTERFACE:SDC>)
endif()

if (Nyx_FORCING)
   target_compile_definitions(${EXE_NAME} PUBLIC $<BUILD_INTERFACE:FORCING>)
endif ()

if (Nyx_SUNDIALS)
   target_compile_definitions(${EXE_NAME}
      PUBLIC $<BUILD_INTERFACE:AMREX_USE_SUNDIALS>)
endif ()

##############################################################################
#                                                                            #
#                          Nyx third party libraries                         #
#                                                                            #
##############################################################################
target_link_libraries(${EXE_NAME} PRIVATE AMReX::amrex)

if (Nyx_SUNDIALS)
   target_link_libraries(${EXE_NAME} PRIVATE SUNDIALS::cvode)
   if (Nyx_OMP)
      target_link_libraries(${EXE_NAME} PRIVATE SUNDIALS::nvecopenmp)
   endif ()
   if (Nyx_GPU_BACKEND STREQUAL CUDA)
      target_link_libraries(${EXE_NAME} PRIVATE SUNDIALS::nveccuda)
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
