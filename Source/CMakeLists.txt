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
      PUBLIC $<BUILD_INTERFACE:AMREX_USE_SUNDIALS_3x4x;AMREX_USE_SUNDIALS>)
   if (Nyx_ARKODE)
      target_compile_definitions(${EXE_NAME}
         PUBLIC $<BUILD_INTERFACE:AMREX_USE_ARKODE_LIBS>)
   endif ()
endif ()

##############################################################################
#                                                                            #
#                          Nyx third party libraries                         #
#                                                                            #
##############################################################################
target_link_libraries(${EXE_NAME} PRIVATE AMReX::amrex)

if (Nyx_SUNDIALS)
   target_link_libraries(${EXE_NAME} PRIVATE SUNDIALS::cvode)
   if (Nyx_ARKODE)
      target_link_libraries(${EXE_NAME} PRIVATE SUNDIALS::arkode)
   endif ()
endif ()

if (DEFAULT_MAIN)
   target_sources( ${EXE_NAME} PRIVATE main.cpp)
endif ()


target_include_directories(${EXE_NAME}
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( ${EXE_NAME}
   PRIVATE
   nyx_main.cpp
   advance_particles.cpp
   Nyx.cpp
   Nyx_advance.cpp
   NyxBld.cpp
   comoving.cpp
   sum_utils.cpp
   Nyx.H
   IndexDefines.H
   )

if (Nyx_REBEER)
   target_sources( ${EXE_NAME} PRIVATE Nyx_halos.cpp )
endif ()


if (Nyx_HYDRO)

   target_sources( ${EXE_NAME}
      PRIVATE
      compute_hydro_sources.cpp
      update_state_with_sources.cpp
      strang_hydro.cpp
      strang_reactions.cpp
      sum_integrated_quantities.cpp
      )

   if (Nyx_SDC)
      target_sources( ${EXE_NAME}
         PRIVATE
         sdc_hydro.cpp
         sdc_reactions.cpp
         )
   endif ()

endif ()

if (Nyx_GIMLET)
   target_sources( ${EXE_NAME}
      PRIVATE
      DoGimletAnalysis.cpp
      postprocess_tau_fields.cpp
      MakePencilBoxes.cpp
      MakeFFTWBoxes.cpp
      temperature_density_pdf2d.cpp
      pdf.cpp
      ps3d.cpp
      )
endif ()



#
#  Add sources and header from current dir and all sub dirs
#
if (Nyx_AGN)
   add_subdirectory(AGN)
endif ()

if (NOT Nyx_CGS)
   add_subdirectory(Constants)
endif ()

add_subdirectory(DerivedQuantities)

add_subdirectory(EOS)

if (Nyx_FORCING)
   add_subdirectory(Forcing)
endif ()

if (Nyx_GRAVITY)
   add_subdirectory(Gravity)
endif ()

add_subdirectory(HeatCool)

add_subdirectory(HydroCpp)

add_subdirectory(Initialization)

add_subdirectory(IO)

if (Nyx_MPI)
   add_subdirectory(Monitors)
endif ()

add_subdirectory(Particle)

if (Nyx_HYDRO)
   add_subdirectory(SourceTerms)
endif ()


add_subdirectory(Src_3d)

add_subdirectory(Tagging)
