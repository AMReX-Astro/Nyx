target_include_directories(nyx
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyx
   PRIVATE
   main.cpp       # should I include these two here???
   nyx_main.cpp   #
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
   target_sources( nyx PRIVATE Nyx_halos.cpp )
endif ()


if (Nyx_HYDRO)

   target_sources( nyx
      PRIVATE
      compute_hydro_sources.cpp
      update_state_with_sources.cpp
      strang_hydro.cpp
      strang_reactions.cpp
      sum_integrated_quantities.cpp
      )

   if (Nyx_SDC)
      target_sources( nyx
         PRIVATE
         sdc_hydro.cpp
         sdc_reactions.cpp
         )
   endif ()

endif ()

if (Nyx_GIMLET)
   target_sources( nyx
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
