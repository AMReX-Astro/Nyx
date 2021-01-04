target_include_directories(${EXE_NAME}
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( ${EXE_NAME}
   PRIVATE
   nyx_main.cpp
   Nyx.cpp
   NyxBld.cpp
   comoving.cpp
   Nyx_error.cpp
   )

target_sources( ${EXE_NAME}
   PRIVATE
   sum_integrated_quantities.cpp
   sum_utils.cpp
   Nyx.H
   Prob.H
   )


if (Nyx_CGS)
   target_sources( ${EXE_NAME}
      PRIVATE
      constants_cosmo.H
      )
endif ()

if (DEFAULT_MAIN)
   target_sources( ${EXE_NAME} PRIVATE main.cpp)
endif ()

if (Nyx_REBEER)
   target_sources( ${EXE_NAME} PRIVATE Nyx_halos.cpp)
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
