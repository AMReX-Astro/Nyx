target_include_directories(nyxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyxcore
   PRIVATE
   nyx_main.cpp
   Nyx.cpp
   NyxBld.cpp
   comoving.cpp
   Nyx_error.cpp
   )

target_sources( nyxcore
   PRIVATE
   sum_integrated_quantities.cpp
   sum_utils.cpp
   Nyx.H
   Prob.H
   )

if (Nyx_CGS)
   target_sources( nyxcore
      PRIVATE
      constants_cosmo.H
      )
endif ()

if (Nyx_REBEER)
   target_sources( nyxcore PRIVATE Nyx_halos.cpp)
endif ()

if (Nyx_GIMLET)
   target_sources( nyxcore
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
# Don't include this file in the source lists for nyxcore
# Final executable will include it if necessary
#
set( NYX_DEFAULT_MAIN ${CMAKE_CURRENT_LIST_DIR}/main.cpp
   CACHE INTERNAL "Default main if none other is provided")
