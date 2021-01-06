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
# At compile time, this will include main.cpp only if
# the consuming target has property DEFAULT_TARGET turned on
# Note that I need to provide the full path to the main file
#
set(_use_default_main $<BOOL:$<TARGET_PROPERTY:DEFAULT_MAIN>>)
target_sources( nyxcore
   PUBLIC
   $<${_use_default_main}:${CMAKE_CURRENT_LIST_DIR}/main.cpp>
   )
unset(_use_default_main)
