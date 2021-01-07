target_include_directories(nyxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyxcore
   PRIVATE
   Godunov.cpp
   Hydro.cpp
   strang_hydro.cpp
   trace_ppm.cpp
   Utilities.cpp
   Godunov.H
   Hydro.H
   IndexDefines.H
   PLM.H
   PPM.H
   Utilities.H
   Riemann.H
   )

if (Nyx_SDC)
   target_sources( nyxcore
      PRIVATE
      sdc_hydro.cpp
      )
endif ()
