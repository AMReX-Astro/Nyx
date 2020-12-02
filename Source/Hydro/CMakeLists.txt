target_include_directories(${EXE_NAME}
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( ${EXE_NAME}
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
   target_sources( ${EXE_NAME}
      PRIVATE
      sdc_hydro.cpp
      )
endif ()
