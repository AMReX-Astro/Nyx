target_include_directories(nyx
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyx
   PRIVATE
   Godunov.cpp
   Hydro.cpp
   trace_ppm.cpp
   Utilities.cpp
   Godunov.H
   Hydro.H
   PLM.H
   PPM.H
   Utilities.H
   Riemann.H
   )
