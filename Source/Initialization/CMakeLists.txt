target_include_directories(${EXE_NAME}
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )


target_sources( ${EXE_NAME}
   PRIVATE
   Nyx_setup.cpp
   Nyx_initdata.cpp
   Nyx_initcosmo.cpp
   bc_fill.cpp
   bc_fill.H
)
