target_include_directories(nyx
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyx
   PRIVATE
   bc_fill.cpp
   bc_fill.H
   )
