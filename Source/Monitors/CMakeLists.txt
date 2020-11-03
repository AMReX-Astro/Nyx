target_include_directories(nyx
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyx
   PRIVATE
   MemInfo.cpp
   MemInfo.H
   )
