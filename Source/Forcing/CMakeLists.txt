target_compile_definitions(nyx PUBLIC $<BUILD_INTERFACE:FORCING>)

target_include_directories(nyx
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyx
   PRIVATE
   MersenneTwister.cpp
   Forcing.cpp
   Forcing_init.cpp
   Forcing.H
   )
