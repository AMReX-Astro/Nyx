target_compile_definitions(nyx PUBLIC $<BUILD_INTERFACE:GRAVITY>)

target_include_directories(nyx
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyx PRIVATE Gravity.cpp Gravity.H )
