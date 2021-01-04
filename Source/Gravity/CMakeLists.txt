target_include_directories(nyxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyxcore PRIVATE Gravity.cpp Gravity.H )
