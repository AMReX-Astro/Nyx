target_compile_definitions(${EXE_NAME} PUBLIC $<BUILD_INTERFACE:FORCING>)

target_include_directories(${EXE_NAME}
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( ${EXE_NAME}
   PRIVATE
   MersenneTwister.cpp
   Forcing.cpp
   Forcing_init.cpp
   Forcing.H
   )
