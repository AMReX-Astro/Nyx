target_include_directories(${EXE_NAME}
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( ${EXE_NAME}
   PRIVATE
   Nyx_output.cpp
   Nyx_output.H
   write_info.cpp
   )
