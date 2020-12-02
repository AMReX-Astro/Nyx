target_include_directories(${EXE_NAME}
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( ${EXE_NAME}
   PRIVATE
   AGN_sources.cpp
   AGN_halos.cpp
   agn_F.H
   )
