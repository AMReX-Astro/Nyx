target_include_directories(nyxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyxcore
   PRIVATE
   AGN_sources.cpp
   AGN_halos.cpp
   agn_F.H
   )
