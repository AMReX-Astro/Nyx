target_include_directories(nyx
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyx
   PRIVATE
   AGN_sources.cpp
   AGN_halos.cpp
   agn_F.H
   )
