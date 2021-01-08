target_include_directories(nyxcore
   PUBLIC
   $<${_has_agn}:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyxcore
   PUBLIC
   $<${_has_agn}:
   ${CMAKE_CURRENT_LIST_DIR}/AGN_sources.cpp
   ${CMAKE_CURRENT_LIST_DIR}/AGN_halos.cpp
   ${CMAKE_CURRENT_LIST_DIR}/agn_F.H
   >
   )
nyx_set_cuda_sources(
   ${CMAKE_CURRENT_LIST_DIR}/AGN_sources.cpp
   ${CMAKE_CURRENT_LIST_DIR}/AGN_halos.cpp
   )
