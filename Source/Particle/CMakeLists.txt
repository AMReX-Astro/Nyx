target_include_directories(nyxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyxcore
   PRIVATE
   NyxParticleContainer.H
   DarkMatterParticleContainer.H
   NyxParticles.cpp
   DarkMatterParticleContainer.cpp
   )


target_sources( nyxcore
   PUBLIC
   $<${_has_agn}:
   ${CMAKE_CURRENT_LIST_DIR}/AGNParticleContainer.H
   ${CMAKE_CURRENT_LIST_DIR}/AGNParticleContainer.cpp
   >
   )


target_sources( nyxcore
   PUBLIC
   $<${_has_neutrino}:
   ${CMAKE_CURRENT_LIST_DIR}/NeutrinoParticleContainer.cpp
   ${CMAKE_CURRENT_LIST_DIR}/NeutrinoParticleContainer.H
   ${CMAKE_CURRENT_LIST_DIR}/NeutrinoParticles_K.H
   >
   )
