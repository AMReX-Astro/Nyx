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

if (Nyx_AGN)
   target_sources( nyxcore
      PRIVATE
      AGNParticleContainer.H
      AGNParticleContainer.cpp
      )
endif ()

if (Nyx_NEUTRINO_PARTICLES)
   target_sources( nyxcore
      PRIVATE
      NeutrinoParticleContainer.cpp
      NeutrinoParticleContainer.H
      NeutrinoParticles_K.H
      )
endif ()
