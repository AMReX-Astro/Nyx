target_include_directories(nyx
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyx
   PRIVATE
   NyxParticleContainer.H
   DarkMatterParticleContainer.H
   NyxParticles.cpp
   DarkMatterParticleContainer.cpp
   )

if (Nyx_AGN)
   target_sources( nyx
      PRIVATE
      AGNParticleContainer.H
      AGNParticleContainer.cpp
      )
endif ()

if (Nyx_NEUTRINO_PARTICLES)
   target_sources( nyx
      PRIVATE
      NeutrinoParticleContainer.cpp
      NeutrinoParticleContainer.H
      NeutrinoParticles_K.H
      )
endif ()
