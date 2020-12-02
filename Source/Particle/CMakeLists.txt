target_include_directories(${EXE_NAME}
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( ${EXE_NAME}
   PRIVATE
   NyxParticleContainer.H
   DarkMatterParticleContainer.H
   NyxParticles.cpp
   DarkMatterParticleContainer.cpp
   )

if (Nyx_AGN)
   target_sources( ${EXE_NAME}
      PRIVATE
      AGNParticleContainer.H
      AGNParticleContainer.cpp
      )
endif ()

if (Nyx_NEUTRINO_PARTICLES)
   target_sources( ${EXE_NAME}
      PRIVATE
      NeutrinoParticleContainer.cpp
      NeutrinoParticleContainer.H
      NeutrinoParticles_K.H
      )
endif ()
