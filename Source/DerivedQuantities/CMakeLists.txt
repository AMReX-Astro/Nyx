target_include_directories(nyxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( nyxcore
   PRIVATE
   ParticleDerive.cpp
   ParticleDerive.H
   Derive_forcing.cpp
   )

if  (Nyx_HYDRO)
   target_sources( nyxcore
      PRIVATE
      Derive.cpp
      Derive.H
      )
endif ()
