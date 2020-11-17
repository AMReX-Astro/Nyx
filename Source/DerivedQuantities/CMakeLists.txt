target_include_directories(${EXE_NAME}
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

target_sources( ${EXE_NAME}
   PRIVATE
   ParticleDerive.cpp
   ParticleDerive.H
   Derive_forcing.cpp
   )

if  (Nyx_HYDRO)
   target_sources( ${EXE_NAME}
      PRIVATE
      Derive.cpp
      Derive.H
      )
endif ()
