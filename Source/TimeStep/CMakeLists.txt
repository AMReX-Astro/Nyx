target_sources( nyxcore
   PRIVATE
   Nyx_advance.cpp
   Nyx_advance_particles.cpp
   )

if (Nyx_HYDRO)
   target_include_directories(nyxcore
      PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
      )

   target_sources( nyxcore
      PRIVATE
      Nyx_enforce_minimum_density.cpp
      Nyx_enforce_minimum_density.H
      Nyx_time_center_sources.cpp
      Nyx_update_state_with_sources.cpp
      )


   if (Nyx_GRAVITY)
      target_sources( nyxcore
         PRIVATE
         Nyx_correct_gsrc.cpp
         )
   endif ()

   if (NOT Nyx_AGN)
      target_sources( nyxcore PRIVATE Nyx_sources.cpp )
   endif ()

endif ()
