target_sources( ${EXE_NAME}
   PRIVATE
   Nyx_advance.cpp
   advance_particles.cpp
   )

if (Nyx_HYDRO)
  target_sources( ${EXE_NAME}
     PRIVATE
     Nyx_time_center_sources.cpp
     update_state_with_sources.cpp
     )

  if (Nyx_GRAVITY)
     target_sources( ${EXE_NAME}
        PRIVATE
        Nyx_correct_gsrc.cpp
        )
  endif ()

  if (NOT Nyx_AGN)
     target_sources( ${EXE_NAME}
        PRIVATE
        Nyx_sources.cpp
        )
  endif ()

endif ()
